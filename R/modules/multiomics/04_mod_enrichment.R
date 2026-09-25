#' Module: Cross-omics pathway enrichment
#'
#' Runs pathway enrichment for each omics layer using DE results, then
#' combines results to identify consistently dysregulated pathways.
#'
#' If pre-computed enrichment results are available from individual pipelines,
#' uses those. Otherwise, runs KEGG enrichment from DE tables.
#'
#' @param enrichment_results Named list of enrichment results per omics (may be NULL/incomplete)
#' @param de_results Named list of DE results per omics
#' @param harmonization_res Harmonization result with MAE and pre-processing data
#' @param config Full config object
#' @param out_dir Output directory
#' @return List with: per_omics, cross_omics, plots
mod_multiomics_enrichment <- function(enrichment_results = NULL,
                                       de_results = NULL,
                                       harmonization_res = NULL,
                                       config, out_dir) {

    message("\n=== Cross-Omics Pathway Enrichment ===\n")

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Cleared here, before any return can be taken -- including the disabled
    # guard immediately below. The report includes this export on file.exists()
    # alone, so a rerun that switches enrichment off would otherwise keep
    # presenting the previous run's scores as current. The invariant is that
    # the file exists only if THIS invocation produced rows, and a guard that
    # returns before the cleanup is exactly how such an invariant is lost.
    write_compound_gsea_export(NULL, out_dir)

    # Check if enrichment is enabled
    if (!isTRUE(config$modes$multiomics$enrichment$run_enrichment)) {
        message("  Cross-omics enrichment disabled in config")
        return(NULL)
    }

    # Build per-omics enrichment results, running enrichment where needed
    per_omics <- build_per_omics_enrichment(
        enrichment_results = enrichment_results,
        de_results = de_results,
        harmonization_res = harmonization_res,
        config = config,
        out_dir = out_dir
    )

    # --- Compound GSEA (metabolomics) ---
    # Scored before the per_omics guard below, and kept out of per_omics
    # entirely. Two separate reasons, both deliberate:
    #
    # Out of per_omics, because per_omics drives the per-layer bar plots, CSVs
    # and the ORA figure, which are ORA's. It reaches the meta-analysis through
    # its own argument instead (rank_tables below), where
    # merge_pathway_pvalues() picks one method per layer -- rank-based where a
    # layer has it -- rather than taking the minimum across methods.
    #
    # Before the guard, because GSEA does not depend on any ORA table. A run
    # where no layer produced an enriched table is precisely a run where ranked
    # evidence may still be informative, and returning early would have thrown
    # it away. It also means the export is cleared on every path, so a run that
    # scores nothing cannot leave the previous run's file behind.
    #
    # It reads the compound-pathway cache the ORA step fills and never fetches
    # one of its own, so it runs after build_per_omics_enrichment().
    compound_gsea <- tryCatch(
        run_compound_gsea_for_contrasts(
            de_results = de_results,
            harmonization_res = harmonization_res,
            config = config,
            out_dir = out_dir
        ),
        error = function(e) {
            message("  Compound GSEA failed: ", e$message)
            NULL
        }
    )
    if (!is.null(write_compound_gsea_export(compound_gsea, out_dir))) {
        message("  Compound GSEA: ", nrow(compound_gsea),
                " scored pathway-contrast rows written")
    }

    if (is.null(per_omics) || length(per_omics) == 0) {
        message("  No omics layers produced enrichment results")

        if (is.null(compound_gsea) || nrow(compound_gsea) == 0) return(NULL)

        # Compound GSEA alone. Reported as what it is -- no ORA tables, no
        # cross-omics analysis, no figures -- rather than dressed up as a
        # partial enrichment result. run_multigsea_plots() reads per_omics and
        # stops on fewer than two layers, so an empty list is the honest shape
        # and the one its own guard already handles.
        message("  Compound GSEA produced results; returning those alone")
        return(list(
            per_omics = list(),
            cross_omics = NULL,
            compound_gsea = compound_gsea,
            plots = list()
        ))
    }

    # Generate per-omics enrichment barplots (always, even with 1 omics)
    per_omics_plots <- list()
    for (om in names(per_omics)) {
        plot_path <- file.path(out_dir, paste0(om, "_top_pathways.png"))
        per_omics_plots[[paste0(om, "_barplot")]] <- plot_path
        png(plot_path, width = 1000, height = 700, res = 120)
        tryCatch({
            plot_per_omics_barplot(per_omics[[om]], om)
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste(om, "barplot failed:", e$message), cex = 1.2)
        })
        dev.off()
        message("  Saved ", om, " enrichment barplot")

        # Save per-omics table
        write.csv(per_omics[[om]],
                  file.path(out_dir, paste0(om, "_enriched_pathways.csv")),
                  row.names = FALSE)
    }

    # Cross-omics comparison requires >= 2 omics
    #
    # Cleared here and not only inside the analysis: this guard means a rerun
    # that drops to one enriched layer never enters that function at all, and
    # the report finds the earlier run's per-collection figures by glob and
    # shows them beside the one-layer results as though they were current.
    .clear_collection_heatmaps(out_dir)
    .clear_cross_lookup_outputs(out_dir)

    # Compound GSEA joins the meta-analysis as the metabolomics layer's
    # rank-based evidence. A metabolomics layer whose ORA found nothing still
    # counts here, which is the case this exists for: no compound passing the
    # DE threshold leaves ORA empty, while GSEA scores every measured compound.
    rank_tables <- if (!is.null(compound_gsea) && nrow(compound_gsea) > 0) {
        list(metabolomics = compound_gsea)
    } else {
        list()
    }
    n_layers <- length(union(names(per_omics), names(rank_tables)))

    cross_omics_enrich <- NULL
    if (n_layers >= 2) {
        cross_omics_enrich <- analyze_cross_omics_enrichment(
            enrichment_results = per_omics,
            config = config,
            out_dir = out_dir,
            rank_tables = rank_tables
        )
        if (!is.null(cross_omics_enrich)) {
            write_cross_omics_enrichment(cross_omics_enrich, out_dir)
        }
    } else {
        message("  Only ", n_layers, " omics with enrichment results; ",
                "skipping cross-omics comparison (need >= 2)")
    }

    # --- Per-contrast enrichment output ---
    # Save per-contrast barplots, tables, and cross-omics comparison in
    # per_contrast/{contrast_name}/ subdirectories so results cannot overwrite.

    # Contrast names are normalized across omics with normalize_contrast_key(),
    # the key the pathview renderer uses too, and discovered from every source
    # that can take part in a per-contrast analysis -- the rank-based tables as
    # well as per_omics -- so a contrast only compound GSEA scored is not
    # silently skipped.
    discovered <- canonicalize_enrichment_contrasts(per_omics, rank_tables)
    per_omics <- discovered$per_omics
    contrast_names <- discovered$contrast_names

    # Lookup files in any existing per-contrast directory, including contrasts
    # this run no longer has: the report renders every directory there, and the
    # loop below only clears the contrasts it visits.
    .clear_stale_contrast_lookups(out_dir)

    if (length(contrast_names) >= 1) {
        message("  Generating per-contrast enrichment output for ",
                length(contrast_names), " contrasts")
        per_contrast_dir <- file.path(out_dir, "per_contrast")
        dir.create(per_contrast_dir, recursive = TRUE, showWarnings = FALSE)

        for (cname in contrast_names) {
            safe_dir <- gsub("[^a-zA-Z0-9._-]", "_", cname)
            contrast_out <- file.path(per_contrast_dir, safe_dir)
            dir.create(contrast_out, recursive = TRUE, showWarnings = FALSE)
            # Same reason as the run-level clear above: the per-contrast call
            # is guarded on >= 2 layers too, so a contrast that drops below it
            # on a rerun would keep showing the previous run's figures.
            .clear_collection_heatmaps(contrast_out)
            .clear_cross_lookup_outputs(contrast_out)

            per_omics_contrast <- list()
            for (om in names(per_omics)) {
                df <- per_omics[[om]]
                if (is.data.frame(df) && "contrast" %in% colnames(df)) {
                    df_c <- df[df$contrast == cname, , drop = FALSE]
                } else {
                    df_c <- df
                }
                if (is.data.frame(df_c) && nrow(df_c) > 0) {
                    per_omics_contrast[[om]] <- df_c

                    plot_path <- file.path(contrast_out, paste0(om, "_top_pathways.png"))
                    png(plot_path, width = 1000, height = 700, res = 120)
                    tryCatch({
                        plot_per_omics_barplot(df_c, paste0(om, " (", gsub("_", " ", cname), ")"))
                    }, error = function(e) {
                        plot.new()
                        text(0.5, 0.5, paste(om, "barplot failed:", e$message), cex = 1.2)
                    })
                    dev.off()

                    write.csv(df_c,
                              file.path(contrast_out, paste0(om, "_enriched_pathways.csv")),
                              row.names = FALSE)
                }
            }

            # The same contrast's compound GSEA rows, matched on the contrast
            # key rather than the raw name: the layers spell one comparison
            # differently, and per_omics was mapped to canonical names above.
            rank_contrast <- lapply(rank_tables, function(df) {
                if (!"contrast" %in% names(df)) return(df)
                df[normalize_contrast_key(df$contrast) == normalize_contrast_key(cname),
                   , drop = FALSE]
            })
            rank_contrast <- Filter(function(df) nrow(df) > 0, rank_contrast)
            contrast_layers <- union(names(per_omics_contrast), names(rank_contrast))

            message("    ", cname, ": ", length(contrast_layers),
                    " omics (", paste(contrast_layers, collapse = ", "), ")")

            if (length(contrast_layers) >= 2) {
                cross_contrast <- tryCatch({
                    analyze_cross_omics_enrichment(
                        enrichment_results = per_omics_contrast,
                        config = config,
                        out_dir = contrast_out,
                        rank_tables = rank_contrast
                    )
                }, error = function(e) {
                    message("    Cross-omics enrichment failed for ", cname, ": ", e$message)
                    NULL
                })
                if (!is.null(cross_contrast)) {
                    write_cross_omics_enrichment(cross_contrast, contrast_out)
                }
            }
        }
        message("  Per-contrast enrichment output saved to: ", per_contrast_dir)
    }

    message("Cross-omics enrichment analysis complete")

    list(
        per_omics = per_omics,
        cross_omics = cross_omics_enrich,
        # Its own slot, beside per_omics rather than inside it.
        compound_gsea = compound_gsea,
        plots = c(per_omics_plots, if (!is.null(cross_omics_enrich)) cross_omics_enrich$plots else list())
    )
}


#' Canonical contrast names across the per-omics and rank-based tables
#'
#' The layers spell one comparison differently, so contrasts are matched on
#' \code{normalize_contrast_key()} and each key is given one canonical name --
#' the first raw spelling met, per_omics first. Discovery reads the usable
#' rank-based tables as well as per_omics: a contrast that only a rank-based
#' source scored must still get its per-contrast analysis, where that source's
#' rows are matched on the same key.
#'
#' @param per_omics Named list of per-omics enrichment data frames.
#' @param rank_tables Named list of rank-based tables, as passed to
#'   \code{analyze_cross_omics_enrichment()}; only their usable rows (see
#'   \code{.usable_rank_tables()}) are read.
#' @return List with \code{per_omics}, its \code{contrast} columns rewritten to
#'   the canonical names, and \code{contrast_names}, the canonical name of
#'   every contrast found in either source.
#' @examples
#' po <- list(proteomics = data.frame(contrast = "A_vs_B", pvalue = 0.01))
#' rk <- list(metabolomics = data.frame(contrast = c("A vs B", "C_vs_D"),
#'                                      method = "fgsea", pval = 0.02))
#' canonicalize_enrichment_contrasts(po, rk)$contrast_names  # "A_vs_B" "C_vs_D"
canonicalize_enrichment_contrasts <- function(per_omics, rank_tables = list()) {
    raw_of <- function(tables) unlist(lapply(tables, function(df) {
        if (is.data.frame(df) && "contrast" %in% colnames(df)) as.character(df$contrast)
        else NULL
    }), use.names = FALSE)

    all_raw_contrasts <- unique(c(raw_of(per_omics),
                                  raw_of(.usable_rank_tables(rank_tables))))
    all_raw_contrasts <- all_raw_contrasts[!is.na(all_raw_contrasts)]
    if (length(all_raw_contrasts) == 0) {
        return(list(per_omics = per_omics, contrast_names = character(0)))
    }

    norm_keys <- normalize_contrast_key(all_raw_contrasts)
    # For each normalized key, the first raw name met is canonical
    canonical <- all_raw_contrasts[!duplicated(norm_keys)]
    names(canonical) <- norm_keys[!duplicated(norm_keys)]
    raw_to_canonical <- stats::setNames(unname(canonical[norm_keys]), all_raw_contrasts)

    for (om in names(per_omics)) {
        df <- per_omics[[om]]
        if (is.data.frame(df) && "contrast" %in% colnames(df)) {
            per_omics[[om]]$contrast <- unname(raw_to_canonical[as.character(df$contrast)])
        }
    }

    list(per_omics = per_omics, contrast_names = unname(canonical))
}
