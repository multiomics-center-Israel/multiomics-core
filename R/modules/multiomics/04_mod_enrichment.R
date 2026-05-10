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

    if (is.null(per_omics) || length(per_omics) == 0) {
        message("  No omics layers produced enrichment results")
        return(NULL)
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
    cross_omics_enrich <- NULL
    if (length(per_omics) >= 2) {
        cross_omics_enrich <- analyze_cross_omics_enrichment(
            enrichment_results = per_omics,
            config = config,
            out_dir = out_dir
        )
        if (!is.null(cross_omics_enrich)) {
            write_cross_omics_enrichment(cross_omics_enrich, out_dir)
        }
    } else {
        message("  Only ", length(per_omics), " omics with enrichment results; ",
                "skipping cross-omics comparison (need >= 2)")
    }

    # --- Per-contrast enrichment output ---
    # Save per-contrast barplots, tables, and cross-omics comparison in
    # per_contrast/{contrast_name}/ subdirectories so results cannot overwrite.
    contrast_names <- unique(unlist(lapply(per_omics, function(df) {
        if (is.data.frame(df) && "contrast" %in% colnames(df)) df$contrast
        else NULL
    })))

    if (length(contrast_names) > 1) {
        message("  Generating per-contrast enrichment output for ",
                length(contrast_names), " contrasts")
        per_contrast_dir <- file.path(out_dir, "per_contrast")
        dir.create(per_contrast_dir, recursive = TRUE, showWarnings = FALSE)

        for (cname in contrast_names) {
            safe_dir <- gsub("[^a-zA-Z0-9._-]", "_", cname)
            contrast_out <- file.path(per_contrast_dir, safe_dir)
            dir.create(contrast_out, recursive = TRUE, showWarnings = FALSE)

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

            message("    ", cname, ": ", length(per_omics_contrast),
                    " omics (", paste(names(per_omics_contrast), collapse = ", "), ")")

            if (length(per_omics_contrast) >= 2) {
                cross_contrast <- tryCatch({
                    analyze_cross_omics_enrichment(
                        enrichment_results = per_omics_contrast,
                        config = config,
                        out_dir = contrast_out
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
        plots = c(per_omics_plots, if (!is.null(cross_omics_enrich)) cross_omics_enrich$plots else list())
    )
}
