# R/modules/lipidomics/11_mod_lipid_links.R
#
# Target wrapper for the lipid-links extension. Glues together:
#   13_links_inputs.R         (load + harmonize external DE tables)
#   14_links_pathway_map.R    (lipid class → pathway curated map)
#   15_links_kegg_fetch.R     (KEGG pathway member fetch w/ cache)
#   16_links_scoring.R        (pathway-based scoring)
#   17_links_pubmed.R         (literature linker)
#   18_links_combine.R        (merge pathway + literature)
#   19_links_report.R         (HTML report)
#
# Returns the list of output file paths so the target can be declared as
# `format = "file"`.


mod_lipid_links <- function(pre, lipid_de_res, config, out_dir) {
    cfg <- config$modes$lipidomics %||% list()
    links_cfg <- cfg$links %||% list()
    if (!isTRUE(links_cfg$enabled)) {
        message("[links] disabled in config — skipping")
        return(character(0))
    }

    # ---- output layout ----
    base <- file.path(out_dir, "links")
    tables_dir  <- file.path(base, "tables")
    cache_dir   <- file.path(base, "cache")
    kegg_cache  <- file.path(cache_dir, "kegg")
    pubmed_cache <- file.path(cache_dir, "pubmed")
    ensure_dir(base); ensure_dir(tables_dir)

    # ---- 1. Harmonized DE inputs ----
    lipid_long <- harmonize_lipid_de(lipid_de_res, pre$row_data)
    ext_long_list <- load_external_omics_for_links(links_cfg)

    if (!is.null(lipid_long)) {
        save_tsv(lipid_long, tables_dir, "harmonized_de.lipidomics.tsv")
    }
    for (lay in names(ext_long_list)) {
        save_tsv(ext_long_list[[lay]], tables_dir,
                  paste0("harmonized_de.", lay, ".tsv"))
    }

    # ---- 2. Lipid class -> pathway map ----
    pmap <- load_lipid_class_pathway_map()
    save_tsv(pmap, tables_dir, "lipid_class_pathways.tsv")
    organism <- links_cfg$pathways$organism %||% "cel"
    fdr_lipid <- links_cfg$output$fdr_cutoff_lipid %||% 0.05
    fdr_layer <- links_cfg$output$fdr_cutoff_layer %||% 0.25
    abs_lfc_min <- links_cfg$output$abs_logfc_min %||% 0.1

    lipid_pw_long <- lipids_to_pathways(lipid_long, pmap,
                                          fdr_cutoff = fdr_lipid,
                                          organism = organism)
    if (!is.null(lipid_pw_long)) {
        save_tsv(lipid_pw_long, tables_dir, "lipid_pathway_expansion.tsv")
    }

    # ---- 3. KEGG pathway membership ----
    pw_members <- NULL
    if (!is.null(lipid_pw_long) && nrow(lipid_pw_long) > 0) {
        pw_members <- fetch_kegg_pathway_members(
            pathway_ids = unique(lipid_pw_long$kegg_pathway_id),
            cache_dir   = kegg_cache,
            organism    = organism
        )
        if (!is.null(pw_members) && nrow(pw_members) > 0) {
            save_tsv(pw_members, tables_dir, "kegg_pathway_members.tsv")
        }
    }

    # ---- 4. Pathway-based scoring per layer ----
    pathway_links_list <- NULL
    if (!is.null(lipid_pw_long) && !is.null(pw_members) &&
        nrow(pw_members) > 0) {
        pathway_links_list <- score_pathway_links(
            lipid_pw_long = lipid_pw_long,
            ext_long_list = ext_long_list,
            pw_members    = pw_members,
            fdr_cutoff_layer = fdr_layer,
            abs_logfc_min    = abs_lfc_min
        )
        for (lay in names(pathway_links_list)) {
            df <- pathway_links_list[[lay]]
            if (!is.null(df) && nrow(df) > 0) {
                save_tsv(df, tables_dir,
                          paste0("pathway_links.", lay, ".tsv"))
            }
        }
    }

    # ---- 5. Literature linker ----
    literature_res <- NULL
    lit_cfg <- links_cfg$literature %||% list()
    if (isTRUE(lit_cfg$enabled) && !is.null(pathway_links_list)) {
        jif_tbl <- load_journal_if_table(min_if = lit_cfg$min_journal_if %||% 15)
        save_tsv(jif_tbl[, c("journal_full","journal_iso","impact_factor")],
                  tables_dir, "journal_if_used.tsv")

        # Build a top-N pair pool per layer per contrast for querying
        top_n <- links_cfg$output$top_n_per_layer %||% 50
        candidate_rows <- list()
        for (lay in names(pathway_links_list)) {
            df <- pathway_links_list[[lay]]
            if (is.null(df) || nrow(df) == 0) next
            # Take top-N pairs by score per contrast
            df_sp <- split(df, df$contrast)
            top_pairs <- do.call(rbind, lapply(df_sp, function(d) {
                d <- d[order(-d$pathway_link_score), , drop = FALSE]
                d[seq_len(min(top_n, nrow(d))), , drop = FALSE]
            }))
            candidate_rows[[lay]] <- unique(
                top_pairs[, c("layer", "contrast", "lipid_class",
                               "feature_name"), drop = FALSE]
            )
        }
        candidates <- do.call(rbind, candidate_rows)
        if (!is.null(candidates) && nrow(candidates) > 0) {
            literature_res <- run_pubmed_literature_links(
                candidate_pairs = candidates,
                journal_if_tbl  = jif_tbl,
                cache_dir       = pubmed_cache,
                year_min        = lit_cfg$year_min %||% 2021,
                api_key         = lit_cfg$ncbi_api_key,
                context_mode    = lit_cfg$context %||% "same_sentence",
                max_papers_per_pair = lit_cfg$max_papers_per_pair %||% 50,
                sleep_seconds   = lit_cfg$sleep_seconds %||% 0.34
            )
            if (!is.null(literature_res$papers)) {
                save_tsv(literature_res$papers, tables_dir,
                          "literature_links.tsv")
            }
            if (!is.null(literature_res$summary)) {
                save_tsv(literature_res$summary, tables_dir,
                          "literature_links_summary.tsv")
            }
        }
    }

    # ---- 6. Combine + final per-layer rankings ----
    combined_list <- combine_lipid_links(pathway_links_list, literature_res)
    if (!is.null(combined_list)) {
        for (lay in names(combined_list)) {
            df <- combined_list[[lay]]
            if (!is.null(df) && nrow(df) > 0) {
                save_tsv(df, tables_dir,
                          paste0("combined_links.", lay, ".tsv"))
            }
        }
    }

    # ---- 7. Enzyme-lipid network figures (LipidONE-style) ----
    fig_files <- tryCatch(
        write_enzyme_lipid_figures(
            lipid_long    = lipid_long,
            ext_long_list = ext_long_list,
            out_dir       = base,
            organism      = organism
        ),
        error = function(e) {
            message("[links/viz] enzyme-lipid figures failed: ",
                    conditionMessage(e))
            character(0)
        }
    )

    # ---- 8. Report ----
    report_files <- render_lipid_links_report(
        combined_list      = combined_list,
        pathway_links_list = pathway_links_list,
        literature_res     = literature_res,
        lipid_long         = lipid_long,
        ext_long_list      = ext_long_list,
        links_cfg          = links_cfg,
        out_dir            = base,
        enzyme_lipid_figs  = fig_files
    )

    written <- list.files(tables_dir, full.names = TRUE)
    fig_dir <- file.path(base, "figures")
    if (dir.exists(fig_dir)) {
        written <- c(written, list.files(fig_dir, full.names = TRUE))
    }
    c(written, report_files)
}
