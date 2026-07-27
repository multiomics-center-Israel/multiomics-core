# R/modules/metabolomics/06_mod_report.R
#
# Metabolomics HTML report module.
#
# Functions:
#   mod_met_qc_summary_report()  — standalone intermediate normalization QC report
#   mod_metabolomics_report()    — full analysis report (includes QC section via child Rmd)
#   mod_metabolome_overview()    — metabolome overview report (detection, pathways, catalog)


#' Render standalone Normalization QC summary report
#'
#' Renders \code{qc_summary_report.Rmd} to a self-contained HTML at
#' \code{\{out_dir\}/qc/normalization_review_report.html}.  This report is
#' intended for intermediate review: open it after \code{tar_make()} completes
#' the QC layer to decide which normalization method to use, before looking at
#' the full analysis report.
#'
#' @param qc_comparison_file Character scalar.  Path to
#'   \code{normalization_qc_benchmark.tsv} (return value of
#'   \code{met_qc_comparison}).
#' @param qc_suite_files     Character vector.  All QC file paths from
#'   \code{met_log_qc}, \code{met_norm_tss_qc}, \code{met_norm_median_qc},
#'   and \code{met_norm_pqn_qc} combined.
#' @param config             Full pipeline config list.
#' @param out_dir            Mode output directory (\code{metab_out_dir}).
#' @return Character scalar: absolute path to the rendered HTML file.
mod_met_qc_summary_report <- function(qc_comparison_file, qc_suite_files,
                                      config, out_dir) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        warning("rmarkdown not available -- skipping QC summary report generation")
        return(character(0))
    }

    template <- file.path("R", "pipeline", "metabolomics", "templates",
                          "qc_summary_report.Rmd")
    if (!file.exists(template)) {
        stop("QC summary report template not found at: ", template)
    }

    qc_dir   <- file.path(out_dir, "qc")
    ensure_dir(qc_dir)
    proj_name <- config$project$analysis_round
    out_file <- file.path(qc_dir, paste0("normalization_review_report_", proj_name, ".html"))

    message("Rendering normalization QC summary report -> ", out_file)

    rmarkdown::render(
        input       = template,
        output_file = basename(out_file),
        output_dir  = dirname(out_file),
        params = list(
            qc_comparison_file = qc_comparison_file,
            qc_suite_files     = qc_suite_files,
            config             = config
        ),
        envir  = new.env(parent = globalenv()),
        quiet  = TRUE
    )

    out_file
}


#' Render metabolomics HTML report
#'
#' Locates the Rmd template and calls rmarkdown::render() with all module
#' results passed as params.
#'
#' @param pre             List from preprocess_metabolomics().
#' @param qc_res          List from mod_metabolomics_qc_pre().
#' @param de_res          List from mod_metabolomics_de().
#' @param clustering_res  List from mod_metabolomics_clustering() (or NULL).
#' @param feature_sel_res List from mod_metabolomics_feature_selection() (or NULL).
#' @param enrichment_res  List from mod_metabolomics_enrichment() (or NULL).
#' @param config          Full pipeline config.
#' @param out_dir         Output directory for this mode.
#' @return Character path to the rendered HTML file.
mod_metabolomics_report <- function(pre, qc_res, de_res,
                                    clustering_res = NULL,
                                    feature_sel_res,
                                    enrichment_res,
                                    network_res    = NULL,
                                    config, out_dir,
                                    qc_comparison_file = NULL,
                                    qc_suite_files     = NULL,
                                    commentary_file    = NULL,
                                    mummichog_pathways = NULL) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        warning("rmarkdown not available -- skipping report generation")
        return(character(0))
    }

    # Locate template relative to project root (targets sets cwd to project root)
    template <- file.path("R", "pipeline", "metabolomics", "templates",
                          "report_metabolomics.Rmd")

    if (!file.exists(template)) {
        stop("Report template not found at: ", template)
    }

    # Render directly to Results root (parent of mode dir) — avoids duplicate
    # inside the metabolomics/ folder.  Matches RNA/proteomics placement.
    results_root <- dirname(out_dir)  # e.g. Results_project_A01/
    ensure_dir(results_root)

    # Copy Rmd to Results root so users can re-render
    dest_rmd <- file.path(results_root, "report_metabolomics.Rmd")
    file.copy(template, dest_rmd, overwrite = TRUE)

    out_file <- file.path(results_root, "report_metabolomics.html")
    message("Rendering metabolomics report -> ", out_file)

    # Attach module-level plots to rf/plsda sub-results so the Rmd template
    # can access them as rf_res$plots$rf_importance, plsda_res$plots$plsda_scores, etc.
    rf_res_out    <- NULL
    plsda_res_out <- NULL
    if (!is.null(feature_sel_res)) {
        fs_plots <- feature_sel_res$plots %||% list()
        if (!is.null(feature_sel_res$rf)) {
            rf_res_out <- feature_sel_res$rf
            rf_res_out$plots <- fs_plots[grep("^rf_", names(fs_plots))]
        }
        if (!is.null(feature_sel_res$plsda)) {
            plsda_res_out <- feature_sel_res$plsda
            plsda_res_out$plots <- fs_plots[grep("^plsda_", names(fs_plots))]
        }
    }

    # Build per-contrast mummichog report sections (06e) from the per-contrast
    # pathway tables. Empty list when mummichog is disabled or no contrast has a
    # result, which keeps the report section hidden.
    mummi_sections <- build_mummichog_report_sections(mummichog_pathways, config)

    # Save standalone presentation exports per contrast (PNG + PDF plot, TSV +
    # CSV table) into mummichog_pinned/, named by the section's de-duplicated
    # slug (so labels that collapse to the same token don't overwrite each
    # other), in addition to embedding them in the HTML. The engine's mcg_*
    # result files are untouched. Paths are returned so the file target tracks them.
    mummi_exports <- character(0)
    for (contrast in names(mummi_sections)) {
        s <- mummi_sections[[contrast]]
        mummi_exports <- c(mummi_exports,
                           save_mummichog_exports(s$plot, s$table, out_dir,
                                                  contrast_label = s$slug))
    }

    render_params <- list(

        pre                = pre,
        qc_res             = qc_res,
        de_res             = de_res,
        clustering_res     = clustering_res,
        rf_res             = rf_res_out,
        plsda_res          = plsda_res_out,
        enrichment_res     = enrichment_res,
        network_res        = network_res,
        config             = config,
        qc_comparison_file = qc_comparison_file,
        qc_suite_files     = qc_suite_files,
        commentary_file    = commentary_file,
        mummichog_sections = mummi_sections
    )


    rmarkdown::render(
        input       = dest_rmd,
        output_file = basename(out_file),
        output_dir  = results_root,
        params      = render_params,
        envir       = new.env(parent = globalenv()),
        quiet       = TRUE
    )

    # Also render the styled report (matching RNA/proteomics look)
    # Place it at the Results root (parent of mode dir), with Rmd alongside
    styled_template <- file.path("R", "pipeline", "metabolomics", "templates",
                                  "report_metabolomics.Rmd")
    if (file.exists(styled_template)) {
        results_root <- dirname(out_dir)  # e.g. Results_project_A01/
        # Copy Rmd to results root (like RNA/proteomics do)
        dest_rmd <- file.path(results_root, "report_metabolomics.Rmd")
        file.copy(styled_template, dest_rmd, overwrite = TRUE)
        message("Rendering styled metabolomics report -> ",
                file.path(results_root, "report_metabolomics.html"))
        tryCatch({
            rmarkdown::render(
                input       = dest_rmd,
                output_file = "report_metabolomics.html",
                output_dir  = results_root,
                params      = render_params,
                envir       = new.env(parent = globalenv()),
                quiet       = TRUE
            )
        }, error = function(e) {
            warning("Styled report rendering failed: ", e$message,
                    "\nOriginal report was generated successfully.")
        })
    }

    unique(c(out_file, mummi_exports))
}


#' Render metabolome overview report
#'
#' Exports pipeline data to temporary files, then renders the standalone
#' \code{metabolome_overview_report.Rmd} into \code{\{out_dir\}/Metabolome_overview/}.
#'
#' @param pre     List from preprocessing (expr_raw, row_data, meta).
#' @param inputs  Raw inputs from \code{load_metabolomics_inputs()}.
#' @param config  Full pipeline config.
#' @param out_dir Mode output directory (metab_out_dir).
#' @return Character path to the rendered HTML file (or \code{character(0)} on skip).
mod_metabolome_overview <- function(pre, inputs, config, out_dir) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        warning("rmarkdown not available -- skipping metabolome overview")
        return(character(0))
    }

    template <- file.path("scripts", "metabolome_overview_report.Rmd")
    if (!file.exists(template)) {
        warning("Metabolome overview template not found at: ", template,
                " -- skipping")
        return(character(0))
    }

    overview_dir <- file.path(out_dir, "Metabolome_overview")
    ensure_dir(overview_dir)

    metab_cfg <- config$modes$metabolomics

    # ---- Export data files the Rmd expects ----

    # 1. Full feature table (annotations + sample intensities)
    #    Use expr_filt (filtered) if available and row_data is filtered;
    #    otherwise fall back to expr_raw with row_data subsetting.
    data_file <- file.path(overview_dir, "all_features.tsv")
    rd <- pre$row_data
    expr_for_export <- if (!is.null(pre$expr_filt) && nrow(rd) == nrow(pre$expr_filt)) {
        pre$expr_filt
    } else if (nrow(rd) == nrow(pre$expr_raw)) {
        pre$expr_raw
    } else {
        # row_data is filtered but expr_raw is not — subset expr_raw
        common <- intersect(rownames(rd), rownames(pre$expr_raw))
        pre$expr_raw[common, , drop = FALSE]
    }
    full_df <- cbind(rd[rownames(expr_for_export), , drop = FALSE],
                     as.data.frame(expr_for_export))
    utils::write.table(full_df, data_file, sep = "\t", row.names = FALSE,
                        quote = FALSE)

    # 2. Level-1 subset (if identification_level column exists)
    level1_file <- NULL
    if ("identification_level" %in% colnames(pre$row_data)) {
        l1_mask <- !is.na(pre$row_data$identification_level) &
                   pre$row_data$identification_level == 1L
        if (any(l1_mask)) {
            level1_file <- file.path(overview_dir, "level1_features.tsv")
            l1_df <- cbind(pre$row_data[l1_mask, ],
                           as.data.frame(pre$expr_raw[l1_mask, ]))
            utils::write.table(l1_df, level1_file, sep = "\t",
                                row.names = FALSE, quote = FALSE)
        }
    }

    # 3. Metadata
    meta_file <- file.path(overview_dir, "metadata.tsv")
    utils::write.table(pre$meta, meta_file, sep = "\t", row.names = FALSE,
                        quote = FALSE)

    # 4. KEGG pathway mapping (build from row_data KEGG column if available)
    kegg_file <- NULL
    if ("KEGG" %in% colnames(pre$row_data)) {
        kegg_ids <- pre$row_data$KEGG[!is.na(pre$row_data$KEGG) &
                                       grepl("^C[0-9]+$", pre$row_data$KEGG)]
        if (length(kegg_ids) > 0) {
            kegg_map <- tryCatch(
                .build_kegg_mapping(kegg_ids),
                error = function(e) {
                    message("KEGG pathway mapping failed: ", e$message)
                    NULL
                }
            )
            if (!is.null(kegg_map)) {
                kegg_file <- file.path(overview_dir, "kegg_mapping.rds")
                saveRDS(kegg_map, kegg_file)
            }
        }
    }

    # ---- Resolve column names from config ----
    col_cfg <- metab_cfg$columns %||% list()

    out_file <- file.path(overview_dir, "metabolome_overview.html")
    message("Rendering metabolome overview -> ", out_file)

    tryCatch({
        rmarkdown::render(
            input       = template,
            output_file = basename(out_file),
            output_dir  = overview_dir,
            params = list(
                data_file        = data_file,
                level1_file      = level1_file,
                metadata_file    = meta_file,
                kegg_mapping_file = kegg_file,
                project_name     = config$project$name %||% "Metabolomics",
                analyst          = config$project$analyst %||% "Bioinformatics Core",
                feature_id_col   = col_cfg$feature_id %||% "feature_id",
                mz_col           = col_cfg$mz %||% "m/z",
                rt_col           = col_cfg$rt %||% "RT [min]",
                hmdb_col         = col_cfg$hmdb %||% "HMDB",
                kegg_col         = col_cfg$kegg %||% "KEGG",
                formula_col      = col_cfg$formula %||% "Formula"
            ),
            envir = new.env(parent = globalenv()),
            quiet = TRUE
        )
    }, error = function(e) {
        warning("Metabolome overview rendering failed: ", e$message)
        return(character(0))
    })
    out_file
}


# Internal: build KEGG pathway mapping from a vector of compound IDs
.build_kegg_mapping <- function(kegg_ids) {
    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        message("KEGGREST not available -- skipping KEGG pathway mapping")
        return(NULL)
    }

    kegg_ids <- unique(kegg_ids)
    compound_pathways <- list()
    pathway_names     <- character(0)

    batch_size <- 10
    for (i in seq(1, length(kegg_ids), by = batch_size)) {
        batch <- kegg_ids[i:min(i + batch_size - 1, length(kegg_ids))]
        tryCatch({
            links <- KEGGREST::keggLink("pathway", paste0("cpd:", batch))
            if (length(links) > 0) {
                cpd_ids <- sub("^cpd:", "", names(links))
                pw_ids  <- sub("^path:", "", unname(links))
                for (j in seq_along(cpd_ids)) {
                    cid <- cpd_ids[j]
                    pid <- pw_ids[j]
                    compound_pathways[[cid]] <- c(compound_pathways[[cid]], pid)
                }
            }
            Sys.sleep(0.3)
        }, error = function(e) {
            message("  KEGG pathway query failed for batch: ", e$message)
        })
    }

    # Deduplicate per-compound pathway lists
    compound_pathways <- lapply(compound_pathways, unique)

    # Count pathways
    all_pws <- unlist(compound_pathways, use.names = FALSE)
    pathway_counts <- sort(table(all_pws), decreasing = TRUE)

    # Fetch pathway names
    pw_unique <- names(pathway_counts)
    for (i in seq(1, length(pw_unique), by = batch_size)) {
        batch <- pw_unique[i:min(i + batch_size - 1, length(pw_unique))]
        tryCatch({
            info <- KEGGREST::keggList("pathway", batch)
            pathway_names <- c(pathway_names, setNames(
                sub(" - .*$", "", unname(info)),
                sub("^path:", "", names(info))
            ))
            Sys.sleep(0.3)
        }, error = function(e) NULL)
    }

    list(
        compound_pathways = compound_pathways,
        pathway_counts    = pathway_counts,
        pathway_names     = pathway_names
    )
}
