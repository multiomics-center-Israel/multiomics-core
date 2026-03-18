# R/modules/metabolomics/06_mod_report.R
#
# Metabolomics HTML report module.
#
# Functions:
#   mod_met_qc_summary_report()  — standalone intermediate normalization QC report
#   mod_metabolomics_report()    — full analysis report (includes QC section via child Rmd)


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
#' @param qc_comparison_file Character scalar. Path to normalization_qc_benchmark.tsv (or NULL).
#' @param qc_suite_files     Character vector. QC file paths (or NULL).
#' @return Character path to the rendered HTML file.
mod_metabolomics_report <- function(pre, qc_res, de_res,
                                    clustering_res = NULL,
                                    feature_sel_res,
                                    enrichment_res,
                                    config, out_dir,
                                    qc_comparison_file = NULL,
                                    qc_suite_files     = NULL) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        warning("rmarkdown not available -- skipping report generation")
        return(character(0))
    }

    # Locate template relative to project root (targets sets cwd to project root)
    template <- file.path("R", "pipeline", "metabolomics", "templates",
                          "metabolomics_report.Rmd")

    if (!file.exists(template)) {
        stop("Report template not found at: ", template)
    }

    out_file <- file.path(out_dir, "metabolomics_report.html")
    ensure_dir(dirname(out_file))

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

    render_params <- list(
        pre                = pre,
        qc_res             = qc_res,
        de_res             = de_res,
        clustering_res     = clustering_res,
        rf_res             = rf_res_out,
        plsda_res          = plsda_res_out,
        enrichment_res     = enrichment_res,
        config             = config,
        qc_comparison_file = qc_comparison_file,
        qc_suite_files     = qc_suite_files
    )

    rmarkdown::render(
        input       = template,
        output_file = basename(out_file),
        output_dir  = dirname(out_file),
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

    out_file
}
