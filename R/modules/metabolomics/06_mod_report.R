# R/modules/metabolomics/06_mod_report.R
#
# Metabolomics HTML report module: renders the parameterized Rmd template
# with all module results and returns the output file path for targets
# format = "file" tracking.


#' Render metabolomics HTML report
#'
#' Locates the Rmd template and calls rmarkdown::render() with all module
#' results passed as params.
#'
#' @param pre             List from preprocess_metabolomics().
#' @param qc_res          List from mod_metabolomics_qc_pre().
#' @param de_res          List from mod_metabolomics_de().
#' @param feature_sel_res List from mod_metabolomics_feature_selection() (or NULL).
#' @param enrichment_res  List from mod_metabolomics_enrichment() (or NULL).
#' @param config          Full pipeline config.
#' @param out_dir         Output directory for this mode.
#' @return Character path to the rendered HTML file.
mod_metabolomics_report <- function(pre, qc_res, de_res, feature_sel_res,
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

    rmarkdown::render(
        input       = template,
        output_file = basename(out_file),
        output_dir  = dirname(out_file),
        params = list(
            pre                = pre,
            qc_res             = qc_res,
            de_res             = de_res,
            rf_res             = rf_res_out,
            plsda_res          = plsda_res_out,
            enrichment_res     = enrichment_res,
            config             = config,
            qc_comparison_file = qc_comparison_file,
            qc_suite_files     = qc_suite_files
        ),
        envir  = new.env(parent = globalenv()),
        quiet  = TRUE
    )

    out_file
}
