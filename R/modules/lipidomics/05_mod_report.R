# R/modules/lipidomics/05_mod_report.R
#
# Lipidomics HTML report module: renders the styled report (matching the
# RNA/proteomics/metabolomics look).


#' Render lipidomics HTML report
#'
#' @param pre             List from preprocess_lipidomics().
#' @param qc_res          List from mod_lipidomics_qc_pre().
#' @param de_res          List from mod_lipidomics_de().
#' @param feature_sel_res List from mod_lipidomics_feature_selection() (or NULL).
#' @param class_res       List from mod_lipidomics_class_analysis() (or NULL).
#' @param config          Full pipeline config.
#' @param out_dir         Output directory for this mode.
#' @return Character path to the rendered HTML file.
mod_lipidomics_report <- function(pre, qc_res, de_res, feature_sel_res,
                                   class_res, config, out_dir,
                                   clustering_res = NULL,
                                   biomarker_res = NULL, pathway_res = NULL,
                                   qc_enhanced = NULL,
                                   pool_cv = NULL,
                                   commentary_file = NULL) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        warning("rmarkdown not available -- skipping report generation")
        return(character(0))
    }

    template <- file.path("R", "pipeline", "lipidomics", "templates",
                          "report_lipidomics.Rmd")

    if (!file.exists(template)) {
        stop("Report template not found at: ", template)
    }

    # Attach module-level plots to sub-results
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
        pre             = pre,
        qc_res          = qc_res,
        de_res          = de_res,
        rf_res          = rf_res_out,
        plsda_res       = plsda_res_out,
        class_res       = class_res,
        clustering_res  = clustering_res,
        biomarker_res   = biomarker_res,
        pathway_res     = pathway_res,
        qc_enhanced     = qc_enhanced,
        pool_cv         = pool_cv,
        config          = config,
        commentary_file = commentary_file
    )

    # Render at the Results root (parent of the mode dir), with the Rmd copied
    # alongside the HTML (matching the other omics' report layout).
    results_root <- dirname(out_dir)
    out_file <- file.path(results_root, "report_lipidomics.html")
    dest_rmd <- file.path(results_root, "report_lipidomics.Rmd")
    ensure_dir(results_root)
    file.copy(template, dest_rmd, overwrite = TRUE)

    message("Rendering lipidomics report -> ", out_file)

    rmarkdown::render(
        input       = dest_rmd,
        output_file = basename(out_file),
        output_dir  = results_root,
        params      = render_params,
        envir       = new.env(parent = globalenv()),
        quiet       = TRUE
    )

    out_file
}
