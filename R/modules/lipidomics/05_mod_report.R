# R/modules/lipidomics/05_mod_report.R
#
# Lipidomics HTML report module: renders both the basic and styled reports.


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
                                   biomarker_res = NULL, pathway_res = NULL,
                                   qc_enhanced = NULL,
                                   commentary_file = NULL) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        warning("rmarkdown not available -- skipping report generation")
        return(character(0))
    }

    template <- file.path("R", "pipeline", "lipidomics", "templates",
                          "lipidomics_report.Rmd")

    if (!file.exists(template)) {
        stop("Report template not found at: ", template)
    }

    out_file <- file.path(out_dir, "lipidomics_report.html")
    ensure_dir(dirname(out_file))

    message("Rendering lipidomics report -> ", out_file)

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
        biomarker_res   = biomarker_res,
        pathway_res     = pathway_res,
        qc_enhanced     = qc_enhanced,
        config          = config,
        commentary_file = commentary_file
    )

    rmarkdown::render(
        input       = template,
        output_file = basename(out_file),
        output_dir  = dirname(out_file),
        params      = render_params,
        envir       = new.env(parent = globalenv()),
        quiet       = TRUE
    )

    # Also render the styled report (matching RNA/proteomics/metabolomics look)
    # Place it at the Results root (parent of mode dir), with Rmd alongside
    styled_template <- file.path("R", "pipeline", "lipidomics", "templates",
                                  "report_lipidomics.Rmd")
    if (file.exists(styled_template)) {
        results_root <- dirname(out_dir)
        dest_rmd <- file.path(results_root, "report_lipidomics.Rmd")
        file.copy(styled_template, dest_rmd, overwrite = TRUE)
        message("Rendering styled lipidomics report -> ",
                file.path(results_root, "report_lipidomics.html"))
        tryCatch({
            rmarkdown::render(
                input       = dest_rmd,
                output_file = "report_lipidomics.html",
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
