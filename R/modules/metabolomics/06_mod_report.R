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
                                    config, out_dir) {
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

    rmarkdown::render(
        input       = template,
        output_file = basename(out_file),
        output_dir  = dirname(out_file),
        params = list(
            pre            = pre,
            qc_res         = qc_res,
            de_res         = de_res,
            rf_res         = if (!is.null(feature_sel_res)) feature_sel_res$rf else NULL,
            plsda_res      = if (!is.null(feature_sel_res)) feature_sel_res$plsda else NULL,
            enrichment_res = enrichment_res,
            config         = config
        ),
        envir  = new.env(parent = globalenv()),
        quiet  = TRUE
    )

    out_file
}
