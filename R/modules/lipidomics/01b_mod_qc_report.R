# R/modules/lipidomics/01b_mod_qc_report.R
#
# Standalone lipidomics QC summary report. Renders an HTML pulling together
# artifacts already produced by Stage-1 / Stage-3 QC targets:
#   - normalization_method_comparison.tsv (from evaluate_methods)
#   - missingness summary (from pre$info$missingness)
#   - intensity boxplots / densities, PCA, sample-distance heatmap
#   - PERMANOVA (when vegan is available)

mod_lipidomics_qc_report <- function(pre, qc_res, qc_enhanced, config, out_dir) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        warning("rmarkdown not available -- skipping lipidomics QC report")
        return(character(0))
    }

    template <- file.path("R", "pipeline", "lipidomics", "templates",
                          "lipidomics_qc_report.Rmd")
    if (!file.exists(template)) {
        stop("Lipidomics QC report template not found at: ", template)
    }

    qc_dir <- file.path(out_dir, "qc")
    ensure_dir(qc_dir)
    proj_round <- config$project$analysis_round %||% "L01"
    out_file <- file.path(qc_dir,
                          paste0("lipidomics_qc_report_", proj_round, ".html"))

    message("Rendering lipidomics QC report -> ", out_file)

    rmarkdown::render(
        input       = template,
        output_file = basename(out_file),
        output_dir  = dirname(out_file),
        params = list(
            pre         = pre,
            qc_res      = qc_res,
            qc_enhanced = qc_enhanced,
            config      = config,
            out_dir     = out_dir
        ),
        envir = new.env(parent = globalenv()),
        quiet = TRUE
    )

    out_file
}
