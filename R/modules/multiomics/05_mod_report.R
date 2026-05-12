#' Module: Multi-omics integration report
#'
#' Generates final HTML report summarizing all multi-omics integration results.
#'
#' @param run_dir Run output directory
#' @param config Full config object
#' @param config_file Path to config file
#' @return Path to rendered report (character)
mod_multiomics_report <- function(run_dir, config, config_file = NULL) {

    message("\n=== Generating Multi-Omics Report ===\n")

    report_path <- render_multiomics_report(
        run_dir = run_dir,
        config = config,
        config_file = config_file
    )

    message("Report generated: ", report_path)

    report_path
}
