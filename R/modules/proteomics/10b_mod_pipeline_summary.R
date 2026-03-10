#' Module wrapper: Proteomics Pipeline Summary
#'
#' Generates a dark-themed HTML summary of the entire proteomics pipeline
#' workflow, including wet lab, upstream proteomics, and downstream analysis
#' steps. Uses Claude Code CLI for polished descriptions when commentary is
#' enabled.
#'
#' @param config      Full pipeline config
#' @param pre         Preprocessed data list
#' @param de_res      DE results list
#' @param pathway_res Pathway results list
#' @param run_dir     Run output directory
#' @return Path to pipeline_summary.html (character, format = "file")
mod_proteomics_pipeline_summary <- function(config, pre, de_res, pathway_res, run_dir) {
    message("[mod_pipeline_summary] Generating proteomics pipeline summary...")
    generate_proteomics_pipeline_summary(
        config      = config,
        pre         = pre,
        de_res      = de_res,
        pathway_res = pathway_res,
        run_dir     = run_dir
    )
}
