#' Module wrapper: RNA-seq Executive Summary
#'
#' Generates an executive summary markdown file with key findings.
#'
#' @param de_res DE results list
#' @param pathway_res Pathway/enrichment results
#' @param qc_pre_obj QC pre-DE results
#' @param pre Preprocessed data list
#' @param config Full config
#' @param out_dir Output directory
#' @return Path to the generated executive_summary.md file
mod_rnaseq_executive_summary <- function(de_res, pathway_res, qc_pre_obj, pre, config, out_dir) {
    message("[mod_executive_summary] Generating executive summary...")
    generate_rnaseq_executive_summary(
        de_res      = de_res,
        pathway_res = pathway_res,
        qc_pre_obj  = qc_pre_obj,
        pre         = pre,
        config      = config,
        out_dir     = out_dir
    )
}
