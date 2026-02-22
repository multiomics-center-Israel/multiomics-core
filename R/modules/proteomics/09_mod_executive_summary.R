#' Proteomics executive summary module
#'
#' Thin wrapper around generate_proteomics_executive_summary().
#'
#' @param de_res       DE results list
#' @param pathway_res  Pathway results (optional)
#' @param qc_pre_obj   QC pre-DE results
#' @param pre          Preprocessed proteomics object
#' @param config       Full config list
#' @param out_dir      Output root directory
#' @param ppi_res      PPI results (optional)
#' @param adv_stats    Advanced stats results (optional)
#' @return Path to executive_summary.md or NULL
mod_proteomics_executive_summary <- function(de_res, pathway_res, qc_pre_obj,
                                               pre, config, out_dir,
                                               ppi_res = NULL, adv_stats = NULL) {
    cfg <- config$modes$proteomics
    # Executive summary is always generated if we have DE results
    if (is.null(de_res)) {
        message("No DE results. Skipping executive summary.")
        return(NULL)
    }

    generate_proteomics_executive_summary(
        de_res      = de_res,
        pathway_res = pathway_res,
        qc_pre_obj  = qc_pre_obj,
        pre         = pre,
        config      = config,
        out_dir     = out_dir,
        ppi_res     = ppi_res,
        adv_stats   = adv_stats
    )
}
