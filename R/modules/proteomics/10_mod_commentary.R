#' Proteomics AI commentary module
#'
#' Thin wrapper around run_proteomics_commentary().
#'
#' @param de_res     DE results list
#' @param qc_pre_obj QC pre-DE results
#' @param config     Full config list
#' @param out_dir    Output root directory
#' @return Path to commentary_all.json or NULL
mod_proteomics_commentary <- function(de_res, qc_pre_obj, config, out_dir) {
    cfg <- config$modes$proteomics
    if (!isTRUE(cfg$commentary$enabled)) {
        message("Commentary disabled.")
        return(NULL)
    }

    run_proteomics_commentary(de_res, qc_pre_obj, config, out_dir)
}
