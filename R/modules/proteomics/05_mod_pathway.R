#' Proteomics pathway enrichment module
#'
#' Thin wrapper around run_proteomics_pathway() that checks config.
#'
#' @param de_res  DE results list
#' @param pre     Preprocessed proteomics object
#' @param config  Full config list
#' @param out_dir Output root directory
#' @return Pathway results list or NULL
mod_proteomics_pathway <- function(de_res, pre, config, out_dir) {
    cfg <- config$modes$proteomics
    if (!isTRUE(cfg$pathway$enabled)) {
        message("Pathway analysis disabled.")
        return(NULL)
    }
    run_proteomics_pathway(de_res, pre, config, out_dir)
}
