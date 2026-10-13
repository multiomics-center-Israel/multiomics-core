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

#' Proteomics enrichment on the Rubin-pooled ranking (module)
#'
#' Runs only when pathway analysis and \code{pathway.pooled_ranking.enabled}
#' are both on and the DE step that actually ran was limma (read from the DE
#' result, not the config: precomputed tables and the other methods carry no
#' per-imputation limma fits to pool).
#'
#' @param de_res  DE results list from mod_proteomics_de()
#' @param pre     Preprocessed proteomics object
#' @param config  Full config list
#' @param out_dir Output root directory
#' @return list(ranking, results, settings) or NULL when skipped
mod_proteomics_pathway_pooled <- function(de_res, pre, config, out_dir) {
    pw_cfg <- config$modes$proteomics$pathway %||% list()
    if (!isTRUE(pw_cfg$enabled) || !isTRUE((pw_cfg$pooled_ranking %||% list())$enabled)) {
        message("Pooled-ranking enrichment disabled.")
        return(NULL)
    }
    if (!identical(de_res$method, "limma")) {
        message("Pooled-ranking enrichment skipped: it needs de.method 'limma' ",
                "(this run used '", de_res$method %||% "unknown", "').")
        return(NULL)
    }
    run_proteomics_pooled_enrichment(de_res, pre, config, out_dir)
}
