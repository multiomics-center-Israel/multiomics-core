#' Proteomics Domain Analysis Module
#'
#' Target-ready wrapper for protein domain enrichment and
#' functional feature annotation.

mod_proteomics_domain_analysis <- function(de_res, config, out_dir) {
    if (!isTRUE(config$modes$proteomics$domain_analysis$enabled %||% TRUE)) {
        message("Domain analysis disabled in config. Skipping.")
        return(NULL)
    }

    tryCatch(
        run_domain_analysis(de_res, config, out_dir),
        error = function(e) {
            message("Domain analysis failed: ", e$message)
            NULL
        }
    )
}
