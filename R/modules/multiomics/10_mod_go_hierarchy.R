#' Module: GO hierarchy of the layers' GSEA results (compute)
#'
#' Builds the GO sub-graphs above each gene-based layer's significant GO fgsea
#' terms. Writes nothing; write_go_hierarchy_outputs() does that from its own
#' target. A failure (e.g. GO.db not installed) is reported and returns a
#' result with no units, so the writer still clears the previous run's files.
#'
#' @param enrichment_results Named list of layer pathway results (the
#'   multiomics_enrichment_results target).
#' @param config Full config object.
#' @return NULL when modes.multiomics.go_hierarchy.enabled is false, else the
#'   output of build_go_hierarchy_dags() (with an \code{error} element and no
#'   units when the build failed).
mod_multiomics_go_hierarchy <- function(enrichment_results, config) {
    gh <- go_hierarchy_config(config)
    if (!gh$enabled) {
        message("GO hierarchy figures disabled (modes.multiomics.go_hierarchy.enabled).")
        return(NULL)
    }
    message("\n=== GO hierarchy of GSEA results ===\n")
    tryCatch(
        build_go_hierarchy_dags(enrichment_results, config),
        error = function(e) {
            message("  GO hierarchy failed: ", conditionMessage(e))
            list(settings = gh, units = list(), error = conditionMessage(e))
        }
    )
}
