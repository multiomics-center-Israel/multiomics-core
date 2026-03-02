#' Module: Cross-omics pathway enrichment
#'
#' Combines pathway enrichment results from multiple omics layers
#' to identify consistently dysregulated pathways.
#'
#' @param enrichment_results Named list of enrichment results per omics
#' @param config Full config object
#' @param out_dir Output directory
#' @return List with: cross_omics_enrichment, plots
mod_multiomics_enrichment <- function(enrichment_results, config, out_dir) {

    message("\n=== Cross-Omics Pathway Enrichment ===\n")

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Check if enrichment is enabled
    if (!isTRUE(config$modes$multiomics$enrichment$run_enrichment)) {
        message("  Cross-omics enrichment disabled in config")
        return(NULL)
    }

    # Skip if insufficient enrichment results
    if (is.null(enrichment_results) || length(enrichment_results) < 2) {
        message("  Skipping cross-omics enrichment: need ≥2 omics with enrichment results")
        return(NULL)
    }

    cross_omics_enrich <- analyze_cross_omics_enrichment(
        enrichment_results = enrichment_results,
        config = config,
        out_dir = out_dir
    )

    if (!is.null(cross_omics_enrich)) {
        write_cross_omics_enrichment(cross_omics_enrich, out_dir)
        message("Cross-omics enrichment analysis complete")
    }

    cross_omics_enrich
}
