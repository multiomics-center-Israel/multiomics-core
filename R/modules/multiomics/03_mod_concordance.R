#' Module: Multi-omics concordance analysis
#'
#' Analyzes agreement/disagreement between omics layers based on
#' differential expression patterns.
#'
#' @param de_results Named list of DE results per omics
#' @param harmonization_res Output from mod_multiomics_harmonization()
#' @param config Full config object
#' @param out_dir Output directory
#' @return List with: concordance_results, plots
mod_multiomics_concordance <- function(de_results, harmonization_res, config, out_dir) {

    message("\n=== Multi-Omics Concordance Analysis ===\n")

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Skip if no DE results provided
    if (is.null(de_results) || length(de_results) < 2) {
        message("  Skipping concordance: need \u22652 omics with DE results")
        return(NULL)
    }

    concordance_res <- analyze_multiomics_concordance(
        de_results = de_results,
        gene_protein_mapping = harmonization_res$gene_protein_mapping,
        mae = harmonization_res$mae,
        config = config,
        out_dir = out_dir
    )

    # Write concordance tables
    if (!is.null(concordance_res) && !is.null(concordance_res$concordance)) {
        for (om_pair in names(concordance_res$concordance)) {
            conc <- concordance_res$concordance[[om_pair]]

            tbl <- conc$concordance_table %||% conc$merged
            if (!is.null(tbl) && nrow(tbl) > 0) {
                write.csv(
                    tbl,
                    file.path(out_dir, paste0("concordance_", om_pair, ".csv")),
                    row.names = FALSE
                )
            }
        }

        message("Concordance analysis complete")
    }

    concordance_res
}
