#' Module: Foundational Cross-Omics Analysis
#'
#' Runs foundational analysis BEFORE advanced integration methods (MOFA, DIABLO, SNF)
#' to establish baseline understanding of cross-omics relationships.
#'
#' Includes:
#'   - Cross-omics feature correlations (pairwise, partial)
#'   - Sample-level concordance across omics
#'   - Pathway overlap analysis
#'   - Cross-omics co-expression modules (WGCNA)
#'   - RNA-protein correlation (if both omics present)
#'
#' @param harmonization_res Output from mod_multiomics_harmonization()
#' @param de_results Named list of DE results per omics
#' @param gene_protein_mapping Gene-protein mapping table (optional)
#' @param config Full config object
#' @param out_dir Output directory
#' @return List with: foundational_results, rna_protein_correlation
mod_multiomics_foundational <- function(harmonization_res, de_results = NULL,
                                         gene_protein_mapping = NULL,
                                         config, out_dir) {

    message("\n=== Foundational Cross-Omics Analysis ===\n")

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    mae <- harmonization_res$mae
    cfg <- config$modes$multiomics$foundational %||% list()

    # Check if foundational analysis is enabled
    if (!(cfg$run_foundational %||% TRUE)) {
        message("Foundational analysis disabled in config. Skipping.")
        return(list(
            foundational_results = NULL,
            rna_protein_correlation = NULL
        ))
    }

    results <- list()

    # 1. Main foundational analysis
    message("\nRunning foundational cross-omics analysis...")
    foundational_dir <- file.path(out_dir, "foundational")
    dir.create(foundational_dir, showWarnings = FALSE)

    results$foundational <- tryCatch({
        run_foundational_analysis(
            mae = mae,
            de_results = de_results,
            gene_protein_mapping = gene_protein_mapping,
            config = config,
            out_dir = foundational_dir
        )
    }, error = function(e) {
        warning("Foundational analysis failed: ", e$message)
        NULL
    })

    # 2. RNA-Protein correlation (if both omics present)
    omics_names <- names(mae@ExperimentList)
    has_rna <- any(grepl("rna|transcriptomics", omics_names, ignore.case = TRUE))
    has_protein <- any(grepl("protein|proteomics", omics_names, ignore.case = TRUE))

    if (has_rna && has_protein) {
        message("\nRunning RNA-Protein correlation analysis...")
        rna_protein_dir <- file.path(out_dir, "rna_protein")
        dir.create(rna_protein_dir, showWarnings = FALSE)

        results$rna_protein <- tryCatch({
            run_rna_protein_correlation(
                mae = mae,
                de_results = de_results,
                gene_protein_mapping = gene_protein_mapping,
                config = config,
                out_dir = rna_protein_dir
            )
        }, error = function(e) {
            warning("RNA-Protein correlation failed: ", e$message)
            NULL
        })
    } else {
        message("\nSkipping RNA-Protein correlation (requires both transcriptomics and proteomics)")
        results$rna_protein <- NULL
    }

    message("\nFoundational analysis complete")

    # Summary
    summary_stats <- list(
        foundational_complete = !is.null(results$foundational),
        rna_protein_complete = !is.null(results$rna_protein)
    )

    if (!is.null(results$foundational$summary)) {
        summary_stats$n_correlation_pairs <- length(names(results$foundational$feature_correlations)) - 1
    }

    if (!is.null(results$rna_protein$rna_protein$summary)) {
        summary_stats$rna_protein_mean_cor <- results$rna_protein$rna_protein$summary$mean_cor
    }

    list(
        foundational_results = results$foundational,
        rna_protein_correlation = results$rna_protein,
        summary = summary_stats
    )
}
