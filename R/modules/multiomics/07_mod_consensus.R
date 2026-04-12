#' Module: Consensus and Stability Analysis
#'
#' Compares results across integration methods and assesses stability.
#' Runs AFTER integration methods (MOFA, DIABLO, SNF) complete.
#'
#' Includes:
#'   - Method comparison (sample clustering, feature importance)
#'   - Consensus patterns (features/clusters robust across methods)
#'   - Method-specific patterns
#'   - Meta-integration of latent factors
#'   - Bootstrap feature stability
#'   - K-fold cross-validation stability
#'   - Cluster stability assessment
#'
#' @param harmonization_res Output from mod_multiomics_harmonization()
#' @param integration_res Output from mod_multiomics_integration()
#' @param config Full config object
#' @param out_dir Output directory
#' @return List with: consensus_results, stability_results, summary
mod_multiomics_consensus <- function(harmonization_res, integration_res,
                                      config, out_dir) {

    message("\n=== Consensus and Stability Analysis ===\n")

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    mae <- harmonization_res$mae
    cfg <- config$modes$multiomics %||% list()

    results <- list()

    # 1. Integration consensus analysis
    if (cfg$consensus$run_consensus %||% TRUE) {
        message("\nRunning integration consensus analysis...")
        consensus_dir <- file.path(out_dir, "consensus")
        dir.create(consensus_dir, showWarnings = FALSE)

        # Build integration results list for consensus
        integration_results <- list(
            mofa = integration_res$mofa_results,
            diablo = integration_res$diablo_results,
            snf = integration_res$snf_results
        )

        results$consensus <- tryCatch({
            run_integration_consensus(
                integration_results = integration_results,
                mae = mae,
                config = config,
                out_dir = consensus_dir
            )
        }, error = function(e) {
            warning("Consensus analysis failed: ", e$message)
            NULL
        })
    } else {
        message("\nConsensus analysis disabled in config. Skipping.")
        results$consensus <- NULL
    }

    # 2. Stability analysis
    if (cfg$stability$run_stability %||% TRUE) {
        message("\nRunning stability analysis...")
        stability_dir <- file.path(out_dir, "stability")
        dir.create(stability_dir, showWarnings = FALSE)

        # Build integration results list
        integration_results <- list(
            mofa = integration_res$mofa_results,
            diablo = integration_res$diablo_results,
            snf = integration_res$snf_results
        )

        # Get feature data from integration module
        feature_data <- integration_res$feature_selection

        results$stability <- tryCatch({
            run_stability_analysis(
                mae = mae,
                feature_data = feature_data,
                integration_results = integration_results,
                config = config,
                out_dir = stability_dir
            )
        }, error = function(e) {
            warning("Stability analysis failed: ", e$message)
            NULL
        })
    } else {
        message("\nStability analysis disabled in config. Skipping.")
        results$stability <- NULL
    }

    message("\nConsensus and stability analysis complete")

    # Summary
    summary_stats <- list(
        consensus_complete = !is.null(results$consensus),
        stability_complete = !is.null(results$stability)
    )

    if (!is.null(results$consensus$summary)) {
        summary_stats$mean_method_ari <- results$consensus$summary$mean_method_ari
        summary_stats$n_robust_features <- results$consensus$summary$n_robust_features
    }

    if (!is.null(results$stability$summary)) {
        summary_stats$n_highly_stable_features <- results$stability$summary$n_highly_stable_features
        summary_stats$mean_cluster_stability <- results$stability$summary$mean_cluster_stability
    }

    list(
        consensus_results = results$consensus,
        stability_results = results$stability,
        summary = summary_stats
    )
}


#' Module: Commentary Generation
#'
#' Generates AI-powered or deterministic commentary for pipeline figures.
#' Runs after all analysis modules complete.
#'
#' @param harmonization_res Output from mod_multiomics_harmonization()
#' @param integration_res Output from mod_multiomics_integration()
#' @param foundational_res Output from mod_multiomics_foundational()
#' @param mechanistic_res Output from mod_multiomics_mechanistic()
#' @param consensus_res Output from mod_multiomics_consensus()
#' @param enrichment_res Output from enrichment analysis
#' @param config Full config object
#' @param out_dir Output directory
#' @return List with: figures_table, commentary_table
mod_multiomics_commentary <- function(harmonization_res,
                                       integration_res = NULL,
                                       foundational_res = NULL,
                                       mechanistic_res = NULL,
                                       consensus_res = NULL,
                                       enrichment_res = NULL,
                                       config,
                                       out_dir) {

    message("\n=== Figure Commentary Generation ===\n")

    commentary_cfg <- config$modes$multiomics$commentary %||% list()

    if (!(commentary_cfg$enabled %||% FALSE)) {
        message("Commentary generation disabled in config. Skipping.")
        return(list(
            figures_table = NULL,
            commentary_table = NULL
        ))
    }

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    mae <- harmonization_res$mae

    # Build integration results
    integration_results <- NULL
    if (!is.null(integration_res)) {
        integration_results <- list(
            mofa = integration_res$mofa_results,
            diablo = integration_res$diablo_results,
            snf = integration_res$snf_results
        )
    }

    # Build concordance results (from foundational)
    concordance_results <- NULL
    if (!is.null(foundational_res$rna_protein_correlation)) {
        concordance_results <- list(
            rna_protein = foundational_res$rna_protein_correlation
        )
    }

    # 1. Build figures table
    message("\nBuilding figures metadata table...")
    figures_tbl <- tryCatch({
        build_figures_table(
            mae = mae,
            foundational_results = foundational_res$foundational_results,
            mechanistic_results = mechanistic_res,
            integration_results = integration_results,
            concordance_results = concordance_results,
            enrichment_results = enrichment_res,
            consensus_results = consensus_res$consensus_results,
            stability_results = consensus_res$stability_results,
            config = config,
            out_dir = out_dir
        )
    }, error = function(e) {
        warning("Failed to build figures table: ", e$message)
        NULL
    })

    if (is.null(figures_tbl) || nrow(figures_tbl) == 0) {
        message("No figures found for commentary generation.")
        return(list(
            figures_table = figures_tbl,
            commentary_table = NULL
        ))
    }

    # 2. Generate commentary
    message("\nGenerating commentary for ", nrow(figures_tbl), " figures...")
    commentary_dir <- file.path(out_dir, "commentary")

    commentary_tbl <- tryCatch({
        generate_all_commentary(
            figures_tbl = figures_tbl,
            config = config,
            mae = mae,
            foundational_results = foundational_res$foundational_results,
            mechanistic_results = mechanistic_res,
            integration_results = integration_results,
            concordance_results = concordance_results,
            consensus_results = consensus_res$consensus_results,
            stability_results = consensus_res$stability_results,
            out_dir = commentary_dir
        )
    }, error = function(e) {
        warning("Commentary generation failed: ", e$message)
        NULL
    })

    message("\nCommentary generation complete")

    list(
        figures_table = figures_tbl,
        commentary_table = commentary_tbl
    )
}
