#' Module: Mechanistic Inference Analysis
#'
#' Runs deep mechanistic analysis on multi-omics data to infer regulatory
#' relationships and causal mechanisms.
#'
#' Includes:
#'   - RNA-Protein regulatory analysis (translation efficiency, protein stability)
#'   - Transcription factor activity inference
#'   - Gene regulatory network inference
#'   - Causal mediation analysis (RNA -> Protein -> Phenotype)
#'
#' @param harmonization_res Output from mod_multiomics_harmonization()
#' @param de_results Named list of DE results per omics
#' @param gene_protein_mapping Gene-protein mapping table (optional)
#' @param foundational_results Results from mod_multiomics_foundational() (optional)
#' @param config Full config object
#' @param out_dir Output directory
#' @return List with: mechanistic_results, summary
mod_multiomics_mechanistic <- function(harmonization_res, de_results = NULL,
                                        gene_protein_mapping = NULL,
                                        foundational_results = NULL,
                                        config, out_dir) {

    message("\n=== Mechanistic Inference Analysis ===\n")

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    mae <- harmonization_res$mae
    cfg <- config$modes$multiomics$mechanistic %||% list()

    # Check if mechanistic analysis is enabled
    if (!(cfg$run_mechanistic %||% TRUE)) {
        message("Mechanistic analysis disabled in config. Skipping.")
        return(list(
            mechanistic_results = NULL,
            summary = list(enabled = FALSE)
        ))
    }

    # Run mechanistic analysis
    mechanistic_dir <- file.path(out_dir, "mechanistic")
    dir.create(mechanistic_dir, showWarnings = FALSE)

    mechanistic_results <- tryCatch({
        run_mechanistic_analysis(
            mae = mae,
            de_results = de_results,
            gene_protein_mapping = gene_protein_mapping,
            foundational_results = foundational_results,
            config = config,
            out_dir = mechanistic_dir
        )
    }, error = function(e) {
        warning("Mechanistic analysis failed: ", e$message)
        NULL
    })

    message("\nMechanistic analysis complete")

    # Build summary
    summary_stats <- list(
        enabled = TRUE,
        complete = !is.null(mechanistic_results)
    )

    if (!is.null(mechanistic_results)) {
        # RNA-Protein regulation summary
        if (!is.null(mechanistic_results$rna_protein)) {
            rp <- mechanistic_results$rna_protein
            summary_stats$rna_protein_mean_cor <- rp$correlations$summary$mean_cor
            summary_stats$n_high_te_genes <- rp$translation_efficiency$summary$n_high_te
            summary_stats$n_ptr_genes <- rp$post_transcriptional$n_ptr
        }

        # TF activity summary
        if (!is.null(mechanistic_results$tf_activity)) {
            tf <- mechanistic_results$tf_activity
            summary_stats$n_tfs_analyzed <- nrow(tf$activity_matrix)
            if (!is.null(tf$associations)) {
                summary_stats$n_sig_tfs <- sum(tf$associations$padj < 0.05, na.rm = TRUE)
            }
        }

        # Regulatory network summary
        if (!is.null(mechanistic_results$regulatory_network)) {
            rn <- mechanistic_results$regulatory_network
            summary_stats$network_method <- rn$method
            summary_stats$n_network_edges <- nrow(rn$edges)
            summary_stats$n_hub_regulators <- sum(rn$hubs$is_hub, na.rm = TRUE)
        }

        # Mediation summary
        if (!is.null(mechanistic_results$mediation)) {
            med <- mechanistic_results$mediation
            summary_stats$n_genes_tested_mediation <- nrow(med$mediation_df)
            summary_stats$n_genes_mediated <- med$n_mediated
        }

        # Multi-layer mediation summary
        if (!is.null(mechanistic_results$multilayer_mediation)) {
            ml <- mechanistic_results$multilayer_mediation
            summary_stats$n_multilayer_paths_tested <- nrow(ml$multilayer_df)
            summary_stats$n_multilayer_significant <- ml$n_significant
        }

        # COSMOS summary
        if (!is.null(mechanistic_results$cosmos)) {
            cosmos <- mechanistic_results$cosmos
            summary_stats$cosmos_n_nodes <- cosmos$n_nodes
            summary_stats$cosmos_n_edges <- cosmos$n_edges
        }
    }

    list(
        mechanistic_results = mechanistic_results,
        summary = summary_stats
    )
}
