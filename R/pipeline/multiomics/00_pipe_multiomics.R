#' Multi-Omics Integration Pipeline
#'
#' Orchestrates multi-omics integration analysis using {targets}.
#' Requires that individual omics pipelines (rna, proteomics, metabolomics)
#' have already run and produced preprocessed data.
#'
#' @return List of target objects
pipe_multiomics <- function() {

    list(
        # Output directory for multi-omics results
        tar_target(
            multiomics_out_dir,
            get_mode_out_dir(run_dir, "multiomics")
        ),

        # Harmonization: Load data, build MAE, gene-protein mapping
        tar_target(
            multiomics_harmonization,
            {
                # Check if at least 2 omics pipelines have run
                available_omics <- character(0)
                rna_data <- if (exists("rna_pre")) rna_pre else NULL
                prot_data <- if (exists("prot_pre")) prot_pre else NULL
                metab_data <- if (exists("metab_pre")) metab_pre else NULL

                if (!is.null(rna_data)) available_omics <- c(available_omics, "transcriptomics")
                if (!is.null(prot_data)) available_omics <- c(available_omics, "proteomics")
                if (!is.null(metab_data)) available_omics <- c(available_omics, "metabolomics")

                if (length(available_omics) < 2) {
                    message(
                        "Multi-omics integration requires ≥2 omics layers. ",
                        "Found: ", paste(available_omics, collapse = ", "), ". ",
                        "Skipping multi-omics integration."
                    )
                    return(NULL)
                }

                mod_multiomics_harmonization(
                    config = config,
                    rna_pre = rna_data,
                    prot_pre = prot_data,
                    metab_pre = metab_data,
                    out_dir = multiomics_out_dir
                )
            }
        ),

        # Collect DE results from individual omics pipelines
        tar_target(
            multiomics_de_results,
            {
                de_list <- list()

                if (exists("rna_de_res") && !is.null(rna_de_res)) {
                    de_list$transcriptomics <- rna_de_res
                }

                if (exists("prot_de_res") && !is.null(prot_de_res)) {
                    de_list$proteomics <- prot_de_res
                }

                if (exists("metab_de_res") && !is.null(metab_de_res)) {
                    de_list$metabolomics <- metab_de_res
                }

                if (length(de_list) == 0) NULL else de_list
            }
        ),

        # Collect enrichment results from individual omics pipelines
        tar_target(
            multiomics_enrichment_results,
            {
                enrich_list <- list()

                if (exists("rna_pathway_res") && !is.null(rna_pathway_res)) {
                    enrich_list$transcriptomics <- rna_pathway_res
                }

                if (exists("prot_pathway_res") && !is.null(prot_pathway_res)) {
                    enrich_list$proteomics <- prot_pathway_res
                }

                if (exists("metab_enrichment_res") && !is.null(metab_enrichment_res)) {
                    enrich_list$metabolomics <- metab_enrichment_res
                }

                if (length(enrich_list) == 0) NULL else enrich_list
            }
        ),

        # Integration: DIABLO, SNF, etc.
        tar_target(
            multiomics_integration,
            {
                if (is.null(multiomics_harmonization)) {
                    message("Skipping integration: harmonization failed or <2 omics available")
                    return(NULL)
                }

                mod_multiomics_integration(
                    harmonization_res = multiomics_harmonization,
                    de_results = multiomics_de_results,
                    config = config,
                    out_dir = file.path(multiomics_out_dir, "integration")
                )
            }
        ),

        # Concordance analysis
        tar_target(
            multiomics_concordance,
            {
                if (is.null(multiomics_harmonization) ||
                    is.null(multiomics_de_results) ||
                    length(multiomics_de_results) < 2) {
                    message("Skipping concordance: insufficient DE results")
                    return(NULL)
                }

                mod_multiomics_concordance(
                    de_results = multiomics_de_results,
                    harmonization_res = multiomics_harmonization,
                    config = config,
                    out_dir = file.path(multiomics_out_dir, "concordance")
                )
            }
        ),

        # Cross-omics enrichment
        tar_target(
            multiomics_cross_enrichment,
            {
                if (is.null(multiomics_enrichment_results) ||
                    length(multiomics_enrichment_results) < 2) {
                    message("Skipping cross-omics enrichment: insufficient enrichment results")
                    return(NULL)
                }

                mod_multiomics_enrichment(
                    enrichment_results = multiomics_enrichment_results,
                    config = config,
                    out_dir = file.path(multiomics_out_dir, "cross_enrichment")
                )
            }
        ),

        # RNA-protein correlation (extended analysis)
        tar_target(
            multiomics_rna_protein_corr,
            {
                if (is.null(multiomics_harmonization)) {
                    message("Skipping RNA-protein correlation: harmonization failed")
                    return(NULL)
                }

                tryCatch({
                    run_rna_protein_correlation(
                        mae = multiomics_harmonization$mae,
                        de_results = multiomics_de_results,
                        gene_protein_mapping = multiomics_harmonization$gene_protein_mapping,
                        config = config,
                        out_dir = file.path(multiomics_out_dir, "rna_protein_correlation")
                    )
                }, error = function(e) {
                    warning("RNA-protein correlation failed: ", e$message)
                    NULL
                })
            }
        ),

        # MultiGSEA plots
        tar_target(
            multiomics_multigsea,
            {
                if (is.null(multiomics_cross_enrichment)) {
                    message("Skipping MultiGSEA plots: no enrichment results")
                    return(NULL)
                }

                tryCatch({
                    run_multigsea_plots(
                        enrichment_results = multiomics_cross_enrichment,
                        config = config,
                        out_dir = file.path(multiomics_out_dir, "multigsea")
                    )
                }, error = function(e) {
                    warning("MultiGSEA plots failed: ", e$message)
                    NULL
                })
            }
        ),

        # Foundational cross-omics analysis (correlations, WGCNA, RNA-protein)
        tar_target(
            multiomics_foundational,
            {
                if (is.null(multiomics_harmonization)) {
                    message("Skipping foundational analysis: harmonization failed")
                    return(NULL)
                }

                mod_multiomics_foundational(
                    harmonization_res = multiomics_harmonization,
                    de_results = multiomics_de_results,
                    gene_protein_mapping = multiomics_harmonization$gene_protein_mapping,
                    config = config,
                    out_dir = file.path(multiomics_out_dir, "foundational")
                )
            }
        ),

        # Mechanistic inference (TF activity, regulatory networks, mediation)
        tar_target(
            multiomics_mechanistic,
            {
                if (is.null(multiomics_harmonization) || is.null(multiomics_foundational)) {
                    message("Skipping mechanistic analysis: dependencies not ready")
                    return(NULL)
                }

                mod_multiomics_mechanistic(
                    harmonization_res = multiomics_harmonization,
                    de_results = multiomics_de_results,
                    gene_protein_mapping = multiomics_harmonization$gene_protein_mapping,
                    foundational_results = multiomics_foundational$foundational_results,
                    config = config,
                    out_dir = file.path(multiomics_out_dir, "mechanistic")
                )
            }
        ),

        # Consensus and stability analysis
        tar_target(
            multiomics_consensus,
            {
                if (is.null(multiomics_harmonization) || is.null(multiomics_integration)) {
                    message("Skipping consensus analysis: dependencies not ready")
                    return(NULL)
                }

                mod_multiomics_consensus(
                    harmonization_res = multiomics_harmonization,
                    integration_res = multiomics_integration,
                    config = config,
                    out_dir = file.path(multiomics_out_dir, "consensus")
                )
            }
        ),

        # Commentary generation
        tar_target(
            multiomics_commentary,
            {
                if (is.null(multiomics_harmonization)) {
                    message("Skipping commentary: harmonization failed")
                    return(NULL)
                }

                mod_multiomics_commentary(
                    harmonization_res = multiomics_harmonization,
                    integration_res = multiomics_integration,
                    foundational_res = multiomics_foundational,
                    mechanistic_res = multiomics_mechanistic,
                    consensus_res = multiomics_consensus,
                    enrichment_res = multiomics_cross_enrichment,
                    config = config,
                    out_dir = file.path(multiomics_out_dir, "commentary")
                )
            }
        ),

        # Final report
        tar_target(
            multiomics_report,
            {
                # Wait for all analyses to complete
                force(multiomics_integration)
                force(multiomics_concordance)
                force(multiomics_cross_enrichment)
                force(multiomics_rna_protein_corr)
                force(multiomics_multigsea)
                force(multiomics_foundational)
                force(multiomics_mechanistic)
                force(multiomics_consensus)
                force(multiomics_commentary)

                if (is.null(multiomics_harmonization)) {
                    message("Skipping multi-omics report: no harmonization results")
                    return(NULL)
                }

                mod_multiomics_report(
                    run_dir = multiomics_out_dir,
                    config = config,
                    config_file = config_file
                )
            },
            format = "file"
        )
    )
}
