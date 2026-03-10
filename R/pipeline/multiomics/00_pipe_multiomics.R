#' Multi-Omics Integration Pipeline
#'
#' Orchestrates multi-omics integration analysis using {targets}.
#' Requires that individual omics pipelines (rna, proteomics, metabolomics)
#' have already run and produced preprocessed data.
#'
#' @param cfg_raw Raw config list (from yaml::read_yaml) used at plan-definition
#'   time to determine which omics modes are active. This allows {targets} to
#'   wire upstream dependencies correctly via direct symbol references.
#' @return List of target objects
pipe_multiomics <- function(cfg_raw) {

    # Determine which omics modes are active at plan-definition time.
    # Gateway targets below use these flags to include a direct symbol reference
    # (so {targets} detects the upstream dependency) or return NULL (no dep).
    has_rna   <- !is.null(cfg_raw$modes$rna)
    has_prot  <- !is.null(cfg_raw$modes$proteomics)
    has_metab <- !is.null(cfg_raw$modes$metabolomics)

    list(
        # Output directory for multi-omics results
        tar_target(
            multiomics_out_dir,
            get_mode_out_dir(run_dir, "multiomics")
        ),

        # ---------------------------------------------------------------
        # Gateway targets: wire upstream single-omics targets into the
        # multiomics DAG only when the corresponding mode is active.
        # When active, the bare symbol (e.g. rna_pre) is a direct reference
        # that {targets} detects for dependency tracking.
        # When inactive, the target resolves to NULL with no upstream dep.
        # ---------------------------------------------------------------

        # Preprocessing inputs
        if (has_rna)   tar_target(mo_rna_input,   rna_pre)   else tar_target(mo_rna_input,   NULL),
        if (has_prot)  tar_target(mo_prot_input,  prot_pre)  else tar_target(mo_prot_input,  NULL),
        if (has_metab) tar_target(mo_metab_input, metab_pre) else tar_target(mo_metab_input, NULL),

        # DE results
        if (has_rna)   tar_target(mo_rna_de,   rna_de_res)   else tar_target(mo_rna_de,   NULL),
        if (has_prot)  tar_target(mo_prot_de,  prot_de_res)  else tar_target(mo_prot_de,  NULL),
        if (has_metab) tar_target(mo_metab_de, metab_de_res) else tar_target(mo_metab_de, NULL),

        # Enrichment results
        if (has_rna)   tar_target(mo_rna_enrich,   rna_pathway_res)      else tar_target(mo_rna_enrich,   NULL),
        if (has_prot)  tar_target(mo_prot_enrich,  prot_pathway_res)     else tar_target(mo_prot_enrich,  NULL),
        if (has_metab) tar_target(mo_metab_enrich, metab_enrichment_res) else tar_target(mo_metab_enrich, NULL),

        # Harmonization: Load data, build MAE, gene-protein mapping
        tar_target(
            multiomics_harmonization,
            {
                input_mode <- config$modes$multiomics$input_mode %||% "pipeline"

                if (input_mode == "outputs") {
                    # Load from pre-computed shiny payload RDS files
                    mod_multiomics_harmonization(
                        config = config,
                        rna_pre = NULL,
                        prot_pre = NULL,
                        metab_pre = NULL,
                        out_dir = multiomics_out_dir
                    )
                } else {
                    rna_data   <- mo_rna_input
                    prot_data  <- mo_prot_input
                    metab_data <- mo_metab_input

                    available_omics <- character(0)
                    if (!is.null(rna_data))   available_omics <- c(available_omics, "transcriptomics")
                    if (!is.null(prot_data))  available_omics <- c(available_omics, "proteomics")
                    if (!is.null(metab_data)) available_omics <- c(available_omics, "metabolomics")

                    if (length(available_omics) < 2) {
                        message(
                            "Multi-omics integration requires \u22652 omics layers. ",
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
            }
        ),

        # Collect DE results from individual omics pipelines
        tar_target(
            multiomics_de_results,
            {
                de_list <- Filter(Negate(is.null), list(
                    transcriptomics = mo_rna_de,
                    proteomics      = mo_prot_de,
                    metabolomics    = mo_metab_de
                ))
                if (length(de_list) == 0) NULL else de_list
            }
        ),

        # Collect enrichment results from individual omics pipelines
        tar_target(
            multiomics_enrichment_results,
            {
                enrich_list <- Filter(Negate(is.null), list(
                    transcriptomics = mo_rna_enrich,
                    proteomics      = mo_prot_enrich,
                    metabolomics    = mo_metab_enrich
                ))
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
                mod_multiomics_enrichment(
                    enrichment_results = multiomics_enrichment_results,
                    de_results = multiomics_de_results,
                    harmonization_res = multiomics_harmonization,
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
