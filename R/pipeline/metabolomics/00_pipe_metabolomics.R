# R/pipeline/metabolomics/00_pipe_metabolomics.R
#
# Targets assembly for metabolomics:
#   Stage 1: inspect input data, normalize, and produce QC diagnostics.
#   Stage 2: differential expression, feature selection, enrichment.


pipe_metabolomics <- function(skip_outputs = FALSE) {
    # ---- Core targets (always needed — multiomics depends on these) ----
    targets <- list(
        # ---- output dir ----
        tar_target(
            metab_out_dir,
            {
                d <- get_mode_out_dir(run_dir, "metabolomics")
                ensure_dir(d)
                d
            },
            format = "file"
        ),

        # ---- declare input files as file targets ----
        tar_target(
            metab_input_files,
            {
                cfg <- config$modes$metabolomics
                paths <- resolve_raw_path(config, cfg$files$data)
                if (!is.null(cfg$files$metadata) && nzchar(cfg$files$metadata)) {
                    paths <- c(paths, resolve_raw_path(config, cfg$files$metadata))
                }
                paths
            },
            format = "file"
        ),

        # ---- load raw inputs ----
        tar_target(
            metab_inputs,
            {
                metab_input_files
                load_metabolomics_inputs(config)
            }
        ),

        # ---- preprocess (parse + normalize) ----
        tar_target(
            metab_pre,
            preprocess_metabolomics(metab_inputs, config)
        ),

        # ---- differential expression (Stage 2) ----
        tar_target(
            metab_de_res,
            mod_metabolomics_de(
                pre     = metab_pre,
                config  = config,
                out_dir = metab_out_dir
            )
        )
    )

    # ---- Single-omics outputs (skipped when multiomics pipeline is active) ----
    if (!skip_outputs) {
        targets <- c(targets, list(
            # ---- QC diagnostics (Stage 1) ----
            tar_target(
                metab_qc_pre_obj,
                mod_metabolomics_qc_pre(
                    pre     = metab_pre,
                    config  = config,
                    out_dir = metab_out_dir
                )
            ),

            # ---- feature selection: RF + PLS-DA (Stage 2, optional) ----
            tar_target(
                metab_feature_sel_res,
                mod_metabolomics_feature_selection(
                    pre     = metab_pre,
                    config  = config,
                    out_dir = metab_out_dir
                )
            ),

            # ---- pathway enrichment: QEA + ssGSEA + ORA + GSEA (Stage 2, optional) ----
            tar_target(
                metab_enrichment_res,
                mod_metabolomics_enrichment(
                    pre     = metab_pre,
                    de_res  = metab_de_res,
                    config  = config,
                    out_dir = metab_out_dir
                )
            ),

            # ---- write standardized outputs ----
            tar_target(
                metab_standard_outputs,
                write_metabolomics_outputs(
                    pre     = metab_pre,
                    config  = config,
                    out_dir = metab_out_dir
                ),
                format = "file"
            ),

            # ---- Shiny payload export (canonical v2.0) ----
            tar_target(
                metab_shiny_payload,
                save_shiny_payload_metabolomics(
                    pre            = metab_pre,
                    de_res         = metab_de_res,
                    inputs         = metab_inputs,
                    config         = config,
                    pca_res        = metab_qc_pre_obj,
                    clustering_res = NULL,
                    rf_res         = if (!is.null(metab_feature_sel_res)) metab_feature_sel_res$rf else NULL,
                    plsda_res      = if (!is.null(metab_feature_sel_res)) metab_feature_sel_res$plsda else NULL,
                    enrichment_res = metab_enrichment_res,
                    include_legacy = TRUE,
                    out_file       = file.path(metab_out_dir,
                                                "shiny_payload_metabolomics.rds")
                ),
                format = "file"
            ),

            # ---- HTML report ----
            tar_target(
                metab_report,
                mod_metabolomics_report(
                    pre             = metab_pre,
                    qc_res          = metab_qc_pre_obj,
                    de_res          = metab_de_res,
                    feature_sel_res = metab_feature_sel_res,
                    enrichment_res  = metab_enrichment_res,
                    config          = config,
                    out_dir         = metab_out_dir
                ),
                format = "file"
            ),

            # ---- pipeline summary HTML ----
            tar_target(
                metab_pipeline_summary,
                generate_metab_pipeline_summary(
                    config          = config,
                    pre             = metab_pre,
                    de_res          = metab_de_res,
                    feature_sel_res = metab_feature_sel_res,
                    enrichment_res  = metab_enrichment_res,
                    run_dir         = metab_out_dir
                ),
                format = "file"
            ),

            # ---- PowerPoint summary presentation ----
            tar_target(
                metab_pptx,
                mod_metabolomics_powerpoint(
                    pre             = metab_pre,
                    qc_res          = metab_qc_pre_obj,
                    de_res          = metab_de_res,
                    feature_sel_res = metab_feature_sel_res,
                    enrichment_res  = metab_enrichment_res,
                    config          = config,
                    out_dir         = metab_out_dir
                ),
                format = "file"
            )
        ))
    }

    targets
}
