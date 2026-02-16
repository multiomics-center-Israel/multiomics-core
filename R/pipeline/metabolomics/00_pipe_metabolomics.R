# R/pipeline/metabolomics/00_pipe_metabolomics.R
#
# Targets assembly for metabolomics Stage 1:
#   inspect input data, normalize, and produce QC diagnostics.


pipe_metabolomics <- function() {
    list(
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

        # ---- QC diagnostics ----
        tar_target(
            metab_qc_pre_obj,
            mod_metabolomics_qc_pre(
                pre     = metab_pre,
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
                pre = metab_pre,
                de_res = NULL,  # No DE analysis yet for metabolomics
                inputs = metab_inputs,
                config = config,
                pca_res = metab_qc_pre_obj,
                clustering_res = NULL,  # No clustering yet for metabolomics
                include_legacy = TRUE,
                out_file = file.path(metab_out_dir, "shiny_payload_metabolomics.rds")
            ),
            format = "file"
        )
    )
}


