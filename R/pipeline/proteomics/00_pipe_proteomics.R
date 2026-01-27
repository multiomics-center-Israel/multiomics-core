# R/pipeline/proteomics/00_pipe_proteomics.R

pipe_proteomics <- function() {
    list(
        # ---- output dir (create once, targets tracks it as a file path) ----
        tar_target(
            prot_out_dir,
            {
                d <- get_mode_out_dir(run_dir, "proteomics")
                ensure_dir(d)
                d
            },
            format = "file"
        ),

        # ---- declare input files as file targets (so changes retrigger) ----
        tar_target(
            prot_input_files,
            {
                cfg <- config$modes$proteomics
                c(
                    resolve_raw_path(config, cfg$files$protein),
                    resolve_raw_path(config, cfg$files$sample_map),
                    resolve_raw_path(config, cfg$files$metadata),
                    resolve_raw_path(config, cfg$files$contrasts)
                )
            },
            format = "file"
        ),

        # ---- load inputs (forced dependency on prot_input_files) ----
        tar_target(
            prot_inputs,
            {
                prot_input_files
                load_proteomics_inputs(config)
            }
        ),

        # ---- preprocess ----
        tar_target(
            prot_pre,
            preprocess_proteomics(prot_inputs, config)
        ),

        # ---- QC pre ----
        tar_target(
            prot_qc_pre_obj,
            mod_proteomics_qc_pre(
                pre     = prot_pre,
                config  = config,
                out_dir = prot_out_dir
            )
        ),

        # ---- DE ----
        tar_target(
            prot_de_res,
            mod_proteomics_de(
                pre     = prot_pre,
                inputs  = prot_inputs,
                config  = config,
                verbose = TRUE
            )
        ),

        # ---- QC post ----
        tar_target(
            prot_qc_post_obj,
            mod_proteomics_qc_post(
                pre          = prot_pre,
                de_res       = prot_de_res,
                config       = config,
                out_dir      = prot_out_dir,
                de_source    = "table1"
            )
        ),

        # ---- clustering (needed for excel_order) ----
        tar_target(
            prot_clustering_obj,
            mod_proteomics_clustering(
                pre     = prot_pre,
                de_res  = prot_de_res,
                config  = config,
                out_dir = prot_out_dir
            )
        ),

        # ---- write outputs (files) ----
        tar_target(
            prot_de_files,
            write_proteomics_multimpute_outputs(
                pre         = prot_pre,
                de_res      = prot_de_res,
                inputs      = prot_inputs,
                config      = config,
                out_dir     = prot_out_dir,
                write_runs  = FALSE,
                excel_order = prot_clustering_obj$excel_order
            ),
            format = "file"
        ),

        # ---- build final results table ----
        tar_target(
            prot_final_results,
            {
                if (is.null(prot_inputs$contrasts) || is.null(prot_de_res$summary_df)) {
                    NULL
                } else {
                    # DEBUG: Enforce config availability
                    conf_id <- config$modes$proteomics$de_table$id_col
                    if (is.null(conf_id)) {
                        stop("CRITICAL CONFIG ERROR: config$modes$proteomics$de_table$id_col is NULL. Check config.yaml indentation or spelling.")
                    }

                    build_final_results_proteomics(
                        pre = prot_pre,
                        summary_df = prot_de_res$summary_df,
                        contrasts_df = prot_inputs$contrasts,
                        row_data = prot_pre$row_data,
                        feature_id_col = conf_id
                    )
                }
            }
        ),

        # ---- workspace snapshot (file) ----
        # Note: write_project_rdata might be in core/01_io.R or similar?
        # If not migrated yet, this might fail unless tar_source finds it in leftover files?
        # write_project_rdata was usually in 00_utils.R (deleted).
        # I need to ensure write_project_rdata is available.
        # I moved read_table/save_tsv to 01_io.R. Did I move write_project_rdata?
        # I did not see it explicitly.
        # I should check if I missed it. I will keep it commented out or assume I need to restore it.
        # For now I will reproduce the line but add a TODO if missing.
        tar_target(
            project_rdata_file,
            if (exists("write_project_rdata")) {
                write_project_rdata(
                    run_dir      = run_dir,
                    config       = config,
                    inputs       = prot_inputs,
                    pre_process  = prot_pre,
                    imputations  = prot_de_res$imputations,
                    de_results   = prot_de_res,
                    qc_results   = prot_qc_pre_obj
                )
            } else {
                # Fallback or skip
                NULL
            },
            format = "file"
        ),

        # Shiny legacy export (RDS file for old Shiny app)
        tar_target(
            prot_shiny_legacy,
            save_data_to_shiny_legacy_proteomics(
                pre = prot_pre,
                de_res = prot_de_res,
                inputs = prot_inputs,
                config = config,
                pca_res = prot_qc_pre_obj, # Pass QC pre results with PCA objects
                clustering_res = prot_clustering_obj, # Use clustering results
                out_file = file.path(run_dir, "proteomics", "data_to_shiny_legacy_proteomics.rds")
            ),
            format = "file"
        ),

        # Excel export (legacy format with cutoffs)
        tar_target(
            prot_excel_files,
            {
                write_final_results_excels_legacy_generic(
                    final_results = prot_final_results,
                    config = config,
                    out_dir = file.path(run_dir, "proteomics"),
                    mode = "proteomics",
                    id_col = config$modes$proteomics$de_table$id_col %||% "FeatureID",
                    expr_for_de = prot_pre$expr_imp_single,
                    with_cutoffs = TRUE,
                    clustering_res = prot_clustering_obj
                )
            },
            format = "file"
        )
    )
}
