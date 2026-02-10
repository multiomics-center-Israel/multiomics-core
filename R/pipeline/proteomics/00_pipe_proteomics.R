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

        # ---- unified exports (TSV, Excel, Shiny payload) ----
        tar_target(
            prot_exports,
            mod_proteomics_exports(
                pre            = prot_pre,
                de_res         = prot_de_res,
                inputs         = prot_inputs,
                config         = config,
                qc_pre_obj     = prot_qc_pre_obj,
                clustering_res = prot_clustering_obj,
                out_dir        = prot_out_dir
            )
        )
    )
}
