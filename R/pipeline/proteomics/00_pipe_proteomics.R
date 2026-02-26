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
                is_preprocessed <- identical(cfg$input$format, "preprocessed")
                paths <- c(
                    resolve_raw_path(config, cfg$files$metadata),
                    resolve_raw_path(config, cfg$files$contrasts)
                )
                if (is_preprocessed) {
                    pp <- cfg$files$preprocessed_protein %||% ""
                    if (nzchar(pp)) paths <- c(paths, resolve_raw_path(config, pp))
                } else {
                    paths <- c(paths, resolve_raw_path(config, cfg$files$protein))
                    if (nzchar(cfg$files$sample_map %||% ""))
                        paths <- c(paths, resolve_raw_path(config, cfg$files$sample_map))
                }
                paths
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
                de_source    = config$modes$proteomics$qc_post$de_source %||% "summary"
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
        ),

# ---- Pathway enrichment ----
        tar_target(
            prot_pathway_res,
            mod_proteomics_pathway(
                de_res  = prot_de_res,
                pre     = prot_pre,
                config  = config,
                out_dir = prot_out_dir
            )
        ),
         
        # ---- PPI network analysis ----
        tar_target(
            prot_ppi_res,
            mod_proteomics_ppi(
                de_res  = prot_de_res,
                pre     = prot_pre,
                config  = config,
                out_dir = prot_out_dir
            )
        ),

        # ---- Advanced statistics ----
        tar_target(
            prot_adv_stats,
            mod_proteomics_advanced_stats(
                de_res  = prot_de_res,
                pre     = prot_pre,
                inputs  = prot_inputs,
                config  = config,
                out_dir = prot_out_dir
            )
        ),

        # ---- Executive summary ----
        tar_target(
            prot_exec_summary,
            mod_proteomics_executive_summary(
                de_res      = prot_de_res,
                pathway_res = prot_pathway_res,
                qc_pre_obj  = prot_qc_pre_obj,
                pre         = prot_pre,
                config      = config,
                out_dir     = prot_out_dir,
                ppi_res     = prot_ppi_res,
                adv_stats   = prot_adv_stats
            )
        ),
         
        # ---- AI commentary ----
        tar_target(
            prot_commentary_file,
            {
                force(prot_qc_post_obj)
                force(prot_pathway_res)
                force(prot_ppi_res)
                force(prot_adv_stats)
                force(prot_exec_summary)
                mod_proteomics_commentary(prot_de_res, prot_qc_pre_obj, config, prot_out_dir)
            },
            format = "file"
        ),

        # User project summary (from user_notes + optional tech report)
        tar_target(
            prot_user_summary,
            {
                summary_html <- generate_project_summary(config)
                out_file <- file.path(get_mode_out_dir(run_dir, "proteomics"),
                                       "project_summary.html")
                if (!is.null(summary_html)) {
                    writeLines(summary_html, out_file)
                    out_file
                } else {
                    NA_character_
                }
            }
        ),

        # Proteomics HTML report
        tar_target(
            prot_report,
            {
                force(prot_pathway_res)
                force(prot_ppi_res)
                force(prot_adv_stats)
                force(prot_commentary_file)
                force(prot_exec_summary)
                force(prot_exports)
                force(prot_user_summary)
                render_proteomics_report(
                    run_dir     = prot_out_dir,
                    config      = config,
                    config_file = config_file
                )
            },
            format = "file"
        ),

#         Pipeline summary — dark-themed workflow overview HTML
        tar_target(
            prot_pipeline_summary,
            {
                force(prot_report)
                mod_proteomics_pipeline_summary(
                    config      = config,
                    pre         = prot_pre,
                    de_res      = prot_de_res,
                    pathway_res = prot_pathway_res,
                    run_dir     = run_dir
                )
            },
            format = "file"
        )
    )
}
