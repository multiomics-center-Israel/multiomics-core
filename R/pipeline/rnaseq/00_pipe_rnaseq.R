pipe_rnaseq <- function(skip_outputs = FALSE) {
  # ---- Core targets — always needed (multiomics depends on these) ----
  targets <- list(
    tar_target(rna_inputs, load_rna_inputs(config)),
    # Optional annotation inputs (NULL if not configured)
    tar_target(rna_annot, load_and_process_annotation(config)),
    tar_target(rna_trinotate_main, load_and_process_trinotate(config)),
    tar_target(rna_pre, preprocess_rna(rna_inputs, config, gene_lengths = NULL, verbose = TRUE)),
    tar_target(rna_out_dir, get_mode_out_dir(run_dir, "rna")),
    tar_target(rna_de_res, mod_rnaseq_de(rna_pre, rna_inputs, config, verbose = TRUE)),
    # Pathway / enrichment analysis (multiomics depends on this)
    tar_target(
      rna_pathway_res,
      mod_rnaseq_pathway(
        de_res  = rna_de_res,
        pre     = rna_pre,
        config  = config,
        out_dir = rna_out_dir
      )
    )
  )
  
  # ---- Single-omics outputs — skipped when multiomics pipeline is active ----
  if (!skip_outputs) {
    targets <- c(targets, list(
      tar_target(
        rna_qc_pre_obj,
        mod_rnaseq_qc_pre(
          pre     = rna_pre,
          config  = config,
          out_dir = rna_out_dir
        )
      ),
      # Clustering (run before outputs to provide excel_order)
      tar_target(
        rna_clustering_obj,
        mod_rnaseq_clustering(
          pre = rna_pre,
          de_res = rna_de_res,
          config = config,
          out_dir = rna_out_dir
        )
      ),
      # Legacy outputs (TSV files) - now receives clustering_res
      tar_target(
        rna_outputs_legacy,
        write_rnaseq_outputs_legacy(
          pre = rna_pre,
          de_res = rna_de_res,
          inputs = rna_inputs,
          config = config,
          out_dir = rna_out_dir,
          clustering_res = rna_clustering_obj
        ),
        format = "file"
      ),
      tar_target(
        rna_qc_post_obj,
        mod_rnaseq_qc_post(
          pre     = rna_pre,
          de_res  = rna_de_res,
          config  = config,
          out_dir = rna_out_dir
        )
      ),
      # Shiny payload export (canonical v2.0 with legacy aliases)
      tar_target(
        rna_shiny_payload,
        save_shiny_payload_rnaseq(
          pre = rna_pre,
          de_res = rna_de_res,
          inputs = rna_inputs,
          config = config,
          pca_res = rna_qc_pre_obj,
          clustering_res = rna_clustering_obj,
          annot = rna_annot,                    # Optional: gene annotation
          trinotate_main = rna_trinotate_main,  # Optional: Trinotate annotation
          out_file = file.path(rna_out_dir, "shiny_payload_rnaseq.rds")
        ),
        format = "file"
      ),
      # Executive summary
      tar_target(
        rna_exec_summary,
        mod_rnaseq_executive_summary(
          de_res      = rna_de_res,
          pathway_res = rna_pathway_res,
          qc_pre_obj  = rna_qc_pre_obj,
          pre         = rna_pre,
          config      = config,
          out_dir     = rna_out_dir
        ),
        format = "file"
      ),
      # AI commentary (runs after all plots are generated)
      tar_target(
        rna_commentary_file,
        {
          force(rna_qc_post_obj)
          force(rna_pathway_res)
          force(rna_exec_summary)
          mod_rnaseq_commentary(
            de_res     = rna_de_res,
            qc_pre_obj = rna_qc_pre_obj,
            config     = config,
            out_dir    = rna_out_dir
          )
        },
        format = "file"
      ),
      # Auto-report (final target — must wait for all analysis)
      tar_target(
        rna_report,
        {
          # Force dependencies so report renders AFTER all results are ready
          force(rna_pathway_res)
          force(rna_qc_post_obj)
          force(rna_outputs_legacy)
          force(rna_commentary_file)
          force(rna_exec_summary)
          render_rnaseq_report(
            run_dir     = rna_out_dir,
            config      = config,
            config_file = config_file
          )
        },
        format = "file"
      ),
      # Pipeline summary — dark-themed workflow overview HTML
      tar_target(
        rna_pipeline_summary,
        {
          force(rna_report)
          mod_rnaseq_pipeline_summary(
            config      = config,
            pre         = rna_pre,
            de_res      = rna_de_res,
            pathway_res = rna_pathway_res,
            run_dir     = run_dir
          )
        },
        format = "file"
      )
    ))
  }
  
  targets
}