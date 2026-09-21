pipe_rnaseq <- function(skip_outputs = FALSE) {
  # ---- Core targets — always needed (multiomics depends on these) ----
  targets <- list(
    # ---- declare input files as file targets (so changes retrigger) ----
    # Mirrors prot_input_files / metab_input_files. Without this the RNA counts,
    # metadata, contrasts and annotation were untracked, so editing any of them
    # left the pipeline believing it was up to date and the run silently kept
    # the previous results.
    #
    # Deliberately not a second opinion about which files are required: this
    # collects the configured paths and load_omics_inputs() remains the one
    # place that decides what RNA actually needs, which differs by input route
    # (counts / txi / preprocessed_counts).
    tar_target(
      rna_input_files,
      {
        files <- config$modes$rna$files
        paths <- character(0)
        for (nm in names(files)) {
          rel <- files[[nm]]
          # Scalar character entries, as load_omics_inputs() takes them: skips
          # flags (is_logtransformed) and multi-value keys, which are picked up
          # separately below where they are read separately too.
          if (is.null(rel) || !is.character(rel) || length(rel) != 1 ||
              !nzchar(rel)) {
            next
          }
          abs <- resolve_raw_path(config, rel)
          if (dir.exists(abs)) next   # directory entries (e.g. data_dir)
          paths <- c(paths, abs)
        }

        # de_table does not come through load_omics_inputs() at all -- that
        # loader skips multi-value keys -- and load_precomputed_rna_de() reads
        # it directly and supports a list of tables, one per contrast. Left out,
        # editing one of a set of pre-computed DE tables would leave this branch
        # looking up to date, which is the staleness this target exists to stop.
        de_tables <- files$de_table
        if (is.list(de_tables)) de_tables <- unlist(de_tables)
        for (rel in de_tables) {
          if (!is.character(rel) || !nzchar(rel)) next
          paths <- c(paths, resolve_raw_path(config, rel))
        }

        unique(paths)
      },
      format = "file"
    ),

    # ---- load inputs (forced dependency on rna_input_files) ----
    tar_target(
      rna_inputs,
      {
        rna_input_files
        load_rna_inputs(config)
      }
    ),
    # Optional annotation inputs (NULL if not configured).
    # These read config$modes$rna$files directly rather than going through
    # rna_inputs, so they need the file dependency of their own -- otherwise
    # editing the annotation or trinotate file rebuilds rna_inputs while these
    # two keep serving the cached objects built from the old file.
    tar_target(
      rna_annot,
      {
        rna_input_files
        load_and_process_annotation(config)
      }
    ),
    tar_target(
      rna_trinotate_main,
      {
        rna_input_files
        load_and_process_trinotate(config)
      }
    ),
    tar_target(rna_pre, preprocess_rna(rna_inputs, config, gene_lengths = NULL, verbose = TRUE)),
    tar_target(rna_out_dir, get_mode_out_dir(run_dir, "rna")),
    tar_target(rna_de_res, mod_rnaseq_de(rna_pre, rna_inputs, config, verbose = TRUE))
  )

  # ---- Pathway / enrichment (multiomics depends on this) ----
  # One stable target; its clustering input depends on the run mode:
  #   - multiomics (skip_outputs = TRUE): clustering is NOT produced, so pass
  #     clustering_res = NULL -> GSEA only, cluster-based ORA skipped with a warning.
  #   - single-omics (skip_outputs = FALSE): clustering_res = rna_clustering_obj
  #     -> GSEA + cluster-based ORA (target defined in the outputs block below).
  if (skip_outputs) {
    targets <- c(targets, list(
      tar_target(
        rna_pathway_res,
        mod_rnaseq_pathway(
          de_res         = rna_de_res,
          pre            = rna_pre,
          config         = config,
          out_dir        = rna_out_dir,
          clustering_res = NULL
        )
      )
    ))
    return(targets)
  }

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
      # Pathway / enrichment with clustering available (GSEA + cluster-based ORA)
      tar_target(
        rna_pathway_res,
        mod_rnaseq_pathway(
          de_res         = rna_de_res,
          pre            = rna_pre,
          config         = config,
          out_dir        = rna_out_dir,
          clustering_res = rna_clustering_obj
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
          xlsx_files = rna_outputs_legacy,      # paths of Final_results_{ALL,DE}_P_*.xlsx
          pathway_res = rna_pathway_res,        # Stage 3C: enrichment -> payload$enrichment
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