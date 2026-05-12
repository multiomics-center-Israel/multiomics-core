# R/pipeline/metabolomics/00_pipe_metabolomics.R
#
# {targets} assembly for metabolomics.
#
# ── Conditional DAG ──────────────────────────────────────────────────────────
#
# The pipeline is controlled by:
#   chosen_norm  → preprocessing path (NULL = review QC mode, else analysis mode)
#   skip_outputs → if TRUE, only core targets needed by multiomics are returned
#
# ── Scale Contract ────────────────────────────────────────────────────────────
#
#   met_raw      (Linear)
#          └─ met_filtered (Linear) ─────────────────────────────┐
#               └─ met_log (Log2)                               │ (linear feed)
#                    ├─ met_norm_median  (Log2, log-shift)      │
#                    │                                          │
#                    │    ┌─ met_norm_tss  (Log2) ←─────────────┤
#                    │    ├─ met_norm_pqn  (Log2) ←─────────────┘
#                    │    │
#                    │    ├─ met_norm_comparison  (file: datasets/)
#                    │    │
#                    └────┴──── [chosen_norm] ──→ met_corrected
#                                                      └─ metab_pre (adapter)
pipe_metabolomics <- function(chosen_norm = NULL, skip_outputs = FALSE) {
  # -- Validate chosen_norm at plan-definition time --------------------------
  valid_norms <- c("tss", "median", "pqn", "eigenms", "eigenms_forced")
  if (!is.null(chosen_norm)) {
    chosen_norm <- tolower(chosen_norm)
    if (!chosen_norm %in% valid_norms) {
      stop(sprintf(
        "pipe_metabolomics: chosen_norm = '%s' is not valid. Must be NULL or one of: %s",
        chosen_norm, paste(valid_norms, collapse = ", ")
      ))
    }
  }
  review_mode <- is.null(chosen_norm)
  
  # ==================================================================
  # BASE TARGETS — always included (both modes)
  # ==================================================================
  base_targets <- list(
    tar_target(
      metab_out_dir,
      {
        d <- get_mode_out_dir(run_dir, "metabolomics")
        ensure_dir(d)
        d
      },
      format = "file"
    ),
    tar_target(
      metab_input_files,
      {
        cfg <- config$modes$metabolomics
        fmt <- cfg$input[["format"]] %||% "cd_raw"
        paths <- if (fmt == "multi_level") {
          dir_path <- resolve_raw_path(config, cfg$files[["data_dir"]])
          pattern  <- cfg$input[["level_pattern"]] %||% "\\.xlsx$"
          sort(list.files(dir_path, pattern = pattern,
                          full.names = TRUE, recursive = FALSE))
        } else {
          resolve_raw_path(config, cfg$files[["data"]])
        }
        meta_path <- cfg$files[["metadata"]]
        if (!is.null(meta_path) && nzchar(meta_path))
          paths <- c(paths, resolve_raw_path(config, meta_path))
        sm_path <- cfg$files[["sample_map"]]
        if (!is.null(sm_path) && nzchar(sm_path))
          paths <- c(paths, resolve_raw_path(config, sm_path))
        paths
      },
      format = "file"
    ),
    tar_target(
      metab_inputs,
      {
        metab_input_files
        load_metabolomics_inputs(config)
      }
    ),
    tar_target(
      met_raw,
      mod_met_raw(metab_inputs, config)
    ),
    tar_target(
      met_dedup_log,
      met_raw$duplicate_log
    ),
    tar_target(
      met_filtered,
      mod_met_filtered(met_raw, config)
    ),
    tar_target(
      met_missingness_stats,
      {
        mat        <- met_filtered$mat
        sample_col <- met_raw$sample_col
        feat_n_miss   <- rowSums(is.na(mat))
        feat_pct_miss <- if (ncol(mat) > 0) feat_n_miss / ncol(mat) else numeric(nrow(mat))
        stats_df <- data.frame(
          feature_id  = rownames(mat),
          n_missing   = feat_n_miss,
          pct_missing = feat_pct_miss,
          stringsAsFactors = FALSE
        )
        samp_n_miss   <- colSums(is.na(mat))
        samp_pct_miss <- if (nrow(mat) > 0) samp_n_miss / nrow(mat) else numeric(ncol(mat))
        samp_miss_df  <- data.frame(
          sample_id   = colnames(mat),
          n_missing   = samp_n_miss,
          pct_missing = samp_pct_miss,
          stringsAsFactors = FALSE
        )
        colnames(samp_miss_df)[1] <- sample_col
        list(
          stats_df     = stats_df,
          samp_miss_df = samp_miss_df
        )
      }
    ),
    tar_target(
      met_log,
      mod_met_log(met_filtered, config)
    ),
    tar_target(
      met_norm_tss,
      mod_met_normalize_linear(met_filtered, method = "tss", config = config)
    ),
    tar_target(
      met_norm_median,
      mod_met_normalize_log(met_log, config = config)
    ),
    tar_target(
      met_norm_pqn,
      mod_met_normalize_linear(met_filtered, method = "pqn", config = config)
    ),
    tar_target(
      met_norm_eigenms,
      mod_met_normalize_eigenms(met_filtered, config = config)
    ),
    tar_target(
      met_norm_eigenms_forced,
      mod_met_normalize_eigenms_forced(met_filtered, config = config)
    ),
    tar_target(
      met_norm_comparison,
      mod_met_norm_comparison(
        norm_tss    = met_norm_tss,
        norm_median = met_norm_median,
        norm_pqn    = met_norm_pqn,
        logged      = met_log,
        out_dir     = metab_out_dir,
        config      = config
      ),
      format = "file"
    )
  )
  
  # ==================================================================
  # REVIEW MODE (chosen_norm = NULL): QC layer + summary report
  # ==================================================================
  if (review_mode) {
    qc_targets <- list(
      tar_target(
        met_qc_cfg,
        config$modes$metabolomics$qc %||% list()
      ),
      tar_target(
        met_log_qc,
        {
          .qc_cfg <- met_qc_cfg
          mod_met_qc_suite(met_log, stage = "log",
                           out_dir = metab_out_dir, config = config)
        },
        format = "file"
      ),
      tar_target(
        met_norm_tss_qc,
        {
          .qc_cfg <- met_qc_cfg
          mod_met_qc_suite(met_norm_tss, stage = "norm_tss",
                           out_dir = metab_out_dir, config = config)
        },
        format = "file"
      ),
      tar_target(
        met_norm_median_qc,
        {
          .qc_cfg <- met_qc_cfg
          mod_met_qc_suite(met_norm_median, stage = "norm_median",
                           out_dir = metab_out_dir, config = config)
        },
        format = "file"
      ),
      tar_target(
        met_norm_pqn_qc,
        {
          .qc_cfg <- met_qc_cfg
          mod_met_qc_suite(met_norm_pqn, stage = "norm_pqn",
                           out_dir = metab_out_dir, config = config)
        },
        format = "file"
      ),
      tar_target(
        met_norm_eigenms_qc,
        {
          .qc_cfg <- met_qc_cfg
          mod_met_qc_suite(met_norm_eigenms, stage = "norm_eigenms",
                           out_dir = metab_out_dir, config = config)
        },
        format = "file"
      ),
      tar_target(
        met_norm_eigenms_forced_qc,
        {
          .qc_cfg <- met_qc_cfg
          mod_met_qc_suite(met_norm_eigenms_forced, stage = "norm_eigenms_forced",
                           out_dir = metab_out_dir, config = config)
        },
        format = "file"
      ),
      tar_target(
        met_qc_comparison,
        mod_met_qc_comparison_table(
          log_qc_files              = met_log_qc,
          tss_qc_files              = met_norm_tss_qc,
          median_qc_files           = met_norm_median_qc,
          pqn_qc_files              = met_norm_pqn_qc,
          eigenms_qc_files          = met_norm_eigenms_qc,
          eigenms_forced_qc_files   = met_norm_eigenms_forced_qc,
          out_dir                   = metab_out_dir,
          config                    = config
        ),
        format = "file"
      ),
      tar_target(
        met_qc_summary_report,
        mod_met_qc_summary_report(
          qc_comparison_file = met_qc_comparison,
          qc_suite_files     = c(met_log_qc,
                                 met_norm_tss_qc,
                                 met_norm_median_qc,
                                 met_norm_pqn_qc,
                                 met_norm_eigenms_qc,
                                 met_norm_eigenms_forced_qc),
          config  = config,
          out_dir = metab_out_dir
        ),
        format = "file"
      )
    )
    return(c(base_targets, qc_targets))
  }
  
  # ==================================================================
  # ANALYSIS MODE — Core targets (always needed; multiomics depends on these)
  # ==================================================================
  analysis_core <- list(
    # met_corrected: select chosen_norm, apply optional LOESS drift correction
    tar_target(
      met_corrected,
      mod_met_corrected(
        norm_tss            = met_norm_tss,
        norm_median         = met_norm_median,
        norm_pqn            = met_norm_pqn,
        logged              = met_log,
        meta                = met_log$meta,
        out_dir             = metab_out_dir,
        config              = config,
        norm_eigenms        = met_norm_eigenms,
        norm_eigenms_forced = met_norm_eigenms_forced
      )
    ),
    # metab_pre ADAPTER — bridges met_* targets → existing contract
    tar_target(
      metab_pre,
      {
        pre_cfg  <- config$modes$metabolomics$preprocessing %||% list()
        norm_cfg <- config$modes$metabolomics$normalization  %||% list()
        miss_stats <- met_missingness_stats$stats_df
        miss_samp  <- met_missingness_stats$samp_miss_df
        info <- list(
          mode            = "metabolomics",
          n_features_raw  = nrow(met_raw$expr_raw),
          n_features_filt = nrow(met_filtered$mat),
          n_samples       = ncol(met_corrected$mat),
          missingness     = list(
            per_sample  = stats::setNames(
              miss_samp$pct_missing,
              miss_samp[[met_raw$sample_col]]
            ),
            per_feature = stats::setNames(
              miss_stats$pct_missing,
              miss_stats$feature_id
            )
          ),
          normalization   = met_corrected$info$normalization,
          norm_comparison_file = met_norm_comparison,
          drift           = met_corrected$info[c("drift_applied", "drift_info")]
        )
        list(
          expr_raw         = met_raw$expr_raw,
          expr_filt        = met_filtered$mat,
          expr_log         = met_log$mat,
          expr_work        = met_corrected$mat,
          meta             = met_corrected$meta,
          row_data         = met_corrected$row_data,
          info             = info,
          normalization_eval = NULL
        )
      }
    ),
    tar_target(
      metab_de_res,
      mod_metabolomics_de(
        pre     = metab_pre,
        config  = config,
        inputs  = metab_inputs,
        out_dir = metab_out_dir
      )
    ),
    tar_target(
      metab_feature_sel_res,
      mod_metabolomics_feature_selection(
        pre     = metab_pre,
        config  = config,
        out_dir = metab_out_dir
      )
    ),
    tar_target(
      metab_enrichment_res,
      mod_metabolomics_enrichment(
        pre     = metab_pre,
        de_res  = metab_de_res,
        config  = config,
        out_dir = metab_out_dir
      )
    )
  )
  
  # ==================================================================
  # ANALYSIS MODE — Output targets (skipped when multiomics is active)
  # ==================================================================
  if (skip_outputs) {
    return(c(base_targets, analysis_core))
  }
  
  analysis_outputs <- list(
    tar_target(
      metab_qc_pre_obj,
      mod_metabolomics_qc_pre(
        pre     = metab_pre,
        config  = config,
        out_dir = metab_out_dir
      )
    ),
    tar_target(
      metab_clustering_obj,
      mod_metabolomics_clustering(
        pre     = metab_pre,
        de_res  = metab_de_res,
        config  = config,
        out_dir = metab_out_dir
      )
    ),
    tar_target(
      metab_network,
      mod_metabolomics_network(
        de_res  = metab_de_res,
        pre     = metab_pre,
        config  = config,
        out_dir = metab_out_dir
      )
    ),
    tar_target(
      metab_standard_outputs,
      write_metabolomics_outputs(
        pre     = metab_pre,
        config  = config,
        out_dir = metab_out_dir
      ),
      format = "file"
    ),
    tar_target(
      metab_final_results,
      write_metabolomics_final_results(
        pre            = metab_pre,
        de_res         = metab_de_res,
        config         = config,
        out_dir        = metab_out_dir,
        clustering_res = metab_clustering_obj
      ),
      format = "file"
    ),
    tar_target(
      metab_shiny_payload,
      save_shiny_payload_metabolomics(
        pre            = metab_pre,
        de_res         = metab_de_res,
        inputs         = metab_inputs,
        config         = config,
        pca_res        = metab_qc_pre_obj,
        clustering_res = metab_clustering_obj,
        rf_res         = if (!is.null(metab_feature_sel_res)) metab_feature_sel_res$rf else NULL,
        plsda_res      = if (!is.null(metab_feature_sel_res)) metab_feature_sel_res$plsda else NULL,
        enrichment_res = metab_enrichment_res,
        include_legacy = TRUE,
        out_file       = file.path(metab_out_dir,
                                   "shiny_payload_metabolomics.rds")
      ),
      format = "file"
    ),
    tar_target(
      metab_commentary,
      mod_metabolomics_commentary(
        de_res     = metab_de_res,
        qc_pre_obj = metab_qc_pre_obj,
        config     = config,
        out_dir    = metab_out_dir
      ),
      format = "file"
    ),
    tar_target(
      metab_report,
      mod_metabolomics_report(
        pre                = metab_pre,
        qc_res             = metab_qc_pre_obj,
        de_res             = metab_de_res,
        clustering_res     = metab_clustering_obj,
        feature_sel_res    = metab_feature_sel_res,
        enrichment_res     = metab_enrichment_res,
        network_res        = metab_network,
        config             = config,
        out_dir            = metab_out_dir,
        qc_comparison_file = NULL,
        qc_suite_files     = NULL,
        commentary_file    = metab_commentary
      ),
      format = "file"
    ),
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
  )
  
  c(base_targets, analysis_core, analysis_outputs)
}