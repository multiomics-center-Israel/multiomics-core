# R/pipeline/metabolomics/00_pipe_metabolomics.R
#
# {targets} assembly for metabolomics:
#
#   Stage 0 (new): met_* preprocessing DAG
#     met_raw → met_missingness_stats → met_missingness_plot
#             → met_filtered → met_imputed → met_log
#                                              ├─ met_norm_tss    → met_norm_tss_qc
#                                              ├─ met_norm_median → met_norm_median_qc
#                                              ├─ met_norm_pqn    → met_norm_pqn_qc
#                                              └─ met_norm_comparison
#                                           [chosen via config] → met_corrected
#
#   metab_pre adapter (bridges new → existing contract)
#
#   Stage 1: metab_qc_pre_obj (QC diagnostics)
#   Stage 2: metab_de_res, metab_feature_sel_res, metab_enrichment_res,
#            metab_standard_outputs, metab_shiny_payload, metab_report


pipe_metabolomics <- function() {
    list(

        # ------------------------------------------------------------------
        # Output directory (shared by all met_* and metab_* targets)
        # ------------------------------------------------------------------
        tar_target(
            metab_out_dir,
            {
                d <- get_mode_out_dir(run_dir, "metabolomics")
                ensure_dir(d)
                d
            },
            format = "file"
        ),

        # ------------------------------------------------------------------
        # Declare input files as a file target (invalidates on file change)
        # ------------------------------------------------------------------
        tar_target(
            metab_input_files,
            {
                cfg   <- config$modes$metabolomics
                paths <- resolve_raw_path(config, cfg$files$data)
                if (!is.null(cfg$files$metadata) && nzchar(cfg$files$metadata)) {
                    paths <- c(paths, resolve_raw_path(config, cfg$files$metadata))
                }
                paths
            },
            format = "file"
        ),

        # ------------------------------------------------------------------
        # metab_inputs: raw loaded inputs object (for Shiny payload)
        # ------------------------------------------------------------------
        tar_target(
            metab_inputs,
            {
                metab_input_files
                load_metabolomics_inputs(config)
            }
        ),

        # ==================================================================
        # met_* PREPROCESSING DAG (new)
        # ==================================================================

        # met_raw: parse + sample filter + align metadata
        tar_target(
            met_raw,
            {
                metab_input_files  # declare file dependency
                mod_met_raw(config)
            }
        ),

        # met_missingness_stats: MNAR/MAR classification on pre-filter matrix
        tar_target(
            met_missingness_stats,
            mod_met_missingness_stats(met_raw, config)
        ),

        # met_missingness_plot: diagnostic binary heatmap (file target)
        tar_target(
            met_missingness_plot,
            mod_met_missingness_plot(met_missingness_stats, met_raw,
                                     metab_out_dir, config),
            format = "file"
        ),

        # met_filtered: threshold-based feature + sample removal
        tar_target(
            met_filtered,
            mod_met_filtered(met_raw, config)
        ),

        # met_imputed: MNAR → min/2, MAR → KNN
        tar_target(
            met_imputed,
            mod_met_imputed(met_filtered, config)
        ),

        # met_log: log2 transform (pseudocount from config)
        tar_target(
            met_log,
            mod_met_log(met_imputed, config)
        ),

        # ---- Three parallel normalization targets ----

        tar_target(
            met_norm_tss,
            mod_met_normalize(met_log, method = "tss",
                              out_dir = metab_out_dir, config = config)
        ),

        tar_target(
            met_norm_median,
            mod_met_normalize(met_log, method = "median",
                              out_dir = metab_out_dir, config = config)
        ),

        tar_target(
            met_norm_pqn,
            mod_met_normalize(met_log, method = "pqn",
                              out_dir = metab_out_dir, config = config)
        ),

        # ---- QC file targets (boxplot + PCA already saved inside mod_met_normalize;
        #      these targets track the saved PNG paths for {targets} caching) ----

        tar_target(
            met_norm_tss_qc,
            {
                files <- c(met_norm_tss$boxplot_file, met_norm_tss$pca_file)
                files[!sapply(files, is.null)]
            },
            format = "file"
        ),

        tar_target(
            met_norm_median_qc,
            {
                files <- c(met_norm_median$boxplot_file, met_norm_median$pca_file)
                files[!sapply(files, is.null)]
            },
            format = "file"
        ),

        tar_target(
            met_norm_pqn_qc,
            {
                files <- c(met_norm_pqn$boxplot_file, met_norm_pqn$pca_file)
                files[!sapply(files, is.null)]
            },
            format = "file"
        ),

        # met_norm_comparison: TSV of RSD + PC1 variance for all three methods
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
        ),

        # met_corrected: select chosen_norm, apply optional LOESS drift correction
        tar_target(
            met_corrected,
            mod_met_corrected(
                norm_tss    = met_norm_tss,
                norm_median = met_norm_median,
                norm_pqn    = met_norm_pqn,
                logged      = met_log,
                meta        = met_log$meta,
                out_dir     = metab_out_dir,
                config      = config
            )
        ),

        # ==================================================================
        # metab_pre ADAPTER — bridges new met_* targets → existing contract
        #
        # Downstream metab_* targets (qc, de, enrichment, report, shiny)
        # depend on metab_pre and remain unchanged.
        # ==================================================================
        tar_target(
            metab_pre,
            {
                pre_cfg  <- config$modes$metabolomics$preprocessing %||% list()
                norm_cfg <- config$modes$metabolomics$normalization  %||% list()

                # Missingness stats for info (pre-filter, on raw matrix)
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
                    normalization_eval = NULL  # superseded by met_norm_comparison
                )
            }
        ),

        # ==================================================================
        # Stage 1: QC diagnostics
        # ==================================================================
        tar_target(
            metab_qc_pre_obj,
            mod_metabolomics_qc_pre(
                pre     = metab_pre,
                config  = config,
                out_dir = metab_out_dir
            )
        ),

        # ==================================================================
        # Stage 2: Differential expression, feature selection, enrichment
        # ==================================================================
        tar_target(
            metab_de_res,
            mod_metabolomics_de(
                pre     = metab_pre,
                config  = config,
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
        ),

        # ==================================================================
        # Standardized outputs, Shiny payload, HTML report
        # ==================================================================
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
        )
    )
}
