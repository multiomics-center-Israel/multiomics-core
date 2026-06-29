# R/pipeline/lipidomics/00_pipe_lipidomics.R
#
# Targets assembly for lipidomics:
#   Stage 1: load inputs, normalize, QC diagnostics
#   Stage 2: DE, feature selection, lipid class analysis
#   Stage 3: biomarker discovery, pathway analysis, enhanced QC
#   Stage 4: HTML report, pipeline summary, PowerPoint


pipe_lipidomics <- function() {
    list(
        # ---- output dir ----
        tar_target(
            lipid_out_dir,
            {
                d <- get_mode_out_dir(run_dir, "lipidomics")
                ensure_dir(d)
                d
            },
            format = "file"
        ),

        # ---- declare input files as file targets ----
        tar_target(
            lipid_input_files,
            {
                cfg <- config$modes$lipidomics
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
            lipid_inputs,
            {
                lipid_input_files
                load_lipidomics_inputs(config)
            }
        ),

        # ---- preprocess (parse + normalize + annotate) ----
        tar_target(
            lipid_pre,
            {
                pre <- preprocess_lipidomics(lipid_inputs, config)
                # Enrich row_data with full lipid names and categories
                if (!is.null(pre$row_data)) {
                    pre$row_data <- annotate_lipid_row_data(pre$row_data)
                }
                pre
            }
        ),

        # ---- standardized outputs (expr matrices, metadata, annotations) ----
        tar_target(
            lipid_standard_outputs,
            write_lipidomics_outputs(
                pre     = lipid_pre,
                config  = config,
                out_dir = lipid_out_dir
            ),
            format = "file"
        ),

        # ---- Pool CV QC (uses pool_matrix snapshot from lipid_pre) ----
        # Returns NULL when no pools were detected pre-filter, which the
        # report skips gracefully.
        tar_target(
            lipid_pool_cv,
            {
                cv_res <- compute_pool_cv(lipid_pre$pool_matrix)
                if (is.null(cv_res)) {
                    NULL
                } else {
                    dirs <- create_legacy_output_dirs(lipid_out_dir)
                    out_qc <- dirs$diagnostic_plots
                    out_ds <- dirs$datasets

                    f_tsv <- save_tsv(cv_res$per_feature, out_ds,
                                      "pool_cv_per_feature.tsv")

                    p <- plot_pool_cv_density(cv_res)
                    f_png <- file.path(out_qc, "pool_cv_density.png")
                    if (!is.null(p)) {
                        ggplot2::ggsave(f_png, p, width = 7, height = 5,
                                        dpi = 300)
                    } else {
                        f_png <- character(0)
                    }

                    list(
                        cv_result = cv_res,
                        plot      = p,
                        files     = c(f_tsv, f_png)
                    )
                }
            }
        ),

        # ---- QC diagnostics (Stage 1) ----
        tar_target(
            lipid_qc_pre_obj,
            mod_lipidomics_qc_pre(
                pre     = lipid_pre,
                config  = config,
                out_dir = lipid_out_dir
            )
        ),

        # ---- differential expression (Stage 2) ----
        tar_target(
            lipid_de_res,
            mod_lipidomics_de(
                pre     = lipid_pre,
                config  = config,
                out_dir = lipid_out_dir
            )
        ),

        # ---- clustering: hierarchical / partition / binary (Stage 2) ----
        tar_target(
            lipid_clustering_obj,
            mod_lipidomics_clustering(
                pre     = lipid_pre,
                de_res  = lipid_de_res,
                config  = config,
                out_dir = lipid_out_dir
            )
        ),

        # ---- feature selection: RF + PLS-DA (Stage 2, optional) ----
        tar_target(
            lipid_feature_sel_res,
            mod_lipidomics_feature_selection(
                pre     = lipid_pre,
                config  = config,
                out_dir = lipid_out_dir
            )
        ),

        # ---- lipid class analysis (lipidomics-specific) ----
        tar_target(
            lipid_class_res,
            mod_lipidomics_class_analysis(
                pre     = lipid_pre,
                de_res  = lipid_de_res,
                config  = config,
                out_dir = lipid_out_dir
            )
        ),

        # ---- biomarker discovery: ROC, AUC, Cohen's d (Stage 3) ----
        tar_target(
            lipid_biomarker_res,
            mod_lipidomics_biomarker(
                pre     = lipid_pre,
                de_res  = lipid_de_res,
                config  = config,
                out_dir = lipid_out_dir
            )
        ),

        # ---- lipid pathway & network analysis (Stage 3) ----
        tar_target(
            lipid_pathway_res,
            mod_lipidomics_pathway(
                pre     = lipid_pre,
                de_res  = lipid_de_res,
                config  = config,
                out_dir = lipid_out_dir
            )
        ),

        # ---- enhanced QC: PCA loading/scree, PLS-DA CV, PERMANOVA (Stage 3) ----
        tar_target(
            lipid_qc_enhanced,
            mod_lipidomics_qc_enhanced(
                pre             = lipid_pre,
                qc_res          = lipid_qc_pre_obj,
                feature_sel_res = lipid_feature_sel_res,
                config          = config,
                out_dir         = lipid_out_dir
            )
        ),

        # ---- standalone QC report (HTML) ----
        tar_target(
            lipid_qc_report,
            mod_lipidomics_qc_report(
                pre         = lipid_pre,
                qc_res      = lipid_qc_pre_obj,
                qc_enhanced = lipid_qc_enhanced,
                config      = config,
                out_dir     = lipid_out_dir
            ),
            format = "file"
        ),

        # ---- AI commentary (depends on all analysis targets) ----
        tar_target(
            lipid_commentary,
            mod_lipidomics_commentary(
                de_res     = lipid_de_res,
                qc_pre_obj = lipid_qc_pre_obj,
                config     = config,
                out_dir    = lipid_out_dir,
                # Force all plot-producing targets to finish before scanning
                # out_dir for figures (see mod_lipidomics_commentary `deps`).
                deps       = list(
                    lipid_feature_sel_res,
                    lipid_class_res,
                    lipid_biomarker_res,
                    lipid_pathway_res,
                    lipid_clustering_obj,
                    lipid_qc_enhanced
                )
            ),
            format = "file"
        ),

        # ---- HTML report ----
        tar_target(
            lipid_report,
            mod_lipidomics_report(
                pre             = lipid_pre,
                qc_res          = lipid_qc_pre_obj,
                de_res          = lipid_de_res,
                feature_sel_res = lipid_feature_sel_res,
                class_res       = lipid_class_res,
                clustering_res  = lipid_clustering_obj,
                config          = config,
                out_dir         = lipid_out_dir,
                biomarker_res   = lipid_biomarker_res,
                pathway_res     = lipid_pathway_res,
                qc_enhanced     = lipid_qc_enhanced,
                pool_cv         = lipid_pool_cv,
                commentary_file = lipid_commentary
            ),
            format = "file"
        ),

        # ---- pipeline summary HTML ----
        tar_target(
            lipid_pipeline_summary,
            generate_lipid_pipeline_summary(
                config          = config,
                pre             = lipid_pre,
                de_res          = lipid_de_res,
                feature_sel_res = lipid_feature_sel_res,
                class_res       = lipid_class_res,
                run_dir         = lipid_out_dir
            ),
            format = "file"
        ),

        # ---- PowerPoint summary presentation ----
        tar_target(
            lipid_pptx,
            mod_lipidomics_powerpoint(
                pre             = lipid_pre,
                qc_res          = lipid_qc_pre_obj,
                de_res          = lipid_de_res,
                feature_sel_res = lipid_feature_sel_res,
                class_res       = lipid_class_res,
                config          = config,
                out_dir         = lipid_out_dir,
                biomarker_res   = lipid_biomarker_res,
                pathway_res     = lipid_pathway_res,
                qc_enhanced     = lipid_qc_enhanced,
                clustering_res  = lipid_clustering_obj,
                pool_cv         = lipid_pool_cv
            ),
            format = "file"
        ),

        # ---- Lipid <-> {gene,protein,metabolite} links (Stage 4, optional) ----
        # Cross-omics linkage via curated lipid-class -> KEGG pathway map +
        # PubMed literature co-mention (IF > 15, recent N years). Gated on
        # cfg$modes$lipidomics$links$enabled — falls through quickly when
        # disabled.
        tar_target(
            lipid_links_res,
            mod_lipid_links(
                pre          = lipid_pre,
                lipid_de_res = lipid_de_res,
                config       = config,
                out_dir      = lipid_out_dir
            ),
            format = "file"
        )
    )
}
