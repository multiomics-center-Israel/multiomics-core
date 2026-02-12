# R/pipeline/metabolomics/00_pipe_metabolomics.R
#
# Targets assembly for metabolomics:
#   Stage 1: inspect input data, normalize, and produce QC diagnostics.
#   Stage 2: differential expression, feature selection, enrichment.


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

        # ---- QC diagnostics (Stage 1) ----
        tar_target(
            metab_qc_pre_obj,
            mod_metabolomics_qc_pre(
                pre     = metab_pre,
                config  = config,
                out_dir = metab_out_dir
            )
        ),

        # ---- differential expression (Stage 2) ----
        tar_target(
            metab_de_res,
            mod_metabolomics_de(
                pre     = metab_pre,
                config  = config,
                out_dir = metab_out_dir
            )
        ),

        # ---- random forest feature importance (Stage 2, optional) ----
        tar_target(
            metab_rf_res,
            mod_metabolomics_rf(
                pre     = metab_pre,
                config  = config,
                out_dir = metab_out_dir
            )
        ),

        # ---- PLS-DA + VIP scores (Stage 2, optional) ----
        tar_target(
            metab_plsda_res,
            mod_metabolomics_plsda(
                pre     = metab_pre,
                config  = config,
                out_dir = metab_out_dir
            )
        ),

        # ---- pathway enrichment: QEA + ssGSEA (Stage 2, optional) ----
        tar_target(
            metab_enrichment_res,
            mod_metabolomics_enrichment(
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
                pre            = metab_pre,
                de_res         = metab_de_res,
                inputs         = metab_inputs,
                config         = config,
                pca_res        = metab_qc_pre_obj,
                clustering_res = NULL,
                rf_res         = metab_rf_res,
                plsda_res      = metab_plsda_res,
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
                pre            = metab_pre,
                qc_res         = metab_qc_pre_obj,
                de_res         = metab_de_res,
                rf_res         = metab_rf_res,
                plsda_res      = metab_plsda_res,
                enrichment_res = metab_enrichment_res,
                config         = config,
                out_dir        = metab_out_dir
            ),
            format = "file"
        )
    )
}


#' Write standardized metabolomics outputs (Stage 1)
#'
#' Saves the internal representation as TSV files for downstream use.
#' @return Character vector of written file paths.
write_metabolomics_outputs <- function(pre, config, out_dir) {
    dirs <- create_legacy_output_dirs(out_dir)
    out_ds <- dirs$datasets

    written <- character(0)

    # Expression matrix (normalized)
    expr_df <- as.data.frame(pre$expr_work)
    expr_df$feature_id <- rownames(pre$expr_work)
    expr_df <- expr_df[, c("feature_id", setdiff(colnames(expr_df), "feature_id")),
                        drop = FALSE]
    written <- c(written, save_tsv(expr_df, out_ds, "expr_normalized.tsv"))

    # Expression matrix (raw)
    raw_df <- as.data.frame(pre$expr_raw)
    raw_df$feature_id <- rownames(pre$expr_raw)
    raw_df <- raw_df[, c("feature_id", setdiff(colnames(raw_df), "feature_id")),
                      drop = FALSE]
    written <- c(written, save_tsv(raw_df, out_ds, "expr_raw.tsv"))

    # Metadata
    written <- c(written, save_tsv(pre$meta, out_ds, "metadata_aligned.tsv"))

    # Row data (feature annotations)
    written <- c(written, save_tsv(pre$row_data, out_ds, "feature_annotations.tsv"))

    # Sample map (if available, from CD raw parsing)
    if (!is.null(pre$sample_map)) {
        written <- c(written, save_tsv(pre$sample_map, out_ds, "sample_map_cd.tsv"))
    }

    # Normalization evaluation (if available)
    if (!is.null(pre$normalization_eval)) {
        written <- c(written, save_tsv(pre$normalization_eval, out_ds,
                                        "normalization_method_comparison.tsv"))
    }

    written
}
