# R/modules/proteomics/05_mod_exports.R
#
# Unified export orchestrator for proteomics pipeline.
# Consolidates all output generation to avoid duplication.


#' Unified proteomics export orchestrator
#'
#' Builds final results ONCE and generates all outputs:
#' - TSV datasets (raw, imputed, limma results)
#' - Final results TSV
#' - Excel files (ALL and DE-only with cutoffs)
#' - Shiny payload
#'
#' @param pre Preprocessing results
#' @param de_res DE analysis results
#' @param inputs Input data (contrasts, etc.)
#' @param config Full pipeline config
#' @param qc_pre_obj QC pre-processing results (with PCA)
#' @param clustering_res Clustering results
#' @param out_dir Output directory
#'
#' @return List with:
#'   - files: character vector of all written file paths
#'   - final_results: the final results data.frame
#'   - shiny_payload_file: path to Shiny RDS file
#'
#' @export
mod_proteomics_exports <- function(
    pre,
    de_res,
    inputs,
    config,
    qc_pre_obj = NULL,
    clustering_res = NULL,
    out_dir
) {
    files <- character(0)
    dirs <- create_legacy_output_dirs(out_dir)

    prot_cfg <- config$modes$proteomics %||% list()
    de_cfg <- prot_cfg$de %||% list()
    id_col <- prot_cfg$de_table$id_col %||% "FeatureID"

    # =========================================================================
    # 1. Write dataset TSVs (raw, imputed matrices)
    # =========================================================================
    files <- c(files, write_proteomics_datasets_legacy(pre, runs = NULL, config, dirs))

    # =========================================================================
    # 2. Write limma summary TSV
    # =========================================================================
    if (!is.null(de_res$summary_df)) {
        files <- c(files, write_limma_multimp_summary_legacy(de_res$summary_df, config, dirs))
    }

    # =========================================================================
    # 3. Write per-contrast limma results (wide format)
    # =========================================================================
    if (!is.null(de_res$runs_de_tables) && length(de_res$runs_de_tables) > 0) {
        contrast_names <- names(de_res$runs_de_tables[[1]])
        for (cn in contrast_names) {
            files <- c(files, write_limma_results_multimp_legacy(
                de_res = de_res,
                contrast_name = cn,
                config = config,
                dirs = dirs
            ))
        }
    }

    # =========================================================================
    # 4. Build final_results ONCE (the key optimization)
    # =========================================================================
    final_results <- NULL

    if (!is.null(inputs$contrasts) && !is.null(de_res$summary_df)) {
        if (is.null(id_col)) {
            stop("config$modes$proteomics$de_table$id_col is NULL. Check config.yaml.")
        }

        final_results <- build_final_results_proteomics(
            pre = pre,
            summary_df = de_res$summary_df,
            contrasts_df = inputs$contrasts,
            row_data = pre$row_data,
            feature_id_col = id_col
        )

        # Write final_results TSV
        files <- c(files, save_tsv(final_results, dirs$datasets, "final_results.tsv"))
    }

    # =========================================================================
    # 5. Write Excel files (using the already-built final_results)
    # =========================================================================
    if (!is.null(final_results)) {
        # Extract Excel config for enriched layout (annotation rows, sample labels)
        excel_cfg <- prot_cfg$excel %||% list()
        prot_sample_id_col <- prot_cfg$effects$samples %||%
            prot_cfg$id_columns$sample_col %||% "SampleID"

        excel_files <- write_final_results_excels_legacy_generic(
            final_results = final_results,
            config = config,
            out_dir = out_dir,
            mode = "proteomics",
            id_col = id_col,
            expr_for_de = pre$expr_imp_single,
            with_cutoffs = TRUE,
            clustering_res = clustering_res,
            sample_meta = pre$meta,
            sample_id_col = prot_sample_id_col,
            annotation_rows = excel_cfg$annotation_rows,
            sample_label_cols = excel_cfg$sample_label_cols
        )
        files <- c(files, excel_files)
    }

    # =========================================================================
    # 6. Build and save Shiny payload
    # =========================================================================
    shiny_payload_file <- file.path(out_dir, "shiny_payload_proteomics.rds")
    final_results = files[grep("Final_results_DE", files)]
    
    shiny_payload <- build_shiny_payload_proteomics(
        pre = pre,
        de_res = de_res,
        inputs = inputs,
        config = config,
        pca_res = qc_pre_obj,
        clustering_res = clustering_res,
        final_results = final_results
    )

    saveRDS(shiny_payload, shiny_payload_file)
    message("Saved proteomics Shiny payload to: ", shiny_payload_file)
    files <- c(files, shiny_payload_file)

    # =========================================================================
    # Return consolidated results
    # =========================================================================
    list(
        files = unique(files),
        shiny_payload_file = shiny_payload_file
    )
}
