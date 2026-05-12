
#' Build Final Results Table for Metabolomics
#'
#' Mirrors the proteomics structure via build_final_results_generic():
#'   feature_id | annotations | expression | DE stats | pass_any_contrast
#'
#' Required contract columns: feature_id, linearFC.<cn>, pvalue.<cn>,
#' padj.<cn>, upDown.<cn>, manual_cutoffs.<cn>, pass_any_contrast.
#'
#' Metabolomics-specific extras (from row_data): compound name, RT, m/z,
#' HMDB ID, molecular formula, etc.
#'
#' @param pre Preprocessing results (requires expr_work, row_data).
#' @param summary_df DE summary with aligned column naming (linearFC., padj., pass.).
#' @param contrast_labels Character vector of contrast labels.
#' @param row_data Optional feature annotation data.frame. Defaults to pre$row_data.
#' @param feature_id_col Name of ID column (default "feature_id").
#' @return A data.frame matching the proteomics final_results structure.
build_final_results_metabolomics <- function(
    pre,
    summary_df,
    contrast_labels,
    row_data       = NULL,
    feature_id_col = "feature_id"
) {
    # Wrap contrast labels into contrasts_df for the generic API
    contrasts_df <- data.frame(
        Contrast_name = contrast_labels,
        stringsAsFactors = FALSE
    )

    # Build annotation column mapping from row_data
    rd <- row_data %||% pre$row_data
    annot_cols <- NULL
    if (!is.null(rd)) {
        available <- setdiff(colnames(rd), feature_id_col)
        if (length(available) > 0) {
            annot_cols <- setNames(available, available)
        }
    }

    build_final_results_generic(
        summary_df     = summary_df,
        expr_df        = pre$expr_work,
        contrasts_df   = contrasts_df,
        feature_id_col = feature_id_col,
        annot_cols     = annot_cols,
        row_data       = rd,
        fc_is_signed   = TRUE,
        mode           = "metabolomics"
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

  written
}


#' Write metabolomics final results (TSV + Excel)
#'
#' Builds the final_results table from DE summary and exports:
#' - Final_results_metabolomics.tsv
#' - Final_results_ALL_P_*.xlsx (all features)
#' - Final_results_DE_P_*.xlsx (DE-only features with clustering order + z-scores)
#'
#' @param pre Preprocessing results.
#' @param de_res DE results (must contain summary_df and de_tables or contrasts).
#' @param config Full config.
#' @param out_dir Output directory.
#' @param clustering_res Optional clustering results for Excel row ordering.
#' @return Character vector of written file paths.
write_metabolomics_final_results <- function(pre, de_res, config, out_dir,
                                              clustering_res = NULL) {
    if (is.null(de_res) || is.null(de_res$summary_df)) {
        message("metabolomics final_results: skipped (no DE results)")
        return(character(0))
    }

    dirs <- create_legacy_output_dirs(out_dir)
    files <- character(0)

    # Extract contrast labels from DE results
    contrast_labels <- names(de_res$de_tables)
    if (is.null(contrast_labels) || length(contrast_labels) == 0) {
        stop("metabolomics final_results: no contrast labels found in de_res$de_tables")
    }

    # Build final_results table
    final_results <- build_final_results_metabolomics(
        pre             = pre,
        summary_df      = de_res$summary_df,
        contrast_labels = contrast_labels,
        row_data        = pre$row_data,
        feature_id_col  = "feature_id"
    )

    # Write TSV
    files <- c(files, save_tsv(final_results, dirs$datasets,
                                "Final_results_metabolomics.tsv"))

    # Extract Excel config
    excel_cfg <- config$modes$metabolomics$excel %||% list()

    # Resolve sample ID column
    sample_id_col <- config$modes$metabolomics$effects$samples %||% "sample_id"

    # Write Excel files (ALL + DE-only)
    excel_files <- tryCatch(
        write_final_results_excels_legacy_generic(
            final_results   = final_results,
            config          = config,
            out_dir         = out_dir,
            mode            = "metabolomics",
            id_col          = "feature_id",
            expr_for_de     = pre$expr_work,
            with_cutoffs    = TRUE,
            clustering_res  = clustering_res,
            sample_meta     = pre$meta,
            sample_id_col   = sample_id_col,
            annotation_rows = excel_cfg$annotation_rows,
            sample_label_cols = excel_cfg$sample_label_cols
        ),
        error = function(e) {
            warning("Excel export failed: ", e$message)
            character(0)
        }
    )
    files <- c(files, excel_files)

    message("metabolomics final_results: wrote ", length(files), " files (",
            nrow(final_results), " features)")
    unique(files)
}