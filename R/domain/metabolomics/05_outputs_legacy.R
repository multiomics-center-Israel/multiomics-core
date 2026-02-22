
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