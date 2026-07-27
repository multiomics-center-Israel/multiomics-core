# R/domain/lipidomics/05_outputs_legacy.R
#
# Write standardized lipidomics outputs (expression matrices, metadata, annotations)
# for downstream integration (MOFA+, PCA, DIABLO) and the dashboard app.


#' Write standardized lipidomics outputs
#'
#' Saves the normalized expression matrix, raw matrix, metadata, and feature
#' annotations as TSV files in the Datasets/ directory — mirroring the
#' metabolomics output contract so that the dashboard integration layer
#' can discover lipidomics data for MOFA+, PCA, and DIABLO.
#'
#' @param pre Preprocessing results (requires expr_work, expr_raw, meta, row_data).
#' @param config Full pipeline config.
#' @param out_dir Output directory (will create Datasets/ inside).
#' @return Character vector of written file paths.
write_lipidomics_outputs <- function(pre, config, out_dir) {
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
    if (!is.null(pre$expr_raw)) {
        raw_df <- as.data.frame(pre$expr_raw)
        raw_df$feature_id <- rownames(pre$expr_raw)
        raw_df <- raw_df[, c("feature_id", setdiff(colnames(raw_df), "feature_id")),
                         drop = FALSE]
        written <- c(written, save_tsv(raw_df, out_ds, "expr_raw.tsv"))
    }

    # Metadata
    written <- c(written, save_tsv(pre$meta, out_ds, "metadata_aligned.tsv"))

    # Feature annotations (lipid class, chains, etc.)
    if (!is.null(pre$row_data)) {
        written <- c(written, save_tsv(pre$row_data, out_ds, "feature_annotations.tsv"))
    }

    message("lipidomics outputs: wrote ", length(written), " files to ", out_ds)
    written
}
