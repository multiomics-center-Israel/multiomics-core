#' Build a standardized proteomics expression object (engine-aware)
#'
#' For now this supports DIANN only, but keeps a stable contract so that
#' adding MaxQuant later does not require changing downstream code.
#'
#' Contract:
#'   - assay_log2   : matrix (features × samples), log2 scale (downstream contract)
#'   - assay_linear : optional linear-scale matrix (NULL unless available/needed)
#'   - row_data     : feature annotations aligned to assay rows
#'   - col_data     : sample metadata aligned to assay columns
#'   - info         : engine/scale metadata
#'
#' @param inputs List from load_proteomics_inputs().
#' @param config Full config list.
#'
#' @return List with fields: assay_log2, assay_linear, row_data, col_data, info.
get_proteomics_expression_matrix <- function(inputs, config) {
  
  cfg <- config$modes$proteomics
  
  # Build matrix depending on engine
  if (cfg$engine == "DIANN") {
    
    # DIA-NN → expression matrix (your current builder already returns log2)
    assay_log2 <- get_measurements_per_sample_diann(
      protein    = inputs$protein,
      sample_map = inputs$sample_map,
      meta       = inputs$metadata,
      cfg        = cfg
    )
    
    assay_linear <- NULL
    
    # Feature annotations (row_data)
    row_data <- inputs$protein[
      ,
      c(cfg$id_columns$protein_id, unlist(cfg$id_columns$protein_annot)),
      drop = FALSE
    ]
    
  } else {
    stop(sprintf("Unsupported proteomics engine: %s", cfg$engine))
  }
  
  # Align col_data to assay columns (SampleID order)
  col_data <- inputs$metadata
  col_data <- col_data[match(colnames(assay_log2), col_data$SampleID), , drop = FALSE]
  stopifnot(all(col_data$SampleID == colnames(assay_log2)))
  
  list(
    assay_log2   = assay_log2,
    assay_linear = assay_linear,
    row_data     = row_data,
    col_data     = col_data,
    info = list(
      mode         = "proteomics",
      engine       = cfg$engine,
      scale_in     = cfg$scale_in %||% "log2",
      target_scale = cfg$transform$target_scale %||% "log2"
    )
  )
}
#' Align sample metadata to an expression matrix by sample ID
#'
#' Fail-fast helper that reorders `meta` rows to match the column order of
#' an expression matrix (`expr_mat`). This prevents silent misalignment bugs
#' in downstream modeling (e.g., limma).
#'
#' The function:
#' - validates that `expr_mat` is a matrix with non-empty colnames
#' - validates that `meta` is a data.frame containing `sample_col`
#' - checks that all matrix samples are present in metadata
#' - returns `meta` reordered to match `colnames(expr_mat)`
#'
#' @param expr_mat A numeric matrix (features × samples) with sample IDs in `colnames(expr_mat)`.
#' @param meta A data.frame with one row per sample.
#' @param sample_col Character. Column name in `meta` holding sample IDs.
#'
#' @return A reordered `meta` data.frame with rows aligned to `colnames(expr_mat)`.
#' @export

align_meta_to_expr <- function(expr_mat, meta, sample_col) {
  if (!is.matrix(expr_mat)) stop("expr_mat must be a matrix.")
  if (!is.data.frame(meta)) stop("meta must be a data.frame.")
  if (!(sample_col %in% colnames(meta))) stop("meta missing sample_col: ", sample_col)
  
  mi <- match(colnames(expr_mat), meta[[sample_col]])
  if (anyNA(mi)) {
    missing <- colnames(expr_mat)[is.na(mi)]
    stop(
      "meta is missing ", sum(is.na(mi)), " samples (sample_col='", sample_col, "'). ",
      "Example: ", paste(head(missing, 5), collapse = ", ")
    )
  }
  meta[mi, , drop = FALSE]
}

#' Align feature annotations to an expression matrix by feature ID
#'
#' Fail-fast helper that subsets and reorders `prot_tbl` rows to match the
#' row order of an expression matrix (`expr_mat`) using a feature ID column.
#'
#' The function:
#' - checks `protein_id_col` exists in `prot_tbl`
#' - checks all requested `annot_cols` exist in `prot_tbl`
#' - ensures every `rownames(expr_mat)` feature is present in `prot_tbl`
#' - returns an annotation data.frame aligned to `rownames(expr_mat)`
#'
#' @param expr_mat A matrix (features × samples) with feature IDs in `rownames(expr_mat)`.
#' @param prot_tbl A data.frame containing feature annotations (e.g. DIA-NN protein table).
#' @param protein_id_col Character. Column in `prot_tbl` containing the feature IDs.
#' @param annot_cols Character vector. Columns to return (should include `protein_id_col`).
#'
#' @return A data.frame of annotations with `nrow == nrow(expr_mat)` and rows aligned to `rownames(expr_mat)`.
#' @export
align_annotations_to_expr <- function(expr_mat, prot_tbl, protein_id_col, annot_cols) {
  if (!is.matrix(expr_mat)) stop("expr_mat must be a matrix.")
  if (!is.data.frame(prot_tbl)) stop("prot_tbl must be a data.frame.")
  if (!(protein_id_col %in% colnames(prot_tbl))) stop("prot_tbl missing protein_id_col: ", protein_id_col)
  
  missing_cols <- setdiff(annot_cols, colnames(prot_tbl))
  if (length(missing_cols) > 0) {
    stop("prot_tbl is missing annotation columns: ", paste(missing_cols, collapse = ", "))
  }
  
  idx <- match(rownames(expr_mat), prot_tbl[[protein_id_col]])
  if (anyNA(idx)) {
    missing_ids <- rownames(expr_mat)[is.na(idx)]
    stop(
      "prot_tbl missing ", sum(is.na(idx)), " proteins (protein_id_col='", protein_id_col, "'). ",
      "Example: ", paste(head(missing_ids, 5), collapse = ", ")
    )
  }
  
  prot_tbl[idx, annot_cols, drop = FALSE]
}
#' Align limma differential expression results to an expression matrix
#'
#' Fail-fast helper that reorders a `topTable()` result to match the feature
#' order of the input expression matrix. This prevents silent mis-annotation
#' when binding DE statistics with feature annotations.
#'
#' The function:
#' - matches `rownames(expr_mat)` against `rownames(de)`
#' - errors if any features are missing from `de`
#' - returns `de` in the same row order as `expr_mat`
#'
#' @param de A data.frame returned by `limma::topTable()` (rownames are feature IDs).
#' @param expr_mat A matrix (features × samples) with feature IDs in `rownames(expr_mat)`.
#' @param contrast_name Optional character. Used only to improve error messages.
#'
#' @return `de` reordered to match `rownames(expr_mat)`.
#' @export
align_de_to_expr <- function(de, expr_mat, contrast_name = NULL) {
  m <- match(rownames(expr_mat), rownames(de))
  if (anyNA(m)) {
    missing <- rownames(expr_mat)[is.na(m)]
    msg <- paste0(
      "topTable missing ", sum(is.na(m)), " features",
      if (!is.null(contrast_name)) paste0(" for contrast '", contrast_name, "'") else "",
      ". Example: ", paste(head(missing, 5), collapse = ", ")
    )
    stop(msg)
  }
  de[m, , drop = FALSE]
}



