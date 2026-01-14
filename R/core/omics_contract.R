# R/core/omics_contract.R

# Minimal contracts for a standard omics object.

assert_meta_contract <- function(meta, sample_col) {
  if (!is.data.frame(meta)) stop("meta must be a data.frame.")
  check_has_cols(meta, sample_col, df_name = "meta")
  samp <- as.character(meta[[sample_col]])
  if (any(!nzchar(samp))) stop("meta$", sample_col, " contains empty sample IDs.")
  if (anyDuplicated(samp) > 0) stop("meta$", sample_col, " contains duplicated sample IDs.")
  invisible(TRUE)
}

assert_expr_meta_alignment <- function(expr_mat, meta, cfg, strict = TRUE) {
  sample_col <- get_sample_col(cfg)
  
  assert_numeric_matrix(expr_mat, "expr_mat")
  assert_meta_contract(meta, sample_col)
  
  expr_ids <- colnames(expr_mat)
  meta_ids <- as.character(meta[[sample_col]])
  
  # meta -> expr (חסר ב-matrix)
  missing_in_expr <- setdiff(meta_ids, expr_ids)
  if (length(missing_in_expr) > 0) {
    stop(
      "expr_mat is missing samples from meta (sample_col='", sample_col, "'): ",
      paste(head(missing_in_expr, 10), collapse = ", "),
      if (length(missing_in_expr) > 10) sprintf(" ... (+%d more)", length(missing_in_expr) - 10) else ""
    )
  }
  
  # expr -> meta (עודפים ב-matrix)
  if (isTRUE(strict)) {
    extra_in_expr <- setdiff(expr_ids, meta_ids)
    if (length(extra_in_expr) > 0) {
      stop(
        "expr_mat has samples not present in meta (sample_col='", sample_col, "'): ",
        paste(head(extra_in_expr, 10), collapse = ", "),
        if (length(extra_in_expr) > 10) sprintf(" ... (+%d more)", length(extra_in_expr) - 10) else ""
      )
    }
  }
  
  invisible(TRUE)
}

assert_omics_obj <- function(obj, stage = c("raw", "work"), sample_col) {
  stage <- match.arg(stage)
  
  if (!is.list(obj)) stop("omics obj must be a list.")
  if (is.null(obj$assay) || !is.list(obj$assay)) stop("omics obj missing $assay list.")
  if (is.null(obj$col_data) || !is.data.frame(obj$col_data)) stop("omics obj missing $col_data (data.frame).")
  
  if (missing(sample_col) || is.null(sample_col) || !nzchar(sample_col)) {
    stop("assert_omics_obj: 'sample_col' must be provided (no default).")
  }
  
  mat <- obj$assay[[stage]]
  if (is.null(mat)) stop("omics obj missing assay$", stage)
  
  assert_numeric_matrix(mat, paste0("assay$", stage))
  assert_meta_contract(obj$col_data, sample_col)
  
  meta_ids <- as.character(obj$col_data[[sample_col]])
  if (!identical(meta_ids, colnames(mat))) {
    # give a helpful diff
    missing_in_mat <- setdiff(meta_ids, colnames(mat))
    extra_in_mat   <- setdiff(colnames(mat), meta_ids)
    
    stop(
      "assay$", stage, " and col_data are not aligned.\n",
      "- Expected colnames(mat) to match col_data[['", sample_col, "']] exactly.\n",
      if (length(missing_in_mat)) paste0("  Missing in mat: ", paste(head(missing_in_mat, 10), collapse = ", "), "\n") else "",
      if (length(extra_in_mat))   paste0("  Extra in mat: ", paste(head(extra_in_mat, 10), collapse = ", "), "\n") else "",
      "Tip: call align_meta_to_expr(mat, col_data, sample_col) upstream."
    )
  }
  
  invisible(TRUE)
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
