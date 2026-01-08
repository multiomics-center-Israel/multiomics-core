# R/core/omics_utils.R

get_mode_cfg <- function(config, mode) {
  cfg <- config$modes[[mode]]
  if (is.null(cfg)) stop("Missing config$modes$", mode)
  cfg
}

get_sample_col <- function(cfg) {
  # source of truth: cfg$effects$samples
  sample_col <- cfg$effects$samples %||% cfg$id_columns$sample_col
  if (is.null(sample_col) || !nzchar(sample_col)) stop("Missing sample column in cfg$effects$samples / cfg$id_columns$sample_col")
  sample_col
}

# Rename matrix columns using sample_map (from -> to).
# Unmatched columns keep original name (warn).
apply_sample_map_to_colnames <- function(mat, sample_map, map_from, map_to) {
  stopifnot(!is.null(colnames(mat)))
  check_has_cols(sample_map, c(map_from, map_to), df_name = "sample_map")
  
  raw_names <- colnames(mat)
  new_names <- sample_map[[map_to]][match(raw_names, sample_map[[map_from]])]
  
  unmatched <- is.na(new_names)
  if (any(unmatched)) {
    warning(
      "These matrix columns did not match any row in sample_map$", map_from, ": ",
      paste(head(raw_names[unmatched], 20), collapse = ", "),
      if (sum(unmatched) > 20) sprintf(" ... (+%d more)", sum(unmatched) - 20) else ""
    )
  }
  
  colnames(mat) <- ifelse(unmatched, raw_names, new_names)
  mat
}

# Reorder columns of mat to match meta sample order.
# If strict=TRUE -> fail if any meta samples are missing in mat.
align_matrix_to_meta <- function(mat, meta, sample_col, strict = TRUE) {
  check_has_cols(meta, sample_col, df_name = "meta")
  samp <- as.character(meta[[sample_col]])
  
  missing_in_mat <- setdiff(samp, colnames(mat))
  if (length(missing_in_mat) > 0) {
    msg <- paste0(
      "meta contains samples missing in matrix columns: ",
      paste(head(missing_in_mat, 10), collapse = ", "),
      if (length(missing_in_mat) > 10) sprintf(" ... (+%d more)", length(missing_in_mat) - 10) else ""
    )
    if (isTRUE(strict)) stop(msg) else warning(msg)
  }
  
  keep <- intersect(samp, colnames(mat))
  mat[, keep, drop = FALSE]
}

#' Align sample metadata to an expression matrix by sample ID
#'
#' Fail-fast helper that reorders `meta` rows to match the column order of
#' an expression matrix (`expr_mat`). This prevents silent misalignment bugs
#' in downstream modeling (e.g., limma).
#'
#' @param expr_mat A numeric matrix (features × samples) with sample IDs in `colnames(expr_mat)`.
#' @param meta A data.frame with one row per sample.
#' @param sample_col Character. Column name in `meta` holding sample IDs.
#'
#' @return A reordered `meta` data.frame with rows aligned to `colnames(expr_mat)`.
#' @export

align_meta_to_expr <- function(expr_mat, meta, cfg) {
  sample_col <- get_sample_col(cfg)
  assert_expr_meta_alignment(expr_mat, meta, cfg, strict = FALSE)
  
  mi <- match(colnames(expr_mat), meta[[sample_col]])
  if (anyNA(mi)) {
    missing <- colnames(expr_mat)[is.na(mi)]
    stop(
      "meta missing ", sum(is.na(mi)), " samples present in expr_mat (sample_col='", sample_col, "'). ",
      "Example: ", paste(head(missing, 5), collapse = ", ")
    )
  }
  
  meta[mi, , drop = FALSE]
}


#' Align differential expression results to an expression matrix
#'
#' Fail-fast helper that reorders a `topTable()` result to match the feature
#' order of the input expression matrix. This prevents silent mis-annotation
#' when binding DE statistics with feature annotations.
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

align_feature_tbl_to_mat <- function(mat, feature_tbl, feature_id_col, annot_cols)
{
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

# Create a standard omics object (generic container)
init_omics_obj <- function(mode, engine = NULL, assay_raw, col_data, row_data = NULL, assay_work = NULL, params = list()) {
  list(
    mode   = mode,
    engine = engine,
    assay  = list(
      raw  = assay_raw,
      work = assay_work
    ),
    row_data = row_data,
    col_data = col_data,
    params   = params
  )
}



