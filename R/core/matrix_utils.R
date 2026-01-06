#' Align a DE table to a reference FeatureID order (fail-fast)
#'
#' @param tab A data.frame containing a FeatureID column.
#' @param ref_ids Character vector of FeatureID in the desired order.
#' @param run_i Integer; imputation run index (for error messages).
#' @param contrast_name Character; contrast name (for error messages).
#' @param id_col Character; name of the ID column (default "FeatureID").
#' @return tab reordered to match ref_ids (same number of rows, same order).
align_de_table_by_feature_id <- function(tab,
                                         ref_ids,
                                         run_i = NA_integer_,
                                         contrast_name = NA_character_,
                                         id_col = "FeatureID") {
  if (is.null(tab)) {
    stop(sprintf("DE table is NULL (run %s, contrast '%s').", run_i, contrast_name))
  }
  if (!is.data.frame(tab)) tab <- as.data.frame(tab)
  
  if (!id_col %in% colnames(tab)) {
    stop(sprintf("DE table missing '%s' column (run %s, contrast '%s').",
                 id_col, run_i, contrast_name))
  }
  
  ids <- tab[[id_col]]
  
  # Duplicates in tab are ambiguous for match()
  if (anyDuplicated(ids)) {
    dup <- unique(ids[duplicated(ids)])
    stop(sprintf(
      "Duplicate %s values in DE table (run %s, contrast '%s'). Examples: %s",
      id_col, run_i, contrast_name, paste(head(dup, 3), collapse = ", ")
    ))
  }
  
  # Reorder by reference ids
  idx <- match(ref_ids, ids)
  
  if (anyNA(idx)) {
    missing_ids <- ref_ids[is.na(idx)]
    stop(sprintf(
      "Run %s: %d %s values missing in DE table for contrast '%s'. Examples: %s",
      run_i, length(missing_ids), id_col, contrast_name,
      paste(head(missing_ids, 3), collapse = ", ")
    ))
  }
  
  tab[idx, , drop = FALSE]
}

# Row-wise z-score of a matrix
zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds <- apply(mat, 1, stats::sd, na.rm = TRUE)
  row_sds[row_sds == 0 | is.na(row_sds)] <- 1
  
  scaled <- sweep(mat, 1, row_means, FUN = "-")
  scaled <- sweep(scaled, 1, row_sds, FUN = "/")
  scaled
}

# compute hierarchical order for rows (כמו pheatmap עם correlation + complete)
get_hclust_row_order <- function(mat, row_distance = "correlation", hclust_method = "complete") {
  mat <- as.matrix(mat)
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < 2) return(rownames(mat))
  
  d <- switch(
    row_distance,
    correlation = stats::as.dist(1 - stats::cor(t(mat), use = "pairwise.complete.obs")),
    euclidean   = stats::dist(mat),
    stop("Unsupported row_distance: ", row_distance)
  )
  
  hc <- stats::hclust(d, method = hclust_method)
  rownames(mat)[hc$order]
}


