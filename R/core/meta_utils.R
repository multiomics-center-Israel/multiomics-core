# Align metadata rows to match matrix column order
align_meta_to_matrix <- function(sample_ids, meta, sample_col) {
  stopifnot(sample_col %in% colnames(meta))
  meta <- as.data.frame(meta)
  
  idx <- match(sample_ids, meta[[sample_col]])
  if (any(is.na(idx))) {
    missing <- sample_ids[is.na(idx)]
    stop("Some samples are missing from metadata[[", sample_col, "]]: ",
         paste(head(missing), collapse = ", "))
  }
  
  meta_sub <- meta[idx, , drop = FALSE]
  rownames(meta_sub) <- sample_ids
  meta_sub
}
#' Apply sample filtering rules to (meta, expr)
#'
#' @return list(meta = meta_filtered, expr = expr_filtered, info = list(n_before, n_after))
apply_sample_filter <- function(sample_col,
                                meta,
                                expr,
                                rules,
                                mode = "omics",
                                strict_cols = FALSE) {
  stopifnot(is.data.frame(meta))
  stopifnot(is.matrix(expr) || is.data.frame(expr))
  expr <- as.matrix(expr)
  
  if (is.null(rownames(meta)) || anyNA(rownames(meta))) {
    stop("[", mode, "] meta must have non-NA rownames (SampleIDs).")
  }
  if (is.null(colnames(expr))) {
    stop("[", mode, "] expr must have colnames (SampleIDs).")
  }
  
  # enforce alignment assumption
  if (!identical(colnames(expr), rownames(meta))) {
    stop("[", mode, "] apply_sample_filter requires identical colnames(expr) == rownames(meta).")
  }
  
  if (is.null(rules) || length(rules) == 0) {
    return(list(meta = meta, expr = expr, info = list(n_before = nrow(meta), n_after = nrow(meta))))
  }
  
  keep <- rep(TRUE, nrow(meta))
  
  for (col in names(rules)) {
    if (!col %in% names(meta)) {
      msg <- paste0("[", mode, "] sample_filter column not found in metadata: ", col, " (skipping)")
      if (isTRUE(strict_cols)) stop(msg) else warning(msg)
      next
    }
    
    vals <- rules[[col]]
    if (is.null(vals)) next
    vals <- unlist(vals, use.names = FALSE)
    
    keep <- keep & (meta[[col]] %in% vals)
  }
  
  n_before <- nrow(meta)
  meta2 <- meta[keep, , drop = FALSE]
  
  if (nrow(meta2) == 0) {
    stop("[", mode, "] sample_filter removed all samples. Check your rules.")
  }
  
  expr2 <- expr[, meta2[[sample_col]], drop = FALSE]
  
  message(sprintf("[%s] sample_filter kept %d/%d samples.", mode, nrow(meta2), n_before))
  
  list(meta = meta2, expr = expr2, info = list(n_before = n_before, n_after = nrow(meta2)))
}

