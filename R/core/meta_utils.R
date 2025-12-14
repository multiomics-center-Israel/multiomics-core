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
