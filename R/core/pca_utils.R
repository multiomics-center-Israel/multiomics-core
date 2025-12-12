# Compute PCA on a features x samples matrix and return scores + variance explained
compute_pca_scores <- function(expr_mat, pcs = c(1, 2), center = TRUE, scale = FALSE) {
  expr_mat <- as.matrix(expr_mat)
  
  pca <- stats::prcomp(t(expr_mat), center = center, scale. = scale)
  
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  
  scores <- as.data.frame(pca$x[, pcs, drop = FALSE])
  colnames(scores) <- c("PCx", "PCy")
  scores$sample <- rownames(scores)
  
  list(
    pca = pca,
    scores = scores,
    var_expl = var_expl
  )
}
