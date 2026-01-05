#' Compute PCA on a features x samples matrix and return scores + variance explained
#'
#' @param expr_mat numeric matrix (features x samples)
#' @param pcs integer vector of PCs to return (e.g. c(1,2) or 1:3)
#' @param center logical; center features (default TRUE)
#' @param scale logical; scale features (default FALSE)
#'
#' @return list with:
#'  - pca: prcomp object
#'  - scores: data.frame with columns PC<k> for requested pcs + sample
#'  - var_expl: numeric vector (length = nPCs) of fraction variance explained
compute_pca_scores <- function(expr_mat, pcs = c(1, 2), center = TRUE, scale = FALSE) {
  expr_mat <- as.matrix(expr_mat)
  
  pca <- stats::prcomp(t(expr_mat), center = center, scale. = scale)
  
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  
  pcs <- as.integer(pcs)
  if (any(is.na(pcs)) || any(pcs < 1)) stop("pcs must be positive integers.")
  if (max(pcs) > ncol(pca$x)) stop("pcs contains PC index larger than available PCs.")
  
  scores <- as.data.frame(pca$x[, pcs, drop = FALSE])
  colnames(scores) <- paste0("PC", pcs)
  scores$sample <- rownames(scores)
  
  list(
    pca = pca,
    scores = scores,
    var_expl = var_expl
  )
}
