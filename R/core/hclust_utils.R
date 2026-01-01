#' Hierarchical clustering helper
#'
#' Computes hierarchical clustering on either the sample or feature axis of an
#' expression matrix. Returns both the distance object and the corresponding
#' `hclust` result so callers can reuse the computation for visualization or
#' downstream analyses (e.g., dendrogram cuts).
#'
#' @param expr_mat Numeric matrix or data frame of expression values
#'        (features x samples).
#' @param axis Which dimension to cluster: "samples" (default) computes
#'        distances between columns; "features" clusters rows.
#' @param dist_method Distance metric passed to [stats::dist()]. Common choices
#'        include "euclidean", "manhattan", or "maximum".
#' @param linkage Linkage method passed to [stats::hclust()], e.g., "complete",
#'        "average", or "ward.D2".
#' @param k Optional number of clusters for cutting the dendrogram via
#'        [stats::cutree()].
#' @param h Optional height at which to cut the dendrogram via
#'        [stats::cutree()].
#'
#' @return List with elements:
#'   - `dist`: distance object
#'   - `hclust`: hierarchical clustering result
#'   - `order`: order of items in the dendrogram
#'   - `labels`: labels used in the clustering
#'   - `clusters`: optional vector of cluster assignments (if `k` or `h` given)
compute_hierarchical_clustering <- function(expr_mat,
                                            axis = c("samples", "features"),
                                            dist_method = "euclidean",
                                            linkage = "complete",
                                            k = NULL,
                                            h = NULL) {
  axis <- match.arg(axis)
  stopifnot(is.matrix(expr_mat) || is.data.frame(expr_mat))

  mat <- as.matrix(expr_mat)
  if (axis == "samples") mat <- t(mat)

  d <- stats::dist(mat, method = dist_method)
  hc <- stats::hclust(d, method = linkage)

  clusters <- NULL
  if (!is.null(k) && !is.null(h)) {
    stop("Specify only one of k or h for cutree")
  }

  if (!is.null(k)) {
    clusters <- stats::cutree(hc, k = k)
  } else if (!is.null(h)) {
    clusters <- stats::cutree(hc, h = h)
  }

  list(
    dist = d,
    hclust = hc,
    order = hc$order,
    labels = hc$labels,
    clusters = clusters
  )
}
