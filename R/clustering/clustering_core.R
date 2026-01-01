# Core clustering utilities for omics-agnostic feature clustering
#
# The `run_clustering` function applies clustering methods to a set of
# differential features after row-wise z-scoring. Supported methods include
# hierarchical clustering and partitioning approaches (k-means, PAM).

#' Run clustering on differential features
#'
#' @param expr_mat Numeric matrix/data.frame with features in rows and samples in columns.
#' @param col_data Sample metadata (not directly used by clustering but retained for context).
#' @param de_features Character vector of feature identifiers to cluster.
#' @param config List with clustering configuration (e.g., method, k, distance).
#'
#' @return A list containing clustering outputs (method, clusters, ordering, details).
run_clustering <- function(expr_mat, col_data, de_features, config) {
  stopifnot(!is.null(expr_mat), !is.null(de_features), !is.null(config))

  method <- tolower(config$method %||% "hierarchical")

  expr_mat <- as.matrix(expr_mat)
  available <- intersect(de_features, rownames(expr_mat))
  if (length(available) == 0) {
    stop("No differential features found in expression matrix")
  }

  expr_sub <- expr_mat[available, , drop = FALSE]
  z_expr <- zscore_rows(expr_sub)

  if (method %in% c("hierarchical", "hclust")) {
    return(run_hierarchical_clustering(z_expr, config))
  }

  if (method %in% c("kmeans", "k-means", "partition", "pam")) {
    return(run_partition_clustering(z_expr, config))
  }

  stop(sprintf("Unsupported clustering method: %s", method))
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

run_hierarchical_clustering <- function(z_expr, config) {
  dist_method <- config$distance %||% "euclidean"
  linkage <- config$linkage %||% "complete"
  k <- config$k

  dist_mat <- stats::dist(z_expr, method = dist_method)
  hc <- stats::hclust(dist_mat, method = linkage)

  clusters <- NULL
  if (!is.null(k)) {
    clusters <- stats::cutree(hc, k = k)
    clusters <- clusters[hc$labels]
  }

  list(
    method = "hierarchical",
    clusters = clusters,
    ordering = hc$labels[hc$order],
    details = hc,
    data = list(z_scores = z_expr, samples = colnames(z_expr))
  )
}

run_partition_clustering <- function(z_expr, config) {
  partition_method <- tolower(config$method)
  k <- config$k
  if (is.null(k) || k < 2) {
    stop("Partition clustering requires a valid 'k' >= 2")
  }

  if (partition_method %in% c("pam", "partition")) {
    res <- cluster::pam(z_expr, k = k, metric = config$distance %||% "euclidean")
    clusters <- res$clustering
    ordering <- names(sort(clusters))
    details <- res
    method <- "pam"
  } else {
    nstart <- config$nstart %||% 10
    res <- stats::kmeans(z_expr, centers = k, nstart = nstart)
    clusters <- res$cluster
    ordering <- names(sort(clusters))
    details <- res
    method <- "kmeans"
  }

  list(
    method = method,
    clusters = clusters,
    ordering = ordering,
    details = details,
    data = list(z_scores = z_expr, samples = colnames(z_expr))
  )
}
