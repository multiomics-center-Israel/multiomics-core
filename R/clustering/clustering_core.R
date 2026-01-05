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

# ---- Clustering guards (effects-driven; no GROUP/GROUP1) ----

#' Count how many distinct groups exist for clustering
#'
#' Groups are derived from pre$meta[[cfg$effects$color]].
#' If the column is missing or all NA -> returns 0.
#'
#' @param pre legacy-style pre object (must contain $meta)
#' @param cfg proteomics mode config (must contain $effects$color)
#' @return integer number of groups (levels)
get_n_groups_from_effects <- function(pre, cfg) {
  stopifnot(!is.null(pre$meta))
  color_col <- cfg$effects$color
  
  if (is.null(color_col) || !nzchar(color_col)) return(0L)
  if (!(color_col %in% colnames(pre$meta))) return(0L)
  
  x <- pre$meta[[color_col]]
  if (all(is.na(x))) return(0L)
  
  # treat as factor levels; if character, make factor
  x <- as.factor(x)
  nlevels(droplevels(x))
}

#' Decide which clustering steps to run (data-driven)
#'
#' Hierarchical can always run when enabled.
#' Partition + Binary patterns run only if n_groups >= min_groups (default 3).
#'
#' @param pre pre object
#' @param cfg proteomics mode config with $clustering
#' @return named logical list: hierarchical/partition/binary_patterns
clustering_run_flags <- function(pre, cfg) {
  cl <- cfg$clustering
  if (is.null(cl) || isFALSE(cl$enabled)) {
    return(list(hierarchical = FALSE, partition = FALSE, binary_patterns = FALSE))
  }
  
  # config defaults (safe)
  min_groups <- cl$min_groups %||% 3L
  n_groups <- get_n_groups_from_effects(pre, cfg)
  
  # step blocks may be missing; treat missing as enabled=FALSE unless explicitly TRUE
  steps <- cl$steps %||% list()
  hier_enabled <- isTRUE(steps$hierarchical$enabled %||% TRUE)   # default TRUE if clustering enabled
  part_enabled <- isTRUE(steps$partition$enabled %||% FALSE)
  bin_enabled  <- isTRUE(steps$binary_patterns$enabled %||% FALSE)
  
  # guards for data suitability
  can_multi_group <- isTRUE(n_groups >= as.integer(min_groups))
  
  list(
    hierarchical    = hier_enabled,
    partition       = isTRUE(part_enabled && can_multi_group),
    binary_patterns = isTRUE(bin_enabled  && can_multi_group)
  )
}

# ---- Partition clustering (legacy-like; effects-driven) ----

#' Build feature x group mean matrix using effects$color + effects$samples
#'
#' @return list(group_means = matrix feature x group,
#'              groups = factor (per sample, aligned),
#'              group_levels = character,
#'              meta_aligned = data.frame)
build_group_means_from_effects <- function(expr_mat, meta, cfg) {
  stopifnot(is.matrix(expr_mat) || is.data.frame(expr_mat))
  stopifnot(is.data.frame(meta))
  expr_mat <- as.matrix(expr_mat)
  
  group_col  <- cfg$effects$color
  sample_col <- cfg$effects$samples
  
  if (is.null(group_col) || !(group_col %in% colnames(meta))) {
    stop(sprintf("Partition clustering: effects$color column '%s' not found in meta", group_col))
  }
  if (is.null(sample_col) || !(sample_col %in% colnames(meta))) {
    stop(sprintf("Partition clustering: effects$samples column '%s' not found in meta", sample_col))
  }
  
  # align meta to expr columns
  samples <- colnames(expr_mat)
  idx <- match(samples, as.character(meta[[sample_col]]))
  if (any(is.na(idx))) stop("Partition clustering: meta missing samples that appear in expr_mat colnames")
  meta2 <- meta[idx, , drop = FALSE]
  
  groups <- droplevels(as.factor(meta2[[group_col]]))
  group_levels <- levels(groups)
  
  group_means <- sapply(group_levels, function(g) {
    cols <- which(groups == g)
    rowMeans(expr_mat[, cols, drop = FALSE], na.rm = TRUE)
  })
  group_means <- as.matrix(group_means) # features x groups
  colnames(group_means) <- group_levels
  
  list(
    group_means = group_means,
    groups = groups,
    group_levels = group_levels,
    meta_aligned = meta2
  )
}

#' Choose k by silhouette (PAM or kmeans) on feature x group matrix
#' @return integer k
choose_k_silhouette <- function(mat_fg, algorithm = c("pam","kmeans"), k_max = 20, nstart = 25) {
  algorithm <- match.arg(algorithm)
  stopifnot(is.matrix(mat_fg))
  n <- nrow(mat_fg)
  if (n < 2) stop("choose_k_silhouette: need at least 2 features")
  k_max <- min(as.integer(k_max), n - 1L)
  if (k_max < 2) return(2L)
  
  # silhouette uses dist on rows (features)
  d <- stats::dist(mat_fg)
  
  best_k <- 2L
  best_s <- -Inf
  
  for (k in 2:k_max) {
    if (algorithm == "pam") {
      cl <- cluster::pam(mat_fg, k = k)$clustering
    } else {
      cl <- stats::kmeans(mat_fg, centers = k, nstart = nstart)$cluster
    }
    
    sil <- cluster::silhouette(cl, d)
    s_mean <- mean(sil[, "sil_width"], na.rm = TRUE)
    
    if (is.finite(s_mean) && s_mean > best_s) {
      best_s <- s_mean
      best_k <- k
    }
  }
  
  best_k
}

#' Perform partition clustering on DE features using group means (legacy-like)
#'
#' @param expr_mat features x samples (imputed)
#' @param meta sample metadata
#' @param cfg mode config (uses effects + clustering$steps$partition)
#' @param de_features feature IDs to include
#' @return list(k, clusters_named, group_means, z_group_means)
perform_partition_clustering_effects <- function(expr_mat, meta, cfg, de_features) {
  stopifnot(is.character(de_features))
  expr_mat <- as.matrix(expr_mat)
  
  feats <- intersect(de_features, rownames(expr_mat))
  if (length(feats) < 2) stop("Partition clustering: need at least 2 DE features")
  
  cl_cfg <- cfg$clustering$steps$partition
  if (is.null(cl_cfg) || isFALSE(cl_cfg$enabled)) {
    stop("Partition clustering requested but config$modes$<mode>$clustering$steps$partition$enabled is FALSE")
  }
  
  # Default algorithm is now hclust to match legacy, but supports others
  alg <- tolower(cl_cfg$algorithm %||% "hclust") 
  if (!(alg %in% c("pam","kmeans", "hclust"))) stop("partition$algorithm must be 'pam', 'kmeans' or 'hclust'")
  
  # feature x sample subset
  x <- expr_mat[feats, , drop = FALSE]
  
  # group means: feature x group
  gm_obj <- build_group_means_from_effects(x, meta, cfg)
  gm <- gm_obj$group_means
  
  # z-score across groups (row-wise)
  # This scales the patterns so we focus on trends, not intensity
  z_gm <- zscore_rows(gm)
  
  # Configuration parameters
  k_fixed <- cl_cfg$k
  k_max   <- cl_cfg$k_max %||% 20
  nstart  <- cl_cfg$nstart %||% 25
  
  clusters <- NULL
  final_k  <- NULL
  
  # --- Logic Split based on Algorithm ---
  
  if (alg == "hclust") {
    # === Legacy Style: Hierarchical Clustering ===
    
    # 1. Distance: Pearson Correlation (1 - cor)
    #    We use t(z_gm) because cor() works on columns
    dist_mat <- as.dist(1 - cor(t(z_gm)))
    
    # 2. Linkage: Ward.D2 (Minimizes variance within clusters)
    hc <- stats::hclust(dist_mat, method = "ward.D2")
    
    # 3. Determine K (if not fixed)
    if (!is.null(k_fixed)) {
      final_k <- as.integer(k_fixed)
    } else {
      # Optimize K using Silhouette on the hierarchical tree cuts
      # This mimics finding the "best cut"
      sil_width <- numeric(k_max)
      
      for (i in 2:k_max) {
        ct <- stats::cutree(hc, k = i)
        # Calculate silhouette for this cut
        sil <- cluster::silhouette(ct, dist_mat)
        sil_width[i] <- mean(sil[, 3])
      }
      
      # Pick K with max silhouette (ignoring k=1 which is 0)
      final_k <- which.max(sil_width[-1]) + 1 
    }
    
    # 4. Cut the tree
    clusters <- stats::cutree(hc, k = final_k)
    
  } else {
    # === Modern Style: Partitioning (PAM / K-means) ===
    
    if (is.null(k_fixed)) {
      final_k <- choose_k_silhouette(z_gm, algorithm = alg, k_max = k_max, nstart = nstart)
    } else {
      final_k <- as.integer(k_fixed)
    }
    
    if (final_k < 2) stop("partition$k must be >= 2")
    
    if (alg == "pam") {
      res <- cluster::pam(z_gm, k = final_k, metric = cl_cfg$distance %||% "euclidean")
      clusters <- res$clustering
    } else { # kmeans
      res <- stats::kmeans(z_gm, centers = final_k, nstart = nstart)
      clusters <- res$cluster
    }
  }
  
  # Ensure names are set correctly
  names(clusters) <- rownames(z_gm)
  
  list(
    algorithm = alg,
    k = final_k,
    clusters = clusters,
    group_means = gm,
    z_group_means = z_gm
  )
}

wrap_clustering_heatmap <- function(expr_mat, meta, cfg, 
                                    feature_ids,           # The DE or Cluster genes
                                    ordering = NULL,       # Optional custom order
                                    out_file = NULL) {
  
  # 1. Filter & Order Matrix
  use_ids <- intersect(feature_ids, rownames(expr_mat))
  
  if (!is.null(ordering)) {
    # If we have a specific clustering order (e.g. from hierarchical cluster object)
    use_ids <- intersect(ordering, use_ids)
  }
  
  mat2plot <- expr_mat[use_ids, , drop = FALSE]
  
  # 2. Prepare Annotation
  annot <- data.frame(
    Condition = meta[[cfg$effects$color]],
    row.names = meta[[cfg$effects$samples]]
  )
  
  # 3. Plot
  ph <- plot_heatmap_core(
    expr_mat = mat2plot,
    annotation_col = annot,
    title = sprintf("Hierarchical Clustering (%d DE features)", nrow(mat2plot)),
    scale_rows = TRUE,
    cluster_rows = TRUE,  # Let pheatmap cluster the subset
    cluster_cols = TRUE,
    max_rows = NULL       # Don't subsample DE results! We want to see them all
  )
  
  # 4. Save & Return
  if (!is.null(out_file)) save_heatmap_to_file(ph, out_file)
  return(ph)
}


