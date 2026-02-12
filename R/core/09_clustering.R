#' Run binary-pattern clustering
#'
#' Binary pattern clustering identifies expression patterns across biological groups.
#' Each feature (gene/protein) is assigned to a binary pattern based on correlation.
#'
#' How it works:
#' 1) Groups are defined by cfg$clustering$steps$binary_patterns$group_col
#'    (e.g., "treatment" column with values: control, drugA, drugB)
#' 2) Calculate mean expression for each feature across groups (feature x group matrix)
#' 3) Generate all binary patterns based on number of groups
#'    - For 3 groups: 000, 001, 010, 011, 100, 101, 110, 111
#'    - Excludes trivial patterns (all-0, all-1) if skip_trivial_patterns=TRUE
#' 4) For each feature, compute Pearson correlation between its group means and each pattern
#' 5) Assign feature to pattern with highest correlation (if above corr_cutoff threshold)
#' 6) Output: heatmaps, gene lists, and statistics per pattern
#'
#' Configuration:
#' - group_col: Metadata column defining biological groups (REQUIRED)
#' - corr_cutoff: Minimum correlation to assign feature to pattern (default 0.8)
#' - counts_cutoff: Minimum count threshold per group (default 0)
#'
#' @param expr_mat Feature x sample expression matrix
#' @param meta Sample metadata
#' @param cfg Full config list (uses clustering$steps$binary_patterns)
#' @param de_features Character vector of feature IDs to cluster
#' @param out_dir Output directory for results
#' @param corr_cutoff Minimum correlation threshold (overrides config)
#' @param counts_cutoff Minimum count threshold (overrides config)
#' @return List with files (paths), plots (ggplot objects), best (pattern assignments)
run_binary_patterns <- function(expr_mat,
                                meta,
                                cfg,
                                de_features,
                                out_dir,
                                corr_cutoff = 0.8,
                                counts_cutoff = 0) {
  stopifnot(is.matrix(expr_mat) || is.data.frame(expr_mat))
  expr_mat <- as.matrix(expr_mat)
  stopifnot(is.data.frame(meta))
  stopifnot(is.character(de_features))

  # Ensure directory exists
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Task 5: Use configured group_col from binary_patterns config
  # If not set, fall back to effects$color (for backward compatibility)
  bin_cfg <- cfg$clustering$steps$binary_patterns %||% list()
  group_col <- bin_cfg$group_col

  if (is.null(group_col)) {
    # Fallback to primary color for backward compatibility
    color_config <- cfg$effects$color
    group_col <- if (!is.null(color_config)) as.character(color_config[[1]]) else NULL
  }

  sample_col <- cfg$effects$samples

  # --- Validations ---
  if (is.null(group_col) || !(group_col %in% colnames(meta))) {
    stop(sprintf("Binary patterns: group_col '%s' not found in meta", group_col))
  }
  if (is.null(sample_col) || !(sample_col %in% colnames(meta))) {
    stop(sprintf("Binary patterns: effects$samples column '%s' not found in meta", sample_col))
  }

  # Align meta order to expression columns
  samples <- colnames(expr_mat)
  m_idx <- match(samples, as.character(meta[[sample_col]]))
  if (any(is.na(m_idx))) {
    stop("Binary patterns: meta is missing some samples present in expr_mat colnames")
  }
  meta2 <- meta[m_idx, , drop = FALSE]

  groups <- droplevels(as.factor(meta2[[group_col]]))
  group_levels <- levels(groups)
  n_groups <- length(group_levels)

  # --- Early Return Fix 1: Return empty structure on failure ---
  if (n_groups < 3) {
    warning("Binary patterns: <3 groups detected; skipping.")
    return(list(files = character(0), plots = list()))
  }

  # Restrict to DE features present
  feats <- intersect(de_features, rownames(expr_mat))

  # --- Early Return Fix 2: Return empty structure on failure ---
  if (length(feats) < 1) {
    warning("Binary patterns: no DE features found in expression matrix rownames")
    return(list(files = character(0), plots = list()))
  }

  x <- expr_mat[feats, , drop = FALSE]

  # 1) Feature x Group means
  group_means <- sapply(group_levels, function(g) {
    cols <- which(groups == g)
    rowMeans(x[, cols, drop = FALSE], na.rm = TRUE)
  })
  group_means <- as.matrix(group_means) # features x groups

  # 2) Patterns (exclude all-0 / all-1)
  patterns <- .get_all_binary_patterns(n_groups)
  patterns <- patterns[patterns != paste(rep("0", n_groups), collapse = "")]
  patterns <- patterns[patterns != paste(rep("1", n_groups), collapse = "")]

  # 3) Counts gate
  pass_counts <- .calc_counts_gate(x, groups, group_levels, patterns, counts_cutoff)

  # 4) Correlation to patterns (Vectorized!)
  cor_mat <- .calc_cor_to_patterns(group_means, patterns)

  # Best pattern selection (Robust Logic)
  best <- .best_pattern_per_feature(cor_mat, patterns, pass_counts, corr_cutoff)

  # --- Writing Outputs ---
  written <- character(0)
  plots <- list()

  # Write summary tables (assuming save_tsv returns the file path)
  written <- c(written, save_tsv(best, out_dir, "corr_results_best_pattern.tsv"))

  # Stats per pattern
  stats <- as.data.frame(table(best$best_pattern), stringsAsFactors = FALSE)
  colnames(stats) <- c("pattern", "n_features")
  written <- c(written, save_tsv(stats, out_dir, "corr_stats_patterns.tsv"))

  # 5) Heatmaps per pattern
  for (pat in patterns) {
    # Filter features belonging to this pattern
    feats_pat <- best$feature_id[which(best$best_pattern == pat)]
    if (length(feats_pat) < 2) next

    mat2plot <- x[feats_pat, , drop = FALSE]
    f_hm <- file.path(out_dir, sprintf("Heatmap_%s.png", pat))

    annot_df <- data.frame(Condition = groups, row.names = samples)

    # Create Object
    p <- plot_heatmap_core(
      expr_mat       = mat2plot,
      annotation_col = annot_df,
      max_rows       = NULL,
      main           = sprintf("Pattern %s (%d features)", pat, length(feats_pat)),
      scale_rows     = TRUE,
      cluster_rows   = TRUE,
      cluster_cols   = TRUE
    )

    # Save File
    save_heatmap_to_file(p, f_hm)

    # Store Object and Path
    plots[[paste0("pattern_", pat)]] <- p
    written <- c(written, f_hm)

    # Gene list per pattern
    gl <- data.frame(feature_id = feats_pat, stringsAsFactors = FALSE)
    written <- c(written, save_tsv(gl, out_dir, sprintf("genes_%s.tsv", pat)))
  }

  # --- Final Return: List containing both files and plots ---
  return(list(
    files = unique(written),
    plots = plots,
    best = best
  ))
}

# ---- Internal Helpers ----

.get_all_binary_patterns <- function(n) {
  grid <- expand.grid(rep(list(c(0, 1)), n), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  apply(grid, 1, paste0, collapse = "")
}

.calc_cor_to_patterns <- function(group_means, patterns) {
  P <- do.call(rbind, strsplit(patterns, split = ""))
  # Ensures numeric matrix even if 1 row or 1 column
  P_mat <- matrix(as.numeric(P), nrow = nrow(P), ncol = ncol(P))

  # Vectorized Correlation: cor(Features_Transposed, Patterns_Transposed)
  # Result: Features x Patterns
  cors <- stats::cor(t(group_means), t(P_mat), use = "pairwise.complete.obs")

  colnames(cors) <- patterns
  rownames(cors) <- rownames(group_means)

  cors
}

.calc_counts_gate <- function(x, groups, group_levels, patterns, counts_cutoff) {
  # returns logical matrix: features x patterns
  out <- matrix(TRUE, nrow = nrow(x), ncol = length(patterns))
  rownames(out) <- rownames(x)
  colnames(out) <- patterns

  for (j in seq_along(patterns)) {
    pat <- patterns[j]
    bits <- as.integer(strsplit(pat, "")[[1]])
    ones <- which(bits == 1)

    if (length(ones) == 0) {
      out[, j] <- FALSE
      next
    }

    cols <- unlist(lapply(group_levels[ones], function(g) which(groups == g)))
    sub <- x[, cols, drop = FALSE]

    # Strict check: ALL samples in "1" groups must be valid and > cutoff
    ok <- apply(sub, 1, function(v) all(is.finite(v) & (v > counts_cutoff)))
    out[, j] <- ok
  }
  out
}

.best_pattern_per_feature <- function(cor_mat, patterns, pass_counts, corr_cutoff) {
  feats <- rownames(cor_mat)

  # Apply gate: set correlation to NA if counts check failed
  cor_mat[!pass_counts] <- NA_real_

  # 1. Calculate max values (suppressWarnings for -Inf when all NA)
  max_vals <- suppressWarnings(apply(cor_mat, 1, max, na.rm = TRUE))

  # 2. Determine validity (Avoids -Inf and NA issues)
  is_valid <- !is.infinite(max_vals) & !is.na(max_vals) & (max_vals >= corr_cutoff)

  best_pat <- rep(NA_character_, length(feats))
  best_cor <- rep(NA_real_, length(feats))

  # 3. Only find indices for valid rows
  if (any(is_valid)) {
    # Subset to valid rows only
    valid_mat <- cor_mat[is_valid, , drop = FALSE]

    # max.col is faster and safer than apply(which.max) for matrices
    # It handles ties deterministically ("first") and doesn't return list
    best_inds <- max.col(valid_mat, ties.method = "first")

    best_pat[is_valid] <- patterns[best_inds]
    best_cor[is_valid] <- max_vals[is_valid]
  }

  data.frame(
    feature_id = feats,
    best_pattern = best_pat,
    best_corr = best_cor,
    stringsAsFactors = FALSE
  )
}


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

  ordering <- hc$labels[hc$order]

  clusters <- NULL
  if (!is.null(k)) {
    clusters <- stats::cutree(hc, k = k)
    clusters <- clusters[hc$labels]
  }

  # ---- NEW: excel_order payload ----
  z_for_excel <- z_expr
  colnames(z_for_excel) <- paste0(colnames(z_for_excel), ".zscore")


  list(
    method = "hierarchical",
    clusters = clusters,
    ordering = ordering,
    excel_order = list(
      ordered_ids = ordering,
      zscore_mat  = z_for_excel
    ),
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

  # Extract primary color (handle array config for multi-color PCA)
  color_config <- cfg$effects$color
  if (is.null(color_config)) {
    return(0L)
  }
  color_col <- as.character(color_config[[1]])

  if (!nzchar(color_col)) {
    return(0L)
  }
  if (!(color_col %in% colnames(pre$meta))) {
    return(0L)
  }

  x <- pre$meta[[color_col]]
  if (all(is.na(x))) {
    return(0L)
  }

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

  hier_enabled <- isTRUE(steps$hierarchical$enabled %||% TRUE)
  part_enabled <- isTRUE(steps$partition$enabled %||% FALSE)

  # Task 5: Binary clustering conditional on group_col
  # If group_col is NULL/missing, don't perform binary clustering
  bin_cfg <- steps$binary_patterns %||% list()
  bin_enabled <- isTRUE(bin_cfg$enabled %||% FALSE)
  bin_group_col <- bin_cfg$group_col
  # Only enable if both enabled flag is TRUE AND group_col is provided (non-NULL)
  bin_enabled <- isTRUE(bin_enabled && !is.null(bin_group_col))

  # guards for data suitability
  can_multi_group <- isTRUE(n_groups >= as.integer(min_groups))

  list(
    hierarchical    = hier_enabled,
    partition       = isTRUE(part_enabled && can_multi_group),
    binary_patterns = isTRUE(bin_enabled && can_multi_group)
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

  # Extract primary color (handle array config for multi-color PCA)
  color_config <- cfg$effects$color
  group_col <- if (!is.null(color_config)) as.character(color_config[[1]]) else NULL
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
choose_k_silhouette <- function(mat_fg, algorithm = c("pam", "kmeans"), k_max = 20, nstart = 25) {
  algorithm <- match.arg(algorithm)
  stopifnot(is.matrix(mat_fg))
  n <- nrow(mat_fg)
  if (n < 2) stop("choose_k_silhouette: need at least 2 features")
  k_max <- min(as.integer(k_max), n - 1L)
  if (k_max < 2) {
    return(2L)
  }

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
  if (!(alg %in% c("pam", "kmeans", "hclust"))) stop("partition$algorithm must be 'pam', 'kmeans' or 'hclust'")

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
  k_max <- cl_cfg$k_max %||% 20
  nstart <- cl_cfg$nstart %||% 25

  clusters <- NULL
  final_k <- NULL

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
#' Write cluster data in exact legacy format
#' Columns: Name (Sample), Group, Exp (Absolute Expression)
#' Summary File Columns: Cluster, Group, Mean, SE, Mean_SE.y, Mean_SE.ymin, Mean_SE.ymax
#'
#' @param expr_mat Original expression matrix (log-intensities), NOT z-scores
#' @param meta Metadata dataframe
#' @param clusters Named integer vector of clusters
#' @param cfg Config list
#' @param out_dir Output directory
write_clustering_legacy_profiles <- function(expr_mat, meta, clusters, cfg, out_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # 1. Prepare Metadata Map (Sample -> Group)
  # Extract primary color (handle array config for multi-color PCA)
  color_config <- cfg$effects$color
  group_col <- if (!is.null(color_config)) as.character(color_config[[1]]) else NULL
  sample_col <- cfg$effects$samples

  meta_map <- meta %>%
    dplyr::select(Name = all_of(sample_col), Group = all_of(group_col))

  # 2. Convert Expression Matrix to Long Format
  # Rows = Genes, Cols = Samples -> Melt
  norm_expr_long <- as.data.frame(expr_mat) %>%
    tibble::rownames_to_column("Gene") %>%
    tidyr::pivot_longer(cols = -Gene, names_to = "Name", values_to = "Exp")

  # 3. Join with Metadata
  df_annotated <- norm_expr_long %>%
    dplyr::inner_join(meta_map, by = "Name")

  # 4. Map Genes to Clusters
  cluster_map <- data.frame(
    Gene = names(clusters),
    Cluster = as.integer(clusters),
    stringsAsFactors = FALSE
  )

  # Final Join: Only keep genes that are in a cluster
  df_final <- df_annotated %>%
    dplyr::inner_join(cluster_map, by = "Gene")

  files_written <- character(0)

  # 5. Write Per-Cluster Files (Raw Data)
  # Keeping raw data as is (or rounding if you want, usually raw is kept precise)
  unique_clusters <- sort(unique(df_final$Cluster))

  for (k in unique_clusters) {
    clus_data <- df_final %>%
      dplyr::filter(Cluster == k) %>%
      dplyr::select(Name, Group, Exp)

    fname <- file.path(out_dir, sprintf("cluster_profiles_cluster%s_data.txt", k))

    write.table(clus_data, fname, sep = "\t", quote = FALSE, row.names = FALSE)
    files_written <- c(files_written, fname)
  }

  # 6. Write Summary File (Calculated Stats)
  fname_all <- file.path(out_dir, "cluster_profiles_data.txt")

  summary_df <- df_final %>%
    dplyr::group_by(Cluster, Group) %>%
    dplyr::summarise(
      Mean = mean(Exp, na.rm = TRUE),
      SE = sd(Exp, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Mean_SE.y    = Mean,
      Mean_SE.ymin = Mean - SE,
      Mean_SE.ymax = Mean + SE
    ) %>%
    # --- Rounding to 4 decimal places ---
    dplyr::mutate(across(where(is.numeric), ~ round(., 4)))

  write.table(summary_df, fname_all, sep = "\t", quote = FALSE, row.names = FALSE)
  files_written <- c(files_written, fname_all)

  return(files_written)
}


#' Row-wise z-score of a matrix
#'
#' @param mat Numeric matrix (features x samples)
#' @return Numeric matrix of same dimensions with rows scaled to mean 0, sd 1
zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds <- apply(mat, 1, stats::sd, na.rm = TRUE)
  row_sds[row_sds == 0 | is.na(row_sds)] <- 1

  scaled <- sweep(mat, 1, row_means, FUN = "-")
  scaled <- sweep(scaled, 1, row_sds, FUN = "/")
  scaled
}

#' Build row annotations for DE heatmap showing up/down patterns per contrast
#'
#' Creates a data frame with genes as rows and contrasts as columns,
#' showing which genes are up-regulated, down-regulated, or not significant
#' in each contrast.
#'
#' @param summary_df DE summary data frame with columns like padj.<contrast>, log2FoldChange.<contrast>
#' @param feature_ids Character vector of feature IDs to include in heatmap
#' @param p_cutoff P-value cutoff (default 0.05)
#' @param log2fc_cutoff log2 fold-change cutoff (default log2(1.5) = 0.585)
#'
#' @return Data frame with genes as rownames, contrasts as columns, values = "up"/"down"/"ns"
#' @export
build_de_row_annotations <- function(summary_df, feature_ids, p_cutoff = 0.05, log2fc_cutoff = 0.585) {
  stopifnot(is.data.frame(summary_df))
  stopifnot("FeatureID" %in% colnames(summary_df))

  # Identify contrast columns from padj.* pattern
  cols <- colnames(summary_df)
  padj_cols <- grep("^padj\\.", cols, value = TRUE)

  if (length(padj_cols) == 0) {
    warning("build_de_row_annotations: No padj.* columns found in summary_df")
    return(NULL)
  }

  # Extract contrast names
  contrasts <- sub("^padj\\.", "", padj_cols)

  # Filter summary_df to only genes in feature_ids
  sumdf_sub <- summary_df[summary_df$FeatureID %in% feature_ids, , drop = FALSE]

  if (nrow(sumdf_sub) == 0) {
    warning("build_de_row_annotations: No features found in summary_df")
    return(NULL)
  }

  # Initialize annotation data frame
  annot_list <- list()

  for (cn in contrasts) {
    padj_col <- paste0("padj.", cn)
    fc_col <- paste0("log2FoldChange.", cn)

    # Check if FC column exists
    if (!fc_col %in% colnames(sumdf_sub)) {
      warning(sprintf("Missing log2FoldChange column for contrast '%s', skipping", cn))
      next
    }

    padj_vals <- sumdf_sub[[padj_col]]
    fc_vals <- sumdf_sub[[fc_col]]

    # Classify each gene
    pattern <- rep("ns", nrow(sumdf_sub))
    pattern[!is.na(padj_vals) & !is.na(fc_vals) &
            padj_vals <= p_cutoff &
            fc_vals > log2fc_cutoff] <- "up"
    pattern[!is.na(padj_vals) & !is.na(fc_vals) &
            padj_vals <= p_cutoff &
            fc_vals < -log2fc_cutoff] <- "down"

    annot_list[[cn]] <- pattern
  }

  if (length(annot_list) == 0) {
    return(NULL)
  }

  # Build data frame
  annot_df <- as.data.frame(annot_list, stringsAsFactors = FALSE)
  rownames(annot_df) <- sumdf_sub$FeatureID

  # Convert to factor for pheatmap
  for (col in colnames(annot_df)) {
    annot_df[[col]] <- factor(annot_df[[col]], levels = c("down", "ns", "up"))
  }

  annot_df
}
