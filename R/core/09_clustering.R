#' Run binary-pattern clustering
#'
#' Binary pattern clustering identifies expression patterns across biological groups.
#' Each feature (gene/protein) is assigned to a binary pattern based on correlation.
#'
#' How it works:
#' 1) Groups are defined by cfg$clustering$steps$binary_patterns$group_col
#'    (e.g., "treatment" column with values: control, drugA, drugB)
#' 2) Generate all binary patterns based on number of groups
#'    - For 3 groups: 000, 001, 010, 011, 100, 101, 110, 111
#' 3) Each pattern is expanded to match the sample replicate structure
#'    (e.g., pattern "101" with 3 reps/group becomes (1,1,1,0,0,0,1,1,1))
#' 4) For each feature, compute Pearson correlation between its full sample-level
#'    expression vector and each expanded pattern (captures within-group variability)
#' 5) Apply dual count gating: "1" positions must exceed counts_cutoff_high,
#'    "0" positions must be below counts_cutoff_low (if set)
#' 6) Assign feature to pattern with highest correlation (if above corr_cutoff threshold)
#' 7) Output: heatmaps, gene lists, and statistics per pattern
#'
#' Configuration:
#' - group_col: Metadata column defining biological groups (REQUIRED)
#' - corr_cutoff: Minimum correlation to assign feature to pattern (default 0.8)
#' - counts_cutoff_high: Minimum count threshold for "on" groups (default 0)
#' - counts_cutoff_low: Maximum count threshold for "off" groups (default NULL = disabled)
#'
#' @param expr_mat_corr Feature x sample expression matrix (log-transformed) for correlations and heatmaps
#' @param expr_mat_counts Feature x sample expression matrix (counts) for gating thresholds.
#'        If NULL, uses expr_mat_corr for both (legacy behavior).
#' @param meta Sample metadata
#' @param cfg Full config list (uses clustering$steps$binary_patterns)
#' @param de_features Character vector of feature IDs to cluster
#' @param out_dir Output directory for results
#' @param corr_cutoff Minimum correlation threshold (overrides config)
#' @param counts_cutoff_high Minimum count for "1" positions (overrides config)
#' @param counts_cutoff_low Maximum count for "0" positions; NULL disables (overrides config)
#' @return List with files (paths), plots (ggplot objects), best (pattern assignments)
run_binary_patterns <- function(expr_mat_corr,
                                expr_mat_counts = NULL,
                                meta,
                                cfg,
                                de_features,
                                out_dir,
                                summary_df = NULL,
                                corr_cutoff = 0.8,
                                counts_cutoff_high = 0,
                                counts_cutoff_low = NULL, 
                                annot_context = NULL) {
  stopifnot(is.matrix(expr_mat_corr) || is.data.frame(expr_mat_corr))
  expr_mat_corr <- as.matrix(expr_mat_corr)
  stopifnot(is.data.frame(meta))
  stopifnot(is.character(de_features))
  
  de_cfg <- cfg$de %||% list()  
  
  # If no separate counts matrix provided, use the corr matrix for gating (legacy behavior)
  if (is.null(expr_mat_counts)) {
    expr_mat_counts <- expr_mat_corr
    if (any(expr_mat_counts < 0, na.rm = TRUE)) {
      warning("[binary_patterns] expr_mat_counts has negative values; counts gating expects raw counts. Consider passing expr_mat_counts separately.")
    }
  } else {
    stopifnot(is.matrix(expr_mat_counts) || is.data.frame(expr_mat_counts))
    expr_mat_counts <- as.matrix(expr_mat_counts)
    # Validate dimensions match
    if (!identical(dim(expr_mat_corr), dim(expr_mat_counts))) {
      stop("[binary_patterns] expr_mat_corr and expr_mat_counts must have identical dimensions")
    }
    if (!identical(rownames(expr_mat_corr), rownames(expr_mat_counts))) {
      stop("[binary_patterns] expr_mat_corr and expr_mat_counts must have identical row names")
    }
    if (!identical(colnames(expr_mat_corr), colnames(expr_mat_counts))) {
      stop("[binary_patterns] expr_mat_corr and expr_mat_counts must have identical column names")
    }
    message("[binary_patterns] Using separate matrices: log data for correlations, counts for gating")
  }
  
  
  # Handle NULL/NA as "disable gating" by using -Inf
  if (is.null(counts_cutoff_high) || is.na(counts_cutoff_high)) {
    counts_cutoff_high <- -Inf
    message("[binary_patterns] counts_cutoff_high=NULL/NA, disabling high counts gate")
  }
  
  # Ensure directory exists
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Use clustering$group_col (strict; errors if missing)
  group_col <- get_clustering_group_col(cfg, meta)
  message(sprintf("[binary_patterns] Using group_col: %s", group_col))
  
  sample_col <- cfg$effects$samples
  
  # --- Validations ---
  if (is.null(group_col) || !(group_col %in% colnames(meta))) {
    stop(sprintf("Binary patterns: group_col '%s' not found in meta", group_col))
  }
  if (is.null(sample_col) || !(sample_col %in% colnames(meta))) {
    stop(sprintf("Binary patterns: effects$samples column '%s' not found in meta", sample_col))
  }
  
  # Align meta order to expression columns
  samples <- colnames(expr_mat_corr)
  m_idx <- match(samples, as.character(meta[[sample_col]]))
  if (any(is.na(m_idx))) {
    stop("Binary patterns: meta is missing some samples present in expr_mat_corr colnames")
  }
  meta2 <- meta[m_idx, , drop = FALSE]
  
  # Legacy-equivalent: preserve group order as they appear in meta2 (like fct_inorder)
  grp_chr <- as.character(meta2[[group_col]])
  group_levels <- unique(grp_chr)                 # order of appearance
  groups <- factor(grp_chr, levels = group_levels)
  
  n_groups <- length(group_levels)
  
  message(sprintf("[binary_patterns] %d groups detected: %s",
                  n_groups, paste(group_levels, collapse = ", ")))
  
  # --- Early Return Fix 1: Return empty structure on failure ---
  if (n_groups < 3) {
    warning("Binary patterns: <3 groups detected; skipping.")
    return(list(files = character(0), plots = list()))
  }
  
  # Restrict to DE features present
  feats <- intersect(de_features, rownames(expr_mat_corr))
  
  # --- Early Return Fix 2: Return empty structure on failure ---
  if (length(feats) < 1) {
    warning("Binary patterns: no DE features found in expression matrix rownames")
    return(list(files = character(0), plots = list()))
  }
  
  # x = log data for correlations and heatmaps
  x <- expr_mat_corr[feats, , drop = FALSE]
  # x_counts = counts data for gating thresholds
  x_counts <- expr_mat_counts[feats, , drop = FALSE]
  
  # 1) Feature x Group means
  group_means <- sapply(group_levels, function(g) {
    cols <- which(groups == g)
    rowMeans(x[, cols, drop = FALSE], na.rm = TRUE)
  })
  group_means <- as.matrix(group_means) # features x groups
  
  # 2) Patterns (optionally exclude all-0 / all-1)
  patterns <- .get_all_binary_patterns(n_groups)
  patterns <- patterns[patterns != paste(rep("0", n_groups), collapse = "")]
  patterns <- patterns[patterns != paste(rep("1", n_groups), collapse = "")]
  patterns <- as.character(unlist(patterns))
  
  # 3) Counts gate (dual cutoff: high for "1" positions, low for "0" positions)
  # Use x_counts (counts matrix) for gating, NOT the log-transformed x
  pass_counts <- .calc_counts_gate(x_counts, groups, group_levels, patterns,
                                   counts_cutoff_high, counts_cutoff_low)
  
  # Diagnostic: How many features pass counts gate for ANY pattern?
  n_pass_any <- sum(rowSums(pass_counts) > 0)
  message(sprintf("[binary_patterns] %d/%d features pass counts gate (cutoff_high=%g, cutoff_low=%s)",
                  n_pass_any, nrow(x_counts), counts_cutoff_high,
                  if (is.null(counts_cutoff_low)) "NULL" else as.character(counts_cutoff_low)))
  
  # 4) Correlation to patterns
  # Legacy equivalence: use sample-level correlation (cor with expanded pattern)
  cor_mat <- .calc_cor_to_patterns(
    expr_mat      = x,
    patterns      = patterns,
    groups        = groups,
    group_levels  = group_levels
  )
  
  # Legacy equivalence: select best pattern based ONLY on correlation, then check counts gate
  best <- .best_pattern_legacy(cor_mat, patterns, pass_counts, corr_cutoff)
  
  # --- Writing Outputs ---
  written <- character(0)
  plots <- list()
  
  # Write summary tables (assuming save_tsv returns the file path)
  written <- c(written, save_tsv(best, out_dir, "corr_results_best_pattern.tsv"))
  
  # Stats per pattern
  stats <- as.data.frame(table(best$best_pattern), stringsAsFactors = FALSE)
  if (ncol(stats) >= 2) {
    colnames(stats) <- c("pattern", "n_features")
  }
  written <- c(written, save_tsv(stats, out_dir, "corr_stats_patterns.tsv"))
  
  # 5) Heatmaps per pattern
  for (pat in patterns) {
    
    # Features belonging to this pattern
    idx_pat <- which(best$best_pattern == pat)
    if (length(idx_pat) < 2) next
    
    feats_pat <- best$feature_id[idx_pat]
    
    # ---- LEGACY-EQUIVALENT ORDERING ----
    # Legacy sorts genes within pattern by correlation to that pattern, descending.
    # Since these genes have best_pattern==pat, best_corr is the correlation to pat.
    ord <- order(best$best_corr[idx_pat], decreasing = TRUE, na.last = TRUE)
    feats_pat <- feats_pat[ord]
    # -----------------------------------
    
    mat2plot <- x[feats_pat, , drop = FALSE]
    f_hm <- file.path(out_dir, sprintf("Heatmap_%s.png", pat))
    
    annot_df <- data.frame(Condition = groups, row.names = samples)
    
    # Pre-compute row clustering for Shiny (Professional approach)
    # Z-score the rows first (same as pheatmap with scale="row")
    mat_z <- zscore_rows(mat2plot)
    
    # Defensive: only cluster if we have enough rows
    if (nrow(mat_z) >= 2) {
      row_dists <- stats::dist(mat_z, method = "euclidean")
      row_hc <- stats::hclust(row_dists, method = "complete")
      mat_ordered <- mat_z[row_hc$order, , drop = FALSE]
    } else {
      row_hc <- NULL
      mat_ordered <- mat_z
    }
    
    # Build DE pattern row annotations (up/down per contrast)
    
    lin_fc_cutoff <- de_cfg$linear_fc_cutoff %||% 1.5
    log2fc_cutoff <- log2(lin_fc_cutoff)
    
    
    p_bin <- wrap_clustering_heatmap(
      expr_mat = mat2plot,
      meta = meta2,
      cfg = cfg,
      feature_ids = feats_pat,
      ordering = NULL,
      annotation_row_builder = TRUE,
      annotation_row_context = annot_context,
      out_file = f_hm,
      title = sprintf("Pattern %s (%d genes)", pat, length(feats_pat)),
      cluster_cols = FALSE
    )
    
    save_heatmap_to_file(p_bin, f_hm)
    
    # Store Object and Path (Professional pre-compute approach for Shiny)
    plots[[pat]] <- list(
      pheatmap = p_bin,
      mat = mat_ordered,                    # Already in clustered order
      row_order = rownames(mat_ordered),    # Ordered row names for Plotly
      col_order = colnames(mat_ordered),    # Column names (sample order)
      tree_row = row_hc                     # hclust object for dendrogram
    )
    written <- c(written, f_hm)
    
    # Gene list per pattern
    gl <- data.frame(feature_id = feats_pat, stringsAsFactors = FALSE)
    written <- c(written, save_tsv(gl, out_dir, sprintf("genes_%s.tsv", pat)))
  }
  
  # --- Final Return: List containing both files and plots ---
  return(list(
    files = unique(written),
    plots = plots,
    best = best,
    bp_pat = names(table(best$best_pattern))
  ))
}

# ---- Internal Helpers ----

# Legacy-equivalent pattern generator (depth-first order: 000, 001, 010, 011, 100, 101, 110, 111)
# This matches the legacy zero_one_sequences_helper() which appends "0" then "1" recursively.
.get_all_binary_patterns <- function(n) {
  .recursive_patterns <- function(prefix, remaining) {
    if (remaining == 0) {
      return(prefix)
    }
    c(
      .recursive_patterns(paste0(prefix, "0"), remaining - 1),
      .recursive_patterns(paste0(prefix, "1"), remaining - 1)
    )
  }
  .recursive_patterns("", n)
}

# Legacy-compatible: correlate GROUP MEANS with binary patterns
# This gives much higher correlations than sample-level because within-group
# variance is eliminated. With 3 groups, we correlate 3 values vs 3-element pattern.
.calc_cor_to_patterns_groupmeans <- function(group_means, patterns) {
  # group_means: features x groups matrix
  # patterns: character vector like c("001", "010", "011", ...)
  n_groups <- ncol(group_means)
  
  # Build pattern matrix: patterns x groups
  P <- do.call(rbind, lapply(patterns, function(pat) {
    as.numeric(strsplit(pat, "")[[1]])
  }))
  colnames(P) <- colnames(group_means)
  rownames(P) <- patterns
  
  # Correlation: features x patterns
  # For each feature (row of group_means), correlate with each pattern (row of P)
  cors <- stats::cor(t(group_means), t(P), use = "pairwise.complete.obs")
  
  colnames(cors) <- patterns
  rownames(cors) <- rownames(group_means)
  cors
}

# Sample-level correlation (kept for reference, but not used by default)
.calc_cor_to_patterns <- function(expr_mat, patterns, groups, group_levels) {
  # Build pattern matrix expanded to sample level
  # Each pattern "101" becomes (1,1,1,0,0,0,1,1,1) matching sample order
  P_expanded <- do.call(rbind, lapply(patterns, function(pat) {
    bits <- as.integer(strsplit(pat, "")[[1]])
    sample_pattern <- numeric(ncol(expr_mat))
    for (g_idx in seq_along(group_levels)) {
      cols <- which(groups == group_levels[g_idx])
      sample_pattern[cols] <- bits[g_idx]
    }
    sample_pattern
  }))
  # P_expanded: patterns x samples
  # Vectorized Pearson correlation: features x patterns
  cors <- stats::cor(t(expr_mat), t(P_expanded), use = "everything")
  
  colnames(cors) <- patterns
  rownames(cors) <- rownames(expr_mat)
  cors
}


.calc_counts_gate <- function(x, groups, group_levels, patterns,
                              counts_cutoff_high, counts_cutoff_low = NULL) {
  # returns logical matrix: features x patterns
  out <- matrix(TRUE, nrow = nrow(x), ncol = length(patterns))
  rownames(out) <- rownames(x)
  colnames(out) <- patterns
  
  for (j in seq_along(patterns)) {
    pat <- patterns[j]
    bits <- as.integer(strsplit(pat, "")[[1]])
    ones <- which(bits == 1)
    zeros <- which(bits == 0)
    
    # Gate "1" groups: ALL samples must have counts > counts_cutoff_high
    if (length(ones) == 0) {
      out[, j] <- FALSE
      next
    }
    
    cols_high <- unlist(lapply(group_levels[ones], function(g) which(groups == g)))
    sub_high <- x[, cols_high, drop = FALSE]
    ok_high <- apply(sub_high, 1, function(v) all(is.finite(v) & (v > counts_cutoff_high)))
    
    # Gate "0" groups: ALL samples must have counts < counts_cutoff_low
    if (!is.null(counts_cutoff_low) && length(zeros) > 0) {
      cols_low <- unlist(lapply(group_levels[zeros], function(g) which(groups == g)))
      sub_low <- x[, cols_low, drop = FALSE]
      ok_low <- apply(sub_low, 1, function(v) all(is.finite(v) & (v < counts_cutoff_low)))
      out[, j] <- ok_high & ok_low
    } else {
      out[, j] <- ok_high
    }
  }
  out
}

# Legacy-equivalent best pattern selection:
# 1. Choose best pattern purely from correlation (above corr_cutoff)
# 2. Then separately check if that pattern passes counts gate
# 3. Does NOT mask correlations with counts gate before selection
.best_pattern_legacy <- function(cor_mat, patterns, pass_counts, corr_cutoff) {
  feats <- rownames(cor_mat)
  n_feats <- length(feats)
  
  best_pat <- rep(NA_character_, n_feats)
  best_cor <- rep(NA_real_, n_feats)
  best_pass_counts <- rep(NA, n_feats)
  
  # For each feature, find best pattern based on correlation only
  for (i in seq_len(n_feats)) {
    row_cors <- cor_mat[i, ]
    
    # Skip if all NA
    if (all(is.na(row_cors))) next
    
    # Find max correlation (legacy: which.max returns first max in pattern order)
    max_idx <- which.max(row_cors)
    max_cor <- row_cors[max_idx]
    
    # Only assign if above cutoff
    if (!is.na(max_cor) && max_cor >= corr_cutoff) {
      best_pat[i] <- patterns[max_idx]
      best_cor[i] <- max_cor
      
      # Legacy: separately check if best pattern passes counts gate
      best_pass_counts[i] <- pass_counts[i, max_idx]
    }
  }
  
  data.frame(
    feature_id = feats,
    best_pattern = best_pat,
    best_corr = best_cor,
    best_pattern_pass_counts_cutoff = best_pass_counts,
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

# ---- Clustering group column ----

#' Get the clustering group column from config, with strict validation
#'
#' Returns \code{cfg$clustering$group_col} after checking it exists in
#' \code{meta}.  Errors if the key is missing or the column is absent.
#'
#' @param cfg  Mode config (e.g. \code{config$modes$proteomics}).
#' @param meta data.frame of sample metadata.
#' @return Character scalar: validated column name in \code{meta}.
#' @export
get_clustering_group_col <- function(cfg, meta) {
  group_col <- cfg$clustering$group_col
  if (is.null(group_col) || !nzchar(group_col)) {
    stop("clustering$group_col is required but missing or empty. ",
         "Set it in the config under clustering: group_col: \"<column_name>\"")
  }
  if (!(group_col %in% colnames(meta))) {
    stop(sprintf(
      "clustering$group_col '%s' not found in metadata columns: %s",
      group_col, paste(colnames(meta), collapse = ", ")
    ))
  }
  group_col
}

#' Build sample-level long data frame annotated with group and cluster
#'
#' Shared data-preparation helper used by both write_clustering_legacy_profiles()
#' (for data export) and save_cluster_profile_outputs() (for plotting).
#' Converts a feature x sample expression matrix into long format, joins with
#' metadata groups and cluster assignments.
#'
#' @param expr_mat Numeric matrix (features x samples)
#' @param meta Sample metadata data.frame
#' @param clusters Named integer vector (feature IDs -> cluster numbers)
#' @param group_col Character: metadata column for group (X-axis)
#' @param sample_col Character: metadata column identifying samples
#' @param color_col Character or NULL: optional metadata column for secondary grouping
#' @return data.frame with columns: Gene, Name, Exp, Group, Cluster,
#'   and optionally ColorGroup (only when color_col is not NULL)
build_cluster_long_df <- function(expr_mat, meta, clusters,
                                  group_col, sample_col,
                                  color_col = NULL) {
  meta_map <- meta |>
    dplyr::select(Name = dplyr::all_of(sample_col), Group = dplyr::all_of(group_col))
  
  if (!is.null(color_col)) {
    meta_map$ColorGroup <- meta[[color_col]]
  }
  
  norm_expr_long <- as.data.frame(expr_mat) |>
    tibble::rownames_to_column("Gene") |>
    tidyr::pivot_longer(cols = -"Gene", names_to = "Name", values_to = "Exp")
  
  df_annotated <- norm_expr_long |>
    dplyr::inner_join(meta_map, by = "Name")
  
  cluster_map <- data.frame(
    Gene = names(clusters),
    Cluster = as.integer(clusters),
    stringsAsFactors = FALSE
  )
  
  df_annotated |>
    dplyr::inner_join(cluster_map, by = "Gene")
}

# ---- Clustering guards ----

#' Count how many distinct groups exist for clustering
#'
#' Groups are derived from \code{cfg$clustering$group_col}.
#' Fails fast if \code{group_col} is missing or invalid when clustering
#' is enabled — the user must fix the config.
#'
#' @param pre pre object (must contain $meta)
#' @param cfg mode config with $clustering$group_col
#' @return integer number of groups (levels)
get_n_groups_from_effects <- function(pre, cfg) {
  stopifnot(!is.null(pre$meta))
  
  group_col <- get_clustering_group_col(cfg, pre$meta)
  
  x <- as.factor(pre$meta[[group_col]])
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
  bin_group_col <- cl$group_col
  # Only enable if both enabled flag is TRUE AND group_col is provided (non-NULL)
  bin_enabled <- isTRUE(bin_enabled && !is.null(bin_group_col))
  
  # guards for data suitability
  can_multi_group <- isTRUE(n_groups >= as.integer(min_groups))
  
  list(
    hierarchical    = hier_enabled,
    partition       = part_enabled,
    binary_patterns = isTRUE(bin_enabled && can_multi_group)
  )
}

# ---- Partition clustering ----

#' Build feature x group mean matrix using clustering$group_col + effects$samples
#'
#' @return list(group_means = matrix feature x group,
#'              groups = factor (per sample, aligned),
#'              group_levels = character,
#'              meta_aligned = data.frame)
build_group_means_from_effects <- function(expr_mat, meta, cfg) {
  stopifnot(is.matrix(expr_mat) || is.data.frame(expr_mat))
  stopifnot(is.data.frame(meta))
  expr_mat <- as.matrix(expr_mat)
  
  group_col <- get_clustering_group_col(cfg, meta)
  sample_col <- cfg$effects$samples
  
  if (is.null(sample_col) || !(sample_col %in% colnames(meta))) {
    stop(sprintf("Partition clustering: effects$samples column '%s' not found in meta", sample_col))
  }
  
  # align meta to expr columns
  samples <- colnames(expr_mat)
  idx <- match(samples, as.character(meta[[sample_col]]))
  if (any(is.na(idx))) stop("Partition clustering: meta missing samples that appear in expr_mat colnames")
  meta2 <- meta[idx, , drop = FALSE]
  
  groups <- droplevels(as.factor(meta2[[group_col]]))
  group_levels <- unique(meta2[[group_col]])
  
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

build_clustering_distance <- function(mat) {
  cmat <- stats::cor(t(mat), method = "pearson")
  cmat[is.na(cmat)] <- 0
  stats::as.dist(1 - cmat)
}
#' Choose k by gap statistic
#'
#' Uses cluster::clusGap() with hclust (Euclidean distance + Ward.D2) to find
#' optimal k. Matches Neat_RNA-Seq behavior: firstSEmax with SE.factor=1.
#'
#' @param mat_fg Feature x group z-scored matrix
#' @param k_max Maximum K to test (default 20)
#' @param B Number of bootstrap samples (default 100)
#' @return integer k
choose_k_gap_statistic <- function(mat_fg, k_max = 20, B = 100) {
  stopifnot(is.matrix(mat_fg))
  n <- nrow(mat_fg)
  if (n < 2) stop("choose_k_gap_statistic: need at least 2 features")
  k_max <- min(as.integer(k_max), n - 1L)
  if (k_max < 2) return(2L)
  
  hclust_func <- function(x, k) {
    d <- build_clustering_distance(x)
    hc <- stats::hclust(d, method = "ward.D2")
    list(cluster = stats::cutree(hc, k = k))
  }
  
  gap <- cluster::clusGap(mat_fg, FUNcluster = hclust_func, K.max = k_max, B = B)
  best_k <- cluster::maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"],
                           method = "firstSEmax", SE.factor = 1)
  
  # Guard: maxSE can return 1; force minimum of 2
  max(as.integer(best_k), 2L)
}

#' Perform partition clustering on DE features using group means (legacy-like)
#'
#' Two normalization scopes are supported via cl_cfg$zscore_scope:
#'   - "groups" (default, legacy): build group means first, then z-score the
#'     gene x group matrix. Backward-compatible.
#'   - "samples" (recommended for K = 2-4 groups): z-score the gene x sample
#'     matrix first, then build group means and feed them directly to
#'     clustering. Per-gene SD estimated from K group means is unstable when K
#'     is small; sample-level z-scoring uses many more observations and gives
#'     a noise-aware weighting.
#'
#' Under "samples", the matrix in $z_group_means is no longer row-scaled to
#' unit variance — noisier genes carry smaller magnitudes by design. Callers
#' that plot $z_group_means assuming row-unit variance (e.g. fixed [-2, +2]
#' color scale) should be reviewed before switching scope.
#'
#' @param expr_mat features x samples (imputed)
#' @param meta sample metadata
#' @param cfg mode config (uses effects + clustering$steps$partition).
#'   Recognized keys under clustering$steps$partition: algorithm, k, k_max,
#'   k_method, gap_B, nstart, distance, min_variance (variance units, default
#'   0 -> only fully constant rows are dropped), zscore_scope ("groups" |
#'   "samples"), seed (overrides the function-arg default).
#' @param de_features feature IDs to include
#' @param seed integer reproducibility seed for clusGap / kmeans (default
#'   1L). cl_cfg$seed wins over this default. Kept for reproducibility even
#'   though cluster::pam() is deterministic in this implementation.
#' @return list(algorithm, k, clusters, group_means, z_group_means, seed,
#'   zscore_scope). $group_means always holds raw group means; $z_group_means
#'   holds whatever matrix was fed to clustering for the chosen scope.
perform_partition_clustering_effects <- function(expr_mat, meta, cfg, de_features, seed = 1L) {
  stopifnot(is.character(de_features))
  expr_mat <- as.matrix(expr_mat)
  
  feats <- intersect(de_features, rownames(expr_mat))
  if (length(feats) < 2) stop("Partition clustering: need at least 2 DE features")
  
  cl_cfg <- cfg$clustering$steps$partition
  if (is.null(cl_cfg) || isFALSE(cl_cfg$enabled)) {
    stop("Partition clustering requested but config$modes$<mode>$clustering$steps$partition$enabled is FALSE")
  }
  
  # Reproducibility: cl_cfg$seed wins over the function-arg default.
  seed <- cl_cfg$seed %||% seed
  set.seed(seed)
  
  # Default algorithm is now hclust to match legacy, but supports others
  alg <- tolower(cl_cfg$algorithm %||% "hclust")
  if (!(alg %in% c("pam", "kmeans", "hclust"))) stop("partition$algorithm must be 'pam', 'kmeans' or 'hclust'")
  
  zscore_scope <- tolower(cl_cfg$zscore_scope %||% "groups")
  if (!(zscore_scope %in% c("groups", "samples"))) {
    stop("partition$zscore_scope must be 'groups' or 'samples'")
  }
  min_var <- cl_cfg$min_variance %||% 0
  
  # feature x sample subset
  x <- expr_mat[feats, , drop = FALSE]
  
  # Drop low-variance rows. Variance, not SD; threshold is also in variance
  # units so the comparison is unit-consistent. zscore_rows()'s row_sds == 0
  # safety net stays intact for other callers.
  filter_low_var <- function(m) {
    vars <- apply(m, 1, stats::var, na.rm = TRUE)
    vars[is.na(vars)] <- 0
    keep <- vars > min_var
    n_dropped <- sum(!keep)
    if (n_dropped > 0) {
      dropped_ids <- rownames(m)[!keep]
      message(sprintf(
        "partition clustering: dropped %d low-variance features (min_variance = %g)",
        n_dropped, min_var
      ))
      if (length(dropped_ids) > 0) {
        message(sprintf("  first dropped: %s",
                        paste(utils::head(dropped_ids, 5), collapse = ", ")))
      }
    }
    m[keep, , drop = FALSE]
  }
  
  if (zscore_scope == "samples") {
    # Recommended for small K: estimate per-gene SD across all samples.
    x <- filter_low_var(x)
    if (nrow(x) < 2) stop("Partition clustering: <2 features after variance filter")
    
    # group_means slot keeps raw expression semantics (un-z-scored), consistent
    # with the "groups" scope. Avoids leaking z-scored values into any future
    # downstream caller that reads $group_means.
    gm_raw_obj <- build_group_means_from_effects(x, meta, cfg)
    gm <- gm_raw_obj$group_means
    
    # z-scored samples -> group means -> fed to clustering. No second z-score
    # on the group-means matrix (that would undo the noise-aware weighting).
    x_z <- zscore_rows(x)
    gm_z_obj <- build_group_means_from_effects(x_z, meta, cfg)
    z_gm <- gm_z_obj$group_means
  } else {
    # Default "groups": legacy behavior — group means first, then z-score.
    gm_obj <- build_group_means_from_effects(x, meta, cfg)
    gm <- gm_obj$group_means
    gm <- filter_low_var(gm)
    if (nrow(gm) < 2) stop("Partition clustering: <2 features after variance filter")
    z_gm <- zscore_rows(gm)
  }
  
  # Configuration parameters
  k_fixed <- cl_cfg$k
  k_max <- cl_cfg$k_max %||% 20
  nstart <- cl_cfg$nstart %||% 25
  
  clusters <- NULL
  final_k <- NULL
  hc <- NULL
  dist_mat <- NULL
  
  # --- Determine final_k (no fitting yet) ---
  
  if (alg == "hclust") {
    # Single source of truth for the distance, paired with Ward.D2 (Euclidean).
    dist_mat <- build_clustering_distance(z_gm)
    hc <- stats::hclust(dist_mat, method = "ward.D2")
    
    if (!is.null(k_fixed)) {
      final_k <- as.integer(k_fixed)
    } else {
      k_method <- tolower(cl_cfg$k_method %||% "gap")
      
      if (k_method == "gap") {
        # Gap statistic (matches Neat_RNA-Seq: firstSEmax, SE.factor=1)
        gap_B <- cl_cfg$gap_B %||% 100
        final_k <- choose_k_gap_statistic(z_gm, k_max = k_max, B = gap_B)
      } else {
        # Silhouette on hierarchical tree cuts: explicit k -> score table,
        # no magic offsets.
        sil_table <- vapply(2:k_max, function(i) {
          ct  <- stats::cutree(hc, k = i)
          sil <- cluster::silhouette(ct, dist_mat)
          mean(sil[, "sil_width"], na.rm = TRUE)
        }, numeric(1))
        names(sil_table) <- as.character(2:k_max)
        final_k <- as.integer(names(sil_table)[which.max(sil_table)])
      }
    }
  } else {
    # PAM / K-means
    if (is.null(k_fixed)) {
      final_k <- choose_k_silhouette(z_gm, algorithm = alg, k_max = k_max, nstart = nstart)
    } else {
      final_k <- as.integer(k_fixed)
    }
  }
  
  # Single point of validation, fires for every branch (incl. k_fixed = 1).
  if (is.null(final_k) || final_k < 2) stop("partition$k must be >= 2")
  
  # --- Fit clusters ---
  
  if (alg == "hclust") {
    # hc and dist_mat are reused from the k-determination block above.
    clusters <- stats::cutree(hc, k = final_k)
  } else if (alg == "pam") {
    res <- cluster::pam(z_gm, k = final_k, metric = cl_cfg$distance %||% "euclidean")
    clusters <- res$clustering
  } else { # kmeans
    res <- stats::kmeans(z_gm, centers = final_k, nstart = nstart)
    clusters <- res$cluster
  }
  
  # Ensure names are set correctly (with defensive check)
  z_gm_rownames <- rownames(z_gm)
  if (length(clusters) != length(z_gm_rownames)) {
    stop(sprintf(
      "Partition clustering: cluster vector length (%d) does not match z_gm rows (%d)",
      length(clusters), length(z_gm_rownames)
    ))
  }
  names(clusters) <- z_gm_rownames
  
  list(
    algorithm = alg,
    k = final_k,
    clusters = clusters,
    group_means = gm,
    z_group_means = z_gm,
    seed = seed,
    zscore_scope = zscore_scope
  )
}

#' Build cluster assignment table and optionally write to TSV
#'
#' @param clusters Named integer vector (names = feature IDs, values = cluster numbers)
#' @param out_file Optional file path; if provided, writes TSV via save_tsv_path
#' @return Data frame with columns: feature_id, cluster
build_clustering_output_table <- function(clusters, out_file = NULL) {
  tbl <- data.frame(
    feature_id = names(clusters),
    cluster = as.integer(clusters),
    stringsAsFactors = FALSE
  )
  if (!is.null(out_file)) save_tsv_path(tbl, out_file)
  tbl
}

#' Compute row gap positions for pheatmap from ordered cluster assignments
#'
#' Given an integer vector of cluster assignments (in heatmap row order),
#' returns cumulative positions where gaps should appear between clusters.
#'
#' @param clusters_ordered Integer vector of cluster assignments in display order
#' @return Integer vector of gap positions (empty for k=1)
#' @export
compute_cluster_gaps <- function(clusters_ordered) {
  stopifnot(is.integer(clusters_ordered) || is.numeric(clusters_ordered))
  rl <- rle(as.integer(clusters_ordered))
  if (length(rl$lengths) <= 1) return(integer(0))
  cumsum(rl$lengths[-length(rl$lengths)])
}

#' Compute column gap positions from a sample annotation column
#'
#' Reads a group label per sample from \code{annot_col} (preserving the
#' existing column order — no reordering) and returns the cumulative
#' positions where pheatmap should draw vertical gaps between contiguous
#' group blocks.
#'
#' @param annot_col data.frame of column annotations (rownames = sample IDs).
#' @param group_col Optional name of the group column in \code{annot_col};
#'   defaults to the first column (per \code{build_heatmap_annotation_col}
#'   convention).
#' @return Integer vector of gap positions (empty if only one block).
compute_column_gaps_by_group <- function(annot_col, group_col = NULL) {
  if (is.null(group_col)) group_col <- colnames(annot_col)[1]
  groups <- as.character(annot_col[, group_col])
  group_int <- as.integer(factor(groups, levels = unique(groups)))
  compute_cluster_gaps(group_int)
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
  
  group_col  <- get_clustering_group_col(cfg, meta)
  sample_col <- cfg$effects$samples
  
  df_final <- build_cluster_long_df(expr_mat, meta, clusters, group_col, sample_col)
  
  files_written <- character(0)
  
  # 5. Write Per-Cluster Files (Raw Data)
  # Keeping raw data as is (or rounding if you want, usually raw is kept precise)
  unique_clusters <- sort(unique(df_final$Cluster))
  
  for (k in unique_clusters) {
    clus_data <- df_final|>
      dplyr::filter(Cluster == k)|>
      dplyr::select(Name, Group, Exp)
    
    fname <- file.path(out_dir, sprintf("cluster_profiles_cluster%s_data.txt", k))
    
    write.table(clus_data, fname, sep = "\t", quote = FALSE, row.names = FALSE)
    files_written <- c(files_written, fname)
  }
  
  # 6. Write Summary File (Calculated Stats)
  fname_all <- file.path(out_dir, "cluster_profiles_data.txt")
  
  summary_df <- df_final|>
    dplyr::group_by(Cluster, Group)|>
    dplyr::summarise(
      Mean = mean(Exp, na.rm = TRUE),
      SE = sd(Exp, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    )|>
    dplyr::mutate(
      Mean_SE.y    = Mean,
      Mean_SE.ymin = Mean - SE,
      Mean_SE.ymax = Mean + SE
    )|>
    # --- Rounding to 4 decimal places ---
    dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
  
  write.table(summary_df, fname_all, sep = "\t", quote = FALSE, row.names = FALSE)
  files_written <- c(files_written, fname_all)
  
  return(files_written)
}

#' Save per-cluster profile PNGs and multi-panel grid PDF
#'
#' Orchestration layer: resolves config, builds sample-level long data via
#' build_cluster_long_df(), generates per-cluster ggplots via
#' build_cluster_profile_plots(), and saves PNGs + grid PDF.
#'
#' @param expr_mat Expression matrix (features x samples)
#' @param meta Sample metadata data.frame
#' @param clusters Named integer vector (feature IDs -> cluster numbers)
#' @param cfg Mode config list
#' @param out_dir Output directory
#' @return list(files = character vector of written paths, plots = named list of ggplots)
#' @export
save_cluster_profile_outputs <- function(expr_mat, meta, clusters, cfg, out_dir) {
  requireNamespace("gridExtra", quietly = TRUE)
  
  group_col  <- get_clustering_group_col(cfg, meta)
  sample_col <- cfg$effects$samples
  color_col  <- cfg$clustering$steps$partition$color_col  # NULL if not set
  x_axis_col <- cfg$clustering$steps$partition$x_axis_col %||% group_col
  
  if (!is.null(color_col) && !(color_col %in% colnames(meta))) {
    warning(sprintf("clustering$steps$partition$color_col '%s' not found in metadata; ignoring.",
                    color_col))
    color_col <- NULL
  }
  if (!(x_axis_col %in% colnames(meta))) {
    warning(sprintf("clustering$steps$partition$x_axis_col '%s' not found in metadata; falling back to group_col.",
                    x_axis_col))
    x_axis_col <- group_col
  }
  
  long_df <- build_cluster_long_df(expr_mat, meta, clusters,
                                   x_axis_col, sample_col, color_col)
  
  color_label <- if (!is.null(color_col)) color_col else NULL
  plot_list <- build_cluster_profile_plots(long_df, x_label = x_axis_col,
                                           color_label = color_label)
  if (length(plot_list) == 0) return(list(files = character(0), plots = list()))
  
  written <- character(0)
  k <- length(plot_list)
  
  # Per-cluster PNGs (source parity: 600 dpi, 3x3 in)
  for (nm in names(plot_list)) {
    f_png <- file.path(out_dir, sprintf("cluster_profiles_cluster%s.png", nm))
    ggplot2::ggsave(f_png, plot = plot_list[[nm]],
                    dpi = 600, width = 3, height = 3, units = "in")
    written <- c(written, f_png)
  }
  
  # Multi-panel grid PDF (source parity layout)
  if (k <= 2) {
    ncol_grid <- k
    nrow_grid <- 1
  } else if (k <= 4) {
    ncol_grid <- 2
    nrow_grid <- 2
  } else {
    ncol_grid <- 3
    nrow_grid <- 2
  }
  
  f_pdf <- file.path(out_dir, "cluster_profiles.pdf")
  ml <- gridExtra::marrangeGrob(
    grobs = plot_list,
    ncol = ncol_grid,
    nrow = nrow_grid,
    top = NULL
  )
  ggplot2::ggsave(f_pdf, ml, width = 6.99, height = 3.99, dpi = 600)
  written <- c(written, f_pdf)
  
  list(files = written, plots = plot_list)
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
#' in each contrast. Non-significant genes are set to NA so pheatmap renders
#' them as white/blank. Columns with all NA (no DE genes) are removed.
#'
#' Auto-detects column style:
#' - RNA-seq:       padj.<contrast>,      linearFC.<contrast>,      <contrast>_pass
#' - Proteomics:    padj.imputs.<contrast>,linearFC.imputs.<contrast>,pass.imputs.<contrast>
#' - Metabolomics:  padj.<contrast>,      linearFC.<contrast>,      pass.<contrast>
#'
#' @param summary_df DE summary data frame
#' @param feature_ids Character vector of feature IDs to include
#' @param p_cutoff P-value cutoff (unused if pass columns exist, kept for API compat)
#' @param log2fc_cutoff log2 fold-change cutoff (unused if pass columns exist)
#'
#' @param id_col Name of the feature ID column in summary_df (default "FeatureID")
#'
#' @return Data frame with genes as rownames, contrasts as columns, values = "up"/"down"/NA
#' @export
build_de_row_annotations <- function(summary_df, feature_ids, p_cutoff, log2fc_cutoff, id_col) {
  
  
  stopifnot(is.data.frame(summary_df))
  stopifnot(id_col %in% colnames(summary_df))
  
  cols <- colnames(summary_df)
  
  # Auto-detect column style
  # TODO: refactor to parameter-driven prefix dispatch (pass prefix from calling module)
  pass_imputs_cols <- grep("^pass\\.imputs\\.", cols, value = TRUE)
  pass_suffix_cols <- grep("_pass$", cols, value = TRUE)
  pass_dot_cols    <- grep("^pass\\.", cols, value = TRUE)
  pass_dot_cols    <- setdiff(pass_dot_cols, c("pass_any_contrast", pass_imputs_cols))
  
  if (length(pass_imputs_cols) > 0) {
    # Proteomics style: pass.imputs.<contrast>, linearFC.imputs.<contrast>
    contrasts <- sub("^pass\\.imputs\\.", "", pass_imputs_cols)
    pass_prefix <- "pass.imputs."
    fc_prefix <- "linearFC.imputs."
  } else if (length(pass_suffix_cols) > 0) {
    # RNA-seq style: <contrast>_pass, linearFC.<contrast>
    contrasts <- sub("_pass$", "", pass_suffix_cols)
    pass_prefix <- NULL # special case: suffix pattern
    fc_prefix <- "linearFC."
  } else if (length(pass_dot_cols) > 0) {
    # Metabolomics style: pass.<contrast>, linearFC.<contrast>
    contrasts <- sub("^pass\\.", "", pass_dot_cols)
    pass_prefix <- "pass."
    fc_prefix <- "linearFC."
  } else {
    warning("build_de_row_annotations: No pass columns found in summary_df")
    return(NULL)
  }
  
  # Filter to features in heatmap
  sumdf_sub <- summary_df[summary_df[[id_col]] %in% feature_ids, , drop = FALSE]
  
  if (nrow(sumdf_sub) == 0) {
    warning("build_de_row_annotations: No features found in summary_df")
    return(NULL)
  }
  
  # Build annotation list
  annot_list <- list()
  
  for (cn in contrasts) {
    # Determine column names
    if (!is.null(pass_prefix)) {
      pass_col <- paste0(pass_prefix, cn)
    } else {
      pass_col <- paste0(cn, "_pass")
    }
    fc_col <- paste0(fc_prefix, cn)
    
    if (!pass_col %in% colnames(sumdf_sub) || !fc_col %in% colnames(sumdf_sub)) {
      warning(sprintf("Missing columns for contrast '%s', skipping", cn))
      next
    }
    
    pass_vals <- sumdf_sub[[pass_col]]
    fc_vals <- sumdf_sub[[fc_col]]
    
    # Determine if gene passes DE threshold
    # Proteomics: pass.imputs = 1 or NA;  RNA-seq: _pass = "up"/"down"/""
    is_de <- !is.na(pass_vals) & (pass_vals == 1 | pass_vals == TRUE | pass_vals %in% c("up", "down"))
    
    # Classify direction: NA = not significant, "up" or "down" = significant
    direction <- rep(NA_character_, nrow(sumdf_sub))
    direction[is_de & !is.na(fc_vals) & fc_vals > 0] <- "up"
    direction[is_de & !is.na(fc_vals) & fc_vals < 0] <- "down"
    
    annot_list[[cn]] <- direction
  }
  
  if (length(annot_list) == 0) {
    return(NULL)
  }
  
  # Build data frame
  annot_df <- as.data.frame(annot_list,
                            stringsAsFactors = FALSE,
                            check.names = FALSE
  )
  rownames(annot_df) <- sumdf_sub[[id_col]]
  
  # Remove columns where all values are NA (no DE genes for that contrast)
  keep <- vapply(annot_df, function(x) any(!is.na(x)), logical(1))
  annot_df <- annot_df[, keep, drop = FALSE]
  
  if (ncol(annot_df) == 0) {
    return(NULL)
  }
  
  annot_df
}

# ---- Module-level step orchestrators -------------------------------------
#
# The three `.run_*_step` helpers are the shared step bodies used by
# mod_rnaseq_clustering / mod_metabolomics_clustering / mod_proteomics_clustering.
# They take prepared inputs (expression matrix, DE feature ids, annotation
# context, etc.) and run a single clustering step end-to-end: clustering call,
# heatmap render, table writes, profile outputs. Per-omic differences live in
# the calling module (input prep + which optional fields the caller attaches
# to its return list).

#' Run the Hierarchical clustering step (shared across omics modules)
#'
#' Writes Hierarchical_DE_heatmap.png and (when k is set) Hierarchical_clusters.tsv
#' under \code{out_dir}. The caller is responsible for attaching the returned
#' fields to its module-level scaffolds (plots, objects, excel_order).
#'
#' @param expr_mat features x samples expression matrix (already coerced).
#' @param meta sample metadata data.frame.
#' @param de_features character vector of DE feature ids.
#' @param cfg mode config (uses cfg$clustering$steps$hierarchical).
#' @param annot_col column annotation data.frame (for the pheatmap payload).
#' @param annot_context list passed to wrap_clustering_heatmap as
#'   \code{annotation_row_context}.
#' @param out_dir destination directory for this step's outputs.
#' @return list(files, plot, payload, clusters, excel_order, table_df)
.run_hierarchical_step <- function(expr_mat, meta, de_features, cfg,
                                   annot_col, annot_context, out_dir) {
  ensure_dir(out_dir)
  hcfg <- cfg$clustering$steps$hierarchical %||% list()
  
  hc_res <- run_clustering(
    expr_mat    = expr_mat,
    col_data    = meta,
    de_features = de_features,
    config      = list(
      method   = "hierarchical",
      distance = hcfg$distance %||% "euclidean",
      linkage  = hcfg$linkage  %||% "complete",
      k        = hcfg$k %||% NULL
    )
  )
  
  mat_de <- expr_mat[intersect(de_features, rownames(expr_mat)), , drop = FALSE]
  z_de   <- zscore_rows(mat_de)
  colnames(z_de) <- paste0(colnames(z_de), ".zscore")
  
  ordered_row_ids <- intersect(hc_res$ordering, rownames(z_de))
  z_de_ordered    <- z_de[ordered_row_ids, , drop = FALSE]
  
  excel_order <- list(
    ordered_ids        = hc_res$ordering,
    zscore_mat         = z_de,
    partition_clusters = NULL,
    partition_k        = NULL,
    binary_best        = NULL
  )
  
  f_hm <- file.path(out_dir, "Hierarchical_DE_heatmap.png")
  p_cluster <- wrap_clustering_heatmap(
    expr_mat               = expr_mat,
    meta                   = meta,
    cfg                    = cfg,
    feature_ids            = de_features,
    ordering               = hc_res$ordering,
    annotation_row_builder = TRUE,
    annotation_row_context = annot_context,
    out_file               = f_hm,
    cluster_cols           = FALSE,
    title = sprintf("Hierarchical Clustering (%d DE features)", length(de_features))
  )
  written <- f_hm
  
  payload <- list(
    pheatmap       = p_cluster,
    mat            = z_de_ordered,
    row_order      = rownames(z_de_ordered),
    col_order      = colnames(z_de_ordered),
    annotation_col = annot_col,
    feature_ids    = de_features,
    is_zscored     = TRUE,
    cluster_cols   = FALSE,
    tree_row       = hc_res$details
  )
  
  table_df <- NULL
  clusters <- NULL
  if (!is.null(hc_res$clusters)) {
    f_tbl    <- file.path(out_dir, "Hierarchical_clusters.tsv")
    table_df <- build_clustering_output_table(hc_res$clusters, f_tbl)
    written  <- c(written, f_tbl)
    clusters <- hc_res$clusters
  }
  
  list(
    files       = written,
    plot        = p_cluster,
    payload     = payload,
    clusters    = clusters,
    excel_order = excel_order,
    table_df    = table_df
  )
}

#' Run the Partition clustering step (shared across omics modules)
#'
#' Fits partition clusters via \code{perform_partition_clustering_effects},
#' writes the partition heatmap, per-cluster heatmaps, profile outputs and
#' legacy profile exports under \code{out_dir/Partition_clustering_<k>_clusters}.
#' Column order in the heatmap is preserved from \code{annot_col} — only
#' \code{gaps_col} is computed from the existing group blocks.
#'
#' @param expr_mat features x samples expression matrix (already coerced).
#' @param meta sample metadata data.frame.
#' @param de_features character vector of DE feature ids.
#' @param cfg mode config.
#' @param annot_col column annotation data.frame.
#' @param out_dir parent directory (e.g. \code{<clustering>/Partition_clustering}).
#' @return list(files, plots, table_df, clusters, k)
.run_partition_step <- function(expr_mat, meta, de_features, cfg,
                                annot_col, out_dir) {
  ensure_dir(out_dir)
  
  part_res <- perform_partition_clustering_effects(
    expr_mat    = expr_mat,
    meta        = meta,
    cfg         = cfg,
    de_features = de_features
  )
  
  part_dir <- file.path(out_dir, sprintf("Partition_clustering_%d_clusters", part_res$k))
  ensure_dir(part_dir)
  
  written <- character(0)
  plots   <- list()
  
  # (1) clusters table
  f_tbl    <- file.path(part_dir, "partition_clusters.tsv")
  table_df <- build_clustering_output_table(part_res$clusters, f_tbl)
  written  <- c(written, f_tbl)
  
  # (2) heatmap (no column reorder; gaps_col from existing meta order)
  feats       <- names(part_res$clusters)
  valid_feats <- intersect(feats, rownames(expr_mat))
  mat_ord     <- expr_mat[valid_feats, ][order(part_res$clusters[valid_feats], valid_feats), ]
  clusters_ordered <- part_res$clusters[rownames(mat_ord)]
  
  annot_row <- data.frame(
    Cluster   = factor(paste0("C", clusters_ordered)),
    row.names = rownames(mat_ord)
  )
  
  gaps_col <- compute_column_gaps_by_group(annot_col)
  
  f_hm <- file.path(part_dir, "Partition_clustering_heatmap.png")
  p_part <- plot_heatmap_core(
    expr_mat       = mat_ord,
    annotation_col = annot_col,
    annotation_row = annot_row,
    title          = sprintf("Partition clustering (k=%d)", part_res$k),
    scale_rows     = TRUE,
    cluster_rows   = FALSE,
    cluster_cols   = FALSE,
    max_rows       = NULL,
    gaps_row       = compute_cluster_gaps(clusters_ordered),
    gaps_col       = gaps_col
  )
  plots$partition_heatmap <- p_part
  save_heatmap_to_file(p_part, f_hm)
  written <- c(written, f_hm)
  
  # (2b) per-cluster heatmaps
  per_clust_hm_files <- save_per_cluster_heatmaps(
    expr_mat       = expr_mat,
    clusters       = part_res$clusters,
    annotation_col = annot_col,
    out_dir        = part_dir
  )
  written <- c(written, per_clust_hm_files)
  
  # (3) per-cluster profile PNGs + multi-panel grid PDF
  prof_out <- save_cluster_profile_outputs(
    expr_mat = expr_mat,
    meta     = meta,
    clusters = part_res$clusters,
    cfg      = cfg,
    out_dir  = part_dir
  )
  written <- c(written, prof_out$files)
  plots$cluster_profiles <- prof_out$plots
  
  # (4) legacy per-cluster data exports
  legacy_files <- write_clustering_legacy_profiles(
    expr_mat = expr_mat,
    meta     = meta,
    clusters = part_res$clusters,
    cfg      = cfg,
    out_dir  = part_dir
  )
  written <- c(written, legacy_files)
  
  list(
    files    = written,
    plots    = plots,
    table_df = table_df,
    clusters = part_res$clusters,
    k        = part_res$k
  )
}

#' Run the Binary patterns step (shared across omics modules)
#'
#' Thin orchestrator that resolves the per-step config, ensures the output
#' directory and forwards to \code{run_binary_patterns}. \code{expr_mat_counts}
#' is supplied by the caller (RNA: counts matrix; metab/prot: NULL).
#'
#' @param expr_mat_corr log/normalised expression matrix used for correlations.
#' @param expr_mat_counts counts matrix for gating (NULL to disable).
#' @param meta sample metadata.
#' @param de_features character vector of DE feature ids.
#' @param summary_df DE summary data.frame.
#' @param cfg mode config.
#' @param annot_context list of summary_df + cutoffs + id_col for row annotations.
#' @param out_dir destination directory for binary-pattern outputs.
#' @return list(files, plots, patterns, patterns_list, binary_best)
.run_binary_patterns_step <- function(expr_mat_corr, expr_mat_counts, meta,
                                      de_features, summary_df, cfg,
                                      annot_context, out_dir) {
  ensure_dir(out_dir)
  bcfg <- cfg$clustering$steps$binary_patterns %||% list()
  
  bp_res <- run_binary_patterns(
    expr_mat_corr      = expr_mat_corr,
    expr_mat_counts    = expr_mat_counts,
    meta               = meta,
    cfg                = cfg,
    de_features        = de_features,
    out_dir            = out_dir,
    summary_df         = summary_df,
    corr_cutoff        = bcfg$corr_cutoff %||% 0.8,
    counts_cutoff_high = bcfg$counts_cutoff_high %||% bcfg$counts_cutoff %||% 0,
    counts_cutoff_low  = bcfg$counts_cutoff_low %||% NULL,
    annot_context      = annot_context
  )
  
  list(
    files         = bp_res$files %||% character(0),
    plots         = bp_res$plots %||% list(),
    patterns      = bp_res$best %||% NULL,
    patterns_list = bp_res$bp_pat %||% NULL,
    binary_best   = bp_res$best %||% NULL
  )
}