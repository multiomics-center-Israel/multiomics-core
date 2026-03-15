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
                                counts_cutoff_low = NULL) {
  stopifnot(is.matrix(expr_mat_corr) || is.data.frame(expr_mat_corr))
  expr_mat_corr <- as.matrix(expr_mat_corr)
  stopifnot(is.data.frame(meta))
  stopifnot(is.character(de_features))
  
  de_cfg <- cfg$modes$rna$de %||% list()
  

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

  # Use clustering$group_col (no fallback to effects$color)
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
    
    
    annot_context <- list(
      summary_df    = summary_df,
      p_cutoff      = de_cfg$p_cutoff %||% 0.05,
      log2fc_cutoff = log2fc_cutoff,
      id_col        = "FeatureID"
    )
    
    
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

# ---- Clustering group column (decoupled from effects$color) ----

#' Get the clustering group column from config, with strict validation
#'
#' Returns \code{cfg$clustering$group_col} after checking it exists in
#' \code{meta}.  No fallback to \code{effects$color} — the config must
#' set \code{clustering$group_col} explicitly.
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

#' Choose k by gap statistic
#'
#' Uses cluster::clusGap() with hclust (Pearson correlation + Ward.D2) to find
#' optimal k. Matches Neat_RNA-Seq behavior: firstSEmax with SE.factor=1.
#'
#' @param mat_fg Feature x group z-scored matrix
#' @param k_max Maximum K to test (default 20)
#' @param B Number of bootstrap samples (default 10)
#' @return integer k
choose_k_gap_statistic <- function(mat_fg, k_max = 20, B = 10) {
  stopifnot(is.matrix(mat_fg))
  n <- nrow(mat_fg)
  if (n < 2) stop("choose_k_gap_statistic: need at least 2 features")
  k_max <- min(as.integer(k_max), n - 1L)
  if (k_max < 2) return(2L)
  
  hclust_func <- function(x, k) {
    cmat <- stats::cor(t(x), use = "pairwise.complete.obs")
    cmat[is.na(cmat)] <- 0
    d <- stats::as.dist(1 - cmat)
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
      k_method <- tolower(cl_cfg$k_method %||% "gap")

      if (k_method == "gap") {
        # Gap statistic (matches Neat_RNA-Seq: firstSEmax, SE.factor=1)
        gap_B <- cl_cfg$gap_B %||% 10
        final_k <- choose_k_gap_statistic(z_gm, k_max = k_max, B = gap_B)
      } else {
        # Default: Silhouette on hierarchical tree cuts
        sil_width <- numeric(k_max)

        for (i in 2:k_max) {
          ct <- stats::cutree(hc, k = i)
          sil <- cluster::silhouette(ct, dist_mat)
          sil_width[i] <- mean(sil[, 3])
        }

        final_k <- which.max(sil_width[-1]) + 1
      }
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
    z_group_means = z_gm
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

#' Build cluster profile data frame from z-scored group means
#'
#' Computes per-cluster mean and SD of z-scored group means. Handles empty
#' clusters and NA group names safely.
#'
#' @param z_group_means Matrix (features x groups) of z-scored group means
#' @param clusters Named integer vector mapping feature IDs to cluster numbers
#' @param k Number of clusters
#' @return Data frame with columns: cluster, group, mean, sd, n_features.
#'         Returns NULL if all clusters are empty.
build_cluster_profiles <- function(z_group_means, clusters, k) {
  clv <- clusters[rownames(z_group_means)]
  groups <- colnames(z_group_means)

  prof_list <- lapply(seq_len(k), function(ci) {
    rows <- which(clv == ci)
    sub <- z_group_means[rows, , drop = FALSE]
    if (nrow(sub) == 0) {
      return(NULL)
    }
    data.frame(
      cluster = ci,
      group = groups,
      mean = colMeans(sub, na.rm = TRUE),
      sd = apply(sub, 2, stats::sd, na.rm = TRUE),
      n_features = nrow(sub),
      stringsAsFactors = FALSE
    )
  })

  prof_list <- Filter(Negate(is.null), prof_list)
  if (length(prof_list) == 0) {
    return(NULL)
  }

  prof <- do.call(rbind, prof_list)
  if (any(is.na(prof$group))) prof$group[is.na(prof$group)] <- "NA"
  prof
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
  group_col <- get_clustering_group_col(cfg, meta)
  sample_col <- cfg$effects$samples

  meta_map <- meta |>
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
build_de_row_annotations <- function(summary_df, feature_ids, p_cutoff = 0.05, log2fc_cutoff = 0.585, id_col = "FeatureID") {
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
