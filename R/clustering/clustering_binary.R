#' Run binary-pattern clustering
#'
#' Logic:
#' 1) Aggregate expression by group -> feature x group means.
#' 2) Enumerate binary patterns.
#' 3) Compute correlation (vectorized).
#' 4) Filter and write results.
#'
#' @return character vector of written file paths
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
  
  group_col  <- cfg$effects$color
  sample_col <- cfg$effects$samples
  
  # --- Validations ---
  if (is.null(group_col) || !(group_col %in% colnames(meta))) {
    stop(sprintf("Binary patterns: effects$color column '%s' not found in meta", group_col))
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
  
  # Write summary tables (assuming do_write returns the file path)
  written <- c(written, do_write(best, out_dir, "corr_results_best_pattern.tsv"))
  
  # Stats per pattern
  stats <- as.data.frame(table(best$best_pattern), stringsAsFactors = FALSE)
  colnames(stats) <- c("pattern", "n_features")
  written <- c(written, do_write(stats, out_dir, "corr_stats_patterns.tsv"))
  
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
    written <- c(written, do_write(gl, out_dir, sprintf("genes_%s.tsv", pat)))
  }
  
  # --- Final Return: List containing both files and plots ---
  return(list(
    files = unique(written),
    plots = plots
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
    feature_id   = feats,
    best_pattern = best_pat,
    best_corr    = best_cor,
    stringsAsFactors = FALSE
  )
}