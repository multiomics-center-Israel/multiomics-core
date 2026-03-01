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
    pca_obj = pca,
    scores = scores,
    var_expl = var_expl
  )
}
#' PCA scatter plot
#'
#' @param expr_mat numeric matrix (features x samples), normalized/imputed.
#' @param meta     metadata table.
#' @param cfg      config$modes$proteomics (effects$samples/color/shape).
#' @param pcs      length-2 vector of PCs to plot, e.g. c(1, 2) or c(1, 3).
#' @param out_file optional path to save plot.
qc_pca_scatter <- function(expr_mat, meta, cfg, pcs = c(1, 2), out_file = NULL) {
  expr_mat <- as.matrix(expr_mat)

  # Basic input validation
  if (is.null(colnames(expr_mat)) || length(colnames(expr_mat)) == 0) {
    stop("qc_pca_scatter(): expr_mat must have sample column names.")
  }
  pcs <- as.integer(pcs)
  if (length(pcs) != 2 || anyNA(pcs) || any(pcs < 1)) {
    stop("qc_pca_scatter(): pcs must be a length-2 vector of positive integers, e.g. c(1,2).")
  }

  eff <- cfg$effects
  sample_col <- eff$samples
  color_col <- eff$color
  shape_col <- eff$shape
  if (is.null(shape_col)) shape_col <- NULL # avoid dependency on %||%

  sample_ids <- colnames(expr_mat)
  meta_sub <- align_meta_to_matrix(sample_ids, meta, sample_col)

  # Ensure color_col exists (fail-fast with clear error)
  if (is.null(color_col) || !color_col %in% colnames(meta_sub)) {
    stop("qc_pca_scatter(): color column '", color_col, "' not found in aligned metadata.")
  }

  # PCA via shared helper
  res <- compute_pca_scores(expr_mat, pcs = pcs, center = TRUE, scale = FALSE)
  scores <- res$scores
  var_expl <- res$var_expl

  # Ensure expected PC columns exist (should always be true now)
  needed_pc_cols <- paste0("PC", pcs)
  if (!all(needed_pc_cols %in% colnames(scores))) {
    stop(
      "qc_pca_scatter(): PCA scores missing expected columns: ",
      paste(setdiff(needed_pc_cols, colnames(scores)), collapse = ", ")
    )
  }

  pc_labels <- sprintf("PC%d: %.2f%% of variance", seq_along(var_expl), 100 * var_expl)

  pc_x <- pcs[1]
  pc_y <- pcs[2]

  # Attach ALL metadata columns to scores (not just color/shape)
  # This ensures Shiny has access to all metadata via pca_scores
  for (col in colnames(meta_sub)) {
    if (!col %in% colnames(scores)) {
      scores[[col]] <- meta_sub[[col]]
    }
  }

  p <- plot_pca_scatter(
    scores    = scores,
    color_col = color_col,
    shape_col = shape_col,
    pc_x      = pc_x,
    pc_y      = pc_y,
    pc_labels = pc_labels
  )

  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file, plot = p, width = 6, height = 5)
  }

  # Attach PCA results as attributes (backward compatible - plot is still returned)
  attr(p, "pca_result") <- res$pca_obj
  attr(p, "scores") <- scores # scores now include all metadata
  attr(p, "var_expl") <- var_expl

  p
}


#' 3D PCA plot (PC1 vs PC2 vs PC3)
#'
#' @param expr_mat numeric matrix (features x samples), normalized/imputed.
#' @param meta     metadata table.
#' @param cfg      config$modes$proteomics (effects$samples/color/shape).
#' @param out_file optional HTML path to save widget (if NULL, does not save).
#'
#' @return plotly object (interactive 3D PCA).
qc_pca_3d <- function(expr_mat, meta, cfg, out_file = NULL) {
  expr_mat <- as.matrix(expr_mat)
  meta <- as.data.frame(meta)

  eff <- cfg$effects
  sample_col <- eff$samples
  color_col <- eff$color
  shape_col <- eff$shape
  label_col <- eff$label
  if (is.null(shape_col)) shape_col <- NULL

  if (!sample_col %in% colnames(meta)) {
    stop("Sample column '", sample_col, "' not found in metadata.")
  }

  sample_ids <- colnames(expr_mat)
  meta_sub <- align_meta_to_matrix(sample_ids, meta, sample_col)

  # PCA via shared helper
  res <- compute_pca_scores(expr_mat, pcs = 1:3, center = TRUE, scale = FALSE)
  scores <- res$scores
  var_expl <- res$var_expl

  pc_labels <- sprintf("PC%d (%.1f%%)", seq_along(var_expl), 100 * var_expl)

  # Attach metadata (order is aligned)
  scores[[color_col]] <- meta_sub[[color_col]]
  if (!is.null(shape_col) && shape_col %in% colnames(meta_sub)) {
    scores[[shape_col]] <- meta_sub[[shape_col]]
  }

  # Attach label column if configured
  if (!is.null(label_col) && label_col %in% colnames(meta_sub)) {
    scores[[label_col]] <- meta_sub[[label_col]]
  }

  # Hover text
  hover_text <- scores$sample
  if (!is.null(label_col) && label_col %in% colnames(scores)) {
    hover_text <- paste0(hover_text, "<br>", label_col, ": ", scores[[label_col]])
  }
  if (!is.null(color_col)) {
    hover_text <- paste0(hover_text, "<br>", color_col, ": ", scores[[color_col]])
  }
  if (!is.null(shape_col) && shape_col %in% colnames(scores)) {
    hover_text <- paste0(hover_text, "<br>", shape_col, ": ", scores[[shape_col]])
  }

  plt <- plotly::plot_ly(
    data = scores,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    type = "scatter3d",
    mode = "markers",
    color = if (!is.null(color_col)) scores[[color_col]] else NULL,
    symbol = if (!is.null(shape_col) && shape_col %in% colnames(scores)) scores[[shape_col]] else NULL,
    text = hover_text,
    hoverinfo = "text"
  )|>
    plotly::layout(
      scene = list(
        xaxis = list(title = pc_labels[1]),
        yaxis = list(title = pc_labels[2]),
        zaxis = list(title = pc_labels[3])
      ),
      title = "3D PCA: PC1 vs PC2 vs PC3"
    )

  if (!is.null(out_file)) {
    # Check if pandoc is available for self-contained HTML
    has_pandoc <- nzchar(Sys.which("pandoc"))
    if (!has_pandoc) {
      warning("Pandoc not found: Saving 3D PCA widget as non-self-contained (creates matching '_files' directory). Install Pandoc to enable self-contained HTML.")
    }
    htmlwidgets::saveWidget(widget = plt, file = out_file, selfcontained = has_pandoc)
  }

  plt
}

#------- HELPERS -------

#' Prepare data for QC plotting
#'
#' Handles common tasks: matrix conversion, metadata matching, and annotation creation.
#' Includes strict validation to fail fast on data inconsistencies (duplicates, missing samples).
#'
#' @param expr Numeric matrix or data.frame.
#' @param meta Metadata data.frame.
#' @param cfg Config object containing effects$samples and effects$color.
#' @return A list containing aligned expr, meta, annot, and column names.
prepare_qc_data <- function(expr, meta, cfg) {
  # 1. Ensure Expression is a Matrix with Names
  expr <- as.matrix(expr)
  sample_ids <- colnames(expr)

  # FIX: Handle empty matrix or NULL colnames
  if (is.null(sample_ids) || length(sample_ids) == 0) {
    stop("Expression matrix must have non-empty column names (sample IDs).")
  }

  if (anyDuplicated(sample_ids)) {
    stop("Expression matrix contains duplicate column names/sample IDs.")
  }

  # 2. Extract Config & Validate Keys
  eff <- cfg$effects
  if (is.null(eff$samples) || is.null(eff$color)) {
    stop("Config must contain cfg$effects$samples and cfg$effects$color")
  }
  sample_col <- eff$samples
  color_col <- eff$color

  # 3. Ensure Metadata is a base data.frame (safe against tibbles)
  meta <- as.data.frame(meta)

  # 4. Validate Columns Existence in Metadata
  if (!sample_col %in% colnames(meta)) {
    stop(sprintf("Sample column '%s' not found in metadata.", sample_col))
  }
  if (!color_col %in% colnames(meta)) {
    stop(sprintf("Color/Condition column '%s' not found in metadata.", color_col))
  }

  # FIX: Critical check for duplicates in metadata ID column
  # match() returns the first hit, so duplicates would be silently ignored without this.
  if (anyDuplicated(meta[[sample_col]])) {
    stop(sprintf("Metadata column '%s' contains duplicate sample IDs.", sample_col))
  }

  # 5. Match Metadata to Expression Columns
  # This guarantees that meta_sub rows are 1:1 with expr columns
  meta_sub <- meta[match(sample_ids, meta[[sample_col]]), , drop = FALSE]

  # 6. Verify Integrity (Missing samples)
  # FIX: Detailed error message with count + examples
  if (any(is.na(meta_sub[[sample_col]]))) {
    missing_mask <- is.na(meta_sub[[sample_col]])
    missing_samples <- sample_ids[missing_mask]

    stop(sprintf(
      "Found %d samples in expression matrix missing from metadata column '%s'. Examples: %s",
      length(missing_samples),
      sample_col,
      paste(head(missing_samples, 3), collapse = ", ")
    ))
  }

  # FIX: Extra sanity check for 1:1 alignment
  if (nrow(meta_sub) != ncol(expr)) {
    stop("Internal error: metadata alignment failed (meta_sub rows != expr columns).")
  }

  # 7. Create Generic Annotation (for pheatmap)
  annot <- data.frame(
    Condition = meta_sub[[color_col]],
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )

  list(
    expr = expr,
    meta = meta_sub, # Perfectly aligned with expr columns
    annot = annot,
    sample_col = sample_col,
    color_col = color_col,
    sample_ids = sample_ids
  )
}
#' Convert Expression Matrix to Long Format
#'
#' Used for ggplot2. Optimized to use 'rep' since data is already aligned.
to_long_format <- function(prep_data) {
  n_features <- nrow(prep_data$expr)

  # Vectorized creation of long vectors
  norm_expr_long <- data.frame(
    sample = rep(prep_data$sample_ids, each = n_features),
    value = as.vector(prep_data$expr),
    stringsAsFactors = FALSE
  )

  # OPTIMIZATION: Use rep() instead of match()
  # Since prep_data$meta is guaranteed to be sorted by sample_ids (from prepare_qc_data),
  # we can simply repeat each condition n_features times.
  cond_values <- prep_data$meta[[prep_data$color_col]]
  norm_expr_long[[prep_data$color_col]] <- rep(cond_values, each = n_features)

  # Filter non-finite values (NA/NaN/Inf) to keep plots clean
  norm_expr_long[is.finite(norm_expr_long$value), , drop = FALSE]
}


#------- MAIN PLOTTING FUNCTIONS -------


#' Boxplots of normalized expression per sample
#' @param expr_norm Expression matrix (genes x samples)
#' @param meta Sample metadata data.frame
#' @param cfg Mode config with effects$color, effects$samples
#' @param out_file Optional output file path
#' @param title Optional custom title (default: "Normalized expression boxplots")
#' @param y_label Optional y-axis label (default: "log2(normalized intensity)")
norm_boxplot <- function(expr_norm, meta, cfg, out_file = NULL, title = NULL,
                         y_label = "log2(normalized intensity)") {
  d <- prepare_qc_data(expr_norm, meta, cfg)
  norm_expr_long <- to_long_format(d)

  if (is.null(title)) title <- "Normalized expression boxplots"

  plot_title <- title %||% "Normalized expression boxplots"

  p <- ggplot2::ggplot(
    norm_expr_long,
    ggplot2::aes(x = sample, y = value, colour = .data[[d$color_col]])
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.4) +
    ggplot2::labs(
      title  = plot_title,
      x      = "Sample",
      y      = y_label,
      colour = d$color_col
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file, plot = p, width = 8, height = 5)
  }

  invisible(p)
}


#' RLE (Relative Log Expression) boxplot per sample
#'
#' Computes per-feature deviation from the row median and draws a boxplot
#' per sample.  Well-normalized data centers around zero with tight IQRs.
#' @param expr_norm Expression matrix (features x samples, log2 scale)
#' @param meta Sample metadata data.frame
#' @param cfg Mode config with effects$color, effects$samples
#' @param out_file Optional output file path
qc_rle_boxplot <- function(expr_norm, meta, cfg, out_file = NULL) {
  d <- prepare_qc_data(expr_norm, meta, cfg)

  row_medians <- apply(d$expr, 1, stats::median, na.rm = TRUE)
  rle_mat <- d$expr - row_medians

  n_features <- nrow(rle_mat)
  rle_long <- data.frame(
    sample = rep(d$sample_ids, each = n_features),
    value  = as.vector(rle_mat),
    stringsAsFactors = FALSE
  )
  cond_values <- d$meta[[d$color_col]]
  rle_long[[d$color_col]] <- rep(cond_values, each = n_features)
  rle_long <- rle_long[is.finite(rle_long$value), , drop = FALSE]

  p <- ggplot2::ggplot(
    rle_long,
    ggplot2::aes(x = sample, y = value, colour = .data[[d$color_col]])
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.4) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::labs(
      title  = "Relative Log Expression (RLE)",
      x      = "Sample",
      y      = "Deviation from feature median",
      colour = d$color_col
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file, plot = p, width = 8, height = 5)
  }

  invisible(p)
}


#' Histogram summary of normalized expression by condition
norm_histogram_summary <- function(expr_norm, meta, cfg, out_file = NULL) {
  d <- prepare_qc_data(expr_norm, meta, cfg)
  norm_expr_long <- to_long_format(d)

  p <- ggplot2::ggplot(
    norm_expr_long,
    ggplot2::aes(x = value, fill = .data[[d$color_col]])
  ) +
    ggplot2::geom_histogram(alpha = 0.6, bins = 60, position = "identity") +
    ggplot2::facet_wrap(as.formula(paste("~", d$color_col)), nrow = 1, scales = "free_y") +
    ggplot2::labs(
      title = "Normalized expression histograms by condition",
      x     = "log2(normalized intensity)",
      y     = "Frequency",
      fill  = d$color_col
    ) +
    ggplot2::theme_minimal()

  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file, plot = p, width = 8, height = 4)
  }

  invisible(p)
}

wrap_qc_heatmap <- function(expr_mat, meta, cfg, stage, out_file = NULL, cluster_cols = TRUE) {
  # 1. Prepare Data
  d <- prepare_qc_data(expr_mat, meta, cfg)

  # Rename column to match config
  col_name <- d$color_col
  annot <- d$annot
  if ("Condition" %in% names(annot)) {
    colnames(annot)[which(names(annot) == "Condition")] <- col_name
  }

  ph <- plot_heatmap_core(
    expr_mat = d$expr,
    annotation_col = annot,
    title = paste0("QC: Sample ", stage, " Expression"),
    max_rows = 2000,
    cluster_rows = TRUE,
    cluster_cols = cluster_cols
  )

  if (!is.null(out_file)) {
    save_heatmap_to_file(ph, out_file)
  }

  return(ph)
}

#' Sample–sample correlation heatmap (proteomics QC)
#' @param annot_cols Character vector of metadata columns for annotations.
#'   If NULL (default), uses cfg$effects$color only (backward compatible).
#' @param show_labels Show sample names on axes (default FALSE to prevent overload)
#' @param cluster_samples Enable hierarchical clustering (default TRUE)
#' @param adjust_scale Adjust color scale for tight correlation ranges (default TRUE)
qc_sample_correlation_heatmap <- function(expr_mat, meta, cfg, out_file,
                                          annot_cols = NULL,
                                          method = "pearson",
                                          fontsize = 10,
                                          show_labels = TRUE,
                                          cluster_samples = TRUE,
                                          adjust_scale = TRUE) {
  d <- prepare_qc_data(expr_mat, meta, cfg)

  # Build annotations: use multi-column if provided, otherwise default
  if (!is.null(annot_cols)) {
    annot <- d$meta[d$sample_ids, annot_cols, drop = FALSE]
    rownames(annot) <- d$sample_ids
  } else {
    annot <- d$annot  # Default single-column annotation
  }

  ph <- plot_sample_correlation_heatmap(
    expr_mat = d$expr,
    method = method,
    annotation_col = annot,
    fontsize = fontsize,
    show_labels = show_labels,
    cluster_rows = cluster_samples,
    cluster_cols = cluster_samples,
    adjust_scale = adjust_scale
  )

  # Issue 3 & 6 FIX: Larger canvas for better legend placement and readability
  save_heatmap_to_file(ph, out_file, width = 2400, height = 2000, res = 150)
  invisible(ph)
}


#' Sample–sample distance heatmap (QC)
#' @param annot_cols Character vector of metadata columns for annotations.
#'   If NULL (default), uses cfg$effects$color only (backward compatible).
#' @param show_labels Show sample names on axes (default FALSE to prevent overload)
#' @param cluster_samples Enable hierarchical clustering (default TRUE)
qc_sample_distance_heatmap <- function(expr_mat, meta, cfg, out_file,
                                       annot_cols = NULL,
                                       with_na = FALSE,
                                       fontsize = 10,
                                       show_labels = TRUE,
                                       cluster_samples = TRUE) {
  # FIX: Ensure matrix conversion happens BEFORE complete.cases logic
  # This prevents data.frame type coercion issues
  expr_mat <- as.matrix(expr_mat)

  if (!with_na) {
    keep <- stats::complete.cases(expr_mat)
    expr_mat <- expr_mat[keep, , drop = FALSE]
  }

  # Validate & Prepare (after filtering NAs)
  d <- prepare_qc_data(expr_mat, meta, cfg)

  # Build annotations: use multi-column if provided, otherwise default
  if (!is.null(annot_cols)) {
    annot <- d$meta[d$sample_ids, annot_cols, drop = FALSE]
    rownames(annot) <- d$sample_ids
  } else {
    annot <- d$annot  # Default single-column annotation
  }

  ph <- plot_sample_distance_heatmap(
    expr_mat = d$expr,
    annotation_col = annot,
    fontsize = fontsize,
    show_labels = show_labels,
    cluster_rows = cluster_samples,
    cluster_cols = cluster_samples
  )

  # Issue 3 & 6 FIX: Larger canvas for better legend placement and readability
  save_heatmap_to_file(ph, out_file, width = 2400, height = 2000, res = 150)
  invisible(ph)
}


#' Density overlay of normalized expression
qc_omic_density <- function(expr_mat, meta, cfg, out_file = NULL,
                                  alpha = 0.3, show_legend = TRUE,
                                  title = "Density plot of normalized intensities") {
  # Keeping prepare_qc_data here intentionally.
  d <- prepare_qc_data(expr_mat, meta, cfg)
  norm_expr_long <- to_long_format(d)

  # Task 2: Density line only (no fill) - remove fill aesthetic, set alpha=1
  p <- ggplot2::ggplot(norm_expr_long, ggplot2::aes(x = value, colour = .data[["sample"]])) +
    ggplot2::geom_density(alpha = 1) +
    ggplot2::labs(title = title, x = "Intensity", colour = d$color_col) +
    ggplot2::theme_minimal()

  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }



  if (!is.null(out_file)) {
    if (!grepl("\\.png$", out_file, ignore.case = TRUE)) {
      out_file <- paste0(out_file, ".png")
    }
    ggplot2::ggsave(out_file, p, width = 10, height = 6, dpi = 150)
  }

  invisible(p)
}

# (Existing function, kept for completeness)
#' Imputation QC Summary (Wrapper)
#'
#' @param imputed Numeric matrix (imputed expression)
#' @param imputed_flag Logical matrix (TRUE if imputed)
#' @param cfg Config object (optional, for consistency but unused logic-wise)
#' @param out_file Output file path
#' @return ggplot object
qc_imputation_summary <- function(imputed, imputed_flag, cfg = NULL, out_file = NULL) {
  imp_width <- cfg$imputation$width %||% NULL
  imp_shift <- cfg$imputation$downshift %||% NULL
  p <- plot_imputation_summary(imputed, imputed_flag, width = imp_width, downshift = imp_shift)

  if (!is.null(out_file)) {
    # Ensure correct extension if missing
    if (!grepl("\\.png$", out_file, ignore.case = TRUE)) {
      out_file <- paste0(out_file, ".png")
    }
    ggplot2::ggsave(out_file, plot = p, width = 12, height = 8, dpi = 150)
  }

  invisible(p)
}

#' QC wrapper: missingness heatmap
#'
#' Calls \code{plot_missingness_heatmap()} and saves the result to a PNG file.
#' Follows the established \code{qc_*} pattern: handles I/O and calls the pure
#' \code{plot_*} function internally.
#'
#' @param mat Numeric matrix (features x samples).
#' @param miss_stats data.frame returned by \code{classify_missingness()}, with
#'   columns \code{feature_id} and \code{mnar_class} ("MNAR", "MAR",
#'   "all_missing", "all_observed"). Alternatively, a named logical vector where
#'   TRUE = MNAR (as returned by \code{classify_mnar_features()}).
#' @param out_file Output PNG file path.
#'
#' @return Invisibly returns the pheatmap object.
#'
#' @seealso plot_missingness_heatmap
#'
qc_missingness_heatmap <- function(mat, miss_stats, out_file) {
  # Accept either a data.frame (from classify_missingness) or a logical vector
  if (is.data.frame(miss_stats)) {
    mnar_class <- stats::setNames(miss_stats$mnar_class, miss_stats$feature_id)
  } else {
    # Named logical vector: TRUE = MNAR
    mnar_class <- miss_stats
  }

  # Restrict to features present in mat
  common <- intersect(names(mnar_class), rownames(mat))
  mnar_class <- mnar_class[common]

  ph <- plot_missingness_heatmap(mat, mnar_class)
  save_heatmap_to_file(ph, out_file, width = 2400, height = 1800)
  invisible(ph)
}

#' Write per-sample imputation histograms (Batch writer)
#' Refactored: Uses ggsave and returns objects
write_imputation_histograms_per_sample <- function(expr_mat, imputed_flag, out_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  sample_ids <- colnames(expr_mat)
  files <- character(0)
  plots <- list()
  for (s in sample_ids) {
    p <- plot_imputation_histogram_one_sample(expr_mat, imputed_flag, s)
    f <- file.path(out_dir, paste0("X", s, ".png"))
    ggplot2::ggsave(filename = f, plot = p, width = 3.33, height = 6, dpi = 150)
    files <- c(files, f)
    plots[[paste0("sample_", s)]] <- p
  }

  return(list(files = files, plots = plots))
}
