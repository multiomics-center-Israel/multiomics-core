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
  # This ensures Shiny has access to all metadata via mat2plot
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
    ggplot2::ggsave(out_file, plot = p, width = 5, height = 5)
  }

  # Attach PCA results as attributes (backward compatible - plot is still returned)
  attr(p, "pca_result") <- res$pca_obj
  attr(p, "scores") <- scores # scores now include all metadata
  attr(p, "var_expl") <- var_expl

  invisible(p)
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

  # Hover text
  hover_text <- scores$sample
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
  ) |>
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
#' Includes strict validation to fail fast on data inconsistencies.
#'
#' @param expr Numeric matrix or data.frame.
#' @param meta Metadata data.frame.
#' @param cfg Config object containing effects$samples and effects$color.
#' @return A list containing aligned expr, meta, annot, and column names.
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
  df_long <- data.frame(
    sample = rep(prep_data$sample_ids, each = n_features),
    value = as.vector(prep_data$expr),
    stringsAsFactors = FALSE
  )

  # OPTIMIZATION: Use rep() instead of match()
  # Since prep_data$meta is guaranteed to be sorted by sample_ids (from prepare_qc_data),
  # we can simply repeat each condition n_features times.
  cond_values <- prep_data$meta[[prep_data$color_col]]
  df_long[[prep_data$color_col]] <- rep(cond_values, each = n_features)

  # Filter non-finite values (NA/NaN/Inf) to keep plots clean
  df_long[is.finite(df_long$value), , drop = FALSE]
}


#------- MAIN PLOTTING FUNCTIONS -------


#' Boxplots of normalized expression per sample
norm_boxplot <- function(expr_norm, meta, cfg, out_file = NULL) {
  d <- prepare_qc_data(expr_norm, meta, cfg)
  df_long <- to_long_format(d)

  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = sample, y = value, colour = .data[[d$color_col]])
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.4) +
    ggplot2::labs(
      title  = "Normalized expression boxplots",
      x      = "Sample",
      y      = "log2(normalized intensity)",
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
  df_long <- to_long_format(d)

  p <- ggplot2::ggplot(
    df_long,
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

wrap_qc_heatmap <- function(expr_mat, meta, cfg, out_file = NULL) {
  # 1. Prepare Data
  # (Assuming prepare_qc_data logic is simple annotation creation)
  annot <- data.frame(
    Condition = meta[[cfg$effects$color]],
    row.names = meta[[cfg$effects$samples]]
  )

  # 2. Plot
  ph <- plot_heatmap_core(
    expr_mat = expr_mat,
    annotation_col = annot,
    title = "QC: Sample Protein Expression",
    max_rows = 2000, # QC usually needs subsampling
    cluster_rows = TRUE,
    cluster_cols = TRUE
  )

  # 3. Save & Return
  if (!is.null(out_file)) save_heatmap_to_file(ph, out_file)
  return(ph)
}

#' Sample–sample correlation heatmap (proteomics QC)
qc_sample_correlation_heatmap <- function(expr_mat, meta, cfg, out_file,
                                          method = "pearson", fontsize = 12) {
  d <- prepare_qc_data(expr_mat, meta, cfg)

  ph <- plot_sample_correlation_heatmap(
    expr_mat = d$expr,
    method = method,
    annotation_col = d$annot,
    fontsize = fontsize
  )

  save_heatmap_to_file(ph, out_file, width = 1600, height = 1200, res = 150)
  invisible(ph)
}


#' Sample–sample distance heatmap (QC)
qc_sample_distance_heatmap <- function(expr_mat, meta, cfg, out_file,
                                       with_na = FALSE, fontsize = 12) {
  # FIX: Ensure matrix conversion happens BEFORE complete.cases logic
  # This prevents data.frame type coercion issues
  expr_mat <- as.matrix(expr_mat)

  if (!with_na) {
    keep <- stats::complete.cases(expr_mat)
    expr_mat <- expr_mat[keep, , drop = FALSE]
  }

  # Validate & Prepare (after filtering NAs)
  d <- prepare_qc_data(expr_mat, meta, cfg)

  ph <- plot_sample_distance_heatmap(
    expr_mat = d$expr,
    annotation_col = d$annot,
    fontsize = fontsize
  )

  save_heatmap_to_file(ph, out_file, width = 1600, height = 1200, res = 150)
  invisible(ph)
}


#' Density overlay of normalized expression
qc_proteomics_density <- function(expr_mat, meta, cfg, out_file = NULL,
                                  alpha = 0.3, show_legend = TRUE,
                                  title = "Density plot of normalized intensities") {
  # Keeping prepare_qc_data here intentionally.
  # Even though density doesn't strictly need metadata to run,
  # we want to ensure the input data is consistent with the pipeline standards.
  d <- prepare_qc_data(expr_mat, meta, cfg)

  p <- plot_density_overlay(
    expr_mat = d$expr,
    title    = title,
    alpha    = alpha
  )

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
imputed_histograms_summary <- function(imputed, imputed_flag, cfg = NULL, out_file = NULL) {
  p <- plot_imputation_summary(imputed, imputed_flag)
  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file, plot = p, width = 12, height = 5, dpi = 150)
  }
  invisible(p)
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
