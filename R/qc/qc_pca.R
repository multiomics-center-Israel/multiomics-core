#' PCA scatter plot
#'
#' @param expr_mat numeric matrix (features x samples), normalized/imputed.
#' @param meta     metadata table.
#' @param cfg      config$modes$proteomics (effects$samples/color/shape).
#' @param pcs      length-2 vector of PCs to plot, e.g. c(1, 2) or c(1, 3).
#' @param out_file optional path to save plot.
qc_pca_scatter <- function(expr_mat, meta, cfg, pcs=c(1, 2), out_file = NULL) {
  # NOTE:
  # Input should be a processed expression matrix with no missing values
  # (e.g., imputed proteomics or normalized/transformed RNA).
  # PCA uses centering only (scale = FALSE) for consistency with legacy analyses.
  
  expr_mat <- as.matrix(expr_mat)
  
  eff         <- cfg$effects
  sample_col  <- eff$samples
  color_col   <- eff$color
  shape_col   <- eff$shape %||% NULL
  
  sample_ids <- colnames(expr_mat)
  meta_sub <- align_meta_to_matrix(sample_ids, meta, sample_col)

  # PCA: samples as rows
  res <- compute_pca_scores(expr_mat, pcs = pcs, center = TRUE, scale = FALSE)
  
  scores <- res$scores
  var_expl <- res$var_expl
  
  pc_labels <- sprintf("PC%d: %.2f%% of variance", seq_along(var_expl), 100 * var_expl)
  
  pc_x <- pcs[1]
  pc_y <- pcs[2]
  
  scores[[color_col]] <- meta_sub[[color_col]]
  
  if (!is.null(shape_col) && shape_col %in% colnames(meta_sub)) {
    scores[[shape_col]] <- meta_sub[[shape_col]]}
  
  p <- build_pca_scatter_plot(
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
  
  invisible(p)
}


#' 3D PCA plot (PC1 vs PC2 vs PC3)
#'
#' @param expr_mat numeric matrix (features x samples), normalized/imputed.
#' @param meta     metadata table.
#' @param cfg      config$modes$proteomics (effects$samples/color/shape).
#'
#' @return plotly object (interactive 3D PCA).
#' 3D PCA plot (PC1 vs PC2 vs PC3)
qc_pca_3d <- function(expr_mat,
                      meta,
                      cfg,
                      out_file = NULL) {
  
  # 2. Basic conversions to prevent Tibble/Matrix errors
  expr_mat <- as.matrix(expr_mat)
  
  # Convert to standard data.frame to avoid warnings about setting row names on a tibble
  meta <- as.data.frame(meta) 
  
  eff         <- cfg$effects
  sample_col  <- eff$samples
  color_col   <- eff$color
  shape_col   <- eff$shape %||% NULL # Assumes %||% is defined, otherwise use standard if-null logic
  if(is.null(eff$shape)) shape_col <- NULL 
  
  # Ensure the sample column exists in metadata
  if (!sample_col %in% colnames(meta)) {
    stop("Sample column '", sample_col, "' not found in metadata.")
  }
  
  sample_ids <- colnames(expr_mat)
  
  # Align metadata order with matrix column order
  # (Important: This ensures consistency with the PCA results later)
  meta_sub <- align_meta_to_matrix(sample_ids, meta, sample_col)
  
  # PCA Calculation
  pca <- stats::prcomp(t(expr_mat), center = TRUE, scale. = FALSE)
  
  
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  pc_labels <- sprintf(
    "PC%d (%.1f%%)",
    seq_along(var_expl),
    100 * var_expl
  )
  
  # Prepare results table (scores)
  scores <- as.data.frame(pca$x[, 1:3, drop = FALSE])
  colnames(scores) <- c("PC1", "PC2", "PC3")
  scores$sample <- rownames(scores)
  
  # Add metadata to Scores table
  # Since meta_sub is sorted by sample_ids, and PCA was run on sample_ids, the order is identical.
  scores[[color_col]] <- meta_sub[[color_col]]
  
  if (!is.null(shape_col) && shape_col %in% colnames(meta_sub)) {
    scores[[shape_col]] <- meta_sub[[shape_col]]
  }
  
  # Build Hover Text
  hover_text <- scores$sample
  if (!is.null(color_col)) {
    hover_text <- paste0(
      hover_text,
      "<br>", color_col, ": ", scores[[color_col]]
    )
  }
  if (!is.null(shape_col) && shape_col %in% colnames(scores)) {
    hover_text <- paste0(
      hover_text,
      "<br>", shape_col, ": ", scores[[shape_col]]
    )
  }
  
  # Create plot (using Native Pipe |> instead of %>%)
  plt <- plotly::plot_ly(
    data = scores,
    x    = ~PC1,
    y    = ~PC2,
    z    = ~PC3,
    type = "scatter3d",
    mode = "markers",
    color = if (!is.null(color_col)) scores[[color_col]] else NULL,
    symbol = if (!is.null(shape_col) && shape_col %in% colnames(scores)) scores[[shape_col]] else NULL,
    text  = hover_text,
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
  htmlwidgets::saveWidget(
    widget = plt,
    file = out_file,
    selfcontained = TRUE
  )
  
  return(plt)
}