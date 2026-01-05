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
  
  eff        <- cfg$effects
  sample_col <- eff$samples
  color_col  <- eff$color
  shape_col  <- eff$shape
  if (is.null(shape_col)) shape_col <- NULL  # avoid dependency on %||%
  
  sample_ids <- colnames(expr_mat)
  meta_sub   <- align_meta_to_matrix(sample_ids, meta, sample_col)
  
  # Ensure color_col exists (fail-fast with clear error)
  if (is.null(color_col) || !color_col %in% colnames(meta_sub)) {
    stop("qc_pca_scatter(): color column '", color_col, "' not found in aligned metadata.")
  }
  
  # PCA via shared helper
  res      <- compute_pca_scores(expr_mat, pcs = pcs, center = TRUE, scale = FALSE)
  scores   <- res$scores
  var_expl <- res$var_expl
  
  # Ensure expected PC columns exist (should always be true now)
  needed_pc_cols <- paste0("PC", pcs)
  if (!all(needed_pc_cols %in% colnames(scores))) {
    stop("qc_pca_scatter(): PCA scores missing expected columns: ",
         paste(setdiff(needed_pc_cols, colnames(scores)), collapse = ", "))
  }
  
  pc_labels <- sprintf("PC%d: %.2f%% of variance", seq_along(var_expl), 100 * var_expl)
  
  pc_x <- pcs[1]
  pc_y <- pcs[2]
  
  # Attach metadata
  scores[[color_col]] <- meta_sub[[color_col]]
  if (!is.null(shape_col) && shape_col %in% colnames(meta_sub)) {
    scores[[shape_col]] <- meta_sub[[shape_col]]
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
  meta     <- as.data.frame(meta)
  
  eff        <- cfg$effects
  sample_col <- eff$samples
  color_col  <- eff$color
  shape_col  <- eff$shape
  if (is.null(shape_col)) shape_col <- NULL
  
  if (!sample_col %in% colnames(meta)) {
    stop("Sample column '", sample_col, "' not found in metadata.")
  }
  
  sample_ids <- colnames(expr_mat)
  meta_sub   <- align_meta_to_matrix(sample_ids, meta, sample_col)
  
  # PCA via shared helper
  res      <- compute_pca_scores(expr_mat, pcs = 1:3, center = TRUE, scale = FALSE)
  scores   <- res$scores
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
    x    = ~PC1,
    y    = ~PC2,
    z    = ~PC3,
    type = "scatter3d",
    mode = "markers",
    color  = if (!is.null(color_col)) scores[[color_col]] else NULL,
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
    htmlwidgets::saveWidget(widget = plt, file = out_file, selfcontained = TRUE)
  }
  
  plt
}
