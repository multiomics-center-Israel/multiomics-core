
# plot_core.R
#
# Core plotting functions used across all omics modalities.
#
# Rules:
# - Functions in this file must be named plot_*.
# - plot_* functions are PURE:
#     * no config access
#     * no file I/O (no ggsave / png / with_png)
#     * no side effects
# - They return plot objects (ggplot / pheatmap / etc.).
#
# QC wrappers (qc_*) are responsible for:
# - metadata alignment
# - reading cfg$effects
# - saving figures to disk
#



#' Density overlay per sample (no saving, just ggplot)
#'
#' @param expr_mat numeric matrix/data.frame (features x samples)
#' @param sample_ids character vector = colnames(expr_mat) (optional)
#' @return ggplot object
plot_density_overlay <- function(expr_mat,
                                 sample_ids = colnames(expr_mat),
                                 alpha      = 0.3,
                                 title      = "Density plot of normalized intensities") {
  stopifnot(is.matrix(expr_mat) || is.data.frame(expr_mat))
  expr_mat <- as.data.frame(expr_mat)
  
  df_long <- expr_mat |>
    tibble::rownames_to_column("feature") |>
    tidyr::pivot_longer(
      cols = -feature,
      names_to = "SampleID",
      values_to = "value"
    )
  
  df_long <- df_long[is.finite(df_long$value), , drop = FALSE]
  
  library(ggplot2)
  ggplot(df_long, aes(x = value, color = SampleID)) +
    geom_density(alpha = alpha, linewidth = 0.7) +
    labs(
      title = title,
      x     = "log2 intensity",
      y     = "Density"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      legend.title    = element_blank()
    )
}

#' Sample–sample distance heatmap (no saving)
#'
#' @param expr_mat numeric matrix (features x samples)
#' @param dist_method distance metric
#' @return invisibly returns pheatmap object
plot_sample_distance_heatmap <- function(expr_mat,
                                         dist_method = "euclidean",
                                         annotation_col = NULL,
                                         main = NULL,
                                         colors = NULL,
                                         out_file = NULL,
                                         fontsize = 12) {
  
  expr_mat <- as.matrix(expr_mat)
  sampleDists <- stats::dist(t(expr_mat), method = dist_method)
  mat <- as.matrix(sampleDists)
  
  if (is.null(colors)) {
    colors <- grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(9, "Blues"))
    )(255)
  }
  if (is.null(main)) main <- sprintf("Sample distance heatmap (%s)", dist_method)
  
  ph <- pheatmap::pheatmap(
    mat,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    annotation_col = annotation_col,
    main = main,
    col = colors,
    fontsize_row = fontsize,
    fontsize_col = fontsize
  )
  
  if (!is.null(out_file)) {
    save_pheatmap_png(ph, out_file, width = 1600, height = 1200, res = 150)
  }
  
  invisible(ph)
}



#' Expression heatmap (features x samples, no saving)
#'
#' @param expr_mat numeric matrix (features x samples)
#' @param max_rows optional subsampling of rows
#' @return pheatmap object
plot_expr_heatmap <- function(expr_mat,
                              annotation_col = NULL,
                              max_rows = 2000,
                              main = NULL,
                              scale_rows = TRUE,
                              cluster_rows = TRUE,
                              cluster_cols = TRUE) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Need pheatmap")
  
  expr_mat <- as.matrix(expr_mat)
  
  if (!is.null(max_rows) && nrow(expr_mat) > max_rows) {
    set.seed(1)
    expr_mat <- expr_mat[sample(seq_len(nrow(expr_mat)), max_rows), , drop = FALSE]
  }
  
  if (is.null(main)) main <- sprintf("Expression heatmap (%d features)", nrow(expr_mat))
  
  pheatmap::pheatmap(
    expr_mat,
    scale = if (scale_rows) "row" else "none",
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = FALSE,
    annotation_col = annotation_col,
    main = main
  )
}

#' Sample–sample correlation heatmap
#'
#' Computes pairwise correlations between samples and visualizes them as a heatmap.
#'
#' This plot reflects **biological similarity** between samples and is typically
#' preferred over distance-based heatmaps for proteomics QC, as it is robust to
#' global intensity shifts and scaling effects.
#'
#' Correlations are computed using pairwise complete observations to tolerate
#' missing values.
#'
#' @param expr_mat Numeric matrix (features x samples), typically log2 normalized.
#' @param method Correlation method; one of "pearson", "spearman", or "kendall".
#' @param annotation_col Optional data frame of sample annotations
#'        (rows named by sample IDs).
#' @param main Optional plot title.
#' @param colors Optional color palette for the heatmap.
#' @param fontsize Numeric font size for row/column labels.
#'
#' @return pheatmap object.
#'
#' @seealso qc_sample_correlation_heatmap
#'
plot_sample_correlation_heatmap <- function(expr_mat,
                                            method = "pearson",
                                            annotation_col = NULL,
                                            main = NULL,
                                            colors = NULL,
                                            fontsize = 12) {
  
  expr_mat <- as.matrix(expr_mat)
  
  cor_mat <- stats::cor(
    expr_mat,
    use = "pairwise.complete.obs",
    method = method
  )
  
  if (is.null(colors)) {
    colors <- grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(9, "Blues"))
    )(255)
  }
  if (is.null(main)) {
    main <- sprintf("Sample correlation heatmap (%s)", method)
  }
  
  pheatmap::pheatmap(
    cor_mat,
    annotation_col = annotation_col,
    annotation_row = annotation_col,
    main = main,
    col = colors,
    fontsize_row = fontsize,
    fontsize_col = fontsize
  )
}


