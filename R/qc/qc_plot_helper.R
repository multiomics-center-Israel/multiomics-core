build_pca_scatter_plot <- function(scores, color_col, shape_col = NULL,
                                   pc_x, pc_y, pc_labels) {
  aes_args <- list(
    x = quote(PCx),
    y = quote(PCy),
    colour = substitute(.data[[color_col]])
  )
  
  if (!is.null(shape_col) && shape_col %in% colnames(scores)) {
    aes_args$shape <- substitute(.data[[shape_col]])
  }
  
  ggplot2::ggplot(scores, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title  = sprintf("PCA: PC%d vs PC%d", pc_x, pc_y),
      x      = pc_labels[pc_x],
      y      = pc_labels[pc_y],
      colour = color_col,
      shape  = if (!is.null(shape_col)) shape_col else NULL
    ) +
    ggplot2::theme_minimal()
}
