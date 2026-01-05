#' Save a pheatmap object to file
save_heatmap_to_file <- function(pheatmap_obj, out_file,
                                 width = 1600, height = 1200, res = 150) {
  if (is.null(out_file)) return(invisible(NULL))
  
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  if (!grepl("\\.png$", out_file, ignore.case = TRUE)) out_file <- paste0(out_file, ".png")
  
  grDevices::png(filename = out_file, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  grid::grid.newpage()
  grid::grid.draw(pheatmap_obj$gtable)
  
  out_file
}
#' Build a long-format table for imputation QC
#' Used by plot_imputation_summary
build_imputation_long_df <- function(expr_mat, imputed_flag) {
  expr_mat <- as.matrix(expr_mat)
  if (is.null(imputed_flag)) stop("imputed_flag is NULL")
  imputed_flag <- as.matrix(imputed_flag)
  
  stopifnot(identical(dim(expr_mat), dim(imputed_flag)))
  
  df_raw <- data.frame(
    sample     = rep(colnames(expr_mat), each = nrow(expr_mat)),
    value      = as.vector(expr_mat),
    is_imputed = as.vector(imputed_flag),
    stringsAsFactors = FALSE
  )
  
  df_plot <- df_raw[is.finite(df_raw$value), , drop = FALSE]
  list(raw = df_raw, plot = df_plot)
}
#' Standard Blue heatmap palette (DRY helper)
get_heatmap_colors <- function(n = 255) {
  grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(n)
}


