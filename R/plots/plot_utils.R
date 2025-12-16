#' Save a plotting expression to PNG (for base plots / pheatmap)
#'
#' @param out_file path to png file (if NULL -> draw to screen)
#' @param expr plotting expression (e.g. pheatmap::pheatmap(...))
#' @param width,height in pixels
#' @param res resolution (dpi)
with_png <- function(out_file, expr, width = 1600, height = 1200, res = 150) {
  # NOTE:
  # We must close the PNG device by device ID, not by dev.off(),
  # because some plotting functions (e.g. pheatmap) may change
  # the active graphics device internally.
  
  # If no output file is provided, just evaluate the expression
  # without opening any graphics device
  if (is.null(out_file)) {
    force(expr)
    return(invisible(NULL))
  }
  
  # Ensure the output file has a .png extension
  if (!grepl("\\.png$", out_file, ignore.case = TRUE)) {
    out_file <- paste0(out_file, ".png")
  }
  
  # Create output directory if it does not exist
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  
  # Open a PNG graphics device
  grDevices::png(
    filename = out_file,
    width    = width,
    height   = height,
    res      = res
  )
  
  # Store the device ID of the PNG device we just opened
  # (important: plotting functions may change the active device internally)
  png_dev <- grDevices::dev.cur()
  
  # Ensure that we close exactly this PNG device on exit,
  # even if the active device changes inside `expr`
  on.exit({
    if (png_dev %in% grDevices::dev.list()) {
      try(grDevices::dev.off(which = png_dev), silent = TRUE)
    }
  }, add = TRUE)
  
  # Evaluate the plotting expression
  force(expr)
  
  invisible(out_file)
  
}
#' Build a long-format table for imputation QC (raw + plot-safe)
#'
#' @param expr_mat Numeric matrix (features x samples), typically log2.
#' @param imputed_flag Logical matrix (features x samples), TRUE where value was imputed.
#' @return list(raw=..., plot=...)
build_imputation_long_df <- function(expr_mat, imputed_flag) {
  expr_mat <- as.matrix(expr_mat)
  
  if (is.null(imputed_flag)) stop("imputed_flag is NULL")
  imputed_flag <- as.matrix(imputed_flag)
  
  stopifnot(all(dim(expr_mat) == dim(imputed_flag)))
  stopifnot(!is.null(colnames(expr_mat)))
  
  df_raw <- data.frame(
    sample     = rep(colnames(expr_mat), each = nrow(expr_mat)),
    value      = as.vector(expr_mat),
    is_imputed = as.vector(imputed_flag),
    stringsAsFactors = FALSE
  )
  
  df_plot <- df_raw[is.finite(df_raw$value), , drop = FALSE]
  
  list(raw = df_raw, plot = df_plot)
}



#' Build an imputed histograms/density summary plot (legacy "imputed_histograms_summary")
#'
#' Produces a single summary figure showing the distribution of observed vs imputed
#' values per sample (faceted), similar in spirit to the legacy pipeline output.
#'
#' @param expr_mat Numeric matrix (features x samples), typically log2.
#' @param flags Logical matrix (features x samples), TRUE where value was imputed.
#' @param cfg Mode config list (used only for labels; safe to pass cfg = NULL).
#' @return A ggplot object.
build_imputed_histograms_summary <- function(expr_mat, imputed_flag, cfg = NULL) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  
  df <- build_imputation_long_df(expr_mat, imputed_flag)
  
  # IMPORTANT: logic check must be on RAW (before filtering finite values)
  if (!any(df$raw$is_imputed, na.rm = TRUE)) {
    return(
      ggplot2::ggplot(df$plot, ggplot2::aes(x = value)) +
        ggplot2::geom_density(na.rm = TRUE) +
        ggplot2::facet_wrap(~ sample, scales = "free_y") +
        ggplot2::labs(
          title = "Imputation QC: no imputed values detected",
          x = "Expression (log2)",
          y = "Density"
        ) +
        ggplot2::theme_bw()
    )
  }
  
  dfp <- df$plot
  dfp$is_imputed <- ifelse(dfp$is_imputed, "Imputed", "Observed")
  
  ggplot2::ggplot(dfp, ggplot2::aes(x = value, fill = is_imputed)) +
    ggplot2::geom_histogram(alpha = 0.6, bins = 60, position = "identity") +
    ggplot2::facet_wrap(~ sample, scales = "free_y") +
    ggplot2::labs(
      title = "Imputation QC: observed vs imputed distributions (per sample)",
      x = "Expression (log2)",
      y = "Count",
      fill = NULL
    ) +
    ggplot2::theme_bw()
}


#' Legacy-style histogram for a single sample (unimputed vs imputed overlay)
#'
#' @param expr_mat Numeric matrix (features x samples), typically log2
#' @param imputed_flag Logical matrix, TRUE where value was imputed
#' @param sample_id Sample column name to plot
#' @param add_x_prefix If TRUE, label sample as "X<sample>" (legacy-like)
#' @return ggplot object
plot_imputation_histogram_one_sample <- function(expr_mat,
                                                 imputed_flag,
                                                 sample_id,
                                                 add_x_prefix = TRUE) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  
  expr_mat <- as.matrix(expr_mat)
  imputed_flag <- as.matrix(imputed_flag)
  
  stopifnot(sample_id %in% colnames(expr_mat))
  stopifnot(all(dim(expr_mat) == dim(imputed_flag)))
  
  v <- expr_mat[, sample_id]
  f <- imputed_flag[, sample_id]
  
  df_raw <- data.frame(
    value = v,
    is_imputed = f,
    stringsAsFactors = FALSE
  )
  
  # Keep only finite values for plotting (do NOT use this for logic decisions)
  df <- df_raw[is.finite(df_raw$value), , drop = FALSE]
  
  df$is_imputed <- ifelse(df$is_imputed, "Imputed", "Observed")
  
  lbl <- if (add_x_prefix) paste0("X", sample_id) else sample_id
  
  ggplot2::ggplot(df, ggplot2::aes(x = value)) +
    # Observed histogram
    ggplot2::geom_histogram(
      data = df[df$is_imputed == "Observed", , drop = FALSE],
      bins = 60,
      alpha = 0.35
    ) +
    # Imputed histogram (overlay)
    ggplot2::geom_histogram(
      data = df[df$is_imputed == "Imputed", , drop = FALSE],
      bins = 60,
      alpha = 0.35
    ) +
    ggplot2::labs(title = NULL, x = NULL, y = NULL) +
    ggplot2::facet_wrap(~ lbl, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "#e9f2ff", colour = NA),
      strip.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )
}
#' Write per-sample imputation histograms (legacy folder style)
#'
#' @param expr_mat Numeric matrix (features x samples)
#' @param imputed_flag Logical matrix (features x samples)
#' @param out_dir Directory to write PNGs into
#' @return Character vector of written file paths
write_imputation_histograms_per_sample <- function(expr_mat, imputed_flag, out_dir) {
  stopifnot(dir.exists(out_dir) || dir.create(out_dir, recursive = TRUE))
  sample_ids <- colnames(expr_mat)
  
  files <- character(0)
  
  for (s in sample_ids) {
    p <- plot_imputation_histogram_one_sample(expr_mat, imputed_flag, s, add_x_prefix = TRUE)
    f <- file.path(out_dir, paste0("X", s, ".png"))
    
    with_png(f, width = 500, height = 900, res = 150, {
      print(p)
    })
    
    files <- c(files, f)
  }
  
  files
}
#' Save a pheatmap object to PNG reliably (grid graphics)
#'
#' @param ph pheatmap object (returned by pheatmap::pheatmap)
#' @param out_file path to png
#' @param width,height,res passed to with_png
save_pheatmap_png <- function(ph, out_file, width = 1600, height = 1200, res = 150) {
  stopifnot(!is.null(ph$gtable))
  with_png(out_file, width = width, height = height, res = res, expr = {
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
  })
  invisible(out_file)
}

