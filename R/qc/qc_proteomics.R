#' Boxplots of normalized expression per sample
#'
#' @param expr_norm numeric matrix/data.frame (features x samples), normalized (log2).
#' @param meta      metadata table with one row per sample.
#' @param cfg       config$modes$proteomics (used to get color/label columns).
#' @param out_file  optional path to save plot (PDF/PNG). If NULL, only prints.
#'
#' @return invisibly returns the ggplot object.
norm_boxplot <- function(expr_norm, meta, cfg, out_file = NULL) {
  
  expr_norm <- as.matrix(expr_norm)
  
  eff        <- cfg$effects
  sample_col <- eff$samples   # e.g. "SampleID"
  color_col  <- eff$color   # e.g. "Condition"
  
  stopifnot(sample_col %in% colnames(meta))
  sample_ids <- colnames(expr_norm)
  
  meta_sub <- meta[match(sample_ids, meta[[sample_col]]), , drop = FALSE]
  if (any(is.na(meta_sub[[sample_col]]))) {
    stop("Some column names in expr_norm were not found in metadata[[", sample_col, "].")
  }
  
  df_long <- data.frame(
    sample = rep(sample_ids, each = nrow(expr_norm)),
    value  = as.vector(expr_norm),
    stringsAsFactors = FALSE
  )
  # drop NA / NaN / Inf / -Inf
  df_long <- df_long[is.finite(df_long$value), , drop = FALSE]
  
  df_long[[color_col]] <- meta_sub[[color_col]][match(df_long$sample, meta_sub[[sample_col]])]
  
  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = sample, y = value, colour = .data[[color_col]])
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.4) +
    ggplot2::labs(
      title  = "Normalized expression boxplots",
      x      = "Sample",
      y      = "log2(normalized intensity)",
      colour = color_col
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
#'
#' @param expr_norm numeric matrix/data.frame (features x samples), normalized (log2).
#' @param meta      metadata table with one row per sample.
#' @param cfg       config$modes$proteomics (uses effects$samples for SampleID and effects$color for condition).
#' @param out_file  optional path to save plot (PDF/PNG).
#'
#' @return invisibly returns the ggplot object.
norm_histogram_summary <- function(expr_norm, meta, cfg, out_file = NULL) {
  
  expr_norm <- as.matrix(expr_norm)
  
  eff        <- cfg$effects
  sample_col <- eff$samples   # e.g. "SampleID"
  cond_col   <- eff$color   # e.g. "Condition"
  
  stopifnot(all(c(sample_col, cond_col) %in% colnames(meta)))
  
  sample_ids <- colnames(expr_norm)
  meta_sub   <- meta[match(sample_ids, meta[[sample_col]]), , drop = FALSE]
  
  if (any(is.na(meta_sub[[sample_col]]))) {
    stop("Some expr_norm column names not found in metadata[[", sample_col, "].")
  }
  
  df_long <- data.frame(
    sample = rep(sample_ids, each = nrow(expr_norm)),
    value  = as.vector(expr_norm),
    stringsAsFactors = FALSE
  )
  df_long[[cond_col]] <- meta_sub[[cond_col]][match(df_long$sample, meta_sub[[sample_col]])]
  
  df_long <- df_long[is.finite(df_long$value), , drop = FALSE]
  
  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = value, fill = .data[[cond_col]])
  ) +
    ggplot2::geom_histogram(alpha = 0.6, bins = 60, position = "identity") +
    ggplot2::facet_wrap(as.formula(paste("~", cond_col)), nrow = 1, scales = "free_y") +
    ggplot2::labs(
      title = "Normalized expression histograms by condition",
      x     = "log2(normalized intensity)",
      y     = "Frequency",
      fill  = cond_col
    ) +
    ggplot2::theme_minimal()
  
  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file, plot = p, width = 8, height = 4)
  }
  
  invisible(p)
}

#' Imputed data histograms per sample (facet by sample)
#'
#' @param imputed      numeric matrix (features x samples), after imputation.
#' @param imputed_flag logical matrix, TRUE where value was imputed.
#' @param cfg          config$modes$proteomics (for title: width/downshift if desired).
#' @param out_file     optional path to save plot.
#'
#' @return invisibly returns the ggplot object.
imputed_histograms_summary <- function(imputed, imputed_flag, cfg = NULL, out_file = NULL) {
  p <- build_imputed_histograms_summary(imputed, imputed_flag, cfg = cfg)
  
  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file, plot = p, width = 12, height = 5, dpi = 150)
  }
  
  invisible(p)
}

#' Expression heatmap with sample annotations
#'
#' @param expr_mat numeric matrix (features x samples), usually normalized/imputed.
#' @param meta     metadata table (one row per sample).
#' @param cfg      config$modes$proteomics (effects$samples + effects$color).
#' @param out_file optional path to save plot (PNG).
#' @param title    plot title.
#' @param max_rows optional: if not NULL, sample this many rows for speed.
#' @param cluster_cols logical: cluster samples (TRUE = legacy default).
#'
qc_expr_heatmap <- function(expr_mat,
                            meta,
                            cfg,
                            out_file = NULL,
                            title = NULL,
                            max_rows = 2000,
                            cluster_cols = TRUE) {
  
  eff <- cfg$effects
  sample_col <- eff$samples
  color_col  <- eff$color
  
  expr_mat <- as.matrix(expr_mat)
  sample_ids <- colnames(expr_mat)
  
  meta <- as.data.frame(meta)
  meta_sub <- meta[match(sample_ids, meta[[sample_col]]), , drop = FALSE]
  
  annot <- data.frame(
    Condition = meta_sub[[color_col]],
    row.names = sample_ids
  )
  
  plot_fun <- function() {
    plot_expr_heatmap(
      expr_mat        = expr_mat,
      annotation_col  = annot,
      max_rows        = max_rows,
      main            = title,
      # NEW: allow disabling column clustering for the "wo_col" version
      cluster_cols    = cluster_cols
    )
  }
  
  if (!is.null(out_file)) {
    ph <- plot_fun()
    save_pheatmap_png(ph, out_file, width = 1600, height = 1200, res = 150)
  }
  
  invisible(plot_fun())
  
}

#' Density overlay of normalized expression for all samples
#'
#' One density curve per sample, all on the same plot.
#'
#' @param expr_mat numeric matrix/data.frame (features x samples).
#' @param meta      metadata table (one row per sample).
#' @param cfg       config$modes$proteomics (effects$samples + effects$color).
#' @param out_file  optional path to save plot (PDF/PNG). If NULL, only prints.
#'
#' @return invisibly returns the ggplot object.
qc_proteomics_density <- function(expr_mat,
                                  meta,
                                  cfg,
                                  out_file   = NULL,
                                  alpha      = 0.3,
                                  show_legend = TRUE,
                                  title      = "Density plot of normalized intensities") {

  
  sample_ids <- colnames(expr_mat)
  
  p <- plot_density_overlay(
    expr_mat = expr_mat,
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

#' Sample–sample correlation heatmap (proteomics QC)
#'
#' Computes a sample correlation matrix (Pearson) and plots it as a heatmap.
#' This is complementary to the sample distance heatmap.
#'
#' @param expr_mat numeric matrix (features x samples), typically log2 and imputed.
#' @param meta     metadata table (one row per sample). Used only for annotation (optional).
#' @param cfg      config$modes$proteomics (effects$samples + effects$color).
#' @param out_file optional path to save PNG.
#' @param method   correlation method ("pearson", "spearman").
#' @return invisibly returns the pheatmap object.
qc_sample_correlation_heatmap <- function(expr_mat,
                                          meta,
                                          cfg,
                                          out_file,
                                          method = "pearson",
                                          fontsize = 12) {
  eff <- cfg$effects
  sample_col <- eff$samples
  color_col  <- eff$color
  
  expr_mat <- as.matrix(expr_mat)
  meta <- as.data.frame(meta)
  sample_ids <- colnames(expr_mat)
  
  meta_sub <- meta[match(sample_ids, meta[[sample_col]]), , drop = FALSE]
  
  annot <- data.frame(
    Condition = meta_sub[[color_col]],
    row.names = sample_ids
  )
  
  ph <- plot_sample_correlation_heatmap(
    expr_mat = expr_mat,
    method = method,
    annotation_col = annot,
    fontsize = fontsize
  )
  
  save_pheatmap_png(ph, out_file, width = 1600, height = 1200, res = 150)
  invisible(ph)
}


#' Sample–sample distance heatmap (QC)
#'
#' Generates a sample–sample distance heatmap from a proteomics expression matrix
#' and saves it to file.
#'
#' Distances are computed between samples based on their expression profiles
#' (features x samples matrix). This plot is primarily a **technical QC tool**,
#' useful for detecting outlier samples, batch effects, or unexpected clustering.
#'
#' Two modes are supported:
#' - with_na = FALSE (default): rows (proteins) containing any NA values are removed
#'   prior to distance computation.
#' - with_na = TRUE: NA values are kept; distances are computed on available values.
#'
#' @param expr_mat Numeric matrix (features x samples), typically log2 normalized
#'        and optionally imputed expression values.
#' @param meta Data frame with one row per sample (sample metadata).
#' @param cfg Proteomics mode configuration (`config$modes$proteomics`),
#'        used to extract sample ID and color annotations.
#' @param out_file Path to output PNG file.
#' @param with_na Logical; whether to keep rows with missing values (default: FALSE).
#' @param fontsize Numeric font size for row/column labels.
#'
#' @return Invisibly returns the pheatmap object.
#'
#' @seealso plot_sample_distance_heatmap
#'
qc_sample_distance_heatmap <- function(expr_mat,
                                       meta,
                                       cfg,
                                       out_file,
                                       with_na = FALSE,
                                       fontsize = 12) {
  
  eff <- cfg$effects
  sample_col <- eff$samples
  color_col  <- eff$color
  
  expr_mat <- as.matrix(expr_mat)
  meta <- as.data.frame(meta)
  sample_ids <- colnames(expr_mat)
  meta_sub <- meta[match(sample_ids, meta[[sample_col]]), , drop = FALSE]
  
  annot <- data.frame(
    Condition = meta_sub[[color_col]],
    row.names = sample_ids
  )
  
  if (!with_na) {
    keep <- stats::complete.cases(expr_mat)
    expr_mat <- expr_mat[keep, , drop = FALSE]
  }
  
  ph <- plot_sample_distance_heatmap(
    expr_mat = expr_mat,
    annotation_col = annot,
    fontsize = fontsize
  )
  
  save_pheatmap_png(ph, out_file, width = 1600, height = 1200, res = 150)
  invisible(ph)
}




