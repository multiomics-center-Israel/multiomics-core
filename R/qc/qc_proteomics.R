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


# helper 
build_imputation_long_df <- function(imputed, imputed_flag) {
  imputed      <- as.matrix(imputed)
  imputed_flag <- as.matrix(imputed_flag)
  
  sample_ids <- colnames(imputed)
  
  df <- data.frame(
    sample = rep(sample_ids, each = nrow(imputed)),
    value  = as.vector(imputed),
    flag   = as.vector(imputed_flag),
    stringsAsFactors = FALSE
  )
  df$type <- ifelse(df$flag, "Imputed", "Observed")
  df <- df[is.finite(df$value), , drop = FALSE]
  df
}


#' Imputed data histograms per sample (facet by sample)
#'
#' @param imputed      numeric matrix (features x samples), after imputation.
#' @param imputed_flag logical matrix, TRUE where value was imputed.
#' @param cfg          config$modes$proteomics (for title: width/downshift if desired).
#' @param out_file     optional path to save plot.
#'
#' @return invisibly returns the ggplot object.
imputed_histograms_summary <- function(imputed,
                                       imputed_flag,
                                       cfg = NULL,
                                       out_file = NULL) {
 
  
  imputed      <- as.matrix(imputed)
  imputed_flag <- as.matrix(imputed_flag)
  
  sample_ids <- colnames(imputed)
  df <- build_imputation_long_df(imputed, imputed_flag)
  

  
  title_txt <- "Imputed data histograms summary"
  if (!is.null(cfg) && !is.null(cfg$imputation)) {
    w  <- cfg$imputation$width     %||% NA
    ds <- cfg$imputation$downshift %||% NA
    if (!is.na(w) && !is.na(ds)) {
      title_txt <- sprintf(
        "Imputed data histograms summary (WIDTH=%.2f, DOWNSHIFT=%.2f)",
        w, ds
      )
    }
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = type)) +
    ggplot2::geom_histogram(alpha = 0.7, bins = 60, position = "identity") +
    ggplot2::facet_wrap(~ sample, nrow = 1) +
    ggplot2::labs(
      title = title_txt,
      x     = "log2 intensity",
      y     = "Frequency",
      fill  = "Value"
    ) +
    ggplot2::theme_minimal()
  
  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file,
                    plot  = p,
                    width = max(7, length(sample_ids) * 1.5),
                    height = 4)
  }
  
  invisible(p)
}


#' Expression heatmap with sample annotations
#'
#' @param expr_mat numeric matrix (features x samples), usually normalized/imputed.
#' @param meta     metadata table (one row per sample).
#' @param cfg      config$modes$proteomics (effects$samples + effects$color).
#' @param title    plot title.
#' @param out_file optional path to save plot (PDF/PNG).
#' @param cluster_rows logical, cluster proteins.
#' @param cluster_cols logical, cluster samples.
#' @param max_rows optional: if not NULL, sample this many rows for speed.
qc_expr_heatmap <- function(expr_mat, meta, cfg, out_file = NULL, title = NULL, max_rows = 2000) {
  
  eff <- cfg$effects
  sample_col <- eff$samples
  color_col  <- eff$color
  
  sample_ids <- colnames(expr_mat)
  meta_sub <- meta[match(sample_ids, meta[[sample_col]]), , drop = FALSE]
  
  annot <- data.frame(
    Condition = meta_sub[[color_col]],
    row.names = sample_ids
  )
  
  # 1) Save to file (no screen)
  if (!is.null(out_file)) {
    with_png(out_file, {
      plot_expr_heatmap(expr_mat, annotation_col = annot, max_rows = max_rows, main = title)
    })
  }
  
  # 2) Always draw to screen and return object
  ph <- plot_expr_heatmap(expr_mat, annotation_col = annot, max_rows = max_rows, main = title)
  invisible(ph)
}



qc_sample_distance_heatmap <- function(expr_mat, meta, cfg, dist_method = "euclidean", out_file = NULL) {
  
  eff <- cfg$effects
  sample_col <- eff$samples
  color_col  <- eff$color
  
  sample_ids <- colnames(expr_mat)
  meta_sub <- meta[match(sample_ids, meta[[sample_col]]), , drop = FALSE]
  
  annot <- data.frame(
    Condition = meta_sub[[color_col]],
    row.names = sample_ids
  )
  
  if (!is.null(out_file)) {
    with_png(out_file, {
      plot_sample_distance_heatmap(expr_mat, dist_method = dist_method, annotation_col = annot)
    })
  }
  
  ph <- plot_sample_distance_heatmap(expr_mat, dist_method = dist_method, annotation_col = annot)
  invisible(ph)
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






