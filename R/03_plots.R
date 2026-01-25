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
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = value, color = SampleID)) +
    ggplot2::geom_density(alpha = alpha, linewidth = 0.7) +
    ggplot2::labs(
      title = title,
      x     = "log2 intensity",
      y     = "Density"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "right",
      legend.title    = ggplot2::element_blank()
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
                                         fontsize = 12) {
  expr_mat <- as.matrix(expr_mat)
  sampleDists <- stats::dist(t(expr_mat), method = dist_method)
  mat <- as.matrix(sampleDists)
  
  if (is.null(colors)) {
    colors <- get_heatmap_colors(255)
  }
  if (is.null(main)) main <- sprintf("Sample distance heatmap (%s)", dist_method)
  
  pheatmap::pheatmap(
    mat,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    annotation_col = annotation_col,
    main = main,
    col = colors,
    fontsize_row = fontsize,
    fontsize_col = fontsize
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
  
  if (is.null(colors)) colors <- get_heatmap_colors(255)
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
#' Core wrapper for pheatmap
#' @return A pheatmap object
plot_heatmap_core <- function(expr_mat, 
                              annotation_col = NULL, 
                              title = NULL, 
                              scale_rows = TRUE,
                              cluster_cols = TRUE,
                              cluster_rows = TRUE,
                              max_rows = NULL,
                              ...) { # ... allows passing extra pheatmap args
  
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Need pheatmap")
  
  # 1. Subsampling if too large (Optimization)
  if (!is.null(max_rows) && nrow(expr_mat) > max_rows) {
    message(sprintf("Subsampling heatmap from %d to %d rows", nrow(expr_mat), max_rows))
    set.seed(42)
    expr_mat <- expr_mat[sample(seq_len(nrow(expr_mat)), max_rows), , drop = FALSE]
  }
  
  # 2. Title default
  if (is.null(title)) title <- sprintf("Heatmap (%d features)", nrow(expr_mat))
  
  # 3. Draw
  pheatmap::pheatmap(
    mat = as.matrix(expr_mat),
    scale = if (scale_rows) "row" else "none",
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = FALSE,
    annotation_col = annotation_col,
    main = title,
    ...
  )
}
# TODO: add width and downshift to the plot title
#' Build an imputed histograms/density summary plot (legacy "imputed_histograms_summary")
#'
#' Produces a single summary figure showing the distribution of observed vs imputed
#' values per sample (faceted), similar in spirit to the legacy pipeline output.
#'
#' @param expr_mat Numeric matrix (features x samples), typically log2.
#' @param flags Logical matrix (features x samples), TRUE where value was imputed.
#' @return A ggplot object.
plot_imputation_summary <- function(expr_mat, imputed_flag) {
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
  df$lbl <- lbl
  
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
plot_pca_scatter <- function(scores, color_col, shape_col = NULL,
                             pc_x, pc_y, pc_labels) {
  # Resolve PC column names (e.g., PC1, PC2)
  x_col <- paste0("PC", pc_x)
  y_col <- paste0("PC", pc_y)
  
  # Validate required columns exist
  missing <- setdiff(c(x_col, y_col, color_col), colnames(scores))
  if (length(missing) > 0) {
    stop("plot_pca_scatter(): missing columns in scores: ",
         paste(missing, collapse = ", "))
  }
  
  aes_args <- list(
    x = rlang::sym(x_col),
    y = rlang::sym(y_col),
    colour = rlang::sym(color_col)
  )
  
  if (!is.null(shape_col) && shape_col %in% colnames(scores)) {
    aes_args$shape <- rlang::sym(shape_col)
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
#' Plot cluster profiles using ggplot2
#' Replaces the manual base-R loop for cluster visualization.
#' @param prof_df Data frame containing: cluster, group, mean, sd, n_features
plot_cluster_profiles <- function(prof_df, x_label = "Group") {
  # Create a clean label for facets
  prof_df$facet_label <- sprintf("Cluster %s (n=%d)", prof_df$cluster, prof_df$n_features)
  
  # Ensure order matches cluster number
  prof_df$facet_label <- factor(prof_df$facet_label, 
                                levels = unique(prof_df$facet_label[order(as.numeric(as.character(prof_df$cluster)))]))
  
  p <- ggplot2::ggplot(prof_df, aes(x = group, y = mean, group = 1)) + 
    # Error bars (SD)
    ggplot2::geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1, color = "grey50") +
    # Line and points
    ggplot2::geom_line(color = "blue") +
    ggplot2::geom_point(size = 2, color = "darkblue") +
    # Zero line
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    # Faceting
    ggplot2::facet_wrap(~ facet_label, scales = "fixed", ncol = 2) + 
    # Styling
    ggplot2::labs(y = "Mean z-score (group means)", x = x_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# R/plots/plot_de.R

#' Volcano plot for a single DE table (one contrast)
#'
#' Expects columns:
#'   - logFC
#'   - P.Value and/or adj.P.Val
#'
#' Uses cfg$de$p_cutoff and cfg$de$linear_fc_cutoff (linear FC) to mark "pass".
#' Converts linear_fc_cutoff to log2FC cutoff.
#'
#' @param de_tbl data.frame with DE results for one contrast
#' @param cfg mode config (e.g., config$modes$proteomics)
#' @param title optional plot title
#' @param use_adj logical; if TRUE uses adj.P.Val, else P.Value
#'
#' @return ggplot object
plot_volcano <- function(de_tbl, cfg, title = NULL, use_adj = TRUE) {
  stopifnot(is.data.frame(de_tbl))
  if (!("logFC" %in% colnames(de_tbl))) stop("plot_volcano: de_tbl missing 'logFC' column.")
  
  p_col <- if (isTRUE(use_adj) && "adj.P.Val" %in% colnames(de_tbl)) "adj.P.Val" else "P.Value"
  if (!(p_col %in% colnames(de_tbl))) {
    stop("plot_volcano: de_tbl missing p-value column. Need 'P.Value' or 'adj.P.Val'.")
  }
  
  # thresholds from config (with safe defaults)
  p_cut <- cfg$de$p_cutoff %||% 0.1
  lin_fc_cut <- cfg$de$linear_fc_cutoff %||% 1.5
  log2fc_cut <- log2(lin_fc_cut)
  
  df <- de_tbl
  df$.p <- as.numeric(df[[p_col]])
  df$.logFC <- as.numeric(df[["logFC"]])
  
  # pass flag (robust to NAs)
  df$.pass <- !is.na(df$.p) & !is.na(df$.logFC) &
    (df$.p <= p_cut) & (abs(df$.logFC) >= log2fc_cut)
  
  df$.neglog10p <- -log10(pmax(df$.p, .Machine$double.xmin))
  
  ggplot2::ggplot(df, ggplot2::aes(x = .logFC, y = .neglog10p)) +
    ggplot2::geom_point(ggplot2::aes(color = .pass), size = 1.5, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "red"), guide = "none") +
    ggplot2::geom_vline(xintercept = c(-log2fc_cut, log2fc_cut), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = "dashed") +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.9), guide = "none") +
    ggplot2::labs(
      title = title %||% "Volcano plot",
      x = "log2 Fold Change",
      y = paste0("-log10(", p_col, ")")
    ) +
    ggplot2::theme_minimal()
}

#' MA plot for a single DE table (one contrast)
#'
#' Expects columns:
#'   - logFC
#'   - AveExpr (preferred) OR .A provided externally
#'   - P.Value and/or adj.P.Val
#'
#' @param de_tbl data.frame with DE results for one contrast
#' @param cfg mode config (e.g., config$modes$proteomics)
#' @param title optional plot title
#' @param use_adj logical; if TRUE uses adj.P.Val, else P.Value
#'
#' @return ggplot object
plot_ma <- function(de_tbl, cfg, title = NULL, use_adj = TRUE) {
  stopifnot(is.data.frame(de_tbl))
  if (!("logFC" %in% colnames(de_tbl))) stop("plot_ma: de_tbl missing 'logFC' column.")
  
  # p column
  p_col <- if (isTRUE(use_adj) && "adj.P.Val" %in% colnames(de_tbl)) "adj.P.Val" else "P.Value"
  if (!(p_col %in% colnames(de_tbl))) stop("plot_ma: need 'P.Value' or 'adj.P.Val'.")
  
  # A column
  if ("AveExpr" %in% colnames(de_tbl)) {
    A <- as.numeric(de_tbl$AveExpr)
  } else if (".A" %in% colnames(de_tbl)) {
    A <- as.numeric(de_tbl$.A)
  } else {
    stop("plot_ma: missing 'AveExpr' (or computed '.A').")
  }
  
  # thresholds
  p_cut <- cfg$de$p_cutoff %||% 0.1
  lin_fc_cut <- cfg$de$linear_fc_cutoff %||% 1.5
  log2fc_cut <- log2(lin_fc_cut)
  
  df <- de_tbl
  df$.A <- A
  df$.logFC <- as.numeric(df$logFC)
  df$.p <- as.numeric(df[[p_col]])
  
  df$.pass <- !is.na(df$.p) & !is.na(df$.logFC) &
    (df$.p <= p_cut) & (abs(df$.logFC) >= log2fc_cut)
  
  ggplot2::ggplot(df, ggplot2::aes(x = .A, y = .logFC)) +
    ggplot2::geom_point(ggplot2::aes(color = .pass), alpha = 0.7, size = 1.2, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = c(-log2fc_cut, log2fc_cut), linetype = "dashed") +
    ggplot2::labs(
      title = title %||% "MA plot",
      x = "A (mean expression)",
      y = "M (log2 Fold Change)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "red"), guide = "none")
}

# ---- Helper ----

add_A_from_expr <- function(de_tbl, expr_mat, id_col = "FeatureID") {
  stopifnot(is.data.frame(de_tbl), is.matrix(expr_mat))
  if (!(id_col %in% colnames(de_tbl))) stop("add_A_from_expr: missing id_col: ", id_col)
  
  ids <- as.character(de_tbl[[id_col]])
  m <- match(ids, rownames(expr_mat))
  A <- rep(NA_real_, length(ids))
  ok <- !is.na(m)
  A[ok] <- rowMeans(expr_mat[m[ok], , drop = FALSE], na.rm = TRUE)
  
  de_tbl$.A <- A
  de_tbl
}

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

wrap_clustering_heatmap <- function(expr_mat, meta, cfg,
                                    feature_ids,
                                    ordering = NULL,
                                    out_file = NULL) {
  
  # 1) choose features present in matrix
  use_ids <- intersect(feature_ids, rownames(expr_mat))
  
  # preserve ordering if provided
  if (!is.null(ordering)) {
    use_ids <- ordering[ordering %in% use_ids]
    cluster_rows_flag <- FALSE
  } else {
    cluster_rows_flag <- TRUE
  }
  
  mat2plot <- expr_mat[use_ids, , drop = FALSE]
  
  # 2) annotation (aligned to matrix columns)
  sample_col <- cfg$effects$samples
  color_col  <- cfg$effects$color
  
  annot <- data.frame(
    Condition = meta[[color_col]],
    row.names = meta[[sample_col]]
  )
  
  # align annotation to columns actually in the matrix
  annot <- annot[colnames(mat2plot), , drop = FALSE]
  
  # 3) plot
  ph <- plot_heatmap_core(
    expr_mat = mat2plot,
    annotation_col = annot,
    title = sprintf("Hierarchical Clustering (%d DE features)", nrow(mat2plot)),
    scale_rows = TRUE,
    cluster_rows = cluster_rows_flag,
    cluster_cols = TRUE,
    max_rows = NULL
  )
  
  # 4) save & return
  if (!is.null(out_file)) save_heatmap_to_file(ph, out_file)
  ph
}





