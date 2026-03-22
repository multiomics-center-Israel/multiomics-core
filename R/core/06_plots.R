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
                                 alpha = 0.3,
                                 title = "Density plot of normalized intensities") {
  stopifnot(is.matrix(expr_mat) || is.data.frame(expr_mat))
  expr_mat <- as.data.frame(expr_mat)

  norm_expr_long <- expr_mat|>
    tibble::rownames_to_column("feature")|>
    tidyr::pivot_longer(
      cols = -feature,
      names_to = "SampleID",
      values_to = "value"
    )

  norm_expr_long <- norm_expr_long[is.finite(norm_expr_long$value), , drop = FALSE]

  ggplot2::ggplot(norm_expr_long, ggplot2::aes(x = value, color = SampleID)) +
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
#' @param show_labels show sample names on axes (default FALSE for readability)
#' @param cluster_rows,cluster_cols enable hierarchical clustering (default TRUE)
#' @return invisibly returns pheatmap object
plot_sample_distance_heatmap <- function(expr_mat,
                                         dist_method = "euclidean",
                                         annotation_col = NULL,
                                         main = NULL,
                                         colors = NULL,
                                         fontsize = 12,
                                         show_labels = TRUE,
                                         cluster_rows = TRUE,
                                         cluster_cols = TRUE) {
  expr_mat <- as.matrix(expr_mat)
  sampleDists <- stats::dist(t(expr_mat), method = dist_method)
  mat <- as.matrix(sampleDists)

  if (is.null(colors)) {
    colors <- get_heatmap_colors(255)
  }
  if (is.null(main)) main <- sprintf("Sample distance heatmap (%s)", dist_method)

  # Issue 1 FIX: Suppress axis labels by default to prevent overload
  # Issue 5 FIX: Make clustering configurable
  pheatmap::pheatmap(
    mat,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    annotation_col = annotation_col,
    annotation_row = annotation_col, # Issue 2 FIX: Mirror annotations for symmetry
    main = main,
    col = colors,
    show_rownames = show_labels, # Issue 1 FIX: Hide by default
    show_colnames = show_labels, # Issue 1 FIX: Hide by default
    fontsize_row = fontsize,
    fontsize_col = fontsize,
    annotation_legend = TRUE,
    legend = TRUE,
    border_color = NA # Issue 2 FIX: Remove grid lines for cleaner look
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
                                            fontsize = 12,
                                            show_labels = TRUE,
                                            cluster_rows = TRUE,
                                            cluster_cols = TRUE,
                                            adjust_scale = FALSE) {
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

  # Issue 4 FIX: Improve contrast for high-correlation matrices
  # When correlations are tight (e.g., 0.8-0.95), use focused scale
  breaks <- NULL
  if (adjust_scale) {
    cor_range <- range(cor_mat[lower.tri(cor_mat)], na.rm = TRUE)
    cor_min <- cor_range[1]
    cor_max <- cor_range[2]

    # If range is narrow (typical for good QC data), adjust color scale
    if ((cor_max - cor_min) < 0.3) {
      # Use quantile-based breaks for better visual separation
      q_breaks <- stats::quantile(cor_mat[lower.tri(cor_mat)],
        probs = seq(0, 1, length.out = 256),
        na.rm = TRUE
      )
      breaks <- unique(q_breaks)

      # Regenerate colors to match breaks
      if (length(breaks) > 2) {
        colors <- get_heatmap_colors(length(breaks) - 1)
      }
    }
  }

  # Issue 1, 2, 5 FIX: Hide labels, add row annotations, make clustering configurable
  pheatmap::pheatmap(
    cor_mat,
    annotation_col = annotation_col,
    annotation_row = annotation_col, # Mirror annotations for symmetry
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    main = main,
    col = colors,
    breaks = breaks, # Issue 4 FIX: Adjusted scale
    show_rownames = show_labels, # Issue 1 FIX
    show_colnames = show_labels, # Issue 1 FIX
    fontsize_row = fontsize,
    fontsize_col = fontsize,
    annotation_legend = TRUE,
    legend = TRUE,
    border_color = NA # Issue 2 FIX: Cleaner appearance
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
  args <- list(...)
  args$mat <- as.matrix(expr_mat)
  args$scale <- if (scale_rows) "row" else "none"
  args$cluster_rows <- cluster_rows
  args$cluster_cols <- cluster_cols
  args$show_rownames <- FALSE
  args$annotation_col <- annotation_col
  args$main <- title

  do.call(pheatmap::pheatmap, args)
}
#' Build an imputed histograms/density summary plot (legacy "imputed_histograms_summary")
#'
#' Produces a single summary figure showing the distribution of observed vs imputed
#' values per sample (faceted), similar in spirit to the legacy pipeline output.
#'
#' @param expr_mat Numeric matrix (features x samples), typically log2.
#' @param imputed_flag Logical matrix (features x samples), TRUE where value was imputed.
#' @param width Imputation width parameter (optional, shown in title).
#' @param downshift Imputation downshift parameter (optional, shown in title).
#' @return A ggplot object.
plot_imputation_summary <- function(expr_mat, imputed_flag, width = NULL, downshift = NULL) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))

  df <- build_imputation_long_df(expr_mat, imputed_flag)

  # IMPORTANT: logic check must be on RAW (before filtering finite values)
  if (!any(df$raw$is_imputed, na.rm = TRUE)) {
    return(
      ggplot2::ggplot(df$plot, ggplot2::aes(x = value)) +
        ggplot2::geom_density(na.rm = TRUE) +
        ggplot2::facet_wrap(~sample, scales = "free_y") +
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
    ggplot2::facet_wrap(~sample, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("Imputed" = "#00BFC4", "Observed" = "#F8766D")) +
    ggplot2::labs(
      title = if (!is.null(width) && !is.null(downshift)) {
        sprintf("Imputation QC: observed vs imputed (width = %s, shift = %s)", width, downshift)
      } else {
        "Imputation QC: observed vs imputed distributions (per sample)"
      },
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
      alpha = 0.35,
      fill = "#F8766D"
    ) +
    # Imputed histogram (overlay)
    ggplot2::geom_histogram(
      data = df[df$is_imputed == "Imputed", , drop = FALSE],
      bins = 60,
      alpha = 0.35,
      fill = "#00BFC4"
    ) +
    ggplot2::labs(title = NULL, x = NULL, y = NULL) +
    ggplot2::facet_wrap(~lbl, scales = "free_y") +
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
    stop(
      "plot_pca_scatter(): missing columns in scores: ",
      paste(missing, collapse = ", ")
    )
  }

  aes_args <- list(
    x = rlang::sym(x_col),
    y = rlang::sym(y_col),
    colour = rlang::sym(color_col)
  )

  if (!is.null(shape_col) && shape_col %in% colnames(scores)) {
    aes_args$shape <- rlang::sym(shape_col)
  }

  p <- ggplot2::ggplot(scores, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title  = sprintf("PCA: PC%d vs PC%d", pc_x, pc_y),
      x      = pc_labels[pc_x],
      y      = pc_labels[pc_y],
      colour = color_col,
      shape  = if (!is.null(shape_col)) shape_col else NULL
    ) +
    ggplot2::theme_minimal()

  p
}

plot_cluster_profiles_legacy_style <- function(group_means, clusters, x_label = "Group") {
  gm <- as.matrix(group_means)
  clv <- clusters[rownames(gm)]
  
  nfeat_df <- data.frame(
    cluster = as.integer(names(table(clv))),
    n_features = as.integer(table(clv))
  )
  
  df <- as.data.frame(gm)
  df$gene <- rownames(gm)
  df$cluster <- unname(clv)
  
  long <- tidyr::pivot_longer(
    df,
    cols = setdiff(colnames(df), c("gene", "cluster")),
    names_to = "group",
    values_to = "EXP"
  )
  

  long$group <- factor(long$group, levels = unique(long$group))
  
  # facet labels with n genes
  long <- merge(long, nfeat_df, by = "cluster")
  
  long$facet_label <- sprintf("Cluster %s (n=%d)", long$cluster, long$n_features)
  
  long$facet_label <- factor(
    long$facet_label,
    levels = unique(long$facet_label[order(as.numeric(as.character(long$cluster)))])
  )
  
  ggplot2::ggplot(long, ggplot2::aes(x = group, y = EXP, group = 1)) +
    ggplot2::stat_summary(fun = mean, geom = "line", linewidth = 1.3) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 0.3) +
    ggplot2::facet_wrap(~facet_label, scales = "fixed", ncol = 2) +
    ggplot2::labs(y = "Expression (group means)", x = x_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot cluster profiles using ggplot2
#' Replaces the manual base-R loop for cluster visualization.
#' @param prof_df Data frame containing: cluster, group, mean, sd, n_features
plot_cluster_profiles <- function(prof_df, x_label = "Group") {
  

  prof_df$group <- factor(prof_df$group, levels = unique(prof_df$group))
  
  
  prof_df$facet_label <- sprintf("Cluster %s (n=%d)", prof_df$cluster, prof_df$n_features)
  prof_df$facet_label <- factor(
    prof_df$facet_label,
    levels = unique(prof_df$facet_label[order(as.numeric(as.character(prof_df$cluster)))])
  )

  
  ggplot2::ggplot(prof_df, ggplot2::aes(x = group, y = mean, group = 1)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sd, ymax = mean + sd),
                           width = 0.1, color = "grey50") +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::facet_wrap(~facet_label, scales = "fixed", ncol = 2) +
    ggplot2::labs(y = "Mean (group means)", x = x_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
# R/plots/plot_de.R

#' Volcano plot for a single DE table (one contrast)
#'
#' Consistent logic:
#' Y-axis: -log10(adj.P.Val)  [Adjusted p-value / FDR]
#' Color:  Up (red), Down (blue), NS (grey) based on padj and logFC thresholds
#' H-line: padj_cutoff (matches the Y-axis and coloring threshold)
#'
#' @param de_tbl Data frame with logFC, P.Value, adj.P.Val
#' @param cfg Config list (sections de$p_cutoff, de$linear_fc_cutoff)
#' @param title Plot title
#' @param ... Ignored
#' @return ggplot object
plot_volcano <- function(de_tbl, cfg, title = NULL, ...) {
  # 1. Flexible Column Mapping
  # Try to find the logFC column
  lfc_col <- intersect(c("log2FoldChange", "logFC"), colnames(de_tbl))[1]
  # Try to find the Adjusted P-value column
  padj_col <- intersect(c("padj", "adj.P.Val", "fdr"), colnames(de_tbl))[1]
  # Try to find the raw P-value column
  pval_col <- intersect(c("pvalue", "P.Value", "p.value"), colnames(de_tbl))[1]
  
  # Validation check
  if (is.na(lfc_col) || is.na(padj_col)) {
    stop(paste0(
      "plot_volcano: Could not find required columns. ",
      "Available: ", paste(colnames(de_tbl), collapse = ", ")
    ))
  }
  
  # 2. Thresholds from config
  p_cut <- cfg$de$p_cutoff %||% 0.05
  lin_fc_cut <- cfg$de$linear_fc_cutoff %||% 1.5
  log2fc_cut <- log2(lin_fc_cut)
  
  # 3. Prepare Data
  # We create a local 'plot_df' so we don't mess with the original de_tbl
  plot_df <- as.data.frame(de_tbl)
  plot_df$.logFC <- as.numeric(plot_df[[lfc_col]])
  plot_df$.padj <- as.numeric(plot_df[[padj_col]])
  
  # Handle NAs (important for DESeq2)
  plot_df$.padj_plot <- ifelse(is.na(plot_df$.padj), 1, plot_df$.padj)
  plot_df$.neglog10p <- -log10(pmax(plot_df$.padj_plot, 1e-300))
  
  # 4. Define Significance & Direction
  is_sig <- !is.na(plot_df$.logFC) & (plot_df$.padj_plot <= p_cut) & (abs(plot_df$.logFC) >= log2fc_cut)
  
  plot_df$.direction <- factor(
    ifelse(!is_sig, "NS", ifelse(plot_df$.logFC > 0, "Up", "Down")),
    levels = c("NS", "Down", "Up")
  )
  
  # 5. Sorting (Significant points on top)
  plot_df <- plot_df[order(plot_df$.direction), ]
  
  # 6. Plotting
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .logFC, y = .neglog10p)) +
    ggplot2::geom_point(ggplot2::aes(color = .direction, alpha = .direction), size = 1.2, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c("NS" = "grey80", "Down" = "#377eb8", "Up" = "#e41a1c"),
      drop = FALSE # Keep all levels in legend even if 0 counts
    ) +
    ggplot2::scale_alpha_manual(values = c("NS" = 0.3, "Down" = 0.8, "Up" = 0.8), guide = "none") +
    ggplot2::geom_vline(xintercept = c(-log2fc_cut, log2fc_cut), linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      title = title %||% "Volcano Plot",
      subtitle = paste("Using:", padj_col, "and", lfc_col),
      x = "log2 Fold Change",
      y = paste0("-log10(", padj_col, ")")
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
#' Color: Up (red), Down (blue), NS (grey) based on padj and logFC thresholds
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

  # A column (already log2 scale from DESeq2/edgeR normalization)
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

  # Categorize: Up, Down, NS (not significant)
  is_sig <- !is.na(df$.p) & !is.na(df$.logFC) &
    (df$.p <= p_cut) & (abs(df$.logFC) >= log2fc_cut)
  is_up <- is_sig & (df$.logFC > 0)
  is_down <- is_sig & (df$.logFC < 0)

  df$.direction <- factor(
    ifelse(is_up, "Up", ifelse(is_down, "Down", "NS")),
    levels = c("NS", "Down", "Up")
  )

  # Count stats
  n_up <- sum(is_up, na.rm = TRUE)
  n_down <- sum(is_down, na.rm = TRUE)

  # Sort so significant points are plotted last (on top)
  df <- df[order(df$.direction), ]

  ggplot2::ggplot(df, ggplot2::aes(x = .A, y = .logFC)) +
    ggplot2::geom_point(ggplot2::aes(color = .direction, alpha = .direction), size = 1.2, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      name = "Regulation",
      values = c("NS" = "grey60", "Down" = "blue", "Up" = "red"),
      labels = c(
        "NS" = sprintf("NS (%d)", nrow(df) - n_up - n_down),
        "Down" = sprintf("Down (%d)", n_down),
        "Up" = sprintf("Up (%d)", n_up)
      )
    ) +
    ggplot2::scale_alpha_manual(values = c("NS" = 0.4, "Down" = 0.8, "Up" = 0.8), guide = "none") +
    ggplot2::geom_hline(yintercept = c(-log2fc_cut, log2fc_cut), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    ggplot2::labs(
      title = title %||% "MA plot",
      x = "log2(Mean Expression)",
      y = "log2 Fold Change"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")
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

#' Binary missingness heatmap for MNAR/MAR diagnostics
#'
#' Visualises a binary (missing/observed) matrix with columns sorted by ascending
#' sample intensity (revealing MNAR patterns on the left) and rows colour-annotated
#' by MNAR/MAR class.
#'
#' PURE function: no config access, no file I/O, no side effects.
#'
#' @param mat Numeric matrix (features x samples). NAs denote missing values.
#' @param mnar_class Named logical or character vector, one element per feature.
#'   Logical: TRUE = MNAR, FALSE = MAR/observed.
#'   Character: "MNAR", "MAR", "all_missing", or "all_observed".
#'   Names must correspond to rownames(mat).
#' @param sample_order Optional character vector of column names specifying column
#'   order. If NULL, columns are sorted by ascending \code{colMeans(mat, na.rm=TRUE)}
#'   (low-intensity first).
#' @param feat_order Optional character vector of row names specifying row order.
#'   If NULL, rows are sorted by descending missingness fraction.
#'
#' @return pheatmap object.
#'
plot_missingness_heatmap <- function(mat, mnar_class,
                                     sample_order = NULL,
                                     feat_order = NULL) {
  mat <- as.matrix(mat)

  # Normalise mnar_class to character
  if (is.logical(mnar_class)) {
    class_vec <- ifelse(mnar_class, "MNAR", "MAR")
  } else {
    class_vec <- as.character(mnar_class)
  }
  if (is.null(names(class_vec))) names(class_vec) <- rownames(mat)

  # Binary missingness matrix (1 = missing, 0 = observed)
  miss_bin <- matrix(
    as.integer(is.na(mat)),
    nrow = nrow(mat), ncol = ncol(mat),
    dimnames = dimnames(mat)
  )

  # Column order: ascending colMeans → low-intensity samples on the left
  if (is.null(sample_order)) {
    col_means    <- colMeans(mat, na.rm = TRUE)
    sample_order <- names(sort(col_means, decreasing = FALSE))
  }
  sample_order <- intersect(sample_order, colnames(miss_bin))
  miss_bin     <- miss_bin[, sample_order, drop = FALSE]

  # Row order: descending missingness fraction → most-missing features at top
  if (is.null(feat_order)) {
    row_miss_pct <- rowMeans(is.na(mat[, sample_order, drop = FALSE]))
    feat_order   <- names(sort(row_miss_pct, decreasing = TRUE))
  }
  feat_order <- intersect(feat_order, rownames(miss_bin))
  miss_bin   <- miss_bin[feat_order, , drop = FALSE]

  # Row annotation aligned to the (possibly reordered) features
  class_aligned <- class_vec[feat_order]
  annot_row     <- data.frame(Class = class_aligned, row.names = feat_order)

  class_colors <- c(
    "MNAR"         = "#d73027",  # red
    "MAR"          = "#4575b4",  # blue
    "all_missing"  = "#313695",  # dark blue
    "all_observed" = "#91cf60"   # green
  )
  present <- intersect(names(class_colors), unique(class_aligned))
  annot_colors <- list(Class = class_colors[present])

  show_rownames <- nrow(miss_bin) <= 80
  show_colnames <- ncol(miss_bin) <= 40

  pheatmap::pheatmap(
    miss_bin,
    color              = c("#f7f7f7", "#252525"),  # white = observed, dark = missing
    breaks             = c(-0.5, 0.5, 1.5),
    cluster_rows       = FALSE,
    cluster_cols       = FALSE,
    annotation_row     = annot_row,
    annotation_colors  = annot_colors,
    main               = "Missingness pattern (dark = missing, sorted by sample intensity)",
    show_rownames      = show_rownames,
    show_colnames      = show_colnames,
    border_color       = NA,
    legend             = TRUE,
    legend_breaks      = c(0, 1),
    legend_labels      = c("Observed", "Missing")
  )
}

#' Save a pheatmap object to file
save_heatmap_to_file <- function(pheatmap_obj, out_file,
                                 width = 2000, height = 1400, res = 150) {
  if (is.null(out_file)) {
    return(invisible(NULL))
  }

  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  if (!grepl("\\.png$", out_file, ignore.case = TRUE)) out_file <- paste0(out_file, ".png")

  # Task 3: Increase size to prevent legend cutoff with multi-column annotations
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
    sample = rep(colnames(expr_mat), each = nrow(expr_mat)),
    value = as.vector(expr_mat),
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
                                    annotation_row_builder = FALSE, # function(use_ids, context) -> df
                                    annotation_row_context = NULL,
                                    out_file = NULL,
                                    title = NULL,
                                    cluster_cols = TRUE,
                                    cluster_rows_default = TRUE,
                                    scale_rows = TRUE) {
  
  use_ids <- intersect(feature_ids, rownames(expr_mat))
  
  if (!is.null(ordering)) {
    use_ids <- ordering[ordering %in% use_ids]
    cluster_rows_flag <- FALSE
  } else {
    cluster_rows_flag <- cluster_rows_default
  }
  
  mat2plot <- expr_mat[use_ids, , drop = FALSE]
  
  sample_col <- cfg$effects$samples
  color_col  <- get_color_config(cfg)
  
  annot_col <- build_heatmap_annotation_col(meta, cfg)
  # annot_col <- data.frame(
  #   Condition = meta[[color_col]],
  #   row.names = meta[[sample_col]]
  # )
  # annot_col <- annot_col[colnames(mat2plot), , drop = FALSE]
  # 
  annot_row <- NULL
  if (annotation_row_builder) {
    annot_row <- build_contrast_row_context(use_ids, annotation_row_context)
  }
  
  if (is.null(title)) title <- sprintf("Heatmap (%d features)", nrow(mat2plot))
  
  ph <- plot_heatmap_core(
    expr_mat       = mat2plot,
    annotation_col = annot_col,
    annotation_row = annot_row,
    title          = title,
    scale_rows     = scale_rows,
    cluster_rows   = cluster_rows_flag,
    cluster_cols   = cluster_cols,
    max_rows       = NULL
  )
  
  if (!is.null(out_file)) save_heatmap_to_file(ph, out_file)
  ph
}


build_contrast_row_context <- function(use_ids, context) {
  summary_df    <- context$summary_df
  p_cutoff      <- context$p_cutoff %||% 0.05
  log2fc_cutoff <- context$log2fc_cutoff %||% 0.585
  id_col        <- context$id_col %||% "FeatureID"
  
  if (is.null(summary_df)) return(NULL)
  
  ar <- build_de_row_annotations(
    summary_df    = summary_df,
    feature_ids   = use_ids,
    p_cutoff      = p_cutoff,
    log2fc_cutoff = log2fc_cutoff,
    id_col        = id_col
  )
  if (is.null(ar)) return(NULL)
  
  # Expand to full set of rows (avoid breaking when some ids missing)
  full <- as.data.frame(
    matrix(NA_character_, nrow = length(use_ids), ncol = ncol(ar)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(full) <- colnames(ar)
  rownames(full) <- use_ids
  
  common <- intersect(use_ids, rownames(ar))
  full[common, ] <- ar[common, , drop = FALSE]
  full
}

build_heatmap_annotation_col <- function(meta, cfg) {
  sample_col <- cfg$effects$samples %||% "sample_id"
  color_col  <- cfg$effects$color %||% NULL
  shape_col  <- cfg$effects$shape %||% NULL
  
  annot_cols <- c(
    if (!is.null(color_col) && all(color_col %in% colnames(meta))) color_col,
    if (!is.null(shape_col) && all(shape_col %in% colnames(meta)) && all(shape_col != color_col)) shape_col
  )
  
  if (length(annot_cols) == 0) return(NULL)
  
  annot <- meta[, annot_cols, drop = FALSE]
  rownames(annot) <- meta[[sample_col]]
  annot
}
