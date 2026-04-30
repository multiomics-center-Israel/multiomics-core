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

  norm_expr_long <- expr_mat %>%
    tibble::rownames_to_column("feature") %>%
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
    # Light-to-dark Blues: small distance = light, large distance = dark
    # (Apr 2026 reviewer feedback — invert the legacy palette).
    colors <- get_heatmap_colors(255, reverse = FALSE)
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
                                            adjust_scale = TRUE) {
  expr_mat <- as.matrix(expr_mat)

  cor_mat <- stats::cor(
    expr_mat,
    use = "pairwise.complete.obs",
    method = method
  )

  # High correlation = darker color (diagonal = 1 = darkest).
  if (is.null(colors)) colors <- get_heatmap_colors(255, reverse = FALSE)
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

      # Regenerate colors to match breaks (high corr = dark)
      if (length(breaks) > 2) {
        colors <- get_heatmap_colors(length(breaks) - 1, reverse = FALSE)
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
#' Plot cluster profiles using ggplot2
#' Replaces the manual base-R loop for cluster visualization.
#' @param prof_df Data frame containing: cluster, group, mean, sd, n_features
plot_cluster_profiles <- function(prof_df, x_label = "Group") {
  # Create a clean label for facets
  prof_df$facet_label <- sprintf("Cluster %s (n=%d)", prof_df$cluster, prof_df$n_features)

  # Ensure order matches cluster number
  prof_df$facet_label <- factor(prof_df$facet_label,
    levels = unique(prof_df$facet_label[order(as.numeric(as.character(prof_df$cluster)))])
  )

  p <- ggplot2::ggplot(prof_df, ggplot2::aes(x = group, y = mean, group = 1)) +
    # Error bars (SD)
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sd, ymax = mean + sd), width = 0.1, color = "grey50") +
    # Line and points
    ggplot2::geom_line(color = "blue") +
    ggplot2::geom_point(size = 2, color = "darkblue") +
    # Zero line
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    # Faceting
    ggplot2::facet_wrap(~facet_label, scales = "fixed", ncol = 2) +
    # Styling
    ggplot2::labs(y = "Mean z-score (group means)", x = x_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  return(p)
}

#' Build per-cluster profile ggplots from sample-level long data
#'
#' Each plot shows the mean expression trend across groups with SE error bars.
#' Supports optional secondary grouping (color/linetype) when a ColorGroup
#' column is present in the input data.
#'
#' @param long_df data.frame with columns: Gene, Name, Exp, Group, Cluster,
#'   and optionally ColorGroup
#' @param x_label Character: X-axis label (group column name)
#' @param color_label Character or NULL: legend label for color/linetype
#'   grouping. NULL means no secondary grouping.
#' @return Named list of ggplot objects, keyed by cluster number (as character)
#' @export
build_cluster_profile_plots <- function(long_df, x_label = "Group",
                                         color_label = NULL) {
  has_color <- "ColorGroup" %in% colnames(long_df) && !is.null(color_label)

  long_df$Group <- factor(long_df$Group, levels = unique(long_df$Group))
  if (has_color) long_df$ColorGroup <- factor(long_df$ColorGroup)

  unique_clusters <- sort(unique(long_df$Cluster))
  plot_list <- list()

  for (ci in unique_clusters) {
    sub <- long_df[long_df$Cluster == ci, , drop = FALSE]
    if (nrow(sub) == 0) next
    n_feat <- length(unique(sub$Gene))

    if (has_color) {
      p <- ggplot2::ggplot(sub, ggplot2::aes(
             x = Group, y = Exp,
             colour = ColorGroup, linetype = ColorGroup,
             group = ColorGroup)) +
        ggplot2::stat_summary(fun = mean, geom = "line", linewidth = 1.3) +
        ggplot2::stat_summary(fun.data = ggplot2::mean_se,
                              geom = "errorbar", width = 0.3) +
        ggplot2::labs(colour = color_label, linetype = color_label)
    } else {
      p <- ggplot2::ggplot(sub, ggplot2::aes(
             x = Group, y = Exp, group = 1)) +
        ggplot2::stat_summary(fun = mean, geom = "line", linewidth = 1.3) +
        ggplot2::stat_summary(fun.data = ggplot2::mean_se,
                              geom = "errorbar", width = 0.3)
    }

    p <- p +
      ggplot2::labs(
        title = sprintf("Cluster %d (n=%d)", ci, n_feat),
        y = "Expression", x = x_label
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    plot_list[[as.character(ci)]] <- p
  }
  plot_list
}

#' Write a separate heatmap PNG per partition cluster
#'
#' For each cluster, subsets features and generates its own heatmap PNG.
#' Within each cluster, rows are hierarchically clustered for readability.
#'
#' @param expr_mat Feature x sample expression matrix
#' @param clusters Named integer vector (feature IDs -> cluster numbers)
#' @param annotation_col Column annotation data.frame for pheatmap
#' @param out_dir Output directory for heatmap PNGs
#' @param ... Additional arguments passed to plot_heatmap_core()
#' @return Character vector of written file paths
#' @export
save_per_cluster_heatmaps <- function(expr_mat, clusters, annotation_col, out_dir, ...) {
  expr_mat <- as.matrix(expr_mat)
  unique_clusters <- sort(unique(clusters))
  written <- character(0)

  for (ci in unique_clusters) {
    feats <- names(clusters[clusters == ci])
    valid_feats <- intersect(feats, rownames(expr_mat))
    if (length(valid_feats) == 0) next

    mat_sub <- expr_mat[valid_feats, , drop = FALSE]
    f_hm <- file.path(out_dir, sprintf("Partition_clustering_heatmap_cluster%d.png", ci))

    # pheatmap's row clustering needs >= 2 features; disable it for singletons.
    cluster_rows_sub <- nrow(mat_sub) >= 2

    p <- plot_heatmap_core(
      expr_mat       = mat_sub,
      annotation_col = annotation_col,
      title          = sprintf("Partition Cluster %d (%d features)", ci, length(valid_feats)),
      scale_rows     = TRUE,
      cluster_rows   = cluster_rows_sub,
      cluster_cols   = FALSE,
      max_rows       = NULL,
      ...
    )
    save_heatmap_to_file(p, f_hm)
    written <- c(written, f_hm)
  }
  written
}

# R/plots/plot_de.R

#' Volcano plot for a single DE table (one contrast)
#'
#' Consistent logic:
#' Y-axis: -log10(P.Value)  [Raw p-value]
#' Color:  Up (red), Down (blue), NS (grey) based on p-value and logFC thresholds
#' H-line: p_cutoff (matches the Y-axis and coloring threshold)
#'
#' @param de_tbl Data frame with logFC, P.Value, adj.P.Val
#' @param cfg Config list (sections de$p_cutoff, de$logfc_cutoff or de$linear_fc_cutoff)
#' @param title Plot title
#' @param pvalue_type Character; "padj" (default, use adj.P.Val) or "pval" (use raw P.Value).
#'   When "padj" is requested but adj.P.Val is absent, falls back to raw P.Value.
#' @param ... Ignored
#' @return ggplot object
plot_volcano <- function(de_tbl, cfg, title = NULL, pvalue_type = "padj", ...) {
  stopifnot(is.data.frame(de_tbl))
  pvalue_type <- match.arg(pvalue_type, c("padj", "pval"))

  # Required columns
  req_cols <- c("logFC", "P.Value")
  missing <- setdiff(req_cols, colnames(de_tbl))
  if (length(missing) > 0) {
    stop("plot_volcano: de_tbl missing columns: ", paste(missing, collapse = ", "))
  }

  # Thresholds
  p_cut <- cfg$de$p_cutoff %||% 0.05
  if (!is.null(cfg$de$logfc_cutoff)) {
    log2fc_cut <- cfg$de$logfc_cutoff
  } else {
    lin_fc_cut <- cfg$de$linear_fc_cutoff %||% 1.5
    log2fc_cut <- log2(lin_fc_cut)
  }

  df <- de_tbl

  # Prepare plotting data
  df$.logFC <- as.numeric(df[["logFC"]])

  # Y-axis: respect pvalue_type, fall back to raw P.Value when adj.P.Val absent
  p_col_vol <- if (pvalue_type == "padj" && "adj.P.Val" %in% colnames(df)) "adj.P.Val" else "P.Value"
  df$.pval <- as.numeric(df[[p_col_vol]])
  df$.pval_plot <- ifelse(is.na(df$.pval), 1, df$.pval)
  df$.neglog10p <- -log10(pmax(df$.pval_plot, 1e-300))

  # Categorize: Up, Down, NS (not significant)
  is_sig <- !is.na(df$.logFC) & (df$.pval_plot <= p_cut) & (abs(df$.logFC) >= log2fc_cut)
  is_up <- is_sig & (df$.logFC > 0)
  is_down <- is_sig & (df$.logFC < 0)

  df$.direction <- factor(
    ifelse(is_up, "Up", ifelse(is_down, "Down", "NS")),
    levels = c("NS", "Down", "Up")
  )

  # Count stats
  n_up <- sum(is_up, na.rm = TRUE)
  n_down <- sum(is_down, na.rm = TRUE)
  n_total <- n_up + n_down

  if (n_total == 0) {
    if (is.null(title)) title <- "Volcano (0 sig)"
    message(sprintf(
      "[plot_volcano] '%s': 0 genes passed (p<=%.2g, |logFC|>=%.2g). Total rows: %d",
      title, p_cut, log2fc_cut, nrow(df)
    ))
  } else {
    # Sort so significant points are plotted last (on top)
    df <- df[order(df$.direction), ]
    message(sprintf(
      "[plot_volcano] '%s': %d up, %d down (p<=%.2g, |logFC|>=%.2g)",
      title, n_up, n_down, p_cut, log2fc_cut
    ))
  }

  ggplot2::ggplot(df, ggplot2::aes(x = .logFC, y = .neglog10p)) +
    ggplot2::geom_point(ggplot2::aes(color = .direction, alpha = .direction), size = 1.5, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      name = "Regulation",
      values = c("NS" = "grey70", "Down" = "blue", "Up" = "red"),
      labels = c(
        "NS" = sprintf("NS (%d)", nrow(df) - n_total),
        "Down" = sprintf("Down (%d)", n_down),
        "Up" = sprintf("Up (%d)", n_up)
      )
    ) +
    ggplot2::scale_alpha_manual(values = c("NS" = 0.25, "Down" = 0.8, "Up" = 0.8), guide = "none") +
    ggplot2::geom_vline(xintercept = c(-log2fc_cut, log2fc_cut), linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "black") +
    ggplot2::labs(
      title = title %||% "Volcano plot",
      x = "log2(FC)",
      y = sprintf("-log10(%s)", if (p_col_vol == "adj.P.Val") "adj.P.Val" else "P.Value")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")
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

  # Use adjusted p-value when available, fall back to raw P.Value
  p_col <- if ("adj.P.Val" %in% colnames(de_tbl)) "adj.P.Val" else "P.Value"
  if (!(p_col %in% colnames(de_tbl))) stop("plot_ma: need 'P.Value' or 'adj.P.Val' column.")

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
  if (!is.null(cfg$de$logfc_cutoff)) {
    log2fc_cut <- cfg$de$logfc_cutoff
  } else {
    lin_fc_cut <- cfg$de$linear_fc_cutoff %||% 1.5
    log2fc_cut <- log2(lin_fc_cut)
  }

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
#'
#' Returns a Blues palette of length \code{n}.
#' \describe{
#'   \item{\code{reverse = TRUE} (default)}{dark-to-light. Use for distance
#'         heatmaps where SMALL values mean MORE similar (dark = similar).}
#'   \item{\code{reverse = FALSE}}{light-to-dark. Use for correlation heatmaps
#'         where LARGE values mean MORE similar (dark = similar, diagonal = darkest).}
#' }
get_heatmap_colors <- function(n = 255, reverse = TRUE) {
  pal <- RColorBrewer::brewer.pal(9, "Blues")
  if (reverse) pal <- rev(pal)
  grDevices::colorRampPalette(pal)(n)
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
  
  annot_col <- data.frame(
    Condition = meta[[color_col]],
    row.names = meta[[sample_col]]
  )
  annot_col <- annot_col[colnames(mat2plot), , drop = FALSE]
  
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
