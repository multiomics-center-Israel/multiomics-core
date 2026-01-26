mod_rnaseq_qc_pre <- function(pre, config, out_dir) {
  assert_pre_contract(pre, stage = "rna")
  
  dirs <- create_legacy_output_dirs(out_dir, create = TRUE)
  
  mat <- pre$expr_work
  meta <- pre$meta
  cfg <- config$modes$rna
  eff_color <- cfg$effects$color %||% "Group"
  eff_shape <- cfg$effects$shape %||% NULL
  
  # --- Sample correlation heatmap ---
  cor_mat <- stats::cor(mat, use = "pairwise.complete.obs", method = "pearson")
  cor_png <- file.path(dirs$diagnostic_plots, "sample_correlation_heatmap.png")
  grDevices::png(cor_png, width = 1200, height = 1000, res = 150)
  pheatmap::pheatmap(
    cor_mat,
    main = "Sample correlation",
    clustering_method = "complete"
  )
  grDevices::dev.off()
  
  # --- PCA plot ---
  pca <- stats::prcomp(t(mat), scale. = TRUE)
  var <- (pca$sdev^2) / sum(pca$sdev^2)
  df <- data.frame(
    SampleID = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  df <- dplyr::left_join(df, meta, by = c("SampleID" = (cfg$id_columns$sample_col %||% "SampleID")))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = .data[[eff_color]])) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = "PCA (expr_work)",
      x = sprintf("PC1 (%.1f%%)", 100 * var[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * var[2])
    ) +
    ggplot2::theme_bw()
  
  if (!is.null(eff_shape) && eff_shape %in% names(df)) {
    p <- p + ggplot2::aes(shape = .data[[eff_shape]])
  }
  
  pca_png <- file.path(dirs$diagnostic_plots, "pca_PC1_PC2.png")
  ggplot2::ggsave(pca_png, p, width = 8, height = 6, dpi = 150)
  
  list(
    files = c(cor_png, pca_png),
    pca = pca,
    cor = cor_mat
  )
}
