#' Proteomics QC module (legacy-style outputs)
#'
#' Generates all diagnostic QC plots for proteomics and writes them
#' to <run_dir>/Diagnostic_plots/, matching legacy naming conventions.
#'
#' @param pre List returned by preprocess_proteomics()
#' @param config Full config list
#' @param run_dir Output root directory
#' @return Character vector of written file paths
mod_proteomics_qc <- function(pre, config, run_dir) {
  
  dirs   <- create_legacy_output_dirs(run_dir)
  out_qc <- dirs$diagnostic_plots
  cfg    <- config$modes$proteomics
  
  files <- character(0)
  
  ## ---------- PCA ----------
  f_pca12 <- file.path(out_qc, "PCA_PC1.vs.PC2.png")
  qc_pca_scatter(pre$expr_imp, pre$meta, cfg, pcs = c(1, 2), out_file = f_pca12)
  files <- c(files, f_pca12)
  
  f_pca13 <- file.path(out_qc, "PCA_PC1.vs.PC3.png")
  qc_pca_scatter(pre$expr_imp, pre$meta, cfg, pcs = c(1, 3), out_file = f_pca13)
  files <- c(files, f_pca13)
  
  f_pca3d <- file.path(out_qc, "PCA_3D.html")
  qc_pca_3d(pre$expr_imp, pre$meta, cfg, out_file = f_pca3d)
  files <- c(files, f_pca3d)
  
  ## ---------- Density ----------
  f_hist <- file.path(out_qc, "protein_histograms_summary.png")
  qc_proteomics_density(pre$expr_imp, pre$meta, cfg, out_file = f_hist)
  files <- c(files, f_hist)
  
  ## ---------- Imputation QC ----------
  if (!is.null(pre$imputation$imputed_flag)) {
    
    f_imp_hist <- file.path(out_qc, "imputed_histograms_samples_summary.png")
    p_imp <- build_imputed_histograms_summary(
      pre$expr_imp,
      pre$imputation$imputed_flag,
      cfg
    )
    with_png(f_imp_hist, { print(p_imp) })
    files <- c(files, f_imp_hist)
    
    f_box <- file.path(out_qc, "imputed_boxplot.png")
    norm_boxplot(pre$expr_imp, pre$meta, cfg, out_file = f_box)
    files <- c(files, f_box)
  }
  
  ## ---------- Sample distance ----------
  annot <- data.frame(
    Condition = pre$meta[[cfg$effects$color]],
    row.names = pre$meta[[cfg$effects$samples]]
  )
  
  f_dist <- file.path(out_qc, "sample_distance_heatmap.png")
  plot_sample_distance_heatmap(pre$expr_imp, annotation_col = annot, out_file = f_dist)
  files <- c(files, f_dist)
  
  f_dist_na <- file.path(out_qc, "sample_distance_heatmap_NA.png")
  plot_sample_distance_heatmap(pre$expr_imp, annotation_col = NULL, out_file = f_dist_na)
  files <- c(files, f_dist_na)
  
  ## ---------- Correlation ----------
  f_cor <- file.path(out_qc, "sample_correlation_heatmap.png")
  qc_sample_correlation_heatmap(pre$expr_imp, pre$meta, cfg, out_file = f_cor)
  files <- c(files, f_cor)
  
  ## ---------- Expression heatmaps ----------
  f_hm <- file.path(out_qc, "samples_protein_heatmap.png")
  qc_expr_heatmap(pre$expr_imp, pre$meta, cfg, out_file = f_hm)
  files <- c(files, f_hm)
  
  f_hm_nocol <- file.path(out_qc, "samples_protein_heatmap_wo_col.png")
  qc_expr_heatmap(pre$expr_imp, pre$meta, cfg, out_file = f_hm_nocol, cluster_cols = FALSE)
  files <- c(files, f_hm_nocol)
  
  files
}
