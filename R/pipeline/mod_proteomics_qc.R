#' Proteomics QC module (legacy-style outputs)
#'
#' Generates all diagnostic QC plots for proteomics and writes them
#' to <run_dir>/Diagnostic_plots/, matching legacy naming conventions.
#'
#' @param pre List returned by preprocess_proteomics()
#' @param config Full config list
#' @param run_dir Output root directory
#' @return list of 1.Character vector of written file paths and 2.list of plots objects
mod_proteomics_qc <- function(pre, config, run_dir) {
  assert_pre_contract(pre, stage = "proteomics")
  
  dirs   <- create_legacy_output_dirs(run_dir)
  out_qc <- dirs$diagnostic_plots
  cfg    <- config$modes$proteomics
  
  plots <- list()        
  files <- character(0)  
  
  ## ---------- PCA ----------
  f_pca12 <- file.path(out_qc, "PCA_PC1.vs.PC2.png")
  p12 <- qc_pca_scatter(pre$expr_imp_single, pre$meta, cfg, pcs = c(1, 2), out_file = f_pca12)
  files <- c(files, f_pca12)
  plots$pca_1_2 <- p12   # 
  
  f_pca13 <- file.path(out_qc, "PCA_PC1.vs.PC3.png")
  p13 <- qc_pca_scatter(pre$expr_imp_single, pre$meta, cfg, pcs = c(1, 3), out_file = f_pca13)
  files <- c(files, f_pca13)
  plots$pca_1_3 <- p13
  
  f_pca3d <- file.path(out_qc, "PCA_3D.html")
  p_3d <- qc_pca_3d(pre$expr_imp_single, pre$meta, cfg, out_file = f_pca3d)
  files <- c(files, f_pca3d)
  if (!is.null(p_3d)) plots$pca_3d <- p_3d
  
  ## ---------- Density ----------
  f_hist <- file.path(out_qc, "protein_histograms_summary.png")
  p_dens <- qc_proteomics_density(pre$expr_imp_single, pre$meta, cfg, out_file = f_hist)
  files <- c(files, f_hist)
  plots$density <- p_dens
  
  ## ---------- Imputation QC ----------
  if (!is.null(pre$imputation) && !is.null(pre$imputation$imputed_flag)) {

    f_imp_hist <- file.path(out_qc, "imputed_histograms_samples_summary.png")
    p_imp <- plot_imputation_summary(
      pre$expr_imp_single,
      pre$imputation$imputed_flag
    )
    ggplot2::ggsave(f_imp_hist, plot = p_imp, width = 12, height = 8, dpi = 150)
    
    files <- c(files, f_imp_hist)
    plots$imputation_hist <- p_imp 
    
    f_box <- file.path(out_qc, "imputed_boxplot.png")
    p_bp <- norm_boxplot(pre$expr_imp_single, pre$meta, cfg, out_file = f_box)
    files <- c(files, f_box)
    plots$boxplot <- p_bp
  }
  
  ## ---------- Sample distance ----------
  annot <- data.frame(
    Condition = pre$meta[[cfg$effects$color]],
    row.names = pre$meta[[cfg$effects$samples]]
  )
  
  f_dist <- file.path(out_qc, "sample_distance_heatmap.png")
  ph <- plot_sample_distance_heatmap(pre$expr_imp_single, annotation_col = annot)
  written <- save_heatmap_to_file(ph, f_dist, width = 1600, height = 1200, res = 150)
  files <- c(files, written)
  plots$dist_heatmap <- ph # pheatmap object
  
  
  f_dist_na <- file.path(out_qc, "sample_distance_heatmap_NA.png")
  ph_na <- plot_sample_distance_heatmap(pre$expr_raw, annotation_col = annot)
  written_na <- save_heatmap_to_file(ph_na, f_dist_na, width = 1600, height = 1200, res = 150)
  files <- c(files, written_na)
  plots$dist_heatmap_na <- ph_na
  
  ## ---------- Correlation ----------
  f_cor <- file.path(out_qc, "sample_correlation_heatmap.png")
  p_cor <- qc_sample_correlation_heatmap(pre$expr_imp_single, pre$meta, cfg, out_file = f_cor)
  files <- c(files, f_cor)
  plots$correlation <- p_cor
  
  ## ---------- Expression heatmaps ----------
  f_hm <- file.path(out_qc, "samples_protein_heatmap.png")
  p_hm <- wrap_qc_heatmap(pre$expr_imp_single, pre$meta, cfg, out_file = f_hm)
  files <- c(files, f_hm)
  plots$heatmap_clusters <- p_hm
  
  f_hm_nocol <- file.path(out_qc, "samples_protein_heatmap_wo_col.png")
  p_hm_nocol <- wrap_qc_heatmap(pre$expr_imp_single, pre$meta, cfg, out_file = f_hm_nocol)
  files <- c(files, f_hm_nocol)
  plots$heatmap_nocol <- p_hm_nocol
  
  return(list(
    plots = plots,  
    files = files  
  ))
}

