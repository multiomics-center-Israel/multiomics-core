#' RNA-seq pre-DE QC module
#'
#' QC plots that depend only on preprocessing outputs.
#' Mirrors Proteomics logic using core plotting functions.
#'
#' @param pre List returned by preprocess_rna()
#' @param config Full config list
#' @param out_dir Output root directory
#' @return list(plots, files, objects, ...)
mod_rnaseq_qc_pre <- function(pre, config, out_dir) {
    assert_pre_contract(pre, stage = "rna")

    dirs <- create_legacy_output_dirs(out_dir, create = TRUE)
    out_qc <- dirs$diagnostic_plots
    cfg <- config$modes$rna

    plots <- list()
    files <- character(0)

    # Use working expression matrix (normalized log-counts) for QC
    mat <- pre$expr_work
    meta <- pre$meta

    # ---------- PCA ----------
    f_pca12 <- file.path(out_qc, "PCA_PC1.vs.PC2.png")
    p12 <- qc_pca_scatter(mat, meta, cfg, pcs = c(1, 2), out_file = f_pca12)
    files <- c(files, f_pca12)
    plots$pca_1_2 <- p12

    # PCA 1 vs 3
    f_pca13 <- file.path(out_qc, "PCA_PC1.vs.PC3.png")
    p13 <- qc_pca_scatter(mat, meta, cfg, pcs = c(1, 3), out_file = f_pca13)
    files <- c(files, f_pca13)
    plots$pca_1_3 <- p13

    # PCA 3D
    f_pca3d <- file.path(out_qc, "PCA_3D.html")
    p_3d <- qc_pca_3d(mat, meta, cfg, out_file = f_pca3d)
    files <- c(files, f_pca3d)
    if (!is.null(p_3d)) plots$pca_3d <- p_3d

    # ---------- Density ----------
    f_hist <- file.path(out_qc, "rna_histograms_summary.png")
    p_dens <- qc_proteomics_density(
        mat, meta, cfg,
        out_file = f_hist,
        title = "Density plot of normalized counts"
    )
    files <- c(files, f_hist)
    plots$density <- p_dens

    # ---------- Boxplots ----------
    f_box <- file.path(out_qc, "norm_boxplot.png")
    p_bp <- norm_boxplot(mat, meta, cfg, out_file = f_box)
    files <- c(files, f_box)
    plots$boxplot <- p_bp

    # ---------- Sample distance heatmap ----------
    # Need to construct annotation_col manually or use helper?
    # qc_sample_distance_heatmap handles prepare_qc_data internally.
    f_dist <- file.path(out_qc, "sample_distance_heatmap.png")
    ph_dist <- qc_sample_distance_heatmap(mat, meta, cfg, out_file = f_dist)
    # Note: qc_sample_distance_heatmap doesn't return file path, just writes it.
    files <- c(files, f_dist)
    plots$dist_heatmap <- ph_dist

    # ---------- Correlation heatmap ----------
    f_cor <- file.path(out_qc, "sample_correlation_heatmap.png")
    p_cor <- qc_sample_correlation_heatmap(mat, meta, cfg, out_file = f_cor)
    files <- c(files, f_cor)
    plots$correlation <- p_cor

    # ---------- Expression heatmaps ----------
    f_hm <- file.path(out_qc, "samples_rna_heatmap.png")
    p_hm <- wrap_qc_heatmap(mat, meta, cfg, out_file = f_hm)
    files <- c(files, f_hm)
    plots$heatmap_clusters <- p_hm

    # heatmap without columns clustering? (Proteomics calls it keys 'heatmap_nocol')
    f_hm_nocol <- file.path(out_qc, "samples_rna_heatmap_wo_col.png")
    # wrap_qc_heatmap defaults to clustering both.
    # To disable col clustering, we need access to plot_heatmap_core args via ...
    # But wrap_qc_heatmap doesn't expose ... in its arg list (based on my read).
    # It hardcodes cluster_cols=TRUE.
    # So I'll skip wo_col for now unless I refactor wrap_qc_heatmap,
    # OR I just call plot_heatmap_core directly.
    # Calling plot_heatmap_core directly is safer.
    annot <- data.frame(
        Condition = meta[[cfg$effects$color]],
        row.names = meta[[cfg$effects$samples]]
    )
    p_hm_nocol <- plot_heatmap_core(
        expr_mat = mat,
        annotation_col = annot,
        title = "QC: Sample RNA Expression",
        max_rows = 2000,
        cluster_rows = TRUE,
        cluster_cols = FALSE # The specific difference
    )
    save_heatmap_to_file(p_hm_nocol, f_hm_nocol)
    files <- c(files, f_hm_nocol)
    plots$heatmap_nocol <- p_hm_nocol

    # ---------- Extract PCA objects from plot attributes (p12) ----------
    pca_obj <- attr(p12, "pca_result")
    scores <- attr(p12, "scores")
    var_expl <- attr(p12, "var_expl")

    eff_color <- cfg$effects$color %||% NULL
    eff_shape <- cfg$effects$shape %||% NULL

    # objects list (Proteomics style)
    objs <- list(
        norm_log_counts_pca = pca_obj,
        mat2plot = scores,
        var_expl = var_expl,
        color = eff_color,
        shape = eff_shape
    )

    list(
        files = unique(files),
        plots = plots,
        objects = objs,
        # Keep these at top level for compatibility with current build_data_to_shiny_legacy_rna
        # (until that function is updated to look inside 'objects')
        norm_log_counts_pca = pca_obj,
        mat2plot = scores
    )
}
