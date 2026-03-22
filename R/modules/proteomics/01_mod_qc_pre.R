#' Proteomics pre-DE QC module
#'
#' QC plots that depend only on preprocessing outputs (no DE needed).
#' @param pre List returned by preprocess_proteomics()
#' @param config Full config list
#' @param out_dir Output root directory
#' @return list(plots, files)
mod_proteomics_qc_pre <- function(pre, config, out_dir) {
    stage <- "proteomics"
    stopifnot(is.character(out_dir), length(out_dir) == 1)
    assert_pre_contract(pre, stage = stage)

    dirs <- create_legacy_output_dirs(out_dir)
    out_qc <- dirs$diagnostic_plots
    cfg <- config$modes$proteomics

    plots <- list()
    files <- character(0)

    # ---------- PCA ----------
    f_pca12 <- file.path(out_qc, "PCA_PC1.vs.PC2.png")
    p12 <- qc_pca_scatter(pre$expr_imp_single, pre$meta, cfg, pcs = c(1, 2), out_file = f_pca12)
    files <- c(files, f_pca12)
    plots$pca_1_2 <- p12

    f_pca13 <- file.path(out_qc, "PCA_PC1.vs.PC3.png")
    p13 <- qc_pca_scatter(pre$expr_imp_single, pre$meta, cfg, pcs = c(1, 3), out_file = f_pca13)
    files <- c(files, f_pca13)
    plots$pca_1_3 <- p13

    # ---------- PCA with top variable proteins (for report dropdown) ----------
    n_top_values <- c(500, 1000, 2000)
    n_features <- nrow(pre$expr_imp_single)
    cfg_temp <- cfg

    for (n_top in n_top_values) {
        if (n_top <= n_features) {
            prot_vars <- apply(pre$expr_imp_single, 1, var, na.rm = TRUE)
            top_idx <- order(prot_vars, decreasing = TRUE)[1:n_top]
            mat_top <- pre$expr_imp_single[top_idx, , drop = FALSE]

            f_pca_top <- file.path(out_qc, sprintf("PCA_top%d.png", n_top))
            tryCatch({
                p_top <- qc_pca_scatter(mat_top, pre$meta, cfg_temp, pcs = c(1, 2), out_file = f_pca_top)
                files <- c(files, f_pca_top)
                plots[[sprintf("pca_top%d", n_top)]] <- p_top
                message(sprintf("  Generated PCA with top %d variable proteins", n_top))
            }, error = function(e) {
                message(sprintf("  Could not generate PCA with top %d proteins: %s", n_top, e$message))
            })
        }
    }

    # ---------- PCA Subsets ----------
    pca_subsets <- cfg$qc$pca_subsets
    if (!is.null(pca_subsets) && length(pca_subsets) > 0) {
        group_col_name <- cfg$effects$color %||% "Condition"
        sample_col_name <- cfg$effects$samples %||% "SampleID"

        for (si in seq_along(pca_subsets)) {
            subset_def <- pca_subsets[[si]]
            subset_name <- subset_def$name %||% paste0("subset_", si)
            conditions <- subset_def$conditions

            if (is.null(conditions)) next  # null = all samples, already covered above

            # Filter to matching conditions
            keep_samples <- pre$meta[[group_col_name]] %in% conditions
            if (sum(keep_samples) < 3) {
                message(sprintf("  PCA subset '%s': fewer than 3 samples, skipping.", subset_name))
                next
            }

            sample_ids <- pre$meta[[sample_col_name]][keep_samples]
            mat_subset <- pre$expr_imp_single[, colnames(pre$expr_imp_single) %in% sample_ids, drop = FALSE]
            meta_subset <- pre$meta[keep_samples, , drop = FALSE]

            safe_name <- gsub("[^a-zA-Z0-9_]", "_", subset_name)
            f_pca_sub <- file.path(out_qc, sprintf("PCA_%s.png", safe_name))
            tryCatch({
                p_sub <- qc_pca_scatter(mat_subset, meta_subset, cfg, pcs = c(1, 2), out_file = f_pca_sub)
                files <- c(files, f_pca_sub)
                plots[[sprintf("pca_subset_%s", safe_name)]] <- p_sub
                message(sprintf("  Generated PCA subset: %s (%d samples)", subset_name, ncol(mat_subset)))
            }, error = function(e) {
                message(sprintf("  Could not generate PCA subset '%s': %s", subset_name, e$message))
            })
        }
    }

    f_pca3d <- file.path(out_qc, "PCA_3D.html")
    p_3d <- qc_pca_3d(pre$expr_imp_single, pre$meta, cfg, out_file = f_pca3d)
    files <- c(files, f_pca3d)
    if (!is.null(p_3d)) plots$pca_3d <- p_3d

    # ---------- Density ----------
    f_hist <- file.path(out_qc, "protein_histograms_summary.png")
    p_dens <- qc_omic_density(pre$expr_imp_single, pre$meta, cfg, out_file = f_hist)
    files <- c(files, f_hist)
    plots$density <- p_dens

    # ---------- Pre-imputation boxplot ----------
    if (!is.null(pre$expr_filt_pre_imp)) {
        f_box_pre <- file.path(out_qc, "pre_imputation_boxplot.png")
        p_bp_pre <- norm_boxplot(pre$expr_filt_pre_imp, pre$meta, cfg,
                                  out_file = f_box_pre,
                                  title = "Expression boxplots (before imputation)")
        files <- c(files, f_box_pre)
        plots$boxplot_pre_imp <- p_bp_pre
    }

    # ---------- Imputation QC ----------
    if (!is.null(pre$imputation_qc) && !is.null(pre$imputation_qc$imputed_flag)) {
        imp_w <- cfg$imputation$width %||% "NA"
        imp_s <- cfg$imputation$downshift %||% "NA"
        f_imp_hist <- file.path(out_qc, sprintf("imputed_histograms_samples_summary_w%s_s%s.png", imp_w, imp_s))
        p_imp <- qc_imputation_summary(
            imputed      = pre$expr_imp_single,
            imputed_flag = pre$imputation_qc$imputed_flag,
            cfg          = cfg,
            out_file     = f_imp_hist
        )
        files <- c(files, f_imp_hist)
        plots$imputation_hist <- p_imp

        f_box <- file.path(out_qc, "imputed_boxplot.png")
        p_bp <- norm_boxplot(pre$expr_imp_single, pre$meta, cfg, out_file = f_box)
        files <- c(files, f_box)
        plots$boxplot <- p_bp
    }

    # ---------- Sample distance ----------
    annot <- data.frame(
        Condition = pre$meta[[cfg$effects$color]],
        row.names = pre$meta[[cfg$effects$samples]]
    )

    f_dist <- file.path(out_qc, "sample_distance_heatmap.png")
    ph <- plot_sample_distance_heatmap(pre$expr_imp_single, annotation_col = annot)
    written <- save_heatmap_to_file(ph, f_dist, width = 1600, height = 1200, res = 150)
    files <- c(files, written)
    plots$dist_heatmap <- ph

    f_dist_na <- file.path(out_qc, "sample_distance_heatmap_NA.png")
    ph_na <- plot_sample_distance_heatmap(pre$expr_raw, annotation_col = annot)
    written_na <- save_heatmap_to_file(ph_na, f_dist_na, width = 1600, height = 1200, res = 150)
    files <- c(files, written_na)
    plots$dist_heatmap_na <- ph_na

    # ---------- Correlation ----------
    f_cor <- file.path(out_qc, "sample_correlation_heatmap.png")
    p_cor <- qc_sample_correlation_heatmap(pre$expr_imp_single, pre$meta, cfg, out_file = f_cor)
    files <- c(files, f_cor)
    plots$correlation <- p_cor

    # ---------- Expression heatmaps ----------
    f_hm <- file.path(out_qc, "samples_protein_heatmap.png")
    hm_clusters <- wrap_qc_heatmap(pre$expr_imp_single, pre$meta, cfg, stage = stage, out_file = f_hm)
    files <- c(files, f_hm)
    plots$heatmap_clusters <- hm_clusters

    f_hm_nocol <- file.path(out_qc, "samples_protein_heatmap_wo_col.png")
    hm_nocol <- wrap_qc_heatmap(pre$expr_imp_single, pre$meta, cfg, stage = stage, out_file = f_hm_nocol, cluster_cols = FALSE)
    files <- c(files, f_hm_nocol)
    plots$heatmap_nocol <- hm_nocol

    # ---------- Extract PCA objects from plot attributes ----------
    # PCA plot p12 has attributes attached by qc_pca_scatter
    pca_obj <- attr(p12, "pca_result")
    scores <- attr(p12, "scores")
    var_expl <- attr(p12, "var_expl")

    # Get color and shape column names from config
    color <- cfg$effects$color %||% NULL
    shape <- cfg$effects$shape %||% NULL

    # ---------- UMAP ----------
    umap_res <- NULL
    if (isTRUE(cfg$qc$run_umap)) {
        umap_res <- tryCatch(
            run_proteomics_umap(pre$expr_imp_single, pre$meta, cfg),
            error = function(e) {
                message("UMAP generation failed: ", e$message)
                NULL
            }
        )
        if (!is.null(umap_res) && !is.null(umap_res$plot)) {
            f_umap <- file.path(out_qc, "UMAP.png")
            ggplot2::ggsave(f_umap, plot = umap_res$plot, width = 8, height = 6, dpi = 150)
            files <- c(files, f_umap)
            plots$umap <- umap_res$plot
        }
    }

    # ---------- Outlier detection ----------
    outlier_res <- NULL
    cor_mat <- cor(pre$expr_imp_single, use = "pairwise.complete.obs")
    outlier_res <- detect_proteomics_outliers(pre$expr_imp_single, pca_obj, cor_mat, cfg)
    if (!is.null(outlier_res) && nrow(outlier_res$flagged_samples_df) > 0) {
        f_outlier <- file.path(out_qc, "outlier_report.csv")
        write.csv(outlier_res$flagged_samples_df, f_outlier, row.names = FALSE)
        files <- c(files, f_outlier)
    }

    # Return plots, files, AND objects for Shiny export
    list(
        plots = plots,
        files = unique(files),
        objects = list(
            norm_log_counts_pca = pca_obj, # prcomp result
            pca_scores = assert_pca_scores(scores, context = "proteomics QC"), # data.frame with PCs + metadata
            var_expl = var_expl, # variance explained
            color = color, # string: column name for color aesthetic
            shape = shape, # string: column name for shape aesthetic
            umap_res = umap_res, # UMAP coordinates + plot
            outlier_res = outlier_res # outlier detection results
        )
    )
}
