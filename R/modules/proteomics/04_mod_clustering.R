#' Proteomics clustering module (legacy-like outputs under out_dir/Clustering)
#'
#' Runs hierarchical clustering on DE features (always, if enabled),
#' and additionally runs partition clustering + binary patterns
#' only when the data has enough groups (>= min_groups), where groups are
#' defined by cfg$effects$color.
#'
#' @param pre preprocessed proteomics object (expects $expr_imp_single and $meta)
#' @param de_res  DE results list (expects $summary_df at least)
#' @param config  full YAML config list
#' @param out_dir output run directory (project folder)
#'
#' @return list(plots, files)
mod_proteomics_clustering <- function(pre, de_res, config, out_dir) {
    stopifnot(is.character(out_dir), length(out_dir) == 1)

    assert_pre_contract(pre, stage = "proteomics")
    assert_de_contract(de_res, stage = "proteomics")

    cfg <- config$modes$proteomics
    cl <- cfg$clustering

    excel_order <- NULL

    # Objects for legacy Shiny export
    pheatmap_payload <- NULL
    patterns_tbl <- NULL
    patterns_list <- NULL
    heatmaps_by_pattern <- NULL
    clusters_vec <- NULL
    excel_order <- NULL
    written <- character(0)
    plots <- list()

    # If clustering is missing/disabled -> no-op
    if (is.null(cl) || isFALSE(cl$enabled)) {
        message("Clustering disabled. Skipping.")
        return(list(plots = list(), files = character(0))) # Return empty list structure
    }

    # ---- output root (UNDER PROJECT RUN DIR) ----
    clustering_dir <- file.path(out_dir, "Clustering")
    ensure_dir(clustering_dir)

    # ---- decide which steps to run ----
    flags <- clustering_run_flags(pre, cfg)
    message(sprintf(
        "Clustering flags: hierarchical=%s, partition=%s, binary=%s",
        flags$hierarchical, flags$partition, flags$binary_patterns
    ))

    # ---- build annotation_col for heatmaps using effects ----
    eff_color <- get_color_config(cfg)
    annot <- NULL
    if (!is.null(eff_color) && eff_color %in% colnames(pre$meta)) {
        annot <- data.frame(
            Condition = pre$meta[[eff_color]],
            row.names = pre$meta[[cfg$effects$samples]]
        )
    }

    # Get DE features
    de_features <- get_de_features(de_res, cfg)

    if (length(de_features) < 2) {
        message(sprintf("Clustering: only %d DE feature(s) found; need >= 2. Skipping.", length(de_features)))
        return(list(
            plots = list(), files = character(0), excel_order = NULL,
            objects = list(hm_hier_de = NULL, patterns = NULL, patterns_list = NULL,
                           heatmaps_by_pattern = NULL, clusters = NULL)
        ))
    }

    # Expression matrix (Imputed)
    expr_mat <- as.matrix(pre$expr_imp_single)

    written <- character(0)
    plots <- list()

    # ------ 1) Hierarchical clustering ---------
    if (isTRUE(flags$hierarchical)) {
        hcfg <- cl$steps$hierarchical %||% list()

        clust_out_dir <- file.path(clustering_dir, "Hierarchical")
        ensure_dir(clust_out_dir)

        # Run clustering
        hc_res <- run_clustering(
            expr_mat = expr_mat,
            col_data = pre$meta,
            de_features = de_features,
            config = list(
                method   = "hierarchical",
                distance = hcfg$distance %||% "euclidean",
                linkage  = hcfg$linkage %||% "complete",
                k        = hcfg$k %||% NULL
            )
        )

        mat_de <- as.matrix(pre$expr_imp_single)[de_features, , drop = FALSE]
        z_de <- zscore_rows(mat_de)
        colnames(z_de) <- paste0(colnames(z_de), ".zscore")

        # Pre-compute ordered matrix for Shiny (Professional approach)
        # Order rows according to clustering result - makes Shiny app robust
        ordered_row_ids <- intersect(hc_res$ordering, rownames(z_de))
        z_de_ordered <- z_de[ordered_row_ids, , drop = FALSE]

        excel_order <- list(
            ordered_ids = hc_res$ordering,
            zscore_mat  = z_de
        )

        # Build DE pattern row annotations (up/down per contrast)
        prot_de_cfg <- cfg$de %||% list()
        prot_p_cutoff <- prot_de_cfg$p_cutoff %||% 0.05
        prot_lin_fc <- prot_de_cfg$linear_fc_cutoff %||% 1.5
        prot_log2fc <- log2(prot_lin_fc)

        # Get ID column from config (proteomics may use different column name)
        prot_id_col <- cfg$de_table$id_col %||% "FeatureID"

        annot_context <- list(
            summary_df    = de_res$summary_df,
            p_cutoff      = prot_p_cutoff,
            log2fc_cutoff = prot_log2fc,
            id_col        = prot_id_col
        )

        # Heatmap setup
        f_hm <- file.path(clust_out_dir, "Hierarchical_DE_heatmap.png")

        # Run wrapper
        p_cluster <- wrap_clustering_heatmap(
            expr_mat = pre$expr_imp_single,
            meta = pre$meta,
            cfg = cfg,
            feature_ids = de_features,
            ordering = hc_res$ordering,
            annotation_row_builder = TRUE,
            annotation_row_context = annot_context,
            out_file = f_hm,
            title = sprintf("Hierarchical Clustering (%d DE features)", length(de_features))
        )
        written <- c(written, f_hm)
        plots$p_cluster_hier <- p_cluster

        # Capture pheatmap payload for Shiny (Professional pre-compute approach)
        # Store matrix in clustered order so Shiny doesn't need to extract from pheatmap
        pheatmap_payload <- list(
            pheatmap = p_cluster,
            mat = z_de_ordered, # Already in clustered order
            row_order = rownames(z_de_ordered), # Ordered row names for Plotly
            col_order = colnames(z_de_ordered), # Column names
            annotation_col = annot,
            feature_ids = de_features,
            is_zscored = TRUE,
            cluster_cols = FALSE,
            tree_row = hc_res$details # hclust object for dendrogram
        )

        # Save cluster assignments
        if (!is.null(hc_res$clusters)) {
            f_tbl <- file.path(clust_out_dir, "Hierarchical_clusters.tsv")
            cl_tbl <- build_clustering_output_table(hc_res$clusters, f_tbl)
            written <- c(written, f_tbl)
            plots$cl_tbl <- cl_tbl

            # Map for Shiny
            clusters_vec <- hc_res$clusters
        }
    }

    # ----- 2) Partition clustering (kmeans / PAM) + heatmap per cluster  ----

    if (isTRUE(flags$partition)) {
        part_dir_name <- "Partition_clustering"
        part_base_dir <- file.path(clustering_dir, part_dir_name)
        ensure_dir(part_base_dir)

        # Fit clusters on group means (using the effects function)
        part_res <- perform_partition_clustering_effects(
            expr_mat = pre$expr_imp_single,
            meta = pre$meta,
            cfg = cfg,
            de_features = de_features
        )

        # Capture clusters for Shiny
        clusters_vec <- part_res$clusters %||% NULL

        # Final output dir includes number of clusters (legacy style)
        part_dir <- file.path(part_base_dir, sprintf("Partition_clustering_%d_clusters", part_res$k))
        ensure_dir(part_dir)

        # (1) write clusters table
        f_tbl <- file.path(part_dir, "partition_clusters.tsv")
        clusters_tbl <- build_clustering_output_table(part_res$clusters, f_tbl)
        written <- c(written, f_tbl)
        plots$pt_tbl <- clusters_tbl

        # (2) heatmap
        feats <- names(part_res$clusters)
        mat <- as.matrix(pre$expr_imp_single)[feats, , drop = FALSE]

        # Order rows: cluster then name
        ord <- order(part_res$clusters, names(part_res$clusters))
        mat_ord <- mat[ord, , drop = FALSE]

        f_hm <- file.path(part_dir, "Partition_clustering_heatmap.png")

        # Using Core Plotter directly for custom ordering
        p_part <- plot_heatmap_core(
            expr_mat       = mat_ord,
            annotation_col = annot,
            title          = sprintf("Partition clustering (k=%d) on DE features (n=%d)", part_res$k, nrow(mat_ord)),
            scale_rows     = TRUE,
            cluster_rows   = FALSE,
            cluster_cols   = TRUE,
            max_rows       = NULL
        )

        save_heatmap_to_file(p_part, f_hm)
        plots$partition_heatmap <- p_part
        written <- c(written, f_hm)

        # (3) cluster profiles pdf
        prof <- build_cluster_profiles(part_res$z_group_means, part_res$clusters, part_res$k)

        if (!is.null(prof)) {
            f_pdf <- file.path(part_dir, "cluster_profiles.pdf")
            p_prof <- plot_cluster_profiles(prof, x_label = eff_color)

            n_clusters <- length(unique(prof$cluster))
            calc_height <- max(6, ceiling(n_clusters / 2) * 3)

            ggplot2::ggsave(f_pdf, plot = p_prof, width = 10, height = calc_height)

            written <- c(written, f_pdf)
            plots$cluster_profiles <- p_prof
        }

        # --- FIX: Export Legacy Data (Moved INSIDE the IF block) ---
        # This must be here because part_res and part_dir are only defined here.
        legacy_files <- write_clustering_legacy_profiles(
            expr_mat = pre$expr_imp_single, # Source of Truth (Absolute values)
            meta     = pre$meta, # Metadata
            clusters = part_res$clusters,
            cfg      = cfg,
            out_dir  = part_dir
        )

        written <- c(written, legacy_files)
    }

    # ---- 3) Binary patterns (only meaningful when >= 3 conditions) ----
    if (isTRUE(flags$binary_patterns)) {
        bcfg <- cl$steps$binary_patterns %||% list()
        clust_out_dir <- file.path(clustering_dir, "Binary_patterns")
        ensure_dir(clust_out_dir)

        # DE cutoffs for row annotations
        de_cfg <- cfg$de %||% list()
        bp_p_cutoff <- de_cfg$p_cutoff %||% 0.05
        bp_lin_fc_cutoff <- de_cfg$linear_fc_cutoff %||% 1.5
        bp_log2fc_cutoff <- log2(bp_lin_fc_cutoff)

        # Get ID column from config (proteomics may use different column name)
        bp_id_col <- cfg$de_table$id_col %||% "FeatureID"

        # Run function and capture result
        # For proteomics: both matrices are log2-transformed (no separate counts matrix)
        # expr_mat_counts = NULL means it will use expr_mat_corr for gating too
        bp_res <- run_binary_patterns(
            expr_mat_corr      = expr_mat,
            expr_mat_counts    = NULL, # proteomics has no separate counts matrix
            meta               = pre$meta,
            cfg                = cfg,
            de_features        = de_features,
            out_dir            = clust_out_dir,
            summary_df         = de_res$summary_df,
            corr_cutoff        = bcfg$corr_cutoff %||% 0.8,
            counts_cutoff_high = bcfg$counts_cutoff_high %||% bcfg$counts_cutoff %||% 0,
            counts_cutoff_low  = bcfg$counts_cutoff_low %||% NULL
        )

        # Append results
        if (!is.null(bp_res$files)) written <- c(written, bp_res$files)
        if (!is.null(bp_res$plots)) plots <- c(plots, bp_res$plots)

        # Capture for Shiny
        patterns_tbl <- bp_res$best %||% NULL
        heatmaps_by_pattern <- bp_res$plots %||% NULL
        patterns_list <- bp_res$bp_pat %||% NULL
    
    }
    return(list(
        plots = plots,
        files = unique(written),
        excel_order = excel_order,
        objects = list(
            hm_hier_de = pheatmap_payload,
            patterns = patterns_tbl,
            patterns_list = patterns_list,
            heatmaps_by_pattern = heatmaps_by_pattern,
            clusters = clusters_vec
        )
      ))
    
}
