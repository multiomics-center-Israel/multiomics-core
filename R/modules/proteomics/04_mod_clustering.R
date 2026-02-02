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
    heatmaps_by_pattern <- NULL
    clusters_vec <- NULL

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
    annot <- NULL
    if (!is.null(cfg$effects$color) && cfg$effects$color %in% colnames(pre$meta)) {
        annot <- data.frame(
            Condition = pre$meta[[cfg$effects$color]],
            row.names = pre$meta[[cfg$effects$samples]]
        )
    }

    # Get DE features
    de_features <- get_de_features(de_res, cfg)

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

        excel_order <- list(
            ordered_ids = hc_res$ordering,
            zscore_mat  = z_de
        )

        # Capture payload for Shiny
        pheatmap_payload <- list(
            mat = z_de,
            annotation_col = annot,
            feature_ids = de_features,
            is_zscored = TRUE
        )

        # Heatmap setup
        f_hm <- file.path(clust_out_dir, "Hierarchical_DE_heatmap.png")

        # Run wrapper
        p_cluster <- wrap_clustering_heatmap(
            expr_mat    = pre$expr_imp_single,
            meta        = pre$meta,
            cfg         = cfg,
            feature_ids = de_features,
            ordering    = hc_res$ordering,
            out_file    = f_hm
        )
        written <- c(written, f_hm)
        plots$p_cluster <- p_cluster

        # Save cluster assignments
        if (!is.null(hc_res$clusters)) {
            f_tbl <- file.path(clust_out_dir, "Hierarchical_clusters.tsv")
            cl_tbl <- data.frame(
                feature_id = names(hc_res$clusters),
                cluster = as.integer(hc_res$clusters),
                stringsAsFactors = FALSE
            )
            save_tsv_path(cl_tbl, f_tbl)
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
        clusters_tbl <- data.frame(
            feature_id = names(part_res$clusters),
            cluster = as.integer(part_res$clusters),
            stringsAsFactors = FALSE
        )
        f_tbl <- file.path(part_dir, "partition_clusters.tsv")
        save_tsv_path(clusters_tbl, f_tbl)
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
        # Prepare data from Z-scored Group Means (returned by the clustering func)
        zgm <- part_res$z_group_means
        clv <- part_res$clusters[rownames(zgm)]

        k <- part_res$k
        groups <- colnames(zgm)

        prof_list <- lapply(1:k, function(ci) {
            rows <- which(clv == ci)
            sub <- zgm[rows, , drop = FALSE]
            data.frame(
                cluster = ci,
                group = groups,
                mean = colMeans(sub, na.rm = TRUE),
                sd = apply(sub, 2, stats::sd, na.rm = TRUE),
                n_features = nrow(sub),
                stringsAsFactors = FALSE
            )
        })

        prof <- do.call(rbind, prof_list)

        f_pdf <- file.path(part_dir, "cluster_profiles.pdf")
        p_prof <- plot_cluster_profiles(prof, x_label = cfg$effects$color)

        # Dynamic height
        n_clusters <- length(unique(prof$cluster))
        calc_height <- max(6, ceiling(n_clusters / 2) * 3)

        ggplot2::ggsave(f_pdf, plot = p_prof, width = 10, height = calc_height)

        written <- c(written, f_pdf)
        plots$cluster_profiles <- p_prof

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

        # Run function and capture result
        bp_res <- run_binary_patterns(
            expr_mat      = expr_mat,
            meta          = pre$meta,
            cfg           = cfg,
            de_features   = de_features,
            out_dir       = clust_out_dir,
            corr_cutoff   = bcfg$corr_cutoff %||% 0.8,
            counts_cutoff = bcfg$counts_cutoff %||% 0
        )

        # Append results
        if (!is.null(bp_res$files)) written <- c(written, bp_res$files)
        if (!is.null(bp_res$plots)) plots <- c(plots, bp_res$plots)

        # Capture for Shiny
        patterns_tbl <- bp_res$best %||% NULL
        heatmaps_by_pattern <- bp_res$plots %||% NULL
    }

    return(list(
        plots = plots,
        files = unique(written),
        excel_order = excel_order,
        objects = list(
            pheatmap_data_DE_genes = pheatmap_payload,
            patterns = patterns_tbl,
            heatmaps_by_pattern = heatmaps_by_pattern,
            clusters = clusters_vec
        )
    ))
}
