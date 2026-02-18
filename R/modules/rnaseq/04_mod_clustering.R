#' RNA-seq clustering module (legacy-like outputs under out_dir/Clustering)
#'
#' Runs hierarchical clustering on DE features (always, if enabled),
#' and additionally runs partition clustering + binary patterns
#' only when the data has enough groups.
#'
#' @param pre preprocessed RNA object (expects $expr_work and $meta)
#' @param de_res DE results list (will build summary_df from this)
#' @param config full YAML config list
#' @param out_dir output run directory
#'
#' @return list(plots, files, objects)
mod_rnaseq_clustering <- function(pre, de_res, config, out_dir) {
    stopifnot(is.character(out_dir), length(out_dir) == 1)

    assert_pre_contract(pre, stage = "rna")

    cfg <- config$modes$rna
    cl <- cfg$clustering

    eff_color <- get_color_config(cfg)

    # Initialize Objects for legacy Shiny export (MUST exist)
    objects <- list(
        patterns = NULL,
        heatmaps = NULL,
        clusters = NULL # New_clusters
    )

    excel_order <- NULL
    plots <- list()
    written <- character(0)

    # 1. Verify module execution conditions
    if (is.null(cl) || isFALSE(cl$enabled)) {
        message("[rnaseq clustering] skipped: Clustering disabled in config.")
        return(list(plots = plots, files = written, objects = objects))
    }

    # 2. Validate prerequisites explicitly
    # Check 1: Effects definition
    if (is.null(eff_color) || is.null(cfg$effects$samples)) {
        warning("[rnaseq clustering] skipped: effects$color or effects$samples missing in config.")
        return(list(plots = plots, files = written, objects = objects))
    }

    # Check 2: Expression matrix dims
    expr_mat <- as.matrix(pre$expr_work)
    if (nrow(expr_mat) < 10) { # Arbitrary small number validation
        warning("[rnaseq clustering] skipped: Too few features in expr_work (< 10).")
        return(list(plots = plots, files = written, objects = objects))
    }

    # Check 3: Metadata alignment
    if (!all(colnames(expr_mat) %in% rownames(pre$meta))) {
        stop("CRITICAL: expr_work columns do not match meta rownames.")
    }

    # ---- output root (UNDER PROJECT RUN DIR) ----
    clustering_dir <- file.path(out_dir, "Clustering")
    ensure_dir(clustering_dir)

    # ---- decide which steps to run ----
    flags <- clustering_run_flags(pre, cfg)
    message(sprintf(
        "[rnaseq clustering] Flags: hierarchical=%s, partition=%s, binary=%s",
        flags$hierarchical, flags$partition, flags$binary_patterns
    ))

    # ---- build annotation_col for heatmaps using effects ----
    annot <- NULL
    if (eff_color %in% colnames(pre$meta)) {
        annot <- data.frame(
            Condition = pre$meta[[eff_color]],
            row.names = pre$meta[[cfg$effects$samples]]
        )
    }

    # Rebuild summary_df to identify DE features
    summary_df <- tryCatch(build_rnaseq_summary_df(de_res$tables, config$modes$rna$de), error = function(e) NULL)
    if (is.null(summary_df)) {
        warning("[rnaseq clustering] skipped: Could not build RNA summary_df.")
        return(list(plots = plots, files = written, objects = objects))
    }

    # Get DE ids
    de_features <- get_rna_de_features(summary_df)
    message(sprintf("[rnaseq clustering] Found %d DE features for clustering.", length(de_features)))

    if (length(de_features) < 2) {
        warning("[rnaseq clustering] skipped: < 2 DE features found.")
        return(list(plots = plots, files = written, objects = objects))
    }

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

        mat_de <- expr_mat[intersect(de_features, rownames(expr_mat)), , drop = FALSE]
        z_de <- zscore_rows(mat_de)
        colnames(z_de) <- paste0(colnames(z_de), ".zscore")

        excel_order <- list(
            ordered_ids = hc_res$ordering,
            zscore_mat  = z_de
        )

        # Build DE pattern row annotations (up/down per contrast)
        de_cfg <- config$modes$rna$de %||% list()
        p_cutoff <- de_cfg$p_cutoff %||% 0.05
        lin_fc_cutoff <- de_cfg$linear_fc_cutoff %||% 1.5
        log2fc_cutoff <- log2(lin_fc_cutoff)

        row_annot <- build_de_row_annotations(
            summary_df    = summary_df,
            feature_ids   = de_features,
            p_cutoff      = p_cutoff,
            log2fc_cutoff = log2fc_cutoff
        )

        # Heatmap
        f_hm <- file.path(clust_out_dir, "Hierarchical_DE_heatmap.png")
        p_cluster <- wrap_clustering_heatmap(
            expr_mat       = expr_mat,
            meta           = pre$meta,
            cfg            = cfg,
            feature_ids    = de_features,
            ordering       = hc_res$ordering,
            annotation_row = row_annot,
            out_file       = f_hm
        )
        written <- c(written, f_hm)
        plots$p_cluster <- p_cluster

        # Save clusters
        if (!is.null(hc_res$clusters)) {
            f_tbl <- file.path(clust_out_dir, "Hierarchical_clusters.tsv")
            build_clustering_output_table(hc_res$clusters, f_tbl)
            written <- c(written, f_tbl)

            # Populate Shiny Object
            objects$clusters <- hc_res$clusters
        }
    }

    # ----- 2) Partition clustering (kmeans / PAM) ----
    if (isTRUE(flags$partition)) {
        part_dir_name <- "Partition_clustering"
        part_base_dir <- file.path(clustering_dir, part_dir_name)
        ensure_dir(part_base_dir)

        # Fit clusters
        part_res <- perform_partition_clustering_effects(
            expr_mat = expr_mat,
            meta = pre$meta,
            cfg = cfg,
            de_features = de_features
        )

        # Overwrite clusters with partition result if available (Proteomics behavior preference)
        if (!is.null(part_res$clusters)) {
            objects$clusters <- part_res$clusters
        }

        part_dir <- file.path(part_base_dir, sprintf("Partition_clustering_%d_clusters", part_res$k))
        ensure_dir(part_dir)

        # (1) write clusters table
        f_tbl <- file.path(part_dir, "partition_clusters.tsv")
        clusters_tbl <- build_clustering_output_table(part_res$clusters, f_tbl)
        written <- c(written, f_tbl)

        # (2) heatmap
        feats <- names(part_res$clusters)
        # Use intersection to be safe
        valid_feats <- intersect(feats, rownames(expr_mat))
        mat <- expr_mat[valid_feats, , drop = FALSE]

        ord <- order(part_res$clusters[valid_feats], valid_feats)
        mat_ord <- mat[ord, , drop = FALSE]

        # Task 6: Multi-column annotations
        annot_cols_config <- cfg$effects$heatmap_annotations %||% NULL
        if (!is.null(annot_cols_config)) {
            annot_col <- pre$meta[, annot_cols_config, drop = FALSE]
            rownames(annot_col) <- pre$meta[[cfg$effects$samples]]
            annot_col <- annot_col[colnames(mat_ord), , drop = FALSE]
        } else {
            # Fallback to single column
            annot_col <- annot[colnames(mat_ord), , drop = FALSE]
        }

        # Task 6: Row annotations showing cluster assignments
        clusters_ordered <- part_res$clusters[rownames(mat_ord)]
        annot_row <- data.frame(
            Cluster = factor(paste0("C", clusters_ordered)),
            row.names = rownames(mat_ord)
        )

        # Task 6: Compute gaps_row for visual cluster separation
        cluster_sizes <- table(clusters_ordered)[unique(clusters_ordered)]
        gaps_row <- cumsum(cluster_sizes)[-length(cluster_sizes)]

        f_hm <- file.path(part_dir, "Partition_clustering_heatmap.png")

        p_part <- plot_heatmap_core(
            expr_mat         = mat_ord,
            annotation_col   = annot_col,
            annotation_row   = annot_row,
            gaps_row         = gaps_row,
            title            = sprintf("Partition clustering (k=%d)", part_res$k),
            scale_rows       = TRUE,
            cluster_rows     = FALSE,
            cluster_cols     = TRUE,
            max_rows         = NULL
        )
        save_heatmap_to_file(p_part, f_hm)
        written <- c(written, f_hm)

        # (3) cluster profiles pdf
        prof <- build_cluster_profiles(part_res$z_group_means, part_res$clusters, part_res$k)

        if (!is.null(prof)) {
            f_pdf <- file.path(part_dir, "cluster_profiles.pdf")
            eff_col_name <- eff_color %||% "Group"

            p_prof <- plot_cluster_profiles(prof, x_label = eff_col_name)

            n_clusters <- length(unique(prof$cluster))
            calc_height <- max(6, ceiling(n_clusters / 2) * 3)

            ggplot2::ggsave(f_pdf, plot = p_prof, width = 10, height = calc_height)
            written <- c(written, f_pdf)
        }
    }

    # ---- 3) Binary patterns ----
    if (isTRUE(flags$binary_patterns)) {
        bcfg <- cl$steps$binary_patterns %||% list()
        clust_out_dir <- file.path(clustering_dir, "Binary_patterns")
        ensure_dir(clust_out_dir)

        bp_res <- run_binary_patterns(
            expr_mat      = expr_mat,
            meta          = pre$meta,
            cfg           = cfg,
            de_features   = de_features,
            out_dir       = clust_out_dir,
            corr_cutoff   = bcfg$corr_cutoff %||% 0.8,
            counts_cutoff = bcfg$counts_cutoff %||% 0,
            summary_df    = summary_df,
            p_cutoff      = p_cutoff,
            log2fc_cutoff = log2fc_cutoff
        )

        if (!is.null(bp_res$files)) written <- c(written, bp_res$files)

        # Populate Shiny Objects
        objects$patterns <- bp_res$best %||% NULL
        objects$heatmaps <- bp_res$plots %||% NULL
    }

    return(list(
        plots = plots,
        files = unique(written),
        excel_order = excel_order,
        objects = objects
    ))
}

# Helper to exact DE features from RNA summary_df
get_rna_de_features <- function(summary_df) {
    if (is.null(summary_df) || nrow(summary_df) == 0) {
        return(character(0))
    }

    # Try generic pass column
    if ("pass_any_contrast" %in% colnames(summary_df)) {
        return(summary_df$FeatureID[which(summary_df$pass_any_contrast == 1)])
    }

    # Fallback: specific pass columns
    pass_cols <- grep("^sum\\.pass\\.", colnames(summary_df), value = TRUE)
    if (length(pass_cols) > 0) {
        row_sums <- rowSums(summary_df[, pass_cols, drop = FALSE], na.rm = TRUE)
        return(summary_df$FeatureID[which(row_sums > 0)])
    }

    character(0)
}
