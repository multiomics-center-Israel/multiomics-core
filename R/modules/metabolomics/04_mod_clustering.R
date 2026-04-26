# R/modules/metabolomics/04_mod_clustering.R
#
# Metabolomics clustering module — thin wrapper over core clustering functions.
# Return contract matches proteomics: list(plots, files, excel_order, objects).
#
# Reuses: clustering_run_flags, get_color_config, run_clustering, zscore_rows,
#         wrap_clustering_heatmap, build_clustering_output_table,
#         perform_partition_clustering_effects, run_binary_patterns,
#         write_clustering_legacy_profiles, build_cluster_profiles,
#         plot_cluster_profiles_legacy_style, plot_heatmap_core,
#         save_heatmap_to_file


#' Metabolomics clustering module
#'
#' @param pre    Preprocessed metabolomics object (expects $expr_work, $meta).
#' @param de_res DE results list (expects $summary_df with pass_any_contrast).
#' @param config Full pipeline config.
#' @param out_dir Output run directory.
#' @return list(plots, files, excel_order, objects)
mod_metabolomics_clustering <- function(pre, de_res, config, out_dir) {
    stopifnot(is.character(out_dir), length(out_dir) == 1)

    assert_pre_contract(pre, stage = "metabolomics")
    assert_de_contract(de_res, stage = "metabolomics")

    cfg <- config$modes$metabolomics
    cl  <- cfg$clustering

    # ---- Filter to biological samples only (exclude QC/blanks) ----
    sample_col    <- cfg$effects$samples %||% "sample_id"
    condition_col <- (cfg$de %||% list())$condition_column %||%
                     cfg$effects$color %||% "sample_type"
    bio <- filter_to_biological(pre$expr_work, pre$meta, condition_col,
                                sample_col, label = "metabolomics clustering")
    pre$expr_work <- bio$mat
    pre$meta      <- bio$meta

    # Return scaffolds
    excel_order        <- NULL
    pheatmap_payload   <- NULL
    patterns_tbl       <- NULL
    patterns_list      <- NULL
    heatmaps_by_pattern <- NULL
    clusters_vec       <- NULL
    written            <- character(0)
    plots              <- list()

    # If clustering is missing/disabled -> no-op
    if (is.null(cl) || isFALSE(cl$enabled)) {
        message("metabolomics clustering: disabled in config — skipping")
        return(list(plots = list(), files = character(0),
                    excel_order = NULL, objects = list()))
    }

    # ---- Extract DE features ----
    de_features <- get_de_features_metabolomics(de_res)

    if (length(de_features) < 2) {
        message(sprintf(
            "metabolomics clustering: skipped — %d DE features (need >= 2)",
            length(de_features)
        ))
        return(list(plots = list(), files = character(0),
                    excel_order = NULL, objects = list()))
    }

    message(sprintf("metabolomics clustering: %d DE features", length(de_features)))

    # ---- output root ----
    clustering_dir <- file.path(out_dir, "Clustering")
    ensure_dir(clustering_dir)

    # ---- decide which steps to run ----
    flags <- clustering_run_flags(pre, cfg)
    message(sprintf(
        "Clustering flags: hierarchical=%s, partition=%s, binary=%s",
        flags$hierarchical, flags$partition, flags$binary_patterns
    ))

    # ---- build annotation_col for heatmaps ----
    annot_col <- build_heatmap_annotation_col(pre$meta, cfg)

    # Expression matrix
    expr_mat <- as.matrix(pre$expr_work)

    # ------ 1) Hierarchical clustering ---------
    if (isTRUE(flags$hierarchical)) {
        hcfg <- cl$steps$hierarchical %||% list()

        clust_out_dir <- file.path(clustering_dir, "Hierarchical")
        ensure_dir(clust_out_dir)

        hc_res <- run_clustering(
            expr_mat    = expr_mat,
            col_data    = pre$meta,
            de_features = de_features,
            config      = list(
                method   = "hierarchical",
                distance = hcfg$distance %||% "euclidean",
                linkage  = hcfg$linkage  %||% "complete",
                k        = hcfg$k %||% NULL
            )
        )

        mat_de <- expr_mat[de_features, , drop = FALSE]
        z_de   <- zscore_rows(mat_de)
        colnames(z_de) <- paste0(colnames(z_de), ".zscore")

        ordered_row_ids <- intersect(hc_res$ordering, rownames(z_de))
        z_de_ordered <- z_de[ordered_row_ids, , drop = FALSE]

        excel_order <- list(
            ordered_ids        = hc_res$ordering,
            zscore_mat         = z_de,
            partition_clusters = NULL,
            partition_k        = NULL,
            binary_best        = NULL
        )

        # DE row annotations
        de_cfg <- cfg$de %||% list()
        met_p_cutoff  <- de_cfg$p_cutoff %||% 0.05
        met_lin_fc    <- de_cfg$linear_fc_cutoff %||% 1.5
        met_log2fc    <- log2(met_lin_fc)

        annot_context <- list(
            summary_df    = de_res$summary_df,
            p_cutoff      = met_p_cutoff,
            log2fc_cutoff = met_log2fc,
            id_col        = "feature_id"
        )

        f_hm <- file.path(clust_out_dir, "Hierarchical_DE_heatmap.png")

        p_cluster <- wrap_clustering_heatmap(
            expr_mat               = expr_mat,
            meta                   = pre$meta,
            cfg                    = cfg,
            feature_ids            = de_features,
            ordering               = hc_res$ordering,
            annotation_row_builder = TRUE,
            annotation_row_context = annot_context,
            out_file               = f_hm,
            cluster_cols           = FALSE,
            title = sprintf("Hierarchical Clustering (%d DE features)", length(de_features))
        )
        written <- c(written, f_hm)
        plots$p_cluster_hier <- p_cluster

        pheatmap_payload <- list(
            pheatmap      = p_cluster,
            mat           = z_de_ordered,
            row_order     = rownames(z_de_ordered),
            col_order     = colnames(z_de_ordered),
            annotation_col = annot_col,
            feature_ids   = de_features,
            is_zscored    = TRUE,
            cluster_cols  = FALSE,
            tree_row      = hc_res$details
        )

        if (!is.null(hc_res$clusters)) {
            f_tbl <- file.path(clust_out_dir, "Hierarchical_clusters.tsv")
            cl_tbl <- build_clustering_output_table(hc_res$clusters, f_tbl)
            written <- c(written, f_tbl)
            plots$cl_tbl <- cl_tbl
            clusters_vec <- hc_res$clusters
        }
    }

    # ----- 2) Partition clustering ----
    if (isTRUE(flags$partition)) {
        part_dir_name <- "Partition_clustering"
        part_base_dir <- file.path(clustering_dir, part_dir_name)
        ensure_dir(part_base_dir)

        part_res <- perform_partition_clustering_effects(
            expr_mat    = expr_mat,
            meta        = pre$meta,
            cfg         = cfg,
            de_features = de_features
        )

        clusters_vec <- part_res$clusters %||% NULL

        part_dir <- file.path(part_base_dir,
                              sprintf("Partition_clustering_%d_clusters", part_res$k))
        ensure_dir(part_dir)

        # Clusters table
        f_tbl <- file.path(part_dir, "partition_clusters.tsv")
        clusters_tbl <- build_clustering_output_table(part_res$clusters, f_tbl)
        written <- c(written, f_tbl)
        plots$pt_tbl <- clusters_tbl

        # Heatmap
        feats <- names(part_res$clusters)
        mat <- expr_mat[feats, , drop = FALSE]
        ord <- order(part_res$clusters, names(part_res$clusters))
        mat_ord <- mat[ord, , drop = FALSE]

        clusters_ordered <- part_res$clusters[rownames(mat_ord)]
        annot_row <- data.frame(
            Cluster = factor(paste0("C", clusters_ordered)),
            row.names = rownames(mat_ord)
        )

        f_hm <- file.path(part_dir, "Partition_clustering_heatmap.png")

        p_part <- plot_heatmap_core(
            expr_mat       = mat_ord,
            annotation_col = annot_col,
            annotation_row = annot_row,
            title = sprintf("Partition clustering (k=%d) on DE features (n=%d)",
                            part_res$k, nrow(mat_ord)),
            scale_rows   = TRUE,
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            max_rows     = NULL,
            gaps_row     = compute_cluster_gaps(clusters_ordered)
        )

        save_heatmap_to_file(p_part, f_hm)
        plots$partition_heatmap <- p_part
        written <- c(written, f_hm)

        # Per-cluster heatmaps (one PNG per cluster)
        per_clust_hm_files <- save_per_cluster_heatmaps(
            expr_mat       = expr_mat,
            clusters       = part_res$clusters,
            annotation_col = annot_col,
            out_dir        = part_dir
        )
        written <- c(written, per_clust_hm_files)

        # Cluster profile outputs (per-cluster PNGs + multi-panel grid PDF)
        prof_out <- save_cluster_profile_outputs(
            expr_mat = expr_mat,
            meta     = pre$meta,
            clusters = part_res$clusters,
            cfg      = cfg,
            out_dir  = part_dir
        )
        written <- c(written, prof_out$files)
        plots$cluster_profiles <- prof_out$plots

        # Legacy profile exports
        legacy_files <- write_clustering_legacy_profiles(
            expr_mat = expr_mat,
            meta     = pre$meta,
            clusters = part_res$clusters,
            cfg      = cfg,
            out_dir  = part_dir
        )
        written <- c(written, legacy_files)

        # Attach partition results to excel_order
        if (!is.null(excel_order)) {
            excel_order$partition_clusters <- part_res$clusters
            excel_order$partition_k        <- part_res$k
        }
    }

    # ---- 3) Binary patterns ----
    if (isTRUE(flags$binary_patterns)) {
        bcfg <- cl$steps$binary_patterns %||% list()
        clust_out_dir <- file.path(clustering_dir, "Binary_patterns")
        ensure_dir(clust_out_dir)

        bp_res <- run_binary_patterns(
            expr_mat_corr      = expr_mat,
            expr_mat_counts    = NULL,
            meta               = pre$meta,
            cfg                = cfg,
            de_features        = de_features,
            out_dir            = clust_out_dir,
            summary_df         = de_res$summary_df,
            corr_cutoff        = bcfg$corr_cutoff %||% 0.8,
            counts_cutoff_high = bcfg$counts_cutoff_high %||%
                                     bcfg$counts_cutoff %||% 0,
            counts_cutoff_low  = bcfg$counts_cutoff_low %||% NULL
        )

        if (!is.null(bp_res$files)) written <- c(written, bp_res$files)
        if (!is.null(bp_res$plots)) plots <- c(plots, bp_res$plots)

        patterns_tbl        <- bp_res$best %||% NULL
        heatmaps_by_pattern <- bp_res$plots %||% NULL
        patterns_list       <- bp_res$bp_pat %||% NULL

        # Attach binary pattern results to excel_order
        if (!is.null(excel_order) && !is.null(bp_res$best)) {
            excel_order$binary_best <- bp_res$best
        }
    }

    list(
        plots       = plots,
        files       = unique(written),
        excel_order = excel_order,
        objects     = list(
            hm_hier_de          = pheatmap_payload,
            patterns            = patterns_tbl,
            patterns_list       = patterns_list,
            heatmaps_by_pattern = heatmaps_by_pattern,
            clusters            = clusters_vec
        )
    )
}


# ---- module-local helper ---------------------------------------------------

#' Extract DE feature IDs from metabolomics DE results
#'
#' @param de_res DE results with summary_df containing feature_id and
#'   pass_any_contrast columns.
#' @return Character vector of feature IDs passing any contrast.
#' @keywords internal
get_de_features_metabolomics <- function(de_res) {
    df <- de_res$summary_df
    stopifnot(is.data.frame(df))
    stopifnot("feature_id" %in% colnames(df))
    stopifnot("pass_any_contrast" %in% colnames(df))

    feats <- df$feature_id[!is.na(df$pass_any_contrast) & df$pass_any_contrast == 1]
    unique(as.character(feats))
}
