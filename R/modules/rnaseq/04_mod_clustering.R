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
#' @return list(plots, files, excel_order, objects)
mod_rnaseq_clustering <- function(pre, de_res, config, out_dir) {
  stopifnot(is.character(out_dir), length(out_dir) == 1)
  
  assert_pre_contract(pre, stage = "rna")
  
  cfg <- config$modes$rna
  cl <- cfg$clustering
  
  # Initialize Scaffolds and Objects (MUST exist for consistent return)
  written <- character(0)
  plots <- list()
  excel_order <- NULL 
  
  objects <- list(
    patterns = NULL,
    heatmaps = NULL,
    clusters = NULL, 
    hm_hier_de = NULL,
    patterns_list = NULL,
    heatmaps_by_pattern = NULL
  )
  
  # 1. Verify module execution conditions
  if (is.null(cl) || isFALSE(cl$enabled)) {
    message("[rnaseq clustering] skipped: Clustering disabled in config.")
    return(list(plots = plots, files = written, excel_order = NULL, objects = objects))
  }
  
  # 2. Validate prerequisites explicitly
  if (is.null(cfg$effects$samples)) {
    stop("[rnaseq clustering] effects$samples is missing in config.")
  }
  
  expr_mat <- as.matrix(pre$expr_work)
  if (nrow(expr_mat) < 10) {
    warning("[rnaseq clustering] skipped: Too few features in expr_work (< 10).")
    return(list(plots = plots, files = written, excel_order = NULL, objects = objects))
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
  annot_col <- build_heatmap_annotation_col(pre$meta, cfg)
  
  # Rebuild summary_df to identify DE features
  summary_df <- tryCatch(build_rnaseq_summary_df(de_res$tables, cfg$de), 
                         error = function(e) NULL)
  
  if (is.null(summary_df)) {
    warning("[rnaseq clustering] skipped: Could not build RNA summary_df.")
    return(list(plots = plots, files = written, excel_order = NULL, objects = objects))
  }
  
  # Get DE ids
  de_features <- get_rna_de_features(summary_df)
  message(sprintf("[rnaseq clustering] Found %d DE features for clustering.", length(de_features)))
  
  if (length(de_features) < 2) {
    warning("[rnaseq clustering] skipped: < 2 DE features found.")
    return(list(plots = plots, files = written, excel_order = NULL, objects = objects))
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
        linkage  = hcfg$linkage  %||% "complete",
        k        = hcfg$k %||% NULL
      )
    )
    
    mat_de <- expr_mat[intersect(de_features, rownames(expr_mat)), , drop = FALSE]
    z_de <- zscore_rows(mat_de)
    colnames(z_de) <- paste0(colnames(z_de), ".zscore")
    
    # Order rows according to clustering result
    ordered_row_ids <- intersect(hc_res$ordering, rownames(z_de))
    z_de_ordered <- z_de[ordered_row_ids, , drop = FALSE]
    
    # Initialize excel_order with Hierarchical results
    excel_order <- list(
      ordered_ids        = hc_res$ordering,
      zscore_mat         = z_de,
      partition_clusters = NULL,
      partition_k        = NULL,
      binary_best        = NULL
    )
    
    # Build DE pattern row annotations
    de_cfg <- cfg$de %||% list()
    annot_context <- list(
      summary_df    = summary_df,
      p_cutoff      = de_cfg$p_cutoff %||% 0.05,
      log2fc_cutoff = log2(de_cfg$linear_fc_cutoff %||% 1.5),
      id_col        = "FeatureID"
    )
    
    # Heatmap
    f_hm <- file.path(clust_out_dir, "Hierarchical_DE_heatmap.png")
    p_cluster <- wrap_clustering_heatmap(
      expr_mat = expr_mat,
      meta = pre$meta,
      cfg = cfg,
      feature_ids = de_features,
      ordering = hc_res$ordering,
      annotation_row_builder = TRUE,
      annotation_row_context = annot_context,
      out_file = f_hm,
      cluster_cols = FALSE,
      title = sprintf("Hierarchical Clustering (%d DE features)", length(de_features))
    )
    
    written <- c(written, f_hm)
    plots$p_cluster_hier <- p_cluster
    
    # Capture pheatmap payload for Shiny
    objects$hm_hier_de <- list(
      pheatmap = p_cluster,
      mat = z_de_ordered, 
      row_order = rownames(z_de_ordered),
      col_order = colnames(z_de_ordered),
      annotation_col = annot_col,
      feature_ids = de_features,
      is_zscored = TRUE,
      tree_row = hc_res$details 
    )
    
    # Save clusters
    if (!is.null(hc_res$clusters)) {
      f_tbl <- file.path(clust_out_dir, "Hierarchical_clusters.tsv")
      build_clustering_output_table(hc_res$clusters, f_tbl)
      written <- c(written, f_tbl)
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
    
    # Update Shiny objects and excel_order
    if (!is.null(part_res$clusters)) {
      objects$clusters <- part_res$clusters
      #  Attach partition results to excel_order
      if (!is.null(excel_order)) {
        excel_order$partition_clusters <- part_res$clusters
        excel_order$partition_k        <- part_res$k
      }
    } else {
      message("[rnaseq clustering] Partition clustering skipped (too few features).")
    }

    if (is.null(part_res$clusters)) {
      # Nothing to plot — skip the rest of the partition block
    } else {

    part_dir <- file.path(part_base_dir, sprintf("Partition_clustering_%d_clusters", part_res$k))
    ensure_dir(part_dir)
    
    # (1) write clusters table
    f_tbl <- file.path(part_dir, "partition_clusters.tsv")
    build_clustering_output_table(part_res$clusters, f_tbl)
    written <- c(written, f_tbl)
    
    # (2) Heatmap
    feats <- names(part_res$clusters)
    valid_feats <- intersect(feats, rownames(expr_mat))
    mat_ord <- expr_mat[order(part_res$clusters[valid_feats], valid_feats), , drop = FALSE]
    
    annot_row <- data.frame(
      Cluster = factor(paste0("C", part_res$clusters[rownames(mat_ord)])),
      row.names = rownames(mat_ord)
    )
    
    f_hm <- file.path(part_dir, "Partition_clustering_heatmap.png")
    p_part <- plot_heatmap_core(
      expr_mat       = mat_ord,
      annotation_col = annot_col,
      annotation_row = annot_row,
      title          = sprintf("Partition clustering (k=%d)", part_res$k),
      scale_rows     = TRUE,
      cluster_rows   = FALSE,
      cluster_cols   = FALSE,
      max_rows       = NULL
    )
    
    plots$partition_heatmap <- p_part
    save_heatmap_to_file(p_part, f_hm)
    written <- c(written, f_hm)
    
    # (3) cluster profiles
    prof <- build_cluster_profiles(part_res$group_means, part_res$clusters, part_res$k)
    if (!is.null(prof)) {
      f_pdf <- file.path(part_dir, "cluster_profiles.pdf")
      p_prof <- plot_cluster_profiles_legacy_style(
        group_means = part_res$group_means,
        clusters = part_res$clusters,
        x_label = cfg$clustering$group_col %||% "Group"
      )
      ggplot2::ggsave(f_pdf, plot = p_prof, width = 10, height = max(6, ceiling(part_res$k / 2) * 3))
      written <- c(written, f_pdf)
    }
    } # end else (partition has clusters)
  }

  # ---- 3) Binary patterns ----
  if (isTRUE(flags$binary_patterns)) {
    bcfg <- cl$steps$binary_patterns %||% list()
    clust_out_dir <- file.path(clustering_dir, "Binary_patterns")
    ensure_dir(clust_out_dir)
    
    bp_res <- run_binary_patterns(
      expr_mat_corr      = expr_mat, 
      expr_mat_counts    = as.matrix(pre$expr_filt), 
      meta               = pre$meta,
      cfg                = cfg,
      de_features        = de_features,
      summary_df         = summary_df,
      out_dir            = clust_out_dir,
      corr_cutoff        = bcfg$corr_cutoff %||% 0.8,
      counts_cutoff_high = bcfg$counts_cutoff_high %||% bcfg$counts_cutoff %||% 0,
      counts_cutoff_low  = bcfg$counts_cutoff_low %||% NULL
    )
    
    if (!is.null(bp_res$files)) written <- c(written, bp_res$files)
    
    # Populate Shiny Objects
    objects$patterns <- bp_res$best %||% NULL
    objects$patterns_list <- bp_res$bp_pat %||% NULL
    objects$heatmaps_by_pattern <- bp_res$plots %||% NULL
    
    #  Attach binary pattern results to excel_order
    if (!is.null(excel_order) && !is.null(bp_res$best)) {
      excel_order$binary_best <- bp_res$best
    }
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
  if (is.null(summary_df) || nrow(summary_df) == 0) return(character(0))
  
  if ("pass_any_contrast" %in% colnames(summary_df)) {
    return(summary_df$FeatureID[which(summary_df$pass_any_contrast == 1)])
  }
  
  pass_cols <- grep("^sum\\.pass\\.", colnames(summary_df), value = TRUE)
  if (length(pass_cols) > 0) {
    row_sums <- rowSums(summary_df[, pass_cols, drop = FALSE], na.rm = TRUE)
    return(summary_df$FeatureID[which(row_sums > 0)])
  }
  
  character(0)
}