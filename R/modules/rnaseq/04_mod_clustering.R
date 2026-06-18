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
  cl  <- cfg$clustering
  
  written     <- character(0)
  plots       <- list()
  excel_order <- NULL
  objects     <- list(
    patterns            = NULL,
    heatmaps            = NULL,
    clusters            = NULL,
    hm_hier_de          = NULL,
    patterns_list       = NULL,
    heatmaps_by_pattern = NULL
  )
  
  if (is.null(cl) || isFALSE(cl$enabled)) {
    message("[rnaseq clustering] skipped: Clustering disabled in config.")
    return(list(plots = plots, files = written, excel_order = NULL, objects = objects))
  }
  
  if (is.null(cfg$effects$samples)) {
    stop("[rnaseq clustering] effects$samples is missing in config.")
  }
  
  expr_mat <- as.matrix(pre$expr_work)
  if (nrow(expr_mat) < 10) {
    warning("[rnaseq clustering] skipped: Too few features in expr_work (< 10).")
    return(list(plots = plots, files = written, excel_order = NULL, objects = objects))
  }
  
  summary_df <- tryCatch(build_rnaseq_summary_df(de_res$tables, cfg$de),
                         error = function(e) NULL)
  if (is.null(summary_df)) {
    warning("[rnaseq clustering] skipped: Could not build RNA summary_df.")
    return(list(plots = plots, files = written, excel_order = NULL, objects = objects))
  }
  
  de_features <- get_rna_de_features(summary_df)
  message(sprintf("[rnaseq clustering] Found %d DE features for clustering.", length(de_features)))
  manual         <- cfg$de$manual_include %||% character(0)
  manual_in_data <- intersect(manual, rownames(expr_mat))
  n_added        <- length(setdiff(manual_in_data, de_features))
  if (n_added > 0) {
    message(sprintf("DE selection: adding %d manually-included features (cfg$de$manual_include)",
                    n_added))
  }
  de_features <- union(de_features, manual_in_data)
  
  if (length(de_features) < 2) {
    warning("[rnaseq clustering] skipped: < 2 DE features found.")
    return(list(plots = plots, files = written, excel_order = NULL, objects = objects))
  }
  
  clustering_dir <- file.path(out_dir, "Clustering")
  ensure_dir(clustering_dir)
  
  flags <- clustering_run_flags(pre, cfg)
  message(sprintf(
    "[rnaseq clustering] Flags: hierarchical=%s, partition=%s, binary=%s",
    flags$hierarchical, flags$partition, flags$binary_patterns
  ))
  
  annot_col <- build_heatmap_annotation_col(pre$meta, cfg)
  
  de_cfg <- cfg$de %||% list()
  annot_context <- list(
    summary_df    = summary_df,
    p_cutoff      = de_cfg$p_cutoff %||% 0.05,
    log2fc_cutoff = log2(de_cfg$linear_fc_cutoff %||% 1.5),
    id_col        = "FeatureID"
  )
  
  if (isTRUE(flags$hierarchical)) {
    h <- .run_hierarchical_step(expr_mat, pre$meta, de_features, cfg,
                                annot_col, annot_context,
                                file.path(clustering_dir, "Hierarchical"))
    written                <- c(written, h$files)
    plots$p_cluster_hier   <- h$plot
    if (!is.null(h$table_df)) plots$cl_tbl <- h$table_df
    objects$hm_hier_de     <- h$payload
    if (!is.null(h$clusters)) objects$clusters <- h$clusters
    excel_order            <- h$excel_order
  }
  
  if (isTRUE(flags$partition)) {
    p <- .run_partition_step(expr_mat, pre$meta, de_features, cfg, annot_col,
                             file.path(clustering_dir, "Partition_clustering"))
    written                  <- c(written, p$files)
    plots                    <- c(plots, p$plots)
    plots$pt_tbl             <- p$table_df
    if (!is.null(p$clusters)) {
      objects$clusters <- p$clusters
      if (!is.null(excel_order)) {
        excel_order$partition_clusters <- p$clusters
        excel_order$partition_k        <- p$k
      }
    }
  }
  
  if (isTRUE(flags$binary_patterns)) {
    b <- .run_binary_patterns_step(
      expr_mat_corr   = expr_mat,
      expr_mat_counts = as.matrix(pre$expr_filt),
      meta            = pre$meta,
      de_features     = de_features,
      summary_df      = summary_df,
      cfg             = cfg,
      annot_context   = annot_context,
      out_dir         = file.path(clustering_dir, "Binary_patterns")
    )
    written                       <- c(written, b$files)
    objects$patterns              <- b$patterns
    objects$patterns_list         <- b$patterns_list
    objects$heatmaps_by_pattern   <- b$plots
    if (!is.null(excel_order) && !is.null(b$binary_best)) {
      excel_order$binary_best <- b$binary_best
    }
  }
  
  list(
    plots       = plots,
    files       = unique(written),
    excel_order = excel_order,
    objects     = objects
  )
}

# Helper to extract DE features from RNA summary_df
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