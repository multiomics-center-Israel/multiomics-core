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
#' @return list(plots, files, excel_order, objects)
mod_proteomics_clustering <- function(pre, de_res, config, out_dir) {
  stopifnot(is.character(out_dir), length(out_dir) == 1)

  assert_pre_contract(pre, stage = "proteomics")
  assert_de_contract(de_res, stage = "proteomics")

  cfg <- config$modes$proteomics
  cl  <- cfg$clustering

  pheatmap_payload    <- NULL
  patterns_tbl        <- NULL
  patterns_list       <- NULL
  heatmaps_by_pattern <- NULL
  clusters_vec        <- NULL
  excel_order         <- NULL
  written             <- character(0)
  plots               <- list()

  if (is.null(cl) || isFALSE(cl$enabled)) {
    message("Clustering disabled. Skipping.")
    return(list(plots = list(), files = character(0),
                excel_order = NULL, objects = list()))
  }

  clustering_dir <- file.path(out_dir, "Clustering")
  ensure_dir(clustering_dir)

  flags <- clustering_run_flags(pre, cfg)
  message(sprintf(
    "Clustering flags: hierarchical=%s, partition=%s, binary=%s",
    flags$hierarchical, flags$partition, flags$binary_patterns
  ))

  annot_col   <- build_heatmap_annotation_col(pre$meta, cfg)
  de_features <- get_de_features(de_res, cfg)
  expr_mat    <- as.matrix(pre$expr_imp_single)

  prot_de_cfg <- cfg$de %||% list()
  annot_context <- list(
    summary_df    = de_res$summary_df,
    p_cutoff      = prot_de_cfg$p_cutoff %||% 0.05,
    log2fc_cutoff = log2(prot_de_cfg$linear_fc_cutoff %||% 1.5),
    id_col        = cfg$de_table$id_col %||% "FeatureID"
  )

  if (isTRUE(flags$hierarchical)) {
    h <- .run_hierarchical_step(expr_mat, pre$meta, de_features, cfg,
                                annot_col, annot_context,
                                file.path(clustering_dir, "Hierarchical"))
    written              <- c(written, h$files)
    plots$p_cluster_hier <- h$plot
    if (!is.null(h$table_df)) plots$cl_tbl <- h$table_df
    pheatmap_payload     <- h$payload
    if (!is.null(h$clusters)) clusters_vec <- h$clusters
    excel_order          <- h$excel_order
  }

  if (isTRUE(flags$partition)) {
    p <- .run_partition_step(expr_mat, pre$meta, de_features, cfg, annot_col,
                             file.path(clustering_dir, "Partition_clustering"))
    written      <- c(written, p$files)
    plots        <- c(plots, p$plots)
    plots$pt_tbl <- p$table_df
    clusters_vec <- p$clusters %||% clusters_vec
    if (!is.null(excel_order)) {
      excel_order$partition_clusters <- p$clusters
      excel_order$partition_k        <- p$k
    }
  }

  if (isTRUE(flags$binary_patterns)) {
    b <- .run_binary_patterns_step(
      expr_mat_corr   = expr_mat,
      expr_mat_counts = NULL,
      meta            = pre$meta,
      de_features     = de_features,
      summary_df      = de_res$summary_df,
      cfg             = cfg,
      annot_context   = annot_context,
      out_dir         = file.path(clustering_dir, "Binary_patterns")
    )
    written             <- c(written, b$files)
    if (length(b$plots)) plots <- c(plots, b$plots)
    patterns_tbl        <- b$patterns
    heatmaps_by_pattern <- b$plots
    patterns_list       <- b$patterns_list
    if (!is.null(excel_order) && !is.null(b$binary_best)) {
      excel_order$binary_best <- b$binary_best
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
