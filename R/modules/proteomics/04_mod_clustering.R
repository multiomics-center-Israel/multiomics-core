#' Proteomics clustering module (legacy-like outputs under out_dir/Clustering)
#'
#' Runs hierarchical clustering on DE features (always, if enabled),
#' and additionally runs partition clustering + binary patterns
#' only when the data has enough groups (>= min_groups), where groups are
#' defined by cfg$effects$color.
#'
#' When \code{de$covariates} was fitted, the hierarchical step clusters the DE
#' matrix with the covariates removed (resolve_proteomics_display_matrix()),
#' so the heatmap shows what the model kept. Its Excel z-score columns stay on
#' the unadjusted matrix, beside the .norm values they summarise; only the row
#' order follows the heatmap. \code{clustering.de_padj_cutoff} /
#' \code{de_log2fc_cutoff} choose the rows by their own rule
#' (get_clustering_de_cutoffs()); without them the rows are pass_any_contrast.
#'
#' @param pre preprocessed proteomics object (expects $expr_imp_single and $meta)
#' @param de_res  DE results list (expects $summary_df at least)
#' @param config  full YAML config list
#' @param out_dir output run directory (project folder)
#' @param display Optional output of resolve_proteomics_display_matrix() (the
#'   prot_display_matrix target); built here when NULL and needed.
#'
#' @return list(plots, files, excel_order, objects)
mod_proteomics_clustering <- function(pre, de_res, config, out_dir, display = NULL) {
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

  # Clear this run's hierarchical heatmap before any path that can skip
  # rewriting it -- the early return below when fewer than two DE features are
  # available, and a run where the hierarchical step is simply off. Otherwise a
  # rerun into the same output directory leaves the previous run's image, which
  # the report both displays and now reads as evidence that clustering ran.
  # Same guard, and the same reason, as PCA_robust.png in 01_mod_qc_pre.R.
  # The hierarchical step recreates the file whenever it actually runs.
  # (Named without its call syntax on purpose: test-clustering-empty-de.R finds
  # the first step by a fixed-string grep, and a comment would shadow it.)
  f_hier_hm <- file.path(clustering_dir, "Hierarchical", "Hierarchical_DE_heatmap.png")
  if (file.exists(f_hier_hm)) file.remove(f_hier_hm)
  # Its PDF and the caption facts the report reads, for the same reason.
  f_hier_extra <- file.path(clustering_dir, "Hierarchical",
                            c("Hierarchical_DE_heatmap.pdf", "Hierarchical_DE_heatmap_info.tsv"))
  unlink(f_hier_extra)

  flags <- clustering_run_flags(pre, cfg)
  message(sprintf(
    "Clustering flags: hierarchical=%s, partition=%s, binary=%s",
    flags$hierarchical, flags$partition, flags$binary_patterns
  ))
  
  annot_col   <- build_heatmap_annotation_col(pre$meta, cfg)
  # Same rows for every step below, chosen by the clustering cut-offs when set.
  de_cutoffs  <- get_clustering_de_cutoffs(cfg)
  summary_cl  <- apply_clustering_de_cutoffs(de_res$summary_df, de_cutoffs)
  de_features <- get_de_features(list(summary_df = summary_cl), cfg)
  expr_mat    <- as.matrix(pre$expr_imp_single)
  matrix_label <- "QC imputation draw (expr_imp_single)"
  if (!is.null(de_res$covariate_design)) {
    display      <- display %||% resolve_proteomics_display_matrix(pre, de_res, config)
    expr_mat     <- display$mat
    matrix_label <- display$label
  }

  # Same guard the metabolomics and lipidomics modules already carry: a DE step
  # that returns nothing is a legitimate result, not a pipeline failure.
  # run_clustering() rightly refuses an empty feature set, but without this the
  # error aborted the whole run -- so an unmoderated or unblocked comparison,
  # expected to return few or no features at this n, could not reach its report.
  n_available <- length(intersect(de_features, rownames(expr_mat)))
  if (n_available < 2) {
    message(sprintf(
      "proteomics clustering: skipped — %d DE features present in the matrix (need >= 2)",
      n_available
    ))
    return(list(plots = list(), files = character(0),
                excel_order = NULL, objects = list()))
  }
  
  prot_de_cfg <- cfg$de %||% list()
  annot_context <- list(
    summary_df    = summary_cl,
    p_cutoff      = prot_de_cfg$p_cutoff %||% 0.05,
    log2fc_cutoff = log2(prot_de_cfg$linear_fc_cutoff %||% 1.5),
    id_col        = cfg$de_table$id_col %||% "FeatureID"
  )
  
  if (isTRUE(flags$hierarchical)) {
    annot_hier <- annotation_for_pheatmap(build_proteomics_sample_annotation(pre$meta, cfg), cfg)
    zlim <- (cl$heatmap %||% list())[["zscore_limit"]]
    h <- .run_hierarchical_step(expr_mat, pre$meta, de_features, cfg,
                                annot_hier, annot_context,
                                file.path(clustering_dir, "Hierarchical"),
                                zscore_limit = zlim,
                                annotation_colors = build_sample_annotation_colors(annot_hier))
    if (!is.null(de_res$covariate_design) && !is.null(h$excel_order)) {
      # The z-score columns sit next to the unadjusted .norm values in the
      # workbook, so they describe the unadjusted matrix, as before.
      ids <- rownames(h$excel_order$zscore_mat)
      z_unadj <- zscore_rows(as.matrix(pre$expr_imp_single)[ids, , drop = FALSE])
      colnames(z_unadj) <- paste0(colnames(z_unadj), ".zscore")
      h$excel_order$zscore_mat <- z_unadj
    }
    f_info <- f_hier_extra[2]
    if (file.exists(f_hier_hm)) {
      info <- data.frame(
        n_features    = length(intersect(de_features, rownames(expr_mat))),
        row_rule      = if (is.null(de_cutoffs)) {
          "pass_any_contrast (DE step's own rule)"
        } else {
          sprintf("padj.imputs < %s and abs(log2FC.imputs) >= %s in at least one contrast",
                  format(de_cutoffs$padj_cutoff), format(de_cutoffs$log2fc_cutoff))
        },
        matrix        = matrix_label,
        distance      = (cfg$clustering$steps$hierarchical %||% list())$distance %||% "euclidean",
        linkage       = (cfg$clustering$steps$hierarchical %||% list())$linkage %||% "complete",
        samples_clustered = isTRUE(cfg$clustering$cluster_samples),
        zscore_limit  = if (is.null(zlim)) NA_real_ else as.numeric(zlim),
        stringsAsFactors = FALSE
      )
      utils::write.table(info, f_info, sep = "\t", quote = FALSE, row.names = FALSE)
      h$files <- c(h$files, f_info)
    }
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
      summary_df      = summary_cl,
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
