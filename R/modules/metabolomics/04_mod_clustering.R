# R/modules/metabolomics/04_mod_clustering.R
#
# Metabolomics clustering module — thin wrapper over core clustering steps.
# Return contract matches proteomics: list(plots, files, excel_order, objects).

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

    excel_order         <- NULL
    pheatmap_payload    <- NULL
    patterns_tbl        <- NULL
    patterns_list       <- NULL
    heatmaps_by_pattern <- NULL
    clusters_vec        <- NULL
    written             <- character(0)
    plots               <- list()

    if (is.null(cl) || isFALSE(cl$enabled)) {
        message("metabolomics clustering: disabled in config — skipping")
        return(list(plots = list(), files = character(0),
                    excel_order = NULL, objects = list()))
    }

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

    clustering_dir <- file.path(out_dir, "Clustering")
    ensure_dir(clustering_dir)

    flags <- clustering_run_flags(pre, cfg)
    message(sprintf(
        "Clustering flags: hierarchical=%s, partition=%s, binary=%s",
        flags$hierarchical, flags$partition, flags$binary_patterns
    ))

    annot_col <- build_heatmap_annotation_col(pre$meta, cfg)
    expr_mat  <- as.matrix(pre$expr_work)

    de_cfg <- cfg$de %||% list()
    annot_context <- list(
      summary_df    = de_res$summary_df,
      p_cutoff      = de_cfg$p_cutoff %||% 0.05,
      log2fc_cutoff = log2(de_cfg$linear_fc_cutoff %||% 1.5),
      id_col        = "feature_id"
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
