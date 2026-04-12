# ============================================================
# Shiny Export — Metabolomics
# ============================================================
# This file contains the canonical builder for metabolomics.
# Since metabolomics is newly added to the pipeline, there is
# no legacy builder - only the canonical contract (v2.0).
#
# New code should use build_shiny_payload_metabolomics().
# ============================================================


# ============================================================
# CANONICAL BUILDER (v2.0)
# ============================================================

#' Build canonical Shiny payload for Metabolomics
#'
#' Creates a Shiny payload conforming to the canonical contract (v2.0).
#' All 26 keys are guaranteed to exist (NULL if not applicable).
#'
#' @param pre Preprocessing results (from preprocess_metabolomics)
#' @param de_res DE results (from metabolomics DE analysis). Can be NULL if DE was skipped.
#' @param inputs Input list (contrasts, metadata, etc.)
#' @param config Full config object
#' @param pca_res Optional: pre-computed PCA results
#' @param clustering_res Optional: pre-computed clustering results
#' @param rf_res Optional: random forest results (from feature_sel_res$rf)
#' @param plsda_res Optional: PLS-DA results (from feature_sel_res$plsda)
#' @param enrichment_res Optional: enrichment results (from mod_metabolomics_enrichment)
#' @param annot Optional: external annotation data.frame
#'
#' @return A named list with 26 canonical keys (+ legacy aliases if requested)
#'
#' @export
build_shiny_payload_metabolomics <- function(
    pre,
    de_res = NULL,
    inputs,
    config,
    pca_res = NULL,
    clustering_res = NULL,
    rf_res = NULL,
    plsda_res = NULL,
    enrichment_res = NULL,
    annot = NULL,
    include_legacy = TRUE,
    out_dir = NULL
) {
    # ============================================================
    # Initialize canonical payload structure
    # ============================================================
    payload <- init_shiny_payload("metabolomics")

    # ============================================================
    # Extract config shortcuts
    # ============================================================
    modes <- config$modes %||% list()
    metab_cfg <- modes$metabolomics %||% list()
    de_cfg <- metab_cfg$de %||% list()
    norm_cfg <- metab_cfg$normalization %||% list()
    effects_cfg <- metab_cfg$effects %||% list()

    # ============================================================
    # METADATA (3 keys)
    # ============================================================

    # sample_meta: Sample metadata with rownames as sample IDs
    payload$sample_meta <- pre$meta
    if (!is.null(payload$sample_meta)) {
        sample_col <- effects_cfg$samples %||% "sample_id"
        if (sample_col %in% colnames(payload$sample_meta)) {
            rownames(payload$sample_meta) <- as.character(payload$sample_meta[[sample_col]])
        }
    }

    # contrasts: Contrast definitions
    payload$contrasts <- inputs$contrasts %||% NULL

    # feature_annot: Feature annotations (metabolite names, m/z, RT, etc.)
    if (!is.null(annot)) {
        payload$feature_annot <- annot
    } else if (!is.null(pre$row_data)) {
        # Use row_data as feature annotations if no external annot provided
        payload$feature_annot <- pre$row_data
    }

    # ============================================================
    # EXPRESSION DATA (2 keys)
    # ============================================================

    # expr_raw: Filtered expression (before normalization, may have NAs)
    payload$expr_raw <- pre$expr_filt

    # expr_norm: Normalized expression (NO NAs allowed)
    # In metabolomics, expr_work is the normalized matrix
    # Check for NA policy - if na_policy="zero", NAs were replaced
    payload$expr_norm <- pre$expr_work

    # expr_long: Long-format expression with metadata
    payload$expr_long <- build_expr_long(payload$expr_norm, payload$sample_meta)

    # Handle NAs in expr_norm if present (warn but don't fail)
    if (!is.null(payload$expr_norm) && anyNA(payload$expr_norm)) {
        na_count <- sum(is.na(payload$expr_norm))
        warning(
            "metabolomics expr_norm contains ", na_count, " NA values. ",
            "Consider using na_policy='zero' in config or imputation. ",
            "Replacing NAs with row medians for contract compliance."
        )
        # Replace NAs with row medians as fallback
        for (i in seq_len(nrow(payload$expr_norm))) {
            row_i <- payload$expr_norm[i, ]
            na_mask <- is.na(row_i)
            if (any(na_mask)) {
                row_median <- median(row_i, na.rm = TRUE)
                if (is.na(row_median)) row_median <- 0
                payload$expr_norm[i, na_mask] <- row_median
            }
        }
    }

    # ============================================================
    # QC/PCA (3 keys)
    # ============================================================

    if (!is.null(pca_res)) {
        # Handle nested structure (objects list) or flat structure
        pca_objects <- pca_res$objects %||% pca_res

        # pca_object: prcomp result
        payload$pca_object <- pca_objects$norm_log_counts_pca %||%
                              pca_objects$pca_object %||%
                              pca_objects$pca %||%
                              NULL

        # pca_scores: PCA scores data.frame with metadata
        payload$pca_scores <- pca_objects$pca_scores %||% NULL

        # pca_3d: 3D PCA plotly widget
        payload$pca_3d <- pca_res$plots$pca_3d %||% NULL
        
        # QC plot
        payload$imp_hist_samp <- pca_res$plots$imputation_hist %||% NULL
        payload$samples_hm_w_qc <- pca_res$plots$dist_heatmap %||% NULL
        payload$samples_hm <- pca_res$plots$dist_heatmap_noQC %||% NULL
    }

    # ============================================================
    # DE RESULTS (5 keys)
    # Note: Keys already initialized to NULL by init_shiny_payload()
    # Do NOT use payload$key <- NULL here as it REMOVES the key!
    # ============================================================

    if (!is.null(de_res)) {
        # de_model: Statistical model object (limma fit, etc.)
        de_model_val <- de_res$de_model %||% de_res$fit
        if (!is.null(de_model_val)) payload$de_model <- de_model_val

        # de_stats: Full DE statistics table
        de_stats_val <- de_res$summary_df %||% de_res$de_table
        if (!is.null(de_stats_val)) payload$de_stats <- de_stats_val

        # Ensure feature_id column exists
        if (!is.null(payload$de_stats)) {
            possible_id_cols <- c("feature_id", "FeatureID", "metabolite", "compound", "name")
            found_col <- intersect(possible_id_cols, names(payload$de_stats))[1]

            if (!is.na(found_col) && !"feature_id" %in% names(payload$de_stats)) {
                payload$de_stats$feature_id <- payload$de_stats[[found_col]]
            }
        }

        # de_sig_stats: Subset of de_stats for significant features
        if (!is.null(payload$de_stats) && "pass_any_contrast" %in% colnames(payload$de_stats)) {
            payload$de_sig_stats <- payload$de_stats[
                !is.na(payload$de_stats$pass_any_contrast) &
                payload$de_stats$pass_any_contrast == 1, ,
                drop = FALSE
            ]
        }

        # de_expr_norm: Expression matrix subset to significant features
        if (!is.null(payload$expr_norm) && !is.null(payload$de_sig_stats) && nrow(payload$de_sig_stats) > 0) {
            sig_ids <- payload$de_sig_stats$feature_id %||%
                       payload$de_sig_stats$FeatureID %||%
                       payload$de_sig_stats$metabolite

            if (!is.null(sig_ids) && length(sig_ids) > 0) {
                matched_ids <- intersect(sig_ids, rownames(payload$expr_norm))
                if (length(matched_ids) > 0) {
                    payload$de_expr_norm <- payload$expr_norm[matched_ids, , drop = FALSE]
                }
            }
        }

        # de_summary: Per-contrast summary counts
        if (!is.null(payload$de_stats)) {
            payload$de_summary <- build_de_summary_counts_metabolomics(payload$de_stats, out_dir = out_dir)
        }

        # de_final_table: DE-significant rows (equivalent to Final_results_DE_P_*.xlsx)
        if (!is.null(payload$de_sig_stats) && nrow(payload$de_sig_stats) > 0) {
            payload$de_final_table <- payload$de_sig_stats
        }
    }

    # ============================================================
    # CLUSTERING (4 keys)
    # Note: Keys already initialized to NULL by init_shiny_payload()
    # Do NOT use payload$key <- NULL here as it REMOVES the key!
    # ============================================================

    if (!is.null(clustering_res)) {
        src <- if (!is.null(clustering_res$objects)) clustering_res$objects else clustering_res

        val <- src$clusters %||% src$New_clusters
        if (!is.null(val)) payload$clust_partition <- val
        if (!is.null(src$patterns)) payload$clust_patterns <- src$patterns
        val <- src$heatmaps %||% src$heatmaps_by_pattern
        if (!is.null(val)) payload$clust_heatmaps_by_pattern <- val

        # clust_heatmap_partition: Partition clustering heatmap (full pheatmap object)
        if (!is.null(clustering_res$plots$partition_heatmap))
            payload$clust_heatmap_partition <- clustering_res$plots$partition_heatmap

        # clust_heatmap_partition_fig: The actual drawable gtable (print to see the plot)
        if (!is.null(clustering_res$plots$partition_heatmap) && !is.null(clustering_res$plots$partition_heatmap$gtable))
            payload$clust_heatmap_partition_fig <- clustering_res$plots$partition_heatmap$gtable
    }

    # clust_heatmap_hier: Hierarchical clustering (pheatmap + dendrogram)
    # Priority: clustering_res$objects$hm_hier_de (full payload with tree_row, row_order, etc.)
    val <- NULL
    if (!is.null(clustering_res)) {
      # First check objects (where the full payload is stored)
      if (!is.null(clustering_res$objects)) {
        val <- clustering_res$objects$hm_hier_de
      }
      # Fallback to plots if needed
      if (is.null(val) && !is.null(clustering_res$plots)) {
        val <- clustering_res$plots$hm_hier_de %||% clustering_res$plots$p_cluster_hier
      }
    }
    if (!is.null(val)) payload$clust_heatmap_hier <- val
    
    # clust_heatmap_hier_fig: The actual drawable gtable from pheatmap (print to see the plot)
    if (!is.null(val) && !is.null(val$pheatmap) && !is.null(val$pheatmap$gtable))
      payload$clust_heatmap_hier_fig <- val$pheatmap$gtable

    # ============================================================
    # CONFIGURATION (6 keys)
    # ============================================================

    # Canonical key names (overwrite init_shiny_payload defaults)
    payload$padj_cutoff <- de_cfg$padj_cutoff %||% de_cfg$p_cutoff %||% 0.05

    linear_fc <- de_cfg$linear_fc_cutoff %||% 1.5
    payload$log_fc_cutoff <- log2(linear_fc)
    payload$fc_cutoff <- linear_fc
    # Build normalization method description
    norm_method <- paste0(
        norm_cfg$sample_norm %||% "none",
        "/",
        norm_cfg$transform %||% "none",
        "/",
        norm_cfg$scaling %||% "none"
    )
    payload$norm_method <- norm_method

    # Group and aesthetic variables (canonical: group, color)
    if (!is.null(effects_cfg$color)) {
        payload$group <- as.character(effects_cfg$color[[1]])
    } else if (!is.null(effects_cfg$group)) {
        payload$group <- as.character(effects_cfg$group[[1]])
    } else {
        payload$group <- "sample_type"
    }

    if (!is.null(effects_cfg$color)) {
        payload$color <- as.character(effects_cfg$color)
    }

    if (!is.null(effects_cfg$shape)) {
        payload$shape <- as.character(effects_cfg$shape)
    }

    # ============================================================
    # VALIDATION
    # ============================================================

    assert_shiny_payload_contract(payload, strict = FALSE, context = "metabolomics")

    # ============================================================
    # LEGACY ALIASES (optional)
    # ============================================================

    if (isTRUE(include_legacy)) {
        payload <- attach_legacy_aliases(payload)

        # Metabolomics-specific legacy keys
        payload$metab_raw <- pre$expr_raw %||% NULL
        payload$metab_filt <- pre$expr_filt %||% NULL
        payload$metab_norm <- pre$expr_work %||% NULL

        # Missingness info (useful for metabolomics QC)
        if (!is.null(pre$info) && !is.null(pre$info$missingness)) {
            payload$missingness <- pre$info$missingness
        }

        # Sample map (CD column to sample_id mapping)
        if (!is.null(pre$sample_map)) {
            payload$sample_map <- pre$sample_map
        }

        # Row data (feature annotations) for compatibility
        payload$row_data <- pre$row_data %||% NULL

        # Random forest results
        if (!is.null(rf_res)) {
            payload$rf_importance <- rf_res$importance_df %||% NULL
            payload$rf_method     <- rf_res$method %||% NULL
        }

        # PLS-DA results
        if (!is.null(plsda_res)) {
            payload$plsda_vip_df           <- plsda_res$vip_df %||% NULL
            payload$plsda_explained_variance <- plsda_res$explained_variance %||% NULL
        }

        # Enrichment results
        if (!is.null(enrichment_res)) {
            if (!is.null(enrichment_res$qea)) {
                payload$enrichment_qea <- enrichment_res$qea$table %||% NULL
            }
            if (!is.null(enrichment_res$ssgsea)) {
                payload$enrichment_ssgsea       <- enrichment_res$ssgsea$table %||% NULL
                payload$enrichment_ssgsea_scores <- enrichment_res$ssgsea$scores %||% NULL
            }
            if (!is.null(enrichment_res$ora)) {
                payload$enrichment_ora <- enrichment_res$ora$table %||% NULL
            }
            if (!is.null(enrichment_res$gsea)) {
                payload$enrichment_gsea <- enrichment_res$gsea$table %||% NULL
            }
        }
    }

    # ============================================================
    # SUMMARY
    # ============================================================

    message("Metabolomics payload created successfully")
    message("  - Canonical keys: ", length(get_canonical_keys()))
    message("  - Total keys: ", length(payload))
    message("  - NULL keys: ", sum(sapply(payload, is.null)))
    message("  - Non-NULL keys: ", sum(!sapply(payload, is.null)))

    payload
}


#' Build DE summary counts for Metabolomics
#'
#' Thin wrapper around \code{\link{build_de_summary_counts_generic}} with
#' metabolomics naming conventions: \code{pass.<contrast>} pass columns and
#' \code{linearFC.<contrast>} fold-change columns.
#'
#' @param de_stats DE statistics data.frame with pass columns
#' @param out_dir Optional: output directory to write TSV file. If provided, writes de_summary.tsv
#' @return data.frame with columns: contrast, up, down, total (invisibly if file written)
#' @keywords internal
build_de_summary_counts_metabolomics <- function(de_stats, out_dir = NULL) {
    result <- build_de_summary_counts_generic(
        de_stats         = de_stats,
        pass_pattern     = "^pass\\.",
        extract_contrast = function(col) sub("^pass\\.", "", col),
        find_fc_col      = function(cn, cols) {
            fc <- paste0("linearFC.", cn)
            if (fc %in% cols) fc else NULL
        }
    )

    if (!is.null(out_dir) && !is.null(result) && nrow(result) > 0) {
        save_tsv(result, out_dir, "de_summary_counts.tsv")
    }

    result
}


#' Save Metabolomics Shiny payload to RDS file
#'
#' @param ... Arguments passed to build_shiny_payload_metabolomics()
#' @param out_file Output file path
#' @return Path to saved file (invisibly)
#' @export
save_shiny_payload_metabolomics <- function(..., out_file = "shiny_payload_metabolomics.rds") {
    out_dir <- dirname(out_file)
    payload <- build_shiny_payload_metabolomics(..., out_dir = out_dir)
    if (nchar(out_dir) > 0 && !dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    saveRDS(payload, out_file)
    message("Saved metabolomics payload to: ", out_file)

    invisible(out_file)
}


#' Load and validate metabolomics Shiny payload
#'
#' @param file Path to RDS file
#' @param verbose Print diagnostic information
#' @return The loaded payload list
#' @export
load_shiny_payload_metabolomics <- function(file, verbose = TRUE) {
    if (!file.exists(file)) {
        stop("File not found: ", file)
    }

    payload <- readRDS(file)

    if (verbose) {
        message("Loaded metabolomics payload from: ", file)
        message("  Version: ", payload$payload_version  %||% "unknown")
        message("  Created: ", payload$payload_created_at  %||% "unknown")
        message("  Source: ", payload$payload_source  %||% "unknown")
        message("  Total keys: ", length(payload))
        message("  NULL keys: ", sum(sapply(payload, is.null)))
        message("  Non-NULL keys: ", sum(!sapply(payload, is.null)))

        # Check critical keys
        critical <- c("sample_meta", "expr_raw", "expr_norm", "de_stats", "pca_object")
        for (key in critical) {
            if (key %in% names(payload)) {
                status <- if (is.null(payload[[key]])) "NULL" else "OK"
                message("  [", status, "] ", key)
            }
        }

        # Report expression dimensions
        if (!is.null(payload$expr_norm)) {
            message("  Expression: ", nrow(payload$expr_norm), " features x ",
                    ncol(payload$expr_norm), " samples")
        }
    }

    payload
}
