# ============================================================
# Shiny Export — RNA-seq
# ============================================================
# This file contains both:
# 1. CANONICAL builder (v2.0): build_shiny_payload_rnaseq()
#
# New code should use build_shiny_payload_rnaseq().
# ============================================================


# ============================================================
# CANONICAL BUILDER (v2.0)
# ============================================================

#' Build canonical Shiny payload for RNA-seq
#'
#' Creates a Shiny payload conforming to the canonical contract (v2.0).
#' All 26 keys are guaranteed to exist (NULL if not applicable).
#'
#' @param pre Preprocessing results (from preprocess_rna)
#' @param de_res DE results (from run_deseq2_de). Can be NULL if DE was skipped.
#' @param inputs Input list (contrasts, metadata, etc.)
#' @param config Full config object
#' @param pca_res Optional: pre-computed PCA results
#' @param clustering_res Optional: pre-computed clustering results
#' @param annot Optional: external annotation data.frame
#' @param trinotate_main Optional: Trinotate annotation data.frame (for legacy compatibility)
#'
#' @return A named list with 26 canonical keys (+ legacy aliases if requested)
#'
#' @export
build_shiny_payload_rnaseq <- function(
    pre,
    de_res = NULL,
    inputs,
    config,
    pca_res = NULL,
    clustering_res = NULL,
    annot = NULL,
    trinotate_main = NULL,
    out_dir = NULL
) {
    # ============================================================
    # Initialize canonical payload structure
    # ============================================================
    payload <- init_shiny_payload("rnaseq")

    # ============================================================
    # Extract config shortcuts
    # ============================================================
    modes <- config$modes %||% list()
    rna_cfg <- modes$rna %||% list()
    de_cfg <- rna_cfg$de %||% list()
    norm_cfg <- rna_cfg$normalization %||% list()
    effects_cfg <- rna_cfg$effects %||% list()

    # ============================================================
    # METADATA (3 keys)
    # ============================================================

    # sample_meta: Sample metadata with rownames as sample IDs
    payload$sample_meta <- pre$meta
    if (!is.null(payload$sample_meta)) {
        # Ensure rownames are set (required by contract)
        sample_col <- effects_cfg$samples %||% NULL
        if (!is.null(sample_col) && sample_col %in% colnames(payload$sample_meta)) {
            rownames(payload$sample_meta) <- as.character(payload$sample_meta[[sample_col]])
        }
    }

    # contrasts: Contrast definitions
    payload$contrasts <- inputs$contrasts

    # feature_annot: Feature annotations (gene symbols, descriptions)
    if (!is.null(annot)) {
        payload$feature_annot <- annot
    }

    # ============================================================
    # EXPRESSION DATA (2 keys)
    # ============================================================

    # expr_raw: Filtered expression (before normalization)
    payload$expr_raw <- pre$expr_filt

    # expr_norm: Normalized expression (NO NAs allowed)
    payload$expr_norm <- pre$expr_work

    # expr_long: Long-format expression with metadata
    payload$expr_long <- build_expr_long(payload$expr_norm, payload$sample_meta)

    # ============================================================
    # QC/PCA (3 keys)
    # ============================================================

    if (!is.null(pca_res)) {
        # pca_object: prcomp result
        pca_obj <- pca_res$norm_log_counts_pca %||% pca_res$pca_object
        if (!is.null(pca_obj)) payload$pca_object <- pca_obj

        # pca_scores: PCA scores data.frame with metadata
        pca_sc <- pca_res$pca_scores %||% pca_res$objects$pca_scores
        if (!is.null(pca_sc)) payload$pca_scores <- pca_sc

        # pca_3d: 3D PCA plotly widget (only set if non-NULL; assigning NULL removes list key)
        if (!is.null(pca_res$plots$pca_3d)) payload$pca_3d <- pca_res$plots$pca_3d

        # QC plot
        if (!is.null(pca_res$plots$dist_heatmap)) payload$samples_hm <- pca_res$plots$dist_heatmap

    }

    # ============================================================
    # DE RESULTS (5 keys)
    # Note: Keys already initialized to NULL by init_shiny_payload()
    # Do NOT use payload$key <- NULL here as it REMOVES the key!
    # ============================================================

    if (!is.null(de_res)) {
        # de_model: DESeq2 dds object
        if (!is.null(de_res$dds)) payload$de_model <- de_res$dds

        # de_stats: Full DE statistics table
        if (!is.null(de_res$tables)) {
            payload$de_stats <- build_rnaseq_summary_df(de_res$tables, de_cfg)

            # Ensure feature_id column exists (standardize from Gene/FeatureID)
            if (!is.null(payload$de_stats)) {
                if ("FeatureID" %in% names(payload$de_stats) && !"feature_id" %in% names(payload$de_stats)) {
                    payload$de_stats$feature_id <- payload$de_stats$FeatureID
                } else if ("Gene" %in% names(payload$de_stats) && !"feature_id" %in% names(payload$de_stats)) {
                    payload$de_stats$feature_id <- payload$de_stats$Gene
                }
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
                       payload$de_sig_stats$Gene %||%
                       payload$de_sig_stats$FeatureID

            if (!is.null(sig_ids) && length(sig_ids) > 0) {
                matched_ids <- intersect(sig_ids, rownames(payload$expr_norm))
                if (length(matched_ids) > 0) {
                    payload$de_expr_norm <- payload$expr_norm[matched_ids, , drop = FALSE]
                }
            }
        }

        # de_summary: Per-contrast summary counts
        if (!is.null(payload$de_stats)) {
            payload$de_summary <- build_de_summary_counts_rnaseq(payload$de_stats, out_dir = out_dir)
        }

        # de_final_table: DE-filtered final results table (richer than de_stats)
        # Use payload$de_stats which was built above from de_res$tables
        if (!is.null(payload$de_stats) && !is.null(inputs$contrasts)) {
            final_results <- tryCatch(
                build_final_results_rnaseq(pre, payload$de_stats, inputs$contrasts, pre$row_data),
                error = function(e) {
                    warning("[shiny_export] de_final_table: ", conditionMessage(e))
                    NULL
                }
            )
            if (!is.null(final_results) && "pass_any_contrast" %in% names(final_results)) {
                is_de <- !is.na(final_results$pass_any_contrast) & final_results$pass_any_contrast == 1
                de_df <- final_results[is_de, , drop = FALSE]
                de_df <- de_df[, !startsWith(names(de_df), "manual_cutoffs") & names(de_df) != "pass_any_contrast", drop = FALSE]

                # Add order and z-score columns from clustering (same as Excel export)
                id_col <- if ("Gene" %in% names(de_df)) "Gene" else "FeatureID"
                excel_ord <- clustering_res$excel_order %||% NULL
                if (!is.null(excel_ord) && !is.null(excel_ord$ordered_ids)) {
                    de_df$order <- match(de_df[[id_col]], excel_ord$ordered_ids)
                    if (!is.null(excel_ord$zscore_mat)) {
                        zmat <- excel_ord$zscore_mat
                        idx <- match(de_df[[id_col]], rownames(zmat))
                        de_df <- cbind(de_df, zmat[idx, , drop = FALSE])
                    }
                }

                payload$de_final_table <- de_df
            }
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
        if (!is.null(src$patterns)) {
          payload$clust_patterns <- src$patterns
          payload$clust_patterns_list <- src$patterns_list}
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

    payload$padj_cutoff <- de_cfg$padj_cutoff %||% de_cfg$p_cutoff %||% 0.05

    linear_fc <- de_cfg$linear_fc_cutoff %||% 1.5
    payload$fc_cutoff <- linear_fc
    payload$log_fc_cutoff <- log2(linear_fc)

    payload$norm_method <- norm_cfg$method %||% "TMMlogCPM"

    if (!is.null(effects_cfg$color)) {
        payload$group <- as.character(effects_cfg$color[[1]])
    } else {
        payload$group <- "unknown"
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

    assert_shiny_payload_contract(payload, strict = FALSE, context = "rnaseq")

    # ============================================================
    # SUMMARY
    # ============================================================

    message("RNA-seq payload created successfully")
    message("  - Canonical keys: ", length(get_canonical_keys()))
    message("  - Total keys: ", length(payload))
    message("  - NULL keys: ", sum(sapply(payload, is.null)))
    message("  - Non-NULL keys: ", sum(!sapply(payload, is.null)))

    payload
}


#' Build DE summary counts for RNA-seq
#'
#' Thin wrapper around \code{\link{build_de_summary_counts_generic}} with
#' RNA-seq naming conventions: \code{<contrast>_pass} pass columns and
#' \code{linearFC.<contrast>} fold-change columns (with grep fallback).
#'
#' @param de_stats DE statistics data.frame with pass columns
#' @param out_dir Optional: output directory to write TSV file. If provided, writes de_summary.tsv
#' @return data.frame with columns: contrast, up, down, total (invisibly if file written)
#' @keywords internal
build_de_summary_counts_rnaseq <- function(de_stats, out_dir = NULL) {
    result <- build_de_summary_counts_generic(
        de_stats         = de_stats,
        pass_pattern     = "_pass$",
        extract_contrast = function(col) sub("_pass$", "", col),
        find_fc_col      = function(cn, cols) {
            fc <- paste0("linearFC.", cn)
            if (fc %in% cols) return(fc)
            # Fallback: grep search
            hits <- grep(paste0("linearFC.*", cn), cols, value = TRUE)
            if (length(hits) > 0) hits[1] else NULL
        }
    )

    if (!is.null(out_dir) && !is.null(result) && nrow(result) > 0) {
        save_tsv(result, out_dir, "de_summary_counts.tsv")
    }

    result
}


#' Save RNA-seq Shiny payload to RDS file
#'
#' @param ... Arguments passed to build_shiny_payload_rnaseq()
#' @param out_file Output file path
#' @return Path to saved file (invisibly)
#' @export
save_shiny_payload_rnaseq <- function(..., out_file = "shiny_payload_rnaseq.rds") {
    out_dir <- dirname(out_file)
    payload <- build_shiny_payload_rnaseq(..., out_dir = out_dir)
    if (nchar(out_dir) > 0 && !dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    saveRDS(payload, out_file)
    message("Saved RNA-seq payload to: ", out_file)

    invisible(out_file)
}


# ============================================================
# LEGACY BUILDER (DEPRECATED)
# ============================================================

#' Build legacy-compatible data structure for Shiny app
#'
#' @description
#' DEPRECATED: Use build_shiny_payload_rnaseq() instead.
#'
#' This function is kept for backward compatibility only.
#' It will be removed in a future version.
#'
#' @param pre Preprocessing results (from preprocess_rna)
#' @param de_res DE results (from run_deseq2_de)
#' @param inputs Input list (contrasts, metadata, etc.)
#' @param config Full config object
#' @param pca_res Optional: pre-computed PCA results
#' @param clustering_res Optional: pre-computed clustering results
#' @param annot Optional: external annotation data.frame
#'
#' @return A flat named list with ALL legacy keys (even if NULL)
#'
#' @export
build_data_to_shiny_legacy_rna <- function(
    pre,
    de_res,
    inputs,
    config,
    pca_res = NULL,
    clustering_res = NULL,
    annot = NULL
) {
    .Deprecated("build_shiny_payload_rnaseq")

    legacy <- build_shiny_legacy_base(
        legacy_source = "RNAseq pipeline",
        pca_basename = "rna_pca_3d"
    )

    # Metadata
    legacy$col_data <- pre$meta
    legacy$contrasts_data <- inputs$contrasts

    # Expression matrices
    legacy$norm_counts <- pre$expr_filt
    legacy$norm_log_counts <- pre$expr_work

    # DESeq2 object
    legacy$dds <- de_res$dds %||% NULL

    # PCA objects
    if (!is.null(pca_res)) {
        legacy$norm_log_counts_pca <- pca_res$norm_log_counts_pca %||% NULL
    }

    # EFFECTS
    modes <- config$modes %||% list()
    rna <- modes$rna %||% list()
    effects <- rna$effects %||% list()

    primary_color <- NULL
    if (!is.null(effects$color)) {
        primary_color <- as.character(effects$color[[1]])
    }

    if (!is.null(effects$color)) {
        # EFFECTS contains all aesthetic variables (shape if present + all colors)
        legacy$EFFECTS <- c(effects$shape, effects$color)
    } else {
        legacy$EFFECTS <- NULL
    }

    # DE summary and stats
    legacy$dds <- de_res$dds

    legacy$stats_df <- build_rnaseq_summary_df(de_res$tables, config$modes$rna$de)
    if (!is.null(legacy$stats_df) && "FeatureID" %in% names(legacy$stats_df) && !"Gene" %in% names(legacy$stats_df)) {
        names(legacy$stats_df)[names(legacy$stats_df) == "FeatureID"] <- "Gene"
    }

    if (!is.null(legacy$stats_df) && "pass_any_contrast" %in% colnames(legacy$stats_df)) {
        legacy$DE_genes_stats <- legacy$stats_df[
            !is.na(legacy$stats_df$pass_any_contrast) & legacy$stats_df$pass_any_contrast == 1, ,
            drop = FALSE
        ]
    }

    # Heatmap data
    legacy["mat2plot"] <- list(NULL)
    legacy["de_expr_norm"] <- list(NULL)
    if (!is.null(legacy$norm_log_counts) && !is.null(legacy$DE_genes_stats)) {
        sig_genes <- legacy$DE_genes_stats$Gene
        if (length(sig_genes) > 0) {
            de_matrix <- legacy$norm_log_counts[rownames(legacy$norm_log_counts) %in% sig_genes, , drop = FALSE]
            legacy$mat2plot <- assert_de_expr_matrix(de_matrix, context = "rnaseq export")
            legacy$de_expr_norm <- legacy$mat2plot
        }
    }

    de_heat <- de_res$pheatmap_data_DE_genes %||% de_res$pheatmap_data %||% NULL
    clust_heat <- if (!is.null(clustering_res)) {
        clustering_res$pheatmap_data_DE_genes %||% clustering_res$pheatmap_data %||% NULL
    } else {
        NULL
    }
    legacy$pheatmap_data_DE_genes <- de_heat %||% clust_heat %||% NULL

    # Clustering results
    if (!is.null(clustering_res)) {
        src <- if (!is.null(clustering_res$objects)) clustering_res$objects else clustering_res

        if (!is.null(src$patterns)) legacy$patterns <- src$patterns
        if (!is.null(src$heatmaps)) legacy$heatmaps_by_pattern <- src$heatmaps
        if (!is.null(src$clusters)) legacy$New_clusters <- src$clusters
    }

    # Annotation
    legacy$annot <- annot

    # Config parameters
    de_cfg <- rna$de %||% list()
    norm_cfg <- rna$normalization %||% list()

    legacy$PADJ_CUTOFF <- de_cfg$padj_cutoff %||% de_cfg$p_cutoff %||% 0.05
    legacy$DESEQ_PADJ_CUTOFF <- legacy$PADJ_CUTOFF

    linear_fc <- de_cfg$linear_fc_cutoff %||% 1.5
    legacy$LINEAR_FC_CUTOFF <- linear_fc
    legacy$LOG_FC_CUTOFF <- log2(linear_fc)

    legacy$NORM_METHOD <- norm_cfg$method %||% "TMMlogCPM"
    legacy$GROUP <- primary_color %||% NULL

    legacy
}


#' Save legacy data structure to RDS file
#'
#' @description DEPRECATED: Use save_shiny_payload_rnaseq() instead.
#'
#' @param ... Arguments passed to build_data_to_shiny_legacy_rna()
#' @param out_file Output file path
#' @return Path to saved file (invisibly)
#' @export
save_data_to_shiny_legacy_rna <- function(..., out_file = "data_to_shiny_legacy.rds") {
    .Deprecated("save_shiny_payload_rnaseq")

    legacy <- build_data_to_shiny_legacy_rna(...)

    out_dir <- dirname(out_file)
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    saveRDS(legacy, out_file)
    message("Saved legacy export to: ", out_file)

    invisible(out_file)
}


#' Load and validate legacy data structure
#'
#' @param file Path to legacy RDS file
#' @param verbose Print diagnostic information
#' @return The loaded legacy list
#' @export
load_legacy_rds <- function(file, verbose = TRUE) {
    if (!file.exists(file)) {
        stop("File not found: ", file)
    }

    legacy <- readRDS(file)

    if (verbose) {
        message("Loaded legacy export from: ", file)
        message("  Version: ", legacy$legacy_version %||% legacy$payload_version %||% "unknown")
        message("  Created: ", legacy$legacy_created_at %||% legacy$payload_created_at %||% "unknown")
        message("  Source: ", legacy$legacy_source %||% legacy$payload_source %||% "unknown")
        message("  Total keys: ", length(legacy))
        message("  NULL keys: ", sum(sapply(legacy, is.null)))
        message("  Non-NULL keys: ", sum(!sapply(legacy, is.null)))

        critical <- c("col_data", "sample_meta", "norm_log_counts", "expr_norm", "stats_df", "de_stats")
        for (key in critical) {
            if (key %in% names(legacy)) {
                status <- if (is.null(legacy[[key]])) "NULL" else "OK"
                message("  [", status, "] ", key)
            }
        }
    }

    legacy
}
