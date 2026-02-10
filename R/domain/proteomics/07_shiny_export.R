# ============================================================
# Shiny Export — Proteomics
# ============================================================
# This file contains both:
# 1. CANONICAL builder (v2.0): build_shiny_payload_proteomics()
# 2. LEGACY builder (deprecated): build_data_to_shiny_legacy_proteomics()
#
# New code should use build_shiny_payload_proteomics().
# Legacy builder is kept for backward compatibility during transition.
# ============================================================


# ============================================================
# CANONICAL BUILDER (v2.0)
# ============================================================

#' Build canonical Shiny payload for Proteomics
#'
#' Creates a Shiny payload conforming to the canonical contract (v2.0).
#' All 26 keys are guaranteed to exist (NULL if not applicable).
#'
#' @param pre Preprocessing results (from preprocess_proteomics)
#' @param de_res DE results (from proteomics DE analysis). Can be NULL if DE was skipped.
#' @param inputs Input list (contrasts, metadata, etc.)
#' @param config Full config object
#' @param pca_res Optional: pre-computed PCA results
#' @param clustering_res Optional: pre-computed clustering results
#' @param annot Optional: external annotation data.frame
#' @param include_legacy Logical. If TRUE (default), attach legacy aliases.
#'
#' @return A named list with 26 canonical keys (+ legacy aliases if requested)
#'
#' @export
build_shiny_payload_proteomics <- function(
    pre,
    de_res = NULL,
    inputs,
    config,
    pca_res = NULL,
    clustering_res = NULL,
    annot = NULL,
    include_legacy = TRUE
) {
    # ============================================================
    # Initialize canonical payload structure
    # ============================================================
    payload <- init_shiny_payload("proteomics")

    # ============================================================
    # Extract config shortcuts
    # ============================================================
    modes <- config$modes %||% list()
    prot_cfg <- modes$proteomics %||% modes$prot %||% list()
    de_cfg <- prot_cfg$de %||% list()
    norm_cfg <- prot_cfg$normalization %||% list()
    effects_cfg <- prot_cfg[["effects"]] %||% list()

    # ============================================================
    # METADATA (3 keys)
    # ============================================================

    # sample_meta: Sample metadata with rownames as sample IDs
    payload$sample_meta <- pre$meta
    if (!is.null(payload$sample_meta)) {
        sample_col <- effects_cfg$samples %||% NULL
        if (!is.null(sample_col) && sample_col %in% colnames(payload$sample_meta)) {
            rownames(payload$sample_meta) <- as.character(payload$sample_meta[[sample_col]])
        }
    }

    # contrasts: Contrast definitions
    payload$contrasts <- inputs$contrasts

    # feature_annot: Feature annotations
    if (!is.null(annot)) {
        payload$feature_annot <- annot
    }

    # ============================================================
    # EXPRESSION DATA (2 keys)
    # ============================================================

    # expr_raw: Filtered expression (before imputation, may have NAs)
    payload$expr_raw <- pre$expr_filt

    # expr_norm: Imputed expression (NO NAs allowed)
    # In proteomics, we use the imputed matrix
    payload$expr_norm <- pre$expr_imp_single %||% pre$expr_work

    # ============================================================
    # QC/PCA (3 keys)
    # ============================================================

    if (!is.null(pca_res)) {
        # Handle nested structure (objects list) or flat structure
        pca_objects <- pca_res$objects %||% pca_res

        # pca_object: prcomp result
        payload$pca_object <- pca_objects$norm_log_counts_pca %||%
                              pca_objects$pca_object %||%
                              NULL

        # pca_scores: PCA scores data.frame with metadata
        payload$pca_scores <- pca_objects$pca_scores %||% NULL
    }

    # qc_plots: Named list of QC plots (proteomics has many)
    if (!is.null(pca_res) && !is.null(pca_res$plots)) {
        payload$qc_plots <- list(
            imputation_hist = pca_res$plots$imputation_hist %||% NULL,
            density = pca_res$plots$density %||% NULL,
            boxplot = pca_res$plots$boxplot %||% NULL,
            heatmap_clusters = pca_res$plots$heatmap_clusters$gtable %||% NULL,
            heatmap_nocol = pca_res$plots$heatmap_nocol$gtable %||% NULL
        )
    } else {
        payload$qc_plots <- list()
    }

    # ============================================================
    # DE RESULTS (5 keys)
    # Note: Keys already initialized to NULL by init_shiny_payload()
    # Do NOT use payload$key <- NULL here as it REMOVES the key!
    # ============================================================

    if (!is.null(de_res)) {
        # de_model: limma fit object (only assign if non-NULL)
        if (!is.null(de_res$de_model)) payload$de_model <- de_res$de_model

        # de_stats: Full DE statistics table (only assign if non-NULL)
        if (!is.null(de_res$summary_df)) {
            payload$de_stats <- de_res$summary_df

            # Ensure feature_id column exists
            possible_id_cols <- c("FeatureID", "Protein", "Protein.Group", "feature_id")
            found_col <- intersect(possible_id_cols, names(payload$de_stats))[1]

            if (!is.na(found_col) && !"feature_id" %in% names(payload$de_stats)) {
                payload$de_stats$feature_id <- payload$de_stats[[found_col]]
            }

            # de_sig_stats: Subset of de_stats for significant features
            if ("pass_any_contrast" %in% colnames(payload$de_stats)) {
                sig_df <- payload$de_stats[
                    !is.na(payload$de_stats$pass_any_contrast) &
                    payload$de_stats$pass_any_contrast == 1, ,
                    drop = FALSE
                ]
                if (nrow(sig_df) > 0) payload$de_sig_stats <- sig_df
            }

            # de_expr_norm: Expression matrix subset to significant features
            if (!is.null(payload$expr_norm) && !is.null(payload$de_sig_stats) && nrow(payload$de_sig_stats) > 0) {
                sig_ids <- payload$de_sig_stats$feature_id %||%
                           payload$de_sig_stats$FeatureID %||%
                           payload$de_sig_stats$Protein

                if (!is.null(sig_ids) && length(sig_ids) > 0) {
                    matched_ids <- intersect(sig_ids, rownames(payload$expr_norm))
                    if (length(matched_ids) > 0) {
                        payload$de_expr_norm <- payload$expr_norm[matched_ids, , drop = FALSE]
                    }
                }
            }

            # de_summary: Per-contrast summary counts
            summary_counts <- build_de_summary_counts_proteomics(payload$de_stats)
            if (!is.null(summary_counts)) payload$de_summary <- summary_counts
        }
    }

    # ============================================================
    # CLUSTERING (4 keys)
    # Note: Keys already initialized to NULL by init_shiny_payload()
    # Do NOT use payload$key <- NULL here as it REMOVES the key!
    # ============================================================

    if (!is.null(clustering_res)) {
        src <- if (!is.null(clustering_res$objects)) clustering_res$objects else clustering_res

        # Only assign if value is non-NULL (to avoid removing keys)
        val <- src$clusters %||% src$New_clusters
        if (!is.null(val)) payload$clust_partition <- val

        if (!is.null(src$patterns)) payload$clust_patterns <- src$patterns

        val <- src$heatmaps %||% src$heatmaps_by_pattern
        if (!is.null(val)) payload$clust_heatmaps_by_pattern <- val
    }

    # clust_heatmap_hier: Hierarchical clustering (pheatmap + dendrogram)
    val <- NULL
    if (!is.null(de_res)) {
        val <- de_res$pheatmap_data_DE_genes %||% de_res$pheatmap_data
    }
    if (is.null(val) && !is.null(clustering_res)) {
        src <- if (!is.null(clustering_res$objects)) clustering_res$objects else clustering_res
        val <- src$pheatmap_data_DE_genes %||% src$pheatmap_data
    }
    if (!is.null(val)) payload$clust_heatmap_hier <- val

    # ============================================================
    # CONFIGURATION (6 keys)
    # ============================================================

    payload$padj_cutoff <- de_cfg$padj_cutoff %||% de_cfg$p_cutoff %||% 0.05

    linear_fc <- de_cfg$linear_fc_cutoff %||% 1.5
    payload$fc_cutoff <- linear_fc
    payload$log_fc_cutoff <- log2(linear_fc)

    payload$norm_method <- norm_cfg$method %||% "none"

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

    assert_shiny_payload_contract(payload, strict = FALSE, context = "proteomics")

    # ============================================================
    # LEGACY ALIASES (optional)
    # ============================================================

    if (isTRUE(include_legacy)) {
        payload <- attach_legacy_aliases(payload)

        # Proteomics-specific legacy keys
        payload$prot_log2_raw <- pre$expr_filt %||% NULL
        payload$prot_log2_imp1 <- pre$expr_imp_single %||% NULL

        # QC plot legacy aliases
        if (!is.null(pca_res) && !is.null(pca_res$plots)) {
            payload$imp_hist_samp <- pca_res$plots$imputation_hist %||% NULL
            payload$prot_hist <- pca_res$plots$density %||% NULL
            payload$imp_box <- pca_res$plots$boxplot %||% NULL
            payload$prot_hm <- pca_res$plots$heatmap_clusters %||% NULL
            payload$prot_hm_noNA <- pca_res$plots$heatmap_nocol %||% NULL
        }

        # DE contrast summary (proteomics-specific)
        if (!is.null(payload$de_stats)) {
            payload$de_contrast_summary <- build_de_contrast_summary(payload$de_stats)
            payload$de_summary_counts <- payload$de_summary
        }

        # Ensure Protein column exists in stats for legacy compatibility
        if (!is.null(payload$stats_df) && "feature_id" %in% names(payload$stats_df)) {
            if (!"Protein" %in% names(payload$stats_df)) {
                payload$stats_df$Protein <- payload$stats_df$feature_id
            }
        }
    }

    # ============================================================
    # SUMMARY
    # ============================================================

    message("Proteomics payload created successfully")
    message("  - Canonical keys: ", length(get_canonical_keys()))
    message("  - Total keys: ", length(payload))
    message("  - NULL keys: ", sum(sapply(payload, is.null)))
    message("  - Non-NULL keys: ", sum(!sapply(payload, is.null)))

    payload
}


#' Build DE summary counts for Proteomics
#'
#' @param de_stats DE statistics data.frame with pass columns
#' @return data.frame with columns: contrast, up, down, total
#' @keywords internal
build_de_summary_counts_proteomics <- function(de_stats) {
    if (is.null(de_stats)) return(NULL)

    pass_cols <- grep("^pass_", names(de_stats), value = TRUE)
    pass_cols <- setdiff(pass_cols, "pass_any_contrast")

    if (length(pass_cols) == 0) return(NULL)

    summaries <- lapply(pass_cols, function(col) {
        contrast_name <- sub("^pass_", "", col)

        # Try multiple FC column name patterns
        fc_patterns <- c(
            paste0("logFC_", contrast_name),
            paste0("log2FoldChange_", contrast_name),
            paste0("logFC.", contrast_name)
        )
        fc_col <- NULL
        for (pat in fc_patterns) {
            if (pat %in% names(de_stats)) {
                fc_col <- pat
                break
            }
        }

        sig_rows <- !is.na(de_stats[[col]]) & de_stats[[col]] == 1

        if (!is.null(fc_col)) {
            up <- sum(sig_rows & de_stats[[fc_col]] > 0, na.rm = TRUE)
            down <- sum(sig_rows & de_stats[[fc_col]] < 0, na.rm = TRUE)
        } else {
            up <- NA
            down <- NA
        }

        data.frame(
            contrast = contrast_name,
            up = up,
            down = down,
            total = sum(sig_rows, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, summaries)
}


#' Build DE contrast summary (proteomics-specific)
#'
#' @param de_stats DE statistics data.frame
#' @return data.frame with contrast summary
#' @keywords internal
build_de_contrast_summary <- function(de_stats) {
    if (is.null(de_stats)) return(NULL)

    # Find pass columns
    pass_cols <- grep("^pass_", names(de_stats), value = TRUE)
    pass_cols <- setdiff(pass_cols, "pass_any_contrast")

    if (length(pass_cols) == 0) return(NULL)

    summaries <- lapply(pass_cols, function(col) {
        contrast_name <- sub("^pass_", "", col)
        sig_count <- sum(!is.na(de_stats[[col]]) & de_stats[[col]] == 1, na.rm = TRUE)

        data.frame(
            Name = contrast_name,
            any = sig_count,
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, summaries)
}


#' Save Proteomics Shiny payload to RDS file
#'
#' @param ... Arguments passed to build_shiny_payload_proteomics()
#' @param out_file Output file path
#' @return Path to saved file (invisibly)
#' @export
save_shiny_payload_proteomics <- function(..., out_file = "shiny_payload_proteomics.rds") {
    payload <- build_shiny_payload_proteomics(...)

    out_dir <- dirname(out_file)
    if (nchar(out_dir) > 0 && !dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    saveRDS(payload, out_file)
    message("Saved proteomics payload to: ", out_file)

    invisible(out_file)
}


# ============================================================
# LEGACY BUILDER (DEPRECATED)
# ============================================================

#' Build legacy-compatible data structure for Shiny app (Proteomics)
#'
#' @description
#' DEPRECATED: Use build_shiny_payload_proteomics() instead.
#'
#' This function is kept for backward compatibility only.
#' It will be removed in a future version.
#'
#' @param pre Preprocessing results (from preprocess_proteomics)
#' @param de_res DE results (from proteomics DE analysis)
#' @param inputs Input list (contrasts, metadata, etc.)
#' @param config Full config object
#' @param pca_res Optional: pre-computed PCA results
#' @param clustering_res Optional: pre-computed clustering results
#' @param annot Optional: external annotation data.frame
#'
#' @return A flat named list with ALL legacy keys (even if NULL)
#'
#' @export
build_data_to_shiny_legacy_proteomics <- function(
    pre,
    de_res,
    inputs,
    config,
    pca_res = NULL,
    clustering_res = NULL,
    annot = NULL
) {
    .Deprecated("build_shiny_payload_proteomics")

    legacy <- build_shiny_legacy_base(
        legacy_source = "Proteomics pipeline",
        pca_basename = "prot_pca_3d"
    )

    # Ensure color/shape keys exist
    if (!("color" %in% names(legacy))) legacy["color"] <- list(NULL)
    if (!("shape" %in% names(legacy))) legacy["shape"] <- list(NULL)

    # Metadata
    legacy$col_data <- pre$meta %||% NULL

    prot_fx <- config$modes$proteomics[["effects"]] %||% list()
    sample_col <- prot_fx$samples %||% NULL
    if (!is.null(legacy$col_data) && !is.null(sample_col) && sample_col %in% colnames(legacy$col_data)) {
        rownames(legacy$col_data) <- as.character(legacy$col_data[[sample_col]])
    }

    legacy$contrasts_data <- inputs$contrasts %||% NULL

    # Expression matrices
    legacy$norm_counts <- pre$expr_filt %||% NULL
    legacy$norm_log_counts <- pre$expr_work %||% NULL

    # PCA objects
    if (!is.null(pca_res) && !is.null(pca_res$objects)) {
        legacy$norm_log_counts_pca <- pca_res$objects$norm_log_counts_pca %||% NULL
        legacy$color <- pca_res$objects$color %||% NULL
        legacy$shape <- pca_res$objects$shape %||% NULL
    }

    # Fallback for color/shape
    val_color <- legacy[["color"]] %||% prot_fx[["color"]]
    val_shape <- legacy[["shape"]] %||% prot_fx[["shape"]]
    legacy["color"] <- list(val_color)
    legacy["shape"] <- list(val_shape)

    # DE summary and stats
    legacy$stats_df <- de_res$summary_df %||% NULL

    legacy["DE_genes_stats"] <- list(NULL)
    if (!is.null(legacy$stats_df) && "pass_any_contrast" %in% colnames(legacy$stats_df)) {
        legacy$DE_genes_stats <- legacy$stats_df[
            !is.na(legacy$stats_df$pass_any_contrast) & legacy$stats_df$pass_any_contrast == 1, ,
            drop = FALSE
        ]
    }

    # Heatmap data matrix
    legacy["mat2plot"] <- list(NULL)
    legacy["de_expr_norm"] <- list(NULL)
    if (!is.null(legacy$norm_log_counts) && !is.null(legacy$DE_genes_stats) && nrow(legacy$DE_genes_stats) > 0) {
        prot_cfg <- config$modes$proteomics %||% list()
        de_table_cfg <- prot_cfg$de_table %||% list()
        id_col <- de_table_cfg$id_col %||% "FeatureID"

        possible_cols <- c(id_col, "Protein", "FeatureID", "Protein.Group")
        id_col_found <- intersect(possible_cols, colnames(legacy$DE_genes_stats))[1]

        if (!is.na(id_col_found)) {
            sig_ids <- legacy$DE_genes_stats[[id_col_found]]
            de_matrix <- legacy$norm_log_counts[rownames(legacy$norm_log_counts) %in% sig_ids, , drop = FALSE]
            legacy$mat2plot <- assert_de_expr_matrix(de_matrix, context = "proteomics export")
            legacy$de_expr_norm <- legacy$mat2plot
        }
    }

    # DE Model
    legacy$de_model <- de_res$de_model %||% NULL

    # DE summaries
    legacy["de_contrast_summary"] <- list(NULL)
    legacy["de_summary_counts"] <- list(NULL)
    if (!is.null(legacy$stats_df)) {
        legacy$de_contrast_summary <- build_de_contrast_summary(legacy$stats_df)
        legacy$de_summary_counts <- build_de_summary_counts_proteomics(legacy$stats_df)
    }

    # Heatmap data
    de_heat <- de_res$pheatmap_data_DE_genes %||% de_res$pheatmap_data %||% NULL
    clust_heat <- if (!is.null(clustering_res)) {
        clustering_res$pheatmap_data_DE_genes %||% clustering_res$pheatmap_data %||% NULL
    } else {
        NULL
    }
    legacy$pheatmap_data_DE_genes <- de_heat %||% clust_heat %||% NULL

    # Clustering results
    if (!is.null(clustering_res)) {
        legacy$patterns <- clustering_res$patterns %||% NULL
        legacy$heatmaps_by_pattern <- clustering_res$heatmaps %||% NULL
        legacy$New_clusters <- clustering_res$clusters %||% NULL
    }

    # Annotation
    legacy$annot <- annot

    # Config parameters
    modes <- config$modes %||% list()
    prot <- modes$proteomics %||% modes$prot %||% list()
    de_cfg <- prot$de %||% list()
    norm_cfg <- prot$normalization %||% list()

    legacy$PADJ_CUTOFF <- de_cfg$padj_cutoff %||% de_cfg$p_cutoff %||% 0.05
    legacy$DESEQ_PADJ_CUTOFF <- legacy$PADJ_CUTOFF

    linear_fc <- de_cfg$linear_fc_cutoff %||% 1.5
    legacy$LINEAR_FC_CUTOFF <- linear_fc
    legacy$LOG_FC_CUTOFF <- log2(linear_fc)

    legacy$NORM_METHOD <- norm_cfg$method %||% "none"

    prot_effects_cfg <- prot[["effects"]] %||% list()
    legacy$GROUP <- prot_effects_cfg$color %||% NULL

    # QC plot objects
    if (!is.null(pca_res) && !is.null(pca_res$plots)) {
        legacy$imp_hist_samp <- pca_res$plots$imputation_hist %||% NULL
        legacy$prot_hist <- pca_res$plots$density %||% NULL
        legacy$imp_box <- pca_res$plots$boxplot %||% NULL
        legacy$prot_hm <- pca_res$plots$heatmap_clusters %||% NULL
        legacy$prot_hm_noNA <- pca_res$plots$heatmap_nocol %||% NULL
    } else {
        legacy["imp_hist_samp"] <- list(NULL)
        legacy["prot_hist"] <- list(NULL)
        legacy["imp_box"] <- list(NULL)
        legacy["prot_hm"] <- list(NULL)
        legacy["prot_hm_noNA"] <- list(NULL)
    }

    # Expression matrices (proteomics-specific)
    legacy$prot_log2_raw <- pre$expr_filt %||% NULL
    legacy$prot_log2_imp1 <- pre$expr_imp_single %||% NULL

    # Ensure expected keys exist
    expected_keys <- c(
        "legacy_version", "legacy_created_at", "legacy_source",
        "col_data", "contrasts_data", "dds", "norm_counts", "norm_log_counts",
        "de_model", "mat2plot", "de_expr_norm",
        "PCA_3D_BASENAME",
        "stats_df", "DE_genes_stats",
        "trinotate_main",
        "PADJ_CUTOFF", "DESEQ_PADJ_CUTOFF", "LOG_FC_CUTOFF",
        "LINEAR_FC_CUTOFF", "NORM_METHOD", "GROUP",
        "color", "shape"
    )

    for (k in expected_keys) {
        if (!k %in% names(legacy)) {
            legacy[k] <- list(NULL)
        }
    }

    legacy
}


#' Save proteomics legacy data structure to RDS file
#'
#' @description DEPRECATED: Use save_shiny_payload_proteomics() instead.
#'
#' @param ... Arguments passed to build_data_to_shiny_legacy_proteomics()
#' @param out_file Output file path
#' @return Path to saved file (invisibly)
#' @export
save_data_to_shiny_legacy_proteomics <- function(..., out_file = "data_to_shiny_legacy_proteomics.rds") {
    .Deprecated("save_shiny_payload_proteomics")

    legacy <- build_data_to_shiny_legacy_proteomics(...)

    out_dir <- dirname(out_file)
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    saveRDS(legacy, out_file)
    message("Saved proteomics legacy export to: ", out_file)

    invisible(out_file)
}
