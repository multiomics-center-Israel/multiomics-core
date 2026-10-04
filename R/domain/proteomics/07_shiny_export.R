# ============================================================
# Shiny Export — Proteomics
# ============================================================
# CANONICAL builder (v2.0): build_shiny_payload_proteomics()
# ============================================================


# ============================================================
# CANONICAL BUILDER (v2.0)
# ============================================================

#' Build canonical Shiny payload for Proteomics
#'
#' Creates a Shiny payload conforming to the canonical contract (v2.0).
#' Every canonical key is guaranteed to exist (NULL if not applicable).
#'
#' @param pre Preprocessing results (from preprocess_proteomics)
#' @param de_res DE results (from proteomics DE analysis). Can be NULL if DE was skipped.
#' @param inputs Input list (contrasts, metadata, etc.)
#' @param config Full config object
#' @param pca_res Optional: pre-computed PCA results
#' @param clustering_res Optional: pre-computed clustering results
#' @param annot Optional: external annotation data.frame
#' @param final_results Optional: the final-results data.frame (as built by
#'   \code{build_final_results_proteomics()}), the source of \code{de_final_table}.
#'   NULL leaves \code{de_final_table} NULL.
#' @return A named list containing every canonical key (NULL where not applicable),
#'   plus any omics-specific extension keys.
#'
#' @export
build_shiny_payload_proteomics <- function(
    pre,
    de_res = NULL,
    inputs,
    config,
    pca_res = NULL,
    clustering_res = NULL,
    final_results = NULL,
    xlsx_files = NULL,
    out_dir = NULL,
    annot = NULL)
  {
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
    # METADATA
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
    payload["contrasts"] <- list(inputs$contrasts)

    # feature_annot: Feature annotations (Protein.Names, Genes, descriptions).
    # An explicit `annot` arg overrides; otherwise derive from pre$row_data via
    # the shared helper so future annotation columns flow through automatically.
    payload["feature_annot"] <- list(if (!is.null(annot)) {
        annot
    } else {
        build_feature_annot(pre$row_data, prot_cfg$id_columns$protein_id)
    })

    # ============================================================
    # EXPRESSION DATA
    # ============================================================

    # expr_raw: Filtered expression (before imputation, may have NAs)
    payload$expr_raw <- pre$expr_filt

    # expr_norm: Imputed expression (NO NAs allowed)
    # In proteomics, we use the imputed matrix
    payload$expr_norm <- pre$expr_imp_single %||% pre$expr_work
    
    # expr_long: Long-format expression with metadata
    payload["expr_long"] <- list(build_expr_long(payload$expr_norm, payload$sample_meta))

    # ============================================================
    # QC/PCA
    # ============================================================

    if (!is.null(pca_res)) {
        # Handle nested structure (objects list) or flat structure
        pca_objects <- pca_res$objects %||% pca_res

        # pca_object: prcomp result
        payload["pca_object"] <- list(pca_objects$norm_log_counts_pca %||%
                                      pca_objects$pca_object)

        # pca_scores: PCA scores data.frame with metadata
        payload["pca_scores"] <- list(pca_objects$pca_scores)

        # pca_3d: 3D PCA plotly widget
        payload["pca_3d"] <- list(pca_res$plots$pca_3d)

        # QC plot
        payload["imp_hist_samp"] <- list(pca_res$plots$imputation_hist)
        payload["samples_hm"] <- list(pca_res$plots$dist_heatmap)
        payload["samples_hm_w_na"] <- list(pca_res$plots$dist_heatmap_na)

    }

    # ============================================================
    # DE RESULTS
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

            # Merge row_data annotation columns into de_stats. on_collision = "skip"
            # is critical: Protein.Names / Genes / First.Protein.Description are already
            # carried in summary_df by summarize_limma_mult_imputation() (and may have
            # been edited by apply_custom_annotation()), so summary_df wins.
            if (!is.null(pre$row_data) && ncol(pre$row_data) > 0) {
                protein_id_col <- prot_cfg$id_columns$protein_id
                payload$de_stats <- annotate_de_stats(
                    payload$de_stats,
                    pre$row_data,
                    id_col_de    = "feature_id",
                    id_col_row   = protein_id_col,
                    on_collision = "skip"
                )
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
            contrasts_vec <- if (!is.null(inputs$contrasts) &&
                                 "Contrast_name" %in% colnames(inputs$contrasts)) {
                as.character(inputs$contrasts$Contrast_name)
            } else {
                character(0)
            }
            summary_counts <- build_de_summary_counts_proteomics(
                payload$de_stats,
                contrasts = contrasts_vec,
                out_dir   = out_dir
            )
            if (!is.null(summary_counts)) payload$de_summary <- summary_counts
        }

        # de_final_table: clean DE table built from the final-results data.frame
        # (not from the workbook, whose Results sheet is a presentation layout).
        # Same construction as the RNA-seq and metabolomics payloads: DE rows
        # only, cutoff/pass helper columns dropped, clustering order + z-scores
        # appended.
        if (!is.null(final_results)) {
            # Same key and default the export module uses to build final_results.
            id_col <- prot_cfg$de_table$id_col %||% "FeatureID"

            if (!all(c(id_col, "pass_any_contrast") %in% names(final_results))) {
                warning("[shiny_export] de_final_table: final_results has no '",
                        id_col, "' or 'pass_any_contrast' column; left NULL.")
            } else {
                is_de <- !is.na(final_results$pass_any_contrast) & final_results$pass_any_contrast == 1
                de_df <- final_results[is_de, , drop = FALSE]
                de_df <- de_df[, !startsWith(names(de_df), "manual_cutoffs") & names(de_df) != "pass_any_contrast", drop = FALSE]

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
    # Embedded xlsx bytes (all_final_xlsx, de_final_xlsx)
    # Independent of the de_final_table path above, which uses the final-results
    # data.frame; this one stores the workbooks' raw bytes.
    # ============================================================
    payload <- attach_final_results_xlsx_bytes(payload, xlsx_files)

    # ============================================================
    # CLUSTERING
    # Note: Keys already initialized to NULL by init_shiny_payload()
    # Do NOT use payload$key <- NULL here as it REMOVES the key!
    # ============================================================

    if (!is.null(clustering_res)) {
        src <- if (!is.null(clustering_res$objects)) clustering_res$objects else clustering_res

        # Only assign if value is non-NULL (to avoid removing keys)
        val <- src$clusters %||% src$New_clusters
        if (!is.null(val)) payload$clust_partition <- val

        if (!is.null(src$patterns)) {
          payload$clust_patterns <- src$patterns
          payload["clust_patterns_list"] <- list(src$patterns_list)}

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
    # CONFIGURATION
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
#' Thin wrapper around \code{\link{build_de_summary_counts_generic}} with
#' proteomics naming conventions: \code{pass.imputs.<contrast>} pass columns
#' and multiple FC column candidates (\code{logFC_}, \code{log2FoldChange_},
#' \code{logFC.}, \code{linearFC.imputs.}).
#'
#' @param de_stats  DE statistics data.frame with pass columns.
#' @param contrasts Character vector of contrast names to count.
#' @param out_dir   Optional: output directory to write TSV file.  If provided,
#'   writes de_summary_counts.tsv.
#' @return data.frame with columns: contrast, up, down, total (invisibly if
#'   file written).
#' @keywords internal
build_de_summary_counts_proteomics <- function(de_stats, contrasts, out_dir = NULL) {
    result <- build_de_summary_counts_generic(
        de_stats     = de_stats,
        contrasts    = contrasts,
        pass_col_for = function(cn) paste0("pass.imputs.", cn),
        fc_col_for   = function(cn, cols) {
            candidates <- c(
                paste0("logFC_", cn),
                paste0("log2FoldChange_", cn),
                paste0("logFC.", cn),
                paste0("linearFC.imputs.", cn)
            )
            matched <- candidates[candidates %in% cols]
            if (length(matched) > 0) matched[1] else NULL
        }
    )

    if (!is.null(out_dir) && !is.null(result) && nrow(result) > 0) {
        save_tsv(result, out_dir, "de_summary_counts.tsv")
    }

    result
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

