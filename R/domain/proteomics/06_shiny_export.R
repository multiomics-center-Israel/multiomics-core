# ============================================================
# Shiny Legacy Export — Proteomics
# ============================================================
# This file exists ONLY to support legacy Shiny applications.
# Do NOT use objects from here inside the pipeline.
# This is a CONSUMER layer, not a PRODUCER.
# ============================================================

#' Build legacy-compatible data structure for Shiny app (Proteomics)
#'
#' This function collects pre-computed objects from the proteomics pipeline
#' and assembles them into a flat named list with exact legacy variable names.
#'
#' CRITICAL: This function is a CONSUMER ONLY. It does NOT compute PCA,
#' clustering, or any transformations. It only collects and renames.
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
    # ============================================================
    # Initialize with ALL keys set to NULL (predictable structure)
    # Using generic base builder with EXACT legacy values
    # ============================================================
    legacy <- build_shiny_legacy_base(
        legacy_source = "Proteomics pipeline", # EXACT pre-refactor value
        pca_basename = "prot_pca_3d"
    )

    # Targets-safe: base object may be old
    if (!("color" %in% names(legacy))) legacy["color"] <- list(NULL)
    if (!("shape" %in% names(legacy))) legacy["shape"] <- list(NULL)

    # ============================================================
    # Populate critical objects
    # ============================================================

    # Metadata
    legacy$col_data <- pre$meta %||% NULL

    # Patch: Set rownames from samples column if specified
    prot_fx <- config$modes$proteomics[["effects"]] %||% list()
    sample_col <- prot_fx$samples %||% NULL
    if (!is.null(legacy$col_data) && !is.null(sample_col) && sample_col %in% colnames(legacy$col_data)) {
        rownames(legacy$col_data) <- as.character(legacy$col_data[[sample_col]])
    }

    legacy$contrasts_data <- inputs$contrasts %||% NULL

    # Expression matrices
    legacy$norm_counts <- pre$expr_filt %||% NULL
    legacy$norm_log_counts <- pre$expr_work %||% NULL

    # ============================================================
    # PCA objects + color/shape aesthetics (from QC pre module)
    # ============================================================
    if (!is.null(pca_res) && !is.null(pca_res$objects)) {
        legacy$norm_log_counts_pca <- pca_res$objects$norm_log_counts_pca %||% NULL
        # legacy$mat2plot was incorrectly assigning PCA scores here.
        # We will calculate it below as the expression matrix of significant features.
        legacy$color <- pca_res$objects$color %||% NULL
        legacy$shape <- pca_res$objects$shape %||% NULL
    } else {
        message("pca_res$objects is NULL; PCA plots will be disabled in Shiny (proteomics)")
    }

    # Fallback: get color/shape from config if not provided by QC
    prot_fx <- config$modes$proteomics[["effects"]] %||% list()
    val_color <- legacy[["color"]] %||% prot_fx[["color"]]
    val_shape <- legacy[["shape"]] %||% prot_fx[["shape"]]

    legacy["color"] <- list(val_color)
    legacy["shape"] <- list(val_shape)

    # Ensure keys always exist (before validation)
    for (k in c("color", "shape")) {
        if (!k %in% names(legacy)) legacy[k] <- list(NULL)
    }

    # ============================================================
    # DE summary and stats
    # ============================================================
    legacy$stats_df <- de_res$summary_df %||% NULL

    # Initialize DE_genes_stats (using list syntax to guarantee key exists)
    legacy["DE_genes_stats"] <- list(NULL)
    if (!is.null(legacy$stats_df) && "pass_any_contrast" %in% colnames(legacy$stats_df)) {
        legacy$DE_genes_stats <- legacy$stats_df[
            !is.na(legacy$stats_df$pass_any_contrast) & legacy$stats_df$pass_any_contrast == 1, ,
            drop = FALSE
        ]
    }

    # ============================================================
    # Heatmap data matrix (mat2plot)
    # ============================================================
    # Initialize mat2plot (using list syntax to guarantee key exists even if NULL)
    legacy["mat2plot"] <- list(NULL)
    # Filter expression matrix by significant DE features (matches RNAseq logic)
    if (!is.null(legacy$norm_log_counts) && !is.null(legacy$DE_genes_stats) && nrow(legacy$DE_genes_stats) > 0) {
        # Use config-defined ID column (fallback to common column names)
        prot_cfg <- config$modes$proteomics %||% list()
        de_table_cfg <- prot_cfg$de_table %||% list()
        id_col <- de_table_cfg$id_col %||% "FeatureID"

        # Try configured column, then fallback to common alternatives
        possible_cols <- c(id_col, "Protein", "FeatureID", "Protein.Group")
        id_col_found <- intersect(possible_cols, colnames(legacy$DE_genes_stats))[1]

        if (!is.na(id_col_found)) {
            sig_ids <- legacy$DE_genes_stats[[id_col_found]]
            legacy$mat2plot <- legacy$norm_log_counts[rownames(legacy$norm_log_counts) %in% sig_ids, , drop = FALSE]
        }
    }

    # ============================================================
    # DE Model (limma fit or DESeq2 dds)
    # ============================================================
    legacy$de_model <- de_res$de_model %||% NULL

    # ============================================================
    # Heatmap data (NO EXTRACTION - only collect)
    # ============================================================
    # Robust: check multiple possible field names
    de_heat <- de_res$pheatmap_data_DE_genes %||% de_res$pheatmap_data %||% NULL
    clust_heat <- if (!is.null(clustering_res)) {
        clustering_res$pheatmap_data_DE_genes %||% clustering_res$pheatmap_data %||% NULL
    } else {
        NULL
    }
    legacy$pheatmap_data_DE_genes <- de_heat %||% clust_heat %||% NULL

    # ============================================================
    # Clustering results (NO COMPUTATION - only collect)
    # ============================================================
    if (!is.null(clustering_res)) {
        legacy$patterns <- clustering_res$patterns %||% NULL
        legacy$heatmaps_by_pattern <- clustering_res$heatmaps %||% NULL
        legacy$New_clusters <- clustering_res$clusters %||% NULL
    }

    # ============================================================
    # Annotation
    # ============================================================
    legacy$annot <- annot

    # ============================================================
    # Config parameters
    # ============================================================
    # Safe navigation: step-by-step to avoid errors on missing intermediate paths
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

    # Get effects safely (avoid conflict with stats::effects function)
    prot_effects_cfg <- prot[["effects"]] %||% list()
    legacy$GROUP <- prot_effects_cfg$color %||% NULL

    # ============================================================
    # Dimension validation (sanity check)
    # ============================================================
    if (!is.null(legacy$norm_log_counts) && !is.null(legacy$col_data)) {
        expr_samples <- colnames(legacy$norm_log_counts)
        meta_samples <- rownames(legacy$col_data)

        if (!is.null(expr_samples) && !is.null(meta_samples)) {
            # Check all expression samples are in metadata
            if (!all(expr_samples %in% meta_samples)) {
                warning(
                    "Dimension mismatch: proteomics norm_log_counts colnames not all found in col_data rownames\n",
                    "  Missing: ", paste(setdiff(expr_samples, meta_samples), collapse = ", ")
                )
            }

            # Check order (informational only)
            if (all(expr_samples %in% meta_samples)) {
                expected <- meta_samples[match(expr_samples, meta_samples)]
                if (!identical(expr_samples, expected)) {
                    message("Note: proteomics col_data rownames not in same order as expression colnames (Shiny may reorder internally)")
                }
            }
        }
    }

    # ============================================================
    # Final validation: ensure all expected keys exist
    # ============================================================
    # Required keys (must always be present, even if NULL)
    expected_keys <- c(
        "legacy_version", "legacy_created_at", "legacy_source",
        "col_data", "contrasts_data", "dds", "norm_counts", "norm_log_counts",
        "de_model", "mat2plot",
        "PCA_3D_BASENAME", # EFFECTS removed
        "stats_df", "DE_genes_stats",
        "trinotate_main",
        "PADJ_CUTOFF", "DESEQ_PADJ_CUTOFF", "LOG_FC_CUTOFF",
        "LINEAR_FC_CUTOFF", "NORM_METHOD", "GROUP",
        "color", "shape" # Shiny aesthetics (column names)
    )

    # Defensive: ensure all expected keys exist (even if NULL)
    for (k in expected_keys) {
        if (!k %in% names(legacy)) {
            legacy[k] <- list(NULL)
        }
    }

    missing_keys <- setdiff(expected_keys, names(legacy))
    if (length(missing_keys) > 0) {
        stop("BUG: Missing expected keys in proteomics legacy export: ", paste(missing_keys, collapse = ", "))
    }

    message("Proteomics legacy export created successfully with ", length(legacy), " keys")
    message("  - ", sum(sapply(legacy, is.null)), " keys are NULL")
    message("  - ", sum(!sapply(legacy, is.null)), " keys have data")

    legacy
}


#' Save proteomics legacy data structure to RDS file
#'
#' Wrapper function that builds and saves the legacy export
#'
#' @param ... Arguments passed to build_data_to_shiny_legacy_proteomics()
#' @param out_file Output file path (default: "data_to_shiny_legacy_proteomics.rds")
#'
#' @return Path to saved file (invisibly)
#'
#' @export
save_data_to_shiny_legacy_proteomics <- function(..., out_file = "data_to_shiny_legacy_proteomics.rds") {
    legacy <- build_data_to_shiny_legacy_proteomics(...)

    # Ensure output directory exists
    out_dir <- dirname(out_file)
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }

    saveRDS(legacy, out_file)
    message("Saved proteomics legacy export to: ", out_file)

    invisible(out_file)
}
