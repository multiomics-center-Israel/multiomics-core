# ============================================================
# Shiny Legacy Export — RNA-seq
# ============================================================
# This file exists ONLY to support legacy Shiny applications.
# Do NOT use objects from here inside the pipeline.
# This is a CONSUMER layer, not a PRODUCER.
# ============================================================

#' Build legacy-compatible data structure for Shiny app
#'
#' This function collects pre-computed objects from the RNA-seq pipeline
#' and assembles them into a flat named list with exact legacy variable names.
#'
#' CRITICAL: This function is a CONSUMER ONLY. It does NOT compute PCA,
#' clustering, or any transformations. It only collects and renames.
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
    # ============================================================
    # Initialize with ALL keys set to NULL (predictable structure)
    # Using generic base builder with EXACT legacy values
    # ============================================================
    legacy <- build_shiny_legacy_base(
        legacy_source = "RNAseq pipeline", # EXACT pre-refactor value (not "RNA-seq")
        pca_basename = "rna_pca_3d"
    )

    # ============================================================
    # Populate critical objects
    # ============================================================

    # Metadata
    legacy$col_data <- pre$meta
    legacy$contrasts_data <- inputs$contrasts

    # Expression matrices
    legacy$norm_counts <- pre$expr_filt
    legacy$norm_log_counts <- pre$expr_work

    # DESeq2 object (conservative: use stored or NULL)
    legacy$dds <- de_res$dds %||% NULL
    if (is.null(legacy$dds)) {
        warning("dds not found in de_res; legacy Shiny DE recomputation will be disabled")
    }

    # ============================================================
    # PCA objects (NO COMPUTATION - only collect)
    # ============================================================
    if (!is.null(pca_res)) {
        legacy$norm_log_counts_pca <- pca_res$norm_log_counts_pca %||% NULL
        # NOTE: In legacy Proteomics, mat2plot IS PCA scores.
        # But in legacy RNA-seq Neat script, mat2plot IS heatmap expression matrix.
        # We will populate mat2plot logic below in the Heatmap section to match RNA legacy.
    } else {
        # Use message instead of warning to reduce noise when PCA is intentionally skipped
        message("pca_res is NULL; PCA plots will be disabled in Shiny")
    }

    # ============================================================
    # EFFECTS (critical for Shiny PCA)
    # ============================================================
    # Safe navigation: step-by-step to avoid errors on missing intermediate paths
    modes <- config$modes %||% list()
    rna <- modes$rna %||% list()
    effects <- rna$effects %||% list()

    if (!is.null(effects$shape) && !is.null(effects$color)) {
        legacy$EFFECTS <- c(effects$shape, effects$color)
    } else {
        warning("EFFECTS not fully defined in config; setting to NULL")
        legacy$EFFECTS <- NULL
    }

    # ============================================================
    # DE summary and stats
    # ============================================================
    # Extract dds object (required by Shiny)
    legacy$dds <- de_res$dds

    legacy$stats_df <- build_rnaseq_summary_df(de_res$tables, config)
    # Ensure ID column is named "Gene" (standard for RNA legacy Shiny)
    if (!is.null(legacy$stats_df) && "FeatureID" %in% names(legacy$stats_df) && !"Gene" %in% names(legacy$stats_df)) {
        names(legacy$stats_df)[names(legacy$stats_df) == "FeatureID"] <- "Gene"
    }

    # Filter for significant genes
    if (!is.null(legacy$stats_df) && "pass_any_contrast" %in% colnames(legacy$stats_df)) {
        legacy$DE_genes_stats <- legacy$stats_df[
            !is.na(legacy$stats_df$pass_any_contrast) & legacy$stats_df$pass_any_contrast == 1, ,
            drop = FALSE
        ]
    }

    # ============================================================
    # Heatmap data (mat2plot)
    # ============================================================
    # In legacy RNA script: mat2plot = filter_expression_matrix_by_gene_list(norm_log_counts, top_DE_genes)
    if (!is.null(legacy$norm_log_counts) && !is.null(legacy$DE_genes_stats)) {
        sig_genes <- legacy$DE_genes_stats$Gene
        if (length(sig_genes) > 0) {
            # Filter expression matrix
            legacy$mat2plot <- legacy$norm_log_counts[rownames(legacy$norm_log_counts) %in% sig_genes, , drop = FALSE]
        } else {
            legacy$mat2plot <- NULL
        }
    } else {
        legacy$mat2plot <- NULL
    }

    # Robust: check multiple possible field names for pre-computed heatmaps
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
        # Check if objects are nested in 'objects' list (standard for modules)
        # or at top level (legacy fallback)
        src <- if (!is.null(clustering_res$objects)) clustering_res$objects else clustering_res

        if (!is.null(src$patterns)) legacy$patterns <- src$patterns
        if (!is.null(src$heatmaps)) legacy$heatmaps_by_pattern <- src$heatmaps
        if (!is.null(src$clusters)) legacy$New_clusters <- src$clusters
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
    rna <- modes$rna %||% list()
    de_cfg <- rna$de %||% list()
    norm_cfg <- rna$normalization %||% list()

    legacy$PADJ_CUTOFF <- de_cfg$padj_cutoff %||% de_cfg$p_cutoff %||% 0.05
    legacy$DESEQ_PADJ_CUTOFF <- legacy$PADJ_CUTOFF

    linear_fc <- de_cfg$linear_fc_cutoff %||% 1.5
    legacy$LINEAR_FC_CUTOFF <- linear_fc
    legacy$LOG_FC_CUTOFF <- log2(linear_fc)

    legacy$NORM_METHOD <- norm_cfg$method %||% "TMMlogCPM"
    legacy$GROUP <- effects$color %||% NULL

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
                    "Dimension mismatch: norm_log_counts colnames not all found in col_data rownames\n",
                    "  Missing: ", paste(setdiff(expr_samples, meta_samples), collapse = ", ")
                )
            }

            # Check order (informational only)
            if (all(expr_samples %in% meta_samples)) {
                # Check if metadata samples in expression order match expression samples
                expected <- meta_samples[match(expr_samples, meta_samples)]
                if (!identical(expr_samples, expected)) {
                    message("Note: col_data rownames are not in the same order as expression colnames (Shiny may reorder internally)")
                }
            }
        }
    }

    # ============================================================
    # Final validation: ensure all expected keys exist
    # ============================================================
    expected_keys <- c(
        "legacy_version", "legacy_created_at", "legacy_source",
        "col_data", "contrasts_data", "dds", "norm_counts", "norm_log_counts",
        "norm_log_counts_pca", "mat2plot", "EFFECTS", "PCA_3D_BASENAME",
        "stats_df", "DE_genes_stats",
        "patterns", "heatmaps_by_pattern", "New_clusters",
        "trinotate_main",
        "PADJ_CUTOFF", "DESEQ_PADJ_CUTOFF", "LOG_FC_CUTOFF",
        "LINEAR_FC_CUTOFF", "NORM_METHOD", "GROUP"
    )

    # Note: pheatmap_data_DE_genes and annot are optional and added separately
    # They're not in expected_keys to avoid validation errors

    missing_keys <- setdiff(expected_keys, names(legacy))
    if (length(missing_keys) > 0) {
        stop("BUG: Missing expected keys in legacy export: ", paste(missing_keys, collapse = ", "))
    }

    message("Legacy export created successfully with ", length(legacy), " keys")
    message("  - ", sum(sapply(legacy, is.null)), " keys are NULL")
    message("  - ", sum(!sapply(legacy, is.null)), " keys have data")

    legacy
}


#' Save legacy data structure to RDS file
#'
#' Wrapper function that builds and saves the legacy export
#'
#' @param ... Arguments passed to build_data_to_shiny_legacy_rna()
#' @param out_file Output file path (default: "data_to_shiny_legacy.rds")
#'
#' @return Path to saved file (invisibly)
#'
#' @export
save_data_to_shiny_legacy_rna <- function(..., out_file = "data_to_shiny_legacy.rds") {
    legacy <- build_data_to_shiny_legacy_rna(...)

    # Ensure output directory exists
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
#' Helper function to load legacy RDS and check structure
#' Useful for debugging and testing
#'
#' @param file Path to legacy RDS file
#' @param verbose Print diagnostic information
#'
#' @return The loaded legacy list
#'
#' @export
load_legacy_rds <- function(file, verbose = TRUE) {
    if (!file.exists(file)) {
        stop("File not found: ", file)
    }

    legacy <- readRDS(file)

    if (verbose) {
        message("Loaded legacy export from: ", file)
        message("  Version: ", legacy$legacy_version %||% "unknown")
        message("  Created: ", legacy$legacy_created_at %||% "unknown")
        message("  Source: ", legacy$legacy_source %||% "unknown")
        message("  Total keys: ", length(legacy))
        message("  NULL keys: ", sum(sapply(legacy, is.null)))
        message("  Non-NULL keys: ", sum(!sapply(legacy, is.null)))

        # Check critical objects
        critical <- c("col_data", "norm_log_counts", "stats_df", "EFFECTS")
        for (key in critical) {
            status <- if (is.null(legacy[[key]])) "MISSING" else "OK"
            message("  [", status, "] ", key)
        }
    }

    legacy
}
