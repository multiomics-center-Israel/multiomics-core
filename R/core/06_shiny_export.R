#' Build base legacy Shiny structure (generic)
#'
#' CRITICAL: This function preserves EXACT legacy values.
#' Do NOT modify legacy_source strings - they must match pre-refactor values exactly.
#'
#' @param legacy_source EXACT legacy source string (e.g., "Proteomics pipeline", "RNAseq pipeline")
#' @param pca_basename PCA 3D basename (e.g., "prot_pca_3d", "rna_pca_3d")
#' @return List with all legacy keys initialized to NULL
build_shiny_legacy_base <- function(legacy_source, pca_basename) {
    # Validate inputs
    if (!is.character(legacy_source) || length(legacy_source) != 1 || nchar(legacy_source) == 0) {
        stop("legacy_source must be a non-empty string")
    }
    if (!is.character(pca_basename) || length(pca_basename) != 1 || nchar(pca_basename) == 0) {
        stop("pca_basename must be a non-empty string")
    }

    # Return structure with EXACT legacy values
    list(
        # Metadata
        legacy_version = "1.0",
        legacy_created_at = Sys.time(),
        legacy_source = legacy_source, # EXACT value, no transformations

        # Data objects (full union of keys from both modes)
        col_data = NULL,
        contrasts_data = NULL,
        dds = NULL,
        norm_counts = NULL,
        norm_log_counts = NULL,
        norm_log_counts_pca = NULL, # PCA prcomp object (optional)
        mat2plot = NULL, # PCA scores data.frame (optional)
        EFFECTS = NULL,
        PCA_3D_BASENAME = pca_basename,
        stats_df = NULL,
        DE_genes_stats = NULL,
        pheatmap_data_DE_genes = NULL,
        patterns = NULL,
        heatmaps_by_pattern = NULL,
        New_clusters = NULL,
        annot = NULL,
        trinotate_main = NULL,

        # Config parameters
        PADJ_CUTOFF = NULL,
        DESEQ_PADJ_CUTOFF = NULL,
        LOG_FC_CUTOFF = NULL,
        LINEAR_FC_CUTOFF = NULL,
        NORM_METHOD = NULL,
        GROUP = NULL,

        # Shiny aesthetics (column names)
        color = NULL,
        shape = NULL
    )
}
