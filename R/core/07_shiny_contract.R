# ============================================================
# Shiny Payload Contract — Canonical Schema v2.0
# ============================================================
# This file defines the UNIFIED contract for all omics Shiny payloads.
# All builders (rnaseq, proteomics, metabolomics) MUST return payloads
# that conform to this contract.
#
# RULES:
# - All keys are snake_case (no uppercase)
# - Same keys across all omics (NULL if not applicable)
# - No redundant/derivable objects
# - Each key has exactly one semantic meaning
# ============================================================

# ------------------------------------------------------------
# Canonical Key Definitions
# ------------------------------------------------------------

#' Get canonical payload key definitions
#'
#' Returns a list describing each key in the canonical contract.
#' Used for documentation and validation.
#'
#' @return Named list with key metadata (type, required, description)
#' @export
get_payload_key_definitions <- function() {
    list(
        # --- Execution Info (3 keys) ---
        payload_version = list(
            type = "character",
            required = TRUE,
            description = "Schema version (e.g., '2.0')"
        ),
        payload_created_at = list(
            type = "POSIXct",
            required = TRUE,
            description = "Timestamp when payload was created"
        ),
        payload_source = list(
            type = "character",
            required = TRUE,
            description = "Omics type identifier ('rnaseq', 'proteomics', 'metabolomics')"
        ),

        # --- Metadata (3 keys) ---
        sample_meta = list(
            type = "data.frame",
            required = TRUE,
            description = "Sample metadata; rownames = sample IDs matching expr colnames"
        ),
        feature_annot = list(
            type = "data.frame",
            required = FALSE,
            description = "Feature annotations (gene symbols, descriptions, etc.)"
        ),

        # --- Expression Data (3 keys) ---
        expr_raw = list(
            type = "matrix",
            required = TRUE,
            description = "Filtered expression matrix (features x samples); may have NAs in proteomics"
        ),
        expr_norm = list(
            type = "matrix",
            required = TRUE,
            description = "Normalized expression matrix (features x samples); NO NAs allowed"
        ),
        expr_long = list(
            type = "data.frame",
            required = FALSE,
            description = "Long-format expression data (feature_id, sample_id, value + metadata)"
        ),

        # --- QC/PCA (4 keys) ---
        pca_object = list(
            type = "prcomp",
            required = FALSE,
            description = "PCA result object (for variance explained, loadings)"
        ),
        pca_scores = list(
            type = "data.frame",
            required = FALSE,
            description = "PCA scores with metadata (PC1, PC2, ..., sample info)"
        ),
        pca_3d = list(
            type = "plotly",
            required = FALSE,
            description = "3-D PCA plotly widget (PC1 x PC2 x PC3)"
        ),
        
        imp_hist_samp = list(
            type = "ggplot",
            required = FALSE,
            description = "Imputation histogram plot (ggplot object)"
        ),
        samples_hm = list(
            type = "pheatmap",
            required = FALSE,
            description = "Samples distance heatmap "
        ),
        samples_hm_w_na = list(
            type = "pheatmap",
            required = FALSE,
            description = "Samples distance heatmap with NA (relevant to proteomics)"
        ),

        # --- DE Results (5 keys) ---
        de_model = list(
            type = "model",
            required = FALSE,
            description = "Fitted DE model object (DESeq2 dds, limma fit, etc.)"
        ),
        de_stats = list(
            type = "data.frame",
            required = FALSE,
            description = "Full DE statistics for all features and contrasts"
        ),
        de_sig_stats = list(
            type = "data.frame",
            required = FALSE,
            description = "Subset of de_stats where pass_any_contrast == TRUE"
        ),
        de_expr_norm = list(
            type = "matrix",
            required = FALSE,
            description = "Expression matrix subset to significant DE features only"
        ),
        de_summary = list(
            type = "data.frame",
            required = FALSE,
            description = "Per-contrast summary counts (up, down, total)"
        ),
        de_final_table = list(
            type = "data.frame",
            required = FALSE,
            description = "DE-significant rows from final results table (equivalent to Final_results_DE_P_*.xlsx content)"
        ),

        # --- Clustering (4 keys) ---
        clust_partition = list(
            type = "data.frame",
            required = FALSE,
            description = "Partition clustering assignments (feature_id, cluster)"
        ),
        clust_patterns = list(
            type = "list",
            required = FALSE,
            description = "Binary clustering patterns"
        ),
        clust_heatmap_hier = list(
            type = "list",
            required = FALSE,
            description = "Hierarchical clustering (pheatmap result + dendrogram)"
        ),
        clust_heatmaps_by_pattern = list(
            type = "list",
            required = FALSE,
            description = "Named list of pheatmap objects per pattern"
        ),
        clust_heatmap_partition = list(
            type = "list",
            required = FALSE,
            description = "Partition clustering heatmap (pheatmap figure)"
        ),
        clust_heatmap_hier_fig = list(
            type = "list",
            required = FALSE,
            description = "Hierarchical clustering heatmap figure (pheatmap object only, extracted from clust_heatmap_hier)"
        ),

        # --- Configuration (6 keys) ---
        padj_cutoff = list(
            type = "numeric",
            required = TRUE,
            description = "Adjusted p-value threshold for significance"
        ),
        log_fc_cutoff = list(
            type = "numeric",
            required = TRUE,
            description = "Log2 fold change threshold"
        ),
        norm_method = list(
            type = "character",
            required = TRUE,
            description = "Normalization method used (e.g., 'TMM', 'VST', 'none')"
        ),
        group = list(
            type = "character",
            required = TRUE,
            description = "Primary grouping variable (column name in sample_meta)"
        ),
        color = list(
            type = "character",
            required = FALSE,
            description = "Color aesthetic variable (column name in sample_meta)"
        ),
        shape = list(
            type = "character",
            required = FALSE,
            description = "Shape aesthetic variable (column name in sample_meta)"
        )
    )
}

#' Get list of all canonical key names
#'
#' @return Character vector of canonical key names (27 keys)
#' @export
get_canonical_keys <- function() {
    names(get_payload_key_definitions())
}

#' Get required canonical keys
#'
#' @return Character vector of required key names
#' @export
get_required_keys <- function() {
    defs <- get_payload_key_definitions()
    names(defs)[vapply(defs, function(x) isTRUE(x$required), logical(1))]
}

#' Get optional canonical keys
#'
#' @return Character vector of optional key names
#' @export
get_optional_keys <- function() {
    setdiff(get_canonical_keys(), get_required_keys())
}


# ------------------------------------------------------------
# Payload Initialization
# ------------------------------------------------------------

#' Initialize empty canonical Shiny payload
#'
#' Creates a payload structure with all 27 canonical keys set to NULL.
#' Use this as the starting point for all omics builders.
#'
#' @param source Character. Omics type: "rnaseq", "proteomics", or "metabolomics"
#' @param version Character. Schema version (default: "2.0")
#'
#' @return Named list with all canonical keys initialized
#' @export
init_shiny_payload <- function(source, version = "2.0") {
    # Validate source
    valid_sources <- c("rnaseq", "proteomics", "metabolomics")
    if (!source %in% valid_sources) {
        stop("source must be one of: ", paste(valid_sources, collapse = ", "))
    }

    # Initialize all keys to NULL
    keys <- get_canonical_keys()
    payload <- setNames(vector("list", length(keys)), keys)

    # Set execution info (always populated)
    payload$payload_version <- version
    payload$payload_created_at <- Sys.time()
    payload$payload_source <- source

    # Set config defaults (required keys need non-NULL defaults)
    payload$padj_cutoff <- 0.05
    payload$log_fc_cutoff <- log2(1.5)
    payload$norm_method <- "unknown"
    payload$group <- "unknown"

    payload
}


# ------------------------------------------------------------
# Payload Validation
# ------------------------------------------------------------

#' Assert Shiny payload conforms to canonical contract
#'
#' Validates that a payload has all required keys, correct types,
#' and passes dimensional consistency checks.
#'
#' @param payload Named list to validate
#' @param strict Logical. If TRUE, fail on unexpected keys (default: FALSE)
#' @param context Character. Context string for error messages
#'
#' @return payload (invisibly) if valid; stops with error otherwise
#' @export
assert_shiny_payload_contract <- function(payload, strict = FALSE, context = "payload") {
  prefix <- sprintf("[%s]", context)
  
  signal <- function(msg) {
    if (strict) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, call. = FALSE)
    }
  }
  
  # 1. Check it's a list
  if (!is.list(payload)) {
    signal(paste0(prefix, " must be a list, got ", class(payload)[1]))
    return(invisible(payload))
  }
  
  # 2. Check all canonical keys exist
  canonical_keys <- get_canonical_keys()
  missing_keys <- setdiff(canonical_keys, names(payload))
  if (length(missing_keys) > 0) {
    signal(paste0(
      prefix, " missing canonical keys: ",
      paste(missing_keys, collapse = ", ")
    ))
  }
  
  # 3. Check required keys are not NULL
  required_keys <- get_required_keys()
  null_required <- required_keys[vapply(payload[required_keys], is.null, logical(1))]
  if (length(null_required) > 0) {
    signal(paste0(
      prefix, " required keys are NULL: ",
      paste(null_required, collapse = ", ")
    ))
  }
  
  # 4. Type checks for critical keys
  
  # payload_source
  if (!is.character(payload$payload_source) || length(payload$payload_source) != 1) {
    signal(paste0(prefix, " payload_source must be a single string"))
  }
  
  # sample_meta
  if (!is.data.frame(payload$sample_meta)) {
    signal(paste0(prefix, " sample_meta must be a data.frame"))
  }
  
  if (!is.null(payload$sample_meta)) {
    if (is.null(rownames(payload$sample_meta)) ||
        any(rownames(payload$sample_meta) == "")) {
      signal(paste0(prefix, " sample_meta must have non-empty rownames (sample IDs)"))
    }
  }
  
  # expr_raw
  if (!is.matrix(payload$expr_raw)) {
    signal(paste0(prefix, " expr_raw must be a matrix"))
  }
  if (!is.numeric(payload$expr_raw)) {
    signal(paste0(prefix, " expr_raw must be numeric"))
  }
  
  # expr_norm
  if (!is.matrix(payload$expr_norm)) {
    signal(paste0(prefix, " expr_norm must be a matrix"))
  }
  if (!is.numeric(payload$expr_norm)) {
    signal(paste0(prefix, " expr_norm must be numeric"))
  }
  if (anyNA(payload$expr_norm)) {
    signal(paste0(prefix, " expr_norm must NOT contain NA values (imputation required)"))
  }
  
  # 5. Dimensional consistency
  expr_samples <- colnames(payload$expr_norm)
  meta_samples <- rownames(payload$sample_meta)
  
  if (is.null(expr_samples)) {
    signal(paste0(prefix, " expr_norm must have colnames (sample IDs)"))
  }
  
  if (!setequal(expr_samples, meta_samples)) {
    missing_in_meta <- setdiff(expr_samples, meta_samples)
    missing_in_expr <- setdiff(meta_samples, expr_samples)
    
    msg <- paste0(prefix, " sample mismatch between expr_norm and sample_meta")
    
    if (length(missing_in_meta) > 0) {
      msg <- paste0(
        msg,
        "\n  In expr but not meta: ",
        paste(head(missing_in_meta, 5), collapse = ", ")
      )
    }
    
    if (length(missing_in_expr) > 0) {
      msg <- paste0(
        msg,
        "\n  In meta but not expr: ",
        paste(head(missing_in_expr, 5), collapse = ", ")
      )
    }
    
    signal(msg)
  }
  
  # 6. Check de_expr_norm subset relationship (if present)
  if (!is.null(payload$de_expr_norm)) {
    if (!is.matrix(payload$de_expr_norm)) {
      signal(paste0(prefix, " de_expr_norm must be a matrix"))
    }
    
    de_features <- rownames(payload$de_expr_norm)
    all_features <- rownames(payload$expr_norm)
    
    if (!all(de_features %in% all_features)) {
      signal(paste0(prefix,
                    " de_expr_norm features must be a subset of expr_norm features"))
    }
  }
  
  # 7. Config variable validation
  if (!is.null(payload$color)) {
    if (!all(payload$color %in% colnames(payload$sample_meta))) {
      signal(paste0(prefix, " color '", payload$color,
                    "' not found in sample_meta columns"))
    }
  }
  
  if (!is.null(payload$shape)) {
    if (!all(payload$shape %in% colnames(payload$sample_meta))) {
      signal(paste0(prefix, " shape '", payload$shape,
                    "' not found in sample_meta columns"))
    }
  }
  
  invisible(payload)
}

#' Validate DE stats data.frame structure
#'
#' Checks that de_stats has expected columns and types.
#'
#' @param de_stats Data.frame to validate
#' @param context Character. Context for error messages
#'
#' @return de_stats (invisibly) if valid
#' @export
assert_de_stats <- function(de_stats, context = "de_stats") {
    prefix <- sprintf("[%s]", context)

    if (is.null(de_stats)) {
        return(invisible(NULL))
    }

    if (!is.data.frame(de_stats)) {
        stop(prefix, " must be a data.frame, got ", class(de_stats)[1], call. = FALSE)
    }

    # Required columns
    required_cols <- c("feature_id", "pass_any_contrast")
    missing <- setdiff(required_cols, colnames(de_stats))
    if (length(missing) > 0) {
        stop(prefix, " missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    # Type checks
    if (!is.character(de_stats$feature_id) && !is.factor(de_stats$feature_id)) {
        stop(prefix, " feature_id must be character or factor", call. = FALSE)
    }

    if (!is.numeric(de_stats$pass_any_contrast) && !is.logical(de_stats$pass_any_contrast)) {
        stop(prefix, " pass_any_contrast must be numeric (0/1) or logical", call. = FALSE)
    }

    invisible(de_stats)
}


# ------------------------------------------------------------
# Legacy Alias Helpers
# ------------------------------------------------------------

#' Attach legacy aliases to a canonical payload (pass-through)
#'
#' @param payload Canonical payload list
#' @return payload unchanged (legacy aliases are added by domain builders)
#' @export
attach_legacy_aliases <- function(payload) {
    payload
}

#' Remove legacy aliases from payload
#'
#' Strips legacy keys, leaving only canonical keys.
#' Use for clean exports or testing.
#'
#' @param payload Payload with legacy aliases
#'
#' @return Payload with only canonical keys
#' @export
# strip_legacy_aliases <- function(payload) {
#     canonical_keys <- get_canonical_keys()
#     payload[canonical_keys]
# }

#' Attach legacy aliases to a Shiny payload
#'
#' Placeholder — returns payload unchanged. Legacy alias mappings
#' can be added here as needed for backward compatibility.
#'
#' @param payload Shiny payload list
#' @return payload (unchanged)
attach_legacy_aliases <- function(payload) {
    payload
}


# ------------------------------------------------------------
# Legacy Base Builder (for backward compatibility during transition)
# ------------------------------------------------------------

#' Build base legacy Shiny structure (generic)
#'
#' DEPRECATED: Use init_shiny_payload() for new code.
#' This function exists only to support existing domain builders during transition.
#'
#' CRITICAL: This function preserves EXACT legacy values.
#' Do NOT modify legacy_source strings - they must match pre-refactor values exactly.
#'
#' @param legacy_source EXACT legacy source string (e.g., "Proteomics pipeline", "RNAseq pipeline")
#' @param pca_basename PCA 3D basename (e.g., "prot_pca_3d", "rna_pca_3d")
#' @return List with all legacy keys initialized to NULL
#' @export
# build_shiny_legacy_base <- function(legacy_source, pca_basename) {
    # Validate inputs
    # if (!is.character(legacy_source) || length(legacy_source) != 1 || nchar(legacy_source) == 0) {
    #     stop("legacy_source must be a non-empty string")
    # }
    # if (!is.character(pca_basename) || length(pca_basename) != 1 || nchar(pca_basename) == 0) {
    #     stop("pca_basename must be a non-empty string")
    # }
    # 
    # # Return structure with EXACT legacy values
    # list(
    #     # Metadata
    #     legacy_version = "1.0",
    #     legacy_created_at = Sys.time(),
    #     legacy_source = legacy_source, # EXACT value, no transformations
    # 
    #     # Data objects (full union of keys from both modes)
    #     col_data = NULL,
    #     contrasts_data = NULL,
        # dds = NULL,
        # norm_counts = NULL,
        # norm_log_counts = NULL,
        # norm_log_counts_pca = NULL, # PCA prcomp object (optional)
        # mat2plot = NULL, # LEGACY: DE expression matrix (features x samples) - use de_expr_norm in new code
        # de_expr_norm = NULL, # DE expression matrix (features x samples) - canonical key
        # EFFECTS = NULL,
        # PCA_3D_BASENAME = pca_basename,
#         stats_df = NULL,
#         DE_genes_stats = NULL,
#         hm_hier_de = NULL,
#         patterns = NULL,
#         heatmaps_by_pattern = NULL,
#         New_clusters = NULL,
#         annot = NULL,
#         trinotate_main = NULL,
# 
#         # Config parameters
#         PADJ_CUTOFF = NULL,
#         DESEQ_PADJ_CUTOFF = NULL,
#         LOG_FC_CUTOFF = NULL,
#         LINEAR_FC_CUTOFF = NULL,
#         NORM_METHOD = NULL,
#         GROUP = NULL,
# 
#         # Shiny aesthetics (column names)
#         color = NULL,
#         shape = NULL
#     )
# }


# ------------------------------------------------------------
# Expression Data Helpers
# ------------------------------------------------------------

#' Build long-format expression data.frame
#'
#' Converts a wide expression matrix to long format and joins with sample metadata.
#' Useful for ggplot2-based visualizations.
#'
#' @param expr_mat Numeric matrix (features x samples) with rownames and colnames
#' @param sample_meta Data.frame with sample metadata; rownames must match expr_mat colnames
#'
#' @return Data.frame with columns: feature_id, sample_id, value, + all sample_meta columns
#' @export
build_expr_long <- function(expr_mat, sample_meta) {
    if (is.null(expr_mat) || is.null(sample_meta)) {
        return(NULL)
    }

    # Convert matrix to long format
    expr_df <- as.data.frame(expr_mat)
    expr_df$feature_id <- rownames(expr_mat)

    expr_long <- tidyr::pivot_longer(
        expr_df,
        cols = -feature_id,
        names_to = "sample_id",
        values_to = "value"
    )

    # Join with sample metadata
    meta_df <- sample_meta
    meta_df$sample_id <- rownames(sample_meta)

    expr_long <- merge(expr_long, meta_df, by = "sample_id", all.x = TRUE)

    expr_long
}
