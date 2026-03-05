# R/modules/metabolomics/00_mod_preprocessing.R
#
# Target-ready wrappers for the new met_* preprocessing pipeline.
#
# Each function is called directly from the {targets} target bodies in
# pipe_metabolomics().  They orchestrate domain functions from:
#   R/domain/metabolomics/00_inputs.R      (load / parse)
#   R/domain/metabolomics/08_missingness.R (classify / filter)
#   R/domain/metabolomics/09_imputation_met.R (impute)
#   R/domain/metabolomics/10_drift_correction.R (LOESS)
#   R/domain/metabolomics/01_normalization.R   (norm_*, transform_metab)
#   R/core/08_qc.R                         (qc_pca_scatter, norm_boxplot)
#   R/core/06_plots.R                      (qc_missingness_heatmap via 08_qc.R)


# ==============================================================================
# mod_met_raw — load, parse, apply sample filter, align metadata
# ==============================================================================

#' Load and parse raw metabolomics inputs
#'
#' Dispatches to the appropriate format parser, applies optional sample filter
#' (QC/blank exclusion), and aligns metadata to the expression matrix.
#'
#' @param inp    Pre-loaded inputs list from \code{load_metabolomics_inputs()}.
#' @param config Full pipeline config list.
#' @return list with: \code{expr_raw}, \code{meta}, \code{row_data},
#'   \code{sample_col}, \code{format}.
#'
mod_met_raw <- function(inp, config) {
  cfg        <- config$modes$metabolomics

  # Use the pre-loaded 'inp' object provided by the previous target.
  # This avoids redundant I/O and ensures data consistency.
  fmt        <- inp$format %||% cfg$input$format %||% "cd_raw"
  sample_col <- cfg$effects$samples %||% "sample_id"
  
  # Dispatch to the appropriate parser based on the defined format.
  # Note: 'inp$data' now contains the processed/merged data from load_metabolomics_inputs.

  fmt        <- inp$format %||% cfg$input$format %||% "cd_raw"
  sample_col <- cfg$effects$samples %||% "sample_id"
  map_from   <- cfg$id_columns$map_from
  map_to     <- cfg$id_columns$map_to

  # For processed_wide the data columns ARE the sample identifiers:
  # parse_processed_wide matches them against metadata IDs, so the rename
  # must happen before the parser sees the data.
  inp_data <- inp$data
  if (!is.null(inp$sample_map) && fmt == "processed_wide") {
    if (is.null(map_from) || is.null(map_to))
      stop("mod_met_raw: id_columns$map_from and id_columns$map_to are required ",
           "when files$sample_map is set")
    inp_data <- apply_sample_map_to_colnames(inp_data, inp$sample_map, map_from, map_to)
  }

  parsed <- switch(fmt,
    cd_raw         = parse_cd_raw(inp_data, cfg),
    processed_wide = parse_processed_wide(inp_data, cfg, inp$metadata),
    multi_level    = parse_multi_level(inp_data, cfg, inp$metadata),
    stop("mod_met_raw: unsupported format: '", fmt, "'")
  )

  expr_raw <- parsed$expr_raw
  row_data <- parsed$row_data

  meta     <- inp$metadata
  
  # Step 2: Apply Sample Mapping (Fix 1: Explicit Map Columns)
  if (!is.null(inp$sample_map)) {
    map_from <- cfg$id_columns$map_from
    map_to   <- cfg$id_columns$map_to
    
    if (is.null(map_from) || is.null(map_to)) {
      stop("mod_met_raw: sample_map is provided but 'map_from' or 'map_to' are missing in config.")
    }
    
    message("mod_met_raw: applying sample map to expression matrix column names.")
    expr_raw <- apply_sample_map_to_colnames(expr_raw, inp$sample_map, map_from, map_to)
  }
  
  # Step 3: Align Metadata
  # Now colnames(expr_raw) are the mapped/final IDs
  meta <- meta %||% build_minimal_meta(colnames(expr_raw))
  meta <- align_meta_to_matrix(colnames(expr_raw), meta, sample_col)
  
  # Step 4: Apply optional sample filters (QC/Blank exclusion)
  rules <- get_sample_filter_rules_metab(cfg)
  if (!is.null(rules)) {
    keep_ids <- apply_sample_filter_metab(colnames(expr_raw), meta, rules, sample_col)
    if (length(keep_ids) < ncol(expr_raw)) {
      message(sprintf(
        "mod_met_raw: sample filter removed %d sample(s); retaining %d.",
        ncol(expr_raw) - length(keep_ids), length(keep_ids)
      ))
      expr_raw <- expr_raw[, keep_ids, drop = FALSE]
      meta     <- meta[meta[[sample_col]] %in% keep_ids, , drop = FALSE]
    }
  }

  # Normalize row_data rownames for downstream feature-wise subsetting

  if (!is.null(row_data) && !is.null(row_data$feature_id)) {
    rownames(row_data) <- row_data$feature_id
  }
  

  # Return a structured list containing all components for the downstream DAG.
  list(
    expr_raw   = expr_raw,
    meta       = meta,
    row_data   = row_data,
    sample_col = sample_col,
    format     = fmt
  )
}

