#' Extract log2-transformed DIA-NN measurements per sample
#'
#' This function:
#' - selects sample columns from the wide protein matrix (excluding annotation columns)
#' - coerces values to numeric (DIA-NN often stores blanks as character)
#' - applies log2 transformation
#' - renames columns from raw sample names to unified SampleID (using sample_map)
#' - reorders columns to match metadata sample order (effects$samples)
#'
#' @param protein    Protein table (wide DIA-NN output)
#' @param sample_map Sample mapping table (e.g. from config$modes$proteomics$files$sample_map)
#' @param meta       Metadata table (e.g. meta_prot.csv)
#' @param cfg        config$modes$proteomics
#'
#' @return A numeric matrix/data.frame of log2-intensities,
#'         columns ordered to match meta[[effects$samples]]
get_measurements_per_sample_diann <- function(protein, sample_map, meta, cfg) {
  id_cols  <- cfg$id_columns
  eff_cols <- cfg$effects
  
  # 1) Identify non-sample columns in the protein table 
  
  annot_cols <- id_cols$protein_annot
  if (is.null(annot_cols)) {
    annot_cols <- character(0)
  } else {
    annot_cols <- unlist(annot_cols)
  }
  
  non_sample_cols <- unique(c(id_cols$protein_id, annot_cols))
  
  # Sample columns = everything else
  sample_cols <- setdiff(colnames(protein), non_sample_cols)
  if (length(sample_cols) == 0) {
    stop("No sample columns detected in protein table after excluding ID/annotation columns.")
  }
  
  df_m <- protein[, sample_cols, drop = FALSE]
  
  # 2) Coerce to numeric (DIA-NN often has blanks) 
  
  df_m <- as.data.frame(
    suppressWarnings(
      sapply(df_m, function(x) as.numeric(as.character(x)))
    ),
    row.names   = row.names(df_m),
    check.names = FALSE
  )
  
  # 3) Log2 transformation 
  
  df_m <- as.data.frame(
    log2(df_m),
    row.names   = row.names(df_m),
    check.names = FALSE
  )
  
  # 4) Rename columns: raw sample names -> SampleID using sample_map 
  
  raw_names <- colnames(df_m)
  
  map_from <- id_cols$map_from
  map_to   <- id_cols$map_to
  
  # Sanity check: mapping columns exist
  check_has_cols(sample_map, c(map_from, map_to), df_name = "sample_map")
  
  new_names <- sample_map[[map_to]][ match(raw_names, sample_map[[map_from]]) ]
  
  # If some didn't match, keep originals but warn
  unmatched <- is.na(new_names)
  if (any(unmatched)) {
    warning(
      "These DIA-NN columns did not match any row in sample_map$", map_from, ": ",
      paste(raw_names[unmatched], collapse = ", ")
    )
  }
  
  colnames(df_m) <- ifelse(unmatched, raw_names, new_names)
  
  # 5) Reorder columns to match metadata sample order
  
  meta_sample_col <- eff_cols$samples
  check_has_cols(meta, meta_sample_col, df_name = "metadata")
  
  meta_sample_ids <- meta[[meta_sample_col]]
  
  missing_in_df <- setdiff(meta_sample_ids, colnames(df_m))
  if (length(missing_in_df) > 0) {
    warning(
      "These samples from metadata were not found in DIA-NN columns after renaming: ",
      paste(missing_in_df, collapse = ", ")
    )
  }
  
  ordered_cols <- intersect(meta_sample_ids, colnames(df_m))
  df_m_ordered <- df_m[, ordered_cols, drop = FALSE]
  
  df_m_ordered
}

#' Convenience wrapper: extract DIA-NN measurements from inputs list
#'
#' @param inputs List returned by load_proteomics_inputs()
#' @param cfg    config$modes$proteomics
#'
#' @return Log2-intensity matrix/data.frame with samples in columns
get_measurements_per_sample_diann_from_inputs <- function(inputs, cfg) {
  get_measurements_per_sample_diann(
    protein    = inputs$protein,
    sample_map = inputs$sample_map,
    meta       = inputs$metadata,
    cfg        = cfg
  )
}

#' Filter proteomics data using min_count per condition
#'
#' Uses the global pass_filter() function (in preprocess_filtering.R).
#' A protein passes the filter if, in at least one condition, the number of
#' non-NA values is >= the configured threshold (default or per-condition).
#'
#' @param expr_mat    Numeric matrix (rows = proteins, cols = samples)
#' @param row_data Annotation table for proteins (same row order as expr_mat)
#' @param meta     Metadata table
#' @param cfg      config$modes$proteomics
#' @param group_col Optional override of condition column name
#'
#' @return A list:
#'   - expr_mat: filtered intensity matrix
#'   - row_data: filtered annotation table
#'   - ids: vector of protein IDs that passed
#'   - keep: logical mask of passing proteins
filter_proteomics_by_min_count <- function(
    expr_mat,
    row_data,
    meta,
    cfg,
    group_col = NULL
) {
  eff <- cfg$effects
  
  # Identify sample + condition columns
  sample_id_col <- eff$samples     # e.g. "SampleID"
  if (is.null(group_col)) {
    group_col <- eff$color       # e.g. "treatment"
  }
  
  # Check metadata contains required columns
  check_has_cols(meta, sample_id_col, df_name = "metadata")
  check_has_cols(meta, group_col,      df_name = "metadata")
  
  # Reorder expr_mat-columns to match metadata sample order
  sample_ids <- meta[[sample_id_col]]
  expr_mat <- expr_mat[, sample_ids, drop = FALSE]
  
  # Group vector
  group <- meta[[group_col]]
  
  # Read min_count from config (default or per-condition)
  min_cfg       <- cfg$filtering$min_count
  min_per_group <- extract_min_count(min_cfg, group)
  
  #  USE pass_filter() FROM THE GLOBAL FILE 
  keep <- pass_filter(
    expr_mat         = expr_mat,
    group         = group,
    min_per_group = min_per_group
  )
  
  # Filter expr_mat + annotations
  expr_mat_filt    <- expr_mat[keep, , drop = FALSE]
  row_data_filt <- row_data[keep, , drop = FALSE]
  
  # Extract IDs
  if (!is.null(rownames(expr_mat))) {
    ids <- rownames(expr_mat)[keep]
  } else if (!is.null(cfg$id_columns$protein_id) &&
             cfg$id_columns$protein_id %in% colnames(row_data)) {
    ids <- row_data[[cfg$id_columns$protein_id]][keep]
  } else {
    ids <- NULL
  }
  
  message(
    "Filtering step: kept ",
    sum(keep), " of ", length(keep), " proteins."
  )
  
  list(
    expr_mat    = expr_mat_filt,
    row_data = row_data_filt,
    ids      = ids,
    keep     = keep
  )
}


#' Impute proteomics expr_mat using Perseus-style method
#'
#' @param expr_mat Numeric matrix (proteins x samples, log2).
#' @param cfg   config$modes$proteomics (uses cfg$imputation$width, cfg$imputation$downshift).
#' @param return_flags logical, if TRUE also returns the imputed_flag matrix.
#'
#' @return If return_flags = FALSE: imputed expr_mat matrix.
#'         If TRUE: list(imputed = matrix, imputed_flag = logical matrix).
impute_proteomics_perseus <- function(expr_mat, cfg, return_flags = FALSE) {
  width     <- cfg$imputation$width     %||% 0.3
  downshift <- cfg$imputation$downshift %||% 1.8
  
  res <- perseus_impute_with_flags(
    expr_mat    = expr_mat,
    width    = width,
    downshift = downshift
  )
  
  if (return_flags) res else res$imputed
}

#' Preprocess proteomics data starting from loaded inputs
#'
#' This function performs all proteomics preprocessing steps required
#' prior to differential expression analysis, starting from already
#' loaded and validated inputs.
#'
#' The preprocessing workflow includes:
#'   1) Validation of proteomics inputs
#'   2) Construction of a sample-by-feature expression matrix from DIA-NN outputs
#'   3) Extraction of protein-level annotation (row data)
#'   4) Feature filtering based on minimum detection count across samples
#'   5) Missing value imputation (currently Perseus-style)
#'
#' The function performs no file I/O and is designed to be used in
#' reproducible pipelines (e.g., targets), where inputs are loaded
#' upstream and passed explicitly.
#'
#' @param inputs List of proteomics inputs as returned by
#'        load_proteomics_inputs(), including protein data, sample mapping,
#'        metadata, and contrasts.
#' @param config Full pipeline configuration list.
#'
#' @return A named list containing:
#'   \item{expr_raw}{Raw expression matrix (samples × proteins).}
#'   \item{expr_filt}{Filtered expression matrix after minimum-count filtering.}
#'   \item{expr_imp}{Imputed expression matrix used for downstream DE analysis.}
#'   \item{row_data}{Protein-level annotation corresponding to the rows of expr_* matrices.}
#'   \item{meta}{Sample metadata aligned to the expression matrices.}
#'
#' @seealso load_proteomics_inputs
#' @seealso filter_proteomics_by_min_count
#' @seealso impute_proteomics_perseus
#'
preprocess_proteomics <- function(inputs, config) {
  cfg <- config$modes$proteomics
  
  # validate inputs
  validate_proteomics_inputs(inputs, cfg)
  
  # Build standardized proteomics object (engine-aware)
  prot_obj <- get_proteomics_expression_matrix(inputs, config)
  
  # Downstream contract: work on log2 assay
  expr_mat  <- prot_obj$assay_log2
  row_data  <- prot_obj$row_data
  col_data  <- prot_obj$col_data
  
  
  # protein annotation
  row_data <- inputs$protein[
    ,
    c(cfg$id_columns$protein_id,
      unlist(cfg$id_columns$protein_annot)),
    drop = FALSE
  ]
  
  # filtering
  filt <- filter_proteomics_by_min_count(
    expr_mat = expr_mat,
    row_data = row_data,
    meta     = col_data,
    cfg      = cfg
  )
  
  expr_filt  <- filt$expr_mat
  row_data_f <- filt$row_data
  
  # imputation (currently: Perseus)
  expr_imp <- impute_proteomics_perseus(
    expr_mat = expr_filt,
    cfg      = cfg
  )
  
  list(
    expr_raw  = expr_mat,
    expr_filt = expr_filt,
    expr_imp  = expr_imp,
    row_data  = row_data_f,
    meta      = col_data
  )
}




