#' Extract log2-transformed DIA-NN measurements per sample (features × samples)
#'
#' This function:
#' - selects sample columns from the wide protein matrix (excluding annotation columns)
#' - coerces values to numeric (DIA-NN often stores blanks as character)
#' - applies log2 transformation
#' - renames columns from raw sample names to unified SampleID (using sample_map)
#' - reorders columns to match metadata sample order (effects$samples)
#' - sets rownames to the configured protein ID column (id_columns$protein_id)
#'
#' @param protein    Protein table (wide DIA-NN output).
#' @param sample_map Sample mapping table.
#' @param meta       Metadata table.
#' @param cfg        config$modes$proteomics.
#'
#' @return A numeric matrix of log2 intensities (features × samples), with:
#'   - rownames = protein IDs (e.g. Protein.Group)\n
#'   - colnames = unified sample IDs (e.g. SampleID), ordered to match metadata.
get_measurements_per_sample_diann <- function(protein, sample_map, meta, cfg) {
  id_cols  <- cfg$id_columns
  eff_cols <- cfg$effects
  
  # 0) Protein ID for rownames
  protein_id_col <- id_cols$protein_id
  check_has_cols(protein, protein_id_col, df_name = "protein")
  feat_ids <- as.character(protein[[protein_id_col]])
  
  if (anyNA(feat_ids) || any(feat_ids == "")) {
    stop(sprintf("Protein ID column '%s' contains NA/empty values.", protein_id_col))
  }
  if (anyDuplicated(feat_ids) > 0) {
    dups <- unique(feat_ids[duplicated(feat_ids)])
    stop(sprintf("Protein ID column '%s' contains duplicated IDs (e.g. %s).",
                 protein_id_col, paste(head(dups, 3), collapse = ", ")))
  }
  
  # 1) Identify non-sample columns in the protein table
  annot_cols <- id_cols$protein_annot
  annot_cols <- if (is.null(annot_cols)) character(0) else unlist(annot_cols)
  non_sample_cols <- unique(c(protein_id_col, annot_cols))
  
  # Sample columns = everything else
  sample_cols <- setdiff(colnames(protein), non_sample_cols)
  if (length(sample_cols) == 0) {
    stop("No sample columns detected in protein table after excluding ID/annotation columns.")
  }
  
  df_m <- protein[, sample_cols, drop = FALSE]
  
  # 2) Coerce to numeric (DIA-NN often has blanks)
  df_m <- as.data.frame(
    suppressWarnings(sapply(df_m, function(x) as.numeric(as.character(x)))),
    check.names = FALSE
  )
  
  # 3) Log2 transformation
  df_m <- as.data.frame(
    log2(df_m),
    check.names = FALSE
  )
  
  # 4) Rename columns: raw sample names -> SampleID using sample_map
  raw_names <- colnames(df_m)
  map_from <- id_cols$map_from
  map_to   <- id_cols$map_to
  
  check_has_cols(sample_map, c(map_from, map_to), df_name = "sample_map")
  
  new_names <- sample_map[[map_to]][match(raw_names, sample_map[[map_from]])]
  
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
      paste(head(missing_in_df, 10), collapse = ", "),
      if (length(missing_in_df) > 10) sprintf(" ... (+%d more)", length(missing_in_df) - 10) else ""
    )
  }
  
  ordered_cols <- intersect(meta_sample_ids, colnames(df_m))
  df_m_ordered <- df_m[, ordered_cols, drop = FALSE]
  
  # 6) Return as numeric matrix with proper rownames
  m <- coerce_df_to_numeric_matrix(df_m_ordered, rownames_vec = feat_ids, name = "diann_measurements")
  
  # Optional: enforce strict contract
  assert_numeric_matrix(m, "diann_measurements")
  
  m
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
  
  # Downstream contract: work on log2 assay (features x samples)
  expr_raw <- prot_obj$assay_log2
  row_data <- prot_obj$row_data
  col_data <- prot_obj$col_data
  
  # ---- enforce matrix + IDs contract early ----
  # Ensure expr_raw is a numeric matrix with proper row/col names.
  # Prefer protein_id column as rownames if needed.
  protein_id_col <- cfg$id_columns$protein_id
  
  # If expr_raw is a data.frame, coerce to matrix
  if (is.data.frame(expr_raw)) {
    rn <- NULL
    if (!is.null(row_data) && protein_id_col %in% colnames(row_data)) {
      rn <- row_data[[protein_id_col]]
    }
    expr_raw <- coerce_df_to_numeric_matrix(expr_raw, rownames_vec = rn, name = "expr_raw")
  }
  
  # If rownames are missing/meaningless, set from row_data protein_id
  if (!is.null(row_data) && protein_id_col %in% colnames(row_data)) {
    rn <- as.character(row_data[[protein_id_col]])
    # replace if empty or looks like 1..N
    if (is.null(rownames(expr_raw)) || all(grepl("^\\d+$", rownames(expr_raw)))) {
      rownames(expr_raw) <- rn
    }
  }
  
  assert_numeric_matrix(expr_raw, "expr_raw")
  
  # Ensure col_data is aligned to expr_raw columns (same sample IDs)
  sample_id_col <- cfg$effects$samples %||% cfg$id_columns$sample_col
  check_has_cols(col_data, sample_id_col, df_name = "col_data")
  
  # Reorder col_data to match expr_raw colnames (fail if missing)
  idx <- match(colnames(expr_raw), col_data[[sample_id_col]])
  if (anyNA(idx)) {
    missing <- colnames(expr_raw)[is.na(idx)]
    stop(sprintf("col_data is missing %d samples present in expr_raw, e.g. %s",
                 length(missing), paste(head(missing, 3), collapse = ", ")))
  }
  col_data <- col_data[idx, , drop = FALSE]
  
  # ---- filtering ----
  filt <- filter_proteomics_by_min_count(
    expr_mat = expr_raw,
    row_data = row_data,
    meta     = col_data,
    cfg      = cfg
  )
  
  expr_filt <- filt$expr_mat
  row_data_f <- filt$row_data
  
  assert_numeric_matrix(expr_filt, "expr_filt")
  
  # ---- optional: single-imputation + QC only (kept for plots) ----
  imp_res <- impute_proteomics_perseus(
    expr_mat     = expr_filt,
    cfg          = cfg,
    return_flags = TRUE
  )
  expr_imp_single <- imp_res$imputed
  assert_numeric_matrix(expr_imp_single, "expr_imp_single")
  
  imputation_qc <- NULL
  if (!is.null(imp_res$imputed_flag)) {
    imputation_qc <- list(
      imputed_flag = imp_res$imputed_flag,
      hist_plot    = plot_imputation_summary(
        expr_mat      = expr_filt,
        imputed_flag  = imp_res$imputed_flag
      )
    )
  }
  
  # Return: matrices only, consistently named
  list(
    expr_raw_mat     = expr_raw,
    expr_filt_mat    = expr_filt,
    expr_imp_single  = expr_imp_single,
    row_data         = row_data_f,
    meta             = col_data,
    imputation_qc    = imputation_qc
  )
}





