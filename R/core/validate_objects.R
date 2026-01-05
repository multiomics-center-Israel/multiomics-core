#' Validate consistency of proteomics imputations and metadata sample IDs
#'
#' Fail-fast validation for runtime objects (not file I/O, not config schema).
#' This is designed to catch silent misalignment bugs before DE/QC steps.
#'
#' Checks:
#' - `imputations` is a non-empty list of numeric matrices
#' - all matrices share identical dim/rownames/colnames (including order)
#' - (optionally) no NA values remain after imputation
#' - no non-finite values (Inf/-Inf/NaN)
#' - metadata contains a sample ID column (cfg$modes$proteomics$id_columns$sample_col)
#' - metadata covers all matrix columns (warn optionally on extra meta samples)
#'
#' @param imputations list of imputed expression matrices (features x samples)
#' @param meta data.frame with one row per sample
#' @param cfg full config list (expects modes$proteomics$id_columns$sample_col)
#' @param warn_extra_meta logical; warn if meta contains samples not present in matrix
#' @param allow_na logical; if TRUE, allow NA values in matrices (useful pre-imputation / NA-QC)
#'
#' @return invisibly TRUE; otherwise stops with informative error

validate_proteomics_imputations <- function(imputations,
                                            meta,
                                            cfg,
                                            warn_extra_meta = TRUE,
                                            allow_na = FALSE,
                                            sample_col = NULL) {
  # ---- basic checks ----
  if (is.null(sample_col) || !nzchar(sample_col)) {
    sample_col <- cfg$modes$proteomics$effects$samples %||%
      cfg$modes$proteomics$id_columns$sample_col %||%
      "SampleID"
  }
  
  ref <- imputations[[1]]
  if (!is.matrix(ref)) stop("imputations[[1]] must be a matrix.")
  
  ref_dim <- dim(ref)
  ref_rn  <- rownames(ref)
  ref_cn  <- colnames(ref)
  
  if (is.null(ref_rn) || is.null(ref_cn)) {
    stop("Imputed matrices must have rownames and colnames (features x samples).")
  }
  
  # uniqueness checks (catch silent misalignment bugs)
  if (anyDuplicated(ref_cn) > 0) {
    dups <- unique(ref_cn[duplicated(ref_cn)])
    stop(sprintf(
      "Expression matrix has duplicated sample IDs in colnames (e.g. %s).",
      paste(head(dups, 3), collapse = ", ")
    ))
  }
  if (anyDuplicated(ref_rn) > 0) {
    dups <- unique(ref_rn[duplicated(ref_rn)])
    stop(sprintf(
      "Expression matrix has duplicated feature IDs in rownames (e.g. %s).",
      paste(head(dups, 3), collapse = ", ")
    ))
  }
  
  # ---- validate each imputation matrix ----
  for (i in seq_along(imputations)) {
    x <- imputations[[i]]
    
    if (!is.matrix(x)) stop(sprintf("Imputation %d is not a matrix.", i))
    
    if (!identical(dim(x), ref_dim)) {
      stop(sprintf(
        "Imputation %d has different dimensions. Expected %s, got %s.",
        i, paste(ref_dim, collapse = "x"), paste(dim(x), collapse = "x")
      ))
    }
    
    # rownames mismatch: report first mismatch
    if (!identical(rownames(x), ref_rn)) {
      j <- which(rownames(x) != ref_rn)[1]
      stop(sprintf(
        "Imputation %d has different rownames/order at position %d: expected '%s' got '%s'.",
        i, j, ref_rn[j], rownames(x)[j]
      ))
    }
    
    # colnames mismatch: report first mismatch
    if (!identical(colnames(x), ref_cn)) {
      j <- which(colnames(x) != ref_cn)[1]
      stop(sprintf(
        "Imputation %d has different colnames/order at position %d: expected '%s' got '%s'.",
        i, j, ref_cn[j], colnames(x)[j]
      ))
    }
    
    if (!is.numeric(x)) {
      stop(sprintf("Imputation %d is not numeric (storage: %s).", i, typeof(x)))
    }
    
    if (!isTRUE(allow_na) && anyNA(x)) {
      stop(sprintf("Imputation %d still contains NA values.", i))
    }
    
    # even if allow_na=TRUE, still guard against Inf/NaN (these break downstream)
    if (!all(is.finite(x[is.na(x) == FALSE]))) {
      stop(sprintf("Imputation %d contains non-finite values (Inf/-Inf/NaN).", i))
    }
  }
  
  # keep the resolved sample_col from the top (or override only if still missing)
  if (is.null(sample_col) || !nzchar(sample_col)) {
    sample_col <- cfg$modes$proteomics$id_columns$sample_col %||%
      cfg$modes$proteomics$id_columns$map_to %||%
      stop("Could not determine sample ID column. Please set cfg$modes$proteomics$id_columns$sample_col.")
  }
  
  # backward-compatible fallback (optional but safe)
  if (is.null(sample_col) || !nzchar(sample_col)) {
    sample_col <- cfg$modes$proteomics$id_columns$map_to
  }
  if (is.null(sample_col) || !nzchar(sample_col)) {
    stop("Could not determine sample ID column. Please set cfg$modes$proteomics$id_columns$sample_col.")
  }
  
  if (!is.data.frame(meta)) stop("`meta` must be a data.frame.")
  if (!(sample_col %in% colnames(meta))) {
    stop(sprintf("`meta` is missing sample ID column '%s' (from config).", sample_col))
  }
  
  meta_samples <- as.character(meta[[sample_col]])
  
  if (anyDuplicated(meta_samples) > 0) {
    dups <- unique(meta_samples[duplicated(meta_samples)])
    stop(sprintf(
      "`meta[[%s]]` contains duplicated sample IDs (e.g. %s).",
      sample_col, paste(head(dups, 3), collapse = ", ")
    ))
  }
  
  missing_in_meta <- setdiff(ref_cn, meta_samples)
  if (length(missing_in_meta) > 0) {
    stop(sprintf(
      "meta is missing %d samples present in expression matrix (sample_col='%s'), e.g. %s",
      length(missing_in_meta), sample_col, paste(head(missing_in_meta, 3), collapse = ", ")
    ))
  }
  
  extra_in_meta <- setdiff(meta_samples, ref_cn)
  if (length(extra_in_meta) > 0 && isTRUE(warn_extra_meta)) {
    warning(sprintf(
      "meta contains %d samples not present in expression matrix (sample_col='%s'), e.g. %s",
      length(extra_in_meta), sample_col, paste(head(extra_in_meta, 3), collapse = ", ")
    ))
  }
  
  invisible(TRUE)
}
#' Assert a strict numeric matrix contract (features × samples)
#'
#' Fail-fast validator for expression matrices used throughout the pipeline.
#' Ensures the object is a numeric matrix with **non-empty, unique** row and
#' column names. This catches silent alignment bugs early (e.g. lost rownames,
#' duplicated sample IDs, character coercions).
#'
#' @param x An object expected to be a numeric matrix.
#' @param name Character. A label used in error messages (helps debugging).
#'
#' @return Invisibly returns TRUE if all checks pass.
#'
#' @examples
#' \dontrun{
#' assert_numeric_matrix(expr_mat, "expr_filt_mat")
#' }
assert_numeric_matrix <- function(x, name = "x") {
  if (!is.matrix(x)) stop(sprintf("`%s` must be a matrix.", name))
  if (!is.numeric(x)) stop(sprintf("`%s` must be numeric.", name))
  if (is.null(rownames(x)) || anyNA(rownames(x)) || any(rownames(x) == "")) {
    stop(sprintf("`%s` must have non-empty rownames (feature IDs).", name))
  }
  if (is.null(colnames(x)) || anyNA(colnames(x)) || any(colnames(x) == "")) {
    stop(sprintf("`%s` must have non-empty colnames (sample IDs).", name))
  }
  if (anyDuplicated(rownames(x)) > 0) stop(sprintf("`%s` has duplicated rownames.", name))
  if (anyDuplicated(colnames(x)) > 0) stop(sprintf("`%s` has duplicated colnames.", name))
  invisible(TRUE)
}
#' Coerce a data.frame to a numeric matrix (optionally set rownames)
#'
#' Converts a data.frame (or already-matrix) to a matrix with storage mode
#' forced to numeric. This is useful for I/O sources (e.g. DIA-NN exports)
#' where columns may be read as character due to blanks or formatting.
#'
#' If `rownames_vec` is provided, it is used to set matrix rownames (feature IDs),
#' and its length must exactly match the number of rows.
#'
#' @param df A data.frame or matrix.
#' @param rownames_vec Optional vector of feature IDs to use as rownames.
#' @param name Character. A label used in error messages (helps debugging).
#'
#' @return A numeric matrix.
#'
#' @examples
#' \dontrun{
#' m <- coerce_df_to_numeric_matrix(df, rownames_vec = row_data$Protein.Group)
#' }
coerce_df_to_numeric_matrix <- function(df, rownames_vec = NULL, name = "df") {
  if (is.matrix(df)) return(df)
  if (!is.data.frame(df)) stop(sprintf("`%s` must be a data.frame or matrix.", name))
  
  m <- as.matrix(df)
  suppressWarnings(storage.mode(m) <- "numeric")
  
  if (!is.null(rownames_vec)) {
    if (length(rownames_vec) != nrow(m)) {
      stop(sprintf("`rownames_vec` length (%d) != nrow(%s) (%d).",
                   length(rownames_vec), name, nrow(m)))
    }
    rownames(m) <- as.character(rownames_vec)
  }
  m
}
