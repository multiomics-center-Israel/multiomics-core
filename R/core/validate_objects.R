#' Validate a list of imputed matrices + metadata alignment
#'
#' @param imputations list of numeric matrices (features x samples)
#' @param meta data.frame with one row per sample
#' @param sample_col column name in meta holding sample IDs
#' @param warn_extra_meta warn if meta contains samples not in matrix
#' @param allow_na allow NA values in matrices
validate_imputations <- function(imputations,
                                 meta,
                                 sample_col,
                                 warn_extra_meta = TRUE,
                                 allow_na = FALSE) {
  if (!is.list(imputations) || length(imputations) < 1) {
    stop("imputations must be a non-empty list of matrices.")
  }
  if (!is.data.frame(meta)) stop("meta must be a data.frame.")
  if (is.null(sample_col) || !nzchar(sample_col)) stop("sample_col must be provided.")
  if (!(sample_col %in% colnames(meta))) stop("meta is missing sample_col: ", sample_col)
  
  ref <- imputations[[1]]
  assert_numeric_matrix(ref, "imputations[[1]]")  # rownames+colnames+unique
  
  ref_dim <- dim(ref)
  ref_rn  <- rownames(ref)
  ref_cn  <- colnames(ref)
  
  # ---- validate each imputation matrix ----
  for (i in seq_along(imputations)) {
    x <- imputations[[i]]
    if (!is.matrix(x)) stop(sprintf("Imputation %d is not a matrix.", i))
    if (!is.numeric(x)) stop(sprintf("Imputation %d is not numeric (storage: %s).", i, typeof(x)))
    
    if (!identical(dim(x), ref_dim)) {
      stop(sprintf(
        "Imputation %d has different dimensions. Expected %s, got %s.",
        i, paste(ref_dim, collapse = "x"), paste(dim(x), collapse = "x")
      ))
    }
    if (!identical(rownames(x), ref_rn)) {
      j <- which(rownames(x) != ref_rn)[1]
      stop(sprintf(
        "Imputation %d has different rownames/order at position %d: expected '%s' got '%s'.",
        i, j, ref_rn[j], rownames(x)[j]
      ))
    }
    if (!identical(colnames(x), ref_cn)) {
      j <- which(colnames(x) != ref_cn)[1]
      stop(sprintf(
        "Imputation %d has different colnames/order at position %d: expected '%s' got '%s'.",
        i, j, ref_cn[j], colnames(x)[j]
      ))
    }
    
    if (!isTRUE(allow_na) && anyNA(x)) {
      stop(sprintf("Imputation %d still contains NA values.", i))
    }
    if (!all(is.finite(x[is.na(x) == FALSE]))) {
      stop(sprintf("Imputation %d contains non-finite values (Inf/-Inf/NaN).", i))
    }
  }
  
  # ---- meta checks vs matrix ----
  meta_ids <- as.character(meta[[sample_col]])
  if (anyDuplicated(meta_ids) > 0) {
    dups <- unique(meta_ids[duplicated(meta_ids)])
    stop("meta[[", sample_col, "]] has duplicated sample IDs (e.g. ", paste(head(dups, 3), collapse = ", "), ").")
  }
  
  missing_in_meta <- setdiff(ref_cn, meta_ids)
  if (length(missing_in_meta) > 0) {
    stop(sprintf(
      "meta is missing %d samples present in expression matrix (sample_col='%s'), e.g. %s",
      length(missing_in_meta), sample_col, paste(head(missing_in_meta, 3), collapse = ", ")
    ))
  }
  
  extra_in_meta <- setdiff(meta_ids, ref_cn)
  if (length(extra_in_meta) > 0 && isTRUE(warn_extra_meta)) {
    warning(sprintf(
      "meta contains %d samples not present in expression matrix (sample_col='%s'), e.g. %s",
      length(extra_in_meta), sample_col, paste(head(extra_in_meta, 3), collapse = ", ")
    ))
  }
  
  invisible(TRUE)
}

#' Validate consistency of proteomics imputations and metadata sample IDs
#'
#' Fail-fast validation for runtime objects (not file I/O, not config schema).
#' This is designed to catch silent misalignment bugs before DE/QC steps.
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
  if (is.null(sample_col) || !nzchar(sample_col)) {
    # בדיוק כמו שעשיתם סטנדרט חדש:
    p <- cfg$modes$proteomics
    sample_col <- p$effects$samples %||% p$id_columns$sample_col %||% "SampleID"
  }
  
  validate_imputations(
    imputations     = imputations,
    meta            = meta,
    sample_col      = sample_col,
    warn_extra_meta = warn_extra_meta,
    allow_na        = allow_na
  )
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
#' assert_numeric_matrix(expr_mat, "expr_filt")
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

assert_pre_contract <- function(pre, stage = "proteomics") {
  stopifnot(is.list(pre))
  
  required <- c("expr_raw", "expr_filt", "expr_imp_single", "meta", "row_data")
  missing <- setdiff(required, names(pre))
  if (length(missing) > 0) {
    stop(sprintf(
      "Preprocess contract failed for %s. Missing fields: %s",
      stage, paste(missing, collapse = ", ")
    ))
  }
  
  invisible(TRUE)
}

assert_de_contract <- function(de_res, stage = "proteomics") {
  stopifnot(is.list(de_res))
  
  required <- c("summary_df")
  missing <- setdiff(required, names(de_res))
  if (length(missing) > 0) {
    stop(sprintf(
      "DE contract failed for %s. Missing fields: %s",
      stage, paste(missing, collapse = ", ")
    ))
  }
  
  invisible(TRUE)
}
