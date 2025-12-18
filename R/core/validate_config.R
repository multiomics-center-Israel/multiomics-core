#' Validate global pipeline configuration
#'
#' Runs fail-fast checks on the YAML-derived config object.
#' This is intended to catch invalid parameter combinations before
#' running long computations (e.g., multi-imputation limma).
#'
#' @param config Full config list (from config.yaml).
#' @return Invisibly TRUE if all checks pass; otherwise stops with an error.
validate_config <- function(config) {
  stopifnot(is.list(config))
  stopifnot(!is.null(config$paths))
  stopifnot(!is.null(config$modes))
  stopifnot(!is.null(config$params))
  
  # Basic required fields
  stopifnot(!is.null(config$paths$out))
  stopifnot(!is.null(config$params$seed))
  
  # Mode-specific validations (add more over time)
  if (!is.null(config$modes$proteomics)) {
    validate_proteomics_config(config$modes$proteomics)
  }
  if (!is.null(config$modes$rna)) {
    validate_rna_config(config$modes$rna)
  }
  if (!is.null(config$modes$metabolomics)) {
    validate_metabolomics_config(config$modes$metabolomics)
  }
  if (!is.null(config$modes$lipidomics)) {
    validate_lipidomics_config(config$modes$lipidomics)
  }
  
  invisible(TRUE)
}


#' Validate proteomics configuration parameters
#'
#' Checks internal consistency of proteomics-related config fields.
#' These are "fail-fast" checks to prevent running long pipelines with invalid parameters.
#'
#' @param cfg Proteomics mode config (config$modes$proteomics).
#' @return Invisibly TRUE if all checks pass; otherwise stops with an error.
validate_proteomics_config <- function(cfg) {
  stopifnot(is.list(cfg))
  
  # Engine + scale contract (future-proof for DIANN/MaxQuant)
  if (!is.null(cfg$engine)) {
    stopifnot(is.character(cfg$engine), length(cfg$engine) == 1)
  }
  if (!is.null(cfg$scale_in)) {
    stopifnot(cfg$scale_in %in% c("linear", "log2"))
  }
  if (!is.null(cfg$transform) && !is.null(cfg$transform$target_scale)) {
    stopifnot(cfg$transform$target_scale %in% c("log2"))
  }
  
  # --- ID columns (NEW: ensure sample_col exists and is consistent) ---
  stopifnot(!is.null(cfg$id_columns))
  stopifnot(!is.null(cfg$id_columns$protein_id))
  
  # You added this:
  stopifnot(!is.null(cfg$id_columns$sample_col))
  stopifnot(is.character(cfg$id_columns$sample_col), length(cfg$id_columns$sample_col) == 1)
  
  # Recommended: keep consistent with map_to if map_to exists
  if (!is.null(cfg$id_columns$map_to)) {
    stopifnot(is.character(cfg$id_columns$map_to), length(cfg$id_columns$map_to) == 1)
    stopifnot(identical(cfg$id_columns$sample_col, cfg$id_columns$map_to))
  }
  
  # Imputation parameters
  if (!is.null(cfg$imputation) && !is.null(cfg$imputation$method)) {
    stopifnot(cfg$imputation$method %in% c("none", "perseus", "dep2"))
    
    # NOTE: currently using your existing keys (no_repetitions/min_no_passed)
    n_rep <- as.integer(cfg$imputation$no_repetitions)
    n_min <- as.integer(cfg$imputation$min_no_passed)
    
    stopifnot(!is.na(n_rep), n_rep >= 1)
    stopifnot(!is.na(n_min), n_min >= 1)
    stopifnot(n_min <= n_rep)
    
    if (identical(cfg$imputation$method, "perseus")) {
      width <- as.numeric(cfg$imputation$width)
      down  <- as.numeric(cfg$imputation$downshift)
      stopifnot(!is.na(width), width > 0)
      stopifnot(!is.na(down),  down > 0)
    }
    
    if (identical(cfg$imputation$method, "dep2")) {
      stopifnot(!is.null(cfg$imputation$dep2_method))
      stopifnot(!is.null(cfg$imputation$dep2_random_seed))
    }
  }
  
  # Limma settings
  if (!is.null(cfg$limma)) {
    
    # Optional sanity checks
    if (!is.null(cfg$limma$p_cutoff)) {
      p <- as.numeric(cfg$limma$p_cutoff)
      stopifnot(!is.na(p), p > 0, p <= 1)
    }
    if (!is.null(cfg$limma$linear_fc_cutoff)) {
      fc <- as.numeric(cfg$limma$linear_fc_cutoff)
      stopifnot(!is.na(fc), fc >= 1)
    }
  }
  
  invisible(TRUE)
}


#' Validate consistency of proteomics imputations and meta/sample IDs
#'
#' Runtime validation (not config): checks that imputations are consistent and
#' that `meta` matches expression sample IDs using
#' `cfg$modes$proteomics$id_columns$sample_col`.
#'
#' Keep this in validate_config.R if you prefer a single file, but call it
#' from analysis code (e.g., DE runner), because it needs runtime objects.
#'
#' @param imputations list of imputed expression matrices
#' @param meta sample metadata data.frame
#' @param cfg full config list (must include modes$proteomics$id_columns$sample_col)
#' @param warn_extra_meta logical; warn if meta contains samples not present in matrix
#' @return invisibly TRUE; otherwise stops with informative error
validate_proteomics_imputations <- function(imputations,
                                            meta,
                                            cfg,
                                            warn_extra_meta = TRUE) {
  
  if (!is.list(imputations) || length(imputations) < 1) {
    stop("`imputations` must be a non-empty list of matrices.")
  }
  
  ref <- imputations[[1]]
  if (!is.matrix(ref)) stop("imputations[[1]] must be a matrix.")
  
  ref_dim <- dim(ref)
  ref_rn  <- rownames(ref)
  ref_cn  <- colnames(ref)
  
  if (is.null(ref_rn) || is.null(ref_cn)) {
    stop("Imputed matrices must have rownames and colnames (features x samples).")
  }
  
  # uniqueness checks
  if (anyDuplicated(ref_cn) > 0) {
    dups <- unique(ref_cn[duplicated(ref_cn)])
    stop(sprintf("Expression matrix has duplicated sample IDs in colnames (e.g. %s).",
                 paste(head(dups, 3), collapse = ", ")))
  }
  if (anyDuplicated(ref_rn) > 0) {
    dups <- unique(ref_rn[duplicated(ref_rn)])
    stop(sprintf("Expression matrix has duplicated feature IDs in rownames (e.g. %s).",
                 paste(head(dups, 3), collapse = ", ")))
  }
  
  for (i in seq_along(imputations)) {
    x <- imputations[[i]]
    
    if (!is.matrix(x)) stop(sprintf("Imputation %d is not a matrix.", i))
    if (!identical(dim(x), ref_dim)) {
      stop(sprintf("Imputation %d has different dimensions. Expected %s, got %s.",
                   i, paste(ref_dim, collapse = "x"), paste(dim(x), collapse = "x")))
    }
    if (!identical(rownames(x), ref_rn)) stop(sprintf("Imputation %d has different rownames.", i))
    if (!identical(colnames(x), ref_cn)) stop(sprintf("Imputation %d has different colnames.", i))
    if (!is.numeric(x)) stop(sprintf("Imputation %d is not numeric (storage: %s).", i, typeof(x)))
    if (anyNA(x)) stop(sprintf("Imputation %d still contains NA values.", i))
    if (any(!is.finite(x))) stop(sprintf("Imputation %d contains non-finite values (Inf/-Inf/NaN).", i))
  }
  
  sample_col <- cfg$modes$proteomics$id_columns$sample_col
  if (is.null(sample_col) || !nzchar(sample_col)) {
    stop("Config missing `cfg$modes$proteomics$id_columns$sample_col` (needed to match meta to samples).")
  }
  
  if (!is.data.frame(meta)) stop("`meta` must be a data.frame.")
  if (!(sample_col %in% colnames(meta))) {
    stop(sprintf("`meta` is missing sample ID column '%s' (from config).", sample_col))
  }
  
  meta_samples <- as.character(meta[[sample_col]])
  
  if (anyDuplicated(meta_samples) > 0) {
    dups <- unique(meta_samples[duplicated(meta_samples)])
    stop(sprintf("`meta[[%s]]` contains duplicated sample IDs (e.g. %s).",
                 sample_col, paste(head(dups, 3), collapse = ", ")))
  }
  
  missing_in_meta <- setdiff(ref_cn, meta_samples)
  if (length(missing_in_meta) > 0) {
    stop(sprintf("meta is missing %d samples present in expression matrix (sample_col='%s'), e.g. %s",
                 length(missing_in_meta), sample_col, paste(head(missing_in_meta, 3), collapse = ", ")))
  }
  
  extra_in_meta <- setdiff(meta_samples, ref_cn)
  if (length(extra_in_meta) > 0 && isTRUE(warn_extra_meta)) {
    warning(sprintf("meta contains %d samples not present in expression matrix (sample_col='%s'), e.g. %s",
                    length(extra_in_meta), sample_col, paste(head(extra_in_meta, 3), collapse = ", ")))
  }
  
  invisible(TRUE)
}


# --- Stubs for future modes (no-op for now) ---
validate_rna_config <- function(cfg) invisible(TRUE)
validate_metabolomics_config <- function(cfg) invisible(TRUE)
validate_lipidomics_config <- function(cfg) invisible(TRUE)
