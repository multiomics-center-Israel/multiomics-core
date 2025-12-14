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
  
  # Imputation parameters
  if (!is.null(cfg$imputation) && !is.null(cfg$imputation$method)) {
    stopifnot(cfg$imputation$method %in% c("none", "perseus", "dep2"))
    
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
  
  # Limma multi-imputation settings (your key check)
  if (!is.null(cfg$limma)) {
    n_rep <- as.integer(cfg$limma$no_repetitions)
    n_min <- as.integer(cfg$limma$min_no_passed)
    stopifnot(!is.na(n_rep), n_rep >= 1)
    stopifnot(!is.na(n_min), n_min >= 1)
    stopifnot(n_min <= n_rep)
    
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


# --- Stubs for future modes (no-op for now) ---
validate_rna_config <- function(cfg) invisible(TRUE)
validate_metabolomics_config <- function(cfg) invisible(TRUE)
validate_lipidomics_config <- function(cfg) invisible(TRUE)
