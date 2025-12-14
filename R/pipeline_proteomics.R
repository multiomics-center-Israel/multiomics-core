#' Full proteomics preprocessing pipeline (DIA-NN)
#'
#' @param config full config list from load_config()
#'
run_proteomics_preprocessing <- function(config) {
  inputs <- load_proteomics_inputs(config)
  preprocess_proteomics(inputs, config)
}


#' Run proteomics differential expression with multiple imputations and stability filtering
#'
#' This function implements the legacy Neat proteomics workflow for differential expression:
#'   1) Perform N independent imputations of the proteomics expression matrix
#'   2) Run limma differential expression analysis on each imputed dataset
#'   3) Aggregate results across imputations using a stability criterion
#'      (features must pass the DE cutoff in >= min_no_passed imputations)
#'
#' The stability summary reproduces the original pipeline behavior:
#'   - "8 out of 10" (or configurable) pass rule
#'   - Average linear fold-change across imputations
#'   - Quantile-based p-value and adjusted p-value summaries
#'
#' This function performs no file I/O and is suitable for use inside
#' a reproducible workflow (e.g., targets).
#'
#' @param expr_mat Numeric expression matrix to be used for DE analysis
#'        (typically filtered, normalized, and log2-transformed).
#' @param inputs List of proteomics inputs as returned by load_proteomics_inputs(),
#'        containing at least protein annotations, sample metadata, and contrasts.
#' @param config Full pipeline configuration list.
#'
#' @return A list with two elements:
#'   \item{runs_de_tables}{List of length N, where each element contains per-contrast
#'        limma DE tables for a single imputation run.}
#'   \item{summary_df}{A data.frame summarizing DE stability across imputations,
#'        containing legacy-compatible columns (sum.pass.*, pass.imputs.*, linearFC.imputs.*, etc.).}
#'
#' @seealso run_limma_mult_imputation_proteomics
#' @seealso summarize_mult_imputation
#'
run_proteomics_de_mult_impute <- function(expr_mat, inputs, config) {
  
  # Extract limma sub-config
  limma_cfg <- config$modes$proteomics$limma
  
  # Fail fast: stability threshold must not exceed number of imputations
  if (as.integer(limma_cfg$min_no_passed) > as.integer(limma_cfg$no_repetitions)) {
    stop("Invalid limma config: min_no_passed must be <= no_repetitions")
  }
  
  # Run N imputations and limma per imputation
  runs <- run_limma_mult_imputation_proteomics(
    expr_mat     = expr_mat,
    meta         = as.data.frame(inputs$metadata),
    contrasts_df = as.data.frame(inputs$contrasts),
    prot_tbl     = as.data.frame(inputs$protein),
    cfg          = config,
    seed_base    = config$params$seed
  )
  
  # Keep only DE tables per run (format expected by summarize_mult_imputation)
  runs_de_tables <- lapply(runs, `[[`, "de_tables")
  
  # Legacy stability summary ("passed in >= min_no_passed imputations")
  summary_df <- summarize_mult_imputation(
    runs_de_tables = runs_de_tables,
    cfg            = limma_cfg
  )
  
  list(
    runs_de_tables = runs_de_tables,
    summary_df     = summary_df
  )
}





