# R/de/build_proteomics_de_results.R

#' Build standardized DE results for proteomics
#'
#' @param de_runs list of DE runs (one per imputation)
#' @param config full config list
#'
#' @return list with:
#'   - runs_de_tables
#'   - runs
#'   - summary_df
build_proteomics_de_results <- function(de_runs, config) {
  
  stopifnot(is.list(de_runs), length(de_runs) > 0)
  
  cfg <- config$modes$proteomics
  de_cfg <- cfg$de
  
  method <- de_cfg$method %||% "limma"
  
  # ---- extract DE tables (engine-agnostic) ----
  runs_de_tables <- lapply(de_runs, function(x) {
    stopifnot(!is.null(x$de_tables))
    x$de_tables
  })
  
  # ---- summarize according to DE method ----
  if (method == "limma") {
    
    summary_df <- summarize_limma_mult_imputation(
      runs_de_tables = runs_de_tables,
      config = config
    )
    
  } else {
    stop(sprintf(
      "DE method '%s' is not implemented for proteomics",
      method
    ))
  }
  
  list(
    runs_de_tables = runs_de_tables,
    runs           = de_runs,
    summary_df     = summary_df
  )
}
