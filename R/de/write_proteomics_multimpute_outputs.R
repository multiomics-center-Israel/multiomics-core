#' Write proteomics multi-imputation DE outputs (legacy-compatible)
#'
#' Writes:
#'   1) A stability summary table (8/10 logic) as a single CSV file.
#'   2) Optionally, per-imputation per-contrast DE tables (CSV).
#'
#' Returns file paths for targets tracking (format = "file").
#'
#' @param de_res Output list returned by run_proteomics_de_mult_impute():
#'        - de_res$summary_df (data.frame)
#'        - de_res$runs_de_tables (list of length N)
#' @param out_dir Base output directory (default: config$paths$out / "proteomics").
#' @param prefix File prefix for outputs.
#' @param write_runs Logical; if TRUE, write per-imputation DE tables.
#'
#' @return Character vector of file paths written.
write_proteomics_multimpute_outputs <- function(de_res,
                                                out_dir,
                                                prefix = "proteomics",
                                                write_runs = FALSE) {
  stopifnot(is.list(de_res))
  stopifnot(is.data.frame(de_res$summary_df))
  stopifnot(is.list(de_res$runs_de_tables))
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  files <- character(0)
  
  # 1) Stability summary (single file)
  summary_file <- file.path(out_dir, paste0(prefix, "_limma_multimp_summary.csv"))
  utils::write.csv(de_res$summary_df, summary_file, row.names = FALSE, quote = FALSE)
  files <- c(files, summary_file)
  
  # 2) Optional: per-imputation DE tables (per contrast)
  if (isTRUE(write_runs)) {
    runs_dir <- file.path(out_dir, paste0(prefix, "_limma_multimp_runs"))
    dir.create(runs_dir, showWarnings = FALSE, recursive = TRUE)
    
    N <- length(de_res$runs_de_tables)
    
    for (i in seq_len(N)) {
      run_i <- de_res$runs_de_tables[[i]]
      stopifnot(is.list(run_i))
      
      for (cn in names(run_i)) {
        df <- run_i[[cn]]
        stopifnot(is.data.frame(df))
        
        f <- file.path(runs_dir, sprintf("%s_run%02d_%s.csv", prefix, i, cn))
        utils::write.csv(df, f, row.names = FALSE, quote = FALSE)
        files <- c(files, f)
      }
    }
  }
  
  files
}
