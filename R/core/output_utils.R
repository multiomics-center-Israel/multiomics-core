#' Build output root directory for a specific project run
#'
#' Creates a stable "run folder" name to separate outputs between analysis rounds.
#' Example: outputs/Results_E_Pick_Analysis_02
#'
#' @param config Full config list.
#' @return Character path to the run output directory.
get_run_out_dir <- function(config) {
  out_base <- config$paths$out %||% "outputs"
  proj <- config$project$name %||% "Project"
  round <- config$project$analysis_round %||% "Analysis"
  file.path(out_base, sprintf("Results_%s_%s", proj, round))
}
#' Create legacy-style output folder structure for a run
#'
#' Mirrors the original Neat proteomics output tree inside the run folder:
#'   Datasets/, Diagnostic_plots/, Clustering/, Enrichment/, GSEA_enrichment/
#'
#' @param run_dir Run output root directory (e.g., outputs/Results_E_Pick_Analysis_02).
#' @return Named list of important subdirectories.
create_legacy_output_dirs <- function(run_dir) {
  dirs <- list(
    run_dir          = run_dir,
    datasets         = file.path(run_dir, "Datasets"),
    diagnostic_plots = file.path(run_dir, "Diagnostic_plots"),
    clustering       = file.path(run_dir, "Clustering"),
    enrichment       = file.path(run_dir, "Enrichment"),
    gsea_enrichment  = file.path(run_dir, "GSEA_enrichment")
  )
  for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)
  dirs
}
#' Write a table twice: legacy .txt (TSV) + improved .tsv
#'
#' Legacy outputs keep the original filenames used by the old pipeline (usually *.txt),
#' but are written as tab-separated text for robustness.
#' Improved outputs use clearer filenames and the .tsv extension.
#'
#' @param x data.frame or matrix to write
#' @param legacy_path Full path to legacy file (typically *.txt)
#' @param improved_path Full path to improved file (typically *.tsv)
#' @return Character vector of written file paths (legacy, improved)
write_dual_tsv <- function(x, legacy_path, improved_path) {
  dir.create(dirname(legacy_path), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(improved_path), showWarnings = FALSE, recursive = TRUE)
  
  write.table(x, legacy_path, sep = "\t", quote = FALSE, col.names = NA)
  write.table(x, improved_path, sep = "\t", quote = FALSE, col.names = NA)
  
  c(legacy_path, improved_path)
}

