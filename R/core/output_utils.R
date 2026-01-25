#' Build output root directory for a specific project run
#'
#' Creates a stable "run folder" name to separate outputs between analysis rounds.
#' Example: outputs/Results_E_Pick_Analysis_02
#'
#' @param config Full config list.
#' @return Character path to the run output directory.
get_run_out_dir <- function(config) {
  proj  <- config$project$name %||% "Project"
  round <- config$project$analysis_round %||% "Analysis"
  resolve_out_path(config, sprintf("Results_%s_%s", proj, round))
}

#' Get output directory for a specific omics mode
#'
#' @param out_dir Base run directory (from get_run_out_dir)
#' @param mode    Omics mode name (e.g. "proteomics", "rnaseq")
#'
#' @return Path to mode-specific output directory
get_mode_out_dir <- function(out_dir, mode) {
  stopifnot(is.character(out_dir), length(out_dir) == 1)
  stopifnot(is.character(mode), length(mode) == 1)
  file.path(out_dir, mode)
}


#' Create legacy-style output folder structure for a run
#'
#' Mirrors the original Neat proteomics output tree inside the run folder:
#'   Datasets/, Diagnostic_plots/, Clustering/, Enrichment/, GSEA_enrichment/
#'
#' @param out_dir Run output root directory (e.g., outputs/Results_E_Pick_Analysis_02).
#' @return Named list of important subdirectories.
create_legacy_output_dirs <- function(out_dir, create = TRUE) {
  stopifnot(is.character(out_dir), length(out_dir) == 1)
  
  dirs <- list(
    datasets         = file.path(out_dir, "Datasets"),
    diagnostic_plots = file.path(out_dir, "Diagnostic_plots"),
    clustering       = file.path(out_dir, "Clustering"),
    enrichment       = file.path(out_dir, "Enrichment"),
    gsea_enrichment  = file.path(out_dir, "GSEA_enrichment")
  )
  
  if (isTRUE(create)) {
    for (d in unique(unname(dirs))) ensure_dir(d)
  }
  
  dirs
}


