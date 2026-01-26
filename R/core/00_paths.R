#' Ensure a directory exists
#'
#' @param dir Character path.
#' @return The path (invisibly).
ensure_dir <- function(dir) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    invisible(dir)
}

#' Resolve a path relative to the project root
#'
#' @param config Config object with `project$dir`.
#' @param ... Path components.
#' @return Absolute path.
resolve_project_path <- function(config, ...) {
    root <- config$project$dir %||% stop("config$project$dir is missing.")
    file.path(root, ...)
}

#' Resolve a path relative to the 'raw' data directory
resolve_raw_path <- function(config, rel_path) {
    base <- config$paths$raw %||% "data"
    resolve_project_path(config, base, rel_path)
}

#' Resolve a path relative to the 'outputs' directory
resolve_out_path <- function(config, rel_path = "") {
    base <- config$paths$out %||% "outputs"
    resolve_project_path(config, base, rel_path)
}

#' Build output root directory for a specific project run
get_run_out_dir <- function(config) {
    proj <- config$project$name %||% "Project"
    round <- config$project$analysis_round %||% "Analysis"
    resolve_out_path(config, sprintf("Results_%s_%s", proj, round))
}

#' Get output directory for a specific omics mode
get_mode_out_dir <- function(out_dir, mode) {
    stopifnot(is.character(out_dir), length(out_dir) == 1)
    stopifnot(is.character(mode), length(mode) == 1)
    file.path(out_dir, mode)
}

#' Create legacy-style output folder structure for a run
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

# Helper for defaults (needed here for %||% usage in resolve_*)
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}
