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

#' Resolve a user-supplied input path (absolute, or relative to raw/)
#'
#' Locates an input file the same way the data files are located: an absolute
#' path is returned unchanged, while a relative path is resolved against the raw
#' data directory (\code{config$paths$raw} under \code{project$dir}). "Absolute"
#' covers Unix (\code{/...}), home (\code{~/...}), Windows drive
#' (\code{C:/...}, \code{C:\\...}) and UNC (\code{\\\\host\\share}) paths, so a
#' config authored on either OS behaves the same. Accepts a vector (each element
#' resolved independently) and returns \code{NULL} for \code{NULL} input, so it
#' can wrap optional keys like \code{gmt_file} / \code{mapping_file} directly.
#'
#' @param config Config object with \code{project$dir} (and optional \code{paths$raw}).
#' @param path Character path(s), or NULL.
#' @return Resolved path(s) as a character vector, or NULL.
resolve_input_path <- function(config, path) {
    if (is.null(path)) return(NULL)
    path <- unlist(path, use.names = FALSE)
    if (length(path) == 0) return(path)
    vapply(path, function(p) {
        if (is.na(p) || !nzchar(p)) return(p)
        # Absolute if it starts with a drive letter ("C:"), or with /, ~ or \.
        if (grepl("^([A-Za-z]:|[/~\\\\])", p)) p else resolve_raw_path(config, p)
    }, character(1), USE.NAMES = FALSE)
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
