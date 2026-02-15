#' Load pipeline configuration from YAML
load_config <- function(path = NULL) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
        stop("Package 'yaml' is required.")
    }
    cfg <- yaml::read_yaml(path)
    cfg$.config_path <- normalizePath(path, mustWork = FALSE)
    cfg$.config_mtime <- file.info(path)$mtime
    cfg
}

#' Validate global pipeline configuration
validate_config <- function(config) {
    # Basic structure checks
    if (is.null(config$paths)) stop("config$paths is missing")
    if (is.null(config$modes)) stop("config$modes is missing")
    if (is.null(config$params)) stop("config$params is missing")

    if (!is.null(config$modes$proteomics)) validate_proteomics_config(config$modes$proteomics)
    if (!is.null(config$modes$rna)) validate_rna_config(config$modes$rna)
    if (!is.null(config$modes$metabolomics)) validate_metabolomics_config(config$modes$metabolomics)

    invisible(TRUE)
}

# (Note: validate_proteomics_config and validate_rna_config are dispatch hooks
# that should be defined in their respective domain files or sourced gracefully.
# We assume they will be available in the environment.)

get_mode_cfg <- function(config, mode) {
    cfg <- config$modes[[mode]]
    if (is.null(cfg)) stop("Missing config$modes$", mode)
    cfg
}

#' Extract primary color column name from config
#'
#' Handles both scalar and array config formats for effects$color.
#' Returns NULL if color config is missing.
#'
#' @param cfg Mode config list (e.g., config$modes$proteomics)
#' @return Character string (column name) or NULL
get_color_config <- function(cfg) {
    color_config <- cfg$effects$color
    if (is.null(color_config)) return(NULL)
    as.character(color_config[[1]])
}

#' Write execution info (snapshot)
write_execution_info <- function(config, run_dir, config_path = NULL, manifest_df = NULL) {
    exec_dir <- file.path(run_dir, "execution_info")
    dir.create(exec_dir, recursive = TRUE, showWarnings = FALSE)

    written_files <- character(0)

    yaml::write_yaml(config, file.path(exec_dir, "config_used.yaml"))
    written_files <- c(written_files, file.path(exec_dir, "config_used.yaml"))

    if (!is.null(config_path)) {
        tryCatch(
            {
                writeLines(normalizePath(config_path, winslash = "/", mustWork = FALSE), file.path(exec_dir, "config_path.txt"))
                written_files <- c(written_files, file.path(exec_dir, "config_path.txt"))
            },
            error = function(e) warning("Could not write config_path.txt")
        )
    }

    # Timestamp
    writeLines(as.character(Sys.time()), file.path(exec_dir, "timestamp.txt"))

    # Session Info
    writeLines(capture.output(sessionInfo()), file.path(exec_dir, "sessionInfo.txt"))

    written_files
}

# Helper for %||%
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}
