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

write_project_rdata <- function(run_dir, ..., filename = "project.Rdata") {
    exec_dir <- file.path(run_dir, "execution_info")
    dir.create(exec_dir, recursive = TRUE, showWarnings = FALSE)
    out <- file.path(exec_dir, filename)

    objects_list <- list(...)
    save(list = names(objects_list), envir = list2env(objects_list), file = out)

    if (!file.exists(out)) stop("save() did not create project.Rdata.")
    out
}

# Helper for %||%
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}
