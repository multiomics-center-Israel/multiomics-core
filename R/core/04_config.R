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
    if (!is.null(config$modes$multiomics)) {
        # validate_multiomics_config() applies defaults and returns the modified
        # section — capture the return value so defaults actually persist.
        config$modes$multiomics <- validate_multiomics_config(config$modes$multiomics)
    }

    invisible(config)
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
    if (is.null(color_config)) {
        return(NULL)
    }
    as.character(color_config[[1]])
}

#' Resolve the canonical group column from a mode config.
#'
#' Walks the standard rnaseq fallback chain to determine which metadata column
#' identifies the experimental group used for filtering, DE summaries,
#' commentary, and design-viability checks. Strict at every level: NULL, NA,
#' empty strings, and non-character values all fall through to the next
#' fallback. Only a length-1, non-NA, non-empty character scalar is accepted.
#'
#' Resolution order:
#' \enumerate{
#'   \item \code{cfg$filtering$group_col}
#'   \item \code{cfg$de_table$group_col}
#'   \item First element of \code{cfg$effects$color}
#'   \item Literal default \code{"Condition"}
#' }
#'
#' @param cfg Mode-level config list (e.g. \code{config$modes$rna}).
#' @return Single non-empty character string. Never \code{NULL}.
#'
#' @note This function is the canonical resolver going forward, but as of
#'   introduction it is NOT yet adopted at every call site. Several rnaseq,
#'   proteomics, and report-template files still inline divergent variants:
#'   the proteomics chain skips \code{filtering$group_col}; some sites use
#'   bare \code{effects$color} (vector) instead of the first element; some
#'   omit the terminal default. Migration is tracked in follow-up issues
#'   filed alongside this change. Do not assume every group-column lookup
#'   in the repo flows through this function.
#'
#' @examples
#' cfg <- list(filtering = list(group_col = "treatment"),
#'             effects   = list(color = c("batch", "treatment")))
#' resolve_group_col(cfg)  # "treatment"
#'
#' cfg <- list(effects = list(color = c("Group", "batch")))
#' resolve_group_col(cfg)  # "Group" (first element only)
#'
#' resolve_group_col(list())  # "Condition" (terminal default)
resolve_group_col <- function(cfg) {
    is_valid <- function(x) {
        is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
    }

    candidates <- list(
        cfg$filtering$group_col,
        cfg$de_table$group_col,
        if (!is.null(cfg$effects$color)) as.character(cfg$effects$color)[1] else NULL
    )

    for (cand in candidates) {
        if (is_valid(cand)) return(cand)
    }
    "Condition"
}

#' Write execution info (snapshot)
write_execution_info <- function(config, run_dir, config_path = NULL, manifest_df = NULL, targets_file = NULL) {
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

    # _targets.R snapshot (copy actual file)
    if (!is.null(targets_file)) {
        tryCatch(
            {
                dest <- file.path(exec_dir, "_targets_used.R")
                file.copy(targets_file, dest, overwrite = TRUE)
                written_files <- c(written_files, dest)
            },
            error = function(e) warning("Could not copy _targets_used.R")
        )
    }

    # Pipeline definition file snapshots (copy actual files)
    tryCatch(
        {
            pipeline_files <- sort(list.files(
                "R/pipeline",
                pattern = "\\.R$",
                full.names = TRUE,
                recursive = TRUE
            ))
            if (length(pipeline_files) > 0) {
                pipe_dest_dir <- file.path(exec_dir, "pipeline_files_used")
                dir.create(pipe_dest_dir, recursive = TRUE, showWarnings = FALSE)
                for (pf in pipeline_files) {
                    file.copy(pf, file.path(pipe_dest_dir, basename(pf)), overwrite = TRUE)
                }
                written_files <- c(written_files, file.path(pipe_dest_dir, basename(pipeline_files)))
            }
        },
        error = function(e) warning("Could not copy pipeline files")
    )

    # Git commit hash
    tryCatch(
        {
            commit_hash <- system("git rev-parse HEAD", intern = TRUE)
            if (length(commit_hash) > 0 && nchar(commit_hash[1]) > 0) {
                writeLines(commit_hash[1], file.path(exec_dir, "git_commit.txt"))
                written_files <- c(written_files, file.path(exec_dir, "git_commit.txt"))
            }
        },
        error = function(e) warning("Could not write git_commit.txt")
    )

    written_files
}

# Helper for %||%
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}
