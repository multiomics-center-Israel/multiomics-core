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

    config <- apply_renamed_config_keys(config)
    warn_retired_config_keys(config)

    if (!is.null(config$modes$proteomics)) validate_proteomics_config(config$modes$proteomics)
    if (!is.null(config$modes$rna)) validate_rna_config(config$modes$rna)
    if (!is.null(config$modes$metabolomics)) validate_metabolomics_config(config$modes$metabolomics)
    if (!is.null(config$modes$multiomics)) {
        # validate_multiomics_config() applies defaults and returns the modified
        # section — capture the return value so defaults actually persist.
        config$modes$multiomics <- validate_multiomics_config(config$modes$multiomics)
    }
    if (!is.null(config$modes$de_integration)) {
        config$modes$de_integration <- validate_de_integration_config(config$modes$de_integration)
    }

    invisible(config)
}

#' Config keys the pipeline no longer reads
#'
#' One row per retired key: its dotted path from the config root, why it is not
#' read, and since when (\code{"never read"} for keys that appear in old configs
#' or templates but have no reader in the history of R/). Checked by
#' validate_config() through warn_retired_config_keys(). Add a row here when a
#' change stops reading a key, instead of leaving it to be silently ignored.
#' @noRd
.retired_config_keys <- data.frame(
    key = c(
        "modes.proteomics.de_table.pass_any_col",
        "modes.rna.de_table.pass_any_col",
        "modes.proteomics.qc_post.outputs.write_de_tables",
        "modes.rna.qc_post.outputs.write_de_tables",
        "modes.proteomics.pathway.gene_set_note",
        "modes.rna.pathway.gene_set_note",
        "modes.proteomics.report.show_ma_imputation",
        "modes.proteomics.report.show_sample_scatter",
        "modes.proteomics.report.show_volcano_nonimputed",
        "modes.proteomics.report.project_id",
        "modes.proteomics.technical_report.references",
        "modes.multiomics.report.sections.diablo_interactive_loadings",
        "modes.multiomics.report.sections.rna_protein_concordance_hist",
        "modes.metabolomics.files.feature_annotations",
        "modes.rna.clustering.steps.partition.x_axis_col",
        "modes.proteomics.clustering.steps.partition.x_axis_col",
        "modes.metabolomics.clustering.steps.partition.x_axis_col",
        "modes.lipidomics.commentary.retry_delay",
        "design.batch_column",
        "output.output_dir",
        "output.report_format"
    ),
    reason = c(
        "the DE pass column is always pass_any_contrast",
        "the RNA pipeline uses the pass_any_contrast column and never read this key",
        "the DE summary tables are always written",
        "the RNA pipeline never read this key",
        "no code reads it",
        "no code reads it",
        "no report section is controlled by it",
        "no report section is controlled by it",
        "no report section is controlled by it",
        "the report does not use it",
        paste("the technical summary uses facility, sample_prep, ms_acquisition,",
              "search_engine and search_parameters only"),
        "the multiomics report has no section of that name",
        "the multiomics report has no section of that name",
        "feature annotations are taken from the columns of the data table",
        "renamed to clustering.steps.partition.profile_x_axis_col",
        "renamed to clustering.steps.partition.profile_x_axis_col",
        "renamed to clustering.steps.partition.profile_x_axis_col",
        "the commentary scripts use a fixed retry delay",
        "no code reads it",
        "outputs are written under paths.out",
        "no code reads it; reports are rendered by their own templates"
    ),
    since = c(
        "22c7cb8 (2026-07-10, #109)",
        "never read",
        "22c7cb8 (2026-07-10, #109)",
        "never read",
        "never read",
        "never read",
        "never read",
        "never read",
        "never read",
        "never read",
        "never read",
        "never read",
        "never read",
        "never read",
        "f56cfa2 (2026-06-29, #102)",
        "f56cfa2 (2026-06-29, #102)",
        "f56cfa2 (2026-06-29, #102)",
        "never read",
        "never read",
        "never read",
        "never read"
    ),
    stringsAsFactors = FALSE
)

#' Move values from renamed config keys to their current names
#'
#' Older wizard and run.R configs wrote \code{modes.rna.id_columns.gene_id_col};
#' the RNA code reads \code{id_columns.gene_id}. Until now the old name worked
#' only through R's partial matching of \code{$}. Its value is copied to the
#' current name (unless that is already set) and the old key is dropped, with
#' one warning.
#'
#' @param config Full pipeline config.
#' @return The config with renamed keys moved.
apply_renamed_config_keys <- function(config) {
    ids <- config$modes$rna$id_columns
    if (is.list(ids) && "gene_id_col" %in% names(ids)) {
        if (is.null(ids[["gene_id", exact = TRUE]])) {
            ids[["gene_id"]] <- ids[["gene_id_col", exact = TRUE]]
        }
        ids[["gene_id_col"]] <- NULL
        config$modes$rna$id_columns <- ids
        warning("modes.rna.id_columns.gene_id_col is the old name of ",
                "modes.rna.id_columns.gene_id; its value was used. Rename it.",
                call. = FALSE)
    }
    config
}

#' Find retired config keys present in a config
#'
#' A key counts as present when its name exists at that path, even with a null
#' value, since a \code{key: null} line in a YAML file is still a stale setting.
#'
#' @param config Full pipeline config (list, as returned by load_config()).
#' @param retired Data frame of retired keys (\code{key}, \code{reason},
#'   \code{since}); defaults to the package list.
#' @return The rows of \code{retired} whose key is present in \code{config}.
find_retired_config_keys <- function(config, retired = .retired_config_keys) {
    has_path <- function(x, parts) {
        for (p in parts) {
            if (!is.list(x) || !p %in% names(x)) return(FALSE)
            x <- x[[p]]
        }
        TRUE
    }
    present <- vapply(strsplit(retired$key, ".", fixed = TRUE),
                      function(parts) has_path(config, parts), logical(1))
    retired[present, , drop = FALSE]
}

#' Warn once about retired config keys present in a config
#'
#' Never an error: a retired key changes nothing, it is only misleading to
#' whoever reads the config.
#'
#' @param config Full pipeline config.
#' @return The present retired keys (data frame), invisibly.
warn_retired_config_keys <- function(config) {
    found <- find_retired_config_keys(config)
    if (nrow(found) > 0) {
        warning("The config contains keys that are no longer read:\n",
                paste0("  - ", found$key, " is no longer read (", found$reason,
                       "); remove it", collapse = "\n"),
                call. = FALSE)
    }
    invisible(found)
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
