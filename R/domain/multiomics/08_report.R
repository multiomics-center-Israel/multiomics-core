#' Multi-Omics Report Rendering
#'
#' Copies the canonical Rmd template into the run directory and renders
#' a self-contained HTML report.

#' Render the multi-omics integration report
#'
#' @param run_dir  The results run directory (e.g. outputs/project/Results_.../multiomics)
#' @param config   Full pipeline config list
#' @param config_file Path to the original YAML config file. Retained for the
#'   pipeline call sites; the config snapshot is serialized from \code{config}
#'   rather than copied from this path, so that it matches what the
#'   \code{execution_info_files} target writes.
#' @return Path to the rendered HTML file (character, format = "file")
#' @export
render_multiomics_report <- function(run_dir, config, config_file = NULL) {

    # Locate the canonical Rmd template shipped with multiomics-core
    template_path <- system.file("report_template_multiomics.Rmd",
                                  package = "multiomics.core",
                                  mustWork = FALSE)

    # Fallback: template lives alongside this source file
    if (!nzchar(template_path) || !file.exists(template_path)) {
        src_dir <- system.file("R", "domain", "multiomics",
                                package = "multiomics.core",
                                mustWork = FALSE)
        if (!nzchar(src_dir)) {
            # When running via targets::tar_source(), try the working directory first,
            # then fall back to the project root from config
            src_dir <- file.path("R", "domain", "multiomics")
            if (!file.exists(file.path(src_dir, "report_template_multiomics.Rmd"))) {
                proj_dir <- config$project$dir %||% "."
                src_dir <- file.path(proj_dir, "R", "domain", "multiomics")
            }
        }
        template_path <- file.path(src_dir, "report_template_multiomics.Rmd")
    }

    if (!file.exists(template_path)) {
        warning("Multi-omics report template not found at: ", template_path,
                ". Skipping report generation.")
        return(NA_character_)
    }

    # Ensure run_dir exists, then copy template into it (the knit directory)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    dest_rmd <- file.path(run_dir, "report_multiomics.Rmd")
    file.copy(template_path, dest_rmd, overwrite = TRUE)

    # Write execution_info/config_used.yaml (needed by the template). Rewritten on
    # every render rather than only when absent: the old guard left the first run's
    # snapshot in place forever, so config edits appeared to do nothing.
    #
    # Always serialized from `config`, never copied from the source YAML. This can
    # be the same file the execution_info_files target owns, and that target writes
    # it with this same yaml::write_yaml(config, ...). Copying the raw YAML here
    # would put different content at a tracked path -- load_config() adds
    # .config_path and .config_mtime, which the source file does not carry --
    # leaving the target outdated the moment the report finished.
    exec_dir <- file.path(run_dir, "execution_info")
    dir.create(exec_dir, recursive = TRUE, showWarnings = FALSE)
    yaml::write_yaml(config, file.path(exec_dir, "config_used.yaml"))

    # Render
    out_html <- file.path(run_dir, "report_multiomics.html")

    # Ensure any stale report from a previous run is gone, as the proteomics
    # renderer already does. The config snapshot above is now rewritten every
    # time, so a render that fails afterwards would otherwise leave last run's
    # HTML sitting beside this run's snapshot, reading as though that report had
    # been produced from these settings.
    if (file.exists(out_html)) try(file.remove(out_html), silent = TRUE)

    message("Rendering multi-omics report to: ", out_html)

    rmarkdown::render(
        input       = dest_rmd,
        output_file = out_html,
        output_dir  = run_dir,
        quiet       = TRUE,
        envir       = new.env(parent = globalenv())
    )

    if (file.exists(out_html)) {
        message("Multi-omics report rendered successfully: ", out_html)
    } else {
        warning("Multi-omics report rendering did not produce output file.")
    }

    out_html
}
