#' Proteomics Report Rendering
#'
#' Copies the canonical Rmd template into the run directory and renders
#' a self-contained HTML report.

#' Render the proteomics analysis report
#'
#' @param run_dir  The results run directory (e.g. outputs/project/Results_...)
#' @param config   Full pipeline config list
#' @param config_file Path to the original YAML config file (for embedding)
#' @return Path to the rendered HTML file (character, format = "file")
#' @export
render_proteomics_report <- function(run_dir, config, config_file = NULL) {

    # Locate the canonical Rmd template shipped with multiomics-core
    template_path <- system.file("report_template_proteomics.Rmd",
                                  package = "multiomics.core",
                                  mustWork = FALSE)

    # Fallback: template lives alongside this source file
    if (!nzchar(template_path) || !file.exists(template_path)) {
        src_dir <- system.file("R", "domain", "proteomics",
                                package = "multiomics.core",
                                mustWork = FALSE)
        if (!nzchar(src_dir)) {
            # When running via targets::tar_source(), try the working directory first,
            # then fall back to the project root from config
            src_dir <- file.path("R", "domain", "proteomics")
            if (!file.exists(file.path(src_dir, "report_template_proteomics.Rmd"))) {
                proj_dir <- config$project$dir %||% "."
                src_dir <- file.path(proj_dir, "R", "domain", "proteomics")
            }
        }
        template_path <- file.path(src_dir, "report_template_proteomics.Rmd")
    }

    if (!file.exists(template_path)) {
        warning("Proteomics report template not found at: ", template_path,
                ". Skipping report generation.")
        return(NA_character_)
    }

    # Ensure run_dir exists, then copy template into it (the knit directory)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    dest_rmd <- file.path(run_dir, "report_proteomics.Rmd")
    file.copy(template_path, dest_rmd, overwrite = TRUE)

    # Ensure execution_info/config_used.yaml exists (needed by the template)
    exec_dir <- file.path(run_dir, "execution_info")
    config_used <- file.path(exec_dir, "config_used.yaml")
    if (!file.exists(config_used)) {
        dir.create(exec_dir, recursive = TRUE, showWarnings = FALSE)
        if (!is.null(config_file) && file.exists(config_file)) {
            file.copy(config_file, config_used, overwrite = TRUE)
        } else {
            yaml::write_yaml(config, config_used)
        }
    }

    # Render
    out_html <- file.path(run_dir, "report_proteomics.html")

    message("Rendering proteomics report to: ", out_html)

    rmarkdown::render(
        input       = dest_rmd,
        output_file = out_html,
        output_dir  = run_dir,
        quiet       = TRUE,
        envir       = new.env(parent = globalenv())
    )

    if (file.exists(out_html)) {
        message("Proteomics report rendered successfully: ", out_html)
    } else {
        warning("Proteomics report rendering did not produce output file.")
    }

    out_html
}
