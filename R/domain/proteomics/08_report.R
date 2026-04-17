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
        if (!nzchar(src_dir) || !file.exists(file.path(src_dir, "report_template_proteomics.Rmd"))) {
            # When running via targets::tar_source(), try the working directory first
            cwd_path <- file.path("R", "domain", "proteomics")
            if (file.exists(file.path(cwd_path, "report_template_proteomics.Rmd"))) {
                src_dir <- cwd_path
            } else {
                # Fallback to project root from config
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

    # Ensure run_dir exists.  The Rmd and its HTML output live in the *parent*
    # results directory (alongside the pptx / pipeline summary), not inside the
    # proteomics/ sub-folder.
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    parent_dir <- if (basename(run_dir) == "proteomics") dirname(run_dir) else run_dir
    dest_rmd <- file.path(parent_dir, "report_proteomics.Rmd")
    file.copy(template_path, dest_rmd, overwrite = TRUE)

    # Ensure execution_info/config_used.yaml exists (needed by the template).
    # The template looks for it relative to its own location, so keep a copy
    # in the parent results dir as well as the proteomics/ subdir.
    for (edir in unique(c(file.path(parent_dir, "execution_info"),
                          file.path(run_dir, "execution_info")))) {
        config_used <- file.path(edir, "config_used.yaml")
        if (!file.exists(config_used)) {
            dir.create(edir, recursive = TRUE, showWarnings = FALSE)
            if (!is.null(config_file) && file.exists(config_file)) {
                file.copy(config_file, config_used, overwrite = TRUE)
            } else {
                yaml::write_yaml(config, config_used)
            }
        }
    }

    # Render into the parent results directory
    out_html <- file.path(parent_dir, "report_proteomics.html")

    message("Rendering proteomics report to: ", out_html)

    render_err      <- NULL
    render_err_full <- NULL
    tryCatch({
        rmarkdown::render(
            input       = dest_rmd,
            output_file = out_html,
            output_dir  = parent_dir,
            quiet       = FALSE,  # surface inner knitr errors to the pipeline log
            envir       = new.env(parent = globalenv())
        )
    }, error = function(e) {
        render_err      <<- e$message
        render_err_full <<- paste(c(e$message,
                                     capture.output(traceback(max.lines = 40))),
                                   collapse = "\n")
        # Print to the live pipeline log so the user can see WHY it failed.
        message("\n==== Proteomics Rmd render ERROR ====")
        message(e$message)
        message("====================================\n")
    })

    # Fallback: if the render did not produce an output file (a setup chunk or
    # the pandoc step crashed), write a minimal HTML stub so downstream targets
    # and the user-facing run flow still see *something* rather than failing.
    if (!file.exists(out_html)) {
        warning("Proteomics report rendering did not produce output file.",
                if (!is.null(render_err)) paste0(" Cause: ", render_err) else "")

        # Drop the full error + traceback next to the stub so the user can send it.
        err_log <- file.path(parent_dir, "report_proteomics_error.log")
        writeLines(render_err_full %||% (render_err %||% "(error not captured)"),
                   err_log)
        message("Full render error + traceback saved to: ", err_log)

        stub <- paste0(
            "<!doctype html><html><head><meta charset='utf-8'>",
            "<title>Proteomics Report (partial)</title>",
            "<style>body{font-family:system-ui,sans-serif;max-width:780px;margin:40px auto;padding:0 20px;color:#222}",
            "h1{color:#b91c1c}code{background:#f3f4f6;padding:2px 6px;border-radius:3px;font-size:0.95em}",
            "pre{background:#0d1220;color:#e7eefc;padding:14px;border-radius:6px;overflow-x:auto;font-size:13px}",
            ".box{background:#fef3c7;border-left:4px solid #d97706;padding:14px 18px;border-radius:4px;margin:18px 0}",
            "</style></head><body>",
            "<h1>Proteomics report did not render</h1>",
            "<div class='box'>The pipeline finished, but the HTML report template failed to render. ",
            "All upstream analysis outputs (CSV tables, plots, PPTX, payloads) are still available under ",
            "the results directory.</div>",
            "<p><strong>Render error:</strong></p>",
            "<pre>", htmltools::htmlEscape(render_err %||% "(not captured)"), "</pre>",
            "<p>Full traceback saved to <code>report_proteomics_error.log</code> in this folder.</p>",
            "</body></html>"
        )
        writeLines(stub, out_html)
        message("Wrote fallback stub report: ", out_html)
    } else {
        message("Proteomics report rendered successfully: ", out_html)
    }

    out_html
}
