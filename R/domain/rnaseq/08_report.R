#' RNA-seq Report Rendering
#'
#' Copies the canonical Rmd template into the run directory and renders
#' a self-contained HTML report.

#' Render the RNA-seq analysis report
#'
#' @param run_dir  The results run directory (e.g. outputs/project/Results_.../rna)
#' @param config   Full pipeline config list
#' @param config_file Path to the original YAML config file (for embedding)
#' @return Path to the rendered HTML file (character, format = "file")
#' @export
render_rnaseq_report <- function(run_dir, config, config_file = NULL) {

    # Locate the canonical Rmd template shipped with multiomics-core.
    template_path <- system.file("report_template.Rmd",
                                  package = "multiomics.core",
                                  mustWork = FALSE)

    if (!nzchar(template_path) || !file.exists(template_path)) {
        template_path <- file.path(getwd(), "R", "domain", "rnaseq",
                                   "report_template.Rmd")
    }

    if (!file.exists(template_path)) {
        warning("Report template not found at: ", template_path,
                ". Skipping report generation.")
        return(NA_character_)
    }

    # The Rmd and its HTML output live in the *parent* results directory
    # (alongside the pptx / pipeline summary), not inside the rna/ sub-folder.
    # The template already references paths as run_dir/rna/... internally.
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    parent_dir <- if (basename(run_dir) %in% c("rna", "rnaseq")) dirname(run_dir) else run_dir
    dest_rmd <- file.path(parent_dir, "report_rnaseq.Rmd")
    file.copy(template_path, dest_rmd, overwrite = TRUE)

    # Ensure execution_info/config_used.yaml exists (needed by the template).
    # Keep a copy in both parent and rna/ subdirectories.
    for (edir in unique(c(file.path(parent_dir, "execution_info"),
                          file.path(run_dir, "execution_info")))) {
        # Rewritten on every render, not only when absent: a stale snapshot from
        # an earlier run silently overrode any config change, so report settings
        # edited in the config appeared to do nothing.
        config_used <- file.path(edir, "config_used.yaml")
        {
            dir.create(edir, recursive = TRUE, showWarnings = FALSE)
            if (!is.null(config_file) && file.exists(config_file)) {
                file.copy(config_file, config_used, overwrite = TRUE)
            } else {
                yaml::write_yaml(config, config_used)
            }
        }
    }

    # Render into the parent results directory
    out_html <- file.path(parent_dir, "report_rnaseq.html")

    message("Rendering RNA-seq report to: ", out_html)

    tryCatch({
        rmarkdown::render(
            input       = dest_rmd,
            output_file = out_html,
            output_dir  = parent_dir,
            quiet       = TRUE,
            envir       = new.env(parent = globalenv())
        )
    }, error = function(e) {
        warning("RNA-seq report rendering failed: ", e$message,
                "\nTemplate: ", dest_rmd,
                "\nCheck the Rmd for errors.")
    })

    rmarkdown::render(
        input       = dest_rmd,
        output_file = out_html,
        output_dir  = run_dir,
        quiet       = TRUE,
        envir       = new.env(parent = globalenv())
    )

    if (file.exists(out_html)) {
        message("Report rendered successfully: ", out_html)
    } else {
        warning("Report rendering did not produce output file.")
    }

    out_html
}
