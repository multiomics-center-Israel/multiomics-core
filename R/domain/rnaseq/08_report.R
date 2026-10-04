#' RNA-seq Report Rendering
#'
#' Copies the canonical Rmd template into the run directory and renders
#' a self-contained HTML report.

#' Render the RNA-seq analysis report
#'
#' @param run_dir  The results run directory (e.g. outputs/project/Results_.../rna)
#' @param config   Full pipeline config list
#' @param config_file Path to the original YAML config file. Retained for the
#'   pipeline call sites; the config snapshot is serialized from \code{config}
#'   rather than copied from this path, so that it matches what the
#'   \code{execution_info_files} target writes.
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

    # Write execution_info/config_used.yaml (needed by the template), keeping a
    # copy in both parent and rna/ subdirectories. Rewritten on every render
    # rather than only when absent: the rna/ copy is written by nothing else, so
    # the old guard left the first run's snapshot in place forever and config
    # edits appeared to do nothing.
    #
    # Always serialized from `config`, never copied from the source YAML. The
    # parent path is the same file the execution_info_files target owns, and that
    # target writes it with this same yaml::write_yaml(config, ...). Copying the
    # raw YAML here would put different content at a tracked path -- load_config()
    # adds .config_path and .config_mtime, which the source file does not carry --
    # leaving the target outdated the moment the report finished.
    for (edir in unique(c(file.path(parent_dir, "execution_info"),
                          file.path(run_dir, "execution_info")))) {
        dir.create(edir, recursive = TRUE, showWarnings = FALSE)
        yaml::write_yaml(config, file.path(edir, "config_used.yaml"))
    }

    # Render into the parent results directory
    out_html <- file.path(parent_dir, "report_rnaseq.html")

    # Ensure any stale report from a previous run is gone, as the proteomics and
    # multiomics renderers do. The config snapshot above is now rewritten every
    # time, so a render that fails afterwards would otherwise leave last run's
    # HTML sitting beside this run's snapshot, reading as though that report had
    # been produced from these settings. It also makes the file.exists(out_html)
    # check below meaningful rather than satisfied by the leftover.
    if (file.exists(out_html)) try(file.remove(out_html), silent = TRUE)

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
