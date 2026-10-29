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
#' @param pre_info The \code{info} list of \code{preprocess_rna()}'s result,
#'   passed to the template as \code{params$rna_info} so the report describes
#'   the gene filter that ran rather than the one the config asked for. NULL
#'   renders the filtering text as not recorded.
#' @return Path to the HTML report (character, format = "file"). If rendering
#'   fails, that file is a stub naming the failing chunk and error, the same
#'   details go to \code{report_rnaseq_error.log}, and a warning is raised.
#' @export
render_rnaseq_report <- function(run_dir, config, config_file = NULL, pre_info = NULL) {

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
    # Likewise a previous failure's log, or a successful re-run would sit beside
    # diagnostics saying it failed.
    err_log <- file.path(parent_dir, "report_rnaseq_error.log")
    if (file.exists(err_log)) try(file.remove(err_log), silent = TRUE)

    message("Rendering RNA-seq report to: ", out_html)

    # knitr's own "Quitting from lines ..." note is lost with quiet = TRUE, so the
    # failing chunk is recorded here while knitr is still inside it. By the time
    # tryCatch's handler runs, the stack has unwound and the label is gone.
    render_err   <- NULL
    render_chunk <- NULL
    tryCatch(
        withCallingHandlers(
            rmarkdown::render(
                input       = dest_rmd,
                output_file = out_html,
                output_dir  = parent_dir,
                quiet       = TRUE,
                params      = list(rna_info = pre_info),
                envir       = new.env(parent = globalenv())
            ),
            error = function(e) {
                if (isTRUE(getOption("knitr.in.progress"))) {
                    render_chunk <<- knitr::opts_current$get("label")
                }
            }
        ),
        error = function(e) render_err <<- conditionMessage(e)
    )

    # rmarkdown can fail after pandoc has already written the HTML, so a file on
    # disk is not proof of success on its own.
    if (is.null(render_err) && file.exists(out_html)) {
        message("Report rendered successfully: ", out_html)
        return(out_html)
    }

    # The report failed, but every analysis output it would have shown is already
    # on disk. Leave the cause where the user will look for the report, plus a
    # stub so the file target still resolves and the run is not marked failed.
    where <- if (!is.null(render_chunk)) {
        lines <- find_rmd_chunk_lines(dest_rmd, render_chunk)
        if (is.null(lines)) sprintf("chunk '%s'", render_chunk)
        else sprintf("chunk '%s' (lines %d-%d)", render_chunk, lines[1], lines[2])
    } else {
        "outside any chunk (setup or pandoc step)"
    }
    cause <- render_err %||% "no error was raised, but no HTML was written"
    writeLines(c(paste("Template:", dest_rmd),
                 paste("Failed in:", where),
                 paste("Error:", cause)),
               err_log)

    writeLines(paste0(
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>RNA-seq Report (not rendered)</title>",
        "<style>body{font-family:system-ui,sans-serif;max-width:780px;margin:40px auto;padding:0 20px;color:#222}",
        "h1{color:#b91c1c}pre{background:#f3f4f6;padding:14px;border-radius:6px;white-space:pre-wrap}</style>",
        "</head><body><h1>RNA-seq report did not render</h1>",
        "<p>The analysis finished; its tables and plots are in this results folder. ",
        "Only the HTML report failed.</p>",
        "<p><strong>Failed in:</strong> ", htmltools::htmlEscape(where), "</p>",
        "<pre>", htmltools::htmlEscape(cause), "</pre>",
        "<p>The same details are in <code>report_rnaseq_error.log</code>.</p>",
        "</body></html>"),
        out_html)
    message("Wrote fallback stub report: ", out_html)

    # Raised last: under options(warn = 2) or an exiting warning handler it ends
    # the call, and the stub must already be on disk by then.
    warning("RNA-seq report rendering failed in ", where, ": ", cause,
            "\nTemplate: ", dest_rmd,
            "\nDetails: ", err_log)

    out_html
}

#' Locate a knitr chunk in an Rmd file
#'
#' Turns the chunk label knitr reports into the line range a person can open,
#' since a rendering error from a 3,000-line template is not much use without it.
#'
#' @param rmd Path to the Rmd file.
#' @param label Chunk label, as \code{knitr::opts_current$get("label")} gives it.
#' @return Integer vector \code{c(start, end)}: the lines of the chunk's opening
#'   and closing fences. \code{NULL} if no chunk carries that label, which is
#'   the case for unlabelled chunks (knitr names them \code{unnamed-chunk-N}).
#' @examples
#' rmd <- tempfile(fileext = ".Rmd")
#' writeLines(c("# Title", "```{r setup, include=FALSE}", "x <- 1", "```"), rmd)
#' find_rmd_chunk_lines(rmd, "setup")  # c(2, 4)
find_rmd_chunk_lines <- function(rmd, label) {
    src <- readLines(rmd, warn = FALSE)
    headers <- grep("^\\s*```\\{r[ ,]", src)
    labels <- sub("^\\s*```\\{r[ ,]+([^ ,}]*).*$", "\\1", src[headers])
    start <- headers[labels == label][1]
    if (is.na(start)) return(NULL)
    fences <- grep("^\\s*```\\s*$", src)
    end <- fences[fences > start][1]
    if (is.na(end)) end <- length(src)
    c(start, end)
}
