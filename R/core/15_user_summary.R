#' User-Provided Project Summary
#'
#' Generates a project overview section for the HTML report from user notes.
#' If Claude CLI is available and a technical report file is provided,
#' blends user notes with extracted technical details into a polished summary.
#' Otherwise, returns user notes as-is wrapped in HTML.

#' Generate project summary HTML from user notes and optional tech report
#'
#' @param config Full pipeline config (reads project$user_notes,
#'               project$technical_report_file)
#' @return HTML string for the Project Overview section, or NULL if no notes
generate_project_summary <- function(config) {

  user_notes <- config$project$user_notes
  tech_report_file <- config$project$technical_report_file
  has_notes <- !is.null(user_notes) && nzchar(trimws(user_notes %||% ""))
  has_tech_report <- !is.null(tech_report_file) && nzchar(tech_report_file) &&
                     file.exists(tech_report_file)

  # Nothing to do if neither notes nor tech report

  if (!has_notes && !has_tech_report) {
    return(NULL)
  }

  user_notes <- if (has_notes) trimws(user_notes) else ""

  # Check if Claude CLI is available
  claude_available <- nzchar(Sys.which("claude"))

  # If no tech report file, return notes as-is
  if (!has_tech_report) {
    message("[user_summary] No tech report file — using notes as-is")
    return(wrap_notes_html(user_notes))
  }

  # If no user notes but tech report exists, read and include raw text
  if (!has_notes && !claude_available) {
    message("[user_summary] No user notes, no Claude — including raw tech report text")
    doc_text <- tryCatch({
      ext <- tolower(tools::file_ext(tech_report_file))
      if (ext %in% c("log", "txt", "tsv", "csv")) {
        paste(readLines(tech_report_file, warn = FALSE), collapse = "\n")
      } else {
        paste(system2("pandoc", c("-t", "plain", shQuote(tech_report_file)),
                       stdout = TRUE, stderr = FALSE), collapse = "\n")
      }
    }, error = function(e) NULL)
    if (!is.null(doc_text) && nchar(doc_text) > 20) {
      # Truncate to first 5000 chars if very large
      if (nchar(doc_text) > 5000) doc_text <- paste0(substr(doc_text, 1, 5000), "\n\n[... truncated]")
      return(wrap_notes_html(paste0(
        "<strong>Technical Report (raw):</strong>\n<pre style='white-space: pre-wrap; font-size: 0.85em; max-height: 400px; overflow-y: auto; background: #f8f9fa; padding: 12px; border-radius: 4px;'>",
        htmltools::htmlEscape(doc_text),
        "</pre>"
      )))
    }
    return(NULL)
  }

  # If no Claude but both notes and tech report, show both
  if (!claude_available) {
    message("[user_summary] No Claude CLI — using notes + raw tech report")
    doc_text <- tryCatch({
      ext <- tolower(tools::file_ext(tech_report_file))
      if (ext %in% c("log", "txt", "tsv", "csv")) {
        paste(readLines(tech_report_file, warn = FALSE), collapse = "\n")
      } else {
        paste(system2("pandoc", c("-t", "plain", shQuote(tech_report_file)),
                       stdout = TRUE, stderr = FALSE), collapse = "\n")
      }
    }, error = function(e) NULL)
    notes_html <- user_notes
    if (!is.null(doc_text) && nchar(doc_text) > 20) {
      if (nchar(doc_text) > 5000) doc_text <- paste0(substr(doc_text, 1, 5000), "\n\n[... truncated]")
      notes_html <- paste0(notes_html,
        "\n\n<details><summary><strong>Technical Report Details</strong></summary>\n",
        "<pre style='white-space: pre-wrap; font-size: 0.85em; max-height: 400px; overflow-y: auto; background: #f8f9fa; padding: 12px; border-radius: 4px;'>",
        htmltools::htmlEscape(doc_text),
        "</pre></details>")
    }
    return(wrap_notes_html(notes_html))
  }

  # Validate tech report file exists
  if (!file.exists(tech_report_file)) {
    message(sprintf("[user_summary] Tech report file not found: %s — using notes as-is",
                    tech_report_file))
    return(wrap_notes_html(user_notes))
  }

  # Extract text from tech report
  message("[user_summary] Extracting text from technical report...")
  doc_text <- tryCatch({
    ext <- tolower(tools::file_ext(tech_report_file))
    if (ext %in% c("log", "txt", "tsv", "csv")) {
      paste(readLines(tech_report_file, warn = FALSE), collapse = "\n")
    } else {
      paste(system2("pandoc", c("-t", "plain", shQuote(tech_report_file)),
                    stdout = TRUE, stderr = FALSE), collapse = "\n")
    }
  }, error = function(e) NULL)

  if (is.null(doc_text) || nchar(doc_text) < 20) {
    message("[user_summary] Could not extract text from tech report — using notes as-is")
    return(wrap_notes_html(user_notes))
  }

  # Blend with Claude
  message(sprintf("[user_summary] Blending notes + tech report (%d chars) via Claude...",
                  nchar(doc_text)))

  prompt <- sprintf(paste0(
    "You are writing a 'Project Overview' section for a scientific analysis report. ",
    "Blend the user's project notes with the technical details from the facility report ",
    "into a cohesive, professional summary. Write 2-4 paragraphs in clear scientific prose. ",
    "Preserve ALL user-provided details and hypotheses. ",
    "Integrate relevant technical details (instrument, method, parameters) naturally. ",
    "Do NOT use markdown headers or bullet points — write flowing paragraphs. ",
    "Return ONLY the summary text, no preamble.\n\n",
    "USER NOTES:\n%s\n\n",
    "TECHNICAL REPORT:\n%s"
  ), user_notes, substr(doc_text, 1, 6000))

  blended <- tryCatch({
    cmd <- sprintf("unset CLAUDECODE; claude --print --model sonnet --no-session-persistence %s",
                    shQuote(prompt))
    raw <- system(cmd, intern = TRUE, timeout = 120)
    paste(raw, collapse = "\n")
  }, error = function(e) {
    message(sprintf("[user_summary] Claude blending failed: %s", e$message))
    NULL
  })

  if (!is.null(blended) && nchar(trimws(blended)) > 50) {
    message("[user_summary] AI-blended summary generated successfully")
    return(wrap_notes_html(blended))
  }

  # Fallback: notes as-is
  message("[user_summary] Claude output insufficient — using notes as-is")
  wrap_notes_html(user_notes)
}


#' Wrap plain text notes into styled HTML
#'
#' @param text Plain text to wrap
#' @return HTML string
wrap_notes_html <- function(text) {
  # Convert newlines to <br> for display
  html_text <- gsub("\n", "<br>\n", htmltools::htmlEscape(text))
  sprintf(
    '<div class="project-overview-box" style="background:#f8f9fa; border-left:4px solid #4a90d9; padding:16px 20px; border-radius:6px; margin-bottom:20px; font-size:14px; line-height:1.7; color:#333;">%s</div>',
    html_text
  )
}
