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

  # Nothing to do if no user notes
  if (is.null(user_notes) || !nzchar(trimws(user_notes))) {
    return(NULL)
  }

  user_notes <- trimws(user_notes)

  # Check if Claude CLI is available
  claude_available <- nzchar(Sys.which("claude"))

  # If no Claude or no tech report file, return notes as-is
  if (!claude_available || is.null(tech_report_file) || !nzchar(tech_report_file)) {
    message("[user_summary] No Claude CLI or no tech report file — using notes as-is")
    return(wrap_notes_html(user_notes))
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
