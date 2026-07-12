# R/core/18_renviron.R
#
# Tiny helper to record a single environment variable in an .Renviron-style
# file, used by `make setup` to pin MUMMICHOG_PYTHON per machine. Kept generic
# and base-R only; it is setup tooling, never called during a pipeline run.


#' Upsert a single VAR=value line in an .Renviron-style file
#'
#' Sets `key` to `value` in the env file at `path`, creating the file if it does
#' not exist. An existing assignment for `key` (optionally written as
#' `export KEY=...` or with surrounding whitespace) is replaced in place and any
#' duplicate assignments are dropped; if `key` is absent the assignment is
#' appended. Every other line is preserved verbatim, so an existing file — which
#' may hold unrelated secrets — is never clobbered, and the file's contents are
#' never printed.
#'
#' Idempotent: calling it again with the same arguments leaves the file in the
#' same state (exactly one `key=value` line, no duplicates).
#'
#' @param key   Environment variable name, e.g. "MUMMICHOG_PYTHON".
#' @param value Value to assign, written verbatim. Quote it yourself if it may
#'              contain spaces or shell metacharacters.
#' @param path  Path to the env file. Defaults to ".Renviron" in the working dir.
#' @return Invisibly one of "added", "updated" or "unchanged".
#' @examples
#' \dontrun{
#' set_renviron_var("MUMMICHOG_PYTHON", "/proj/envs/mummichog/bin/python")
#' }
set_renviron_var <- function(key, value, path = ".Renviron") {
  if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", key)) {
    stop("Invalid environment variable name: '", key, "'.", call. = FALSE)
  }
  new_line <- paste0(key, "=", value)

  existing <- if (file.exists(path)) readLines(path, warn = FALSE) else character()

  # Match an assignment for `key`, tolerating leading whitespace and an
  # `export ` prefix (both valid in .Renviron / shell env files).
  pattern <- paste0("^\\s*(export\\s+)?", key, "\\s*=")
  hits    <- grepl(pattern, existing)

  if (any(hits)) {
    unchanged <- sum(hits) == 1L && identical(existing[hits], new_line)
    first     <- which(hits)[1]
    existing[first] <- new_line
    # Drop any further duplicate assignments of the same key.
    keep <- rep(TRUE, length(existing))
    keep[which(hits)[-1]] <- FALSE
    out    <- existing[keep]
    action <- if (unchanged) "unchanged" else "updated"
  } else {
    out    <- c(existing, new_line)
    action <- "added"
  }

  writeLines(out, path)
  invisible(action)
}
