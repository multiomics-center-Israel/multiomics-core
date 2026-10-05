#' Building blocks for the per-run results fact sheet
#'
#' The fact sheet is a flat table of the headline numbers a run produced, each
#' one paired with the file it can be checked against. Its purpose is to let a
#' reader recompute any figure that ends up in a report or a slide without
#' having to work out which artifact it came from.
#'
#' Rows are assembled per omics mode (see \code{build_rnaseq_fact_sheet}) from
#' the files already written to the results directory, rather than from objects
#' held in memory. That is deliberate: the sheet is a statement about what is on
#' disk, so reading the artifacts is what makes it true.

#' Make one fact-sheet row
#'
#' @param claim Short description of the quantity, in lower case.
#' @param value The value, coerced to a single string. Free text is fine
#'   (\code{"1.68 to 351, median 2.61"}); the sheet is read by people.
#' @param source_file Path of the artifact the value came from, relative to the
#'   mode's results directory.
#' @return A one-row data.frame, or NULL when \code{value} is NULL or NA, so
#'   that a missing input drops the row instead of writing an empty one.
fact <- function(claim, value, source_file) {
    if (is.null(value) || length(value) == 0) return(NULL)
    if (length(value) == 1 && is.na(value)) return(NULL)
    data.frame(
        claim       = as.character(claim),
        value       = paste(as.character(value), collapse = ", "),
        source_file = as.character(source_file),
        stringsAsFactors = FALSE
    )
}

#' Combine fact rows into a sheet
#'
#' @param ... One-row data.frames from \code{\link{fact}}, or lists of them.
#'   NULLs are dropped, so a caller can pass a row that may not apply.
#' @return A data.frame with columns claim, value, source_file; zero rows if
#'   nothing was supplied.
bind_facts <- function(...) {
    rows <- unlist(list(...), recursive = FALSE, use.names = FALSE)
    rows <- Filter(function(r) is.data.frame(r) && nrow(r) > 0, rows)
    if (length(rows) == 0) {
        return(data.frame(claim = character(0), value = character(0),
                          source_file = character(0), stringsAsFactors = FALSE))
    }
    do.call(rbind, rows)
}

#' Write a fact sheet to the results directory
#'
#' @param sheet Data frame from \code{\link{bind_facts}}.
#' @param out_dir The mode's results directory.
#' @param filename File name to write.
#' @return Path to the written file, or character(0) if the sheet was empty.
write_results_fact_sheet <- function(sheet, out_dir, filename = "results_fact_sheet.tsv") {
    if (is.null(sheet) || nrow(sheet) == 0) {
        warning("results fact sheet: nothing to write; skipping.")
        return(character(0))
    }
    save_tsv(sheet, out_dir, filename)
}

#' Read a delimited results artifact, returning NULL if it is absent
#'
#' The fact sheet is assembled opportunistically: a run without enrichment, or
#' one stopped early, should still produce a sheet covering what it did make.
#'
#' @param path File path.
#' @param sep Field separator.
#' @return A data.frame, or NULL when the file is missing or unreadable.
read_artifact_or_null <- function(path, sep = "\t") {
    if (length(path) != 1 || is.na(path) || !file.exists(path)) return(NULL)
    tryCatch(
        utils::read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE),
        error = function(e) {
            warning(sprintf("results fact sheet: could not read '%s' (%s); skipping its rows.",
                            basename(path), conditionMessage(e)))
            NULL
        }
    )
}

#' Format a numeric range compactly for a fact-sheet value
#'
#' @param x Numeric vector; non-finite values are dropped.
#' @param digits Significant digits.
#' @return A string such as "1.68 to 351, median 2.61", or NULL if x is empty.
fmt_range <- function(x, digits = 3) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NULL)
    sprintf("%s to %s, median %s",
            signif(min(x), digits), signif(max(x), digits),
            signif(stats::median(x), digits))
}
