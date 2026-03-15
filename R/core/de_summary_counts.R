#' Generic DE summary counts builder (internal)
#'
#' Counts up-/down-/total significant features per contrast from a DE
#' statistics data.frame.  Domain-specific wrappers supply the naming
#' conventions and significance logic via callback arguments.
#'
#' @param de_stats  DE statistics data.frame with pass and FC columns.
#' @param pass_pattern  Regex that matches per-contrast pass columns
#'   (e.g. \code{"^pass\\\\.imputs\\\\."} or \code{"_pass$"}).
#' @param extract_contrast  \code{function(col_name)} returning the contrast
#'   name derived from the pass column name.
#' @param find_fc_col  \code{function(contrast_name, col_names)} returning the
#'   matching fold-change column name, or \code{NULL} if none found.
#' @param is_significant  \code{function(x)} returning a logical vector
#'   indicating which rows are significant.  Receives the raw pass-column
#'   values (may contain \code{NA}).  Default: \code{!is.na(x) & x == 1}.
#' @return A \code{data.frame} with columns \code{contrast}, \code{up},
#'   \code{down}, \code{total}; or \code{NULL} when no pass columns exist.
#'   If \code{pass_any_contrast} is present in \code{de_stats}, an extra
#'   \code{"any"} row is appended with the total number of features
#'   significant in any contrast (up/down are 0).
#' @keywords internal
build_de_summary_counts_generic <- function(de_stats,
                                            pass_pattern,
                                            extract_contrast,
                                            find_fc_col,
                                            is_significant = function(x) !is.na(x) & x == 1) {
    if (is.null(de_stats)) return(NULL)

    pass_cols <- grep(pass_pattern, names(de_stats), value = TRUE)
    pass_cols <- setdiff(pass_cols, "pass_any_contrast")

    if (length(pass_cols) == 0) return(NULL)

    summaries <- lapply(pass_cols, function(col) {
        contrast_name <- extract_contrast(col)
        fc_col <- find_fc_col(contrast_name, names(de_stats))

        sig_rows <- is_significant(de_stats[[col]])

        if (!is.null(fc_col)) {
            up   <- sum(sig_rows & de_stats[[fc_col]] > 0, na.rm = TRUE)
            down <- sum(sig_rows & de_stats[[fc_col]] < 0, na.rm = TRUE)
        } else {
            up   <- NA
            down <- NA
        }

        data.frame(
            contrast = contrast_name,
            up       = up,
            down     = down,
            total    = sum(sig_rows, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    })

    result <- do.call(rbind, summaries)

    # Append "any" row if pass_any_contrast exists
    if ("pass_any_contrast" %in% names(de_stats)) {
        any_total <- sum(
            is_significant(de_stats$pass_any_contrast)
        )
        any_row <- data.frame(
            contrast = "any",
            up       = 0,
            down     = 0,
            total    = any_total,
            stringsAsFactors = FALSE
        )
        result <- rbind(result, any_row)
    }

    result
}
