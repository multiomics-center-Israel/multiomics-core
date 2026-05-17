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


#' Merge row_data annotation columns into a DE statistics data.frame
#'
#' Left-joins annotation columns from \code{row_data} onto \code{de_stats} by
#' matching values of \code{id_col_de} to values of \code{id_col_row}.  The row
#' count of \code{de_stats} is preserved.  Columns that already exist in
#' \code{de_stats} are skipped by default, treating the upstream
#' \code{summary_df} as the source of truth (override with
#' \code{on_collision = "overwrite"}).
#'
#' Final column order is: \code{id_col_de} -> annotation columns (in
#' \code{row_data} order, minus \code{id_col_row}) -> remaining \code{de_stats}
#' columns.  When nothing is merged (NULL/empty \code{row_data}, no annotation
#' columns, or all candidates collided under \code{"skip"}), \code{de_stats} is
#' returned unchanged.
#'
#' @param de_stats   data.frame of DE statistics.  Must contain
#'   \code{id_col_de}.
#' @param row_data   data.frame of feature annotations, or \code{NULL}.  Must
#'   contain \code{id_col_row} when non-NULL.
#' @param id_col_de  Join column in \code{de_stats}.  Default
#'   \code{"feature_id"}.
#' @param id_col_row Join column in \code{row_data}.  Default
#'   \code{"feature_id"}.
#' @param on_collision What to do when an annotation column already exists in
#'   \code{de_stats}: \code{"skip"} (default) keeps the existing column;
#'   \code{"overwrite"} replaces it with the \code{row_data} value.
#' @return A data.frame with the same number of rows as \code{de_stats}, with
#'   annotation columns merged in and reordered.
#' @keywords internal
annotate_de_stats <- function(de_stats,
                              row_data,
                              id_col_de    = "feature_id",
                              id_col_row   = "feature_id",
                              on_collision = c("skip", "overwrite")) {
    on_collision <- match.arg(on_collision)

    if (is.null(de_stats)) return(de_stats)
    if (!is.data.frame(de_stats)) {
        stop("annotate_de_stats: de_stats must be a data.frame, got ",
             class(de_stats)[1], call. = FALSE)
    }
    if (!(id_col_de %in% colnames(de_stats))) {
        stop("annotate_de_stats: id_col_de '", id_col_de,
             "' not found in de_stats. Available: ",
             paste(head(colnames(de_stats), 10), collapse = ", "),
             call. = FALSE)
    }

    if (is.null(row_data)) {
        message("annotate_de_stats: row_data is NULL or empty - no annotations merged.")
        return(de_stats)
    }
    if (!is.data.frame(row_data)) {
        stop("annotate_de_stats: row_data must be a data.frame, got ",
             class(row_data)[1], call. = FALSE)
    }
    if (nrow(row_data) == 0) {
        message("annotate_de_stats: row_data is NULL or empty - no annotations merged.")
        return(de_stats)
    }
    if (!(id_col_row %in% colnames(row_data))) {
        stop("annotate_de_stats: id_col_row '", id_col_row,
             "' not found in row_data. Available: ",
             paste(head(colnames(row_data), 10), collapse = ", "),
             call. = FALSE)
    }

    annot_cols <- setdiff(colnames(row_data), id_col_row)
    if (length(annot_cols) == 0) {
        message("annotate_de_stats: row_data has no annotation columns - nothing to merge.")
        return(de_stats)
    }

    row_keys <- as.character(row_data[[id_col_row]])
    nonNA_keys <- row_keys[!is.na(row_keys)]
    if (anyDuplicated(nonNA_keys) > 0) {
        dups <- unique(nonNA_keys[duplicated(nonNA_keys)])
        warning("annotate_de_stats: row_data has ", length(dups),
                " duplicate key(s) in '", id_col_row,
                "'; first match wins. Examples: ",
                paste(head(dups, 5), collapse = ", "), call. = FALSE)
    }

    collisions <- intersect(annot_cols, colnames(de_stats))
    if (length(collisions) > 0 && identical(on_collision, "skip")) {
        message("annotate_de_stats: skipping ", length(collisions),
                " column(s) already present in de_stats: ",
                paste(head(collisions, 10), collapse = ", "))
        annot_cols <- setdiff(annot_cols, collisions)
    }
    if (length(annot_cols) == 0) {
        return(de_stats)
    }

    de_keys <- as.character(de_stats[[id_col_de]])
    idx <- match(de_keys, row_keys)

    match_rate <- if (length(idx) > 0) sum(!is.na(idx)) / length(idx) else 1
    if (match_rate < 0.95) {
        warning(sprintf(
            "annotate_de_stats: low match rate %.1f%% (%d/%d features matched between de_stats$%s and row_data$%s).",
            match_rate * 100, sum(!is.na(idx)), length(idx),
            id_col_de, id_col_row
        ), call. = FALSE)
    }

    annot_block <- row_data[idx, annot_cols, drop = FALSE]
    rownames(annot_block) <- NULL

    if (identical(on_collision, "overwrite") && length(collisions) > 0) {
        de_stats <- de_stats[, setdiff(colnames(de_stats), collisions),
                             drop = FALSE]
    }

    out <- cbind(de_stats, annot_block)

    other <- setdiff(colnames(out), c(id_col_de, annot_cols))
    out[, c(id_col_de, annot_cols, other), drop = FALSE]
}
