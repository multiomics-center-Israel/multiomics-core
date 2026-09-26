# R/domain/de_integration/02_observed_counts.R
#
# How many values each group actually measured, per feature and contrast, so
# later steps can set aside fold changes that rest mostly on imputation.

#' Numerator and denominator groups of a contrast
#'
#' From the layer's contrasts file when it has one, matched on
#' \code{normalize_contrast_key()}. Otherwise from the groups the table itself
#' counts (its \code{N.observed.<group>} columns): the one ordered pair whose
#' "<numerator>_vs_<denominator>" has the contrast's key. That reads "A_vs_B"
#' and the space-stripped "SvsNS" the proteomics export writes for "S vs NS"
#' alike. When several pairs fit, nothing is guessed.
#'
#' @param contrast Contrast label.
#' @param contrasts_df Contrasts table (\code{Contrast_name}, \code{Factor},
#'   \code{Numerator}, \code{Denominator}), or NULL.
#' @param groups Group names the table carries observed counts for, or NULL.
#' @param layer Layer name, for the warning.
#' @return List with \code{numerator}, \code{denominator} and \code{factor}
#'   (NA where unknown), or NULL.
#' @examples
#' contrast_groups("SvsNS", groups = c("S", "NS"))$numerator   # "S"
contrast_groups <- function(contrast, contrasts_df = NULL, groups = NULL, layer = "") {
    key <- normalize_contrast_key(contrast)
    if (is.data.frame(contrasts_df) && "Contrast_name" %in% names(contrasts_df)) {
        i <- which(normalize_contrast_key(contrasts_df$Contrast_name) == key)
        if (length(i) == 1) {
            return(list(numerator = as.character(contrasts_df$Numerator[i]),
                        denominator = as.character(contrasts_df$Denominator[i]),
                        factor = as.character(contrasts_df$Factor[i] %||% NA)))
        }
    }

    groups <- unique(groups[!is.na(groups) & nzchar(groups)])
    if (length(groups) < 2) return(NULL)
    pairs <- expand.grid(num = groups, den = groups, stringsAsFactors = FALSE)
    pairs <- pairs[pairs$num != pairs$den, , drop = FALSE]
    fit <- pairs[normalize_contrast_key(paste0(pairs$num, "_vs_", pairs$den)) == key, ,
                 drop = FALSE]
    if (nrow(fit) == 1) {
        return(list(numerator = fit$num, denominator = fit$den, factor = NA_character_))
    }
    if (nrow(fit) > 1) {
        warning("Layer '", layer, "', contrast '", contrast, "': its groups are ",
                "ambiguous -- ", paste(fit$num, fit$den, sep = " vs ", collapse = "; "),
                " all fit. Observed counts are left unavailable; give the layer an ",
                "observed block with a contrasts_file to name the groups.", call. = FALSE)
    }
    NULL
}


#' Observed (non-imputed) values per group for one contrast
#'
#' A fold change resting on one measured value per group is mostly imputation,
#' and later steps flag it through \code{well_observed}. Counts come from the
#' table's own columns when it has them -- \code{N.observed.<group>} in the
#' proteomics export, the mapped \code{n_obs_num}/\code{n_obs_den} in a generic
#' table -- otherwise from the layer's unimputed matrix and sample sheet
#' (\code{observed} block), using the grouping
#' \code{compute_group_stat_columns()} applies to the exported means. Without
#' either, they are NA.
#'
#' @param ly The layer's config.
#' @param df The layer's table.
#' @param ids Feature ids, one per row of \code{df}.
#' @param contrast Contrast label.
#' @param config Full config, for path resolution.
#' @param observed Pre-read observed inputs from \code{read_observed_inputs()},
#'   or NULL.
#' @return List with \code{num}, \code{den} (integer vectors, one per row) and
#'   \code{source}.
layer_observed_counts <- function(ly, df, ids, contrast, config, observed = NULL) {
    none <- list(num = rep(NA_integer_, nrow(df)), den = rep(NA_integer_, nrow(df)),
                 source = "not available")
    where <- sprintf("Layer '%s' (%s)", ly$name, basename(ly$path))
    count_col <- function(col) {
        x <- .dei_numeric_col(df, col, where, 0, Inf, finite = TRUE,
                              what = "counts of observed values")
        if (any(!is.na(x) & x != round(x))) {
            stop(where, ": column '", col, "' must hold whole counts of observed ",
                 "values.", call. = FALSE)
        }
        as.integer(x)
    }
    from_table <- function(num_col, den_col) {
        list(num = count_col(num_col), den = count_col(den_col),
             source = paste(num_col, den_col, sep = ", "))
    }

    if (identical(ly$format, "generic")) {
        cols <- ly$columns %||% list()
        # read_de_layer() has already checked that mapped columns exist.
        if (!is.null(cols$n_obs_num)) return(from_table(cols$n_obs_num, cols$n_obs_den))
    } else {
        prefix <- de_layer_columns(ly$format, contrast)$n_obs_prefix
        if (!is.null(prefix)) {
            counted <- substring(names(df)[startsWith(names(df), prefix)], nchar(prefix) + 1L)
            grp <- contrast_groups(contrast, observed$contrasts, counted, ly$name)
            if (!is.null(grp)) {
                num_col <- paste0(prefix, grp$numerator)
                den_col <- paste0(prefix, grp$denominator)
                if (all(c(num_col, den_col) %in% names(df))) {
                    return(from_table(num_col, den_col))
                }
            }
        }
    }

    if (is.null(observed)) return(none)
    grp <- contrast_groups(contrast, observed$contrasts, layer = ly$name)
    if (is.null(grp)) {
        warning("Layer '", ly$name, "': contrast '", contrast, "' is not in ",
                "observed.contrasts_file (", basename(ly$observed$contrasts_file),
                "); its observed counts are unavailable.", call. = FALSE)
        return(none)
    }
    # A group no counted sample carries (a typo, or a group missing from this
    # matrix) would come back as an all-NA column, read as "not available".
    sheet <- observed$samplesheet
    in_matrix <- as.character(sheet[[ly$observed$sample_col]]) %in% colnames(observed$matrix)
    present <- unique(as.character(sheet[[grp$factor]][in_matrix]))
    absent <- setdiff(c(grp$numerator, grp$denominator), present)
    if (length(absent) > 0) {
        stop("Layer '", ly$name, "': observed.contrasts_file (",
             basename(ly$observed$contrasts_file), ") names group(s) ",
             paste(absent, collapse = ", "), " for contrast '", contrast, "', but no ",
             "sample of observed.matrix has that value in the sample sheet's '",
             grp$factor, "' column.", call. = FALSE)
    }
    counts <- compute_group_stat_columns(
        expr = observed$matrix, sample_meta = observed$samplesheet,
        sample_id_col = ly$observed$sample_col,
        contrasts_df = observed$contrasts[
            normalize_contrast_key(observed$contrasts$Contrast_name) ==
                normalize_contrast_key(contrast), , drop = FALSE],
        stat_fn = function(m) rowSums(!is.na(m)),
        prefix = "N.observed.")
    num_col <- paste0("N.observed.", grp$numerator)
    den_col <- paste0("N.observed.", grp$denominator)
    if (is.null(counts) || !all(c(num_col, den_col) %in% names(counts))) return(none)
    idx <- match(ids, rownames(counts))
    list(num = as.integer(counts[[num_col]][idx]), den = as.integer(counts[[den_col]][idx]),
         source = paste0("unimputed matrix ", basename(ly$observed$matrix)))
}


#' Read a layer's observed-count inputs once, and check them
#'
#' A configured \code{observed} block is a request for counts, so an input
#' that cannot be read, or lacks a column the counting needs, stops the run
#' here, naming the layer and the key -- rather than leaving every
#' \code{well_observed} NA with nothing said. Only the matrix columns the
#' sample sheet names are kept, so annotation columns beside the samples are
#' not read as samples.
#'
#' @param ly The layer's config (with an \code{observed} block).
#' @param config Full config, for path resolution.
#' @return List with \code{matrix} (features x samples, NA where unobserved),
#'   \code{samplesheet} and \code{contrasts}; NULL when the layer has no
#'   \code{observed} block.
read_observed_inputs <- function(ly, config) {
    obs <- ly$observed
    if (is.null(obs)) return(NULL)
    fail <- function(key, path, why) {
        stop("Layer '", ly$name, "': observed.", key, " (", path, ") ", why, ".",
             call. = FALSE)
    }
    read_or_fail <- function(key, reader) {
        path <- resolve_input_path(config, obs[[key]])
        out <- tryCatch(reader(path), error = function(e) {
            fail(key, path, paste("cannot be read:", conditionMessage(e)))
        })
        if (!is.data.frame(out) || nrow(out) == 0) fail(key, path, "cannot be read, or is empty")
        .dei_check_header(out, sprintf("Layer '%s': observed.%s (%s)", ly$name, key, path))
        list(df = out, path = path)
    }

    sheet <- read_or_fail("samplesheet", read_samplesheet)
    if (!obs$sample_col %in% names(sheet$df)) {
        fail("samplesheet", sheet$path,
             paste0("has no column '", obs$sample_col, "' (observed.sample_col)"))
    }
    # Samples are matched to their group by id; a repeated or blank id would
    # silently take the first row's group.
    sid <- as.character(sheet$df[[obs$sample_col]])
    if (any(is.na(sid) | !nzchar(trimws(sid)))) {
        fail("samplesheet", sheet$path, paste0("has blank sample ids in '", obs$sample_col, "'"))
    }
    if (anyDuplicated(sid)) {
        fail("samplesheet", sheet$path,
             paste0("repeats sample id(s) ", paste(unique(sid[duplicated(sid)]), collapse = ", "),
                    " in '", obs$sample_col, "'"))
    }

    contr <- read_or_fail("contrasts_file", read_table_auto)
    need <- c("Contrast_name", "Factor", "Numerator", "Denominator")
    gap <- setdiff(need, names(contr$df))
    if (length(gap) > 0) {
        fail("contrasts_file", contr$path, paste("lacks column(s):", paste(gap, collapse = ", ")))
    }
    keys <- normalize_contrast_key(as.character(contr$df$Contrast_name))
    if (any(is.na(keys) | !nzchar(keys))) {
        fail("contrasts_file", contr$path, "has a blank or missing Contrast_name")
    }
    if (anyDuplicated(keys)) {
        fail("contrasts_file", contr$path,
             paste0("names one contrast more than once (",
                    paste(unique(contr$df$Contrast_name[keys %in% keys[duplicated(keys)]]),
                          collapse = ", "), ")"))
    }
    num <- trimws(as.character(contr$df$Numerator))
    den <- trimws(as.character(contr$df$Denominator))
    if (any(is.na(num) | !nzchar(num) | is.na(den) | !nzchar(den))) {
        fail("contrasts_file", contr$path, "has a blank Numerator or Denominator")
    }
    # One group on both sides would count it twice and call the fold change
    # well observed with no second group looked at.
    same <- num == den
    if (any(same)) {
        fail("contrasts_file", contr$path,
             paste0("gives the same group as Numerator and Denominator for ",
                    paste(contr$df$Contrast_name[same], collapse = ", ")))
    }
    factors <- setdiff(unique(as.character(contr$df$Factor)), names(sheet$df))
    if (length(factors) > 0) {
        fail("contrasts_file", contr$path,
             paste0("groups by ", paste(factors, collapse = ", "),
                    ", which the sample sheet has no column for"))
    }

    mat <- read_or_fail("matrix", read_table_auto)
    id_col <- obs$id_col %||% names(mat$df)[1]
    if (!id_col %in% names(mat$df)) {
        fail("matrix", mat$path, paste0("has no column '", id_col, "' (observed.id_col)"))
    }
    samples <- intersect(setdiff(names(mat$df), id_col),
                         as.character(sheet$df[[obs$sample_col]]))
    if (length(samples) == 0) {
        fail("matrix", mat$path, paste0("has no column named after a sample in the ",
                                        "sample sheet's '", obs$sample_col, "'"))
    }
    # A sample with nothing measured reads back as an all-NA logical column.
    not_num <- samples[!vapply(mat$df[samples], function(x) is.numeric(x) || all(is.na(x)),
                               logical(1))]
    if (length(not_num) > 0) {
        fail("matrix", mat$path, paste("has non-numeric sample column(s):",
                                       paste(not_num, collapse = ", ")))
    }
    # A count is of values present; Inf (a log of zero, say) is not a
    # measurement, and !is.na() would count it as one.
    not_finite <- samples[vapply(mat$df[samples], function(x) {
        any(!is.na(x) & !is.finite(as.numeric(x)))
    }, logical(1))]
    if (length(not_finite) > 0) {
        fail("matrix", mat$path, paste0("has infinite values in sample column(s) ",
                                        paste(not_finite, collapse = ", "),
                                        "; write unmeasured values as NA"))
    }
    # Counts are looked up by feature id; a repeated id would silently give
    # the first row's counts to a feature whose DE row came from another.
    fid <- as.character(mat$df[[id_col]])
    if (any(is.na(fid) | !nzchar(trimws(fid)))) {
        fail("matrix", mat$path, paste0("has blank feature ids in '", id_col, "'"))
    }
    if (anyDuplicated(fid)) {
        fail("matrix", mat$path,
             paste0("repeats feature id(s) ",
                    paste(utils::head(unique(fid[duplicated(fid)]), 5), collapse = ", "),
                    " in '", id_col, "'"))
    }
    m <- as.matrix(mat$df[, samples, drop = FALSE])
    rownames(m) <- fid

    list(matrix = m, samplesheet = sheet$df, contrasts = contr$df)
}
