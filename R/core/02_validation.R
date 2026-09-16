# Error/warning prefix convention in this file (and validators called from
# here): use [<mode>] for mode-scoped messages (e.g. [rna]) and [<scope>] for
# cross-cutting ones (e.g. [design], [counts]). Do not use [<function_name>].

#' Check that required columns exist in a data frame
check_has_cols <- function(df, required, df_name = deparse(substitute(df))) {
    required <- unlist(required)
    missing <- setdiff(required, colnames(df))
    if (length(missing) > 0) {
        stop(
            "In table '", df_name, "': missing columns: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
}

#' Warn when a tab-separated sample sheet has commas inside its values
#'
#' The pipeline reads metadata by extension (see \code{read_table_auto}), so a
#' tab-separated sheet with commas in its cells parses correctly here. Readers
#' that assume CSV do not: they see more fields in the data rows than in the
#' header, silently decide column 1 holds row names, and abort with
#' "duplicate 'row.names' are not allowed". That is how a comma in a sample
#' sheet takes down the report render long after the analysis itself succeeded,
#' with no traceback pointing back at the sheet. Flag it at load time instead.
#'
#' Only column names and affected-row counts are reported — never cell values,
#' which may carry identifying information.
#'
#' @param df Parsed metadata data frame.
#' @param path Path the metadata was read from; its extension decides whether
#'   the file is tab-separated and therefore at risk.
#' @param mode Omics mode or scope, used as the message prefix.
#' @return Invisibly \code{TRUE}. Emits a warning when commas are found.
#' @examples
#' df <- data.frame(SampleName = c("S1", "S2"),
#'                  Source = c("run-1,run-2", "run-3,run-4"))
#' # check_metadata_delimiter_safety(df, "samples.txt", "proteomics")
check_metadata_delimiter_safety <- function(df, path, mode = "design") {
    if (!is.data.frame(df) || nrow(df) == 0) return(invisible(TRUE))
    # Only tab-separated sheets are at risk; a real CSV quotes its commas.
    if (!tolower(tools::file_ext(path)) %in% c("tsv", "txt")) return(invisible(TRUE))

    chr_cols <- names(df)[vapply(df, is.character, logical(1))]
    n_affected <- vapply(
        chr_cols,
        function(cn) sum(grepl(",", df[[cn]], fixed = TRUE), na.rm = TRUE),
        integer(1)
    )
    hits <- chr_cols[n_affected > 0]
    if (length(hits) == 0) return(invisible(TRUE))

    warning(
        sprintf(
            paste0(
                "[%s] Metadata '%s' is tab-separated but %d column(s) contain commas: %s.\n",
                "  This pipeline parses it correctly, but any reader that assumes CSV — ",
                "including the report templates — will mis-parse it and abort with\n",
                "  \"duplicate 'row.names' are not allowed\".\n",
                "  Fix: remove the commas from those values (e.g. rename 'A,B' to 'A_B'), ",
                "or drop the column if the analysis does not use it."
            ),
            mode, basename(path), length(hits),
            paste(sprintf("%s (%d/%d rows)", hits, n_affected[hits], nrow(df)),
                  collapse = ", ")
        ),
        call. = FALSE
    )
    invisible(TRUE)
}

#' Check that all values in x are present in y
check_all_in <- function(x, y, label_x = "x", label_y = "y") {
    missing <- setdiff(x, y)
    if (length(missing) > 0) {
        stop(
            "Values in ", label_x, " not found in ", label_y, ": ",
            paste(head(missing, 10), collapse = ", "),
            if (length(missing) > 10) " ... (truncated)" else "",
            call. = FALSE
        )
    }
}

#' Helper to extract sample column from config
get_sample_col <- function(cfg) {
    sample_col <- cfg$effects$samples %||% cfg$id_columns$sample_col
    if (is.null(sample_col) || !nzchar(sample_col)) stop("Missing sample column in cfg$effects$samples / cfg$id_columns$sample_col")
    sample_col
}

# Helper for %||% (duplicated here or relies on core loading order, better strictly define/import if possible, but for tar_source likely fine if defined once available project-wide.
# However, to be safe, I'll define it here if not already defined or assume it's available.
# Actually, since these are independent files, I should define it if used inside functions?
# Only get_sample_col uses it. I'll include it here locally or assume common utils. Using local define for safety.)
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

#' Assert a strict numeric matrix contract
assert_numeric_matrix <- function(x, name = "x") {
    if (!is.matrix(x)) stop(sprintf("`%s` must be a matrix.", name))
    if (!is.numeric(x)) stop(sprintf("`%s` must be numeric.", name))
    if (is.null(rownames(x)) || anyNA(rownames(x)) || any(rownames(x) == "")) {
        stop(sprintf("`%s` must have non-empty rownames (feature IDs).", name))
    }
    if (is.null(colnames(x)) || anyNA(colnames(x)) || any(colnames(x) == "")) {
        stop(sprintf("`%s` must have non-empty colnames (sample IDs).", name))
    }
    if (anyDuplicated(rownames(x)) > 0) stop(sprintf("`%s` has duplicated rownames.", name))
    if (anyDuplicated(colnames(x)) > 0) stop(sprintf("`%s` has duplicated colnames.", name))
    invisible(TRUE)
}

assert_meta_contract <- function(meta, sample_col) {
    if (!is.data.frame(meta)) stop("meta must be a data.frame.")
    check_has_cols(meta, sample_col, df_name = "meta")
    samp <- as.character(meta[[sample_col]])
    if (any(!nzchar(samp))) stop("meta$", sample_col, " contains empty sample IDs.")
    if (anyDuplicated(samp) > 0) stop("meta$", sample_col, " contains duplicated sample IDs.")
    invisible(TRUE)
}

assert_expr_meta_alignment <- function(expr_mat, meta, cfg, strict = FALSE) {
    sample_col <- get_sample_col(cfg)

    assert_numeric_matrix(expr_mat, "expr_mat")
    assert_meta_contract(meta, sample_col)

    expr_ids <- colnames(expr_mat)
    meta_ids <- as.character(meta[[sample_col]])

    missing_in_expr <- setdiff(meta_ids, expr_ids)
    if (length(missing_in_expr) > 0) {
        msg <- paste0(
            "expr_mat is missing samples from meta (sample_col='", sample_col, "'): ",
            paste(head(missing_in_expr, 10), collapse = ", "),
            if (length(missing_in_expr) > 10) sprintf(" ... (+%d more)", length(missing_in_expr) - 10) else ""
        )
        if (isTRUE(strict)) stop(msg) else warning(msg)
    }

    if (isTRUE(strict)) {
        extra_in_expr <- setdiff(expr_ids, meta_ids)
        if (length(extra_in_expr) > 0) {
            stop(
                "expr_mat has samples not present in meta (sample_col='", sample_col, "'): ",
                paste(head(extra_in_expr, 10), collapse = ", "),
                if (length(extra_in_expr) > 10) sprintf(" ... (+%d more)", length(extra_in_expr) - 10) else ""
            )
        }
    }

    invisible(TRUE)
}

assert_omics_obj <- function(obj, stage = c("raw", "work"), sample_col) {
    stage <- match.arg(stage)

    if (!is.list(obj)) stop("omics obj must be a list.")
    if (is.null(obj$assay) || !is.list(obj$assay)) stop("omics obj missing $assay list.")
    if (is.null(obj$col_data) || !is.data.frame(obj$col_data)) stop("omics obj missing $col_data (data.frame).")

    if (missing(sample_col) || is.null(sample_col) || !nzchar(sample_col)) {
        stop("assert_omics_obj: 'sample_col' must be provided (no default).")
    }

    mat <- obj$assay[[stage]]
    if (is.null(mat)) stop("omics obj missing assay$", stage)

    assert_numeric_matrix(mat, paste0("assay$", stage))
    assert_meta_contract(obj$col_data, sample_col)

    meta_ids <- as.character(obj$col_data[[sample_col]])
    if (!identical(meta_ids, colnames(mat))) {
        missing_in_mat <- setdiff(meta_ids, colnames(mat))
        extra_in_mat <- setdiff(colnames(mat), meta_ids)

        stop(
            "assay$", stage, " and col_data are not aligned.\n",
            "- Expected colnames(mat) to match col_data[['", sample_col, "']] exactly.\n",
            if (length(missing_in_mat)) paste0("  Missing in mat: ", paste(head(missing_in_mat, 10), collapse = ", "), "\n") else "",
            if (length(extra_in_mat)) paste0("  Extra in mat: ", paste(head(extra_in_mat, 10), collapse = ", "), "\n") else "",
            "Tip: call align_meta_to_expr(mat, col_data, sample_col) upstream."
        )
    }
    invisible(TRUE)
}

assert_pre_contract <- function(pre, stage = c("proteomics", "rna", "metabolomics", "lipidomics")) {
    stage <- match.arg(stage)
    stopifnot(is.list(pre))

    base_fields <- c("expr_raw", "expr_filt", "expr_work", "meta", "row_data")
    stage_fields <- switch(stage,
        proteomics = c("expr_imp_single"),
        rna = character(0),
        metabolomics = character(0),
        lipidomics = character(0)
    )

    required <- c(base_fields, stage_fields)
    missing <- setdiff(required, names(pre))

    if (length(missing) > 0) {
        stop(
            sprintf(
                "Preprocess contract failed for %s. Missing fields: %s",
                stage, paste(missing, collapse = ", ")
            ),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

assert_de_contract <- function(de_res, stage = c("proteomics", "rna", "metabolomics", "lipidomics")) {
    stage <- match.arg(stage)
    stopifnot(is.list(de_res))

    base_fields <- c("summary_df")
    stage_fields <- switch(stage,
        proteomics = c("method", "imputations"),
        rna = character(0),
        metabolomics = character(0),
        lipidomics = character(0)
    )

    required <- c(base_fields, stage_fields)
    missing <- setdiff(required, names(de_res))

    if (length(missing) > 0) {
        stop(
            sprintf(
                "DE contract failed for %s. Missing fields: %s",
                stage, paste(missing, collapse = ", ")
            ),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

coerce_df_to_numeric_matrix <- function(df, rownames_vec = NULL, name = "df") {
    if (is.matrix(df)) {
        return(df)
    }
    if (!is.data.frame(df)) stop(sprintf("`%s` must be a data.frame or matrix.", name))

    m <- as.matrix(df)
    suppressWarnings(storage.mode(m) <- "numeric")

    if (!is.null(rownames_vec)) {
        if (length(rownames_vec) != nrow(m)) {
            stop(sprintf(
                "`rownames_vec` length (%d) != nrow(%s) (%d).",
                length(rownames_vec), name, nrow(m)
            ))
        }
        rownames(m) <- as.character(rownames_vec)
    }
    m
}

assert_named_list <- function(x, name) {
    if (is.null(x) || !is.list(x)) stop(sprintf("'%s' must be a list", name))
    invisible(TRUE)
}

assert_scalar_bool <- function(x, name, allow_null = FALSE) {
    if (allow_null && is.null(x)) {
        return(invisible(TRUE))
    }
    if (!is.logical(x) || length(x) != 1 || is.na(x)) {
        stop(sprintf("'%s' must be TRUE/FALSE", name))
    }
    invisible(TRUE)
}

assert_scalar_chr <- function(x, name, allow_null = FALSE) {
    if (allow_null && is.null(x)) {
        return(invisible(TRUE))
    }
    if (!is.character(x) || length(x) != 1 || is.na(x) || x == "") {
        stop(sprintf("'%s' must be a non-empty string", name))
    }
    invisible(TRUE)
}

assert_scalar_num <- function(x, name, allow_null = FALSE, min_val = -Inf, max_val = Inf) {
    if (allow_null && is.null(x)) {
        return(invisible(TRUE))
    }
    if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
        stop(sprintf("'%s' must be a number", name))
    }
    if (x < min_val || x > max_val) {
        stop(sprintf("'%s' must be between %s and %s", name, min_val, max_val))
    }
    invisible(TRUE)
}

#' Convert signed linear fold change to log2 fold change
#'
#' Proteomics linear FC convention: positive values = up, negative = down
#' (e.g. 2 = 2x up, -1.5 = 1.5x down). Converts to log2 space preserving sign.
#' @param fc Numeric vector of signed linear fold changes
#' @return Numeric vector of log2 fold changes
signed_fc_to_log2 <- function(fc) {
    fc <- as.numeric(fc)
    ifelse(is.na(fc) | fc == 0, NA_real_, log2(abs(fc)) * sign(fc))
}

assert_one_of <- function(x, name, choices, allow_null = FALSE) {
    if (allow_null && is.null(x)) {
        return(invisible(TRUE))
    }
    assert_scalar_chr(x, name)
    if (!(tolower(x) %in% tolower(choices))) {
        stop(sprintf("'%s' must be one of: %s", name, paste(choices, collapse = ", ")))
    }
    invisible(TRUE)
}

#' Assert PCA Scores Contract
#'
#' Validates that an object conforms to the PCA scores data.frame contract.
#' Use at: QC module return points, when storing PCA scores in objects list.
#'
#' @param x Object to validate (expected: data.frame with PC1, PC2, ... columns)
#' @param context Character string for error messages (e.g., "proteomics QC")
#' @return x (invisibly) if valid, otherwise stops with informative error
assert_pca_scores <- function(x, context = "PCA scores") {
    prefix <- sprintf("[%s]", context)

    if (is.null(x)) {
        return(invisible(NULL))
    }

    if (!is.data.frame(x)) {
        stop(prefix, " Expected data.frame, got ", class(x)[1], call. = FALSE)
    }

    required_cols <- c("PC1", "PC2")
    missing <- setdiff(required_cols, colnames(x))
    if (length(missing) > 0) {
        stop(prefix, " Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    if (!is.numeric(x$PC1) || !is.numeric(x$PC2)) {
        stop(prefix, " PC1 and PC2 must be numeric", call. = FALSE)
    }

    if (nrow(x) == 0) {
        warning(prefix, " PCA scores has 0 rows (no samples)")
    }

    invisible(x)
}

#' Assert DE Expression Matrix Contract
#'
#' Validates that an object conforms to the DE expression matrix contract.
#' Use at: Shiny export layer, before assigning to legacy$mat2plot or legacy$de_expr_norm.
#'
#' @param x Object to validate (expected: numeric matrix with rownames and colnames)
#' @param context Character string for error messages
#' @return x (invisibly) if valid, otherwise stops with informative error
assert_de_expr_matrix <- function(x, context = "DE expression matrix") {
    prefix <- sprintf("[%s]", context)

    if (is.null(x)) {
        return(invisible(NULL))
    }

    if (!is.matrix(x)) {
        stop(prefix, " Expected matrix, got ", class(x)[1], call. = FALSE)
    }

    if (!is.numeric(x)) {
        stop(prefix, " Matrix must be numeric, got ", typeof(x), call. = FALSE)
    }

    if (is.null(rownames(x))) {
        stop(prefix, " Matrix must have rownames (feature IDs)", call. = FALSE)
    }

    if (is.null(colnames(x))) {
        stop(prefix, " Matrix must have colnames (sample IDs)", call. = FALSE)
    }

    if (nrow(x) == 0) {
        warning(prefix, " Matrix has 0 rows (no DE features)")
    }

    invisible(x)
}

#' Assert that an experimental design is statistically viable for DE.
#'
#' Validates group-replicate counts in the sample metadata against thresholds
#' for differential-expression feasibility. Halts the pipeline on
#' configurations that make DE impossible; warns on configurations that are
#' technically runnable but statistically unreliable.
#'
#' Thresholds (per resolved group column):
#' \itemize{
#'   \item Any group with n < 2 : \code{stop()} — DE cannot be performed.
#'   \item Group with n == 2    : \code{warning()} — severe; one composite
#'         warning that subsumes the low-power case.
#'   \item Group with 3 <= n < 5: \code{warning()} — low statistical power.
#'   \item max(n) / min(n) > 5  : \code{warning()} — heavily imbalanced
#'         design (emitted once globally, separate from per-group warnings).
#' }
#'
#' @param meta Sample metadata as a \code{data.frame}.
#' @param cfg  Mode-level config list (e.g. \code{config$modes$rna}). The
#'   group column is resolved via \code{\link{resolve_group_col}}.
#' @return \code{invisible(TRUE)} on viable designs. Stops with an
#'   informative error on infeasible ones; emits at most one warning
#'   per underpowered group, plus an optional global imbalance warning.
#'
#' @note Caller contract: \code{meta} MUST be the post-\code{apply_sample_filter}
#'   analysis cohort — i.e. the metadata that will actually enter DE — not
#'   the raw input metadata. Design viability is a property of the cohort
#'   that survives filtering; calling this function on unfiltered metadata
#'   produces false-positive failures for configs that intentionally drop
#'   singleton or NA-labeled groups via \code{cfg$sample_filter$rules}.
#'   Canonical call site: \code{R/domain/rnaseq/03_preprocess.R}, immediately
#'   after the sample-filter block.
#' @note If the resolved group column is absent from \code{meta}, this
#'   function \code{stop()}s with the list of available columns. This makes
#'   \code{assert_design_viable} the canonical gate for group-column
#'   presence in the rnaseq path; no upstream validator currently performs
#'   this check.
#' @note NA values in the resolved group column also cause a \code{stop()}
#'   — they are not silently dropped and not bucketed as a virtual NA group.
#'   The user must assign a label or remove the offending samples upstream.
#'
#' @examples
#' meta <- data.frame(sample = paste0("S", 1:6),
#'                    Condition = rep(c("ctrl", "trt"), each = 3))
#' cfg <- list(filtering = list(group_col = "Condition"))
#' assert_design_viable(meta, cfg)  # invisible(TRUE) + one low-power warning
assert_design_viable <- function(meta, cfg) {
    if (!is.data.frame(meta)) {
        stop("[design] meta must be a data.frame.", call. = FALSE)
    }

    group_col <- resolve_group_col(cfg)

    if (!(group_col %in% names(meta))) {
        stop(sprintf(
            "[design] Resolved group column '%s' not found in metadata. Available columns: %s.",
            group_col, paste(names(meta), collapse = ", ")
        ), call. = FALSE)
    }

    grp <- as.character(meta[[group_col]])
    na_idx <- which(is.na(grp) | !nzchar(grp))
    if (length(na_idx) > 0) {
        sample_col <- tryCatch(get_sample_col(cfg), error = function(e) NULL)
        sample_ids <- if (!is.null(sample_col) && sample_col %in% names(meta)) {
            as.character(meta[[sample_col]])[na_idx]
        } else {
            as.character(na_idx)
        }
        shown <- head(sample_ids, 10)
        more <- length(sample_ids) - length(shown)
        stop(sprintf(
            "[design] Group column '%s' has NA values for %d sample(s): %s%s. Assign a group label or remove these samples from the metadata before running DE.",
            group_col, length(sample_ids), paste(shown, collapse = ", "),
            if (more > 0) sprintf(" ... (+%d more)", more) else ""
        ), call. = FALSE)
    }

    group_sizes <- table(grp)

    critical <- names(group_sizes)[group_sizes < 2]
    if (length(critical) > 0) {
        if (length(critical) == 1L) {
            stop(sprintf(
                "[design] Group '%s' has n=%d; differential expression requires n >= 2 per group.",
                critical, as.integer(group_sizes[[critical]])
            ), call. = FALSE)
        } else {
            details <- paste(sprintf("'%s' (n=%d)", critical, as.integer(group_sizes[critical])), collapse = ", ")
            stop(sprintf(
                "[design] Groups with n < 2: %s. DE cannot be performed.",
                details
            ), call. = FALSE)
        }
    }

    for (g in names(group_sizes)) {
        n <- as.integer(group_sizes[[g]])
        if (n == 2L) {
            warning(sprintf(
                "[design] Group '%s' has n=2: severe — within-group variance cannot be reliably estimated; DE results will be unreliable.",
                g
            ), call. = FALSE)
        } else if (n < 5L) {
            warning(sprintf(
                "[design] Group '%s' has n=%d: low statistical power (recommended n >= 5).",
                g, n
            ), call. = FALSE)
        }
    }

    if (length(group_sizes) >= 2L) {
        max_idx <- which.max(group_sizes)
        min_idx <- which.min(group_sizes)
        ratio <- as.numeric(group_sizes[[max_idx]]) / as.numeric(group_sizes[[min_idx]])
        if (ratio > 5) {
            warning(sprintf(
                "[design] Group sizes imbalanced: max/min ratio = %.1f ('%s' n=%d vs '%s' n=%d). DE estimates may be biased toward the larger group.",
                ratio,
                names(group_sizes)[max_idx], as.integer(group_sizes[[max_idx]]),
                names(group_sizes)[min_idx], as.integer(group_sizes[[min_idx]])
            ), call. = FALSE)
        }
    }

    invisible(TRUE)
}
