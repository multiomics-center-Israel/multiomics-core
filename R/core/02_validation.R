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

assert_expr_meta_alignment <- function(expr_mat, meta, cfg, strict = TRUE) {
    sample_col <- get_sample_col(cfg)

    assert_numeric_matrix(expr_mat, "expr_mat")
    assert_meta_contract(meta, sample_col)

    expr_ids <- colnames(expr_mat)
    meta_ids <- as.character(meta[[sample_col]])

    missing_in_expr <- setdiff(meta_ids, expr_ids)
    if (length(missing_in_expr) > 0) {
        stop(
            "expr_mat is missing samples from meta (sample_col='", sample_col, "'): ",
            paste(head(missing_in_expr, 10), collapse = ", "),
            if (length(missing_in_expr) > 10) sprintf(" ... (+%d more)", length(missing_in_expr) - 10) else ""
        )
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

assert_pre_contract <- function(pre, stage = c("proteomics", "rna", "metabolomics")) {
    stage <- match.arg(stage)
    stopifnot(is.list(pre))

    base_fields <- c("expr_raw", "expr_filt", "expr_work", "meta", "row_data")
    stage_fields <- switch(stage,
        proteomics = c("expr_imp_single"),
        rna = character(0),
        metabolomics = character(0)
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

assert_de_contract <- function(de_res, stage = c("proteomics", "rna", "metabolomics")) {
    stage <- match.arg(stage)
    stopifnot(is.list(de_res))

    base_fields <- c("summary_df")
    stage_fields <- switch(stage,
        proteomics = c("method", "imputations"),
        rna = character(0),
        metabolomics = character(0)
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
