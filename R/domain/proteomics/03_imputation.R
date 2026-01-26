#' Impute proteomics expr_mat using Perseus-style method
impute_proteomics_perseus <- function(expr_mat, cfg, return_flags = FALSE) {
    width <- cfg$imputation$width %||% 0.3
    downshift <- cfg$imputation$downshift %||% 1.8

    res <- perseus_impute_with_flags(
        expr_mat = expr_mat,
        width = width,
        downshift = downshift
    )

    if (return_flags) res else res$imputed
}

#' Wrapper to run multiple imputations (with seed increments)
make_imputations_proteomics <- function(expr_mat, cfg, verbose = FALSE) {
    imp_cfg <- cfg$modes$proteomics$imputation
    n_imputations <- as.integer(imp_cfg$no_repetitions)
    seed_base <- cfg$params$seed

    stopifnot(is.matrix(expr_mat))

    imps <- vector("list", n_imputations)
    for (i in seq_len(n_imputations)) {
        if (isTRUE(verbose)) message(sprintf("Imputation: %d / %d", i, n_imputations))
        set.seed(as.integer(seed_base) + i)

        expr_imp_i <- impute_proteomics_perseus(
            expr_mat,
            cfg = cfg$modes$proteomics,
            return_flags = FALSE
        )
        imps[[i]] <- expr_imp_i
    }
    imps
}

#' Perseus-style imputation (downshifted & narrowed normal distribution)
perseus_impute_with_flags <- function(expr_mat, width = 0.3, downshift = 1.8) {
    expr_mat <- as.matrix(expr_mat)
    sample_cols <- colnames(expr_mat)
    imputed_flag <- is.na(expr_mat)
    imputed <- expr_mat

    for (j in seq_len(ncol(imputed))) {
        x <- imputed[, j]
        if (all(is.na(x))) stop("Imputation failed: sample '", colnames(imputed)[j], "' is all-NA.")

        obs <- x[!is.na(x)]
        s <- stats::sd(obs)
        if (!is.finite(s) || s == 0) s <- 1e-8
        m <- mean(obs)

        imp_sd <- width * s
        imp_mean <- m - downshift * s

        n_missing <- sum(is.na(x))
        if (n_missing > 0) {
            x[is.na(x)] <- stats::rnorm(n_missing, mean = imp_mean, sd = imp_sd)
        }
        imputed[, j] <- x
    }

    list(imputed = imputed, imputed_flag = imputed_flag)
}

validate_proteomics_imputations <- function(imputations, meta, cfg, warn_extra_meta = TRUE, allow_na = FALSE, sample_col = NULL) {
    if (is.null(sample_col) || !nzchar(sample_col)) {
        p <- cfg$modes$proteomics
        sample_col <- p$effects$samples %||% p$id_columns$sample_col %||% "SampleID"
    }

    if (!is.list(imputations) || length(imputations) < 1) stop("imputations must be a non-empty list.")

    # Validate ref
    ref <- imputations[[1]]
    assert_numeric_matrix(ref, "imputations[[1]]")

    # Check alignment of all matrices
    ref_rn <- rownames(ref)
    ref_cn <- colnames(ref)

    for (i in seq_along(imputations)) {
        x <- imputations[[i]]
        if (!identical(dim(x), dim(ref))) stop(sprintf("Imputation %d dim mismatch.", i))
        if (!identical(rownames(x), ref_rn)) stop(sprintf("Imputation %d rowname mismatch.", i))
        if (!identical(colnames(x), ref_cn)) stop(sprintf("Imputation %d colname mismatch.", i))
        if (!isTRUE(allow_na) && anyNA(x)) stop(sprintf("Imputation %d has NAs.", i))
    }

    # Check meta alignment
    meta_ids <- as.character(meta[[sample_col]])
    missing_in_meta <- setdiff(ref_cn, meta_ids)
    if (length(missing_in_meta) > 0) stop("Meta missing samples present in imputations.")

    invisible(TRUE)
}
