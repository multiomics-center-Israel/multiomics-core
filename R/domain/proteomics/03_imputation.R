#' Dispatch imputation based on config method
#'
#' Routes to the correct imputation function based on cfg$imputation$method.
#' @param expr_mat numeric matrix (proteins x samples), may contain NAs
#' @param cfg      mode-level config (config$modes$proteomics)
#' @param return_flags if TRUE, return list(imputed, imputed_flag)
#' @return imputed matrix or list(imputed, imputed_flag)
impute_proteomics <- function(expr_mat, cfg, return_flags = FALSE) {
    method <- cfg$imputation$method %||% "perseus_like"
    if (method == "perseus") method <- "perseus_like"

    if (method == "none") {
        # No imputation: return as-is (NAs remain)
        if (return_flags) {
            return(list(imputed = expr_mat, imputed_flag = is.na(expr_mat)))
        }
        return(expr_mat)
    }

    if (method == "perseus_like") {
        return(impute_proteomics_perseus_like(expr_mat, cfg, return_flags))
    }

    if (method == "dep2") {
        return(impute_proteomics_dep2(expr_mat, cfg, return_flags))
    }

    if (method == "qrilc") {
        return(impute_proteomics_qrilc(expr_mat, cfg, return_flags))
    }

    if (method == "minval") {
        return(impute_proteomics_minval(expr_mat, cfg, return_flags))
    }

    stop(sprintf("Unknown imputation method: '%s'. Supported: 'none', 'perseus_like', 'dep2', 'qrilc', 'minval'.", method))
}

#' Impute proteomics expr_mat using Perseus-like method
impute_proteomics_perseus_like <- function(expr_mat, cfg, return_flags = FALSE) {
    width <- cfg$imputation$width %||% 0.3
    downshift <- cfg$imputation$downshift %||% 1.8

    res <- perseus_like_impute_with_flags(
        expr_mat = expr_mat,
        width = width,
        downshift = downshift
    )

    if (return_flags) res else res$imputed
}

#' Impute proteomics expr_mat using DEP2 / MinDet method
#'
#' Deterministic minimum-based imputation: replaces each missing value with
#' a small value drawn from the lower tail of the observed distribution per sample.
#' MinDet uses the q-th quantile of observed values per sample as the imputed value.
#' MinProb draws from a Gaussian centered at that quantile (stochastic variant).
#'
#' @param expr_mat numeric matrix (proteins x samples), may contain NAs
#' @param cfg      mode-level config (config$modes$proteomics)
#' @param return_flags if TRUE, return list(imputed, imputed_flag)
#' @return imputed matrix or list(imputed, imputed_flag)
impute_proteomics_dep2 <- function(expr_mat, cfg, return_flags = FALSE) {
    dep2_method <- cfg$imputation$dep2_method %||% "MinDet"
    dep2_seed <- as.integer(cfg$imputation$dep2_random_seed %||% 1)

    expr_mat <- as.matrix(expr_mat)
    imputed_flag <- is.na(expr_mat)
    imputed <- expr_mat

    set.seed(dep2_seed)

    if (dep2_method == "MinDet") {
        # Deterministic: replace NAs with 1st percentile of observed values per sample
        for (j in seq_len(ncol(imputed))) {
            x <- imputed[, j]
            if (all(is.na(x))) stop("Imputation failed: sample '", colnames(imputed)[j], "' is all-NA.")
            obs <- x[!is.na(x)]
            q01 <- quantile(obs, probs = 0.01, na.rm = TRUE)
            x[is.na(x)] <- q01
            imputed[, j] <- x
        }
    } else if (dep2_method == "MinProb") {
        # Stochastic: draw from Gaussian at the 1st percentile, narrow width
        for (j in seq_len(ncol(imputed))) {
            x <- imputed[, j]
            if (all(is.na(x))) stop("Imputation failed: sample '", colnames(imputed)[j], "' is all-NA.")
            obs <- x[!is.na(x)]
            q01 <- quantile(obs, probs = 0.01, na.rm = TRUE)
            s <- sd(obs)
            if (!is.finite(s) || s == 0) s <- 1e-8
            n_miss <- sum(is.na(x))
            if (n_miss > 0) {
                x[is.na(x)] <- rnorm(n_miss, mean = q01, sd = s * 0.3)
            }
            imputed[, j] <- x
        }
    } else {
        stop(sprintf("Unknown dep2_method: '%s'. Supported: 'MinDet', 'MinProb'.", dep2_method))
    }

    if (return_flags) list(imputed = imputed, imputed_flag = imputed_flag) else imputed
}

#' Impute proteomics expr_mat using QRILC method
#'
#' Quantile Regression Imputation of Left-Censored data (QRILC).
#' Uses the imputeLCMD package implementation when available (matching DEP behavior),
#' falling back to a faithful reimplementation of the same algorithm.
#' The algorithm uses quantile regression on a Q-Q plot to estimate the
#' complete data distribution (CDD) parameters, then samples from a truncated
#' normal with upper bound at the pNAs-th quantile of the CDD.
#'
#' @param expr_mat numeric matrix (proteins x samples), may contain NAs
#' @param cfg      mode-level config (config$modes$proteomics)
#' @param return_flags if TRUE, return list(imputed, imputed_flag)
#' @return imputed matrix or list(imputed, imputed_flag)
impute_proteomics_qrilc <- function(expr_mat, cfg, return_flags = FALSE) {
    qrilc_seed <- as.integer(cfg$imputation$qrilc_random_seed %||% 1)

    expr_mat <- as.matrix(expr_mat)
    imputed_flag <- is.na(expr_mat)

    set.seed(qrilc_seed)

    # Use imputeLCMD::impute.QRILC when available (exact DEP match)
    if (requireNamespace("imputeLCMD", quietly = TRUE)) {
        result <- imputeLCMD::impute.QRILC(expr_mat, tune.sigma = 1)
        imputed <- as.matrix(result[[1]])
        rownames(imputed) <- rownames(expr_mat)
        colnames(imputed) <- colnames(expr_mat)
    } else {
        # Fallback: faithful reimplementation of imputeLCMD::impute.QRILC algorithm
        message("imputeLCMD package not available; using built-in QRILC (quantile regression).")
        imputed <- qrilc_impute_builtin(expr_mat)
    }

    if (return_flags) list(imputed = imputed, imputed_flag = imputed_flag) else imputed
}

#' Built-in QRILC implementation matching imputeLCMD algorithm
#'
#' Reimplements the quantile regression approach from imputeLCMD::impute.QRILC
#' for environments where the package is not available.
#' @param expr_mat numeric matrix with NAs
#' @return imputed numeric matrix
qrilc_impute_builtin <- function(expr_mat) {
    nFeatures <- nrow(expr_mat)
    nSamples <- ncol(expr_mat)
    imputed <- expr_mat

    for (j in seq_len(nSamples)) {
        curr_sample <- expr_mat[, j]
        if (all(is.na(curr_sample))) stop("Imputation failed: sample '", colnames(expr_mat)[j], "' is all-NA.")

        n_miss <- sum(is.na(curr_sample))
        if (n_miss == 0) next

        # Fraction of missing values
        pNAs <- n_miss / length(curr_sample)
        upper_q <- 0.99

        # Build theoretical normal quantiles at observed probability points
        q_normal <- stats::qnorm(
            seq(pNAs + 0.001, upper_q + 0.001, (upper_q - pNAs) / (upper_q * 100)),
            mean = 0, sd = 1
        )

        # Build empirical quantiles of observed values
        q_curr <- stats::quantile(
            curr_sample,
            probs = seq(0.001, upper_q + 0.001, 0.01),
            na.rm = TRUE
        )

        # Quantile regression: empirical ~ theoretical => intercept = mean_CDD, slope = sd_CDD
        temp_qr <- stats::lm(q_curr ~ q_normal)
        mean_cdd <- as.numeric(temp_qr$coefficients[1])
        sd_cdd <- as.numeric(temp_qr$coefficients[2])
        if (!is.finite(sd_cdd) || sd_cdd <= 0) sd_cdd <- 1e-8

        # Upper bound: pNAs-th quantile of the complete data distribution
        upper_bound <- stats::qnorm(pNAs + 0.001, mean = mean_cdd, sd = sd_cdd)

        # Sample from truncated normal via inverse CDF
        u <- stats::runif(nFeatures, min = 0, max = stats::pnorm(upper_bound, mean = mean_cdd, sd = sd_cdd))
        data_to_imp <- stats::qnorm(u, mean = mean_cdd, sd = sd_cdd)

        imputed[is.na(curr_sample), j] <- data_to_imp[is.na(curr_sample)]
    }

    imputed
}

#' Impute proteomics expr_mat using MinVal method (floor of minimum observed value)
#'
#' Deterministic imputation: replaces each missing value with the floor of the
#' minimum observed value in that sample column. Simple LOD-style approach.
#'
#' @param expr_mat numeric matrix (proteins x samples), may contain NAs
#' @param cfg      mode-level config (config$modes$proteomics)
#' @param return_flags if TRUE, return list(imputed, imputed_flag)
#' @return imputed matrix or list(imputed, imputed_flag)
impute_proteomics_minval <- function(expr_mat, cfg, return_flags = FALSE) {
    expr_mat <- as.matrix(expr_mat)
    imputed_flag <- is.na(expr_mat)
    imputed <- expr_mat

    for (j in seq_len(ncol(imputed))) {
        x <- imputed[, j]
        if (all(is.na(x))) stop("Imputation failed: sample '", colnames(imputed)[j], "' is all-NA.")
        obs <- x[!is.na(x)]
        fill_val <- floor(min(obs))
        x[is.na(x)] <- fill_val
        imputed[, j] <- x
    }

    if (return_flags) list(imputed = imputed, imputed_flag = imputed_flag) else imputed
}

#' Wrapper to run multiple imputations (with seed increments)
make_imputations_proteomics <- function(expr_mat, cfg, verbose = FALSE) {
    imp_cfg <- cfg$modes$proteomics$imputation
    method <- imp_cfg$method %||% "perseus_like"
    if (method == "perseus") method <- "perseus_like"
    multi_imp <- imp_cfg$multi_imputation %||% TRUE
    n_imputations <- as.integer(imp_cfg$no_repetitions)
    seed_base <- cfg$params$seed

    if (isFALSE(multi_imp)) {
        n_imputations <- 1L
        if (isTRUE(verbose)) message("Single imputation mode (multi_imputation: false)")
    }

    stopifnot(is.matrix(expr_mat))

    # For deterministic methods (none, MinDet, minval), all runs are identical —
    # still produce n_imputations copies for pipeline compatibility
    is_deterministic <- method == "none" ||
        method == "minval" ||
        (method == "dep2" && (imp_cfg$dep2_method %||% "MinDet") == "MinDet")

    imps <- vector("list", n_imputations)
    for (i in seq_len(n_imputations)) {
        if (isTRUE(verbose)) message(sprintf("Imputation [%s]: %d / %d", method, i, n_imputations))
        set.seed(as.integer(seed_base) + i)

        expr_imp_i <- impute_proteomics(
            expr_mat,
            cfg = cfg$modes$proteomics,
            return_flags = FALSE
        )
        imps[[i]] <- expr_imp_i
    }
    imps
}

#' Perseus-like imputation (downshifted & narrowed normal distribution)
perseus_like_impute_with_flags <- function(expr_mat, width = 0.3, downshift = 1.8) {
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
