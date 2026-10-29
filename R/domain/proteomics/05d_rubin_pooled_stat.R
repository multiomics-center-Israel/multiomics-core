#' Rubin's-rules pooling of per-imputation limma fits
#'
#' The proteomics limma path fits one moderated linear model per imputation
#' draw (run_limma_multimp()). The DE summary combines those draws by voting and
#' quantiles; this file combines them into one per-protein statistic instead, for
#' ranking proteins in the pooled-ranking enrichment. It never reads or changes
#' the DE summary.
#'
#' The pipeline method pools the eBayes-moderated variances ("option A"), as
#' RNAseqCovarImpute does (Baker et al. 2024, Genome Biol 25:236): per-draw
#' coefficients are averaged, and the within-draw moderated variance is combined
#' with (1 + 1/m) times the between-draw variance (Rubin 1987), with
#' Barnard-Rubin degrees of freedom (Barnard & Rubin 1999, Biometrika 86:948).
#' Every choice of variance and complete-data df is made in one place,
#' rubin_pool_estimates(), so a variant (such as the mi4p-style "option B" in
#' scripts/validate_pooled_enrichment.R) only changes what it passes in.

#' Complete-data degrees of freedom rule used by pool_limma_imputations()
#'
#' Baker et al. (2024) take the average of the per-draw moderated df as the
#' complete-data df in the Barnard-Rubin formula; this follows them. Kept as a
#' constant so the rule recorded in the target settings is the rule applied.
#'
#' @return Character scalar naming the rule.
pooled_df_com_rule <- function() "mean_df_total"

#' Combine per-draw estimates and variances with Rubin's rules
#'
#' @param Q Numeric matrix (features x m) of per-draw estimates.
#' @param U Numeric matrix (features x m) of per-draw squared standard errors.
#' @param df_com Numeric vector (one per feature) of complete-data degrees of
#'   freedom for the Barnard-Rubin adjustment.
#' @return data.frame with one row per feature: \code{estimate} (mean of
#'   \code{Q}), \code{W} (mean within-draw variance), \code{B} (between-draw
#'   variance), \code{T} (\code{W + (1 + 1/m) B}), \code{t}, \code{df}
#'   (Barnard-Rubin), \code{riv} (relative increase in variance),
#'   \code{lambda}, \code{fmi} (fraction of missing information) and
#'   \code{pooled}. Rows with one draw, or with a between-draw variance that is
#'   numerically zero, are not pooled: \code{t} is the single-fit statistic on
#'   \code{df_com}, and \code{riv}, \code{lambda} and \code{fmi} are 0.
rubin_pool_estimates <- function(Q, U, df_com) {
    Q <- as.matrix(Q)
    U <- as.matrix(U)
    if (!identical(dim(Q), dim(U))) {
        stop("rubin_pool_estimates(): Q and U must have the same dimensions ",
             "(got ", paste(dim(Q), collapse = "x"), " and ",
             paste(dim(U), collapse = "x"), ").", call. = FALSE)
    }
    if (length(df_com) != nrow(Q)) {
        stop("rubin_pool_estimates(): df_com needs one value per row of Q.", call. = FALSE)
    }
    m <- ncol(Q)

    estimate <- rowMeans(Q)
    W <- rowMeans(U)
    B <- if (m > 1L) apply(Q, 1, stats::var) else rep(0, nrow(Q))
    # Identical draws (deterministic imputation, or a protein with nothing
    # imputed) give a between-draw variance that is zero up to rounding. Treat
    # it as zero so those rows reproduce the single-fit statistic exactly.
    pooled <- is.finite(B) & is.finite(W) & B > .Machine$double.eps * W
    T_total <- W + (1 + 1 / m) * B
    T_total[!pooled] <- W[!pooled]

    lambda <- ifelse(pooled, (1 + 1 / m) * B / T_total, 0)
    riv    <- ifelse(pooled, (1 + 1 / m) * B / W, 0)
    df_old <- (m - 1) / lambda^2
    df_obs <- (df_com + 1) / (df_com + 3) * df_com * (1 - lambda)
    df_br  <- ifelse(pooled, 1 / (1 / df_old + 1 / df_obs), df_com)
    fmi    <- ifelse(pooled, (riv + 2 / (df_br + 3)) / (1 + riv), 0)

    data.frame(
        estimate = estimate,
        W        = W,
        B        = B,
        T        = T_total,
        t        = estimate / sqrt(T_total),
        df       = df_br,
        riv      = riv,
        lambda   = lambda,
        fmi      = fmi,
        pooled   = pooled,
        row.names = rownames(Q),
        stringsAsFactors = FALSE
    )
}

#' Pool per-imputation limma fits for one contrast
#'
#' Option A: per draw, the squared standard error is the eBayes posterior
#' variance times the squared unscaled standard deviation
#' (\code{s2.post * stdev.unscaled^2}), so the moderated t of a single draw is
#' reproduced exactly when all draws agree. fdrtool, when enabled, rewrites only
#' \code{fit2$p.value}; this reads coefficients and variances, so it is
#' unaffected.
#'
#' @param runs List of per-imputation results from run_limma_multimp(); each
#'   needs \code{fit2} (an eBayes MArrayLM). Per-contrast limma, t-test and
#'   ANOVA runs carry no \code{fit2} and are rejected.
#' @param contrast Contrast name (a column of \code{fit2$coefficients}).
#' @param observed Optional logical matrix (proteins x samples), TRUE where the
#'   pre-imputation value was observed. Drives \code{n_observed} and
#'   \code{frac_imputed}; NULL leaves them NA.
#' @return data.frame, one row per protein in fit order: \code{FeatureID},
#'   \code{contrast}, \code{m}, \code{pooled_logFC}, \code{W}, \code{B},
#'   \code{T}, \code{se_pooled}, \code{t_pooled}, \code{df_com}, \code{df_BR},
#'   \code{riv}, \code{fmi}, \code{z_pooled}, \code{AveExpr},
#'   \code{n_observed_numerator}, \code{n_observed_denominator},
#'   \code{n_observed}, \code{n_samples}, \code{frac_imputed}, \code{pooled}.
pool_limma_imputations <- function(runs, contrast, observed = NULL) {
    if (length(runs) == 0) stop("pool_limma_imputations(): no imputation runs given.", call. = FALSE)
    fits <- lapply(runs, function(r) r$fit2)
    if (any(vapply(fits, is.null, logical(1)))) {
        stop("pool_limma_imputations(): the runs carry no limma fit (fit2). Pooling ",
             "needs de.method 'limma'; 'limma_percontrast', 'ttest', 'welch' and ",
             "'anova' fit per-contrast or non-limma models and are not supported.",
             call. = FALSE)
    }
    ref_ids <- rownames(fits[[1]]$coefficients)
    for (i in seq_along(fits)) {
        if (!contrast %in% colnames(fits[[i]]$coefficients)) {
            stop(sprintf("pool_limma_imputations(): contrast '%s' is not in fit %d (has: %s).",
                         contrast, i, paste(colnames(fits[[i]]$coefficients), collapse = ", ")),
                 call. = FALSE)
        }
        if (!identical(rownames(fits[[i]]$coefficients), ref_ids)) {
            stop(sprintf(paste0("pool_limma_imputations(): fit %d has different protein rows ",
                                "from fit 1. Every imputation must be fitted on the same ",
                                "proteins in the same order."), i), call. = FALSE)
        }
    }

    Q  <- vapply(fits, function(f) f$coefficients[, contrast], numeric(length(ref_ids)))
    U  <- vapply(fits, function(f) f$s2.post * f$stdev.unscaled[, contrast]^2,
                 numeric(length(ref_ids)))
    df_total <- vapply(fits, function(f) as.numeric(f$df.total), numeric(length(ref_ids)))
    Q <- matrix(Q, nrow = length(ref_ids), dimnames = list(ref_ids, NULL))
    U <- matrix(U, nrow = length(ref_ids), dimnames = list(ref_ids, NULL))
    df_total <- matrix(df_total, nrow = length(ref_ids))
    df_com <- switch(pooled_df_com_rule(),
                     mean_df_total = rowMeans(df_total),
                     min_df_total  = apply(df_total, 1, min))

    rp <- rubin_pool_estimates(Q, U, df_com)
    ave <- rowMeans(matrix(vapply(fits, function(f) as.numeric(f$Amean),
                                  numeric(length(ref_ids))), nrow = length(ref_ids)))

    counts <- count_observed_in_contrast(runs[[1]], contrast, observed, ref_ids)

    data.frame(
        FeatureID    = ref_ids,
        contrast     = contrast,
        m            = length(fits),
        pooled_logFC = rp$estimate,
        W            = rp$W,
        B            = rp$B,
        T            = rp$T,
        se_pooled    = sqrt(rp$T),
        t_pooled     = rp$t,
        df_com       = df_com,
        df_BR        = rp$df,
        riv          = rp$riv,
        fmi          = rp$fmi,
        z_pooled     = limma::zscoreT(rp$t, df = rp$df),
        AveExpr      = ave,
        counts,
        pooled       = rp$pooled,
        row.names    = NULL,
        stringsAsFactors = FALSE
    )
}

#' Count observed (pre-imputation) values in the two sides of a contrast
#'
#' Samples on the numerator side are those in groups with a positive contrast
#' coefficient, the denominator side those with a negative one.
#'
#' @param run One element of run_limma_multimp() (needs \code{design},
#'   \code{contrast_matrix} and \code{expr_imp}).
#' @param contrast Contrast name.
#' @param observed Logical matrix (proteins x samples) or NULL.
#' @param ids Protein IDs in fit order.
#' @return data.frame with \code{n_observed_numerator},
#'   \code{n_observed_denominator}, \code{n_observed}, \code{n_samples} and
#'   \code{frac_imputed}, one row per protein.
count_observed_in_contrast <- function(run, contrast, observed, ids) {
    na_out <- data.frame(n_observed_numerator = NA_integer_, n_observed_denominator = NA_integer_,
                         n_observed = NA_integer_, n_samples = NA_integer_,
                         frac_imputed = NA_real_)[rep(1L, length(ids)), , drop = FALSE]
    rownames(na_out) <- NULL
    if (is.null(observed) || is.null(run$design) || is.null(run$contrast_matrix)) return(na_out)

    cvec <- run$contrast_matrix[, contrast]
    samples <- colnames(run$expr_imp)
    side <- function(lvls) {
        if (length(lvls) == 0) return(character(0))
        samples[rowSums(run$design[, lvls, drop = FALSE]) > 0]
    }
    num_s <- side(names(cvec)[cvec > 0])
    den_s <- side(names(cvec)[cvec < 0])

    missing_rows <- setdiff(ids, rownames(observed))
    missing_cols <- setdiff(c(num_s, den_s), colnames(observed))
    if (length(missing_rows) > 0 || length(missing_cols) > 0) {
        stop(sprintf(paste0("count_observed_in_contrast(): the observed-value mask does not ",
                            "cover the fitted data (%d proteins, %d samples missing)."),
                     length(missing_rows), length(missing_cols)), call. = FALSE)
    }
    obs <- observed[ids, , drop = FALSE]
    n_num <- as.integer(rowSums(obs[, num_s, drop = FALSE]))
    n_den <- as.integer(rowSums(obs[, den_s, drop = FALSE]))
    n_tot <- length(num_s) + length(den_s)
    data.frame(
        n_observed_numerator   = n_num,
        n_observed_denominator = n_den,
        n_observed             = n_num + n_den,
        n_samples              = n_tot,
        frac_imputed           = 1 - (n_num + n_den) / n_tot,
        stringsAsFactors = FALSE
    )
}
