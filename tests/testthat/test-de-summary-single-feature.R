# tests/testthat/test-de-summary-single-feature.R
#
# summarize_limma_mult_imputation() on a dataset with exactly one feature.
#
# The per-run blocks were built with sapply(), which simplifies to a plain
# vector whenever each run returns a single value -- that is, whenever there is
# one feature. rowSums() on the pass matrix then aborted with "'x' must be an
# array of at least two dimensions", so the whole summary failed. One surviving
# feature after missingness filtering is enough to reach it.
#
# These tests drive the production summariser rather than a helper, so they
# fail against the sapply() version and pass against the reshaped one.
#
# All fixtures are synthetic (p1..p3, contrast S_vs_NS).

sf_config <- function(n_reps = 3, min_passed = 2, multi = TRUE) {
    list(modes = list(proteomics = list(
        de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5, use_adj_for_pass1 = TRUE),
        imputation = list(multi_imputation = multi, no_repetitions = n_reps,
                          min_no_passed = min_passed),
        de_table = list(id_col = "FeatureID")
    )))
}

# One limma-shaped table per run. `logfc`, `pval` and `padj` are given per run,
# so a test can state exactly what each run reported.
sf_runs <- function(logfc, pval = NULL, padj = NULL,
                    features = "p1", contrast = "S_vs_NS") {
    logfc <- as.matrix(logfc)                       # features x runs
    n_feat <- nrow(logfc); n_runs <- ncol(logfc)
    pval <- if (is.null(pval)) matrix(1e-4, n_feat, n_runs) else as.matrix(pval)
    padj <- if (is.null(padj)) matrix(1e-3, n_feat, n_runs) else as.matrix(padj)

    lapply(seq_len(n_runs), function(i) {
        per_contrast <- list(data.frame(
            FeatureID = features,
            logFC     = logfc[, i],
            P.Value   = pval[, i],
            adj.P.Val = padj[, i],
            stringsAsFactors = FALSE
        ))
        names(per_contrast) <- contrast
        per_contrast
    })
}


# =============================================================================
# One feature
# =============================================================================

test_that("one feature with multi-imputation summarises without error", {
    # The regression case. Under sapply() this aborted in rowSums().
    cfg <- sf_config(n_reps = 3, min_passed = 2)
    runs <- sf_runs(logfc = matrix(c(2.0, 2.2, 2.4), nrow = 1))

    expect_no_error(sdf <- summarize_limma_mult_imputation(runs, cfg))
    expect_equal(nrow(sdf), 1L)
    expect_equal(sdf$FeatureID, "p1")
})

test_that("one feature pools exactly as the documented rule says", {
    cfg <- sf_config(n_reps = 3, min_passed = 2)
    lfc <- c(2.0, 2.2, 2.4)
    sdf <- summarize_limma_mult_imputation(sf_runs(matrix(lfc, nrow = 1)), cfg)

    # log2( mean( 2^logFC ) ), written a second way so the production formula
    # cannot make this agree with itself by construction.
    expected_ratio <- mean(exp(log(2) * lfc))
    expect_equal(sdf$linearRatio.imputs.S_vs_NS, expected_ratio)
    expect_equal(sdf$log2FC.imputs.S_vs_NS, log2(expected_ratio))
    expect_equal(sdf$linearFC.imputs.S_vs_NS, signif(expected_ratio, 3))
})

test_that("one feature summarises p-values by the same quantile rule", {
    cfg <- sf_config(n_reps = 3, min_passed = 2)
    pv <- c(1e-5, 1e-4, 1e-3)
    pa <- c(1e-4, 1e-3, 1e-2)
    sdf <- summarize_limma_mult_imputation(
        sf_runs(matrix(c(2, 2, 2), nrow = 1),
                pval = matrix(pv, nrow = 1), padj = matrix(pa, nrow = 1)), cfg)

    q <- 2 / 3   # min_no_passed / no_repetitions
    expect_equal(unname(sdf$pvalue.imputs.S_vs_NS), unname(quantile(pv, probs = q)))
    expect_equal(unname(sdf$padj.imputs.S_vs_NS), unname(quantile(pa, probs = q)))
})

test_that("one feature votes with the same min_no_passed rule", {
    # Two of three runs clear both cutoffs, min_no_passed is 2, and the summary
    # padj must also clear p_cutoff for the consensus call to stand.
    cfg <- sf_config(n_reps = 3, min_passed = 2)

    passing <- summarize_limma_mult_imputation(
        sf_runs(matrix(c(2, 2, 0.1), nrow = 1),
                padj = matrix(c(1e-3, 1e-3, 1e-3), nrow = 1)), cfg)
    expect_equal(passing$sum.pass.S_vs_NS, 2)
    expect_equal(passing$pass.imputs.S_vs_NS, 1)

    # One of three passes: below the vote, so no consensus call.
    failing <- summarize_limma_mult_imputation(
        sf_runs(matrix(c(2, 0.1, 0.1), nrow = 1),
                padj = matrix(c(1e-3, 1e-3, 1e-3), nrow = 1)), cfg)
    expect_equal(failing$sum.pass.S_vs_NS, 1)
    expect_true(is.na(failing$pass.imputs.S_vs_NS))
})

test_that("one feature with a single imputation reports that one run", {
    cfg <- sf_config(n_reps = 1, min_passed = 1, multi = FALSE)
    sdf <- summarize_limma_mult_imputation(sf_runs(matrix(2.5, nrow = 1)), cfg)

    expect_equal(nrow(sdf), 1L)
    # One run: the consensus is that run's own value, with no pooling to do.
    expect_equal(sdf$log2FC.imputs.S_vs_NS, 2.5)
    expect_equal(sdf$linearRatio.imputs.S_vs_NS, 2^2.5)
    expect_equal(unname(sdf$pvalue.imputs.S_vs_NS), 1e-4)
    expect_equal(sdf$sum.pass.S_vs_NS, 1)
    expect_equal(sdf$pass.imputs.S_vs_NS, 1)
})

test_that("one feature and several contrasts keeps one row per feature", {
    cfg <- sf_config(n_reps = 2, min_passed = 1)
    runs <- lapply(1:2, function(i) {
        one <- function(lfc) data.frame(FeatureID = "p1", logFC = lfc,
                                        P.Value = 1e-4, adj.P.Val = 1e-3,
                                        stringsAsFactors = FALSE)
        list(S_vs_NS = one(2 + 0.1 * i), S_vs_REF = one(-1 - 0.1 * i))
    })

    sdf <- summarize_limma_mult_imputation(runs, cfg)
    expect_equal(nrow(sdf), 1L)
    expect_true(all(c("log2FC.imputs.S_vs_NS", "log2FC.imputs.S_vs_REF") %in% names(sdf)))
    expect_gt(sdf$log2FC.imputs.S_vs_NS, 0)
    expect_lt(sdf$log2FC.imputs.S_vs_REF, 0)
})


# =============================================================================
# The multi-feature path is unchanged
# =============================================================================

test_that("several features summarise exactly as before", {
    cfg <- sf_config(n_reps = 3, min_passed = 2)
    lfc <- matrix(c(2.0, 2.2, 2.4,      # p1
                   -1.0, -1.1, -1.2,    # p2
                    0.05, 0.05, 0.05),  # p3, below the fold-change cutoff
                  nrow = 3, byrow = TRUE)
    sdf <- summarize_limma_mult_imputation(
        sf_runs(lfc, features = c("p1", "p2", "p3")), cfg)

    expect_equal(nrow(sdf), 3L)
    expect_equal(sdf$FeatureID, c("p1", "p2", "p3"))
    expect_equal(sdf$linearRatio.imputs.S_vs_NS, rowMeans(2^lfc))
    expect_equal(sdf$log2FC.imputs.S_vs_NS, log2(rowMeans(2^lfc)))
    # p1 and p2 clear the cutoff in all three runs; p3 in none.
    expect_equal(sdf$sum.pass.S_vs_NS, c(3, 3, 0))
    expect_equal(sdf$pass.imputs.S_vs_NS, c(1, 1, NA))
})

test_that("several features with a single run are unchanged too", {
    # sapply() already returned a features x 1 matrix here, so this pins that
    # the reshape did not disturb the case it was not aimed at.
    cfg <- sf_config(n_reps = 1, min_passed = 1, multi = FALSE)
    lfc <- matrix(c(2.0, -1.0), nrow = 2)
    sdf <- summarize_limma_mult_imputation(
        sf_runs(lfc, features = c("p1", "p2")), cfg)

    expect_equal(nrow(sdf), 2L)
    expect_equal(sdf$log2FC.imputs.S_vs_NS, as.numeric(lfc))
})


# =============================================================================
# The reshape itself: same numbers, whatever the shape
# =============================================================================

test_that("a one-feature run gives the same answer as that row of a larger one", {
    # The strongest statement of the contract: summarising a feature alone and
    # summarising it alongside others must produce identical numbers. Only the
    # sapply() simplification could have made those differ.
    cfg <- sf_config(n_reps = 3, min_passed = 2)
    lfc_all <- matrix(c(2.0, 2.2, 2.4,
                        -1.0, -1.1, -1.2),
                      nrow = 2, byrow = TRUE)

    many <- summarize_limma_mult_imputation(
        sf_runs(lfc_all, features = c("p1", "p2")), cfg)
    alone <- summarize_limma_mult_imputation(
        sf_runs(lfc_all[1, , drop = FALSE], features = "p1"), cfg)

    for (col in c("log2FC.imputs.S_vs_NS", "linearRatio.imputs.S_vs_NS",
                  "linearFC.imputs.S_vs_NS", "pvalue.imputs.S_vs_NS",
                  "padj.imputs.S_vs_NS", "sum.pass.S_vs_NS")) {
        expect_equal(unname(alone[[col]]), unname(many[[col]][1]), info = col)
    }
})
