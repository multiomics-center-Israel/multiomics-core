# tests/testthat/test-export-raw-group-stats.R
#
# Tests for the pre-imputation per-group export columns (R/core/05_export_excel.R).
#
# The results table already carried Mean.<group> and log2FC_from_means.<contrast>
# computed on the matrix the model was fitted on (imputed for proteomics,
# library-size-normalized for RNA-seq). Those answer "what did the model see".
# Mean.raw.<group>, N.observed.<group> and log2FC_from_raw.<contrast> answer
# "what was actually measured", which is what a reader needs to check a fold
# change by hand.

mk_meta <- function() {
    data.frame(SampleID = c("A1", "A2", "A3", "B1", "B2", "B3"),
               Group = rep(c("trt", "ctl"), each = 3),
               stringsAsFactors = FALSE)
}

mk_contrasts <- function() {
    data.frame(Contrast_name = "trt_vs_ctl", Factor = "Group",
               Numerator = "trt", Denominator = "ctl",
               stringsAsFactors = FALSE)
}

mk_expr <- function(na_positions = NULL) {
    m <- matrix(c(10, 10, 10,  8,  8,  8,
                  20, 20, 20, 20, 20, 20),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("F1", "F2"), c("A1", "A2", "A3", "B1", "B2", "B3")))
    for (p in na_positions) m[p[1], p[2]] <- NA
    m
}

test_that("Mean.raw columns are emitted under the requested prefix", {
    out <- compute_group_mean_columns(mk_expr(), mk_meta(), "SampleID",
                                      mk_contrasts(), prefix = "Mean.raw.")
    expect_setequal(colnames(out), c("Mean.raw.trt", "Mean.raw.ctl"))
    expect_equal(out[["Mean.raw.trt"]], c(10, 20))
    expect_equal(out[["Mean.raw.ctl"]], c(8, 20))
})

test_that("the default prefix is unchanged, so existing callers keep Mean.", {
    out <- compute_group_mean_columns(mk_expr(), mk_meta(), "SampleID", mk_contrasts())
    expect_setequal(colnames(out), c("Mean.trt", "Mean.ctl"))
})

test_that("raw means ignore NAs rather than propagating them", {
    e <- mk_expr(list(c(1, 1)))   # F1 / A1 unobserved
    out <- compute_group_mean_columns(e, mk_meta(), "SampleID",
                                      mk_contrasts(), prefix = "Mean.raw.")
    expect_equal(out[["Mean.raw.trt"]][1], 10)  # mean of the two observed
    expect_false(is.na(out[["Mean.raw.trt"]][1]))
})

test_that("N.observed counts measured values per group", {
    e <- mk_expr(list(c(1, 1), c(1, 2)))  # F1 loses two of three trt samples
    out <- compute_group_observed_columns(e, mk_meta(), "SampleID", mk_contrasts())
    expect_setequal(colnames(out), c("N.observed.trt", "N.observed.ctl"))
    expect_equal(out[["N.observed.trt"]], c(1, 3))
    expect_equal(out[["N.observed.ctl"]], c(3, 3))
})

test_that("N.observed is skipped entirely when nothing is missing", {
    # RNA-seq counts case: a constant column would only widen the workbook.
    expect_null(compute_group_observed_columns(mk_expr(), mk_meta(), "SampleID",
                                               mk_contrasts()))
})

test_that("log2FC_from_raw reads the Mean.raw prefix on a log2 scale", {
    means <- compute_group_mean_columns(mk_expr(), mk_meta(), "SampleID",
                                        mk_contrasts(), prefix = "Mean.raw.")
    fc <- compute_naive_log2fc_columns(means, mk_contrasts(), scale = "log2",
                                       prefix = "Mean.raw.")
    expect_equal(fc[["trt_vs_ctl"]], c(2, 0))   # 10-8 and 20-20
})

test_that("log2FC_from_raw takes a ratio on a linear scale", {
    means <- compute_group_mean_columns(mk_expr(), mk_meta(), "SampleID",
                                        mk_contrasts(), prefix = "Mean.raw.")
    fc <- compute_naive_log2fc_columns(means, mk_contrasts(), scale = "linear",
                                       prefix = "Mean.raw.")
    expect_equal(fc[["trt_vs_ctl"]], c(log2(10 / 8), 0))
})

test_that("a mismatched prefix yields no columns rather than wrong ones", {
    means <- compute_group_mean_columns(mk_expr(), mk_meta(), "SampleID",
                                        mk_contrasts(), prefix = "Mean.raw.")
    # Asking for the default prefix against Mean.raw. columns must not silently
    # fall through to whatever happens to be there.
    expect_warning(fc <- compute_naive_log2fc_columns(means, mk_contrasts(),
                                                      scale = "log2"),
                   "no group means")
    expect_null(fc)
})

test_that("the shared mean helper still returns NaN for an all-missing group", {
    # This is the GENERIC helper's behaviour and it is unchanged: rowMeans(na.rm
    # = TRUE) over nothing is NaN. The proteomics raw path deliberately
    # overrides it (see the NA test below); no other mode can reach this case,
    # because their Mean. columns are taken on a complete matrix.
    e <- mk_expr(list(c(1, 1), c(1, 2), c(1, 3)))
    out <- compute_group_mean_columns(e, mk_meta(), "SampleID",
                                      mk_contrasts(), prefix = "Mean.raw.")
    expect_true(is.nan(out[["Mean.raw.trt"]][1]))
    n <- compute_group_observed_columns(e, mk_meta(), "SampleID", mk_contrasts())
    expect_equal(n[["N.observed.trt"]][1], 0)
})


# =============================================================================
# The proteomics raw path: NA for an unmeasured arm, and the signed linear
# presentation of log2FC_from_raw.
# =============================================================================

mk_raw_cfg <- function() {
    list(modes = list(proteomics = list(excel = list(group_cv = TRUE))))
}

mk_pre <- function(expr) {
    list(expr_filt = expr, meta = mk_meta())
}

test_that("the proteomics raw means give NA for a group with nothing measured", {
    # NaN reads as a failed calculation; NA reads as "never measured", which is
    # what actually happened. Excel renders the two differently too.
    e <- mk_expr(list(c(1, 1), c(1, 2), c(1, 3)))   # all three trt samples of F1
    out <- build_group_raw_stats_proteomics(mk_pre(e), mk_contrasts(), mk_raw_cfg())

    expect_true(is.na(out$means[["Mean.raw.trt"]][1]))
    expect_false(is.nan(out$means[["Mean.raw.trt"]][1]))
    # The measured arm is untouched, and so is the fully measured feature.
    expect_equal(out$means[["Mean.raw.ctl"]][1], 8)
    expect_equal(out$means[["Mean.raw.trt"]][2], 20)
})

test_that("the proteomics raw means still ignore NAs where something was measured", {
    e <- mk_expr(list(c(1, 1)))   # F1 loses one of three trt samples
    out <- build_group_raw_stats_proteomics(mk_pre(e), mk_contrasts(), mk_raw_cfg())
    # The mean of the two observed values, not a value dragged towards zero by
    # counting the blank as a measurement.
    expect_equal(out$means[["Mean.raw.trt"]][1], 10)
})

test_that("an unmeasured arm propagates NA into log2FC_from_raw", {
    e <- mk_expr(list(c(1, 1), c(1, 2), c(1, 3)))
    out <- build_group_raw_stats_proteomics(mk_pre(e), mk_contrasts(), mk_raw_cfg())
    fc <- compute_naive_log2fc_columns(out$means, mk_contrasts(), scale = "log2",
                                       prefix = "Mean.raw.")
    expect_true(is.na(fc[["trt_vs_ctl"]][1]))
    expect_equal(fc[["trt_vs_ctl"]][2], 0)
})
