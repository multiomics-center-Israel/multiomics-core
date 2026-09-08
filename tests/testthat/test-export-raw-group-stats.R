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

test_that("a group whose values are all missing gives NaN, not a wrong mean", {
    e <- mk_expr(list(c(1, 1), c(1, 2), c(1, 3)))
    out <- compute_group_mean_columns(e, mk_meta(), "SampleID",
                                      mk_contrasts(), prefix = "Mean.raw.")
    expect_true(is.nan(out[["Mean.raw.trt"]][1]))
    n <- compute_group_observed_columns(e, mk_meta(), "SampleID", mk_contrasts())
    expect_equal(n[["N.observed.trt"]][1], 0)
})
