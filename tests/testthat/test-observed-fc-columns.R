# Tests for compute_observed_fc_columns(): the measured-only fold change that
# sits beside the pipeline's imputation-based linearFC in the final table.

make_meta <- function() {
    data.frame(
        SampleID = paste0("S", 1:6),
        Group    = rep(c("A", "B"), each = 3),
        stringsAsFactors = FALSE
    )
}

make_contrasts <- function(name = "A_vs_B") {
    data.frame(
        Contrast_name = name,
        Factor        = "Group",
        Numerator     = "A",
        Denominator   = "B",
        stringsAsFactors = FALSE
    )
}

make_expr <- function(values, features = "F1") {
    matrix(values,
           nrow = length(features), byrow = TRUE,
           dimnames = list(features, paste0("S", 1:6)))
}

test_that("observed log2FC is the difference of complete group means", {
    expr <- make_expr(c(10, 12, 14, 4, 6, 8))  # A mean 12, B mean 6
    res <- compute_observed_fc_columns(expr, make_meta(), "SampleID", make_contrasts())

    expect_equal(res[["obs.log2FC.A_vs_B"]], 6)
    expect_equal(res[["obs.linearFC.A_vs_B"]], 64)
    expect_equal(res[["n_obs.A"]], 3L)
    expect_equal(res[["n_obs.B"]], 3L)
})

test_that("missing values are dropped from the group mean, not treated as zero", {
    # A has one missing value; the mean must use the two observed values only.
    expr <- make_expr(c(10, 12, NA, 4, 6, 8))
    res <- compute_observed_fc_columns(expr, make_meta(), "SampleID", make_contrasts())

    expect_equal(res[["obs.log2FC.A_vs_B"]], 11 - 6)
    expect_equal(res[["n_obs.A"]], 2L)
    expect_equal(res[["n_obs.B"]], 3L)
})

test_that("a group with nothing observed yields NA, not a one-armed difference", {
    expr <- make_expr(c(NA, NA, NA, 4, 6, 8))
    res <- compute_observed_fc_columns(expr, make_meta(), "SampleID", make_contrasts())

    expect_true(is.na(res[["obs.log2FC.A_vs_B"]]))
    expect_true(is.na(res[["obs.linearFC.A_vs_B"]]))
    expect_equal(res[["n_obs.A"]], 0L)
})

test_that("linearFC is signed: ratios below 1 are written as -1/ratio", {
    expr <- make_expr(c(4, 4, 4, 6, 6, 6))  # log2FC = -2, ratio = 0.25
    res <- compute_observed_fc_columns(expr, make_meta(), "SampleID", make_contrasts())

    expect_equal(res[["obs.log2FC.A_vs_B"]], -2)
    expect_equal(res[["obs.linearFC.A_vs_B"]], -4)
})

test_that("contrast names are normalized for proteomics but left alone elsewhere", {
    expr <- make_expr(c(10, 12, 14, 4, 6, 8))
    ctr  <- make_contrasts("A vs B")

    prot <- compute_observed_fc_columns(expr, make_meta(), "SampleID", ctr,
                                        mode = "proteomics")
    metab <- compute_observed_fc_columns(expr, make_meta(), "SampleID", ctr,
                                         mode = "metabolomics")

    expect_true(paste0("obs.log2FC.", normalize_contrast_name("A vs B")) %in% names(prot))
    expect_true("obs.log2FC.A vs B" %in% names(metab))
})

test_that("multiple features keep their identity and row order", {
    expr <- matrix(
        c(10, 12, 14, 4, 6, 8,
          1, 1, 1, 1, 1, 1),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("F1", "F2"), paste0("S", 1:6))
    )
    res <- compute_observed_fc_columns(expr, make_meta(), "SampleID", make_contrasts())

    expect_equal(rownames(res), c("F1", "F2"))
    expect_equal(res[["obs.log2FC.A_vs_B"]], c(6, 0))
})

test_that("insufficient inputs return NULL rather than a malformed table", {
    expr <- make_expr(c(10, 12, 14, 4, 6, 8))

    expect_null(compute_observed_fc_columns(NULL, make_meta(), "SampleID", make_contrasts()))
    expect_null(compute_observed_fc_columns(expr, make_meta(), "NotAColumn", make_contrasts()))
    expect_null(compute_observed_fc_columns(expr, make_meta(), "SampleID",
                                            data.frame(x = 1)))
})

test_that("an unresolvable grouping column warns and skips instead of guessing", {
    expr <- make_expr(c(10, 12, 14, 4, 6, 8))
    ctr  <- make_contrasts()
    ctr$Factor <- "NoSuchColumn"

    expect_warning(
        res <- compute_observed_fc_columns(expr, make_meta(), "SampleID", ctr),
        "could not resolve a single grouping column"
    )
    expect_null(res)
})
