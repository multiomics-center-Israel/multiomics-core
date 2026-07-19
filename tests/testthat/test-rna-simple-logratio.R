# tests/testthat/test-rna-simple-logratio.R
#
# Unit tests for the GLM-free simple log2 fold change helper used by
# run_deseq2_de() to add the `simple_log2FC` column:
#   - compute_simple_log2fc()   (R/domain/rnaseq/04_de_summary.R)
#
# Fully deterministic: values are fixed normalized counts, no RNG.

make_norm_counts <- function() {
    m <- rbind(
        up_in_test   = c(10, 20,  4,  6),   # higher in the test group
        down_in_test = c( 0,  0, 100, 100),  # near-zero in test, high in control
        flat         = c(50, 50, 50, 50)     # no difference
    )
    colnames(m) <- c("t1", "t2", "c1", "c2")
    m
}

test_that("compute_simple_log2fc matches log2((mean_num + 1)/(mean_den + 1))", {
    m <- make_norm_counts()
    grp <- c("Test", "Test", "Ctrl", "Ctrl")

    out <- compute_simple_log2fc(m, grp, num = "Test", den = "Ctrl")

    expect_equal(out[["up_in_test"]],   log2((15 + 1) / (5 + 1)))
    expect_equal(out[["down_in_test"]], log2((0 + 1) / (100 + 1)))
    expect_equal(out[["flat"]],         0)
})

test_that("orientation is num - den (positive means higher in num)", {
    m <- make_norm_counts()
    grp <- c("Test", "Test", "Ctrl", "Ctrl")

    out <- compute_simple_log2fc(m, grp, num = "Test", den = "Ctrl")

    expect_gt(out[["up_in_test"]], 0)
    expect_lt(out[["down_in_test"]], 0)
    # Swapping num/den flips the sign.
    swapped <- compute_simple_log2fc(m, grp, num = "Ctrl", den = "Test")
    expect_equal(swapped[["up_in_test"]], -out[["up_in_test"]])
})

test_that("names are preserved and length equals number of genes", {
    m <- make_norm_counts()
    grp <- c("Test", "Test", "Ctrl", "Ctrl")

    out <- compute_simple_log2fc(m, grp, num = "Test", den = "Ctrl")

    expect_named(out, rownames(m))
    expect_length(out, nrow(m))
})

test_that("an absent group yields a named all-NA vector rather than an error", {
    m <- make_norm_counts()
    grp <- c("Test", "Test", "Ctrl", "Ctrl")

    out <- compute_simple_log2fc(m, grp, num = "Missing", den = "Ctrl")

    expect_true(all(is.na(out)))
    expect_named(out, rownames(m))
})

test_that("pseudocount is configurable", {
    m <- make_norm_counts()
    grp <- c("Test", "Test", "Ctrl", "Ctrl")

    out <- compute_simple_log2fc(m, grp, num = "Test", den = "Ctrl", pseudocount = 0.5)
    expect_equal(out[["up_in_test"]], log2((15 + 0.5) / (5 + 0.5)))
})
