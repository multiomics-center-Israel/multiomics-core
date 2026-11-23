# tests/testthat/test-proteomics-sample-figures.R
#
# Inputs of the DE-matrix sample figures:
#   - select_top_variable_measured()        (07g_sample_figures.R)
#
# Synthetic data only.

test_that("top-N ranks by variance among fully measured proteins only", {
    x <- rbind(
        steady   = c(10, 10.1, 9.9, 10, 10.1, 9.9),
        wide     = c(5, 15, 5, 15, 5, 15),        # largest variance, one value imputed
        middle   = c(8, 12, 8, 12, 8, 12),
        small    = c(9, 11, 9, 11, 9, 11),
        flat     = c(7, 7, 7, 7, 7, 7)
    )
    colnames(x) <- paste0("S", 1:6)
    observed <- matrix(TRUE, nrow(x), ncol(x), dimnames = dimnames(x))
    observed["wide", "S3"] <- FALSE

    expect_equal(select_top_variable_measured(x, observed, n = 2), c("middle", "small"))
    expect_equal(select_top_variable_measured(x, observed, n = 10),
                 c("middle", "small", "steady", "flat"))
    # The full set keeps matrix order and drops the protein with an imputed value.
    expect_equal(select_top_variable_measured(x, observed, n = NULL),
                 c("steady", "middle", "small", "flat"))
    expect_error(select_top_variable_measured(x, observed[, 1:5], n = 2), "same shape")
})

test_that("ties keep their input order", {
    x <- rbind(a = c(1, 3, 1, 3), b = c(1, 3, 1, 3), c = c(2, 2, 2, 2.1))
    colnames(x) <- paste0("S", 1:4)
    observed <- matrix(TRUE, 3, 4, dimnames = dimnames(x))
    expect_equal(select_top_variable_measured(x, observed, n = 2), c("a", "b"))
})
