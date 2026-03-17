# tests/testthat/test-filter-to-biological.R
#
# Unit tests for filter_to_biological() (QC/blank sample exclusion)

test_that("filter_to_biological removes QC and blank by condition_column", {
  mat <- matrix(1:24, nrow = 3, ncol = 8,
                dimnames = list(paste0("feat_", 1:3),
                                paste0("s", 1:8)))
  meta <- data.frame(
    sample_id = paste0("s", 1:8),
    treatment = c("A", "B", "A", "B", "QC", "qc", "Blank", "blank"),
    stringsAsFactors = FALSE
  )

  res <- filter_to_biological(mat, meta, "treatment", "sample_id",
                              label = "test")

  expect_equal(ncol(res$mat), 4L)
  expect_equal(nrow(res$meta), 4L)
  expect_equal(levels(res$condition), c("A", "B"))
  expect_true(all(res$meta$sample_id %in% c("s1", "s2", "s3", "s4")))
})

test_that("filter_to_biological also checks is_QC and is_blank columns", {
  mat <- matrix(1:15, nrow = 3, ncol = 5,
                dimnames = list(paste0("f", 1:3), paste0("s", 1:5)))
  meta <- data.frame(
    sample_id = paste0("s", 1:5),
    treatment = c("A", "B", "A", "B", "A"),
    is_QC     = c(FALSE, FALSE, TRUE, FALSE, FALSE),
    is_blank  = c(FALSE, FALSE, FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  res <- filter_to_biological(mat, meta, "treatment", "sample_id",
                              label = "test")

  expect_equal(ncol(res$mat), 3L)
  expect_equal(nrow(res$meta), 3L)
  expect_true(all(res$meta$sample_id %in% c("s1", "s2", "s5")))
})

test_that("filter_to_biological is no-op when all samples are biological", {
  mat <- matrix(1:9, nrow = 3, ncol = 3,
                dimnames = list(paste0("f", 1:3), paste0("s", 1:3)))
  meta <- data.frame(
    sample_id = paste0("s", 1:3),
    treatment = c("A", "B", "A"),
    stringsAsFactors = FALSE
  )

  res <- filter_to_biological(mat, meta, "treatment", "sample_id",
                              label = "test")

  expect_equal(ncol(res$mat), 3L)
  expect_equal(nrow(res$meta), 3L)
})

test_that("isTRUE_vec handles NA and non-logical inputs", {
  expect_equal(isTRUE_vec(c(TRUE, FALSE, NA)), c(TRUE, FALSE, FALSE))
  expect_equal(isTRUE_vec(c(1, 0, NA)), c(TRUE, FALSE, FALSE))
  expect_equal(isTRUE_vec(c("TRUE", "FALSE", NA)), c(TRUE, FALSE, FALSE))
})
