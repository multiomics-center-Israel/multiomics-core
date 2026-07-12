# tests/testthat/test-mummichog-pinned-input.R
#
# Unit tests for the pinned mummichog v2 input writer + validation (06c). These
# are pure R (no Python) and run anywhere testthat is available.

make_stats <- function(n = 3, mz = NULL) {
  if (is.null(mz)) mz <- seq(150, by = 20, length.out = n)
  data.frame(
    feature_id     = paste0("feat_", seq_len(n)),
    mz             = mz,
    retention_time = seq(10, by = 5, length.out = n),
    p_value        = seq(0.001, 0.05, length.out = n),
    statistic      = seq(-2, 2, length.out = n),
    stringsAsFactors = FALSE
  )
}

test_that("missing required columns are rejected", {
  tbl <- make_stats()[, c("feature_id", "mz", "retention_time", "p_value")]  # no statistic
  expect_error(
    write_mummichog_input(tbl, file.path(withr::local_tempdir(), "in.tsv")),
    "Missing required columns"
  )
})

test_that("non-numeric values are rejected", {
  tbl <- make_stats()
  tbl$mz <- c("abc", "150", "170")
  expect_error(
    write_mummichog_input(tbl, file.path(withr::local_tempdir(), "in.tsv")),
    "non-numeric"
  )
})

test_that("p-values outside [0, 1] are rejected", {
  tbl <- make_stats()
  tbl$p_value <- c(0.01, 1.5, 0.2)
  expect_error(
    write_mummichog_input(tbl, file.path(withr::local_tempdir(), "in.tsv")),
    "p-values"
  )
})

test_that("non-positive / non-finite m/z is rejected", {
  tbl <- make_stats()
  tbl$mz <- c(0, 150, 170)
  expect_error(
    write_mummichog_input(tbl, file.path(withr::local_tempdir(), "in.tsv")),
    "m/z"
  )
})

test_that("out-of-range m/z is warned about and pre-filtered", {
  tbl <- make_stats(n = 3, mz = c(10, 300, 3000))  # 10 < 50 and 3000 > 2000 dropped
  path <- file.path(withr::local_tempdir(), "in.tsv")
  expect_warning(res <- write_mummichog_input(tbl, path, id_col = "feature_id"),
                 "m/z outside")

  out <- readr::read_tsv(res, col_types = readr::cols(.default = readr::col_character()))
  expect_equal(colnames(out), c("mz", "rtime", "p-value", "t-score", "CompoundID"))
  expect_equal(out$CompoundID, "feat_2")            # only the in-range feature survives
  expect_equal(nrow(out), 1L)
})

test_that("all-out-of-range m/z is a hard error", {
  tbl <- make_stats(n = 2, mz = c(10, 5000))
  expect_error(
    suppressWarnings(write_mummichog_input(tbl, file.path(withr::local_tempdir(), "in.tsv"))),
    "No features remain"
  )
})

test_that("feature id is written as the 5th column and mirrored in the id-map", {
  tbl  <- make_stats(n = 3)
  path <- file.path(withr::local_tempdir(), "in.tsv")
  res  <- write_mummichog_input(tbl, path, id_col = "feature_id")

  out <- readr::read_tsv(res, col_types = readr::cols(.default = readr::col_character()))
  expect_equal(out$CompoundID, c("feat_1", "feat_2", "feat_3"))

  idmap <- readr::read_tsv(paste0(res, ".idmap.tsv"),
                           col_types = readr::cols(.default = readr::col_character()))
  expect_equal(colnames(idmap), c("written_row", "feature_id"))
  expect_equal(idmap$feature_id, c("feat_1", "feat_2", "feat_3"))
})

test_that("missing feature id falls back to row numbers in the id-map", {
  tbl <- make_stats(n = 3)
  tbl$feature_id <- NULL
  path <- file.path(withr::local_tempdir(), "in.tsv")
  res  <- write_mummichog_input(tbl, path)            # id_col defaults to NULL

  idmap <- readr::read_tsv(paste0(res, ".idmap.tsv"),
                           col_types = readr::cols(.default = readr::col_character()))
  expect_equal(idmap$feature_id, c("1", "2", "3"))

  out <- readr::read_tsv(res, col_types = readr::cols(.default = readr::col_character()))
  expect_equal(out$CompoundID, c("1", "2", "3"))
})
