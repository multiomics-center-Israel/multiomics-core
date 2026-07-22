# tests/testthat/test-metab-bio-factor-norm.R
#
# Tests for biological-factor ("bio_factor") metabolomics normalization:
#   - core divide-by-covariate math                     (norm_biological_factor)
#   - sample alignment via the configured sample column (effects$samples)
#   - the pipeline-safety guard: unconfigured -> module returns NULL, never
#     stop()s (so non-bio_factor runs don't break at tar_make())
#   - input validation (missing column, unknown samples, NA/<=0 covariate)
#
# All fixtures are synthetic (no filesystem, no real data).

# --- fixture: build a small metadata table with named sample + covariate cols
make_meta <- function(samples, factor_vals,
                      sample_col = "SampleID", factor_col = "total_protein") {
  m <- data.frame(a = samples, b = factor_vals, stringsAsFactors = FALSE)
  names(m) <- c(sample_col, factor_col)
  m
}

test_that("norm_biological_factor divides each sample column by its covariate", {
  mat <- matrix(c(10, 20, 30, 40,
                  100, 200, 300, 400,
                  1,   2,   3,   4),
                nrow = 3, byrow = TRUE,
                dimnames = list(c("f1", "f2", "f3"), c("S1", "S2", "S3", "S4")))
  meta <- make_meta(c("S1", "S2", "S3", "S4"), c(2, 4, 5, 10))

  out <- norm_biological_factor(mat, meta, factor_col = "total_protein",
                                sample_col = "SampleID")

  expect_equal(dim(out), dim(mat))
  expect_equal(out[, "S1"], mat[, "S1"] / 2)
  expect_equal(out[, "S4"], mat[, "S4"] / 10)
})

test_that("norm_biological_factor aligns via sample_col when meta rows are shuffled", {
  mat <- matrix(1:8, nrow = 2,
                dimnames = list(c("f1", "f2"), c("S1", "S2", "S3", "S4")))
  # meta rows deliberately in a different order than mat columns
  meta <- make_meta(c("S3", "S1", "S4", "S2"), c(3, 1, 4, 2))

  out <- norm_biological_factor(mat, meta, factor_col = "total_protein",
                                sample_col = "SampleID")

  expect_equal(out[, "S1"], mat[, "S1"] / 1)
  expect_equal(out[, "S2"], mat[, "S2"] / 2)
  expect_equal(out[, "S3"], mat[, "S3"] / 3)
  expect_equal(out[, "S4"], mat[, "S4"] / 4)
})

test_that("norm_biological_factor falls back to conventional sample columns", {
  mat <- matrix(1:4, nrow = 1, dimnames = list("f1", c("S1", "S2", "S3", "S4")))
  meta <- make_meta(c("S1", "S2", "S3", "S4"), c(1, 2, 4, 4), sample_col = "sample_id")

  # sample_col NULL -> should still find the conventional "sample_id" column
  out <- norm_biological_factor(mat, meta, factor_col = "total_protein")
  expect_equal(out[, "S3"], mat[, "S3"] / 4)
})

test_that("norm_biological_factor validates column, samples, and covariate values", {
  mat  <- matrix(1:4, nrow = 1, dimnames = list("f1", c("S1", "S2", "S3", "S4")))
  meta <- make_meta(c("S1", "S2", "S3", "S4"), c(1, 2, 3, 4))

  expect_error(norm_biological_factor(mat, meta, factor_col = "nope"), "not found")

  meta_missing <- make_meta(c("S1", "S2", "S3"), c(1, 2, 3))
  expect_error(
    norm_biological_factor(mat, meta_missing, "total_protein", sample_col = "SampleID"),
    "missing from meta"
  )

  meta_bad <- make_meta(c("S1", "S2", "S3", "S4"), c(1, 0, 3, 4))  # 0 is invalid
  expect_error(
    norm_biological_factor(mat, meta_bad, "total_protein", sample_col = "SampleID"),
    "invalid"
  )
})

test_that("normalize_samples dispatches bio_factor and threads sample_col", {
  mat  <- matrix(c(4, 8, 2, 6), nrow = 1, dimnames = list("f1", c("A", "B", "C", "D")))
  meta <- make_meta(c("A", "B", "C", "D"), c(2, 2, 1, 3), sample_col = "run_id")

  out <- normalize_samples(mat, method = "bio_factor", meta = meta,
                           bio_factor_col = "total_protein", sample_col = "run_id")

  expect_equal(out[, "A"], mat[, "A"] / 2)
  expect_equal(out[, "C"], mat[, "C"] / 1)
})

test_that("mod_met_normalize_bio_factor returns NULL (not stop) when unconfigured", {
  data <- list(
    mat      = matrix(1:4, nrow = 1, dimnames = list("f1", c("S1", "S2", "S3", "S4"))),
    meta     = make_meta(c("S1", "S2", "S3", "S4"), c(1, 2, 3, 4)),
    row_data = data.frame(feature = "f1")
  )
  cfg <- list(modes = list(metabolomics = list(preprocessing = list())))  # no biological_factor_col

  expect_null(suppressMessages(mod_met_normalize_bio_factor(data, cfg)))
})

test_that("mod_met_normalize_bio_factor normalizes then log2s, honoring effects$samples", {
  mat  <- matrix(c(2, 4, 8, 16), nrow = 1, dimnames = list("f1", c("S1", "S2", "S3", "S4")))
  meta <- make_meta(c("S1", "S2", "S3", "S4"), c(1, 2, 4, 8), sample_col = "MySample")
  data <- list(mat = mat, meta = meta, row_data = data.frame(feature = "f1"))
  cfg  <- list(modes = list(metabolomics = list(
    preprocessing = list(chosen_norm = "bio_factor",
                         biological_factor_col = "total_protein",
                         pseudocount = 1),
    effects = list(samples = "MySample")
  )))

  out <- suppressMessages(mod_met_normalize_bio_factor(data, cfg))

  expect_false(is.null(out))
  # 2/1, 4/2, 8/4, 16/8 all == 2; then log2(2 + 1)
  expect_equal(unname(out$mat["f1", ]), rep(log2(2 + 1), 4))
})
