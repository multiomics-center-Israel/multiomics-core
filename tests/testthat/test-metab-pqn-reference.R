# tests/testthat/test-metab-pqn-reference.R
#
# The PQN reference is an explicit choice (#367):
#   - preprocessing.pqn_reference: "all_samples" (default) uses every sample,
#     even when the metadata flags QC samples;
#   - "pooled_qc" uses the samples flagged in drift_correction.qc_flag_col and
#     stops, rather than falling back, when fewer than 2 are available;
#   - the reference actually used is recorded in info$normalization.
# norm_pqn() itself is unchanged; these tests only check which reference it is
# given. All data are synthetic.

# 12 features x 8 samples on a linear scale; S7 and S8 are pooled QCs.
make_pqn_fixture <- function(qc_values = c(rep(FALSE, 6), TRUE, TRUE)) {
  ids <- sprintf("S%d", 1:8)
  mat <- withr::with_seed(367, matrix(
    exp(stats::rnorm(12 * 8, mean = 8, sd = 1)), nrow = 12,
    dimnames = list(sprintf("f%d", 1:12), ids)
  ))
  meta <- data.frame(sample_id = ids, stringsAsFactors = FALSE)
  if (!is.null(qc_values)) meta$is_QC <- qc_values
  list(mat = mat, meta = meta, row_data = data.frame(feature_id = rownames(mat)))
}

pqn_cfg <- function(pqn_reference = NULL, chosen_norm = "pqn") {
  pre <- list(chosen_norm = chosen_norm, scaling = "none",
              drift_correction = list(enabled = FALSE))
  pre["pqn_reference"] <- list(pqn_reference)
  list(modes = list(metabolomics = list(
    preprocessing = pre,
    effects = list(samples = "sample_id")
  )))
}

# What mod_met_normalize_linear() should return for a given reference, from
# the unchanged norm_pqn().
expected_pqn <- function(mat, qc_idx = NULL) {
  transform_metab(norm_pqn(mat, qc_idx = qc_idx), method = "log2", pseudocount = 1)
}

run_pqn <- function(fx, cfg) {
  suppressMessages(mod_met_normalize_linear(fx, method = "pqn", config = cfg))
}


test_that("the fixture separates the two references", {
  fx <- make_pqn_fixture()
  expect_false(isTRUE(all.equal(expected_pqn(fx$mat),
                                expected_pqn(fx$mat, qc_idx = 7:8))))
})

test_that("the default reference is all samples, even when QC samples are flagged", {
  fx  <- make_pqn_fixture()
  out <- run_pqn(fx, pqn_cfg())
  expect_identical(out$mat, expected_pqn(fx$mat))
  expect_identical(out$pqn_reference, "all_samples")
  expect_identical(out$pqn_reference_n_samples, 8L)

  explicit <- run_pqn(fx, pqn_cfg("all_samples"))
  expect_identical(explicit$mat, out$mat)
})

test_that("all_samples needs no QC column and does not warn", {
  fx <- make_pqn_fixture(qc_values = NULL)
  expect_no_warning(out <- run_pqn(fx, pqn_cfg()))
  expect_identical(out$mat, expected_pqn(fx$mat))
})

test_that("pooled_qc uses the flagged QC samples as the reference", {
  fx  <- make_pqn_fixture()
  out <- run_pqn(fx, pqn_cfg("pooled_qc"))
  expect_identical(out$mat, expected_pqn(fx$mat, qc_idx = 7:8))
  expect_identical(out$pqn_reference, "pooled_qc")
  expect_identical(out$pqn_reference_n_samples, 2L)
})

test_that("pooled_qc stops instead of falling back when QC samples are missing", {
  # No QC column at all.
  expect_error(run_pqn(make_pqn_fixture(qc_values = NULL), pqn_cfg("pooled_qc")),
               "needs the QC flag column 'is_QC'")
  # Only one QC sample.
  expect_error(run_pqn(make_pqn_fixture(qc_values = c(rep(FALSE, 7), TRUE)),
                       pqn_cfg("pooled_qc")),
               "needs at least 2 pooled QC samples")
  # Values that as.logical() does not read as TRUE.
  expect_error(run_pqn(make_pqn_fixture(qc_values = c(rep("no", 6), "yes", "yes")),
                       pqn_cfg("pooled_qc")),
               "marks 0")
})

test_that("an unsupported pqn_reference value is rejected", {
  expect_error(run_pqn(make_pqn_fixture(), pqn_cfg("pools")),
               "preprocessing\\$pqn_reference' must be one of")
})

test_that("TSS is unaffected by pqn_reference", {
  fx  <- make_pqn_fixture(qc_values = NULL)
  out <- suppressMessages(mod_met_normalize_linear(fx, method = "tss",
                                                   config = pqn_cfg("pooled_qc")))
  expect_null(out$pqn_reference)
})

test_that("the reference used is recorded in info$normalization only for PQN", {
  fx  <- make_pqn_fixture()
  run_corrected <- function(cfg) {
    pqn <- run_pqn(fx, cfg)
    suppressMessages(mod_met_corrected(
      norm_tss = pqn, norm_median = pqn, norm_pqn = pqn, logged = pqn,
      meta = fx$meta, out_dir = tempdir(), config = cfg
    ))
  }

  ni <- run_corrected(pqn_cfg("pooled_qc"))$info$normalization
  expect_identical(ni$pqn_reference, "pooled_qc")
  expect_identical(ni$pqn_reference_n_samples, 2L)

  ni_default <- run_corrected(pqn_cfg())$info$normalization
  expect_identical(ni_default$pqn_reference, "all_samples")

  ni_tss <- run_corrected(pqn_cfg(chosen_norm = "tss"))$info$normalization
  expect_null(ni_tss$pqn_reference)
})
