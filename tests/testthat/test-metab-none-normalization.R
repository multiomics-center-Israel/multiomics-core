# tests/testthat/test-metab-none-normalization.R
#
# Tests for chosen_norm = "none" — the "already normalized upstream" path.
# "none" must skip sample normalization and select the transform-only matrix
# (met_log), while transform/scaling stay governed by their own config keys.
#
# All fixtures are synthetic (no filesystem, no real data). Drift correction is
# disabled in the config so the selected matrix is compared verbatim.

# --- fixture: a norm-stage list as produced by the mod_met_normalize_* wrappers
make_stage <- function(mat, meta) {
  list(mat = mat, meta = meta, row_data = data.frame(feature_id = rownames(mat)))
}

# --- fixture: minimal metabolomics config with drift correction turned off
make_cfg <- function(chosen_norm = "none", scaling = "none") {
  list(modes = list(metabolomics = list(
    preprocessing = list(
      chosen_norm       = chosen_norm,
      scaling           = scaling,
      transform         = "log2",
      drift_correction  = list(enabled = FALSE)
    ),
    effects = list(samples = "sample_id")
  )))
}

# Distinct matrices per stage so a wrong switch branch would fail loudly.
make_fixtures <- function() {
  cols <- c("S1", "S2", "S3", "S4")
  meta <- data.frame(sample_id = cols, stringsAsFactors = FALSE)
  logged <- make_stage(matrix(c(1, 2, 3, 4,
                                5, 6, 7, 8),
                              nrow = 2, byrow = TRUE,
                              dimnames = list(c("f1", "f2"), cols)), meta)
  # Deliberately different values in the other stages.
  tss    <- make_stage(logged$mat + 100, meta)
  median <- make_stage(logged$mat + 200, meta)
  pqn    <- make_stage(logged$mat + 300, meta)
  list(logged = logged, tss = tss, median = median, pqn = pqn, meta = meta)
}

test_that("mod_met_corrected with chosen_norm='none' selects the log-only matrix", {
  fx  <- make_fixtures()
  cfg <- make_cfg(chosen_norm = "none", scaling = "none")

  out <- suppressMessages(mod_met_corrected(
    norm_tss    = fx$tss,
    norm_median = fx$median,
    norm_pqn    = fx$pqn,
    logged      = fx$logged,
    meta        = fx$meta,
    out_dir     = tempdir(),
    config      = cfg
  ))

  # No sample normalization applied -> identical to met_log's matrix.
  expect_equal(out$mat, fx$logged$mat)
  expect_equal(out$info$chosen_norm, "none")
  expect_false(out$info$drift_applied)
})

test_that("chosen_norm='none' still honors the scaling config key", {
  fx  <- make_fixtures()
  cfg <- make_cfg(chosen_norm = "none", scaling = "center")

  out <- suppressMessages(mod_met_corrected(
    norm_tss    = fx$tss,
    norm_median = fx$median,
    norm_pqn    = fx$pqn,
    logged      = fx$logged,
    meta        = fx$meta,
    out_dir     = tempdir(),
    config      = cfg
  ))

  # "none" skips sample normalization but scaling still applies on top.
  expect_equal(out$mat, scale_metab(fx$logged$mat, method = "center"))
})

test_that("pipe_metabolomics accepts 'none' as a valid chosen_norm", {
  skip_if_not_installed("targets")
  library(targets)
  # "none" is analysis mode (not QC-review), so the plan must include
  # met_corrected. NULL would instead return the QC-review targets.
  plan <- pipe_metabolomics(chosen_norm = "none", skip_outputs = TRUE)
  target_names <- vapply(plan, function(t) t$settings$name, character(1))
  expect_true("met_corrected" %in% target_names)
})
