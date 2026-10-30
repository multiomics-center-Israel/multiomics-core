# tests/testthat/test-metab-drift-opt-in.R
#
# LOESS drift correction is opt-in: apply_drift_correction() runs it only for an
# explicit `preprocessing.drift_correction.enabled: true`. A missing block, an
# empty block, a block without `enabled`, and `enabled: false` all leave the
# matrix untouched. Before this contract a missing `enabled` meant "on".
#
# Every fixture here MEETS the LOESS preconditions (injection order, QC flag,
# >= 4 QC samples), so a test passes only because the setting is respected, not
# because the metadata happened to make drift correction skip itself.
#
# This file checks activation and validation only. Whether the LOESS correction
# itself is mathematically right is a separate question (P3 issue).
#
# All data are synthetic.

# --- fixture: a log-scale matrix with a linear injection-order trend --------
make_drift_fixture <- function() {
  withr::with_seed(11, {
    n_feat <- 6
    n_samp <- 16
    ids    <- sprintf("S%02d", seq_len(n_samp))
    trend  <- 0.08 * seq_len(n_samp)
    mat <- matrix(stats::rnorm(n_feat * n_samp, mean = 20, sd = 0.3),
                  nrow = n_feat,
                  dimnames = list(sprintf("f%d", seq_len(n_feat)), ids))
    mat <- sweep(mat, 2, trend, FUN = "+")
  })
  # 8 QC samples, including the first and last injection. loess() cannot
  # extrapolate, so a sample injected outside the QC range would make every
  # feature skip. With 8 QCs the span stays at 0.75 (~6 points per local
  # quadratic fit); with the 4-QC minimum it is forced to 1 and the fit has no
  # residual degrees of freedom left (#362).
  # The fixture guarantees every feature is corrected: no missing QC values,
  # QCs bracket the run, and values (~20) are far above epsilon.
  meta <- data.frame(
    sample_id       = ids,
    injection_order = seq_len(n_samp),
    is_QC           = seq_len(n_samp) %in% c(1, 3, 5, 7, 10, 12, 14, 16),
    stringsAsFactors = FALSE
  )
  list(mat = mat, meta = meta)
}

# Asserts that LOESS really corrected every feature of the fixture: counts,
# finiteness, and a change in every row. Identical-but-skipped results (or two
# identical invalid ones) cannot pass this.
expect_every_feature_corrected <- function(corrected, original, info) {
  expect_identical(info$n_corrected, nrow(original))
  expect_identical(info$n_skipped, 0L)
  expect_true(all(is.finite(corrected)))
  expect_true(all(rowSums(corrected != original) > 0))
}

# Mode config (config$modes$metabolomics) parsed from YAML text, the same way a
# config file reaches the pipeline.
cfg_from_yaml <- function(txt) {
  cfg <- yaml::yaml.load(txt)
  cfg$effects <- cfg$effects %||% list(samples = "sample_id")
  cfg
}

# The preprocessing block the run.R wizard writes (wizard_metabolomics() emits
# only chosen_norm and the missingness thresholds, never drift_correction).
# The wizard's separate `normalization:` block is left out: validation rejects
# it (#280), so a config carrying it never reaches drift correction at all.
wizard_preprocessing_yaml <- paste(
  "preprocessing:",
  "  chosen_norm: \"pqn\"",
  "  feat_missing_threshold: 0.80",
  "  sample_missing_threshold: 0.99",
  sep = "\n"
)

off_configs <- list(
  "no preprocessing block"         = "effects:\n  samples: \"sample_id\"",
  "no drift_correction block"      = "preprocessing:\n  scaling: \"none\"",
  "drift_correction with no value" = "preprocessing:\n  drift_correction:",
  "empty drift_correction block"   = "preprocessing:\n  drift_correction: {}",
  "block without enabled"          = "preprocessing:\n  drift_correction:\n    loess_span: 0.75",
  "enabled: false"                 = "preprocessing:\n  drift_correction:\n    enabled: false",
  "wizard-generated preprocessing" = wizard_preprocessing_yaml
)

test_that("drift correction stays off unless enabled is explicitly true", {
  fx <- make_drift_fixture()
  for (label in names(off_configs)) {
    cfg <- cfg_from_yaml(off_configs[[label]])
    out <- suppressMessages(apply_drift_correction(fx$mat, fx$meta, cfg))
    expect_false(out$applied, info = label)
    expect_identical(out$mat, fx$mat, info = label)
  }
})

test_that("enabled: true runs the existing LOESS correction unchanged", {
  fx  <- make_drift_fixture()
  cfg <- cfg_from_yaml("preprocessing:\n  drift_correction:\n    enabled: true")

  out <- suppressMessages(apply_drift_correction(fx$mat, fx$meta, cfg))
  expect_true(out$applied)

  # Same call the dispatcher makes with its defaults: activation is the only
  # thing this contract changes, not the correction.
  direct <- correct_drift_loess(fx$mat, fx$meta,
                                sample_col          = "sample_id",
                                injection_order_col = "injection_order",
                                qc_flag_col         = "is_QC",
                                loess_span          = 0.75,
                                epsilon             = 1e-8)
  expect_identical(out$mat, direct$mat_corrected)
  expect_every_feature_corrected(out$mat, fx$mat, out$info)
})

test_that("mod_met_corrected does not drift-correct a wizard-style config", {
  fx  <- make_drift_fixture()
  stage <- list(mat = fx$mat, meta = fx$meta,
                row_data = data.frame(feature_id = rownames(fx$mat)))
  cfg <- list(modes = list(metabolomics = cfg_from_yaml(wizard_preprocessing_yaml)))

  out <- suppressMessages(mod_met_corrected(
    norm_tss = stage, norm_median = stage, norm_pqn = stage, logged = stage,
    meta = fx$meta, out_dir = tempdir(), config = cfg
  ))
  expect_false(out$info$drift_applied)
  expect_identical(out$mat, fx$mat)
})

# --- validation --------------------------------------------------------------

metab_cfg_with_drift <- function(dc) {
  cfg <- list(
    input = list(format = "processed_wide"),
    files = list(data = "data.csv"),
    preprocessing = list(chosen_norm = "pqn")
  )
  cfg$preprocessing["drift_correction"] <- list(dc)
  cfg
}

test_that("validation accepts true, false, missing and empty drift settings", {
  for (dc in list(NULL, list(), list(enabled = TRUE), list(enabled = FALSE),
                  list(loess_span = 0.5))) {
    expect_silent(validate_metabolomics_config(metab_cfg_with_drift(dc)))
  }
})

test_that("validation rejects a non-logical enabled and a non-block drift_correction", {
  expect_error(validate_metabolomics_config(metab_cfg_with_drift(list(enabled = "true"))),
               "drift_correction\\$enabled")
  expect_error(validate_metabolomics_config(metab_cfg_with_drift(list(enabled = c(TRUE, FALSE)))),
               "drift_correction\\$enabled")
  expect_error(validate_metabolomics_config(metab_cfg_with_drift(TRUE)),
               "must be a block")
})

test_that("the shipped metabolomics template keeps drift correction off and validates", {
  root <- normalizePath(if (dir.exists("R")) "." else "../..")
  tmpl <- yaml::read_yaml(file.path(root, "config/templates/metabolomics_template.yaml"))
  metab <- tmpl$modes$metabolomics
  expect_false(isTRUE(metab$preprocessing$drift_correction$enabled))
  expect_silent(validate_metabolomics_config(metab))
})

# --- safeguard: LOESS on feature-scaled data (temporary, see #362) -----------

run_corrected <- function(fx, scaling, drift_enabled, meta = fx$meta) {
  stage <- list(mat = fx$mat, meta = meta,
                row_data = data.frame(feature_id = rownames(fx$mat)))
  cfg <- list(modes = list(metabolomics = list(
    preprocessing = list(chosen_norm = "pqn", scaling = scaling,
                         drift_correction = list(enabled = drift_enabled)),
    effects = list(samples = "sample_id")
  )))
  suppressMessages(mod_met_corrected(
    norm_tss = stage, norm_median = stage, norm_pqn = stage, logged = stage,
    meta = meta, out_dir = tempdir(), config = cfg
  ))
}

# On centered data the existing LOESS skips most features and warns once per
# skipped feature. That is exactly the #362 defect the safeguard stops; these
# tests assert the stop, so only that one warning is muffled, nothing else.
muffle_loess_skip_warnings <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (startsWith(conditionMessage(w),
                   "correct_drift_loess: predicted values < epsilon")) {
      invokeRestart("muffleWarning")
    }
  })
}

test_that("explicit LOESS with any feature scaling stops with an actionable error", {
  fx <- make_drift_fixture()
  for (sc in c("center", "auto", "pareto", "range")) {
    expect_error(
      muffle_loess_skip_warnings(run_corrected(fx, scaling = sc, drift_enabled = TRUE)),
      "cannot run on feature-scaled data", info = sc)
  }
})

test_that("explicit LOESS with scaling 'none' still runs", {
  fx  <- make_drift_fixture()
  out <- run_corrected(fx, scaling = "none", drift_enabled = TRUE)
  expect_true(out$info$drift_applied)
  expect_every_feature_corrected(out$mat, fx$mat, out$info$drift_info)
})

test_that("the safeguard does not fire when LOESS could not run anyway", {
  fx   <- make_drift_fixture()
  # No injection-order column: the preconditions fail, nothing is corrected,
  # so there is nothing invalid to refuse.
  meta <- fx$meta[, c("sample_id", "is_QC")]
  out  <- run_corrected(fx, scaling = "auto", drift_enabled = TRUE, meta = meta)
  expect_false(out$info$drift_applied)
  expect_equal(out$mat, scale_metab(fx$mat, method = "auto"))
})
