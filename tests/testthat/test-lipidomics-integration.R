# tests/testthat/test-lipidomics-integration.R
#
# Integration test for the lipidomics pipeline on the committed synthetic dataset
# (fixtures/test_lipidomics_config.yaml -> test_data/synthetic_lipidomics.csv).
# Exercises load -> preprocess -> QC -> DE, plus the downstream stages (feature
# selection, lipid-class analysis, HTML report) when their optional deps exist.
#
# Rewritten from a narrative cat()/tryCatch() script that asserted nothing and
# wrote outputs into the repo tree (test_outputs/): it now uses real expect_*()
# so breakage fails the suite, and keeps every artifact inside a per-test temp
# dir. R/ is already sourced by helper.R, so no manual sourcing here.

#' TRUE when the fixture and the lipidomics entry-point functions are available.
lipid_fixture_available <- function() {
  cfg_path <- testthat::test_path("fixtures", "test_lipidomics_config.yaml")
  file.exists(cfg_path) &&
    exists("load_lipidomics_inputs") &&
    exists("preprocess_lipidomics")
}

#' Load the lipidomics fixture config, anchored to the repo root so the
#' fixture's relative data paths resolve regardless of testthat's working dir.
lipid_config <- function() {
  config <- load_config(testthat::test_path("fixtures", "test_lipidomics_config.yaml"))
  config$project$dir <- normalizePath(testthat::test_path("..", ".."))
  config
}

test_that("lipidomics: config validates and inputs load", {
  skip_if_not(lipid_fixture_available(), "lipidomics fixture/functions unavailable")

  config <- lipid_config()
  expect_error(validate_lipidomics_config(config$modes$lipidomics), NA)

  inputs <- load_lipidomics_inputs(config)
  expect_type(inputs, "list")
  expect_gt(nrow(inputs$data), 0)
  expect_gt(ncol(inputs$data), 0)
})

test_that("lipidomics: preprocessing produces the `pre` contract", {
  skip_if_not(lipid_fixture_available(), "lipidomics fixture/functions unavailable")

  config <- lipid_config()
  inputs <- load_lipidomics_inputs(config)
  pre    <- preprocess_lipidomics(inputs, config)

  expect_true(all(c("expr_raw", "expr_filt", "expr_work", "meta", "row_data", "info")
                  %in% names(pre)))
  expect_gt(nrow(pre$expr_raw), 0)
  expect_gt(nrow(pre$expr_filt), 0)
  expect_equal(ncol(pre$expr_work), nrow(pre$meta))
  expect_true("lipid_class" %in% names(pre$row_data))
})

test_that("lipidomics: QC and DE modules run and the DE result is well-formed", {
  skip_if_not(lipid_fixture_available(), "lipidomics fixture/functions unavailable")
  skip_if_not_installed("ggplot2")

  config  <- lipid_config()
  out_dir <- withr::local_tempdir()
  inputs  <- load_lipidomics_inputs(config)
  pre     <- preprocess_lipidomics(inputs, config)

  qc_res <- mod_lipidomics_qc_pre(pre, config, out_dir)
  expect_type(qc_res, "list")
  expect_gte(length(qc_res$plots), 1L)

  de_res <- mod_lipidomics_de(pre, config, out_dir)
  expect_error(assert_de_contract(de_res, stage = "lipidomics"), NA)
  expect_s3_class(de_res$summary_df, "data.frame")
  expect_gt(nrow(de_res$summary_df), 0)
  expect_true("pass_any_contrast" %in% names(de_res$summary_df))
  expect_length(de_res$de_tables, 1L)
})

test_that("lipidomics: downstream stages run when their deps are available", {
  skip_if_not(lipid_fixture_available(), "lipidomics fixture/functions unavailable")
  skip_if_not_installed("ggplot2")

  config  <- lipid_config()
  # Trim the RF forest for test speed; the fixture uses 500 in production.
  config$modes$lipidomics$rf$n_trees <- 100
  out_dir <- withr::local_tempdir()
  inputs  <- load_lipidomics_inputs(config)
  pre     <- preprocess_lipidomics(inputs, config)
  qc_res  <- mod_lipidomics_qc_pre(pre, config, out_dir)
  de_res  <- mod_lipidomics_de(pre, config, out_dir)

  # Feature selection (RF / PLS-DA). The module swallows backend errors and
  # returns NULL, so gate on backend availability and then require real results
  # — otherwise a broken learner would leave this test green.
  have_rf    <- requireNamespace("ranger", quietly = TRUE) ||
                requireNamespace("randomForest", quietly = TRUE)
  have_plsda <- requireNamespace("mixOmics", quietly = TRUE)
  if (!have_rf && !have_plsda) {
    skip("no RF/PLS-DA backend installed (ranger/randomForest/mixOmics)")
  }
  fs_res <- mod_lipidomics_feature_selection(pre, config, out_dir)
  expect_false(is.null(fs_res))                       # at least one backend ran
  if (have_rf) {
    expect_s3_class(fs_res$rf$importance_df, "data.frame")
    expect_gt(nrow(fs_res$rf$importance_df), 0)
  }
  if (have_plsda) {
    expect_s3_class(fs_res$plsda$vip_df, "data.frame")
    expect_gt(nrow(fs_res$plsda$vip_df), 0)
  }

  # Lipid-class analysis: class composition is computed from the lipid_class
  # column parsed at preprocessing, so it must be present for this fixture (the
  # module otherwise swallows per-computation errors and still returns a list).
  class_res <- mod_lipidomics_class_analysis(pre, de_res, config, out_dir)
  expect_type(class_res, "list")
  expect_false(is.null(class_res$class_comp))
  expect_false(is.null(class_res$class_comp$class_norm))

  # HTML report — needs a working pandoc plus the Rmd template's own packages
  # (DT, etc.). Skip cleanly when pandoc or any such package is missing rather
  # than failing the suite; a real report bug still surfaces as an error.
  # Report templates resolve relative to the repo root, so render from there.
  skip_if_not_installed("rmarkdown")
  if (!rmarkdown::pandoc_available()) skip("pandoc not available")
  report_path <- tryCatch(
    withr::with_dir(
      config$project$dir,
      mod_lipidomics_report(pre, qc_res, de_res, fs_res, class_res, config, out_dir)
    ),
    error = function(e) {
      if (grepl("there is no package called", conditionMessage(e), fixed = TRUE)) {
        skip(paste("report dependency missing:", conditionMessage(e)))
      }
      stop(e)
    }
  )
  expect_true(file.exists(report_path))
  expect_gt(file.info(report_path)$size, 0)
})
