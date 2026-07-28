# tests/testthat/test-e2e-metabolomics.R
#
# End-to-end regression test: run the metabolomics pipeline for real via
# tar_make() on the repo's shipped synthetic data and assert the DE result.
# Heavy + opt-in: set RUN_E2E_OMICS=1 (and RUN_E2E_OMICS_FULL=1 for the
# file-output checks). Scaffolding lives in helper-e2e.R.

METAB_STAGE  <- e2e_stage_shipped("example_metabolomics/metabolomics", "metabolomics")
METAB_CONFIG <- "e2e_metabolomics_config.yaml"

test_that("metabolomics runs end-to-end and recovers the injected DE signal", {
  if (!e2e_should_run()) skip("Set RUN_E2E_OMICS=1 to run the heavy end-to-end omic tests.")
  skip_if_not_installed("limma")

  res <- run_omic_e2e(METAB_CONFIG, METAB_STAGE, "metab_de_res")
  de  <- res$values[["metab_de_res"]]

  # --- structure: the DE contract holds ---
  expect_true(is.list(de))
  expect_error(assert_de_contract(de, stage = "metabolomics"), NA)
  expect_s3_class(de$summary_df, "data.frame")
  expect_gt(nrow(de$summary_df), 0)
  expect_true("pass_any_contrast" %in% names(de$summary_df))
  expect_length(de$de_tables, 1L)   # exactly the one configured contrast

  # --- numbers: significant count in a sane band + injected signal recovered ---
  # The generator injects strong DE into METAB_001..METAB_025 (25 of 150).
  # Bands are intentionally wide; tighten to a snug range after the first green run.
  sig <- e2e_significant_ids(de$summary_df)
  expect_gte(length(sig), 15L)
  expect_lte(length(sig), 90L)

  truth  <- sprintf("METAB_%03d", 1:25)
  recall <- mean(truth %in% sig)
  expect_gte(recall, 0.6)
})

test_that("metabolomics DE is reproducible across runs (same seed)", {
  if (!e2e_should_run()) skip("Set RUN_E2E_OMICS=1 to run the heavy end-to-end omic tests.")
  skip_if_not_installed("limma")

  r1 <- run_omic_e2e(METAB_CONFIG, METAB_STAGE, "metab_de_res")
  r2 <- run_omic_e2e(METAB_CONFIG, METAB_STAGE, "metab_de_res")

  expect_identical(
    sort(e2e_significant_ids(r1$values[["metab_de_res"]]$summary_df)),
    sort(e2e_significant_ids(r2$values[["metab_de_res"]]$summary_df))
  )
})

test_that("metabolomics writes final-results files (full run)", {
  if (!e2e_should_run_full()) skip("Set RUN_E2E_OMICS_FULL=1 to run the heavier file-output checks.")
  skip_if_not_installed("limma")
  skip_if_not_installed("openxlsx")

  res   <- run_omic_e2e(METAB_CONFIG, METAB_STAGE,
                        c("metab_standard_outputs", "metab_final_results"))
  files <- list.files(res$out_dir, recursive = TRUE, full.names = TRUE)

  xlsx <- files[grepl("\\.xlsx$", files)]
  tsv  <- files[grepl("\\.tsv$", files)]
  expect_gte(length(xlsx), 1L)
  expect_true(all(file.size(xlsx) > 0))
  expect_gte(length(tsv), 1L)
})
