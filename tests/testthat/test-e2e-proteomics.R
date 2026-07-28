# tests/testthat/test-e2e-proteomics.R
#
# End-to-end regression test: run the proteomics pipeline for real via
# tar_make() on the repo's shipped synthetic data and assert the DE result.
# Heavy + opt-in: set RUN_E2E_OMICS=1 (and RUN_E2E_OMICS_FULL=1 for the
# file-output checks). Scaffolding lives in helper-e2e.R.

PROT_STAGE  <- e2e_stage_shipped("example_proteomics/proteomics", "proteomics")
PROT_CONFIG <- "e2e_proteomics_config.yaml"

test_that("proteomics runs end-to-end and produces a well-formed DE result", {
  if (!e2e_should_run()) skip("Set RUN_E2E_OMICS=1 to run the heavy end-to-end omic tests.")
  skip_if_not_installed("limma")

  res <- run_omic_e2e(PROT_CONFIG, PROT_STAGE, "prot_de_res")
  de  <- res$values[["prot_de_res"]]

  # --- structure: the DE contract holds (proteomics also needs method + imputations) ---
  expect_true(is.list(de))
  expect_error(assert_de_contract(de, stage = "proteomics"), NA)
  expect_s3_class(de$summary_df, "data.frame")
  expect_gt(nrow(de$summary_df), 0)
  expect_true("pass_any_contrast" %in% names(de$summary_df))
  expect_length(de$de_tables, 1L)   # exactly the one configured contrast

  # The 3 cRAP- contaminants must be filtered out (77 real proteins remain).
  ids <- e2e_feature_ids(de$summary_df)
  expect_false(any(grepl("^cRAP-", ids)))

  # --- numbers: significant count in a sane band ---
  # The generator injects strong DE into the first 20 real proteins.
  # Wide band; tighten after the first green run.
  sig <- e2e_significant_ids(de$summary_df)
  expect_gte(length(sig), 8L)
  expect_lte(length(sig), 70L)
})

test_that("proteomics DE is reproducible across runs (same seed)", {
  if (!e2e_should_run()) skip("Set RUN_E2E_OMICS=1 to run the heavy end-to-end omic tests.")
  skip_if_not_installed("limma")

  r1 <- run_omic_e2e(PROT_CONFIG, PROT_STAGE, "prot_de_res")
  r2 <- run_omic_e2e(PROT_CONFIG, PROT_STAGE, "prot_de_res")

  expect_identical(
    sort(e2e_significant_ids(r1$values[["prot_de_res"]]$summary_df)),
    sort(e2e_significant_ids(r2$values[["prot_de_res"]]$summary_df))
  )
})

test_that("proteomics writes export files (full run)", {
  if (!e2e_should_run_full()) skip("Set RUN_E2E_OMICS_FULL=1 to run the heavier file-output checks.")
  skip_if_not_installed("limma")
  skip_if_not_installed("openxlsx")

  res   <- run_omic_e2e(PROT_CONFIG, PROT_STAGE, "prot_exports")
  files <- list.files(res$out_dir, recursive = TRUE, full.names = TRUE)

  xlsx <- files[grepl("\\.xlsx$", files)]
  tsv  <- files[grepl("\\.tsv$", files)]
  expect_gte(length(xlsx), 1L)
  expect_true(all(file.size(xlsx) > 0))
  expect_gte(length(tsv), 1L)
})
