# tests/testthat/test-e2e-rnaseq-bundle.R
#
# End-to-end check that a single-omics RNA-seq run writes its omics bundle
# through {targets}: rna_bundle is built by tar_make() in an isolated temp
# project and store (run_omic_e2e(), helper-e2e.R), then validated and read
# back. Uses the synthetic RNA fixture as it is: no enrichment block, so the
# pathway step returns its empty result without network access.
#
# Heavy + opt-in: set RUN_E2E_OMICS=1. Kept in its own file so it can be run
# alone, e.g. on Windows.

test_that("a single-omics RNA run writes a valid bundle through tar_make()", {
  if (!e2e_should_run()) skip("Set RUN_E2E_OMICS=1 to run the heavy end-to-end omic tests.")
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("edgeR")

  res  <- run_omic_e2e("e2e_rnaseq_config.yaml", e2e_stage_rna(seed = 1L), "rna_bundle")
  path <- res$values[["rna_bundle"]]

  # Stored under <paths.out>/bundles/rnaseq/, beside the Results_ folder and
  # not inside it, so a rerun of the round cannot overwrite it. Paths are
  # normalised so the comparison holds with Windows separators too.
  norm <- function(p) normalizePath(p, winslash = "/", mustWork = TRUE)
  expect_true(dir.exists(path))
  expect_identical(norm(dirname(path)),
                   norm(file.path(res$proj, res$config$paths$out, "bundles", "rnaseq")))
  expect_false(grepl("/Results_", norm(path), fixed = TRUE))

  # Every checksum is verified on read.
  b <- read_omics_bundle(path)
  expect_identical(b$manifest$omics, "rnaseq")
  expect_identical(b$manifest$project, res$config$project$name)
  expect_identical(b$manifest$analysis_round, res$config$project$analysis_round)

  # The run's objects are there; the DESeq2 model and load bookkeeping are not.
  expect_true(all(c("pre", "de", "pathway", "config") %in% names(b$parts)))
  expect_true(all(c("expr_raw", "expr_filt", "expr_work", "meta", "row_data")
                  %in% names(b$parts$pre)))
  expect_length(b$parts$de$tables, 1L)
  expect_false("dds" %in% names(b$parts$de))
  expect_false(".config_path" %in% names(b$parts$config))
  expect_identical(nrow(list_omics_bundles(file.path(res$proj, res$config$paths$out),
                                           omics = "rnaseq")), 1L)
})
