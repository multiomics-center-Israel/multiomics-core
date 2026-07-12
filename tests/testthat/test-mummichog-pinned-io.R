# tests/testthat/test-mummichog-pinned-io.R
#
# Unit tests for output discovery, file classification, result parsing and the
# feature-id join (06c). A minimal but faithful v2.7.0 result tree is built in a
# temp dir (no Python needed).

# Build a fake mummichog v2 result tree; returns its root directory.
build_fake_tree <- function() {
  root   <- withr::local_tempdir(.local_envir = parent.frame())
  resdir <- file.path(root, "1699999999.99.smoke_test")
  tables <- file.path(resdir, "tables")
  figs   <- file.path(resdir, "figures")
  jsd    <- file.path(resdir, "js")
  for (d in c(tables, figs, jsd)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  writeLines("<html></html>", file.path(resdir, "result.html"))

  writeLines(c(
    "pathway\toverlap_size\tpathway_size\tp-value\toverlap_EmpiricalCompounds (id)\toverlap_features (id)\toverlap_features (name)",
    "PathA\t1\t10\t0.01\tE1\tC1\tname1",
    "PathB\t2\t20\t0.02\tE1,E2\tC1,C2\tnames"
  ), file.path(tables, "mcg_pathwayanalysis_smoke_test.tsv"))
  writeLines("binary-placeholder", file.path(tables, "mcg_pathwayanalysis_smoke_test.xlsx"))

  writeLines(c(
    "MODULE\tp-value\tsize\tmembers (name)\tmembers (id)\tThis module overlaps with",
    "module_1\t0.03\t3\tn\tE1\t-"
  ), file.path(tables, "mcg_modularanalysis_smoke_test.tsv"))

  writeLines(c(
    "EID\tmassfeature_rows\tstr_row_ion\tcompounds\tcompound_names",
    "E1\trow1\trow1_M+H[1+]\tC1\tn1"
  ), file.path(tables, "ListOfEmpiricalCompounds.tsv"))

  writeLines(c(
    "massfeature_rows\tm/z\tretention_time\tp_value\tstatistic\tCompoundID_from_user",
    "row1\t100\t10\t0.01\t2\tfeat_1",
    "row2\t200\t20\t0.02\t-1\tfeat_2",
    "row3\t201\t20\t0.03\t1\tfeat_3"
  ), file.path(tables, "userInputData.txt"))

  # Note the real file repeats the "input_row" column; the reader must tolerate it.
  writeLines(c(
    "input_row\tEID\tstr_row_ion\tcompounds\tcompound_names\tinput_row\tm/z\tretention_time\tp_value\tstatistic\tCompoundID_from_user",
    "row1\tE1\trow1_M+H[1+]\tC1\tn1\trow1\t100\t10\t0.01\t2\tfeat_1",
    "row2\tE2\trow2_M+H[1+]\tC2\tn2\trow2\t200\t20\t0.02\t-1\tfeat_2",
    "row3\tE2\trow3_M+H[1+]\tC3\tn3\trow3\t201\t20\t0.03\t1\tfeat_3"
  ), file.path(tables, "userInput_to_EmpiricalCompounds.tsv"))

  writeLines("{}", file.path(tables, "exported_Compounds.json"))
  writeLines("%PDF-1.4", file.path(figs, "plot_pathwayModel_smoke_test.pdf"))
  writeLines("// js", file.path(jsd, "result.js"))

  root
}

test_that("files are discovered and classified by kind", {
  root  <- build_fake_tree()
  files <- list_mummichog_files(root)
  kinds <- vapply(files, classify_mummichog_file, character(1))

  for (k in c("v2_result_html", "pathway_tsv", "pathway_xlsx", "module_tsv",
              "empirical_compounds", "user_input_to_ec", "user_input_copy",
              "exported_compounds_json", "figure_pdf", "js")) {
    expect_true(k %in% kinds, info = paste("missing kind:", k))
  }
  expect_false(any(kinds == "other"))
})

test_that("pathway and module tables parse", {
  root  <- build_fake_tree()
  files <- list_mummichog_files(root)

  pw <- read_mummichog_pathways(files)
  expect_equal(nrow(pw), 2L)
  expect_true("overlap_EmpiricalCompounds (id)" %in% colnames(pw))

  md <- read_mummichog_modules(files)
  expect_equal(nrow(md), 1L)
})

test_that("empirical map reader tolerates the duplicated input_row column", {
  root  <- build_fake_tree()
  files <- list_mummichog_files(root)
  emap  <- read_mummichog_empirical_map(files)
  expect_setequal(emap$EID, c("E1", "E2"))
  expect_setequal(emap$feature_id, c("feat_1", "feat_2", "feat_3"))
})

test_that("join maps pathways back to feature ids via carried ids (not rowN)", {
  root  <- build_fake_tree()
  files <- list_mummichog_files(root)
  pw    <- read_mummichog_pathways(files)

  joined <- join_features_to_results(pw, files)
  expect_true("feature_ids" %in% colnames(joined))
  expect_equal(joined$feature_ids[joined$pathway == "PathA"], "feat_1")
  expect_equal(joined$feature_ids[joined$pathway == "PathB"], "feat_1;feat_2;feat_3")
})

test_that("manifest records every file with a kind and byte size", {
  root  <- build_fake_tree()
  files <- list_mummichog_files(root)
  mpath <- file.path(root, "mummichog_manifest.tsv")
  mf    <- write_mummichog_manifest(files, mpath)

  man <- readr::read_tsv(mf, show_col_types = FALSE)
  expect_setequal(colnames(man), c("path", "file", "directory", "kind", "bytes"))
  expect_true("pathway_tsv" %in% man$kind)
  expect_true(all(man$bytes >= 0))
})

test_that("readers fail loudly when an expected table is absent", {
  empty <- withr::local_tempdir()
  expect_error(read_mummichog_pathways(list_mummichog_files(empty)), "No pathway")
  expect_error(read_mummichog_modules(list_mummichog_files(empty)), "No modular")
})

test_that("an input placed inside out_dir is staged out so it survives the wipe", {
  root    <- withr::local_tempdir()
  out_dir <- file.path(root, "v2")
  dir.create(out_dir, recursive = TRUE)
  infile  <- file.path(out_dir, "input.tsv")           # nested INSIDE out_dir
  lines   <- c("mz\trtime\tp-value\tt-score\tCompoundID",
               "180.06\t50\t0.01\t2\tfeat_1")
  writeLines(lines, infile)

  staged <- .mmc_stage_infile(infile, out_dir)
  # the staged copy must live OUTSIDE out_dir (compare with a consistent separator)
  expect_false(startsWith(
    normalizePath(staged, winslash = "/"),
    paste0(normalizePath(out_dir, winslash = "/"), "/")
  ))

  unlink(out_dir, recursive = TRUE, force = TRUE)       # simulate run_mummichog's wipe
  expect_true(file.exists(staged))                      # copy survives the wipe
  expect_identical(readLines(staged), lines)            # content preserved
})

test_that("an input outside out_dir is left untouched", {
  root    <- withr::local_tempdir()
  out_dir <- file.path(root, "v2")
  infile  <- file.path(root, "input.tsv")               # OUTSIDE out_dir
  writeLines("x", infile)
  # .mmc_stage_infile returns a winslash = "/" path; compare on the same convention
  expect_identical(.mmc_stage_infile(infile, out_dir),
                   normalizePath(infile, winslash = "/", mustWork = TRUE))
})

test_that(".mmc_default_python resolves to a platform-appropriate venv path", {
  dp <- .mmc_default_python()
  if (.Platform$OS.type == "windows") {
    expect_match(dp, "Scripts/python\\.exe$")
  } else {
    expect_match(dp, "bin/python$")
  }
})

test_that("mod_mummichog_pinned is a no-op when mummichog is disabled or omitted", {
  out_dir <- withr::local_tempdir()

  cfg_off <- list(modes = list(metabolomics = list(
    enrichment = list(mummichog = list(enabled = FALSE)))))
  expect_null(mod_mummichog_pinned(pre = NULL, de_res = NULL,
                                   config = cfg_off, out_dir = out_dir))

  cfg_omit <- list(modes = list(metabolomics = list(enrichment = list())))
  expect_null(mod_mummichog_pinned(pre = NULL, de_res = NULL,
                                   config = cfg_omit, out_dir = out_dir))

  # nothing written, no venv invoked -> the output subdir is never created
  expect_false(dir.exists(file.path(out_dir, "mummichog_pinned")))
})

test_that("pipe_metabolomics adds the pinned mummichog targets only when enabled", {
  skip_if_not_installed("targets")
  library(targets)
  on  <- pipe_metabolomics(chosen_norm = "pqn", mummichog_enabled = TRUE)
  off <- pipe_metabolomics(chosen_norm = "pqn", mummichog_enabled = FALSE)
  # enabling adds exactly the three metab_mummichog_pinned_* targets
  expect_equal(length(on) - length(off), 3L)

  # the stage lives in analysis_core, so it is also present in multiomics runs
  # (skip_outputs = TRUE), matching the pre-switch behaviour
  mo_on  <- pipe_metabolomics(chosen_norm = "pqn", skip_outputs = TRUE, mummichog_enabled = TRUE)
  mo_off <- pipe_metabolomics(chosen_norm = "pqn", skip_outputs = TRUE, mummichog_enabled = FALSE)
  expect_equal(length(mo_on) - length(mo_off), 3L)
})

test_that("mod_mummichog_pinned refuses a non-human organism without a model", {
  cfg <- list(modes = list(metabolomics = list(
    organism   = "Mus musculus",
    enrichment = list(mummichog = list(enabled = TRUE)))))
  expect_error(
    mod_mummichog_pinned(pre = NULL, de_res = NULL, config = cfg,
                         out_dir = withr::local_tempdir()),
    "non-human"
  )
})

test_that("a configured model bypasses the organism guard (human-only) check", {
  # model_json set -> network != human_mfn -> guard skipped; the call then fails
  # later on the missing DE table, NOT with the organism error.
  cfg <- list(modes = list(metabolomics = list(
    organism   = "Mus musculus",
    enrichment = list(mummichog = list(enabled = TRUE, model_json = "/tmp/model.json")))))
  expect_error(
    mod_mummichog_pinned(pre = NULL, de_res = NULL, config = cfg,
                         out_dir = withr::local_tempdir()),
    "DE table"
  )
})

test_that("a human organism passes the guard (fails later, not on organism)", {
  cfg <- list(modes = list(metabolomics = list(
    organism   = "Homo sapiens",
    enrichment = list(mummichog = list(enabled = TRUE)))))
  err <- tryCatch(
    mod_mummichog_pinned(pre = NULL, de_res = NULL, config = cfg,
                         out_dir = withr::local_tempdir()),
    error = function(e) conditionMessage(e)
  )
  expect_false(grepl("non-human", err))   # cleared the organism guard
  expect_match(err, "DE table")           # stopped later, as expected
})
