# tests/testthat/test-mummichog-pinned-smoke.R
#
# End-to-end smoke test for the pinned mummichog v2 stage. Skips gracefully when
# no venv is available. By default it uses $MUMMICHOG_PYTHON; set
# MUMMICHOG_SMOKE_INSTALL=1 to let the test build a throwaway venv and
# pip install mummichog==2.7.0 itself (needs network + python3).
#
# Output is permutation-based and stochastic, so we assert only that the
# expected files appear — never exact values.

test_that("pinned mummichog v2 runs end-to-end and emits the expected tree", {
  skip_on_cran()

  py <- Sys.getenv("MUMMICHOG_PYTHON", "")
  if (!nzchar(py) || !file.exists(py)) {
    if (Sys.getenv("MUMMICHOG_SMOKE_INSTALL") != "1") {
      skip("MUMMICHOG_PYTHON not set (set it, or MUMMICHOG_SMOKE_INSTALL=1 to build a temp venv)")
    }
    py3 <- Sys.which("python3")
    if (!nzchar(py3)) skip("python3 not available to build a temp venv")
    venv <- withr::local_tempdir()
    system2(py3, c("-m", "venv", venv), stdout = FALSE, stderr = FALSE)
    py <- file.path(venv, "bin", "python")
    if (!file.exists(py)) skip("could not create a temp venv")
    pip <- suppressWarnings(system2(
      py, c("-m", "pip", "install", "-q", "mummichog==2.7.0"),
      stdout = FALSE, stderr = FALSE
    ))
    if (!identical(as.integer(pip), 0L)) skip("could not pip install mummichog==2.7.0 (offline?)")
  }

  # Deterministic synthetic table (a few metabolite-like M+H masses so there is
  # something to match; no RNG, so no seeding needed).
  tbl <- data.frame(
    feature_id     = paste0("feat_", 1:12),
    mz             = c(176.1263, 181.0703, 147.0768, 205.0972, 167.0932,
                       134.0448, 269.0877, 348.0704, 666.1322, 260.0295,
                       300.1000, 411.2000),
    retention_time = c(199, 82.5, 55, 120, 210, 33, 145, 260, 90, 175, 60, 240),
    p_value        = c(0.004, 0.005, 0.70, 0.90, 0.002, 0.60,
                       0.80, 0.001, 0.95, 0.003, 0.40, 0.85),
    statistic      = c(-2.25, 1.77, 0.10, -0.30, 2.40, 0.05,
                       -0.20, 3.10, -0.10, 1.90, 0.20, -0.40),
    stringsAsFactors = FALSE
  )

  work    <- withr::local_tempdir()
  input   <- write_mummichog_input(tbl, file.path(work, "input.tsv"),
                                   id_col = "feature_id")
  files <- run_mummichog_v2(
    infile       = input,
    out_dir      = file.path(work, "v2"),
    project      = "smoke_test",
    python       = py,
    permutations = 20,          # small for speed; results are stochastic anyway
    cutoff       = 0.05
  )

  base <- basename(files)
  expect_true(any(base == "result.html"))
  expect_true(any(grepl("^mcg_pathwayanalysis.*\\.tsv$", base)))
  expect_true(any(base == "mummichog_manifest.tsv"))
  expect_s3_class(read_mummichog_pathways(files), "data.frame")
})
