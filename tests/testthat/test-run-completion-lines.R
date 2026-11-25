# tests/testthat/test-run-completion-lines.R
#
# run.R must not end on an unqualified "Pipeline complete." when a target
# failed without stopping tar_make() -- rna_bundle has error = "continue", and
# tar_make() returns normally when it fails. Checked without running a pipeline:
# pipeline_completion_lines() is evaluated out of run.R, which cannot be
# sourced because it ends by calling main().

run_r_path <- file.path(normalizePath(if (dir.exists("R")) "." else "../..",
                                      mustWork = FALSE), "run.R")

run_r_function <- function(name) {
  exprs <- parse(run_r_path, keep.source = FALSE, encoding = "UTF-8")
  env <- new.env(parent = globalenv())
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
        identical(as.character(e[[2]]), name)) {
      eval(e, env)
    }
  }
  get(name, envir = env, inherits = FALSE)
}


test_that("a clean run ends on 'Pipeline complete.'", {
  lines <- run_r_function("pipeline_completion_lines")
  expect_identical(lines(character(0)), "  Pipeline complete.")
})

test_that("a failed RNA bundle is reported even when the reports were written", {
  lines <- run_r_function("pipeline_completion_lines")
  out <- paste(lines("rna_bundle"), collapse = "\n")
  expect_match(out, "WARNING: RNA bundle export failed", fixed = TRUE)
  expect_match(out, "No RNA bundle is", fixed = TRUE)
  expect_match(out, "errored targets: rna_bundle", fixed = TRUE)
  expect_false(grepl("Pipeline complete.", out, fixed = TRUE))
})

test_that("other errored targets are listed without the bundle warning", {
  lines <- run_r_function("pipeline_completion_lines")
  out <- paste(lines(c("b_target", "a_target", "a_target")), collapse = "\n")
  expect_match(out, "errored targets: a_target, b_target", fixed = TRUE)
  expect_false(grepl("RNA bundle", out, fixed = TRUE))
  expect_false(grepl("Pipeline complete.", out, fixed = TRUE))
})

test_that("run_pipeline() ends on the errored-target check, not a fixed message", {
  body_text <- paste(deparse(body(run_r_function("run_pipeline"))), collapse = "\n")
  expect_match(body_text, "tar_errored(store = store)", fixed = TRUE)
  expect_match(body_text, "pipeline_completion_lines(errored)", fixed = TRUE)
  expect_false(grepl("Pipeline complete.", body_text, fixed = TRUE))
})
