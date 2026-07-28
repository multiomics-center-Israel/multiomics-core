# tests/testthat/test-mummichog-force-primary-ion.R
#
# Unit tests for the -z / force_primary_ion CLI wiring (06c). Exercises the pure
# arg builder .mmc_build_cli_args() directly — no live mummichog run needed.
#
# Ground truth (mummichog 2.7.0, get_user_data.py): force_primary_ion defaults
# to TRUE and `-z` takes a VALUE (getopt "z:"), mapped via a booleandict; the
# long --force_primary_ion is broken. So we emit `-z True` / `-z False` (never a
# bare -z), and only when explicitly set.

base_args <- function(...) {
  .mmc_build_cli_args(
    infile = "in.tsv", project = "proj", network = "human_mfn",
    mode = "pos_default", instrument_ppm = 10, permutations = 100, ...)
}

# Index of the value following the first "-z" flag, or NA if absent.
z_value <- function(args) {
  i <- which(args == "-z")
  if (length(i) == 0) NA_character_ else args[i[[1]] + 1L]
}

test_that("no -z is emitted when force_primary_ion is unset (NULL)", {
  args <- base_args()                       # force_primary_ion defaults to NULL
  expect_false("-z" %in% args)
  # sanity: the core flags are all present
  expect_true(all(c("-f", "in.tsv", "-n", "human_mfn", "-m", "pos_default") %in% args))
})

test_that("force_primary_ion = FALSE emits `-z False`", {
  args <- base_args(force_primary_ion = FALSE)
  expect_true("-z" %in% args)
  expect_identical(z_value(args), "False")
})

test_that("force_primary_ion = TRUE emits `-z True`", {
  args <- base_args(force_primary_ion = TRUE)
  expect_identical(z_value(args), "True")
})

test_that("quoted / string boolean values coerce (YAML may pass a string)", {
  expect_identical(z_value(base_args(force_primary_ion = "true")),  "True")
  expect_identical(z_value(base_args(force_primary_ion = "false")), "False")
  expect_identical(z_value(base_args(force_primary_ion = "no")),    "False")
})

test_that("a non-logical force_primary_ion errors loudly (never silent -z False)", {
  expect_error(base_args(force_primary_ion = "ture"), "force_primary_ion")
  expect_error(base_args(force_primary_ion = "maybe"), "force_primary_ion")
})

test_that("cutoff and force_primary_ion coexist, extra_args stay last", {
  args <- base_args(cutoff = 0.05, force_primary_ion = FALSE,
                    extra_args = c("--flag", "x"))
  # -c <cutoff> present
  ci <- which(args == "-c")
  expect_identical(args[ci + 1L], "0.05")
  # -z False present
  expect_identical(z_value(args), "False")
  # extra_args appended verbatim at the very end
  expect_identical(tail(args, 2), c("--flag", "x"))
})
