# tests/testthat/test-renviron.R
#
# Unit tests for set_renviron_var() (R/core/18_renviron.R): the idempotent
# .Renviron upsert used by `make setup`.

test_that("it appends the variable to a fresh file", {
  f <- withr::local_tempfile()
  expect_identical(set_renviron_var("MUMMICHOG_PYTHON", "/x/python", f), "added")
  expect_identical(readLines(f), "MUMMICHOG_PYTHON=/x/python")
})

test_that("it is idempotent — a second identical call changes nothing", {
  f <- withr::local_tempfile()
  set_renviron_var("MUMMICHOG_PYTHON", "/x/python", f)
  expect_identical(set_renviron_var("MUMMICHOG_PYTHON", "/x/python", f), "unchanged")
  # exactly one assignment, no duplicate
  expect_identical(readLines(f), "MUMMICHOG_PYTHON=/x/python")
})

test_that("it updates the value in place without duplicating the line", {
  f <- withr::local_tempfile()
  set_renviron_var("MUMMICHOG_PYTHON", "/old/python", f)
  expect_identical(set_renviron_var("MUMMICHOG_PYTHON", "/new/python", f), "updated")
  expect_identical(readLines(f), "MUMMICHOG_PYTHON=/new/python")
})

test_that("it preserves unrelated lines (never clobbers other secrets)", {
  f <- withr::local_tempfile()
  writeLines(c("SECRET_KEY=abc123", "OTHER=1"), f)
  set_renviron_var("MUMMICHOG_PYTHON", "/x/python", f)
  expect_identical(
    readLines(f),
    c("SECRET_KEY=abc123", "OTHER=1", "MUMMICHOG_PYTHON=/x/python")
  )
})

test_that("it replaces an existing assignment among other lines, keeping order", {
  f <- withr::local_tempfile()
  writeLines(c("SECRET_KEY=abc123", "MUMMICHOG_PYTHON=/old/python", "OTHER=1"), f)
  set_renviron_var("MUMMICHOG_PYTHON", "/new/python", f)
  expect_identical(
    readLines(f),
    c("SECRET_KEY=abc123", "MUMMICHOG_PYTHON=/new/python", "OTHER=1")
  )
})

test_that("it tolerates an `export ` prefix and collapses duplicates", {
  f <- withr::local_tempfile()
  writeLines(c("export MUMMICHOG_PYTHON=/a", "KEEP=1", "MUMMICHOG_PYTHON=/b"), f)
  set_renviron_var("MUMMICHOG_PYTHON", "/new/python", f)
  # first (export) line replaced, later duplicate dropped, KEEP preserved
  expect_identical(
    readLines(f),
    c("MUMMICHOG_PYTHON=/new/python", "KEEP=1")
  )
})

test_that("it rejects an invalid variable name", {
  f <- withr::local_tempfile()
  expect_error(set_renviron_var("BAD NAME", "x", f), "Invalid")
  expect_false(file.exists(f))   # nothing written on a rejected key
})
