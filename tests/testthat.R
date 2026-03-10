# Standard testthat runner for multiomics-core
# Run from project root: Rscript tests/testthat.R

if (!requireNamespace("testthat", quietly = TRUE)) {
    stop("testthat is required to run tests. Install with: install.packages('testthat')")
}

library(testthat)

test_dir(
    "tests/testthat",
    reporter = "summary",
    stop_on_failure = TRUE
)
