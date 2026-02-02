#!/usr/bin/env Rscript
# Run metadata profiler unit tests
#
# Usage:
#   Rscript tests/run_profiler_tests.R
#
# Or from R console:
#   source("tests/run_profiler_tests.R")

library(testthat)

# Set working directory to project root if running as script
if (basename(getwd()) == "tests") {
  setwd("..")
}

# Run only profiler tests
cat("\n========================================\n")
cat("Running Metadata Profiler Unit Tests\n")
cat("========================================\n\n")

test_results <- test_file(
  "tests/testthat/test-metadata_profile.R",
  reporter = "progress"
)

# Print summary
cat("\n========================================\n")
cat("Test Summary:\n")
cat("========================================\n")
print(test_results)

# Exit with appropriate code
if (any(as.data.frame(test_results)$failed > 0)) {
  cat("\n❌ TESTS FAILED\n\n")
  quit(status = 1)
} else {
  cat("\n✅ ALL TESTS PASSED\n\n")
  quit(status = 0)
}
