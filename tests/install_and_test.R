# Install testthat if needed and run tests

# Check if testthat is available
if (!requireNamespace("testthat", quietly = TRUE)) {
  cat("Installing testthat...\n")
  install.packages("testthat", repos = "https://cloud.r-project.org", quiet = TRUE)
}

library(testthat)

# Source the profiler code directly
cat("Sourcing profiler code...\n")
source("R/core/08_metadata_profile.R")

# Run tests
cat("\n========================================\n")
cat("Running Metadata Profiler Unit Tests\n")
cat("========================================\n\n")

results <- test_file(
  "tests/testthat/test-metadata_profile.R",
  reporter = "progress"
)

# Summary
cat("\n\n========================================\n")
cat("TEST SUMMARY\n")
cat("========================================\n")
print(as.data.frame(results))

# Check if all passed
n_failed <- sum(as.data.frame(results)$failed)
n_passed <- sum(as.data.frame(results)$passed)

cat(sprintf("\nPassed: %d\n", n_passed))
cat(sprintf("Failed: %d\n", n_failed))

if (n_failed > 0) {
  cat("\n❌ SOME TESTS FAILED\n")
  quit(status = 1)
} else {
  cat("\n✅ ALL TESTS PASSED\n")
  quit(status = 0)
}
