# Simple test runner
library(testthat)

# Run profiler tests
cat("Running metadata profiler tests...\n\n")

results <- test_file(
  "tests/testthat/test-metadata_profile.R",
  reporter = "progress"
)

cat("\n\n")
cat("========================================\n")
cat("SUMMARY\n")
cat("========================================\n")
summary(results)
