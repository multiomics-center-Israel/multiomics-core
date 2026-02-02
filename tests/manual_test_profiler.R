#!/usr/bin/env Rscript
# Manual test of metadata profiler (no testthat dependency)

# Source the profiler
cat("Loading profiler code...\n")
source("R/core/08_metadata_profile.R")

# Track test results
tests_passed <- 0
tests_failed <- 0
failed_tests <- character(0)

# Helper function to run a test
run_test <- function(test_name, condition, expected = TRUE) {
  if (identical(condition, expected)) {
    tests_passed <<- tests_passed + 1
    cat(sprintf("✓ %s\n", test_name))
    return(TRUE)
  } else {
    tests_failed <<- tests_failed + 1
    failed_tests <<- c(failed_tests, test_name)
    cat(sprintf("✗ %s (Expected: %s, Got: %s)\n",
                test_name, expected, condition))
    return(FALSE)
  }
}

cat("\n========================================\n")
cat("MANUAL PROFILER TESTS\n")
cat("========================================\n\n")

# ============================================================================
# Test 1: Simple 2-group design
# ============================================================================
cat("Test Group 1: Simple 2-group design\n")
cat("────────────────────────────────────\n")

meta_simple <- data.frame(
  SampleID = paste0("S", 1:20),
  Treatment = rep(c("Control", "Treated"), each = 10),
  Age = rnorm(20, mean = 45, sd = 10),
  row.names = paste0("S", 1:20),
  stringsAsFactors = FALSE
)
effects_simple <- list(color = "Treatment", samples = "SampleID")
prof_simple <- assess_metadata_complexity(meta_simple, effects_simple)

run_test("Returns list", is.list(prof_simple))
run_test("Has n_samples field", "n_samples" %in% names(prof_simple))
run_test("Correct sample count", prof_simple$n_samples == 20)
run_test("Identifies as simple", prof_simple$complexity_label == "simple")
run_test("Detects 2 groups", prof_simple$n_groups_primary == 2)
run_test("Flags as 2-group", prof_simple$is_simple_2group == TRUE)
run_test("Primary column correct", prof_simple$primary_col == "Treatment")
run_test("Core plots enabled", prof_simple$plots$pca_1v2 == TRUE)
run_test("Complex plots disabled", prof_simple$plots$pca_1v3 == FALSE)

cat("\n")

# ============================================================================
# Test 2: Complex design with batches
# ============================================================================
cat("Test Group 2: Complex design with batches\n")
cat("──────────────────────────────────────────\n")

meta_complex <- data.frame(
  SampleID = paste0("S", 1:60),
  CellLine = rep(paste0("CL", 1:6), each = 10),
  Treatment = rep(c("DMSO", "Drug1", "Drug2", "Drug3"), length.out = 60),
  Batch = rep(paste0("Batch", 1:3), length.out = 60),
  PlateID = rep(paste0("P", 1:10), length.out = 60),
  row.names = paste0("S", 1:60),
  stringsAsFactors = FALSE
)
effects_complex <- list(color = "CellLine", samples = "SampleID")
prof_complex <- assess_metadata_complexity(meta_complex, effects_complex)

run_test("Complex: 60 samples", prof_complex$n_samples == 60)
run_test("Complex: 6 groups", prof_complex$n_groups_primary == 6)
run_test("Complex: labeled correctly", prof_complex$complexity_label == "complex")
run_test("Complex: detects batches", prof_complex$has_batch_like == TRUE)
run_test("Complex: Batch in batch_cols", "Batch" %in% prof_complex$batch_cols)
run_test("Complex: PlateID in batch_cols", "PlateID" %in% prof_complex$batch_cols)
run_test("Complex: enables pca_1v3", prof_complex$plots$pca_1v3 == TRUE)
run_test("Complex: suggests shape", !is.null(prof_complex$suggested_mappings$shape))
run_test("Complex: primary in annotations", prof_complex$annotation_cols[1] == "CellLine")
run_test("Complex: batch in annotations", "Batch" %in% prof_complex$annotation_cols)

cat("\n")

# ============================================================================
# Test 3: ID exclusion
# ============================================================================
cat("Test Group 3: ID column exclusion\n")
cat("──────────────────────────────────\n")

meta_ids <- data.frame(
  SampleID = paste0("S", 1:20),
  BarcodeID = paste0("BC-", sprintf("%05d", 1:20)),  # Unique
  Condition = rep(c("A", "B"), each = 10),
  PatientID = paste0("PT", 1:20),  # Unique
  row.names = paste0("S", 1:20),
  stringsAsFactors = FALSE
)
effects_ids <- list(color = "Condition", samples = "SampleID")
prof_ids <- assess_metadata_complexity(meta_ids, effects_ids)

run_test("IDs: BarcodeID excluded", !("BarcodeID" %in% prof_ids$annotation_cols))
run_test("IDs: PatientID excluded", !("PatientID" %in% prof_ids$annotation_cols))
run_test("IDs: SampleID excluded", !("SampleID" %in% prof_ids$annotation_cols))
run_test("IDs: Condition included", "Condition" %in% prof_ids$annotation_cols)

cat("\n")

# ============================================================================
# Test 4: Plot gating by sample count
# ============================================================================
cat("Test Group 4: Plot gating thresholds\n")
cat("─────────────────────────────────────\n")

# Small dataset (< 10 samples)
meta_small <- data.frame(
  SampleID = paste0("S", 1:8),
  Condition = rep(c("A", "B"), each = 4),
  row.names = paste0("S", 1:8),
  stringsAsFactors = FALSE
)
effects_small <- list(color = "Condition", samples = "SampleID")
prof_small <- assess_metadata_complexity(meta_small, effects_small)

run_test("Small: 8 samples", prof_small$n_samples == 8)
run_test("Small: skips 3D PCA", prof_small$plots$pca_3d == FALSE)
run_test("Small: enables core plots", prof_small$plots$pca_1v2 == TRUE)

# Large dataset (> 120 samples)
meta_large <- data.frame(
  SampleID = paste0("S", 1:150),
  Condition = rep(c("A", "B"), length.out = 150),
  row.names = paste0("S", 1:150),
  stringsAsFactors = FALSE
)
effects_large <- list(color = "Condition", samples = "SampleID")
prof_large <- assess_metadata_complexity(meta_large, effects_large)

run_test("Large: 150 samples", prof_large$n_samples == 150)
run_test("Large: skips distance heatmap", prof_large$plots$distance_heatmap == FALSE)
run_test("Large: skips correlation heatmap", prof_large$plots$correlation_heatmap == FALSE)
run_test("Large: skips expr heatmap", prof_large$plots$expr_heatmap == FALSE)

cat("\n")

# ============================================================================
# Test 5: Threshold overrides
# ============================================================================
cat("Test Group 5: Custom thresholds\n")
cat("────────────────────────────────\n")

meta_thresh <- data.frame(
  SampleID = paste0("S", 1:50),
  Condition = rep(c("A", "B"), each = 25),
  row.names = paste0("S", 1:50),
  stringsAsFactors = FALSE
)
effects_thresh <- list(color = "Condition", samples = "SampleID")

# Default thresholds
prof_default <- assess_metadata_complexity(meta_thresh, effects_thresh)
run_test("Default: heatmap enabled", prof_default$plots$distance_heatmap == TRUE)

# Custom strict thresholds
thresholds_strict <- list(max_samples_for_heatmaps = 30)
prof_strict <- assess_metadata_complexity(meta_thresh, effects_thresh, thresholds_strict)
run_test("Strict: heatmap disabled", prof_strict$plots$distance_heatmap == FALSE)
run_test("Strict: threshold stored", prof_strict$thresholds$max_samples_for_heatmaps == 30)

cat("\n")

# ============================================================================
# Test 6: Edge cases
# ============================================================================
cat("Test Group 6: Edge cases\n")
cat("─────────────────────────\n")

# No categorical columns
meta_nocats <- data.frame(
  SampleID = paste0("S", 1:15),
  Value1 = rnorm(15),
  Value2 = rnorm(15),
  row.names = paste0("S", 1:15),
  stringsAsFactors = FALSE
)
effects_nocats <- list(samples = "SampleID")
prof_nocats <- assess_metadata_complexity(meta_nocats, effects_nocats)

run_test("NoCats: returns list", is.list(prof_nocats))
run_test("NoCats: n_groups is NA", is.na(prof_nocats$n_groups_primary))
run_test("NoCats: not simple 2-group", prof_nocats$is_simple_2group == FALSE)

cat("\n")

# ============================================================================
# Test 7: Print function
# ============================================================================
cat("Test Group 7: Print function\n")
cat("─────────────────────────────\n")

# Should not crash
print_result <- tryCatch({
  capture.output(print_metadata_profile(prof_simple))
  TRUE
}, error = function(e) {
  FALSE
})

run_test("Print: no errors", print_result == TRUE)

cat("\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat("========================================\n")
cat("TEST SUMMARY\n")
cat("========================================\n\n")

total_tests <- tests_passed + tests_failed
cat(sprintf("Total tests:  %d\n", total_tests))
cat(sprintf("Passed:       %d\n", tests_passed))
cat(sprintf("Failed:       %d\n", tests_failed))

if (tests_failed > 0) {
  cat("\n❌ Failed tests:\n")
  for (test in failed_tests) {
    cat(sprintf("   - %s\n", test))
  }
  cat("\n")
  quit(status = 1)
} else {
  cat("\n✅ ALL TESTS PASSED!\n\n")

  # Show example output
  cat("Example profiler output:\n")
  cat("────────────────────────────────────\n")
  print_metadata_profile(prof_simple)

  quit(status = 0)
}
