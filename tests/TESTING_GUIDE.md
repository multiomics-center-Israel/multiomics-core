# Testing Guide: Metadata Profiler

## Quick Start (RStudio or R Console)

Since you're on Windows and R may not be in your system PATH, the easiest way to run tests is from RStudio or R console:

### Option 1: Run all profiler tests
```r
# From R console (working directory = project root)
testthat::test_file("tests/testthat/test-metadata_profile.R")
```

### Option 2: Use the test runner script
```r
# From R console
source("tests/run_profiler_tests.R")
```

### Option 3: Run tests with detailed output
```r
# See exactly which tests pass/fail
testthat::test_file(
  "tests/testthat/test-metadata_profile.R",
  reporter = "progress"
)
```

---

## Running Specific Test Suites

### Only run tests matching a pattern
```r
# Run only "simple" tests
testthat::test_file(
  "tests/testthat/test-metadata_profile.R",
  filter = "simple"
)

# Run only batch detection tests
testthat::test_file(
  "tests/testthat/test-metadata_profile.R",
  filter = "batch"
)

# Run only edge case tests
testthat::test_file(
  "tests/testthat/test-metadata_profile.R",
  filter = "edge"
)
```

---

## Interactive Testing (Recommended for Development)

When developing or debugging, test interactively:

```r
# 1. Source the profiler code
source("R/core/08_metadata_profile.R")

# 2. Create test metadata
meta <- data.frame(
  SampleID = paste0("S", 1:20),
  Condition = rep(c("Control", "Treated"), each = 10),
  Batch = rep(paste0("B", 1:2), 10),
  Age = rnorm(20, 45, 10),
  row.names = paste0("S", 1:20),
  stringsAsFactors = FALSE
)

# 3. Set up effects config
effects <- list(
  color = "Condition",
  samples = "SampleID"
)

# 4. Run profiler
prof <- assess_metadata_complexity(meta, effects)

# 5. Inspect results
print_metadata_profile(prof)

# 6. Check specific properties
prof$complexity_label        # "simple"
prof$n_groups_primary        # 2
prof$is_simple_2group        # TRUE
prof$has_batch_like          # TRUE
prof$annotation_cols         # "Condition", "Batch"
prof$plots$pca_1v3           # FALSE (simple design)
prof$plots$distance_heatmap  # TRUE (20 samples)
```

---

## Testing with Your Own Metadata

### Load your actual metadata
```r
# Source profiler
source("R/core/08_metadata_profile.R")

# Load your metadata (adjust path)
meta <- read.csv("data/rna/meta_rna.csv", row.names = 1)

# Define effects from your config
effects <- list(
  color = "YourConditionColumn",  # Adjust to your column name
  shape = "YourBatchColumn",      # Optional
  samples = "SampleID"
)

# Profile your metadata
prof <- assess_metadata_complexity(meta, effects)

# See what plots would be generated
print_metadata_profile(prof)

# Check specific decisions
prof$complexity_label
prof$annotation_cols
prof$plots
```

---

## Expected Test Results

When all tests pass, you should see:

```
test-metadata_profile.R
  ✓ assess_metadata_complexity returns expected structure
  ✓ profiler correctly identifies simple 2-group design
  ✓ profiler correctly identifies moderate complexity
  ✓ profiler correctly identifies complex design
  ✓ profiler detects batch-like columns
  ✓ profiler does not falsely detect batches in simple metadata
  ✓ profiler excludes ID-like columns from annotations
  ✓ annotation columns prioritize primary grouping variable
  ✓ annotation columns include batch columns
  ✓ annotation columns respect max_annotation_cols threshold
  ✓ annotation columns include categorical before numeric
  ✓ profiler suggests primary column for color
  ✓ profiler suggests shape for batch when appropriate
  ✓ profiler does not suggest shape for high-cardinality batch
  ✓ profiler respects configured shape column
  ✓ plot gating respects sample count thresholds
  ✓ heatmaps are skipped for large datasets
  ✓ core plots are always enabled
  ✓ profiler handles metadata with no categorical columns
  ✓ profiler handles missing color column gracefully
  ✓ profiler handles metadata with all ID columns
  ✓ profiler handles single-sample metadata
  ✓ custom thresholds override defaults
  ✓ partial threshold specification merges with defaults
  ✓ print_metadata_profile does not crash
  ✓ typical case-control RNA-seq experiment (SIMPLE)
  ✓ time-course experiment with batch effects (COMPLEX)
  ✓ multi-factor drug screen (VERY COMPLEX)

[ FAIL 0 | WARN 0 | SKIP 0 | PASS 28 ]
```

---

## Debugging Failed Tests

### If a test fails:

1. **Run with detailed reporter:**
```r
testthat::test_file(
  "tests/testthat/test-metadata_profile.R",
  reporter = "location"
)
```

2. **Run the specific failing test interactively:**
```r
# Look at test file, find the failing test code
# Copy the code and run it line by line

# Example:
meta <- data.frame(
  SampleID = paste0("S", 1:20),
  Condition = rep(c("A", "B"), each = 10),
  row.names = paste0("S", 1:20)
)
effects <- list(color = "Condition", samples = "SampleID")
prof <- assess_metadata_complexity(meta, effects)

# Now inspect what went wrong
str(prof)
prof$complexity_label  # What did it return vs. what was expected?
```

3. **Check profiler logic:**
```r
# Step through the profiler function
debug(assess_metadata_complexity)
prof <- assess_metadata_complexity(meta, effects)
# Use 'n' to step, 'c' to continue, 'Q' to quit
```

---

## Testing Different Threshold Configurations

```r
# Test with custom thresholds
meta <- make_simple_meta(50)  # 50 samples
effects <- list(color = "Condition", samples = "SampleID")

# Default thresholds (heatmaps enabled for ≤120 samples)
prof_default <- assess_metadata_complexity(meta, effects)
prof_default$plots$distance_heatmap  # TRUE

# Strict thresholds (heatmaps only for ≤30 samples)
thresholds_strict <- list(max_samples_for_heatmaps = 30)
prof_strict <- assess_metadata_complexity(meta, effects, thresholds_strict)
prof_strict$plots$distance_heatmap  # FALSE (50 > 30)

# Permissive thresholds
thresholds_permissive <- list(
  min_samples_for_pca3d = 5,
  max_samples_for_heatmaps = 200
)
prof_permissive <- assess_metadata_complexity(meta, effects, thresholds_permissive)
```

---

## Common Issues

### Issue: "could not find function 'assess_metadata_complexity'"
**Solution:** Source the profiler first
```r
source("R/core/08_metadata_profile.R")
```

### Issue: Tests fail because of missing packages
**Solution:** Install testthat
```r
install.packages("testthat")
```

### Issue: "helper.R sources all R files and causes errors"
**Solution:** Run tests from project root, not from tests/ directory
```r
setwd("path/to/multiomics-core")  # Project root
testthat::test_file("tests/testthat/test-metadata_profile.R")
```

---

## Next Steps After Tests Pass

Once all tests pass, you can:

1. **Test with your real data** (see "Testing with Your Own Metadata" above)
2. **Proceed to integration** (modify preprocessing and QC modules)
3. **Add custom thresholds** to your config files
4. **Add more tests** for your specific use cases

---

## Adding Your Own Tests

```r
# Add to tests/testthat/test-metadata_profile.R

test_that("my specific use case works", {

  # Create metadata matching your experiment
  meta <- data.frame(
    SampleID = paste0("S", 1:48),
    Treatment = rep(c("Vehicle", "LowDose", "HighDose"), each = 16),
    Timepoint = rep(c("0h", "6h", "24h", "48h"), each = 4, times = 3),
    Batch = rep(paste0("Batch", 1:3), 16),
    row.names = paste0("S", 1:48),
    stringsAsFactors = FALSE
  )

  effects <- list(color = "Treatment", samples = "SampleID")

  prof <- assess_metadata_complexity(meta, effects)

  # Assert your expectations
  expect_equal(prof$n_groups_primary, 3)
  expect_equal(prof$complexity_label, "complex")
  expect_true("Batch" %in% prof$annotation_cols)
  expect_true(prof$plots$pca_1v3)  # Complex design should enable this
})
```

Then run your new test:
```r
testthat::test_file(
  "tests/testthat/test-metadata_profile.R",
  filter = "my specific"
)
```

---

## Test Coverage Summary

- **28 unit tests** covering:
  - Basic functionality (5 tests)
  - Batch detection (2 tests)
  - ID exclusion (1 test)
  - Annotation selection (5 tests)
  - Aesthetic mapping (4 tests)
  - Plot gating (3 tests)
  - Edge cases (5 tests)
  - Threshold overrides (2 tests)
  - Print function (1 test)
  - Realistic scenarios (3 tests)

- **Test data generators** for quick mock metadata creation
- **Comprehensive documentation** in [PROFILER_TESTS.md](PROFILER_TESTS.md)
