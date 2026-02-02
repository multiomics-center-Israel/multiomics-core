# Metadata Profiler Unit Tests

Comprehensive test suite for the metadata complexity profiler (`R/core/08_metadata_profile.R`).

## Quick Start

### Run all profiler tests
```r
# From R console
testthat::test_file("tests/testthat/test-metadata_profile.R")

# Or use the runner script
source("tests/run_profiler_tests.R")
```

### Run from command line
```bash
Rscript tests/run_profiler_tests.R
```

### Run specific test
```r
testthat::test_file(
  "tests/testthat/test-metadata_profile.R",
  filter = "simple 2-group"  # runs only tests matching this pattern
)
```

---

## Test Coverage

### 1. **Basic Functionality** (5 tests)
- ✅ Returns expected data structure
- ✅ Identifies simple 2-group designs
- ✅ Identifies moderate complexity (3-5 groups)
- ✅ Identifies complex designs (>5 groups, multi-factor)
- ✅ Validates all required output fields

### 2. **Batch Detection** (2 tests)
- ✅ Detects batch/lane/run/plate columns via pattern matching
- ✅ Does not falsely flag non-batch columns

### 3. **ID Column Exclusion** (1 test)
- ✅ Excludes high-uniqueness columns (>90% unique values)
- ✅ Preserves low-cardinality categorical columns

### 4. **Annotation Column Selection** (5 tests)
- ✅ Prioritizes primary grouping variable
- ✅ Includes batch columns
- ✅ Respects `max_annotation_cols` threshold
- ✅ Orders categorical before numeric
- ✅ Excludes sample ID and other ID-like columns

### 5. **Aesthetic Mapping Suggestions** (4 tests)
- ✅ Suggests primary column for color
- ✅ Suggests batch for shape when low cardinality
- ✅ Skips shape suggestion for high-cardinality columns
- ✅ Respects user-configured shape column

### 6. **Plot Gating Decisions** (3 tests)
- ✅ Respects sample count thresholds (3D PCA, heatmaps)
- ✅ Skips expensive heatmaps for large datasets (>120 samples)
- ✅ Always enables core plots (PCA 1v2, density, boxplot)

### 7. **Edge Cases** (5 tests)
- ✅ Handles metadata with no categorical columns
- ✅ Handles missing color column (falls back to first categorical)
- ✅ Handles all-ID columns gracefully
- ✅ Handles single-sample metadata
- ✅ Handles empty or minimal metadata

### 8. **Threshold Overrides** (2 tests)
- ✅ Custom thresholds override defaults
- ✅ Partial threshold specification merges with defaults

### 9. **Print Function** (1 test)
- ✅ `print_metadata_profile()` executes without errors

### 10. **Realistic RNA-seq Scenarios** (3 tests)
- ✅ **Simple**: Case-control study (12 samples, 2 groups)
- ✅ **Moderate**: Time-course with batch effects (24 samples, 2×3 factorial)
- ✅ **Complex**: Multi-factor drug screen (96 samples, 8 cell lines × 4 drugs)

---

## Test Data Generators

The test suite includes 5 mock metadata generators:

| Generator | Use Case | Complexity |
|-----------|----------|------------|
| `make_simple_meta()` | 2-group comparison | Simple |
| `make_moderate_meta()` | 3-5 groups, 2-3 factors | Moderate |
| `make_complex_meta()` | >5 groups, batch structure | Complex |
| `make_meta_with_ids()` | Tests ID column exclusion | N/A |
| `make_minimal_meta()` | Numeric-only metadata | Edge case |

---

## Expected Test Output

```
Running Metadata Profiler Unit Tests
========================================

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

✅ ALL TESTS PASSED (28 tests)
```

---

## Adding New Tests

### Template for new test
```r
test_that("descriptive test name", {

  # 1. Setup mock metadata
  meta <- data.frame(
    SampleID = paste0("S", 1:20),
    YourColumn = ...,
    row.names = paste0("S", 1:20)
  )
  effects <- list(color = "YourColumn", samples = "SampleID")

  # 2. Run profiler
  prof <- assess_metadata_complexity(meta, effects)

  # 3. Assert expectations
  expect_equal(prof$n_samples, 20)
  expect_true(prof$some_property)
  expect_false(prof$other_property)
})
```

### Common assertions
```r
# Structure checks
expect_type(prof, "list")
expect_true("plots" %in% names(prof))

# Complexity checks
expect_equal(prof$complexity_label, "simple")
expect_true(prof$is_simple_2group)
expect_equal(prof$n_groups_primary, 2)

# Annotation checks
expect_true("Condition" %in% prof$annotation_cols)
expect_length(prof$annotation_cols, 3)

# Plot gating checks
expect_true(prof$plots$pca_1v2)
expect_false(prof$plots$expr_heatmap)

# Threshold checks
expect_equal(prof$thresholds$max_samples_for_heatmaps, 120)
```

---

## Debugging Failed Tests

### Run with verbose output
```r
testthat::test_file(
  "tests/testthat/test-metadata_profile.R",
  reporter = "location"  # Shows exact line numbers
)
```

### Test interactively
```r
# Source the profiler
source("R/core/08_metadata_profile.R")

# Create test data
meta <- data.frame(
  SampleID = paste0("S", 1:20),
  Condition = rep(c("A", "B"), each = 10),
  row.names = paste0("S", 1:20)
)
effects <- list(color = "Condition", samples = "SampleID")

# Run and inspect
prof <- assess_metadata_complexity(meta, effects)
str(prof)  # Inspect structure
print_metadata_profile(prof)  # Pretty print
```

### Check specific component
```r
# Test annotation selection
prof$annotation_cols
# Expected: "Condition"

# Test plot gating
prof$plots$distance_heatmap
# Expected: TRUE (20 samples < 120 threshold)
```

---

## Integration with CI/CD

### GitHub Actions example
```yaml
- name: Run profiler tests
  run: Rscript tests/run_profiler_tests.R
```

### Local pre-commit hook
```bash
#!/bin/bash
# .git/hooks/pre-commit
Rscript tests/run_profiler_tests.R || exit 1
```

---

## Maintenance Notes

- **Add tests when**: New complexity rules, new thresholds, new edge cases discovered
- **Update tests when**: Profiler logic changes, new metadata patterns emerge
- **Remove tests when**: Deprecated functionality removed (update this doc!)

---

## Questions or Issues?

If tests fail unexpectedly:
1. Check if profiler logic changed
2. Verify test assumptions are still valid
3. Inspect failed test output for clues
4. Run test interactively to debug

If you need to modify profiler behavior:
1. Write failing test first (TDD)
2. Update profiler code
3. Verify all tests pass
4. Update this documentation
