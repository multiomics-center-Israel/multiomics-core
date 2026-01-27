# Force Complete Reload - Run this in R console

# Step 1: Completely destroy targets cache
targets::tar_destroy()

# Step 2: Clear R environment
rm(list = ls())

# Step 3: Manually source the updated files to verify they load
source("R/core/06_shiny_export.R")
source("R/domain/proteomics/06_shiny_export.R")

# Step 4: Verify the function is correct
# This should show "color", "shape" in expected_keys, NOT "norm_log_counts_pca", "mat2plot"
cat("Checking expected_keys in build_data_to_shiny_legacy_proteomics...\n")
func_body <- deparse(body(build_data_to_shiny_legacy_proteomics))
if (any(grepl('"color", "shape"', func_body, fixed = TRUE))) {
    cat("✓ CORRECT: Found 'color', 'shape' in expected_keys\n")
} else {
    cat("✗ ERROR: 'color', 'shape' NOT found in expected_keys\n")
}

if (any(grepl('"norm_log_counts_pca", "mat2plot"', func_body, fixed = TRUE))) {
    cat("✗ ERROR: OLD CODE STILL LOADED - found 'norm_log_counts_pca', 'mat2plot' in expected_keys\n")
} else {
    cat("✓ CORRECT: 'norm_log_counts_pca', 'mat2plot' NOT in expected_keys (they are optional)\n")
}

# Step 5: Run targets
cat("\nRunning tar_make()...\n")
targets::tar_make()
