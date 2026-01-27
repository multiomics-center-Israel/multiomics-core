# Diagnostic Script - Print Actual Function Content

cat("=== DIAGNOSTIC: Checking what code is actually loaded ===\n\n")

# Step 1: Source the files
cat("Step 1: Sourcing files...\n")
source("R/core/06_shiny_export.R")
source("R/domain/proteomics/06_shiny_export.R")

# Step 2: Get the function body
cat("\nStep 2: Extracting function body...\n")
func_body <- deparse(body(build_data_to_shiny_legacy_proteomics))

# Step 3: Find the expected_keys section
cat("\nStep 3: Looking for expected_keys definition...\n")
expected_keys_start <- which(grepl("expected_keys.*<-.*c\\(", func_body))

if (length(expected_keys_start) > 0) {
    cat("Found expected_keys at line", expected_keys_start, "of function body\n")
    cat("\nShowing 20 lines starting from expected_keys definition:\n")
    cat("---\n")
    for (i in expected_keys_start:(min(expected_keys_start + 19, length(func_body)))) {
        cat(sprintf("%3d: %s\n", i, func_body[i]))
    }
    cat("---\n")
} else {
    cat("ERROR: Could not find expected_keys definition!\n")
    cat("\nShowing first 50 lines of function:\n")
    cat(paste(head(func_body, 50), collapse = "\n"), "\n")
}

# Step 4: Check specific patterns
cat("\nStep 4: Checking for specific patterns...\n")
has_color_shape <- any(grepl("color", func_body) & grepl("shape", func_body))
has_pca_fields <- any(grepl("norm_log_counts_pca", func_body) & grepl("mat2plot", func_body))

cat("- Contains 'color' AND 'shape':", has_color_shape, "\n")
cat("- Contains 'norm_log_counts_pca' AND 'mat2plot':", has_pca_fields, "\n")

cat("\n=== END DIAGNOSTIC ===\n")
