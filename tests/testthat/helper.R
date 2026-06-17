# tests/testthat/helper.R

# This script runs before tests.
# Since this project is not a formal package with a DESCRIPTION file,
# we need to manually source the R files to make functions available to tests.

# Load libraries commonly used
suppressPackageStartupMessages({
    library(testthat)
    library(dplyr)
    library(limma)
    library(edgeR)
})

# Source all R files (similar to _targets.R approach)
# We assume the working directory is the project root or tests/testthat
# standardized way to find root:
root_dir <- normalizePath(if (dir.exists("R")) "." else "../..")

r_files <- list.files(file.path(root_dir, "R"), pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE)

# Filter out pipeline files if they conflict or are not needed for unit tests
# (Optional: exclusion logic)

for (f in r_files) {
    # We wrap in tryCatch to avoid stopping if one file fails (e.g. syntax error or missing dep)
    tryCatch(source(f), error = function(e) {
        message("Error sourcing ", f, ": ", e$message)
    })
}
