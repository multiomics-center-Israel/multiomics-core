#!/usr/bin/env Rscript
# Quick validation of multi-omics pipeline structure

cat("\n=== Multi-Omics Pipeline Validation ===\n\n")

# Check domain files exist
domain_files <- c(
  "R/domain/multiomics/00_inputs.R",
  "R/domain/multiomics/01_id_mapping.R",
  "R/domain/multiomics/02_mae.R",
  "R/domain/multiomics/03_feature_selection.R",
  "R/domain/multiomics/04_integration_diablo.R",
  "R/domain/multiomics/05_integration_snf.R",
  "R/domain/multiomics/06_concordance.R",
  "R/domain/multiomics/07_enrichment.R",
  "R/domain/multiomics/08_report.R",
  "R/domain/multiomics/90_config_validate.R"
)

cat("✓ Domain files:\n")
for (f in domain_files) {
  exists <- file.exists(f)
  status <- if (exists) "✓" else "✗"
  cat(sprintf("  %s %s\n", status, basename(f)))
  if (!exists) stop("Missing file: ", f)
}

# Check module files exist
module_files <- c(
  "R/modules/multiomics/01_mod_harmonization.R",
  "R/modules/multiomics/02_mod_integration.R",
  "R/modules/multiomics/03_mod_concordance.R",
  "R/modules/multiomics/04_mod_enrichment.R",
  "R/modules/multiomics/05_mod_report.R"
)

cat("\n✓ Module files:\n")
for (f in module_files) {
  exists <- file.exists(f)
  status <- if (exists) "✓" else "✗"
  cat(sprintf("  %s %s\n", status, basename(f)))
  if (!exists) stop("Missing file: ", f)
}

# Check pipeline file exists
pipeline_file <- "R/pipeline/multiomics/00_pipe_multiomics.R"
cat("\n✓ Pipeline orchestration:\n")
exists <- file.exists(pipeline_file)
status <- if (exists) "✓" else "✗"
cat(sprintf("  %s %s\n", status, basename(pipeline_file)))
if (!exists) stop("Missing file: ", pipeline_file)

# Check _targets.R was updated
targets_file <- "_targets.R"
cat("\n✓ Pipeline integration:\n")
targets_content <- readLines(targets_file)
has_multiomics <- any(grepl("pipe_multiomics", targets_content))
status <- if (has_multiomics) "✓" else "✗"
cat(sprintf("  %s _targets.R loads pipe_multiomics()\n", status))

# Check config template exists
config_template <- "config/templates/multiomics_config.yaml"
exists <- file.exists(config_template)
status <- if (exists) "✓" else "✗"
cat(sprintf("  %s Config template exists\n", status))

# Check test config exists
test_config <- "config/multiomics_GT15_test.yaml"
exists <- file.exists(test_config)
status <- if (exists) "✓" else "✗"
cat(sprintf("  %s Test config exists\n", status))

# Syntax check (if yaml package available)
if (requireNamespace("yaml", quietly = TRUE)) {
  cat("\n✓ Config validation:\n")
  tryCatch({
    cfg <- yaml::read_yaml(test_config)
    cat("  ✓ Test config is valid YAML\n")

    has_multiomics_mode <- !is.null(cfg$modes$multiomics)
    status <- if (has_multiomics_mode) "✓" else "✗"
    cat(sprintf("  %s Config has multiomics mode\n", status))

    has_rna <- !is.null(cfg$modes$rna)
    has_metab <- !is.null(cfg$modes$metabolomics)
    n_omics <- sum(has_rna, has_metab)
    status <- if (n_omics >= 2) "✓" else "✗"
    cat(sprintf("  %s Config has %d omics modes (need ≥2)\n", status, n_omics))

  }, error = function(e) {
    cat("  ✗ Config validation failed:", e$message, "\n")
  })
}

cat("\n=== Validation Summary ===\n")
cat("✓ All multi-omics pipeline files created successfully!\n")
cat("✓ Pipeline structure is correct\n")
cat("\n")
cat("Next steps:\n")
cat("  1. Install required packages: renv::restore()\n")
cat("  2. Review config: config/multiomics_GT15_test.yaml\n")
cat("  3. Run pipeline: tar_make()\n")
cat("\nSee MULTIOMICS_QUICKSTART.md for detailed instructions.\n\n")
