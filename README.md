# multiomics-core

A modular and reproducible R framework for single-omics and multi-omics analyses\
(RNA-seq, proteomics, metabolomics).

This repository provides: - standardized data loading and validation - modular preprocessing (filtering, normalization, imputation) - unified QC (including PCA and heatmaps) - a configuration-driven workflow - reproducibility via `renv` and (future) `{targets}` pipelines

------------------------------------------------------------------------

## Repository structure

R/ ├── core/ \# Small reusable helpers (PCA utils, metadata alignment, etc.) ├── plots/ \# plot\_\* functions (pure plotting, no saving, no config) ├── qc/ \# qc\_\* wrappers (QC logic + saving) ├── preprocess\_\* \# Omics-specific preprocessing ├── pipeline\_\* \# High-level pipeline entry points config/ ├── config.yml \# Central configuration file outputs/ \# Analysis outputs (ignored by git)

yaml Copy code

------------------------------------------------------------------------

## Requirements

-   R (\>= 4.3 recommended)
-   `renv`
-   System libraries required by Bioconductor (if using RNA/proteomics)

------------------------------------------------------------------------

## Setup

Clone the repository and restore the R environment:

``` r
install.packages("renv") 
renv::restore() 
```

This will install all required CRAN and Bioconductor packages. Configuration

Edit the configuration file:

```         
config/config.yml
```

This file controls:

file paths

omics-specific parameters

filtering / normalization / imputation settings

QC aesthetics (color, shape, sample ID columns)

# Running proteomics QC (example)

``` r
# Load all functions
r_files <- list.files("R", pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE)
invisible(lapply(sort(r_files), source))

# Load config
config <- load_config("config/config.yml")

# Run proteomics pipeline
res <- run_proteomics_preprocessing(config)

# Example QC: PCA
qc_pca_scatter(
  expr_mat = res$expr_imp,
  meta     = res$meta,
  cfg      = config$modes$proteomics,
  out_file = "outputs/proteomics/qc/pca_pc1_pc2.png"
)
```

# Reproducibility

All package versions are locked in renv.lock

Outputs, caches, and intermediate files are excluded from git

Future versions will use {targets} for fully declarative pipelines

# Status

Version: v0.1.0 This release focuses on:

proteomics preprocessing

QC refactoring and PCA unification

RNA and multi-omics integration will be added in later versions.
