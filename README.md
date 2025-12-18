# multiomics-core

A modular and reproducible R framework for single-omics and multi-omics analyses (RNA-seq, proteomics, metabolomics, lipidomics).

### The project emphasizes:

-   clear separation of concerns (I/O, preprocessing, DE, QC, plotting)
-   configuration-driven workflows - reproducibility via `renv`
-   scalable orchestration via `{targets}`

------------------------------------------------------------------------

### What this repository provides

-   Standardized data loading and validation
-   Omics-specific preprocessing (filtering, normalization, imputation)
-   Proteomics differential expression using **limma** with legacy-style **multiple imputations + stability filtering**
-   Unified QC utilities (PCA, heatmaps, sample distance)
-   A central YAML configuration file controlling all parameters
-   A `{targets}` pipeline for reproducible, dependency-aware execution

------------------------------------------------------------------------

## Getting started (new users)

If you are new to **multiomics-core**, start here:

-   📘 **Onboarding guide:** `docs/onboarding.md`

The onboarding guide explains:

-   How to open the project in RStudio.
-   How to restore the R environment with `renv`.
-   How to run analyses interactively or via `{targets}`.
-   How to reproduce previous runs.

------------------------------------------------------------------------

## Repository structure

```         
R/
├── core/ # Core utilities (contracts, validation, matrix/meta helpers)
├── de/ # Differential expression logic (limma, multi-imputation)
├── preprocess_* # Omics-specific preprocessing pipelines
├── pipeline_* # High-level pipeline entry points
├── plots/ # Pure plotting functions (no I/O)
├── qc/ # QC wrappers (logic + saving)
├── utils.R
config/
├── config.yaml # Central configuration file
docs/ # Design notes and architecture decisions
outputs/ # Analysis outputs (git-ignored)
_targets.R # {targets} pipeline definition
```

------------------------------------------------------------------------

# Requirements

-   R (≥ 4.3 recommended)
-   `renv`
-   System libraries required by Bioconductor (for RNA / proteomics workflows)

------------------------------------------------------------------------

# Setup

Clone the repository and restore the R environment:

``` r
install.packages("renv") 
renv::restore() 
```

## Configuration

#### Edit the central configuration file:

`config/config.yaml`

***The configuration controls***:

-   input and output file paths

-   omics-specific parameters

-   filtering, normalization, and imputation settings

-   differential expression thresholds

-   QC aesthetics (color, shape, sample ID columns)

------------------------------------------------------------------------

## Running the proteomics pipeline (via `{targets}`)

The recommended way to run analyses is via `{targets}`.

From an R session in the project root:

``` r
library(targets)
tar_make()
```

**This will run, in order**:

1.  configuration validation

2.  proteomics input loading

3.  proteomics preprocessing

4.  differential expression using limma with multiple imputations

5.  writing result tables to outputs/

Targets ensures that only steps affected by changes are recomputed.

## Running preprocessing interactively (example)

For exploratory work or debugging, preprocessing can also be run interactively:

``` r
# Load all functions

r_files <- list.files("R", pattern = "\\.[Rr]\$", full.names = TRUE, recursive = TRUE)
invisible(lapply(r_files, source))

# Load config

config <- load_config("config/config.yaml")

# Load inputs

inputs <- load_proteomics_inputs(config)

# Run preprocessing

res <- preprocess_proteomics(inputs, config)

# Example QC: PCA

qc_pca_scatter( expr_mat = res$expr_imp,
  meta     = res$meta, cfg = config$modes$proteomics, out_file = "outputs/proteomics/qc/pca_pc1_pc2.png" )
```

## Starting a new analysis run

To start a new analysis with different parameters:

1.  Copy a config template from `config/templates/`

    ``` bash
    cp config/templates/proteins_config.yaml config/<PROJECT>_<ROUND>.yaml
    ```

2.  Update project identifiers, input paths, and analysis parameters

3.  Place your data under data/<mode>/

4.  Restore the environment if needed:

``` r
renv::restore()
```

5.  Run the pipeline:

``` r
library(targets)
tar_make()
```

Each run produces a separate output folder, allowing multiple analyses to coexist without overwriting results.

## Proteomics DE details

Proteomics differential expression is implemented using:

-   log2-scale expression matrices

-   Perseus-style missing value imputation

-   multiple imputation runs (N repetitions)

-   limma run on each imputed dataset

-   stability filtering (features must pass DE cutoffs in ≥ K/N runs)

All thresholds are controlled via `config.yaml`.

## Reproducibility

-   All package versions are locked in `renv.lock`

-   Outputs and caches are excluded from git

-   `{targets}` provides deterministic, restartable pipelines

## Outputs

All analysis outputs are written to the `outputs/` directory.

-   This directory is intentionally excluded from git
-   Results should be shared by zipping the relevant output folder
-   Each run is isolated by its configuration parameters

This design keeps the repository clean while allowing full reproducibility.

# Status

**Current version**: v0.2.0

Implemented:

-   Proteomics preprocessing

-   Proteomics DE (limma + multi-imputation)

-   Unified config validation

-   `{targets}` pipeline integration

# Planned:

-   RNA-seq DE integration

-   Metabolomics / lipidomics DE

-   Extended QC and reporting
