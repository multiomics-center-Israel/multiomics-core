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
├── config.yml # Central configuration file
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

`config/config.yml`

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

config <- load_config("config/config.yml")

# Load inputs

inputs <- load_proteomics_inputs(config)

# Run preprocessing

res <- preprocess_proteomics(inputs, config)

# Example QC: PCA

qc_pca_scatter( expr_mat = res$expr_imp,
  meta     = res$meta, cfg = config$modes$proteomics, out_file = "outputs/proteomics/qc/pca_pc1_pc2.png" )
```

## Proteomics DE details

Proteomics differential expression is implemented using:

-   log2-scale expression matrices

-   Perseus-style missing value imputation

-   multiple imputation runs (N repetitions)

-   limma run on each imputed dataset

-   stability filtering (features must pass DE cutoffs in ≥ K/N runs)

All thresholds are controlled via `config.yml`.

## Reproducibility

-   All package versions are locked in `renv.lock`

-   Outputs and caches are excluded from git

-   `{targets}` provides deterministic, restartable pipelines

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
