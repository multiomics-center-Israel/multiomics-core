# multiomics-core

A modular and reproducible R framework for single-omics and multi-omics analyses (RNA-seq, proteomics, metabolomics, lipidomics).

### The project emphasizes

-   clear separation of concerns (I/O, preprocessing, DE, QC, plotting)
-   configuration-driven workflows with reproducibility via `renv`
-   scalable orchestration via `{targets}`

#### This project relies heavily on {targets} for reproducible pipeline orchestration.

For in-depth documentation and tutorials, see the official targets book: <https://books.ropensci.org/targets/>

------------------------------------------------------------------------

## What this repository provides

-   Standardized data loading and validation
-   Omics-specific preprocessing (filtering, normalization, imputation)
-   Proteomics differential expression via a method-based interface (currently **limma**) with legacy-style multiple imputations and stability filtering
-   Unified QC utilities (PCA, heatmaps, sample distance)
-   A central YAML configuration file controlling all parameters
-   A `{targets}` pipeline for reproducible, dependency-aware execution

------------------------------------------------------------------------

## Getting started (new users)

If you are new to **multiomics-core**, start here:

-   📘 **Onboarding guide:** `docs/onboarding.md`
-   📘 **Developer guide:** `docs/developer_guide.md`

The onboarding guide explains:

-   How to open the project in RStudio
-   How to restore the R environment with `renv`
-   How to run analyses interactively or via `{targets}`
-   How to reproduce previous runs

------------------------------------------------------------------------

## Repository structure

```         
R/
├── core/        # Core utilities (contracts, validation, matrix/meta helpers)
├── io/          # Data loading, config parsing, and I/O helpers
├── preprocess/  # Omics-specific preprocessing (filtering, normalization, imputation)
├── de/          # Differential expression logic (engines, summarization, builders)
├── qc/          # QC computations (PCA, density, correlation)
├── plots/       # Pure plotting functions (no I/O)
├── clustering/  # Clustering algorithms and legacy exporters
├── pipeline/    # Pipeline modules and {targets} factories
config/
├── config.yaml              # Central configuration file
├── templates/               # Analysis config templates
docs/                         # Onboarding and developer documentation
outputs/                      # Analysis outputs (git-ignored)
_targets.R                    # {targets} pipeline definition
```

------------------------------------------------------------------------

## Requirements

-   **R ≥ 4.3** (tested with R 4.5.x)
-   **RStudio** (recommended)
-   **`renv`** (for reproducible environments)

### Windows users (IMPORTANT)

On Windows, some CRAN / Bioconductor packages may need to be **compiled from source** (e.g. `SparseArray`, `IRanges`, `Biobase`). Therefore, **Rtools is required**.

**Install Rtools (matching your R version):** 👉 <https://cran.r-project.org/bin/windows/Rtools/>

After installation, **restart R / RStudio**, then verify:

``` r
Sys.which("gcc")
Sys.which("make")
```

Both commands should return a valid path. If they return `""`, Rtools is not correctly installed or not on `PATH`.

Missing Rtools may cause `renv::restore()` to fail.

### Linux / macOS users

-   A standard compiler toolchain is required
-   System libraries commonly needed by Bioconductor (e.g. `libxml2`, `curl`, `openssl`)

------------------------------------------------------------------------

## Package repositories (CRAN + Bioconductor)

This project relies on **CRAN** and **Bioconductor** packages.

We recommend using **Posit Package Manager (PPM)** for CRAN together with standard Bioconductor repositories.

Recommended setup:

``` r
options(repos = c(
  CRAN = "https://packagemanager.posit.co/cran/latest"
))

BiocManager::repositories()
```

The Bioconductor version must match the one recorded in `renv.lock` (e.g. `Bioconductor 3.22`).

To check:

``` r
BiocManager::version()
```

If needed:

``` r
BiocManager::install(version = "3.22")
```

------------------------------------------------------------------------

## Setup

Clone the repository and restore the R environment:

``` r
install.packages("renv")
renv::restore()
```

------------------------------------------------------------------------

## Configuration

Edit the central configuration file:

```         
config/config.yaml
```

Or start from a template:

``` bash
cp config/templates/proteins_config.yaml config/<PROJECT>_<ROUND>.yaml
```

The configuration controls:

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

This will run, in order:

1.  configuration validation
2.  proteomics input loading
3.  proteomics preprocessing
4.  differential expression (multi-imputation, method-based)
5.  QC analysis and diagnostic plots (if enabled)
6.  clustering analysis (if enabled)
7.  writing result tables to `outputs/`

`{targets}` ensures that only steps affected by changes are recomputed.

### Learning more about `{targets}`

This project uses `{targets}` for reproducible, dependency-aware pipeline orchestration.

For a detailed introduction, tutorials, and best practices, see the official **targets** book: <https://books.ropensci.org/targets/>

------------------------------------------------------------------------

## Running preprocessing interactively (example)

For exploratory work or debugging:

``` r
# Load all functions
r_files <- list.files("R", pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE)
invisible(lapply(r_files, source))

# Load config
config <- load_config("config/config.yaml")

# Load inputs
inputs <- load_proteomics_inputs(config)

# Run preprocessing
res <- preprocess_proteomics(inputs, config)

# Example QC: PCA
qc_pca_scatter(
  expr_mat = res$expr_imp_single,
  meta     = res$meta,
  cfg      = config$modes$proteomics,
  out_file = "outputs/proteomics/qc/pca_pc1_pc2.png"
)
```

------------------------------------------------------------------------

## Reproducibility

-   All package versions are locked in `renv.lock`
-   Outputs and caches are excluded from git
-   `{targets}` provides deterministic, restartable pipelines
-   Each run records execution metadata (config snapshot, git commit, session info)

------------------------------------------------------------------------

## Outputs

All analysis outputs are written to the `outputs/` directory.

-   This directory is intentionally excluded from git
-   Results should be shared by zipping the relevant output folder
-   Each run is isolated by its configuration parameters

------------------------------------------------------------------------

## Developer notes

If you want to extend, modify, or maintain **multiomics-core**, see:

-   📘 **Developer guide:** `docs/developer_guide.md`

------------------------------------------------------------------------

## Status

**Current version:** v0.2.0

### Implemented

-   Proteomics preprocessing
-   Proteomics DE (multi-imputation, method-based; currently limma)
-   QC module
-   Clustering module (hierarchical, partition, binary patterns)
-   Unified config validation
-   `{targets}` pipeline integration

### Planned

-   MA / Volcano plots
-   RNA-seq DE integration
-   Metabolomics / lipidomics DE
-   Multi-omics integration and reporting
