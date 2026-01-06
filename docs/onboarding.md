# Onboarding — multiomics-core

This document is intended for a new team member who wants to **run an existing analysis or start a new project** using **multiomics-core**.

The framework is **configuration-driven**, fully reproducible, and designed to be used primarily through **RStudio**.

------------------------------------------------------------------------

## 1. Working with multiomics-core in RStudio

This project is fully compatible with **RStudio**, and new users are encouraged to work through RStudio rather than the command line.

------------------------------------------------------------------------

### 1.1 Clone and open the project in RStudio

**Option A — via terminal (recommended if you use git regularly):**

``` bash
git clone <REPO_URL>
```

Then open RStudio and choose:

-   *File → Open Project…*
-   Select `multiomics-core` (the folder containing the `.Rproj` file)

**Option B — via RStudio UI:**

-   *File → New Project → Version Control → Git*
-   Paste the repository URL
-   Choose a local directory

> Always work inside the `.Rproj`. This ensures correct relative paths and a clean R session.

------------------------------------------------------------------------

### 1.2 Restore the R environment (renv)

Once the project is open in RStudio, run in the **Console**:

``` r
install.packages("renv")   # once per machine
renv::restore()
```

Verify that the environment is consistent:

``` r
renv::status()
sessionInfo()
```

------------------------------------------------------------------------

## 2. Expected data structure

The pipeline expects input data under:

```         
data/<mode>/...
```

Example — proteomics:

```         
data/
  proteomics/
    protein_matrix.csv
    sample_map.csv
    meta_prot.csv
    contrasts_prot.csv
```

> File names and paths are fully controlled via the configuration file.\
> The repository does **not** enforce file names — only consistency with the config.

------------------------------------------------------------------------

## 3. Starting a new project (configuration-driven)

### 3.1 Create a new config file from a template

**Do not edit `config/config.yaml` directly.**

Always start from a template:

``` bash
cp config/templates/proteins_config.yaml config/<PROJECT>_<ROUND>.yaml
```

Example:

``` bash
cp config/templates/proteomics.yaml config/E_Pick_A02.yaml
```

------------------------------------------------------------------------

### 3.2 Fields you **must** update

-   `project.name`
-   `project.analysis_round`
-   `paths.out` (if different from default)
-   `modes.proteomics.files.*` (paths to your data files)
-   `modes.proteomics.id_columns.*`
-   `params.seed`

------------------------------------------------------------------------

### 3.3 Fields you may / should update

-   filtering parameters (`filtering.min_count`)
-   normalization method
-   imputation parameters (method, repetitions, thresholds)
-   de thresholds
-   QC aesthetics (`effects`)

> Starting a new project should **never require modifying R code**.

------------------------------------------------------------------------

## 4. Running the pipeline with `{targets}`

The recommended way to run analyses is via `{targets}`.

From an R session in the project root:

``` r
library(targets)
tar_make()
```

Useful helpers:

``` r
tar_visnetwork()    # visualize dependency graph
tar_progress()      # execution status
tar_read(prot_de_res)  # read DE summary results
tar_destroy()       # clear cache (use with care)
```

------------------------------------------------------------------------

## 5. Interactive execution (for debugging & exploration)

For exploratory work or debugging, steps can be run interactively.

``` r
# Load all functions
r_files <- list.files("R", pattern = "\\.[Rr]\$", full.names = TRUE, recursive = TRUE)
invisible(lapply(r_files, source))
```

``` r
config <- load_config("config/proj1_02.yaml")
validate_config(config)
```

``` r
prot_inputs <- load_proteomics_inputs(config)
prot_pre    <- preprocess_proteomics(prot_inputs, config)
```

### 5.1 Multiple imputation

``` r
prot_imps <- make_imputations_proteomics(
  expr_mat  = prot_pre$expr_filt,
  cfg       = config,
  verbose   = TRUE
)

validate_proteomics_imputations(
  imputations = prot_imps,
  meta        = prot_pre$meta,
  cfg         = config
)
```

### 5.2 Differential expression (limma)

``` r
limma_runs <- run_limma_on_imputations_proteomics(
  imputations  = prot_imps,
  meta         = prot_pre$meta,
  contrasts_df = prot_inputs$contrasts,
  prot_tbl     = prot_inputs$protein,
  cfg          = config
)

runs_de_tables <- lapply(limma_runs, function(x) x$de_tables)

summary_df <- summarize_limma_mult_imputation(
  runs_de_tables = runs_de_tables,
  config         = config
)
```

------------------------------------------------------------------------

## 6. Reproducing a previous run

To reproduce an analysis exactly:

1.  Use the same git commit (hash) of the repository
2.  Run `renv::restore()` with the same `renv.lock`
3.  Use the exact same config file and `params.seed`

Each pipeline run automatically generates an **`execution_info/`** directory containing:

-   config snapshot
-   git commit hash
-   execution timestamp

These files should be archived together with the outputs.

------------------------------------------------------------------------

## Working principles

-   **Config = analysis decisions**
-   **Code = infrastructure**
-   No copy–paste of legacy scripts
-   Every run must be reproducible

------------------------------------------------------------------------

For questions, contact the **multiomics-core** maintainer.
