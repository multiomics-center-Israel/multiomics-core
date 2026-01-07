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

**Do not edit `config/config.yaml` manually unless you know what you are doing.**

Instead, start from a template:

``` bash
cp config/templates/proteomics.yaml config/<PROJECT>_<ROUND>.yaml
```

Then point the pipeline to your new config by updating the config_file target in `_targets.R`:

``` r
tar_target(config_file, "config/<PROJECT>_<ROUND>.yaml", format = "file")
```

Example:

``` r
tar_target(config_file, "config/E_Pick_A02.yaml", format = "file")
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

### Learning more about `{targets}`

This project uses `{targets}` for reproducible, dependency-aware pipeline orchestration.

For a detailed introduction, tutorials, and best practices, see the official **targets** book: <https://books.ropensci.org/targets/>

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

------------------------------------------------------------------------

### 5.2 Differential expression (method-based)

Differential expression is executed via a **method-agnostic builder**, controlled entirely by the configuration file (`cfg$de$method`).

At this stage, the proteomics pipeline supports:

-   `method: "limma"` — limma with multiple imputations + stability filtering

Additional DE methods may be added in the future without changing the pipeline interface.

### Run DE interactively

``` r
prot_de_res <- build_proteomics_de_results(
  pre          = prot_pre,
  inputs       = prot_inputs,
  config       = config,
  verbose      = TRUE
)
```

This function performs, internally:

1.  generation of multiple imputed datasets (if configured)
2.  execution of the selected DE method
3.  summarization across imputations
4.  validation of the resulting DE tables

### Returned object

`prot_de_res` is a list with the following fields:

-   `runs` — per-imputation DE results (method-specific)
-   `runs_de_tables` — extracted DE tables per contrast
-   `summary_df` — unified DE summary table used downstream
-   `method` — DE method used (e.g. `"limma"`)

Downstream steps (clustering, QC, writing outputs) **must not assume a specific DE method**, and should rely only on the standardized fields above.

------------------------------------------------------------------------

### Configuration example

``` yaml
modes:
  proteomics:
    de:
      method: "limma"
      use_adj_for_pass1: true
      p_cutoff: 0.1
      linear_fc_cutoff: 1.5
```

Changing the DE method requires **only updating the config**, not modifying R code.

------------------------------------------------------------------------

### Design note (for developers)

-   `build_proteomics_de_results()` is the **single entry point** for DE
-   Method-specific logic lives under `R/de/`
-   No downstream code should call `run_limma_*()` directly

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
