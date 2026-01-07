# Developer Guide — multiomics-core

This guide is intended for developers who want to **extend**, **modify**, or **maintain** the **multiomics-core** framework.

The project is **configuration-driven**, **modular**, and **fully reproducible**, with orchestration handled by `{targets}`. The goal of this guide is to protect the infrastructure from degrading into legacy code as the team grows.

------------------------------------------------------------------------

## Core Principles

1.  **Config = analysis decisions, Code = infrastructure** Any decision that may vary between runs must live in YAML configuration, not hard-coded in R.

2.  **Clear separation of concerns**

    -   **I/O**: loading, saving, paths, legacy-compatible outputs
    -   **Logic**: preprocessing, QC, DE, clustering
    -   **Plotting**: pure plotting functions (no file writing)
    -   **Orchestration**: `_targets.R` only

3.  **Pure functions by default** Logic functions should not read/write files, rely on global state, or mutate external objects.

4.  **Fail-fast validation** Critical objects must be validated immediately to prevent silent downstream errors.

5.  **Backward compatibility matters**

    -   API changes must be explicit and intentional
    -   Existing outputs must remain usable
    -   Prefer additive changes over breaking ones

------------------------------------------------------------------------

## High-level Architecture Overview

### Main directories

-   `R/core/` — contracts, validations, matrix/meta helpers, utilities
-   `R/io/`— config parsing, loaders, I/O helpers, Excel legacy writers
-   `R/preprocess/` — omics preprocessing (filtering, normalization, imputation)
-   `R/de/` — DE engines + summarization + builders (method-based)
-   `R/qc/` — QC computations (PCA, density, correlation)
-   `R/plots/` — pure plotting (no I/O)
-   `R/clustering/` — clustering algorithms + legacy exporters
-   `R/pipeline/` — pipeline modules + {targets} factories
-   `_targets.R` — orchestration only

### Typical pipeline flow (proteomics example)

```         
load_inputs → preprocess → DE (impute+run+summarize) → QC/clustering → write_outputs
```

------------------------------------------------------------------------

## Contracts: What Each Stage Must Accept and Return

To keep modules composable and stable, each **main stage returns a predictable structure**.

#### Example: `preprocess_proteomics()`

Recommended return structure:

-   `expr_raw`, `expr_filt`, `expr_imp_single` — numeric matrices (features × samples)
-   `meta` — data.frame (samples × covariates)
-   `row_data` — feature annotations aligned to rows of expr matrices
-   `qc` / `*_qc` (optional)

**Golden rule:** If a downstream step consumes `expr_mat + meta`, it should not care *how* they were loaded or named originally.

------------------------------------------------------------------------

## Validation: When and What to Validate

### When to add validation

-   After `load_*` (files, columns, uniqueness)
-   After sample ID mapping (`sample_map`)
-   After matrix construction (alignment, numeric type, NA policy)
-   After imputation generation (number of runs, dimensions, alignment)
-   After DE results (required columns, duplicated features)

### Minimum checks

-   `is.matrix(expr)` and numeric storage mode
-   `nrow(expr) > 0`, `ncol(expr) > 1`
-   Column/sample alignment with metadata
-   No duplicated critical IDs
-   No silent dimension recycling

> If something is wrong — fail early with a clear error: *what was expected*, *what was found*, and *how to fix it*.

------------------------------------------------------------------------

## When to Add a Function vs. When to Add a Target

### Add a **function** when:

-   It represents reusable logic
-   It should be unit-testable or contract-validated
-   It is (mostly) pure
-   It will be used both interactively and in `{targets}`

### Add a **target** when:

-   You want caching and dependency tracking
-   It is an expensive computation
-   Multiple downstream steps depend on it
-   It produces files/artifacts

**Rule of thumb:** Function = *what is computed* Target = *when and with which dependencies*

------------------------------------------------------------------------

## How to Add a New Pipeline Step

### 1. Design decision: where does it belong?

Choose the category:

-   `preprocess` — transformations before modeling
-   `qc` — quality metrics or visualizations
-   `de` — statistical models / differential analysis
-   `clustering` — clustering algorithms and summaries (writing handled by modules/writers)

Create a new directory only for **families of functionality**, not single functions.

------------------------------------------------------------------------

### 2. Define a clear API

Prefer consistent signatures:

-   Pass `cfg` instead of many parameters
-   Use `config$modes$<mode>` when possible

Example:

``` r
run_clustering_proteomics <- function(expr_mat, meta, cfg, verbose = FALSE) {
  # validations
  # algorithm
  # return list(assignments = ..., embeddings = ..., metrics = ...)
}
```

------------------------------------------------------------------------

### 3. Define the output contract

Clearly document what the function returns:

-   `assignments` — cluster labels
-   `embeddings` — PCA/UMAP coordinates
-   `params` — parameter snapshot
-   **No file writing**

File writing belongs to wrappers/writers.

------------------------------------------------------------------------

### 4. Integrate with `{targets}`

Add a new target **after** its dependencies:

``` r
tar_target(
  prot_clusters,
  run_clustering_proteomics(
    expr_mat = prot_pre$expr_filt,
    meta     = prot_pre$meta,
    cfg      = config$modes$proteomics,
    verbose  = TRUE
  )
)
```

For file outputs:

-   Wrap with `write_*`
-   Use `format = "file"`

------------------------------------------------------------------------

## Writing to Disk: Writers and Wrappers

### Hard rule

-   `R/plots/` → plotting only, no saving
-   `R/qc/` and `write_*` → responsible for disk output

### Recommended writer template

``` r
write_proteomics_cluster_outputs <- function(clust_res, run_dir, cfg, ...) {
  # create directories
  # write tables
  # save plots
  # return character vector of file paths
}
```

Target definition:

``` r
tar_target(..., format = "file")
```

------------------------------------------------------------------------

## Naming Conventions

### Targets

-   Prefix by mode: `prot_*`, `rna_*`, `met_*`, `lip_*`

-   Common patterns:

    -   `*_inputs`
    -   `*_pre`
    -   `*_imputations`
    -   `*_de_res`
    -   `*_qc_files`
    -   `*_de_files`

### Functions

-   Loaders: `load_<mode>_inputs()`
-   Preprocessing: `preprocess_<mode>()`
-   Imputation: `make_imputations_<mode>()`
-   Modeling: `run_<model>_<mode>()`
-   Summaries: `summarize_<model>_<mode>()`
-   Writers: `write_<mode>_<artifact>()`
-   Validation: `validate_<object>()`

### Internal objects

-   Matrices: `expr_*_mat`
-   Metadata: always `meta`
-   Config: `cfg` (usually `config$modes$<mode>`)

------------------------------------------------------------------------

## How Not to Break Backward Compatibility

1.  **Avoid renaming return fields once they are part of the public contract. If a rename is necessary, do it in a single refactor PR + update docs + bump version.**

    -   Add new fields instead
    -   Deprecate old ones gradually with warnings

2.  **Preserve legacy output formats**

    -   If external tools depend on them, changes must be opt-in via config

3.  **Config schema**

    -   New fields must have defaults
    -   Semantic changes require documentation + version bump

4.  **Versioning**

    -   Bump version on meaningful API or behavior changes
    -   Document changes in `docs/`

------------------------------------------------------------------------

## Working Correctly with `{targets}`

### Principles

-   Targets must be deterministic given inputs + config
-   Expensive computations should run once and be reused
-   No hidden file I/O outside tracked targets

### Seeds and reproducibility

-   Use `params.seed` as the global seed anchor
-   Derive imputation seeds systematically
-   Avoid ad-hoc `set.seed()` calls unrelated to config
-   All stochastic steps (e.g., imputations) must derive seeds from params.seed and record them in outputs/metadata.

### Debugging

-   Use the same functions interactively as in `{targets}`
-   Do not maintain parallel “debug-only” logic

------------------------------------------------------------------------

## Pre-PR / Pre-Merge Checklist

-   [ ] Validations added where needed (fail-fast)
-   [ ] Function follows single-responsibility principle
-   [ ] No hidden I/O in logic
-   [ ] New target added for expensive or reusable steps
-   [ ] Naming conventions followed
-   [ ] No breaking changes to existing contracts
-   [ ] New config fields have defaults + documentation
-   [ ] Full pipeline runs (`tar_make()`) with sensible caching

------------------------------------------------------------------------

## Recommended Design Patterns

### “Thin orchestration, thick modules”

`_targets.R` should remain short and readable. All logic lives in `R/`.

### “Data in, data out”

Modules consume objects and return objects; only writers touch disk.

### “Mode-first thinking”

When adding functionality, consider:

-   What is shared across omics?
-   What is mode-specific?
-   Can this scale to RNA-seq, metabolomics, lipidomics?

------------------------------------------------------------------------
