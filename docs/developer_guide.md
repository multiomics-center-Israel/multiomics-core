## Developer Guide — multiomics-core

This guide is intended for developers who want to **extend**, **modify**, or **maintain** the **multiomics-core** framework.

The project is **configuration-driven**, **modular**, and **reproducible**, with orchestration handled by `{targets}`. The goal of this guide is to protect the infrastructure from drifting into legacy code as the team grows.

------------------------------------------------------------------------

## Core Principles

1.  **Config = analysis decisions, Code = infrastructure** Any decision that may vary between runs must live in YAML configuration, not hard-coded in R.

2.  **Clear separation of concerns**

    -   **I/O**: loading, saving, paths, legacy-compatible outputs
    -   **Logic**: preprocessing, QC computation, DE, clustering
    -   **Plotting**: pure plotting functions (no file writing)
    -   **Orchestration**: `_targets.R` only

3.  **Pure functions by default** Logic functions should not read/write files, rely on global state, or mutate external objects.

4.  **Fail-fast validation** Validate critical objects immediately to prevent silent downstream errors.

5.  **Backward compatibility matters**

    -   API changes must be explicit and intentional
    -   Existing outputs must remain usable
    -   Prefer additive changes over breaking ones

------------------------------------------------------------------------

## High-level Architecture Overview

### Main directories

-   `R/core/` — contracts, validations, matrix/meta helpers, execution snapshots
-   `R/io/` — config loading + input loaders
-   `R/preprocess/` — omics preprocessing logic
-   `R/de/` — differential expression logic (e.g. limma, multi-imputation)
-   `R/qc/` — QC wrappers (may include file writing)
-   `R/plots/` — pure plotting functions (no I/O)
-   `_targets.R` — orchestration only
-   `docs/` — architecture and contribution guides

### Typical pipeline flow (proteomics example)

```         
load_inputs → preprocess → make_imputations → run_DE → summarize → write_outputs + QC
```

------------------------------------------------------------------------

## Contracts: What Each Stage Must Accept and Return

To keep modules composable and stable, each stage must return a predictable structure.

### Example: `preprocess_proteomics()`

Recommended return structure:

-   `expr_*_mat` — numeric matrix (features × samples), colnames = unified `SampleID`
-   `meta` — data.frame with one row per sample, aligned to matrix columns
-   `row_data` / `feature_tbl` — annotations for features, aligned to matrix rows
-   `qc` (optional) — computed QC metrics (no files)

**Golden rule:** If a downstream step consumes `expr_mat + meta`, it should not care how they were loaded or named originally.

------------------------------------------------------------------------

## Validation: When and What to Validate

### When to add validation

-   After `load_*` (files exist, required columns, uniqueness)
-   After sample ID mapping (`sample_map`)
-   After matrix construction (numeric type, alignment, duplicates)
-   After filtering (dimensions + alignment preserved)
-   After imputation generation (number of runs, dims, same sample order)
-   After DE results (required columns, feature alignment, duplicates)

### Minimum checks

-   `is.matrix(expr)` and numeric storage mode
-   `nrow(expr) > 0`, `ncol(expr) > 1`
-   `identical(colnames(expr), meta[[SampleID_col]])`
-   no duplicated critical IDs
-   no silent recycling or partial matching

If something is wrong — fail early with a clear error: **what was expected**, **what was found**, and **how to fix it**.

------------------------------------------------------------------------

## When to Add a Function vs. When to Add a Target

### Add a **function** when:

-   it represents reusable logic
-   it should be unit-testable / contract-validated
-   it is (mostly) pure
-   it will be used both interactively and in `{targets}`

### Add a **target** when:

-   you want caching and dependency tracking
-   it is expensive computation
-   multiple downstream steps depend on it
-   it produces files / artifacts

**Rule of thumb:** Function = *what is computed* Target = *when it runs + what it depends on*

------------------------------------------------------------------------

## How to Add a New Pipeline Step (e.g. clustering module)

### 1) Decide where it belongs

Choose the category:

-   `preprocess/` — transformations before modeling
-   `qc/` — QC metrics and visualization wrappers
-   `de/` — statistical modeling / differential analysis
-   `clustering/` (new) — clustering, embeddings, feature selection

Create a new directory only for a **family of functionality**, not one-off helpers.

------------------------------------------------------------------------

### 2) Define a clear API

Prefer consistent signatures:

-   pass `cfg` instead of many parameters
-   use `config$modes$<mode>` (mode-scoped config)

Example:

``` r
run_clustering_proteomics <- function(expr_mat, meta, cfg, verbose = FALSE) {
  # validations
  # compute clusters + embeddings
  # return results object (no I/O)
}
```

------------------------------------------------------------------------

### 3) Define the output contract

Return an object that downstream steps can rely on:

-   `assignments` — cluster labels per sample / feature
-   `embeddings` — PCA/UMAP coordinates
-   `metrics` — silhouette / withinSS / etc.
-   `params` — parameters snapshot (derived from cfg)

**No file writing here.** Writers handle disk.

------------------------------------------------------------------------

### 4) Integrate with `{targets}`

Add a new target after its dependencies:

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

-   wrap in `write_*`
-   use `format = "file"`

------------------------------------------------------------------------

## Writing to Disk: Writers and Wrappers

### Hard rule

-   `R/plots/` → plotting only, no saving
-   writers (`write_*`) or `R/qc/` wrappers → responsible for disk output

Recommended writer template:

``` r
write_proteomics_cluster_outputs <- function(clust_res, run_dir, cfg, ...) {
  # create directories
  # write tables
  # save plots
  # return character vector of file paths
}
```

Target:

``` r
tar_target(..., format = "file")
```

------------------------------------------------------------------------

## Naming Conventions

### Targets (mode-prefixed)

Use:

-   `prot_*`, `rna_*`, `met_*`, `lip_*`

Common patterns:

-   `*_inputs`
-   `*_pre`
-   `*_imputations`
-   `*_de_res`
-   `*_qc_files`
-   `*_de_files`
-   `*_clusters` / `*_cluster_files`

### Functions

-   Loaders: `load_<mode>_inputs()`
-   Preprocessing: `preprocess_<mode>()`
-   Imputation: `make_imputations_<mode>()`
-   Modeling: `run_<model>_<mode>()`
-   Summaries: `summarize_<model>_<mode>()`
-   Writers: `write_<mode>_<artifact>()`
-   Validation: `validate_<object>()`

### Internal objects

-   matrices: `expr_*_mat`
-   metadata: always `meta`
-   config: `cfg` (`config$modes$<mode>`)

------------------------------------------------------------------------

## How Not to Break Backward Compatibility

1.  **Do not rename existing return fields**

    -   add new fields instead
    -   deprecate old ones gradually with warnings

2.  **Preserve legacy output formats**

    -   if external tools depend on them, changes must be opt-in via config

3.  **Config schema**

    -   new fields must have defaults
    -   semantic changes require docs + version bump

4.  **Versioning**

    -   bump version on meaningful API or behavior change
    -   document changes under `docs/`

------------------------------------------------------------------------

## Working Correctly with `{targets}`

### Core rules

-   targets must be deterministic given inputs + config
-   expensive steps should run once and be reused
-   no hidden I/O outside tracked file targets

### Seeds & reproducibility

-   use `params.seed` as the global anchor
-   derive all random seeds systematically from it
-   avoid ad-hoc `set.seed()` not tied to config

### Debugging

-   the same functions must work interactively and in `{targets}`
-   do not maintain parallel “debug-only” code paths

------------------------------------------------------------------------

## Pre-PR / Pre-Merge Checklist

-   [ ] validations added where needed (fail-fast)
-   [ ] single responsibility per function
-   [ ] no hidden I/O in logic
-   [ ] new target added for expensive/reusable steps
-   [ ] naming conventions followed
-   [ ] no breaking changes to existing contracts
-   [ ] new config fields have defaults + docs
-   [ ] pipeline runs end-to-end (`tar_make()`), caching behaves as expected

------------------------------------------------------------------------

## Recommended Design Patterns

### “Thin orchestration, thick modules”

`_targets.R` stays short; all logic lives in `R/`.

### “Data in, data out”

Modules consume objects and return objects; only writers touch disk.

### “Mode-first thinking”

When adding functionality, consider:

-   what is shared across omics?
-   what is mode-specific?
-   can this scale to RNA-seq / metabolomics / lipidomics?
