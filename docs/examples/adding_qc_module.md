---
editor_options: 
  markdown: 
    wrap: 72
---

# Example — Adding a New QC Module (Worked Example)

**Goal:** Add a **sample-level clustering QC** step for proteomics
without breaking the pipeline.

We will:

-   add a pure QC computation function
-   add a plotting function
-   add a writer wrapper
-   integrate with `{targets}`

------------------------------------------------------------------------

## Step 1 — Define the QC goal clearly

We want:

-   hierarchical clustering of samples
-   based on preprocessed expression matrix
-   configurable distance + linkage
-   reproducible
-   saved as a PNG

**Key decision:** This is *QC*, not preprocessing → lives under `R/qc/`
and `R/plots/`.

------------------------------------------------------------------------

## Step 2 — Add config fields (with defaults)

In `config.yaml` (or template):

``` yaml
modes:
  proteomics:
    qc:
      clustering:
        distance: "euclidean"
        linkage: "complete"
        scale: true
```

Rules applied:

-   Nested under mode
-   Descriptive names
-   Sensible defaults

------------------------------------------------------------------------

## Step 3 — Pure QC computation function

**File:** `R/qc/qc_clustering_proteomics.R`

``` r
qc_cluster_samples_proteomics <- function(expr_mat, meta, cfg) {
  validate_expr_meta_alignment(expr_mat, meta)

  mat <- if (isTRUE(cfg$qc$clustering$scale)) {
    scale(t(expr_mat))
  } else {
    t(expr_mat)
  }

  dist_mat <- dist(mat, method = cfg$qc$clustering$distance)
  hc <- hclust(dist_mat, method = cfg$qc$clustering$linkage)

  list(
    hclust = hc,
    params = cfg$qc$clustering
  )
}
```

✔ pure ✔ validated ✔ config-driven ✔ returns objects, not files

------------------------------------------------------------------------

## Step 4 — Plotting function (pure)

**File:** `R/plots/plot_sample_dendrogram.R`

``` r
plot_sample_dendrogram <- function(hclust_obj, meta) {
  plot(
    hclust_obj,
    main = "Sample clustering",
    xlab = "",
    sub  = ""
  )
}
```

✔ no paths ✔ no saving ✔ reusable

------------------------------------------------------------------------

## Step 5 — Writer wrapper

**File:** `R/qc/write_proteomics_qc_clustering.R`

``` r
write_proteomics_qc_clustering <- function(qc_res, meta, run_dir) {
  out_dir <- file.path(run_dir, "qc")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_file <- file.path(out_dir, "sample_clustering.png")

  png(out_file, width = 1200, height = 1000)
  plot_sample_dendrogram(qc_res$hclust, meta)
  dev.off()

  out_file
}
```

✔ isolated I/O ✔ returns file path ✔ deterministic

------------------------------------------------------------------------

## Step 6 — Integrate into `{targets}`

In `_targets.R`:

``` r
tar_target(
  prot_qc_clustering,
  qc_cluster_samples_proteomics(
    expr_mat = prot_pre$expr_filt_mat,
    meta     = prot_pre$meta,
    cfg      = config$modes$proteomics
  )
),

tar_target(
  prot_qc_clustering_file,
  write_proteomics_qc_clustering(
    qc_res  = prot_qc_clustering,
    meta    = prot_pre$meta,
    run_dir = prot_run_dir
  ),
  format = "file"
)
```

✔ computation cached ✔ output tracked ✔ clean dependency graph

------------------------------------------------------------------------

## Step 7 — Sanity check

-   `tar_visnetwork()` shows QC depends on preprocessing

-   Changing clustering params invalidates only QC targets

-   Interactive usage works:

    ``` r
    qc <- qc_cluster_samples_proteomics(expr, meta, cfg)
    plot_sample_dendrogram(qc$hclust, meta)
    ```

------------------------------------------------------------------------

## Why This Pattern Works

-   QC logic is reusable
-   Plotting is testable
-   I/O is explicit
-   Config controls behavior
-   No legacy entanglement

This is the **canonical pattern** for adding:

-   PCA QC
-   batch effect diagnostics
-   clustering
-   embeddings
-   feature-level QC summaries
