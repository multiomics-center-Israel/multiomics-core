# Metabolomics Preprocessing Pipeline

Production-grade preprocessing for metabolomics data in the multiomics-core
{targets} pipeline.

---

## Pipeline Overview

The metabolomics pipeline runs in two conceptual stages:

| Stage | Targets prefix | Description |
|---|---|---|
| **0 — Preprocessing** | `met_*` | New preprocessing DAG (this document) |
| **1 — QC** | `metab_qc_pre_obj` | PCA, heatmaps, normalization diagnostics |
| **2 — Analysis** | `metab_de_*`, `metab_enrichment_*` | DE, feature selection, enrichment |

Stages 1 and 2 consume `metab_pre`, an adapter target that converts the new
`met_*` outputs into the established pre-processing contract.

---

## DAG (Stage 0)

```
met_input_files ──► met_raw
                        ├──► met_missingness_stats
                        └──► met_filtered
                                  └──► met_imputed
                                            └──► met_log
                                                      ├──► met_norm_tss    ──► met_norm_tss_qc
                                                      ├──► met_norm_median ──► met_norm_median_qc
                                                      ├──► met_norm_pqn    ──► met_norm_pqn_qc
                                                      └──► met_norm_comparison (TSV)
                                                 [chosen via config]
                                                      └──► met_corrected
                                                                └──► metab_pre (adapter)
                                                                          └──► (Stage 1 & 2)
```

---

## Processing Steps

### 1. Input loading (`met_raw`)

- Supports two input formats:
  - `cd_raw` — Compound Discoverer raw export (`Area: SAMPLEID.raw (Fn)` columns)
  - `processed_wide` — pre-processed wide table with separate metadata
- Optional sample filter: excludes QC/blank samples before preprocessing
  (`sample_filter.enabled: true` in config)
- Metadata is built automatically if not provided

### 2. Missingness classification (`met_missingness_stats`)

Operated on the **pre-filter** matrix to give the complete picture.

- Per-feature Spearman correlation between missingness indicator and sample
  intensity rank
- Negative correlation → **MNAR** (missing in low-intensity samples → min/2)
- Non-negative correlation → **MAR** (missing at random → KNN)
- Special classes: `all_observed`, `all_missing`

Threshold: `preprocessing.mnar_threshold` (default 0.3, sign-agnostic)
Saved to: `diagnostic_plots/missingness_heatmap.png`

### 4. Missingness filtering (`met_filtered`)

Applied sequentially:

1. **Samples first**: drop samples where `% missing > sample_missing_threshold`
2. **Features second**: drop features where `% missing > feat_missing_threshold`
   (evaluated on the retained sample set)

### 5. Imputation (`met_imputed`)

Re-classifies missingness on the **filtered** matrix, then:

- **MNAR features** → `min(observed) / 2` per feature
  (global min/2 fallback for all-NA rows)
- **MAR features** → `impute::impute.knn` (Bioconductor)
  - `k` clamped to `min(knn_k, ncol - 1)`; column-mean fallback if `k < 1`
  - All-NA MAR rows are logged and skipped

### 6. Log transformation (`met_log`)

`log2(x + pseudocount)` applied to the imputed matrix.

Config: `normalization.transform` (default `"log2"`), `normalization.pseudocount`
(default `1`).

### 7. Parallel normalization (`met_norm_tss`, `met_norm_median`, `met_norm_pqn`)

All three methods always run independently so {targets} can cache each:

| Target | Method | Description |
|---|---|---|
| `met_norm_tss` | TSS | Scale column sums to median column sum |
| `met_norm_median` | Median | Scale by sample median |
| `met_norm_pqn` | PQN | Probabilistic Quotient Normalization |

Each target saves a boxplot and PCA scatter to `diagnostic_plots/`.

**Note**: `normalization.scaling` (auto/pareto/range) is **not applied** in
the new pipeline — sample normalization replaces scaling for `expr_work`.

### 8. Normalization comparison (`met_norm_comparison`)

TSV of QC metrics (median RSD, PC1 variance) for all three methods plus the
log-only baseline. Saved to `datasets/normalization_comparison.tsv`.

### 9. Drift correction + final selection (`met_corrected`)

1. Picks the normalization matrix from `preprocessing.chosen_norm`
   (`"tss"`, `"median"`, or `"pqn"`)
2. Applies LOESS drift correction **only if**:
   - `preprocessing.drift_correction.enabled: true`
   - Metadata has columns `injection_order_col` (numeric, no NAs)
     and `qc_flag_col` (logical)
   - ≥ 4 QC samples (`qc_flag_col == TRUE`)
3. If preconditions fail, returns the chosen norm matrix unchanged (no error)
4. Saves before/after PCA coloured by injection order if drift was applied

---

## Configuration Reference

Add this block under `modes: metabolomics:` in your `config.yaml`.

```yaml
modes:
  metabolomics:

    input:
      format: "cd_raw"          # "cd_raw" | "processed_wide"

    files:
      data: "metabolomics/cd_export.xlsx"
      metadata: "metabolomics/sample_metadata.csv"  # optional

    preprocessing:
      # Missingness filtering
      feat_missing_threshold: 0.20     # drop features missing in > 20% samples
      sample_missing_threshold: 0.30   # drop samples missing in > 30% features
      mnar_threshold: 0.3              # |Spearman corr| for MNAR classification

      # Imputation
      knn_k: 10                        # KNN neighbours for MAR features
      epsilon: 1.0e-8                  # near-zero floor for imputed values

      # Normalization selection (TSS/Median/PQN always run; this picks the one
      # forwarded to drift correction and all downstream targets)
      chosen_norm: "pqn"               # "tss" | "median" | "pqn"

      # LOESS drift correction (skipped silently if preconditions absent)
      drift_correction:
        enabled: true
        injection_order_col: "injection_order"
        qc_flag_col: "is_QC"
        loess_span: 0.75
        epsilon: 1.0e-8

    normalization:
      transform: "log2"          # used by met_log target
      pseudocount: 1
      # sample_norm, scaling, na_policy, evaluate_methods are legacy fields
      # (used by preprocess_metabolomics() only, not the new pipeline)
```

---

## Output Structure

```
outputs/<run_id>/metabolomics/
├── diagnostic_plots/
│   ├── missingness_heatmap.png         # binary MNAR/MAR heatmap
│   ├── norm_tss_boxplot.png
│   ├── norm_tss_pca.png
│   ├── norm_median_boxplot.png
│   ├── norm_median_pca.png
│   ├── norm_pqn_boxplot.png
│   ├── norm_pqn_pca.png
│   ├── drift_correction_pca.png        # before/after (if drift applied)
│   ├── PCA_PC1.vs.PC2.png
│   └── ...
└── datasets/
    ├── normalization_comparison.tsv    # RSD + PC1 var for all three methods
    ├── missingness_per_sample.tsv
    ├── missingness_per_feature.tsv
    ├── normalization_applied.tsv
    ├── expr_normalized.tsv
    ├── expr_raw.tsv
    ├── metadata_aligned.tsv
    └── feature_annotations.tsv
```

---

## How to Run

```r
library(targets)

# Set config path (or use default config.yaml in project root)
Sys.setenv(MULTIOMICS_CONFIG = "config/my_config.yaml")

# Inspect the DAG
tar_manifest(fields = c("name", "command"))
tar_visnetwork()

# Run only preprocessing
tar_make(names = starts_with("met_"))

# Run everything
tar_make()
```

### Invalidation behaviour

| Config change | Targets invalidated |
|---|---|
| `chosen_norm` | `met_corrected`, `metab_pre`, Stage 1+2 |
| `feat_missing_threshold` | `met_filtered` and all downstream |
| `mnar_threshold` | `met_missingness_stats`, `met_imputed` and downstream |
| `drift_correction.*` | `met_corrected`, `metab_pre`, Stage 1+2 |
| Input data file | `met_input_files`, all `met_*`, all `metab_*` |
| Normalization method params | Only the affected `met_norm_*` and downstream |

TSS, Median, and PQN targets remain **cached** when only `chosen_norm` changes.

---

## Required Packages

```r
# CRAN
install.packages(c("targets", "yaml", "ggplot2", "pheatmap",
                   "patchwork", "readr", "readxl"))

# Bioconductor
BiocManager::install("impute")  # KNN imputation for MAR features
```
