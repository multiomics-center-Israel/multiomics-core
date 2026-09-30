# Shiny payload: contract and handoff (v2.0)

> **Status: for review by the app owner.**
> This document describes what the `multiomics-core` pipeline **emits**. It was
> written without access to the app code, so nothing here has been checked
> against what the app actually reads. Section 8 lists the questions that only
> the app owner can answer.

## Quick start

1. Load the file for your omics and check `payload$payload_source` and
   `payload$payload_version` (`"2.0"`); file names are in section 1.
2. **Canonical keys** (section 3) are always in the list, `NULL` when they do not
   apply. Use `is.null()`. **Extension keys** (section 5) exist only for
   metabolomics; treat missing as `NULL`. Payloads from older runs lack the keys
   added later (section 7).
3. Three things differ by omics: what `expr_norm` is (4.1), what `de_final_table`
   contains (4.2), and how per-contrast columns in `de_stats` are named (4.3).
   Use `contrasts` and `de_summary` for the contrast list instead of parsing
   column names.
4. `de_final_table` is a clean data.frame in all three omics.
5. What changed compared with earlier payloads: section 7. What we need from
   you: section 8.

## 1. What you receive

Each omics run writes one RDS file next to its other outputs, named
`shiny_payload_<omics>_<run_name>.rds`, where `<omics>` is `rnaseq`,
`proteomics` or `metabolomics` and `<run_name>` is the run directory name
without its `Results_` prefix (project name and analysis round):

```r
payload <- readRDS("shiny_payload_metabolomics_Noam_Ziv_A07.rds")
payload$payload_source   # "rnaseq" | "proteomics" | "metabolomics"
payload$payload_version  # "2.0"
```

The payload is a named list. The contract is defined in code by
`get_payload_key_definitions()` in `R/core/07_shiny_contract.R`; the builders are
`build_shiny_payload_{rnaseq,proteomics,metabolomics}()`.

**R packages needed to load and draw the objects:** `ggplot2` (plus `ggrepel` for
the mummichog plots), `plotly` (`pca_3d`), `pheatmap`/`grid` (heatmaps), and
`limma` or `DESeq2` if you use `de_model`. Plots are stored as R objects, so it
is worth testing that they print in the app's environment (see section 8).

## 2. Conventions

- All keys are `snake_case`.
- **Canonical keys** (section 3) are always present in the list. When a key does
  not apply, or the step that produces it did not run, its value is `NULL`.
  `is.null(payload$key)` is the test to use.
- **Extension keys** (section 5) are specific to one omics. Treat a missing
  extension key and an extension key that is `NULL` the same way: not applicable.
- Payloads written by older pipeline versions do not contain keys added later
  (see section 7). Read every optional key defensively.
- The pipeline validates payloads with `assert_shiny_payload_contract()`, but
  only as warnings, not errors. A payload can therefore be written even when it
  deviates from the contract.
- Feature IDs are row names of `expr_raw`, `expr_norm` and `feature_annot`.
  Sample IDs are column names of the matrices and row names of `sample_meta`.

## 3. Canonical keys (all omics)

`Req` = the pipeline treats a `NULL` value as a contract violation.

### Execution info and metadata

| Key | Type | Req | Notes |
|---|---|---|---|
| `payload_version` | character | yes | `"2.0"` |
| `payload_created_at` | POSIXct | yes | |
| `payload_source` | character | yes | omics identifier |
| `sample_meta` | data.frame | yes | row names = sample IDs, matching `colnames(expr_norm)` |
| `feature_annot` | data.frame | no | row names = feature IDs; the ID column is **not** in the body. Every annotation column of the feature table flows through. Index by `rownames(feature_annot)`. |
| `contrasts` | data.frame | no | one row per contrast, with a `Contrast_name` column. `NULL` if no contrasts were supplied. |

### Expression

| Key | Type | Req | Notes |
|---|---|---|---|
| `expr_raw` | matrix | yes | filtered features x samples, **before** normalization; may contain `NA` |
| `expr_norm` | matrix | yes | features x samples, no `NA`. **What it is differs by omics**, see section 4. |
| `expr_long` | data.frame | no | long format of `expr_norm` joined with `sample_meta` (`feature_id`, `sample_id`, `value`, plus metadata columns) |

### QC / PCA

| Key | Type | Notes |
|---|---|---|
| `pca_object` | prcomp | variance explained, loadings |
| `pca_scores` | data.frame | scores with sample metadata |
| `pca_3d` | plotly | 3-D PCA widget |
| `imp_hist_samp` | ggplot | imputation histogram (proteomics; not expected in metabolomics) |
| `samples_hm` | pheatmap | sample-distance heatmap. In metabolomics this is the version **without** QC samples. |
| `samples_hm_w_na` | pheatmap | sample-distance heatmap with NA (proteomics) |

### Differential analysis

| Key | Type | Notes |
|---|---|---|
| `de_model` | model | fitted model (DESeq2 `dds` for RNA-seq, limma fit for proteomics/metabolomics). Can be large. |
| `de_stats` | data.frame | full statistics, all features and contrasts. Contains `feature_id` and `pass_any_contrast`. Per-contrast column names differ by omics, see section 4. |
| `de_sig_stats` | data.frame | rows of `de_stats` with `pass_any_contrast == 1`. If nothing passes it is a zero-row data.frame in RNA-seq and metabolomics, and `NULL` in proteomics. |
| `de_expr_norm` | matrix | rows of `expr_norm` for the significant features (same scale as `expr_norm`) |
| `de_summary` | data.frame | per-contrast counts: `contrast`, `up`, `down`, `total` |
| `de_final_table` | data.frame | DE table for display, see section 4 |
| `all_final_xlsx`, `de_final_xlsx` | raw | bytes of `Final_results_ALL_P_*.xlsx` / `Final_results_DE_P_*.xlsx`. Write with `writeBin()` to round-trip. The workbook is a presentation layout (sample and annotation rows above the table), not a plain table. |

### Clustering

| Key | Type | Notes |
|---|---|---|
| `clust_partition` | data.frame | partition-clustering assignment per feature |
| `clust_patterns` | list | binary clustering patterns |
| `clust_patterns_list` | character | labels of the binary patterns that have at least one feature (e.g. `"010"`) |
| `clust_heatmaps_by_pattern` | list | pheatmap objects per pattern |
| `clust_heatmap_hier` | list | hierarchical clustering: pheatmap result plus dendrogram data |
| `clust_heatmap_hier_fig` | gtable | drawable figure extracted from `clust_heatmap_hier` |
| `clust_heatmap_partition` | list | partition clustering heatmap |
| `clust_heatmap_partition_fig` | gtable | drawable figure extracted from `clust_heatmap_partition` |

### Configuration

| Key | Type | Notes |
|---|---|---|
| `padj_cutoff` | numeric | adjusted p-value threshold |
| `log_fc_cutoff` | numeric | log2 fold-change threshold (`log2` of the linear cutoff) |
| `norm_method` | character | RNA-seq/proteomics: the configured method. Metabolomics: `"<chosen_norm>/<transform>/<scaling>"`. |
| `group` | character | primary grouping variable (a column of `sample_meta`) |
| `color`, `shape` | character | aesthetic variables (columns of `sample_meta`); `color` can have several entries |

### Enrichment (RNA-seq only)

`enrichment` is a compact list (`available`, `gene_index`, `config`, `manifest`,
`ora`, `gsea`, `gsea_leading_edge`, `gsea_rankings`, `pathway_membership`) or
`NULL` when enrichment was not run. Genes are referenced by index into
`gene_index`; expression, `de_stats` and `feature_annot` are not duplicated in it.
The block is validated by `validate_enrichment_payload()`.

## 4. Where the omics differ

### 4.1 What `expr_norm` is

- **RNA-seq:** the normalized working matrix (`pre$expr_work`).
- **Proteomics:** a single imputation draw (`pre$expr_imp_single`), seeded so it
  is reproducible. The QC/PCA objects and the clustering objects in the payload
  are based on this matrix. The differential-abundance statistics are **not**
  computed on it: they come from separate multiple imputations that are pooled,
  so no single matrix in the payload is "the matrix the DE used". When
  clustering ran, the `<sample>.zscore` columns of `de_final_table` come from the
  clustering's z-scored matrix, which is derived from the same
  `pre$expr_imp_single`, so they are consistent with `expr_norm`. The
  `<sample>.norm` columns of `de_final_table` are different: they hold the
  imputation the DE model was fitted on (its first draw), and can differ from
  `expr_norm` for cells that were imputed.
- **Metabolomics:** the normalized working matrix (`pre$expr_work`) in which any
  `NA` left after preprocessing is filled with the median of that feature.
  **This fill is for display in the app only.** The differential analysis uses
  `pre$expr_work` with its missing values, not the filled matrix. The PCA in the
  payload uses the same row-median principle. The clustering handles `NA` in its
  own code path, so `expr_norm` should not be assumed to be the exact matrix the
  clustering ran on.
  The unfilled matrix is available as `metab_norm` (section 6), so
  `is.na(payload$metab_norm)` gives the positions of the original missing values.

### 4.2 `de_final_table`

A plain data.frame meant for use in the app, not a copy of the Excel layout.
It is built the same way in all three omics from the final-results data.frame
that also feeds the workbook: DE rows only (`pass_any_contrast == 1`), without
the `manual_cutoffs*` and `pass_any_contrast` columns. When the clustering order
is available it also contains `order` (rank of the feature in the clustering
order) and one `<sample>.zscore` column per sample. Rows are **not** sorted by
`order`. If no feature passes, it is a zero-row data.frame; if it could not be
built, or no final-results table was available, it is `NULL` (with a warning in
the first case). The workbook-only columns `Hierarchical_Order`, `Partition_*`
and `Binary_*`, and the sheet layout (sample and annotation rows above the
table), are not part of it.

- **RNA-seq:** the ID column is `Gene` or `FeatureID`.
- **Metabolomics:** the ID column is `feature_id`, followed by `original_id` (when
  available) and the other annotation columns; then one column per sample
  (values from the working matrix **before** the display fill, so `NA` where
  the measurement was missing); optional per-group `CV.<group>` columns
  (present only when enabled and valid for the chosen normalization); then the
  per-contrast statistics (`linearFC.`, `pvalue.`, `padj.`, `upDown.` plus the
  contrast name).
- **Proteomics:** the ID column is the one named by
  `modes$proteomics$de_table$id_col` in the config (default `FeatureID`). Besides
  annotation columns and the statistics per contrast, it holds one column per
  sample with the measured values (`NA` where not observed), the
  `<sample>.norm` block described in 4.1, and group mean, CV and pre-imputation
  summary columns.
- In all omics the ordering column is called `order` (in the workbook:
  `Hierarchical_Order`).

### 4.3 Per-contrast column names in `de_stats`

The column names that carry contrast statistics are not identical across omics.
As emitted by the builders:

| Omics | Pass column for a contrast | Fold-change column |
|---|---|---|
| RNA-seq | `<contrast>_pass` | `linearFC.<contrast>` |
| Proteomics | `pass.imputs.<contrast>` (contrast names have spaces removed) | one of `logFC_<c>`, `log2FoldChange_<c>`, `logFC.<c>`, `linearFC.imputs.<c>` |
| Metabolomics | `pass.<contrast>` | `linearFC.<contrast>` |

All three also have `feature_id` and `pass_any_contrast`. To avoid hard-coding
these patterns, `de_summary` (contrast, up, down, total) and `contrasts` give the
contrast list and the counts.

## 5. Metabolomics-specific extension keys

These are **not** canonical. They exist only in `payload_source == "metabolomics"`.

| Key | Content | Notes |
|---|---|---|
| `mummichog` | `mummichog[[contrast]] = list(title, subtitle, plot, table, slug)` | see 5.1 |
| `chosen_norm` | selected sample normalization (`NULL` in QC-review mode) | also encoded in `norm_method` |
| `samples_hm_w_qc` | sample-distance heatmap **with** QC samples | `samples_hm` is the version without QC |
| `rf_importance`, `rf_method` | random-forest importance table and method | |
| `plsda_vip_df`, `plsda_explained_variance` | PLS-DA VIP table and explained variance | |
| `enrichment_qea`, `enrichment_ssgsea`, `enrichment_ssgsea_scores`, `enrichment_ora`, `enrichment_gsea` | metabolomics enrichment tables | pending migration to the canonical `enrichment` block; names may change |
| `missingness`, `sample_map` | missingness summary; sample-ID mapping from the raw file | |

`rf_*`, `plsda_*`, `enrichment_*`, `missingness`, `sample_map` and the keys in
section 6 are written when the pipeline runs with `include_legacy = TRUE` (its
current setting; the name is historical). `mummichog` and `chosen_norm` are
always present (possibly `NULL`). `samples_hm_w_qc` does not depend on
`include_legacy`, but is present only when the corresponding QC heatmap is
available.

### 5.1 `mummichog`

- One entry per DE contrast **that has a mummichog result**; contrasts without a
  result are absent. The entry name is the original DE contrast label (the
  sanitised directory name is used only if the pipeline's mapping file is missing).
- The whole key is `NULL` when mummichog was not run or no contrast has a result.
- Fields of each entry:
  - `title`, `subtitle`: character, as used in the report.
  - `plot`: a `ggplot` (enrichment ratio vs. -log10 p, bubble size = pathway
    size, dashed line at the p cutoff, top pathways labelled when `ggrepel` is
    available where the plot is built).
  - `table`: data.frame with columns `Pathway`, `Overlap`, `Pathway size`,
    `Enrichment ratio`, `p.value`, sorted by `p.value`.
  - `slug`: filesystem-safe, de-duplicated token for the contrast (used for
    export file names).
- The plot and table are produced by the same presentation function as the
  metabolomics HTML report (`build_mummichog_report_sections()`), so they follow
  the same presentation logic and output. They are separately built objects, not
  the report's own object.
- Each ggplot carries its plotting environment, which increases the size of the
  payload file. If size becomes a concern, a more compact representation
  (data only) is possible; this has not been needed so far.

## 6. Compatibility and legacy keys (metabolomics)

Kept so existing consumers do not break. Not part of the contract; may be removed
after the app owner confirms they are unused.

| Key | Content | Note |
|---|---|---|
| `fc_cutoff` | the **linear** fold-change cutoff. Exists in all three omics. | derived: `log_fc_cutoff = log2(fc_cutoff)`. Retained, not canonical. |
| `metab_raw` | matrix before missingness filtering | **not** the same as `expr_raw` |
| `metab_filt` | filtered matrix | same content as `expr_raw` |
| `metab_norm` | working matrix **without** the display fill | the only place where the original `NA` positions of the working matrix are visible |
| `row_data` | feature annotation with the ID as a column | duplicate of `feature_annot` |

## 7. Changes compared with earlier payloads

If the app was written against an earlier payload, these are the differences.

1. **New canonical keys:** `contrasts`, `clust_patterns_list`,
   `clust_heatmap_partition_fig`. They were already emitted by the builders but
   were not part of the contract. Canonical keys now stay in the list as `NULL`
   when there is no value; previously some (for example `feature_annot`,
   `de_summary` or the QC plot keys) could be missing from the list entirely.
2. **Metabolomics `de_final_table`** was a copy of `de_sig_stats`. It is now the
   richer table described in 4.2. When no feature passes it is a zero-row
   data.frame instead of `NULL`.
3. **Proteomics `de_final_table`** was read from the `Results` sheet of the DE
   workbook, which is a presentation layout, so it was not a clean table. It is
   now built from the final-results data.frame like the other omics (4.2); the
   workbook-only columns `Hierarchical_Order`, `Partition_*` and `Binary_*` are
   no longer in it.
4. **Metabolomics `clust_patterns_list`** is now populated (it was missing).
5. **New metabolomics extension key `mummichog`.**
6. `payload_version` stays `"2.0"`: all changes are additive except items 2 and 3
   (the content of `de_final_table`).

## 8. Open items and questions for the app owner

Not verified against the app:

1. **Which keys does the app read?** A list would let unused extension and
   compatibility keys (section 6) be removed.
2. **Does the app hard-code contrast column patterns** such as `pass_<contrast>`
   or `logFC_<contrast>`? The three schemas in 4.3 differ, and the planned
   rename of internal differential-expression names to "differential abundance"
   for metabolomics (`docs/DE_to_DA_rename_map.md`) is blocked until this is known.
3. **Compatibility testing:** please check that `imp_hist_samp` and the
   `mummichog` plots print in the app's environment (they depend on the
   installed `ggplot2` and `ggrepel`).
4. **Empty results:** `de_sig_stats` is `NULL` for proteomics but a zero-row
   data.frame for RNA-seq and metabolomics when nothing passes. Which should the
   app expect? (To be aligned once known.)

Open on the pipeline side (not blocking this handoff):

- Enrichment for metabolomics and proteomics in the canonical `enrichment` key.
- Renaming the `include_legacy` switch.
