# Lipidomics pipeline — pending enhancements plan

Branch: `Lipidomics_Oz` (created from `Lipidomics_dor`, 2026-05-12)
Source request: extend the R `{targets}` lipidomics pipeline with a group-aware
feature filter, pool-CV QC, DE-aware structural plots, and a top-DE bar — and
wire everything into `R/pipeline/lipidomics/templates/report_lipidomics.Rmd`.

This file freezes the plan we agreed on before implementation, so the work can
be picked up later. See `lipidomics_amir_sapir_run_memo.md` for L01 run state
and known issues.

## Audit — what already exists (do NOT duplicate)

| Step | Already in code | Section in `report_lipidomics.Rmd` |
|---|---|---|
| Volcano per contrast | `mod_lipidomics_de` → `plot_volcano` | Differential Expression (L911) |
| Top-DE heatmap (z-scored, hierarchical) | clustering module | Top DE Heatmaps (L1727) |
| Class pie / stacked bars | `plot_lipid_category_donut`, `plot_lipid_class_barplot` | Category Donut (L830), Class Composition (L1980) |
| Chain length / saturation (class-wise) | `compute_chain_saturation`, `compute_chain_length_distribution`, `plot_chain_*` in `05_lipid_class_analysis.R` | Chain Saturation/Length (L2002, L2026) |
| Class enrichment (Fisher ORA) | `lipid_class_ora` | Class Enrichment (L2052) |
| PCA PC1×PC2 + scree + loading | `01_mod_qc_pre` + `plot_pca_scree/loading` (`08_qc_enhancements.R`) | QC Diagnostics → PCA (L597), Scree (L810), Loading (L820) |
| PERMANOVA | `run_permanova` | PERMANOVA (L858) |
| Density / boxplots / distance heatmap | `01_mod_qc_pre` | Density / Boxplots / Heatmaps (L698–L761) |

## What is genuinely missing — the work to do

| # | Item | Where it goes |
|---|---|---|
| 1 | **Group-aware missingness filter** — drop lipids with >`max_group_missing_frac` NA/zero in **every** biological group; keep if robust in ≥1 group. Runs *after* sample filter, *before* normalization. Config key: `modes.lipidomics.feature_filter.max_group_missing_frac` (default `0.30`). | `R/domain/lipidomics/02_preprocess.R`: new `filter_features_by_group_missingness(expr, meta, group_col, max_frac)`; call inside `preprocess_lipidomics()` after the sample-filter block. |
| 2 | **Pool CV** — per-feature CV across `^Pool` samples computed *before* `exclude_pools` strips them. Median pool CV reported (expect <20–30%). | New file `R/domain/lipidomics/01c_pool_qc.R`: `compute_pool_cv(pool_matrix, ...)` → tibble + `median_cv`. Modify `preprocess_lipidomics()` to stash the pool sub-matrix on the result list **before** the pool filter step, so this stays a pure function. |
| 3 | **Pool CV density plot** | Same `01c_pool_qc.R`: `plot_pool_cv_density(cv_df)`. New target `lipid_pool_cv` between `lipid_pre` and `lipid_qc`. Skip silently if no pools detected. |
| 5 | **DE-aware structural bubble** — chain length × DB count, one point per significant lipid, size = `-log10(p)`, colour = logFC sign. Distinct from the existing class-wise chain distributions. | `R/domain/lipidomics/05_lipid_class_analysis.R`: `plot_de_chain_bubble(de_tbl, row_data, sig_cutoff)`. Reuse the existing `FAKey` parsing in `extract_individual_chains()` / chain-saturation helpers. One plot per contrast. |
| 6 | **PCA PC1 vs PC3 tab** — `pca_obj` already carries the components; just a second render. | `R/modules/lipidomics/01_mod_qc_pre.R`: add `pca_13` plot alongside existing PC1×PC2. Surface as sibling tab under QC → PCA in the Rmd. |
| 7 | **Report wiring** — single Rmd pass after all data targets are populated. | `R/pipeline/lipidomics/templates/report_lipidomics.Rmd`: <br/>• **QC → Pool CV** sub-tab after PCA scree (skip if `lipid_pool_cv` is NULL)<br/>• **Structural (DE) Bubble** as sub-tab under "Lipid Class Analysis"<br/>• **PC1 vs PC3** sibling tab under PCA |

## Implementation order (dependency-driven)

1. Feature-missingness filter (data shape change — must come first; everything downstream re-runs)
2. Pool-CV stash in `preprocess_lipidomics()` + `compute_pool_cv()` + `plot_pool_cv_density()` + `lipid_pool_cv` target
3. `plot_de_chain_bubble()` (depends on `de_res` + `row_data`; no preprocess changes)
4. PC1×PC3 in QC module
5. Report Rmd wiring (single pass after all data targets exist)

(Step 4 "Top-N DE bar chart per contrast" was dropped from scope on user request, 2026-05-12.)

## Decisions confirmed before coding

1. **Pool-CV gating.** Compute pool CV *upstream* of the existing `exclude_pools: true` filter rather than introducing a "keep pools through QC" config flag. Less surface area; the existing exclude-pools behavior is preserved.
2. **Default `max_group_missing_frac = 0.30`** per the user's request. Configurable under `modes.lipidomics.feature_filter`.

## Out of scope (explicitly NOT in this batch)

- Switching to a Python pipeline (rejected — stay in R, consistent with the rest of the project per team coding policy).
- Re-writing the QC report Rmd (`lipidomics_qc_report.Rmd`). New plots go to the **main** `report_lipidomics.Rmd`.
- Touching the open `plsda_confusion_matrix.png` issue — that fix already landed on this branch (`R/domain/lipidomics/08_qc_enhancements.R:255`, switched to explicit S3 method lookup).

## Related changes already on `Lipidomics_Oz`

These are committable on this branch but separate from the plan above:

| File | Change | Status |
|---|---|---|
| `R/domain/lipidomics/08_qc_enhancements.R` | mixOmics S3 `predict` dispatch fix (memo "open issue") | Applied |
| `_targets.R` | `config_file` target now bakes resolved `MULTIOMICS_CONFIG` into its command via `tar_target_raw` + `bquote`, so switching configs invalidates cache | Applied |
| `config/Lipidomics_Amir_Sapir2.yaml` | `claude_model: claude-sonnet-4.5-20250514` → `claude-sonnet-4-5-20250929` (matches code fallback in `R/services/12_commentary.R:749`) | Applied |
