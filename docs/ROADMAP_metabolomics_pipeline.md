# Metabolomics Pipeline Roadmap

> Last updated: 2026-03-10

---

## Overview

The metabolomics pipeline is a production-grade, modular R framework built on `{targets}` for reproducible multi-omics analysis. It follows a two-stage architecture:

- **Stage 0 (`met_*` targets):** New preprocessing DAG — input loading, missingness, filtering, imputation, log transform, parallel normalization (TSS/Median/PQN), drift correction
- **Stage 1–2 (`metab_*` targets):** Analysis — differential expression, feature selection, pathway enrichment, QC reporting, Shiny export

**Key source locations:**

| Layer | Path | Files |
|-------|------|-------|
| Domain logic | `R/domain/metabolomics/` | 11 files (~4,600 lines) |
| Module wrappers | `R/modules/metabolomics/` | 7 files |
| Pipeline DAG | `R/pipeline/metabolomics/00_pipe_metabolomics.R` | 1 file (~458 lines) |
| Config template | `config/templates/metabolomics_template.yaml` | 1 file |
| Documentation | `metabolomics/README.md` | 1 file |

---

## Data Flow Diagram

```
met_input_files → met_raw (load & parse)
                      ├→ met_missingness_stats 
                      └→ met_filtered (sample + feature filtering)
                                  └→ met_log (log2 + pseudocount)
                                        ├→ met_norm_tss    → met_norm_tss_qc
                                        ├→ met_norm_median → met_norm_median_qc
                                        ├→ met_norm_pqn    → met_norm_pqn_qc
                                        └→ met_norm_comparison (TSV)
                                             └→ met_qc_comparison
                                                  └→ met_qc_summary_report (HTML)

                                   [chosen_norm selects one method]
                                        └→ met_corrected (drift correction)
                                              └→ metab_pre (adapter → legacy contract)
                                                    ├→ metab_qc_pre_obj (pre-DE diagnostics)
                                                    ├→ metab_de_res (limma/ttest/wilcoxon)
                                                    ├→ metab_feature_sel_res (RF + PLS-DA)
                                                    └→ metab_enrichment_res (fGSEA/ORA/QEA/ssGSEA)
                                                          └→ metab_standard_outputs
                                                          └→ metab_shiny_payload
                                                          └→ metab_report (HTML)
```

---

## Phase 1: Core Pipeline (Complete)

All capabilities below are implemented and operational.

### Input Loading (`00_inputs.R`)
- 3 input formats: `cd_raw` (Compound Discoverer), `processed_wide`, `multi_level`
- Optional sample filtering (QC/blank exclusion)
- Metadata alignment and validation

### Missingness Analysis (`08_missingness.R`)
- Spearman correlation between missingness & sample intensity
- MNAR/MAR classification (threshold: 0.3, configurable)
- Missingness heatmap visualization

### Filtering & Imputation (`09_imputation_met.R`)
- Sequential: sample missing threshold (30%) → feature missing threshold (20%)
- MNAR: min(observed)/2
- MAR: KNN imputation (k clamped to min(knn_k, ncol-1))

### Normalization (`01_normalization.R`)
- **TSS:** scale to median column sum (applied on linear scale → log2)
- **Median:** log-shift by sample median (applied on log scale)
- **PQN:** Probabilistic Quotient Normalization (applied on linear scale → log2)
- All 3 methods always computed in parallel; one selected via `chosen_norm`

### Drift Correction (`10_drift_correction.R`)
- LOESS correction for injection order bias
- Requires: `injection_order_col` + `is_QC` flag + ≥4 QC samples
- Silently skipped if preconditions not met

### Differential Expression (`03_differential.R`)
- Methods: limma, t-test, Wilcoxon
- Configurable p-value and fold-change cutoffs
- Per-contrast results with volcano/MA plots

### Feature Selection (`04_feature_selection.R`)
- Random Forest variable importance
- PLS-DA VIP scores

### Pathway Enrichment (`06_enrichment.R`)
- fGSEA, ORA, QEA, ssGSEA
- Pathway database support (KEGG, custom GMT)

### Reporting & Export
- HTML report (`mod_metabolomics_report`)
- CSV/TSV dataset outputs
- Shiny payload export (`07_shiny_export.R`)

---

## Phase 2: Two-Pass Normalization Review Workflow (Planned)

**Goal:** Allow researchers to review normalization QC before committing to a method.

### Design

Two execution modes controlled by `preprocessing.chosen_norm`:

| `chosen_norm` value | Mode | Behavior |
|---------------------|------|----------|
| `null` | **Review mode** | Run all 3 norms + QC suites + summary report → **stop** before `met_corrected` |
| `"pqn"` / `"tss"` / `"median"` | **Execution mode** | Skip QC comparison layer → run full pipeline |

### Workflow

```
Pass 1 (chosen_norm: null)
─────────────────────────
  met_raw → met_filtered → met_log
  ├→ met_norm_tss    → met_norm_tss_qc    (PCA, density, heatmaps)
  ├→ met_norm_median → met_norm_median_qc
  ├→ met_norm_pqn    → met_norm_pqn_qc
  └→ met_norm_comparison → met_qc_comparison → met_qc_summary_report (HTML)

  met_corrected → STOP with message:
    "chosen_norm is NULL — review mode. Set chosen_norm in config
     after reviewing qc/normalization_review_report.html, then re-run."

  (All downstream targets do not execute.)

Pass 2 (chosen_norm: "pqn")
───────────────────────────
  All norm targets cached from Pass 1 (no recomputation).

  SKIP: met_norm_comparison, met_*_qc, met_qc_comparison, met_qc_summary_report
  LOG:  "chosen_norm = 'pqn'; skipping normalization review QC."

  met_corrected → metab_pre → DE → feature selection → enrichment → report
```

### Implementation

| File | Change |
|------|--------|
| `R/modules/metabolomics/00_mod_preprocessing.R` | `mod_met_corrected()`: remove `%||% "pqn"` default; if `chosen_norm` is NULL, `stop()` with review-mode message |
| `R/modules/metabolomics/00_mod_preprocessing.R` | `mod_met_norm_comparison()`: if `chosen_norm` is set, return early with `cli::cli_inform()` |
| `R/modules/metabolomics/02_mod_qc_suite.R` | `mod_met_qc_suite()`, `mod_met_qc_comparison_table()`, `mod_met_qc_summary_report()`: check `chosen_norm`; if set, return early |
| `config/templates/metabolomics_template.yaml` | Change `chosen_norm` default to `null`; document two-pass workflow in comments |

---

## Phase 3: Legacy Normalization Evaluation Cleanup (Planned)

**Goal:** Remove `evaluate_normalization_methods()` and all references — fully superseded by `met_norm_comparison` + QC suite.

| File | What to remove |
|------|----------------|
| `R/domain/metabolomics/01_normalization.R` | `evaluate_normalization_methods()` function (~70 lines) |
| `R/domain/metabolomics/02_preprocess.R` | `normalization_eval` field in return list + call to `evaluate_normalization_methods()` |
| `R/domain/metabolomics/05_outputs_legacy.R` | `normalization_eval` TSV export block (lines ~37-41) |
| `R/domain/metabolomics/07_shiny_export.R` | `normalization_eval` payload block (lines ~301-303) |
| `config/templates/metabolomics_template.yaml` | `evaluate_methods` section (lines ~136-152) |

---

## Phase 4: Diagnostic Plots Directory Alignment (Planned)

**Goal:** Match proteomics/RNA-seq directory structure — separate pre-DE and post-DE QC outputs.

### Current state (metabolomics)

```
{out_dir}/Diagnostic_plots/     ← pre-DE QC + post-DE plots mixed together
{out_dir}/qc/                   ← normalization comparison (Phase 2 review mode)
```

### Target state (matching proteomics/RNA-seq)

```
{out_dir}/Diagnostic_plots/          ← pre-DE QC (PCA, density, heatmaps)
{out_dir}/Diagnostic_plots/QC_post/  ← post-DE (volcano, MA, top-DE heatmaps)
{out_dir}/qc/                        ← normalization comparison (review mode only)
```

### Implementation

| File | Change |
|------|--------|
| `R/modules/metabolomics/02_mod_de.R` | Change `out_qc <- dirs$diagnostic_plots` to `out_qc <- file.path(dirs$diagnostic_plots, "QC_post")` + `ensure_dir(out_qc)` |

This is a single-line change. No new files needed — reuses the existing `QC_post` pattern from:
- `R/modules/proteomics/03_mod_qc_post.R` (line 26)
- `R/modules/rnaseq/03_mod_qc_post.R` (line 14)

---

## Phase 5: Advanced Multi-Omics Integration (Planned)

> These are **multi-omics level** features in `R/domain/multiomics/`, not metabolomics-specific.

| Feature | Domain file | Status |
|---------|-------------|--------|
| MOFA2 integration | `04b_integration_mofa.R` | Code exists; needs config validation + pipeline wiring |
| Foundational cross-omics correlations | `09_foundational_correlations.R` | Code exists; needs pipeline wiring |
| Integration consensus | `10_integration_consensus.R` | Code exists; needs pipeline wiring |
| Stability / bootstrap validation | `11_stability_analysis.R` | Code exists; needs pipeline wiring |
| Extended RNA-protein correlation | `06b_rna_protein_correlation.R` | Code exists; needs pipeline wiring |

---

## Phase 6: Interpretation & Automation (Future)

| Feature | Domain file | Status |
|---------|-------------|--------|
| AI figure commentary | `13_commentary.R` | Core services partially implemented |
| KEGG Pathview overlays | (multiomics domain) | Domain code exists |
| Mechanistic inference | (multiomics domain) | Domain code exists |

---

## Configuration Reference

Key config knobs under `modes.metabolomics`:

| Key | Default | Purpose |
|-----|---------|---------|
| `input.format` | `"cd_raw"` | Input format: `cd_raw`, `processed_wide`, `multi_level` |
| `preprocessing.feat_missing_threshold` | `0.20` | Max fraction of missing values per feature |
| `preprocessing.sample_missing_threshold` | `0.30` | Max fraction of missing values per sample |
| `preprocessing.mnar_threshold` | `0.3` | Spearman correlation threshold for MNAR classification |
| `preprocessing.knn_k` | `10` | KNN imputation k parameter |
| `preprocessing.chosen_norm` | `"pqn"` (→ `null` after Phase 2) | Normalization method selection; `null` = review mode |
| `normalization.transform` | `"log2"` | Log transformation method |
| `normalization.pseudocount` | `1` | Pseudocount for log transform |
| `preprocessing.drift_correction.enabled` | `true` | LOESS drift correction toggle |
| `de.method` | `"limma"` | DE method: `limma`, `ttest`, `wilcoxon` |
| `de.p_cutoff` | `0.05` | Significance threshold |
| `de.linear_fc_cutoff` | `1.5` | Fold-change cutoff (linear scale) |

---

## File Index

### Domain (`R/domain/metabolomics/`)

| File | Lines | Purpose |
|------|-------|---------|
| `00_inputs.R` | 714 | Load & parse metabolomics data |
| `01_normalization.R` | 336 | TSS, Median, PQN normalization |
| `02_preprocess.R` | 236 | Legacy preprocessing orchestrator |
| `03_differential.R` | 421 | DE analysis (limma, t-test, Wilcoxon) |
| `04_feature_selection.R` | 360 | RF & PLS-DA feature selection |
| `05_outputs_legacy.R` | 43 | Legacy output compatibility |
| `06_enrichment.R` | 1,279 | Pathway enrichment (fGSEA, ORA, QEA, ssGSEA) |
| `07_shiny_export.R` | 471 | Shiny app export |
| `08_missingness.R` | 214 | MNAR/MAR classification |
| `09_imputation_met.R` | 242 | KNN + min/2 imputation |
| `10_drift_correction.R` | 284 | LOESS drift correction |

### Modules (`R/modules/metabolomics/`)

| File | Purpose |
|------|---------|
| `00_mod_preprocessing.R` | New met_* preprocessing DAG orchestration |
| `01_mod_qc_pre.R` | PCA, heatmaps, sample distance QC |
| `02_mod_de.R` | Differential expression wrapper |
| `02_mod_qc_suite.R` | Comprehensive QC suite (PCA, density, outliers) |
| `03_mod_feature_selection.R` | RF & PLS-DA wrappers |
| `05_mod_enrichment.R` | Pathway enrichment orchestration |
| `06_mod_report.R` | HTML report generation |

### Pipeline

| File | Purpose |
|------|---------|
| `R/pipeline/metabolomics/00_pipe_metabolomics.R` | {targets} DAG definition |

### Config

| File | Purpose |
|------|---------|
| `config/templates/metabolomics_template.yaml` | Full config reference |
| `config/config_GT4.yaml` | Genotype 4 analysis |
| `config/config_GT6.yaml` | Genotype 6 analysis |
| `config/config_GT15.yaml` | Genotype 15 analysis |
