# Migration Plan: multiomics_pipeline → multiomics-core

**Source:** `/home/ozsol/multiomics/src/pipelines/multiomics_pipeline/`  
**Destination:** `/home/ozsol/multiomics-core/`  
**Date:** 2026-02-28  
**Goal:** Transfer remaining multiomics pipeline functionality with **minimum code changes**.

---

## 1. Architecture Comparison

### Source Pipeline (standalone, flat structure)
```
R/
  00_utils.R → 14_rna_protein_correlation.R   (22 files, ~14,630 lines)
_targets.R                                     (monolithic pipeline definition)
config.yml                                     (flat config structure)
scripts/                                       (Python helpers)
reports/                                       (Rmd/Qmd templates)
```

### Target: multiomics-core (layered 3-tier architecture)
```
R/core/       → shared utilities (io, config, validation, annotation, etc.)
R/domain/     → pure analytical functions (stateless, per-omics subdirs)
R/modules/    → orchestration wrappers (combine domain functions + I/O)
R/pipeline/   → {targets} DAG definitions
R/services/   → external services (AI commentary)
```

### Key Architectural Differences

| Aspect | Source Pipeline | multiomics-core |
|--------|---------------|-----------------|
| Config path | `config$integration$methods` | `config$modes$multiomics$integration$methods` |
| Config path | `config$design$condition_column` | `config$design$condition_column` (same) |
| Config path | `config$output$output_dir` | `config$paths$out` + `run_dir` |
| Config path | `config$feature_selection$transcriptomics` | `config$modes$multiomics$feature_selection` |
| Config path | `config$enrichment$*` | `config$modes$multiomics$enrichment` |
| Config loader | `load_config()` in `00_utils.R` | `load_config()` in `core/04_config.R` |
| Null coalesce | `%||%` in `00_utils.R` | `%||%` in `core/04_config.R` |
| Logging | `log_message()` in `00_utils.R` | `message()` (standard R) |
| Plot saving | `save_plot()` in `00_utils.R` | `ggsave()` directly |
| Data structure | `mae_data` (custom list w/ `harmonized_omics$*`) | `mae` (formal `MultiAssayExperiment`) |
| Pipeline def | Monolithic `_targets.R` | `pipe_multiomics()` function |
| Sourcing | Explicit `source("R/11b_*.R")` calls | Auto-sourced by `tar_source()` |

---

## 2. What Already Exists in multiomics-core

| Capability | Source File(s) | Core File(s) | Status |
|-----------|---------------|-------------- |--------|
| Config loading/validation | `00_utils.R` | `core/04_config.R`, `domain/multiomics/90_config_validate.R` | ✅ Done |
| Data ingestion | `01_ingestion.R` | Per-omics `domain/*/00_inputs.R` | ✅ Done (per-omics pipelines handle this) |
| Per-omics preprocessing | `02_preprocessing.R` | Per-omics `domain/*/02_preprocess.R` etc. | ✅ Done |
| ID mapping (gene↔protein) | `03_mapping.R` | `domain/multiomics/01_id_mapping.R` | ✅ Done |
| MAE construction | `04_harmonize.R` | `domain/multiomics/02_mae.R` | ✅ Done |
| Feature selection | `05_feature_selection.R` | `domain/multiomics/03_feature_selection.R` | ✅ Done |
| DIABLO integration | `07_diablo.R` | `domain/multiomics/04_integration_diablo.R` | ✅ Done |
| SNF integration | `08_snf.R` | `domain/multiomics/05_integration_snf.R` | ✅ Done |
| Concordance (basic) | `09_concordance.R` | `domain/multiomics/06_concordance.R` | ✅ Done |
| Cross-omics enrichment | `10_enrichment.R` (partial) | `domain/multiomics/07_enrichment.R` | ✅ Done (Fisher meta-analysis) |
| Report | Rmd template | `domain/multiomics/08_report.R` | ✅ Done (placeholder) |
| Pipeline DAG | `_targets.R` | `pipeline/multiomics/00_pipe_multiomics.R` | ✅ Done |
| Module wrappers | N/A | `modules/multiomics/01-05_mod_*.R` | ✅ Done |

---

## 3. What Is MISSING from multiomics-core

These are substantial capabilities in the source pipeline with **no equivalent** in multiomics-core:

| # | Capability | Source File | Lines | Priority |
|---|-----------|-------------|-------|----------|
| **M1** | MOFA2 integration | `06_mofa.R` + `mofa_wrapper.R` | 736 | **High** |
| **M2** | Foundational cross-omics analysis | `05b_foundational_correlations.R` | 1,880 | **Medium** |
| **M3** | Mechanistic inference | `05c_mechanistic_inference.R` | 1,403 | **Low** |
| **M4** | Integration consensus | `09b_integration_consensus.R` | 1,041 | **Medium** |
| **M5** | Stability analysis | `09c_stability_analysis.R` | 1,096 | **Low** |
| **M6** | AI figure commentary | `11_commentary.R` + `11b_commentary_fallbacks.R` | 1,182 | **Low** (already in `services/12_commentary.R` for other omics) |
| **M7** | MultiGSEA plots | `13_multigsea_plots.R` | 401 | **Medium** |
| **M8** | RNA-protein correlation (extended) | `14_rna_protein_correlation.R` | 343 | **Medium** |
| **M9** | KEGG Pathview overlays | `15_kegg_pathview.R` | 770 | **Medium** |
| **M10** | Per-omics enrichment (full w/ collections) | `10_enrichment.R` (per-contrast, multi-collection) | 1,593 | **Medium** |
| **M11** | Rmd report template | `reports/analysis_report.Rmd` | ~400 | **Low** |
| **M12** | Python MOFA helper | `scripts/run_mofa.py` | ~200 | **High** (if M1) |
| **M13** | Python commentary scripts | `scripts/figure_commentary_*.py` | ~200 | **Low** |

**Total new code: ~9,445 lines** to port.

---

## 4. Detailed Migration Plan

### Phase 1: MOFA2 Integration (High Priority)

**Files to create:**
- `R/domain/multiomics/04b_integration_mofa.R` — Domain logic
- `scripts/run_mofa.py` — Python helper (copy)

**Source files to adapt:**
- `06_mofa.R` → `R/domain/multiomics/04b_integration_mofa.R`
- `mofa_wrapper.R` → Merge into `04b_integration_mofa.R` (or separate `R/core/15_mofa_wrapper.R`)

**Required changes (minimal):**

1. **Config path adaptation:** Replace all `config$integration$mofa2$*` with `config$modes$multiomics$integration$mofa2$*`

2. **Data structure adaptation:** Source uses `feature_data$filtered_mae` as a custom list with `extract_matrices_for_integration()`. Core uses formal `MultiAssayExperiment`. Need a thin adapter:
   ```r
   # Adapter: extract matrices from formal MAE
   extract_matrices_from_mae <- function(mae) {
       lapply(experiments(mae), function(x) as.matrix(assay(x)))
   }
   ```

3. **Replace `log_message()` → `message()`** (or alias at top of file)

4. **Replace `save_plot()`/`save_table()`** → standard `ggsave()`/`write.csv()` with explicit path

5. **Wire into pipeline targets:**
   ```r
   # In R/modules/multiomics/02_mod_integration.R — add MOFA2 block
   if ("MOFA2" %in% methods) {
       results$mofa <- run_mofa2_integration(mae = mae_subset, config = config, out_dir = mofa_dir)
   }
   ```

6. **Update config validator** in `90_config_validate.R` to accept MOFA2 params.

7. **Copy** `scripts/run_mofa.py` to `scripts/run_mofa.py`.

**Estimated changes: ~50 lines modified, ~750 lines copied as-is.**

---

### Phase 2: RNA-Protein Correlation & Concordance Enhancement (Medium Priority)

**Files to create:**
- `R/domain/multiomics/06b_rna_protein_correlation.R`

**Source file to adapt:**
- `14_rna_protein_correlation.R` → `R/domain/multiomics/06b_rna_protein_correlation.R`

**Required changes:**

1. **Data access:** Replace `mae_data$harmonized_omics$transcriptomics$normalized_matrix` with `assay(mae[["transcriptomics"]])` (formal MAE access)

2. **DE table access:** Replace `mae_data$harmonized_omics$*$de_table` with the separate `de_results` object already available in the pipeline

3. **Config paths:** Replace `config$output$output_dir` with explicit `out_dir` parameter

4. **Wire into pipeline:** Add a new target in `pipe_multiomics()` after concordance

**Estimated changes: ~30 lines modified, ~310 lines copied as-is.**

---

### Phase 3: Foundational Cross-Omics Analysis (Medium Priority)

**Files to create:**
- `R/domain/multiomics/09_foundational_correlations.R`
- `R/modules/multiomics/06_mod_foundational.R`

**Source file to adapt:**
- `05b_foundational_correlations.R` → `R/domain/multiomics/09_foundational_correlations.R`

**Required changes:**

1. **Data structure:** Replace all `mae_data$harmonized_omics$*$normalized_matrix` with formal MAE access patterns. Create a helper to extract harmonized omics:
   ```r
   extract_harmonized_omics <- function(mae) {
       lapply(names(experiments(mae)), function(nm) {
           list(normalized_matrix = as.matrix(assay(mae[[nm]])))
       }) |> setNames(names(experiments(mae)))
   }
   ```

2. **Replace `gene_mapping`** references with `gene_protein_mapping` from harmonization

3. **Config:** Replace `config$foundational$*` → `config$modes$multiomics$foundational$*`

4. **Replace `save_plot()`/`save_table()`/`log_message()`** with standard equivalents

5. **Wire into pipeline:** Add `multiomics_foundational` target after harmonization, before integration

**Estimated changes: ~80 lines modified, ~1,800 lines copied as-is. Most of the code is self-contained analytical logic.**

---

### Phase 4: Integration Consensus & Stability (Medium-Low Priority)

**Files to create:**
- `R/domain/multiomics/10_integration_consensus.R`
- `R/domain/multiomics/11_stability_analysis.R`
- `R/modules/multiomics/07_mod_consensus.R`

**Source files to adapt:**
- `09b_integration_consensus.R` → `R/domain/multiomics/10_integration_consensus.R`
- `09c_stability_analysis.R` → `R/domain/multiomics/11_stability_analysis.R`

**Required changes:**

1. **Integration results format:** Source uses `integration_results$mofa$results$*`, `integration_results$diablo$model$*`. Core's `mod_multiomics_integration` returns `diablo_results`/`snf_results`/`mofa_results`. Need wrapper:
   ```r
   integration_results <- list(
       mofa = integration$mofa_results,
       diablo = integration$diablo_results,
       snf = integration$snf_results
   )
   ```

2. **Config paths:** `config$consensus$*` → `config$modes$multiomics$consensus$*` and `config$stability$*` → `config$modes$multiomics$stability$*`

3. **Standard replacements:** `log_message()` → `message()`, `save_plot()` → `ggsave()`, `save_table()` → `write.csv()`

4. **Wire into pipeline:** Add `multiomics_consensus` and `multiomics_stability` targets after integration

**Estimated changes: ~60 lines modified, ~2,000 lines copied as-is.**

---

### Phase 5: Enhanced Enrichment & Pathview (Medium Priority)

**Files to create:**
- `R/domain/multiomics/07b_multigsea_plots.R`
- `R/domain/multiomics/07c_kegg_pathview.R`

**Source files to adapt:**
- `13_multigsea_plots.R` → `R/domain/multiomics/07b_multigsea_plots.R`
- `15_kegg_pathview.R` → `R/domain/multiomics/07c_kegg_pathview.R`

**Required changes:**

1. **Config:** `config$enrichment$multigsea$*` → `config$modes$multiomics$enrichment$multigsea$*`

2. **Data access:** Adapt to formal MAE + separate `enrichment_results`

3. **KEGG Pathview:** The `15_kegg_pathview.R` is largely standalone. Only needs:
   - Remove `install_load_packages()` (packages managed externally)
   - Replace hard-coded `data/HMDB2kegg_cpd.Jan2026.v2.txt` with config-driven path

4. **Copy** `data/HMDB2kegg_cpd.Jan2026.v2.txt` if not already present (it already exists in multiomics-core `data/`)

**Estimated changes: ~40 lines modified, ~1,100 lines copied as-is.**

---

### Phase 6: Mechanistic Inference (Low Priority)

**Files to create:**
- `R/domain/multiomics/12_mechanistic_inference.R`
- `R/modules/multiomics/08_mod_mechanistic.R`

**Source file to adapt:**
- `05c_mechanistic_inference.R` → `R/domain/multiomics/12_mechanistic_inference.R`

**Required changes:** Same pattern as Phase 3 (data access, config path, logging). Code is highly self-contained.

**Estimated changes: ~50 lines modified, ~1,350 lines copied as-is.**

---

### Phase 7: AI Commentary & Report Template (Low Priority)

**Files to create:**
- `R/domain/multiomics/13_commentary.R`
- `R/domain/multiomics/13b_commentary_fallbacks.R`
- `scripts/figure_commentary_claude.py` (if not present)
- `scripts/figure_commentary_openai.py` (if not present)

**Note:** multiomics-core already has `R/services/12_commentary.R` for per-omics commentary. The multiomics-specific commentary from the source adds multi-omics figure types (MOFA, DIABLO, concordance, etc.). This should extend the existing service rather than duplicate it.

---

## 5. Config Migration Cheatsheet

All source `config$X` references need remapping to `config$modes$multiomics$X`:

| Source Config Path | multiomics-core Equivalent |
|---|---|
| `config$integration$methods` | `config$modes$multiomics$integration$methods` |
| `config$integration$mofa2$*` | `config$modes$multiomics$integration$mofa2$*` |
| `config$integration$diablo$*` | `config$modes$multiomics$integration$diablo$*` |
| `config$integration$snf$*` | `config$modes$multiomics$integration$snf$*` |
| `config$feature_selection$transcriptomics$*` | `config$modes$multiomics$feature_selection$*` (unified) |
| `config$enrichment$*` | `config$modes$multiomics$enrichment$*` |
| `config$foundational$*` | `config$modes$multiomics$foundational$*` (new) |
| `config$mechanistic$*` | `config$modes$multiomics$mechanistic$*` (new) |
| `config$consensus$*` | `config$modes$multiomics$consensus$*` (new) |
| `config$stability$*` | `config$modes$multiomics$stability$*` (new) |
| `config$commentary$*` | `config$modes$multiomics$commentary$*` (new) |
| `config$output$output_dir` | `out_dir` parameter (passed in) |
| `config$design$condition_column` | `config$design$condition_column` (same!) |
| `config$design$contrasts` | `config$design$contrasts` (same!) |
| `config$global$organism` | `config$modes$multiomics$organism` or detect from sub-pipelines |

---

## 6. Shared Utility Deduplication

These functions exist in both codebases and should **NOT** be re-ported:

| Function | Source | Already in Core |
|----------|--------|-----------------|
| `load_config()` | `00_utils.R` | `core/04_config.R` |
| `%\|\|%` | `00_utils.R` | `core/04_config.R` |
| `read_gmt()` | `00_utils.R` | `core/13_gmt_utils.R` |
| Gene annotation/mapping | `03_mapping.R` | `core/11_annotation.R`, `core/10_organism_detection.R` |
| `build_gene_protein_mapping()` | `03_mapping.R` | `domain/multiomics/01_id_mapping.R` |
| `create_multiassay_experiment()` | `04_harmonize.R` | `domain/multiomics/02_mae.R` |
| `select_features_for_integration()` | `05_feature_selection.R` | `domain/multiomics/03_feature_selection.R` |
| `run_diablo_integration()` | `07_diablo.R` | `domain/multiomics/04_integration_diablo.R` |
| `run_snf_integration()` | `08_snf.R` | `domain/multiomics/05_integration_snf.R` |
| Concordance (basic) | `09_concordance.R` | `domain/multiomics/06_concordance.R` |
| Fisher combined p-values | `00_utils.R` | `domain/multiomics/07_enrichment.R` |

---

## 7. Data Compatibility Bridge

The source pipeline uses a custom `mae_data` list structure, while multiomics-core uses a formal `MultiAssayExperiment`. To minimize changes, create a **thin adapter** in `R/core/` or at the top of migrated files:

```r
#' Convert formal MAE to the source pipeline's mae_data list format
#' This allows porting source code with minimal changes to data access patterns.
mae_to_legacy_list <- function(mae, gene_protein_mapping = NULL, metadata = NULL) {
    harmonized_omics <- lapply(names(experiments(mae)), function(nm) {
        exp <- experiments(mae)[[nm]]
        list(
            normalized_matrix = as.matrix(assay(exp)),
            de_table = if (!is.null(rowData(exp)$logFC)) as.data.frame(rowData(exp)) else NULL,
            feature_annotation = as.data.frame(rowData(exp))
        )
    })
    names(harmonized_omics) <- names(experiments(mae))

    list(
        mae = mae,
        harmonized_omics = harmonized_omics,
        metadata = metadata %||% as.data.frame(colData(mae)),
        common_samples = colnames(mae),
        gene_mapping = gene_protein_mapping
    )
}

#' Extract matrices for integration (formal MAE → named list of matrices)
extract_matrices_for_integration <- function(mae_data, omics_to_use = NULL) {
    if (inherits(mae_data, "MultiAssayExperiment")) {
        nm <- omics_to_use %||% names(experiments(mae_data))
        lapply(nm, function(n) as.matrix(assay(mae_data[[n]]))) |> setNames(nm)
    } else {
        # Legacy list format
        mae_data$harmonized_omics
    }
}
```

**This avoids rewriting every `mae_data$harmonized_omics$*$normalized_matrix` access.** Instead, convert once at the module boundary.

---

## 8. Recommended Execution Order

```
Phase 1 (MOFA2)               →  ~1 day   |  Unlocks full 3-method integration
Phase 2 (RNA-Protein corr.)   →  ~0.5 day |  Independent, quick win
Phase 5 (Enrichment/Pathview) →  ~0.5 day |  Independent, extends enrichment
Phase 3 (Foundational)        →  ~1 day   |  Large file, mostly copy
Phase 4 (Consensus/Stability) →  ~1 day   |  Depends on Phase 1
Phase 6 (Mechanistic)         →  ~0.5 day |  Independent, low priority
Phase 7 (Commentary/Report)   →  ~0.5 day |  Lowest priority, extending existing service
```

**Total estimated effort: ~5 days**

---

## 9. Config Template Addition

Add these new sections to `config/multiomics_GT15_test.yaml` (or equivalent):

```yaml
  multiomics:
    # ... existing fields ...

    integration:
      methods:
        - "DIABLO"
        - "SNF"
        - "MOFA2"          # NEW

      mofa2:               # NEW
        num_factors: 10
        convergence_mode: "fast"
        seed: 42

      diablo:
        ncomp: 3
        design_matrix: "full"
        cv_folds: 3
        cv_repeats: 1

      snf:
        K: 15
        alpha: 0.5
        T: 20
        n_clusters: 2

    # NEW sections
    foundational:
      run_foundational: true
      correlation_method: "spearman"
      top_variable_features: 500

    consensus:
      compare_methods: true

    stability:
      run_stability: false    # Expensive, opt-in
      n_bootstrap: 100

    enrichment:
      run_enrichment: true
      methods: "gsea"
      collections:
        - name: "GO_BP"
          type: "GO"
          ont: "BP"
        - name: "KEGG"
          type: "KEGG"
      multigsea:
        run_multigsea: true
        pvalue_threshold: 0.05
      pathview:
        run_pathview: true
        top_n: 5

    commentary:
      enabled: false
      backend: "none"
```

---

## 10. Files to Copy As-Is (No Changes)

| Source | Destination | Notes |
|--------|-------------|-------|
| `scripts/run_mofa.py` | `scripts/run_mofa.py` | Python MOFA helper |
| `scripts/figure_commentary_claude.py` | Already exists | Verify identical |
| `scripts/figure_commentary_openai.py` | Already exists | Verify identical |
| `data/HMDB2kegg_cpd.Jan2026.v2.txt` | Already exists | Verify identical |
| `reports/analysis_report.Rmd` | `reports/multiomics_report.Rmd` | Adapt for core |

---

## 11. Testing Strategy

1. **Unit tests:** For each migrated domain function, add test in `tests/testthat/test-multiomics-<module>.R`
2. **Integration test:** Run full pipeline with `config/multiomics_GT15_test.yaml`
3. **Regression check:** Compare outputs from source pipeline vs. multiomics-core on same data
4. **Validation script:** Extend `validate_multiomics.R` to cover new targets

---

## 12. Summary of Changes Per File

### New Files to Create (8-14 files)

| File | Source | Lines (approx) |
|------|--------|----------------|
| `R/domain/multiomics/04b_integration_mofa.R` | `06_mofa.R` + `mofa_wrapper.R` | ~700 |
| `R/domain/multiomics/06b_rna_protein_correlation.R` | `14_rna_protein_correlation.R` | ~340 |
| `R/domain/multiomics/07b_multigsea_plots.R` | `13_multigsea_plots.R` | ~400 |
| `R/domain/multiomics/07c_kegg_pathview.R` | `15_kegg_pathview.R` | ~770 |
| `R/domain/multiomics/09_foundational_correlations.R` | `05b_foundational_correlations.R` | ~1,880 |
| `R/domain/multiomics/10_integration_consensus.R` | `09b_integration_consensus.R` | ~1,040 |
| `R/domain/multiomics/11_stability_analysis.R` | `09c_stability_analysis.R` | ~1,096 |
| `R/domain/multiomics/12_mechanistic_inference.R` | `05c_mechanistic_inference.R` | ~1,400 |
| `R/modules/multiomics/06_mod_foundational.R` | New | ~40 |
| `R/modules/multiomics/07_mod_consensus.R` | New | ~40 |
| `scripts/run_mofa.py` | Copy from source | ~200 |

### Existing Files to Modify (4 files, ~100 lines total)

| File | Change |
|------|--------|
| `R/modules/multiomics/02_mod_integration.R` | Add MOFA2 integration block (~15 lines) |
| `R/pipeline/multiomics/00_pipe_multiomics.R` | Add targets for foundational, MOFA, consensus, stability, pathview (~60 lines) |
| `R/domain/multiomics/90_config_validate.R` | Add MOFA2/foundational/consensus/stability config validation (~25 lines) |
| `config/multiomics_GT15_test.yaml` | Add new config sections (~40 lines) |

**Total: ~8,000 lines of new code (mostly copied), ~100 lines of modifications to existing files.**
