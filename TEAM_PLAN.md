# Multiomics-Core: Team Roadmap & Development Plan

> **Version:** 1.0
> **Date:** 2026-02-12
> **Status:** v0.2.1 — Proteomics complete, RNA-seq functional, Metabolomics scaffolded
> **Architecture:** R / {targets} pipeline, contract-based, config-driven

---

## Executive Summary

The multiomics-core platform provides a reproducible, configuration-driven pipeline for multi-omics analysis built on the {targets} framework. Proteomics (DIA-NN / Limma) and RNA-seq (counts/tximport / DESeq2) are functional. This plan outlines 6 milestones to reach a production-grade, multi-omics integration platform.

---

## Current Architecture Snapshot

```
core (validation, I/O, plots, contracts)
  └─ domain (pure omics logic: proteomics | rnaseq | metabolomics)
       └─ modules (targets-ready wrappers: QC → DE → clustering → export)
            └─ pipeline (target factories consumed by _targets.R)
                 └─ Shiny payload (canonical v2.0 schema)
```

**Strengths:** Strict layering, contract validation at every boundary, config-not-code philosophy, {targets} caching.
**Gaps:** Metabolomics is a stub, no enrichment analysis, no multi-omics integration, minimal test coverage, no CI/CD, no Shiny app.

---

## Milestone 1: Foundation Hardening

**Priority:** Critical
**Estimated effort:** 5–7 weeks (1–2 developers)
**Goal:** Make the existing proteomics and RNA-seq pipelines production-reliable; add a QC-only workflow and DEP2 imputation.

### 1.1 Expand Test Suite

| Item | Detail |
|------|--------|
| **What** | Increase unit test coverage from ~5% to ≥60% across core/, domain/, and modules/ |
| **Technical approach** | Add testthat tests per layer: (a) core contracts — feed valid/invalid inputs to every `assert_*` function; (b) domain logic — test filtering, imputation, normalization, DE summarization with small synthetic datasets; (c) modules — integration tests using `tar_test()` from {targets} |
| **Files to create** | `tests/testthat/test-validation.R`, `test-alignment.R`, `test-io.R`, `test-config.R`, `test-plots.R`, `test-shiny-contract.R`, `test-proteomics-filtering.R`, `test-proteomics-imputation.R`, `test-proteomics-de.R`, `test-rna-preprocess.R`, `test-rna-de.R`, `test-clustering.R` |
| **Dependencies** | None — uses existing code |
| **Risks** | Discovering latent bugs that need fixing; time may exceed estimate if refactoring is needed |
| **Acceptance criteria** | (1) `devtools::test()` runs green; (2) ≥60% line coverage measured by `covr`; (3) every `assert_*` function has at least one positive and one negative test; (4) synthetic data fixtures live in `tests/fixtures/` |

### 1.2 CI/CD Pipeline

| Item | Detail |
|------|--------|
| **What** | GitHub Actions workflow for automated testing on every push/PR |
| **Technical approach** | `.github/workflows/check.yaml`: (1) install R + system deps; (2) `renv::restore()`; (3) `devtools::test()`; (4) `targets::tar_make()` on a small test config; (5) coverage report via `covr` + Codecov badge |
| **Dependencies** | 1.1 (test suite must exist first) |
| **Risks** | Bioconductor packages have long install times — mitigate with caching (`actions/cache` on renv library). Large renv.lock (~200 packages) may cause flaky builds |
| **Acceptance criteria** | (1) PRs show green/red check status; (2) coverage badge in README; (3) build completes in <15 minutes with cache hit |

### 1.3 Remove Hardcoded Paths & Improve Config Bootstrapping

| Item | Detail |
|------|--------|
| **What** | Replace the hardcoded Windows path in `_targets.R` (`C:/Users/sharabmi/...`) with a project-relative or environment-variable-driven config path |
| **Technical approach** | (a) Default to `config/config.yaml` relative to project root; (b) allow override via `OMICS_CONFIG` environment variable; (c) add `scripts/init_project.R` that copies a template config to the project and sets the path; (d) update onboarding docs |
| **Files to modify** | `_targets.R`, `R/core/04_config.R`, `docs/onboarding.md` |
| **Dependencies** | None |
| **Risks** | Breaking existing user workflows — mitigate by checking the old path as a fallback with a deprecation warning |
| **Acceptance criteria** | (1) `_targets.R` contains no absolute user-specific paths; (2) pipeline runs with config at default location; (3) `OMICS_CONFIG` env var override works; (4) deprecation warning if old path is used |

### 1.4 QC-Only Pipeline Mode

| Item | Detail |
|------|--------|
| **What** | A `run_mode: "qc_only"` option that executes only up to the QC pre-stage, produces a self-contained QC report, and stops — so the bioinformatician can review sample quality, clustering, and distributions before committing to DE parameters |
| **Technical approach** | **(a) Config flag.** Add a top-level `run_mode` field (`"full"` default, `"qc_only"`) to the config schema; validate in `validate_config()`. **(b) Pipeline gating.** In each `pipe_*()` factory, wrap DE / clustering / export targets in an `if (run_mode != "qc_only")` guard so {targets} simply does not register them. Alternatively, use `tar_cue(mode = "never")` — but conditional list building is cleaner and avoids stale cache confusion. **(c) QC report target.** Add a new terminal target (`*_qc_report`) that renders a parameterised Quarto/RMarkdown template (`reports/qc_report.qmd`) receiving the QC pre output object. The report should include: (1) sample count & group balance table, (2) missingness heatmap (proteomics) or library-size bar chart (RNA-seq), (3) PCA 1v2 + 1v3 colored by every configured color variable, (4) sample-distance & correlation heatmaps, (5) density/boxplot overlays, (6) filtering impact summary (features before → after), (7) a "Recommended next steps" checklist stub the bioinformatician fills in. **(d) {targets} caching.** When the user later switches to `run_mode: "full"`, all QC pre targets are already cached — DE starts from warm cache with zero rework. **(e) Per-mode support.** Works for proteomics, RNA-seq, and (once implemented) metabolomics, since all share the same `mod_*_qc_pre()` interface and return the same structure (plots, files, objects). |
| **Files to create/modify** | `R/core/04_config.R` (validate `run_mode`), `R/pipeline/proteomics/00_pipe_proteomics.R`, `R/pipeline/rnaseq/00_pipe_rnaseq.R`, `R/pipeline/metabolomics/00_pipe_metabolomics.R` (gating logic), `reports/qc_report.qmd` (new), `R/modules/*/01_mod_qc_pre.R` (ensure filtering-impact stats are returned) |
| **Config example** | `run_mode: "qc_only"  # "qc_only" or "full"` at the top level of the YAML |
| **Dependencies** | None — uses existing QC pre modules |
| **Risks** | (a) Bioinformaticians may forget to switch back to `"full"` — mitigate with a clear log message: *"QC-only mode: pipeline stopped after QC. Set run_mode: full to continue."* (b) Report template maintenance — keep it minimal and data-driven so it works across omics modes |
| **Acceptance criteria** | (1) `tar_make()` with `run_mode: "qc_only"` completes in <30 s on a typical dataset (no DE computation); (2) an HTML QC report is written to the output directory; (3) report contains all 7 sections listed above; (4) switching to `run_mode: "full"` reuses all cached QC targets; (5) works for both proteomics and RNA-seq pipelines; (6) config validation rejects unknown `run_mode` values |

### 1.5 DEP2 Imputation for Proteomics

| Item | Detail |
|------|--------|
| **What** | Add DEP2-based imputation as an alternative to the existing Perseus-style method, selectable via the existing config field `imputation.method: "dep2"` |
| **Technical approach** | **(a) New imputation function.** Create `impute_proteomics_dep2()` in `R/domain/proteomics/03_imputation.R`. This wraps `DEP::impute()` (which itself delegates to `MSnbase::impute()`), supporting the methods already anticipated in the config: MinDet, MinProb, KNN, MLE, QRILC, man, min, zero, mixed, bpca, and others provided by MSnbase. The function accepts `expr_mat` (log2, features × samples with NAs) and `cfg` (which contains `dep2_method` and `dep2_random_seed`), converts the matrix to a `SummarizedExperiment`, calls the imputation, and extracts back to a plain matrix. Returns `list(imputed, imputed_flag)` — same contract as `perseus_impute_with_flags()`. **(b) Multiple-imputation wrapper.** Create `make_imputations_dep2()` analogous to `make_imputations_proteomics()`. For stochastic methods (MinDet, MinProb, QRILC, bpca), run N repetitions with incremented seeds, exactly like Perseus. For deterministic methods (KNN, MLE, SVD), run once and replicate (with a log warning that multiple imputation is not meaningful for deterministic methods). **(c) Dispatch in preprocessing.** Modify `preprocess_proteomics()` and `mod_proteomics_de()` to dispatch on `cfg$imputation$method`: `"perseus"` → current path, `"dep2"` → new path. The downstream contract (matrix dimensions, rownames, colnames, no NAs) is identical, so DE and everything downstream is unchanged. **(d) Config integration.** The config template already has the fields (`dep2_method: "MinDet"`, `dep2_random_seed: 1`). Add validation: if `method == "dep2"`, require `dep2_method` to be one of the supported values; warn if `dep2_method` is deterministic and `no_repetitions > 1`. **(e) QC comparison.** In the QC pre module, if imputation method is DEP2, generate the same imputation histogram and boxplot using the existing `qc_imputation_summary()` (it only needs the imputed matrix and flag — method-agnostic). |
| **New dependency** | `DEP` (Bioconductor — already depends on `SummarizedExperiment` and `limma` which are in renv). Transitively pulls in `MSnbase` for the actual imputation engines |
| **Files to modify** | `R/domain/proteomics/03_imputation.R` (add `impute_proteomics_dep2()`, `make_imputations_dep2()`), `R/domain/proteomics/04_preprocess.R` (dispatch), `R/modules/proteomics/02_mod_de.R` (dispatch), `R/domain/proteomics/90_config_validate.R` (validate dep2 fields), `config/templates/proteins_config.yaml` (document dep2 options) |
| **Dependencies** | None — proteomics pipeline is already stable |
| **Risks** | (a) `DEP`/`MSnbase` add ~15 transitive dependencies — test that `renv::restore()` still works cleanly; (b) some MSnbase imputation methods require complete columns or have minimum-sample-count requirements — add pre-flight checks; (c) KNN imputation can be slow for >10k proteins — document performance expectations |
| **Acceptance criteria** | (1) `imputation.method: "dep2"` with `dep2_method: "MinDet"` produces a fully imputed matrix passing `assert_numeric_matrix()`; (2) multiple imputations with stochastic methods produce N distinct matrices; (3) deterministic methods produce 1 matrix with a warning; (4) DE results downstream are structurally identical to Perseus path; (5) imputation QC plots generate correctly; (6) existing Perseus path is completely unaffected (no regression); (7) unit tests cover at least MinDet, KNN, and MLE methods with synthetic data |

### 1.6 Standardize Error Messages & Logging

| Item | Detail |
|------|--------|
| **What** | Replace ad-hoc `message()`/`warning()`/`stop()` calls with a lightweight structured logging approach |
| **Technical approach** | Add `R/core/00a_logging.R` with `log_info()`, `log_warn()`, `log_error()` wrappers that include timestamp and calling function. Use `cli` package for colored terminal output. Ensure all domain functions use these instead of raw `message()` |
| **Dependencies** | None. `cli` is already an indirect dependency |
| **Risks** | Low — purely additive |
| **Acceptance criteria** | (1) All pipeline stages produce timestamped log output; (2) errors include the originating function name; (3) log level can be set via config (`verbose: true/false`) |

---

## Milestone 2: Metabolomics Pipeline

**Priority:** High
**Estimated effort:** 4–5 weeks (1–2 developers)
**Goal:** Bring metabolomics from stub to feature-parity with proteomics.

### 2.1 Metabolomics Input & Normalization

| Item | Detail |
|------|--------|
| **What** | Implement metabolomics data loading, validation, and normalization |
| **Technical approach** | (a) Support common input formats: peak intensity tables (CSV/TSV) from MetaboAnalyst, Compound Discoverer, MS-DIAL, mzMine; (b) implement normalization methods in `R/domain/metabolomics/01_normalization.R`: log2 transform, median normalization, PQN (probabilistic quotient normalization), quantile normalization; (c) add metabolomics-specific validation in `90_config_validate.R` |
| **Files to modify** | `R/domain/metabolomics/00_inputs.R`, `01_normalization.R`, `02_preprocess.R`, `90_config_validate.R` |
| **Dependencies** | None |
| **Risks** | Diverse vendor formats may require iterative format support. Start with a generic intensity matrix and add vendor-specific parsers as needed |
| **Acceptance criteria** | (1) Can load peak-intensity CSV with metabolite IDs as rows; (2) at least 3 normalization methods work and produce `assert_numeric_matrix()`-passing output; (3) config validation catches missing/malformed metabolomics config; (4) unit tests for each normalization method |

### 2.2 Metabolomics DE & QC

| Item | Detail |
|------|--------|
| **What** | Differential abundance analysis and QC for metabolomics |
| **Technical approach** | (a) Limma-based DE (same as proteomics, metabolites behave similarly to proteins post-normalization); (b) optional: t-test / Wilcoxon fallback for small sample sizes; (c) QC module: PCA, sample-distance heatmap, density overlay — reuse existing `R/core/06_plots.R` and `08_qc.R`; (d) volcano and MA plots via existing core functions |
| **Files to create/modify** | `R/domain/metabolomics/03_de.R`, `R/modules/metabolomics/02_mod_de.R`, `03_mod_qc_post.R` |
| **Dependencies** | 2.1 |
| **Risks** | Metabolomics data often has high missingness — may need imputation similar to proteomics. Consider reusing `impute_proteomics_perseus()` or generalizing it |
| **Acceptance criteria** | (1) DE results pass `assert_de_contract()`; (2) QC plots generate for metabolomics data; (3) Shiny payload passes `assert_shiny_payload_contract()` |

### 2.3 Metabolomics Pipeline & Shiny Export

| Item | Detail |
|------|--------|
| **What** | Wire metabolomics into the {targets} pipeline and produce Shiny-compatible output |
| **Technical approach** | (a) Complete `R/pipeline/00_pipe_metabolomics.R` following the proteomics template; (b) implement `R/domain/metabolomics/07_shiny_export.R` to build canonical v2.0 payload; (c) add clustering support (reuse binary patterns from `09_clustering.R`); (d) uncomment metabolomics in `_targets.R` |
| **Dependencies** | 2.1, 2.2 |
| **Risks** | Config template (`metabolomics_template.yaml`) needs to be fleshed out; may discover missing contract keys |
| **Acceptance criteria** | (1) `tar_make()` runs metabolomics end-to-end on test data; (2) Shiny payload passes full contract validation; (3) Excel export works via existing `05_export_excel.R`; (4) metabolomics config template is documented |

---

## Milestone 3: Enrichment & Pathway Analysis

**Priority:** High
**Estimated effort:** 3–4 weeks (1 developer)
**Goal:** Implement gene set / pathway enrichment for all omics types.

### 3.1 Core Enrichment Engine

| Item | Detail |
|------|--------|
| **What** | Replace the 1-line stub in `R/core/10_enrichment.R` with a functional enrichment module |
| **Technical approach** | (a) ORA (over-representation analysis) using `clusterProfiler::enrichGO()` / `enrichKEGG()` for gene-level data; (b) GSEA (gene set enrichment analysis) using `fgsea` for ranked gene lists; (c) metabolite set enrichment via `MetaboAnalystR` or a simple hypergeometric test against KEGG compound sets; (d) keep the engine pure: input = feature list + universe + gene sets → output = enrichment table |
| **New dependencies** | `clusterProfiler`, `fgsea`, `org.Hs.eg.db` (and other organism DBs as needed), `msigdbr` for MSigDB gene sets |
| **Files to create/modify** | `R/core/10_enrichment.R` (core logic), `R/domain/rnaseq/05_enrichment.R`, `R/domain/proteomics/05_enrichment.R`, `R/domain/metabolomics/04_enrichment.R`, corresponding module files |
| **Dependencies** | Milestone 2 for metabolomics; RNA-seq and proteomics can start immediately |
| **Risks** | (a) Organism-specific annotation databases are large — need clear config for organism selection; (b) ID mapping (Ensembl ↔ Entrez ↔ UniProt) is error-prone — build a robust `map_feature_ids()` utility; (c) additional Bioconductor dependencies increase build time |
| **Acceptance criteria** | (1) ORA returns a data.frame with term, p-value, adjusted p-value, gene list, gene ratio; (2) GSEA returns ranked enrichment with NES and leading edge; (3) results are added to Shiny payload under new canonical keys (`enrichment_ora`, `enrichment_gsea`); (4) at least GO and KEGG are supported; (5) unit tests with synthetic gene lists |

### 3.2 Enrichment Visualization

| Item | Detail |
|------|--------|
| **What** | Dot plots, bar plots, enrichment maps, and network plots for enrichment results |
| **Technical approach** | (a) Use `clusterProfiler`'s built-in plotting where possible; (b) for custom needs, add `plot_enrichment_dotplot()` and `plot_enrichment_barplot()` to `R/core/06_plots.R`; (c) optional: `enrichplot::cnetplot()` for gene-concept networks |
| **Dependencies** | 3.1 |
| **Risks** | Low — visualization is additive |
| **Acceptance criteria** | (1) At least dot plot and bar plot generated per contrast; (2) plots saved to output directory; (3) plot objects included in Shiny payload |

---

## Milestone 4: Multi-Omics Integration

**Priority:** Medium-High
**Estimated effort:** 5–7 weeks (1–2 developers)
**Goal:** Enable cross-omics analysis for samples profiled by multiple platforms.

### 4.1 Multi-Omics Data Alignment Layer

| Item | Detail |
|------|--------|
| **What** | A new module that aligns samples across omics layers and produces a unified multi-omics object |
| **Technical approach** | (a) Create `R/core/11_integration.R` with `align_omics_layers()`: takes a list of omics objects (each with expression matrix + metadata), intersects samples, and returns aligned matrices; (b) support both inner join (only shared samples) and full join (with NA padding); (c) validate sample IDs across layers using existing `assert_*` functions |
| **Files to create** | `R/core/11_integration.R`, `R/domain/integration/00_alignment.R` |
| **Dependencies** | Milestones 1–2 (all single-omics pipelines must be stable) |
| **Risks** | Sample ID mismatches between platforms (different naming conventions) — mitigate with a sample mapping table in config |
| **Acceptance criteria** | (1) Given 2+ omics layers, returns aligned matrices with shared samples; (2) warns about dropped samples; (3) handles edge case of zero overlap gracefully |

### 4.2 MOFA / Factor Analysis Integration

| Item | Detail |
|------|--------|
| **What** | Multi-Omics Factor Analysis (MOFA+) integration for unsupervised multi-omics dimensionality reduction |
| **Technical approach** | (a) Add `R/domain/integration/01_mofa.R`: prepares input for `MOFA2::create_mofa_from_matrix()`, runs training, extracts factors and weights; (b) visualization: factor correlation heatmap, variance explained per omics layer, top features per factor; (c) config: `integration.mofa.n_factors`, `integration.mofa.convergence_mode` |
| **New dependencies** | `MOFA2` (Bioconductor) |
| **Dependencies** | 4.1 |
| **Risks** | MOFA2 requires Python backend (mofapy2) — adds complexity. Mitigate by documenting setup clearly and providing a Docker option (Milestone 6) |
| **Acceptance criteria** | (1) MOFA model trains successfully on ≥2 omics layers; (2) variance decomposition plot generated; (3) factor scores added to Shiny payload; (4) top features per factor exported |

### 4.3 Correlation & Co-expression Networks

| Item | Detail |
|------|--------|
| **What** | Cross-omics correlation analysis and weighted co-expression networks |
| **Technical approach** | (a) Pairwise feature correlation across omics (e.g., protein–mRNA Pearson/Spearman correlation for matched gene products); (b) optionally: WGCNA-style module detection within and across omics; (c) network visualization with `igraph` for high-correlation subnetworks |
| **New dependencies** | `WGCNA` (optional), `igraph` |
| **Dependencies** | 4.1 |
| **Risks** | Computational cost: all-vs-all correlation is O(n*m) and can be large. Mitigate by restricting to DE features or top-variable features |
| **Acceptance criteria** | (1) Correlation matrix produced for matched features across 2 omics; (2) scatter plot of protein vs mRNA log2FC per gene; (3) results exportable to Shiny payload |

---

## Milestone 5: Shiny Application

**Priority:** Medium
**Estimated effort:** 5–6 weeks (1–2 developers)
**Goal:** Build an interactive Shiny dashboard that consumes the canonical payload.

### 5.1 Core Shiny App Skeleton

| Item | Detail |
|------|--------|
| **What** | A modular Shiny app with tab-based navigation, consuming the v2.0 payload |
| **Technical approach** | (a) Use `{golem}` or `{rhino}` framework for production Shiny; (b) tabs: Overview, QC, Differential Expression, Clustering, Enrichment, Integration; (c) load payload RDS on startup; (d) use `{bslib}` for modern Bootstrap 5 UI |
| **New dependencies** | `shiny`, `bslib`, `DT`, `plotly`, `golem` (or `rhino`) |
| **Files to create** | New `app/` directory (or separate repository — team decision) |
| **Dependencies** | All prior milestones (payload must be complete) |
| **Risks** | Scope creep — Shiny apps grow quickly. Mitigate by strictly tying features to payload keys (if it's not in the payload, it's not in the app) |
| **Acceptance criteria** | (1) App loads any valid v2.0 payload; (2) QC tab shows PCA, heatmaps; (3) DE tab shows interactive volcano/MA plots with gene hover; (4) app passes `shinytest2` smoke tests |

### 5.2 Interactive DE Explorer

| Item | Detail |
|------|--------|
| **What** | Interactive volcano/MA plots with click-to-select, gene search, and dynamic filtering |
| **Technical approach** | (a) `plotly` for interactive volcano with hover info (gene name, FC, p-value); (b) click-to-highlight → update linked expression boxplot; (c) DT table with server-side filtering (search by gene, filter by FC/p-value thresholds); (d) download selected gene lists as CSV |
| **Dependencies** | 5.1 |
| **Risks** | Performance with >20k features in plotly — mitigate with WebGL rendering (`toWebGL()`) and server-side data tables |
| **Acceptance criteria** | (1) Volcano plot renders in <2s for 20k features; (2) clicking a point shows expression across groups; (3) table filters update in real-time; (4) CSV download works |

### 5.3 Dynamic Report Generation

| Item | Detail |
|------|--------|
| **What** | Generate a downloadable HTML/PDF report from the Shiny app or pipeline |
| **Technical approach** | (a) Quarto (`.qmd`) or RMarkdown template parameterized by payload path; (b) sections: experimental design, QC summary, DE results (top genes table, volcano), clustering, enrichment; (c) callable from Shiny via `rmarkdown::render()` or from the pipeline as a final target |
| **Files to create** | `reports/analysis_report.qmd` |
| **Dependencies** | 5.1 (for Shiny integration) or standalone |
| **Risks** | Rendering large reports with many plots is slow — mitigate by pre-rendering plot PNGs in the pipeline |
| **Acceptance criteria** | (1) Report renders successfully from payload; (2) includes at least: sample table, PCA, top-20 DE genes per contrast, volcano plots; (3) takes <60s for a typical dataset |

---

## Milestone 6: Deployment & Reproducibility

**Priority:** Medium
**Estimated effort:** 2–3 weeks (1 developer)
**Goal:** Containerize the pipeline and make it deployable anywhere.

### 6.1 Docker Image

| Item | Detail |
|------|--------|
| **What** | A Docker image with all R, Bioconductor, and system dependencies pre-installed |
| **Technical approach** | (a) Base image: `bioconductor/bioconductor_docker:3.22`; (b) `COPY renv.lock` + `renv::restore()`; (c) entrypoint: `Rscript scripts/run_pipeline.R`; (d) mount data and config as volumes; (e) multi-stage build to reduce image size |
| **Files to create** | `Dockerfile`, `docker-compose.yml`, `.dockerignore` |
| **Dependencies** | 1.2 (CI/CD should exist first) |
| **Risks** | Large image size (~5–8 GB with Bioconductor). Mitigate with layer caching and `.dockerignore` |
| **Acceptance criteria** | (1) `docker build` succeeds; (2) `docker run` with mounted test data completes pipeline; (3) CI builds and pushes image to registry; (4) documented in README |

### 6.2 Additional Proteomics Engine Support

| Item | Detail |
|------|--------|
| **What** | Extend proteomics input beyond DIA-NN to support MaxQuant, FragPipe, Spectronaut |
| **Technical approach** | (a) Add engine-specific parsers in `R/domain/proteomics/01_expression.R`; (b) config field `engine: DIANN | MAXQUANT | FRAGPIPE | SPECTRONAUT`; (c) each parser normalizes to the same internal matrix format; (d) use a factory/dispatch pattern (already in place for DIA-NN) |
| **Dependencies** | None — can start anytime |
| **Risks** | Vendor format changes between software versions — pin tested versions in docs. Need test data for each engine |
| **Acceptance criteria** | (1) Each engine loads data into `assert_numeric_matrix()`-passing format; (2) downstream pipeline produces identical structure regardless of engine; (3) at least one test dataset per engine in `tests/fixtures/` |

### 6.3 Data Versioning & Provenance

| Item | Detail |
|------|--------|
| **What** | Track input data versions and link results to specific data snapshots |
| **Technical approach** | (a) Compute SHA-256 checksums of all input files at pipeline start; (b) store in `execution_info` alongside config snapshot; (c) optionally integrate with DVC (Data Version Control) for large file tracking; (d) add `data_checksums` to Shiny payload for traceability |
| **Dependencies** | 1.3 (config improvements) |
| **Risks** | DVC adds tool complexity. Start with checksums-only, add DVC later if needed |
| **Acceptance criteria** | (1) Every pipeline run records input file checksums; (2) checksums are included in execution metadata; (3) re-running with different data produces different checksums |

---

## Milestone Summary & Timeline

```
          Month 1–2          Month 3–4          Month 5–6         Month 7–8
         ┌─────────────────────────────┐
    M1   │ Foundation Hardening        │
         │ Tests, CI, QC Mode, DEP2    │
         └─────────────────────────────┘
         ┌──────────┐      ┌──────────┐
    M2   │Metabolom. │──────│──────────│
         │ Pipeline  │      │          │      ┌──────────┐      ┌──────────┐
         └──────────┘      └──────────┘      │          │      │          │
                           ┌──────────┐      │          │      │          │
    M3   ·················>│Enrichment│──────│          │      │          │
                           │& Pathways│      │          │      │          │
                           └──────────┘      │          │      │          │
                                             ┌──────────┐      │          │
    M4   ···································>│Multi-Omics│─────│          │
                                             │Integration│     │          │
                                             └──────────┘      │          │
                                             ┌──────────┐      ┌──────────┐
    M5   ···································>│  Shiny    │─────│  Shiny   │
                                             │  App v1   │     │  Polish  │
                                             └──────────┘      └──────────┘
                                                               ┌──────────┐
    M6   ·····················································>│ Deploy & │
                                                               │ Docker   │
                                                               └──────────┘
```

---

## Cross-Cutting Concerns

### ID Mapping Strategy
Many features (enrichment, integration, annotation) require mapping between identifier systems (Ensembl gene IDs, Entrez IDs, UniProt accessions, KEGG compound IDs). Build a central `R/core/12_id_mapping.R` utility early (as part of Milestone 3) that:
- Uses `AnnotationDbi` + organism-specific OrgDb packages
- Caches mappings per session
- Logs unmapped IDs with warnings

### Batch Effect Correction
Consider adding `limma::removeBatchEffect()` or `sva::ComBat()` as an optional preprocessing step (configurable via YAML). This is relevant for all omics types when samples come from multiple batches. Add to `R/core/` as a generic utility, callable from any domain's preprocessing.

### Config Schema Validation
As configuration grows more complex (integration, enrichment, Shiny), consider adopting a formal schema (JSON Schema or a custom R-based validator). This catches config errors at load time rather than mid-pipeline.

---

## Risk Register

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Bioconductor dependency conflicts | Medium | High | Pin versions via renv; test upgrades in a separate branch |
| MOFA2 Python dependency complexity | Medium | Medium | Docker container; document conda/pip setup; make MOFA2 optional |
| Shiny app scope creep | High | Medium | Strict rule: only visualize what the payload contains |
| Test data availability (multi-engine, multi-omics) | Medium | Medium | Create synthetic fixtures; partner with lab for anonymized real data |
| Single developer bus factor | Medium | High | Document everything; pair-program on critical modules; enforce PR reviews |
| DEP/MSnbase transitive dependency bloat | Medium | Medium | Pin versions via renv; test `renv::restore()` after adding DEP; fallback to Perseus if install fails |
| QC-only mode cache invalidation confusion | Low | Low | Clear log messages when mode switches; document {targets} caching behavior in onboarding |
| Large dataset performance (>50k features) | Low | Medium | Profile with large data early; use sparse matrices where applicable |

---

## Definition of Done (Global)

Every feature delivered under this plan must meet:

1. **Code:** Follows existing architecture (core → domain → modules → pipeline)
2. **Contracts:** All inputs/outputs validated by `assert_*` functions
3. **Config:** New parameters documented in YAML templates with comments
4. **Tests:** Unit tests with ≥80% coverage for new code; integration test via `tar_test()`
5. **Docs:** Updated developer guide if architecture changes; updated onboarding if user-facing
6. **Review:** PR reviewed by at least one other team member
7. **CI:** All tests pass in CI before merge
