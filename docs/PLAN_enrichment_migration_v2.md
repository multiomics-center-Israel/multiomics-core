# Current Status

Status: **Phase 1 COMPLETE and fully validated** (local & uncommitted). All Phase 1 validation tasks passed — see §4.4 for today's results. Implementation smoke-validated end-to-end and green; two latent downstream-consumer bugs were found, fixed, and verified: (1) `rna_pipeline_summary` integer/double (§4.2) and (2) `extract_enrichment_df()` `rbind` schema mismatch on heterogeneous ORA+GSEA tables (§4.3). Both multiomics (`clustering_res = NULL`) behavior and online-regression routing (`annotation_dir: ""`) are validated. Both fixes are in downstream consumers; the enrichment producer and its Phase 1 contract are unchanged. Parallel enrichment execution (§12) is implemented. **Phase 1 is committed locally** as `aaf8d8b` (not pushed/merged), followed by an **isolated RNG-reproducibility commit** (§12.1: `params$seed` → `future.seed`, making GSEA identical across worker counts and independent rebuilds; one-time GSEA-baseline shift, ORA unchanged). Phase 2 (legacy plots) is next, resuming from this fully reproducible baseline.

Current Branch:
feature/enrichment-migration-v2

Current Approved Phase:
Phase 1 – RNAseq Local Enrichment Core — **implementation + validation COMPLETE** (uncommitted)

Next Planned Phase:
Parallel enrichment execution (§12) — **IMPLEMENTED & validated** (workers=1 ≡ workers=4, identical results). **Next: Phase 2** – Legacy enrichment plots (gated; §13 preparation)

Implementation Started:
Yes — Phase 1 complete in the working tree, NOT committed/pushed. A reduced "smoke" config drove a full RNA single-omics run end-to-end (green). Two latent downstream bugs were fixed (`R/domain/rnaseq/10_pipeline_summary.R` §4.2; `R/domain/multiomics/07_enrichment.R` §4.3). Today's validation (§4.4) closed every remaining Phase 1 item: unit tests, `tar_validate()`, full smoke run, offline ORA + GSEA, exec summary / Shiny payload / HTML report, both downstream fixes verified, multiomics `clustering_res = NULL` behavior (GSEA + DE-derived ORA; partition & binary ORA omitted), and online-regression routing (`annotation_dir: ""` → legacy online path, local path bypassed, expected non-model-organism behavior). Remaining before any code tomorrow: the commit decision only.

Merged to Main:
No

---

This document is the **authoritative, self-contained source of truth** for the enrichment migration project. Someone resuming this work months from now should be able to continue from this document alone, without the originating chat history. Before making any implementation decision, read this document first.

> _Maintenance note: keep the **Current Status** block above current — update `Status`, `Current Approved Phase`, `Next Planned Phase`, `Implementation Started`, and `Merged to Main` whenever a phase is approved, implementation begins, or anything is committed/merged._

**Branch discipline:** work only on `feature/enrichment-migration-v2`. **No merge to main, no PR, no rebase/cherry-pick, no push** until explicitly approved.
**Reference branch (read-only):** `claude/enrichment-migration-plan-cuduc` — never merged; consulted only to re-apply/validate code.
**Companion docs (same branch):** `docs/SPEC_enrichment_migration.md`, `docs/HANDOVER_enrichment_migration.md`, `docs/audit-enrichment-layer.md`.
**Legacy source of truth:** `Projects/Uri_Gat/` — `Neat_RNA-Seq_1.0.Rmd` (section "Enrichment - preparations" onward), `Functions_for_clustering_and_enrichment_1.0.R`, `Functions_for_Neat_RNA-Seq_1.0.R`, `Global_parameters_wo_samples_02_07_09.R`.

---

# 1. Project Goal

**Why this migration exists.** The current pipeline's RNA-seq enrichment relies on *online/dynamic* gene-set loading (OrgDb, KEGG REST, biomaRt). The team's established, trusted enrichment workflow (the legacy "Neat RNA-Seq" Rmd) is **offline and table-driven**. This project replaces the online RNA-seq enrichment with a faithful reproduction of that legacy workflow, while fitting cleanly into the current `{targets}` architecture.

**What "behavioral compatibility with the legacy workflow" means.** Reproduce the legacy *analysis intent and results*, not its line-by-line code: same enrichment inputs (ORA collections), same statistical parameters, same ranking methods, same offline data source, same output organization. Modernization of *implementation* is acceptable only when the *behavior/results* are preserved.

**Which architecture it integrates into.** The current repo uses a 5-layer structure (`core → services → domain → modules → pipeline`) and a `{targets}` DAG. RNA-seq enrichment is one `tar_target` (`rna_pathway_res`) produced by `mod_rnaseq_pathway()` (module layer), which calls shared functions in `R/core/09_enrichment.R`. The pipeline factory `pipe_rnaseq(skip_outputs=)` splits **core** targets (always run; multiomics depends on them) from **single-omics output** targets (skipped when multiomics is active).

**Why incrementally.** The legacy workflow is large (offline loading, ranked lists, GSEA, cluster-based ORA, plots, per-pathway artifacts) and touches owner-gated areas (the Shiny contract, the HTML report). Phasing keeps each change reviewable, lets behavior be verified against legacy before expanding scope, and avoids destabilizing the live online path. Each phase is separately approved.

---

# 2. Migration Strategy

**Guiding principles**
- **Behavior first.** The legacy workflow defines correct behavior; code may be modernized only if results are unchanged.
- **Offline-first.** Runtime enrichment uses only local tab files. No KEGG REST, no biomaRt, no OrgDb-driven pathway loading, no GMT as primary source, no runtime downloads.
- **Backward compatibility.** The existing online pathway path must keep working unchanged for current users; the new offline path is opt-in.
- **Fail soft.** Enrichment must never crash the pipeline — missing data / no overlap degrade to warnings + valid (possibly empty) results.
- **Phased + gated.** Implement, verify against legacy, then expand. No commits/pushes/merges without explicit approval.

**Architecture decisions (high level — details in §5)**
- One stable target `rna_pathway_res`, with its clustering input wired conditionally on `skip_outputs` (single target, no duplicates).
- A single config switch (`enrichment.annotation_dir`) selects offline vs the existing online path.
- New functions are appended to the already-de-duplicated `R/core/09_enrichment.R`; the module is *extended surgically* (online path preserved); the pipeline is *rewired minimally*.
- Sequential execution only (no `future.apply`); no `renv.lock` change.
- Preserve the legacy output directory layout.

**Compatibility goals:** same ORA collections, same ORA/GSEA parameters, same ranking methods, same fold-enrichment math, same offline data contract, downstream consumers (`collect_pipeline_stats`, `extract_enrichment_df`, `get_pathway_highlights`) keep working. **Contract clarification (A12):** `pathway_results` may legitimately contain both ORA and GSEA tables with *different* column schemas; consumers must combine them by aligning columns and NA-filling (not assume a single schema). Two consumers were hardened for this during validation — `collect_pipeline_stats` (§4.2) and `extract_enrichment_df` (§4.3).

---

# 3. Legacy Behavioral Knowledge (what the legacy workflow actually does)

Audited from `Projects/Uri_Gat/` (Rmd "Enrichment - preparations" onward + the two function files + global params).

**Offline annotation.** Annotation is read only from local tab files in `functional_annot_dir`:
`KEGG_pathway2gene.tab`, `KEGG_pathway2name.tab`, and `GO2gene_{BP,MF,CC}.tab` / `GO2name_{BP,MF,CC}.tab`. `library(KEGGREST)` is loaded but **never called** in the enrichment body (only a KEGG hyperlink string is built for display); no biomaRt; no runtime downloads.

**TERM2GENE / TERM2NAME usage.** Both ORA (`clusterProfiler::enricher`) and GSEA (`clusterProfiler::GSEA`) are driven by the local two-column tables: `Pathway2gene`/`GO2gene` → TERM2GENE, `Pathway2name`/`GO2name` → TERM2NAME.

**GO Expansion is preprocessing, not runtime.** Legacy default `functional_annot_dir = "Func_annot_data_expanded"` — GO assignments are expanded ahead of time (e.g. `clusterProfiler::buildGOmap`) during annotation preparation. Runtime enrichment reads the **already-expanded** tables and performs no expansion itself. "With GO Expansion" = supply pre-expanded GO tables.

**GO databases.** BP, MF, CC — all three, explicit loop.

**ORA is cluster-based**, run over many gene-list *collections* (`gene_lists[[method]][[round]]`), each via `Clusters_Enrichment_Test()`:
- `all_DE` (union of DE genes, single `"all"` cluster), `contrasts` (up/down as separate clusters), `contrasts_wo_direction` (all DE per contrast), `partition_clustering`, `binary_patterns`, **and** `manual_clustering` (×2) + `unified_clusters`.
- ORA params: `enricher(minGSSize=0, maxGSSize=10000, qvalueCutoff=1, pvalueCutoff=0.05, pAdjustMethod="fdr")`; per-cluster results merged via `merge_result()` (Cluster column).
- Fold enrichment: `signif((in_cluster_in_term/in_cluster)/(in_term/in_genome), 2)`, with GeneRatio/BgRatio expanded into numeric columns.

**GO simplify behavior.** For GO ORA the legacy **always** runs `simplify(allRes, cutoff=0.7, by="p.adjust", select_fun=min)` (after the regular result), producing a separate "Simplify_*" output. It passes **no explicit `semData`/OrgDb** — it relied on clusterProfiler/GOSemSim defaults (i.e. on a locally available OrgDb / GO.db in the legacy environment). KEGG ORA is not simplified.

**GSEA.** Run **over contrasts only** (iterates ranked gene lists, not clusters), for KEGG and GO BP/MF/CC. Params: `GSEA(minGSSize=4, maxGSSize=length(unique(TERM2GENE genes)), pvalueCutoff=0.05, pAdjustMethod="fdr")`.

**GSEA ranking methods.** All of these are computed (none is singled out as "the" default):
- `pval_wo_direction` = `-log10(pvalue)`;
- `pval_with_direction` = `sign(linearFC) * -log10(pvalue)`;
- `fc` = `log2(ifelse(linearFC>0, linearFC, -1/linearFC))` then `signif(digits=4)` (operates on the signed **linear** FC);
- plus `any_contrast` (cross-contrast row-wise `min(pvalue)`), stored under `pval_wo_direction`.

**Per-pathway GSEA artifacts** (gseaplot2 PNG + leading-edge / core-gene CSV) were produced **only on demand** in a manual "Explore a specific pathway" section — not automatically for every significant pathway.

**Defaults (Global params):** `ENRICHMENT_PVAL_CUTOFF=0.05`, `ENRICHMENT_PADJ_METHOD='fdr'`, `MAX_TERMS_IN_DOTPLOT=20`, `LINEAR_FC_CUTOFF=1.5`, `PADJ_CUTOFF=0.05`, `FILTER_KEGG_PATHWAYS_BY_TAXON=NA`.

**Key assumptions discovered:** the legacy had no `Gene:`-prefix handling (it was internally consistent); GO simplify depended on an OrgDb being present in the environment though the code never passed one explicitly; manual clustering is off by default (`PERFORM_MANUAL_CLUSTERING=F`).

---

# 4. Implementation Progress (Phase 1, in working tree)

What functionality now exists (all local/uncommitted on `feature/enrichment-migration-v2`):

**Offline data loading** — `load_local_pathway_tables()` reads KEGG + GO_{BP,MF,CC} tab files into TERM2GENE/TERM2NAME pairs; header auto-detected; **missing DB files skipped with a warning** (others still load); **<5% gene-overlap warning**; returns a valid empty structure (never crashes).

**Ranked gene lists** — the four legacy rankers (`rank_by_pval_wo_direction`, `rank_by_pval_with_direction`, `rank_by_fc`, `rank_by_min_pval_any_contrast`) and `build_ranked_gene_lists()`. `fc` recovers linear FC from the pipeline's `log2FoldChange` (`ifelse(lfc>=0, 2^lfc, -(2^-lfc))`) then applies the legacy transform + `signif(4)`; `any_contrast` reconstructed from per-contrast tables (row-wise min p) — both verified equivalent to legacy.

**GSEA** — `run_gsea_local()` (`minGSSize=4`, `maxGSSize=length(unique(TERM2GENE[,2]))`) and `run_gsea_all()` over ranking × contrast × DB, **sequential** (parallel branch present but gated off; `future.apply` not used). Writes `Enrichment/GSEA/<db>/ranking_by_<method>/<contrast>/GSEA_results_<contrast>.csv` and dotplots.

**Cluster-based ORA** — `run_cluster_ora()` (per-cluster `enricher(minGSSize=0, maxGSSize=10000, qvalueCutoff=1)` → `merge_result()` → `process_enrichment_table()` fold-enrichment), GO simplify optional/gated (see §5). `build_gene_lists()` assembles the ORA collections: **`all_DE`**, `contrasts`, `contrasts_wo_direction`, `partition`, `binary_patterns`. Writes `Enrichment/ORA/<db>/...` CSVs (+ `Simplify_*` when enabled).

**Module routing** — `mod_rnaseq_pathway(de_res, pre, config, out_dir, clustering_res = NULL)`: if `enrichment.annotation_dir` is set → `.run_local_enrichment()` (offline); otherwise the **existing online path runs unchanged**. `.store_ora_result()` adds `padj`/`pathway` columns for downstream compatibility. Return shape preserved: `list(annotation, pathway_results, plot_files)` with `pathway_results` always a (possibly empty) named list of data.frames carrying `padj`.

**Pipeline wiring** — single `rna_pathway_res`, `pre = rna_pre`; clustering wired conditionally by the `pipe_rnaseq()` factory: `skip_outputs=TRUE` (multiomics) → `clustering_res = NULL` (GSEA-only + ORA-skip warning); `skip_outputs=FALSE` (single-omics) → `clustering_res = rna_clustering_obj` (GSEA + ORA). Verified by deparsing both target commands.

**Config** — `modes.rna.enrichment:` block added to `config/templates/rna_config.yaml` (see §7).

**Tests** — `tests/testthat/test-enrichment-local.R` + synthetic fixtures under `tests/testthat/fixtures/enrichment_local/` (no real data): loader (incl. missing-DB + low overlap), the four rankers (exact values), `build_ranked_gene_lists`, `process_enrichment_table`, `build_gene_lists` (incl. `all_DE`, partition disambiguation, prefix stripping, binary data.frame), and clusterProfiler-gated ORA/GSEA guards.

**Untouched (intentionally):** `R/domain/rnaseq/07_pathway.R` (only `build_pathway_volcano_data()`), the online path, `renv.lock`, the HTML report, Shiny.

**Files touched (Phase 1):** `R/core/09_enrichment.R` (append + edits), `R/modules/rnaseq/05_mod_pathway.R` (extend), `R/pipeline/rnaseq/00_pipe_rnaseq.R` (rewire), `config/templates/rna_config.yaml` (add), `tests/testthat/test-enrichment-local.R` + fixtures (new), **`R/domain/rnaseq/10_pipeline_summary.R`** (one-line downstream-compat fix — see §4.2), and **`R/domain/multiomics/07_enrichment.R`** (one-line downstream-compat fix — see §4.3).

### 4.1 Smoke validation completed (this session)
A reduced **smoke configuration** (`Projects/Uri_Gat/config_smoke.yaml`, `analysis_round: Smoke_01`) drove the **full RNA single-omics pipeline end-to-end**. Validated as PASS:
- Unit tests (`testthat`, filter `enrichment-local`) and `tar_manifest()` / `tar_validate()`.
- `config_smoke.yaml` loaded; **offline local enrichment activated** (GO_CC tables loaded); **ORA** and **GSEA** completed.
- RNA legacy outputs (`Datasets/`, `Final_results_*.xlsx`), **executive summary**, **Shiny payload**, and the **rendered HTML report** all generated.
- Output inspection (`Projects/Uri_Gat/outputs/.../rna/Enrichment/`): all 5 ORA collections (`all_DE`, `contrasts`, `contrasts_wo_direction`, `partition`, `binary_patterns`) + GSEA per method×contrast (+`any_contrast`); correct columns; no `Simplify_*`, no per-pathway artifacts, no KEGG, no empty files — i.e. the M1/M2 fixes and config-gated behaviors are confirmed at runtime.

The **only** failure was the final target `rna_pipeline_summary` (see §4.2). **After the §4.2 fix the smoke pipeline was rebuilt and is now green** — verified read-only against the smoke `_targets` store: `tar_meta()` reports no errored targets, `rna_pipeline_summary` completed (0.46 s, output `Results_U_Gat_Smoke_01/rna_pipeline_summary.html`), and `tar_outdated()` returns `character(0)` (the entire DAG, incl. `rna_exec_summary` / `rna_shiny_payload` / `rna_report`, is in sync with the current source).

### 4.2 Downstream-compatibility fix (`10_pipeline_summary.R`)
`rna_pipeline_summary` failed with `values must be type 'integer', but FUN(X[[1]]) result is type 'double'`. Root cause (full analysis in §10): in `collect_pipeline_stats()` the helper `count_from_df()` returned `up`/`down` as **double `0`** for tables lacking a `NES` column, while the aggregation uses `vapply(..., FUN.VALUE = 0L, ...)` (integer). The **new ORA tables intentionally have no `NES`** (a GSEA-only field), so they were the first NES-less input this function ever received. Fix: change the two fallbacks from `else 0` to **`else 0L`** (lines ~134–135). Minimal, logic-preserving, one file; does not touch enrichment objects, the Phase 1 contract, ORA/GSEA outputs, or reports (other than letting the pipeline complete).

### 4.3 Downstream-compatibility fix (`extract_enrichment_df()` in multiomics)
`extract_enrichment_df(rna_pathway_res)` failed with `numbers of columns of arguments do not match` (call `rbind(deparse.level, ...)`).

**Exact location:** `R/domain/multiomics/07_enrichment.R`, function `extract_enrichment_df()`, at `do.call(rbind, dfs)` (the `$pathway_results` branch). This file is **byte-identical to `origin/main`** — untouched by the migration.

**Root cause (downstream consumer, not enrichment).** The function walks `pathway_results`, harvests **every** leaf data.frame (both the per-contrast GSEA tables and the `cluster_ora` list of ORA tables) into one flat list, and row-binds them with base `rbind`, which requires identical column schemas. ORA and GSEA legitimately have **different columns** — verified at runtime: ORA tables 20 cols (`…, Fold_enrichment, Count, geneID, …`), GSEA tables 16 cols (`…, NES, enrichmentScore, core_enrichment, …`). The native multiomics producer (`run_kegg_enrichment_for_omics`) only ever emits **one** method per run (ORA *or* GSEA), so its results were always schema-homogeneous and the `rbind` assumption was never exercised. The **new offline RNA path is the first producer to place both ORA and GSEA tables under a single `pathway_results`** — a valid, contract-compliant shape the consumer never anticipated. Same category as §4.2: a latent assumption in unmodified downstream code, exposed by a new-but-valid object shape.

**Fix (minimal, consumer-only).** Replace `do.call(rbind, dfs)` with **`dplyr::bind_rows(dfs)`**, which aligns columns by name and NA-fills the missing method-specific columns (`dplyr` is already a core dependency). One line changed; enrichment producers, ORA outputs, GSEA outputs, and the Phase 1 contract are untouched. **Verified** against the live smoke target: `extract_enrichment_df(rna_pathway_res)` now returns a valid `data.frame` of **611 × 29** (union of the 20 ORA + 16 GSEA columns; 7 shared), with **348 ORA rows** (non-NA `Fold_enrichment`) + **263 GSEA rows** (non-NA `NES`) = 611, and the method-specific columns correctly NA-filled on the other method's rows (`NES` NA on all ORA rows; `Fold_enrichment`/`Count` NA on all GSEA rows).

### 4.4 Phase 1 validation — COMPLETE (results & conclusions)
All Phase 1 validation tasks are now done. Summary of what was confirmed:

**Core functional validation (via the smoke config — see §4.1):**
- ✅ Unit tests pass (`testthat`, `enrichment-local`).
- ✅ `tar_validate()` passes (pipeline parses; DAG well-formed).
- ✅ Full RNA single-omics smoke pipeline completes end-to-end and is **green** (`tar_meta()` no errored targets; `tar_outdated()` empty).
- ✅ Offline local (table-driven) enrichment activates and runs (GO_CC tables loaded).
- ✅ **ORA** works (all 5 collections in single-omics: `all_DE`, `contrasts`, `contrasts_wo_direction`, `partition`, `binary_patterns`; correct columns; fold-enrichment computed).
- ✅ **GSEA** works (per ranking method × contrast, incl. `any_contrast`).
- ✅ Executive summary, Shiny payload, and HTML report all generate.

**Downstream-consumer fixes (both verified):**
- ✅ `rna_pipeline_summary` integer/double bug — fixed (§4.2) and the summary target completes.
- ✅ `extract_enrichment_df()` `rbind` schema-mismatch bug — fixed (§4.3); reader returns a valid 611 × 29 data.frame with both methods and correct NA-fill.

**Multiomics (`clustering_res = NULL`) behavior — validated:**
- ✅ GSEA is produced.
- ✅ DE-derived ORA is produced (`all_DE`, `contrasts`, `contrasts_wo_direction`).
- ✅ Partition ORA is correctly omitted.
- ✅ Binary-pattern ORA is correctly omitted.
- (Matches the §9.6 project decision: GSEA + DE-derived ORA when clustering is unavailable, not GSEA-only.)

**Online regression (`annotation_dir: ""`) — validated:**
- ✅ Empty `annotation_dir` correctly routes to the legacy **online** implementation ([05_mod_pathway.R:57](../R/modules/rnaseq/05_mod_pathway.R#L57) gate).
- ✅ The local (offline) enrichment path is correctly bypassed (local banner absent).
- ✅ Expected legacy behavior for a **non-model organism** observed (organism auto-detect degrades/limits, as `origin/main` does) — backward compatibility preserved (online branch is byte-identical to `origin/main`).

**Conclusion:** Phase 1 (RNA-seq local enrichment core) is functionally complete and validated against the legacy behavioral intent. The only outstanding Phase 1 action is the **commit decision** (everything is intentionally uncommitted). The next engineering task before Phase 2 is **parallel enrichment execution (§12)**.

---

# 5. Architectural Decisions (decision · rationale · legacy relationship)

| # | Decision | Rationale | vs legacy |
|---|---|---|---|
| A1 | **Offline, local-table-driven**; online path kept only as fallback when `annotation_dir` unset | Reproduce legacy; no network at runtime | **Matches** |
| A2 | **Single `rna_pathway_res` target**, clustering wired conditionally on `skip_outputs` (`NULL` in multiomics, `rna_clustering_obj` in single-omics). When clustering is `NULL` the run produces **GSEA + DE-derived ORA** (`all_DE`/`contrasts`/`contrasts_wo_direction`); only the clustering-derived collections (`partition`/`binary`) are omitted — see §9.6 | Stable name for multiomics + consumers; no duplicate GSEA; reuses the existing factory pattern. DE-derived ORA is useful biology even without clustering, so it is retained rather than reducing to GSEA-only | **Intentional divergence** from the frozen SPEC's "clustering-NULL → GSEA only" (project decision, §9.6) |
| A3 | **`pre = rna_pre`** | Current main's DE/pathway input; the reference branch's `rna_batch_corr` predates the current structure | Adapts to current arch |
| A4 | **Sequential execution; no `future.apply`; no `renv.lock` change** | Avoid a dependency/lockfile change; determinism | Neutral (perf only) |
| A5 | **Preserve legacy output directory layout** (`Enrichment/GSEA/...`, `Enrichment/ORA/...`) | Downstream/report path expectations | **Matches** |
| A6 | **Fail-soft**: missing DB → skip+warn; low/empty overlap → valid empty result; never crash | Robustness in a pipeline; legacy would error | **Adds** robustness (intentional) |
| A7 | **GO simplify default OFF**, opt-in via explicit pre-installed `enrichment.orgdb` (no auto-detect, no online); **regular GO ORA always produced** | Avoid implicit organism-specific GOSemSim/OrgDb dependency and runtime failures | **Intentional divergence** (legacy simplified by default) |
| A8 | **Per-pathway GSEA artifacts configurable, default OFF** (`gsea_per_pathway_artifacts`) | Legacy emitted them only on demand; avoid thousands of files by default | **Matches** legacy default |
| A9 | **`max_terms_in_dotplot` configurable, default 20** | Legacy display default; display-only | **Matches** |
| A10 | **GO expansion treated as preprocessing**; runtime expects pre-expanded tables (documented, not coded) | Legacy expands during annotation prep, not runtime | **Matches** |
| A11 | **Phased + gated; no commit/push/merge without approval; reference branch read-only** | Reviewability; protect the live pipeline | Process |
| A12 | **Heterogeneous ORA + GSEA tables are an expected, supported part of the `pathway_results` contract** — downstream consumers must tolerate differing column schemas (align by name, NA-fill), not assume a single enrichment schema | The offline path legitimately emits both methods in one result; ORA and GSEA have genuinely different columns. Enforcing one schema would drop method-specific information | **New** (current arch); enabled the §4.3 fix |

**GO simplify — detailed findings (supporting A7):** the OrgDb requirement comes only from `clusterProfiler::simplify()` (ORA via `enricher()` over local tables needs no OrgDb). `simplify()` reduces redundant GO terms via GOSemSim semantic similarity, which needs `godata(OrgDb, ont)`; an `enricher()` result lacks `@organism/@ont/@keytype`, so we build and pass `semData = GOSemSim::godata(<orgdb>, ont)` explicitly when enabled. No purely local-table alternative reproduces semantic simplify (the repo's `cluster_enrichment_terms()` Jaccard grouping is a *different*, non-equivalent method).

---

# 6. Audit-driven Corrections (compatibility fixes; issue · why it mattered · resolution)

These were introduced specifically to better match legacy after the behavioral + clustering audits.

- **`all_DE` ORA collection** — *Issue:* `build_gene_lists()` lacked the legacy single "all DE genes" list. *Why:* it is a required default ORA input. *Resolution:* added a union of DE genes across all contrasts in a single `"all"` cluster.
- **GO Expansion documentation** — *Issue:* nothing stated that GO expansion is a preprocessing step. *Why:* a user supplying *un*-expanded GO tables would silently diverge from legacy. *Resolution:* config documents that supplied GO tables must already be expanded (e.g. via `buildGOmap`); runtime does not expand.
- **Configurable GSEA per-pathway artifacts** — *Issue:* artifacts were generated automatically for every significant pathway. *Why:* legacy produced them only on demand; auto-generation can create huge numbers of files. *Resolution:* `gsea_per_pathway_artifacts` (default `false`); regular GSEA tables/dotplots still produced.
- **Configurable dotplot limit** — *Issue:* dotplot term count was hardcoded (1000). *Why:* legacy default is 20; calculations must be unaffected. *Resolution:* `max_terms_in_dotplot` (default 20), threaded into ORA (`run_cluster_ora` `maxCategory`) and GSEA dotplots; display-only.
- **M1 — `Gene:` prefix alignment** — *Issue:* the module strips a leading `"Gene:"` from DE-derived IDs but not from clustering-derived IDs. *Why:* for prefixed datasets (e.g. Trinity), partition/binary genes wouldn't intersect `TERM2GENE` → silent empty ORA. *Resolution:* `build_gene_lists()` strips the same prefix from partition/binary IDs (no-op when absent).
- **M2 — partition source correction** — *Issue:* `objects$clusters` is overwritten (hierarchical → partition); a hierarchical-only run could be mislabeled as `partition`. *Why:* legacy never ran ORA on hierarchical cuts; mislabeling produces a wrong collection. *Resolution:* prefer `excel_order$partition_clusters` (set only when partition ran); fall back to `objects$clusters` only when no `excel_order` exists (partition-only run); never treat hierarchical cuts as partition.
- **M3 — robustness cleanup** — *Issue:* binary guard used `length()` on a data.frame (counts columns); misleading "named vector" comment. *Why:* correctness/clarity. *Resolution:* `is.data.frame(...) && nrow(...) > 0`, drop NA patterns, corrected comment. (Partition round-key left as `"k"` to avoid changing output basenames.)

All audit corrections are **resolved** and covered by unit tests.

---

# 7. Configuration Reference (`modes.rna.enrichment` in `rna_config.yaml`)

| Option | Purpose | Default | Legacy compatibility | Behavior |
|---|---|---|---|---|
| `annotation_dir` | Activation switch + path to local `.tab` tables (abs or relative to `project.dir`) | `""` | n/a (current-arch switch) | Empty → existing **online** pathway workflow runs unchanged; set → **offline** local path activates |
| `databases` | Subset of DBs to load | all four: `KEGG, GO_BP, GO_MF, GO_CC` | Matches (legacy uses KEGG + GO BP/MF/CC) | Missing files for a DB → skip that DB + warn; others still run |
| `go_simplify` | Enable GO term simplify | `false` | **Intentional divergence** (legacy always simplified GO) | When `true` **and** `orgdb` is installed → also emit `Simplify_*`; regular GO ORA always produced regardless |
| `orgdb` | Local, pre-installed OrgDb package for simplify only (e.g. `"org.Hs.eg.db"`) | unset | Legacy relied on an implicit OrgDb | Used only to build `GOSemSim::godata` for simplify; no auto-detect, no online |
| `gsea_per_pathway_artifacts` | Emit per-pathway gseaplot2 PNGs + core-gene CSVs for every significant pathway | `false` | Matches legacy default (on-demand only) | `true` → emit for all significant pathways; regular GSEA tables/dotplots always produced |
| `max_terms_in_dotplot` | Pathways shown in dotplots | `20` | Matches legacy `MAX_TERMS_IN_DOTPLOT` | **Display-only**; does not affect calculations; all results stay in output tables |
| `pvalue_cutoff` | ORA/GSEA adjusted-p cutoff | `0.05` | Matches `ENRICHMENT_PVAL_CUTOFF` | Passed to enricher/GSEA |
| `padj_method` | p-adjust method | `fdr` | Matches `ENRICHMENT_PADJ_METHOD` | Passed to enricher/GSEA |
| `gsea_pvalue_cutoff` / `gsea_padj_method` | GSEA-specific overrides | fall back to the two above | Matches | Optional overrides |
| `workers` | Parallel GSEA workers | `1` (sequential) | Neutral | `>1` needs `future` + `future.apply` (not enabled; do not set without adding deps) |

> **GO Expansion (preprocessing):** the supplied `GO2gene_*.tab` tables must already contain expanded GO assignments (e.g. produced with `clusterProfiler::buildGOmap` during annotation prep). Runtime enrichment does **not** expand GO. Expected files in `annotation_dir`: `KEGG_pathway2gene.tab`, `KEGG_pathway2name.tab`, `GO2gene_{BP,MF,CC}.tab`, `GO2name_{BP,MF,CC}.tab` (gene-based; two-column, header auto-detected).

---

# 8. Remaining Work

### Phase 1 — finish criteria

**✅ Phase 1 validation is COMPLETE** — full results and conclusions are recorded in **§4.4**. In brief: unit tests, `tar_validate()`, full green smoke run, offline ORA + GSEA, exec summary / Shiny payload / HTML report, both downstream-consumer fixes (§4.2, §4.3), multiomics `clustering_res = NULL` behavior (GSEA + DE-derived ORA; partition & binary ORA omitted), and online-regression routing (`annotation_dir: ""` → legacy online path bypassing local, expected non-model-organism behavior) all passed. The §9.6 scope question is resolved by project decision (§9, item 6).

**⏳ Only outstanding Phase 1 action:**
- **Commit decision** — all Phase 1 work is intentionally uncommitted; decide when to commit (and whether to push the branch). This is the sole remaining item before the §12 parallel-execution task and Phase 2.

### Phase 2 — Legacy enrichment plots *(gated)*
- **Objectives:** port/adapt the legacy plotting helpers — `ridgeplot_edited()` / `ridgeplot_edited1()` (GSEA ridgeplots), `plot_shared_genes()` (gene-term + term-term Jaccard heatmaps), `plot_expression_heatmap_for_precomputed_clusters` (core-gene / all-gene heatmaps, incl. z-score / manual-cluster joins).
- **Dependencies:** these helpers exist **only in the external legacy Rmd** (`Functions_for_clustering_and_enrichment_1.0.R`) — absent from the reference branch and current main; that script is the source.
- **Investigations done:** confirmed which helpers exist and where; identified reusable infra already present — `save_gsea_per_pathway_artifacts()` (reference branch), `generate_pathway_plots()`, `generate_clustered_dotplots()` (core/09), heatmap helpers in `R/core/06_plots.R`.
- **Prerequisites:** Phase 1 complete ✅ (§4.4) + commit decision; the §12 parallel-execution task recommended first; new plotting helpers added (likely `R/core/09_enrichment.R` or `R/core/06_plots.R`) and called from `.run_local_enrichment()`; outputs written into the existing legacy dirs; smoke/dimension tests + visual review. **Full preparation in §13.**

### Phase 3 — Proteomics adaptation *(gated; investigation-first)*
- **Objectives:** offer the same offline local enrichment for proteomics.
- **Investigations done (`R/domain/proteomics/07b_pathway.R`):** proteomics enrichment runs on **gene symbols mapped from proteins** — `map_proteins_to_gene_symbols()` (first symbol from the `Genes` column; dedups keeping lowest `padj`; stores `ProteinID`) and `extract_de_table_for_pathway()` (wide `*.imputs.<cn>` → standard `FeatureID/log2FoldChange/pvalue/padj/stat`). Local KEGG/GO `.tab` are **gene-based** → compatible *after* mapping. `run_proteomics_pathway()` already returns the same nested `pathway_results` (+`padj`).
- **Dependencies / open questions:** does proteomics have `prot_clustering_obj`, and does it cluster on protein IDs or mapped gene symbols (cluster-ORA needs the gene-symbol space)? `Genes` population/NA rate → overlap threshold? Reuse the RNA local-enrichment core (factor omics-agnostic) vs a proteomics-specific wrapper? `fc` ranking parity for `linearFC.imputs`? Wiring in `R/modules/proteomics/05_mod_pathway.R` + the proteomics pathway target (mirror the single-target/conditional-clustering decision).
- **Prerequisites:** Phase 1 complete; the above questions answered + approach approved.

### Phase 4 — HTML report integration *(planning only; deferred)*
- **Objectives:** surface enrichment in the generated HTML report (not Shiny): GSEA tables (method × contrast × DB), cluster-ORA tables (incl. `Cluster`, `Fold_enrichment`), GO simplified vs unsimplified, dotplots, and (later) ridgeplots/heatmaps/shared-gene plots.
- **Dependencies:** likely `R/domain/rnaseq/report_template*.Rmd` and `render_rnaseq_report()` (consumed by the `rna_report` target).
- **Prerequisites:** Phases 1–2 outputs stable; keep report blocks additive/guarded so absent outputs never break existing reports. **No template edits until Phase 4 is approved.**

### Out of scope (legacy has, we deliberately omit for now)
`manual_clustering` (×2) and `unified_clusters` ORA collections — legacy includes them, but manual clustering is off by default (`PERFORM_MANUAL_CLUSTERING=F`) and they are not in the agreed Phase-1 defaults.

---

# 9. Open Decisions (genuine, still pending)

1. **Commit timing** — Phase 1 changes are intentionally uncommitted; decide when to commit (and whether to push the branch).
2. **Phase 2 go/no-go** — approve starting the plots port (requires fetching helpers from the external legacy Rmd).
3. **Phase 3 proteomics approach** — answer the §8 Phase-3 questions and approve reuse-core-vs-wrapper before implementing.
4. **`manual_clustering` / `unified_clusters`** — decide whether these legacy ORA collections should ever be added (currently out of scope).
5. **Parallel enrichment — ✅ RESOLVED & IMPLEMENTED (§12).** `future` + `future.apply` (+ `globals`, `listenv`) were added to `renv.lock` via `renv::install()` + `renv::record()` (approved). `enrichment.workers` now controls a generic, method-agnostic parallel orchestration layer for both ORA and GSEA; `workers = 1` and `workers = 4` produce identical results.
6. **Multiomics ORA scope vs SPEC (§9.6) — ✅ RESOLVED (project decision, Phase 1).** The frozen SPEC said `clustering_res = NULL` (multiomics/`skip_outputs=TRUE`) → **GSEA only, all ORA skipped + warning**. **Decision: keep the current behavior — GSEA + DE-derived ORA** (`all_DE`, `contrasts`, `contrasts_wo_direction`) when clustering results are unavailable; only the clustering-derived collections (`partition`, `binary`) are omitted. We are **not** reducing to GSEA-only. **Rationale:** (1) DE-derived ORA provides useful biological information even without clustering; (2) this behavior is acceptable for the current architecture; (3) the correct place to accommodate mixed results is the **downstream consumers**, which should support heterogeneous ORA + GSEA tables (align by column, NA-fill — see A12 and §4.3) rather than assuming a single enrichment schema. This chooses option (a) from the earlier framing and supersedes the SPEC's "GSEA-only" wording for Phase 1. No producer code changed; the §4.3 consumer fix makes this behavior fully consumable. *(Originally discovered while planning the fast validation.)*

(All earlier planning-phase decisions — offline switch, single-target wiring, `databases` default, overlap threshold, config coexistence, test fixtures, GO-simplify default — are **resolved and implemented**; see §5–§7.)

---

# 10. Session Knowledge (discoveries worth preserving)

- **Legacy enrichment is genuinely offline** and entirely local-table-driven; the `KEGGREST` import is a red herring (unused at runtime). Any "online" behavior is purely the *current* pipeline's, not legacy's.
- **GO Expansion lives in annotation preprocessing**, not runtime — the legacy `Func_annot_data_expanded` directory holds `buildGOmap`-expanded tables. This is easy to miss and silently changes results if the supplied tables aren't expanded.
- **`clusterProfiler::simplify()` is the sole OrgDb dependency** in the GO path; the ORA itself never needs an OrgDb. Legacy simplified GO by default but never passed an OrgDb explicitly — it leaned on an implicit local OrgDb.
- **Clustering object shapes (current pipeline):** `rna_clustering_obj$objects$clusters` is a **named integer vector** (feature → cluster); `…$objects$patterns` is the `run_binary_patterns()$best` **data.frame** (`feature_id, best_pattern, best_corr, …`); `…$excel_order$partition_clusters` is set **only** when partition ran. `objects$clusters` is overwritten hierarchical→partition (the source of the M2 ambiguity).
- **ID-space asymmetry:** the pathway module strips a `Gene:` prefix from DE-derived IDs but not from clustering-derived IDs — the root of M1. Only affects prefixed datasets (e.g. Trinity).
- **GSEA `fc`/`any_contrast` math** reproduces legacy exactly once linear FC is recovered from `log2FoldChange` — verified numerically.
- **Environment caveat:** in this workspace `testthat` segfaults and the renv library is incomplete; validate with base-R assertions + `parse()`/`tar_target` introspection here, and run the full `tar_validate()`/`testthat` suite only in a restored renv.
- **Smoke-validation strategy works.** GSEA runtime ∝ (ranking methods × contrasts × databases); ORA is cheap. A reduced config — **fewer contrasts + a single small GO DB (GO_CC)**, everything else identical — cut a ~5-hour run to minutes while still exercising nearly every Phase-1 code path (loader/skip, all 4 rankers, GSEA plumbing, all 5 ORA collections, module routing, pipeline wiring, all downstream consumers). It is a **validation-only** config (`analysis_round: Smoke_01`, separate output dir), **not** for biological analysis. This is the recommended way to re-validate without repeating the full 5-hour run.
- **Latent downstream bugs surface only under a full pipeline run.** The `rna_pipeline_summary` integer/double failure (§4.2, §10 root cause) lived in unmodified code (`collect_pipeline_stats()`), was invisible to unit tests and to every prior online run (which always emitted `NES`), and only appeared once the smoke run fed contract-valid **NES-less ORA tables** through the whole pipeline. Lesson: exercise the *complete* downstream chain, not just the enrichment functions — a new-but-valid object shape can break a consumer's hidden assumptions. Root cause detail: `count_from_df()` used `else 0` (double) for NES-less tables while the aggregator used `vapply(..., 0L)` (integer); the enrichment implementation was **not** at fault (ORA legitimately has no `NES`); fix = `else 0L`.
- **A second latent consumer assumed one enrichment schema.** `extract_enrichment_df()` (multiomics; unmodified vs `origin/main`) row-bound *all* harvested tables with base `rbind`, which needs identical columns. It never failed before because the native multiomics producer emits a **single** method per run (all-ORA or all-GSEA). The offline RNA path is the first to return **both ORA and GSEA under one `pathway_results`** — with legitimately different columns (ORA 20 cols incl. `Fold_enrichment/Count`; GSEA 16 cols incl. `NES/core_enrichment`) — so `rbind` aborted with "numbers of columns of arguments do not match". Fix = `dplyr::bind_rows()` (align by name, NA-fill). Same lesson as above; also drove the A12 contract clarification: **heterogeneous ORA+GSEA tables are expected and supported** — consumers must not assume a single schema. See §4.3.

---

# 11. Next Recommended Task

**Phase 1 is complete and fully validated (§4.4).** When development resumes:

- **Branch:** stay on `feature/enrichment-migration-v2` (do not work on `main`). Confirm the Phase-1 working-tree changes are present (uncommitted): `R/core/09_enrichment.R`, `R/modules/rnaseq/05_mod_pathway.R`, `R/pipeline/rnaseq/00_pipe_rnaseq.R`, `R/domain/rnaseq/10_pipeline_summary.R` (the §4.2 fix), `R/domain/multiomics/07_enrichment.R` (the §4.3 fix), `config/templates/rna_config.yaml`, `tests/testthat/test-enrichment-local.R` (+ fixtures).
- **Sole remaining Phase 1 action — the commit decision:** all Phase 1 work is intentionally uncommitted; decide whether/when to commit and push the branch.
- **Next engineering task (before Phase 2): parallel enrichment execution — see §12.** This comes first because it changes the enrichment engine itself; doing it before Phase 2 avoids re-validating plots against a later engine change.
- **Then Phase 2 (gated): legacy enrichment plots — see the §13 Phase 2 Preparation** and the Phase 2 entry in §8.
- **Do not touch yet:** HTML report / `report_template*.Rmd` (Phase 4), Shiny, proteomics (Phase 3), `renv.lock` (unless the §12 decision explicitly adds a dependency), reference branch `claude/enrichment-migration-plan-cuduc` (read-only). No commit/push/merge without explicit approval.

---

# 12. Parallel Enrichment Execution (IMPLEMENTED — before Phase 2)

**Status: IMPLEMENTED and validated (local & uncommitted). Byte-identical results for `workers = 1` vs `workers = 4`.**

**Goal (met).** Replace the sequential enrichment execution with configurable parallelism driven by `enrichment.workers`, while producing **identical biological results** regardless of worker count.

### What was implemented
A single **generic, method-agnostic orchestration layer** — `run_enrichment_jobs(jobs, fun, workers)` in `R/core/09_enrichment.R`. A "job" is any independent unit of enrichment work `(method, input-unit, database)`; the parallel logic lives in this one function, and every method funnels its job list through it. This is architecture-driven, not method-driven: GO, KEGG, ORA, GSEA, and any future method (Reactome, WikiPathways, MSigDB, …) parallelize simply by contributing jobs to the same mechanism — no new parallel plumbing per method.

Both existing methods now route through it:
- **GSEA:** `run_gsea_all()` already built a flat `(ranking method × contrast × database)` job list; its inline parallel block was replaced by a call to `run_enrichment_jobs()`.
- **ORA:** the module's previously-sequential triple-nested loop was refactored into a flat `(collection × round × database)` job list dispatched through `run_enrichment_jobs()`. To keep I/O serial, `run_cluster_ora()` was split into `run_cluster_ora_compute()` (pure compute — enricher/merge/fold-enrichment/optional simplify; no I/O, worker-safe) and `write_cluster_ora_outputs()` (serial CSV/PDF writer). `run_cluster_ora()` remains a thin wrapper (compute + write) so its public behavior/signature are unchanged.

### What is parallelized vs sequential
- **Parallel (compute only):** the independent per-job enrichment computations (`clusterProfiler::GSEA()` / `enricher()`), dispatched through the generic layer.
- **Sequential (always):** result assembly into the nested `pathway_results`; **all file writing** (CSV / PDF / PNG); dotplot/plot generation; per-pathway GSEA artifacts; table loading (`load_local_pathway_tables`); gene-list building. Results are assembled in **input-job order**, so object/list/data-frame ordering and file content are deterministic and independent of scheduling. **No plotting, CSV/PDF/PNG writing, Shiny-payload creation, or report generation ever happens inside a worker.**

### How `workers` controls execution (single control for ORA + GSEA)
`enrichment.workers` (read once in `mod_rnaseq_pathway`) is the sole knob for the whole enrichment engine. It selects only the *backend*, never the results:
- `workers <= 1` → `future::sequential` plan (one job at a time, in-process, no worker spawn).
- `workers > 1` → `future::multisession` with that many workers (separate processes; Windows-safe).
- **Both paths go through `future.apply::future_lapply(..., future.seed = seed)`** where `seed` is the project's `params$seed` (see §12.1). An **explicit integer** `future.seed` derives each job's L'Ecuyer-CMRG RNG stream from that fixed seed + job position — **identical regardless of backend, worker count, or ambient RNG state** — so permutation-based GSEA returns identical results for `workers = 1`, `4`, or any N, and across independent pipeline rebuilds. If `future`/`future.apply` are unavailable, it degrades to plain `lapply()` (sequential; no per-job streams) with a message.

**Reproducibility decision (important).** Early testing confirmed the RCA's concern: with `workers = 1` using plain `lapply` (ambient global RNG) and `workers = 4` using `future.seed` streams, **all 5 GSEA tables differed** (ORA was already identical — it has no RNG). Because `clusterProfiler::GSEA()` is permutation-based (`seed = FALSE`), the two RNG mechanisms diverge. The fix routes **every** worker count through `future_lapply()` with a shared RNG mechanism (sequential plan for `workers <= 1`). Initially this used `future.seed = TRUE`, which fixed worker-count invariance but still keyed off ambient RNG (non-reproducible across rebuilds); it was subsequently replaced by an **explicit integer `future.seed = params$seed`** (see §12.1) for full reproducibility. This intentionally supersedes "`workers = 1` == plain-`lapply` behavior": the pre-parallel path used un-seeded ambient RNG and was itself non-reproducible run-to-run, so it was not a stable target; the seeded path is a strict correctness improvement (and aligns with the repo's reproducibility mandate).

### Validation results
- **Unit tests:** 51 passed, 0 failed (`tests/testthat/test-enrichment-local.R`); the one warning is the intended invalid-cluster warning and its test passes — `run_cluster_ora()` still returns `list()` for invalid input (backward-compatible).
- **`tar_validate()`:** PASS.
- **`workers = 1` vs `workers = 4` (real `mod_rnaseq_pathway`, single-omics, same seed):** **identical** — ORA 7/7 leaves byte-identical, GSEA 5/5 leaves byte-identical, Enrichment CSVs 12/12 identical, same file sets. Verdict PASS. (Before the reproducibility fix: GSEA 5/5 differed — see above.)
- **Speedup:** `workers = 4` ≈ 190 s vs `workers = 1` ≈ 2210 s on the smoke config (~11×; GSEA-dominated).
- **Multiomics (`clustering_res = NULL`, workers = 4):** GSEA present ✓; DE-derived ORA present (`all_DE`, `contrasts`, `contrasts_wo_direction`) ✓; partition ORA absent ✓; binary-pattern ORA absent ✓ — matches the §9.6 decision.
- **Online routing (`annotation_dir = ""`):** local offline banner absent ✓; online path entered ✓ (organism auto-detect degrades for the non-model organism, exactly as `origin/main` does) — backward compatibility preserved.
- **Full smoke `tar_make` (post-change rebuild):** GREEN. `rna_pathway_res` rebuilt in the `{targets}` context (all 7 ORA jobs + GSEA jobs), then downstream `rna_pipeline_summary`, `rna_shiny_payload`, `rna_report`, `rna_exec_summary`, `rna_commentary_file` all built — **no errored targets**, `tar_outdated()` = `character(0)`. (One transient failure during validation was an environment-only issue — `pandoc` not on the headless `Rscript` PATH for `rna_report`; resolved by pointing `RSTUDIO_PANDOC` at RStudio's bundled pandoc. Not a code issue; in a normal RStudio session pandoc is already present.)
- **Downstream unchanged across worker counts:** because `rna_pathway_res` is byte-identical for `workers = 1` vs `4`, every downstream consumer (`rna_pipeline_summary`, `rna_shiny_payload`, `rna_report`, `rna_exec_summary`) is a deterministic function of it and is therefore identical between worker counts (differences vs the pre-parallel baseline are limited to GSEA's now-seeded p-values and runtime metadata; see limitations).

### Files changed
- `R/core/09_enrichment.R` — add `run_enrichment_jobs()`; refactor `run_gsea_all()` to use it; split `run_cluster_ora()` into `run_cluster_ora_compute()` + `write_cluster_ora_outputs()` + thin wrapper.
- `R/modules/rnaseq/05_mod_pathway.R` — hoist `workers`; refactor the ORA block into a job list dispatched through `run_enrichment_jobs()`, with serial writing/assembly; add the `.make_ora_worker()` factory that builds the compute worker with a **minimal captured environment** (only `gene_lists`, `local_tables`, and scalar params — ~5 MiB), so `future` serializes just the data the worker uses, not unrelated large objects (e.g. the expression matrix in `pre` / `de_res`) from the module frame.
- `config/templates/rna_config.yaml` — expand the `workers` doc comment (now the single ORA+GSEA control).
- `renv.lock` — added `future`, `future.apply`, `globals`, `listenv` **via `renv::install()` + `renv::record()`** (4 records only; the lockfile was not hand-edited, and an unrelated `renv::snapshot()` sweep was reverted to keep the change minimal).

### Known limitations
- **`workers = 1` is slower than the pre-parallel plain-`lapply` path** (future's sequential backend adds per-job overhead). This is the deliberate cost of routing all worker counts through one RNG-reproducible mechanism.
- **GSEA numeric values differ from the pre-parallel cached baseline** because GSEA is now reproducibly seeded (`future.seed`) instead of using ambient RNG. **ORA is unchanged**, and the `pathway_results` **structure** (keys/schemas) is unchanged, so all downstream consumers work unmodified. Significant-pathway *counts* may shift slightly — a reproducibility improvement, not a parallelization artifact.
- **Memory** scales with worker count on `multisession` (each worker holds a copy of the gene-set tables); set `workers` to fit RAM.
- **`future.globals.maxSize` left at its 500 MiB default (by design).** The `.make_ora_worker()` factory keeps exported globals to ~5 MiB (measured: `local_tables` 2.4 + `gene_lists` 2.1 + scalars), verified to run under the default limit — so no bump is needed. The guard is deliberately *not* raised: it is a useful early warning if a future method ever starts broadcasting large objects. If a real dataset ever legitimately approached 500 MiB, the right response is to reconsider the data-passing architecture (don't broadcast large tables to every worker), not to silence the guard.
- **Worker startup cost** (~1–2 s/worker) makes tiny runs not benefit; the win grows with GSEA job count.
- If `future`/`future.apply` are absent at runtime, execution silently degrades to plain `lapply` (identical structure, no per-job RNG streams, no speedup).
- **RNG scope of the guarantee — fully resolved (§12.1).** The orchestration layer now passes an **explicit integer** `future.seed = params$seed`, so GSEA results are invariant to worker count **and** independent of the ambient RNG state — identical across separate pipeline rebuilds, not merely within one run. (The original `future.seed = TRUE` keyed off ambient RNG, which drifted between builds; see §12.1 for the fix.)
- **Rendering tooling (environment, not code):** `rna_report` needs `pandoc` (via `rmarkdown`). A headless `Rscript` without `RSTUDIO_PANDOC` set will fail that target; a normal RStudio session (or setting `RSTUDIO_PANDOC` to RStudio's bundled pandoc) renders it fine. Unrelated to enrichment/parallelism.

**Untouched.** Enrichment parameters, collections, ranking math; the Phase 1 downstream contract (incl. A12 heterogeneous tables); the online path; the report/Shiny consumers (structurally). This is a performance + reproducibility change with a zero-behavior-change guarantee *between worker counts*.

## 12.1 Enrichment RNG architecture (`params$seed` → deterministic GSEA)

**RNG architecture (single source of truth).** The project's one seed governs enrichment reproducibility:

```
config$params$seed
        ↓  (mod_rnaseq_pathway: enr_seed <- config$params$seed %||% 1L)
run_enrichment_jobs(jobs, fun, workers, seed = enr_seed)      # ORA and GSEA
        ↓
future.apply::future_lapply(jobs, fun, future.seed = seed)   # explicit integer
        ↓
clusterProfiler::GSEA(...)  # fgsea permutations use the future-assigned stream
```

`params$seed` is threaded through `mod_rnaseq_pathway()` into **both** enrichment dispatches (the ORA job list and, via `run_gsea_all(..., seed=)`, the GSEA job list). ORA (`enricher`) has no RNG so it is unaffected; GSEA (`fgsea`, `seed = FALSE`) uses the per-job L'Ecuyer stream that `future.seed` assigns.

**Old behavior.** `future_lapply(..., future.seed = TRUE)` derived per-job RNG streams from the **ambient `.Random.seed` at call time**. That ambient state drifts between builds (targets' per-target seed, RNG consumed by the earlier ORA dispatch, RNG-kind switches, etc.), so GSEA permutation p-values — and thus which borderline pathways cross `padj < 0.05` — varied across independent pipeline rebuilds. It was invariant *within* a run (hence `workers = 1 ≡ 4` held), but **not** across rebuilds. Compounding this, `params$seed` was **not wired to enrichment at all** (only proteomics imputation read it), so the project's declared seed had no effect on GSEA.

**New behavior.** `future_lapply(..., future.seed = seed)` with an **explicit integer** derives streams from that fixed seed + job position — **independent of ambient RNG, backend, and worker count**. GSEA is now identical across worker counts *and* across independent rebuilds, and `params$seed` is the single control.

**Why `future.seed = TRUE` was replaced.** `TRUE` = "seed from current RNG state" (ambient-dependent, non-reproducible across builds); an integer = "seed from this fixed value" (ambient-independent, fully reproducible). Verified directly: two `future_lapply` calls under *different* ambient states give identical results with `future.seed = <int>` but different results with `future.seed = TRUE`.

**Why this improves reproducibility.** Enrichment results now depend only on inputs + `params$seed`, nothing implicit. Re-running or rebuilding the pipeline yields byte-identical GSEA. Method-agnostic: any future RNG-using enrichment method routed through `run_enrichment_jobs()` inherits the same guarantee.

**One-time baseline change (expected).** Because GSEA is now seeded deterministically from `params$seed` (previously ambient/unseeded), the **committed Phase 1 GSEA baseline shifts once** to this stable, reproducible set. **ORA is unchanged** (no RNG); `pathway_results` **structure** and all downstream contracts are unchanged. After this one-time shift, subsequent rebuilds are byte-identical.

**Validation (this change) — PASS.** Three enrichment runs on the smoke config, each with a *different* ambient RNG state (`set.seed(999/1/424242)` before the call): (1) `workers = 1` (ambient 999) ≡ `workers = 4` (ambient 1) — all leaves byte-identical; (2) **ambient-independent** — `workers = 4` ambient 1 ≡ ambient 424242 (the "independent rebuild" property); (3) `workers = 1`/ambient 999 ≡ `workers = 4`/ambient 424242 — identical across *both* axes; (4) ORA identical across all three (deterministic). GSEA leaf counts were **stable at 5/5/5** (vs the 1-vs-2 drift observed under the old `future.seed = TRUE`). Conclusion: enrichment results now depend only on inputs + `params$seed`.

## 12.2 No nested parallelism — `fgsea` runs serially inside each GSEA job

**Root cause (found during RStudio validation).** `clusterProfiler::GSEA()` has no `BPPARAM` argument, so `fgsea` uses `BiocParallel::bpparam()`, which on **Windows defaults to `SnowParam` with (here) 14 SOCK workers**. Every GSEA job therefore spawned — and tore down — a 14-worker PSOCK cluster, *independent of* our `future` plan or `workers` setting. That is **nested parallelism**: `future.apply` fans out GSEA jobs, and each job then spun up its own SOCK cluster.

**Why it was problematic in our architecture.** The nested SOCK cluster is fragile and environment-dependent: in headless `Rscript` it usually worked (which is why earlier validation "passed"), but in an **RStudio session** the master→worker `serialize(data, node$con)` send failed with *"error writing to connection"* (`sendData.SOCKnode`). Because GSEA dispatch is fail-soft, every GSEA job was silently swallowed and `rna_pathway_res` completed with **only `cluster_ora`** (ORA has no `BiocParallel`, so all ORA jobs succeeded). Even when it didn't fail outright, the flaky cluster **silently dropped jobs** (a smoke run showed 3 GSEA leaves where 5 were expected). With `workers > 1` it would be strictly worse — `future` multisession workers each spawning a 14-worker SOCK cluster. The repeated "project is out-of-sync" renv messages were child SOCK workers each autoloading renv on startup.

**Fix.** Pass `BPPARAM = BiocParallel::SerialParam()` to **both** `clusterProfiler::GSEA()` call sites (`run_gsea_all()`'s per-job worker, and `run_gsea_local()`). Each GSEA job now runs `fgsea` **in a single process**; the *only* parallelism is the outer `future.apply` job fan-out controlled by `workers`. Execution model:

```
future.apply distributes independent GSEA jobs   (controlled by `workers`)
        ↓
each job calls clusterProfiler::GSEA()
        ↓
fgsea runs with BiocParallel::SerialParam()       (no inner SOCK cluster)
```

No new config option; no new dependency (`BiocParallel` is already a `fgsea`/`clusterProfiler` dependency); ORA untouched; the §12.1 RNG architecture (`params$seed` → integer `future.seed`) is unchanged and still the sole GSEA RNG source.

**Biological results unchanged (no new baseline shift).** `SerialParam` and the old `SnowParam` produce **byte-identical** GSEA results per job (verified: same result digest); the fix only makes all jobs run *reliably*. It is not a numeric baseline shift — it recovers jobs the flaky cluster was dropping.

**Validation — PASS.** `parse` + `tar_validate` + unit tests (51/0) clean. `tar_invalidate(rna_pathway_res)` + rebuild: all 7 ORA and all GSEA jobs completed; **no "error writing to connection"; no trailing `SOCKnode` serialize error**; renv child-worker messages dropped from ~14×7 to just the main process; `names(pathway_results)` = `cluster_ora, G.vs.Sy, any_contrast, Mud.vs.Sw` with **5 GSEA leaves** (was silently ORA-only in RStudio before). `workers = 1` ≡ `workers = 4` — 12 leaves each, all byte-identical (and `workers = 4` no longer nests SOCK clusters). Exposed during independent **RStudio** validation; committed as an isolated robustness fix.

---

# 13. Phase 2 Preparation (start tomorrow)

**What Phase 2 accomplishes.** Port/adapt the **legacy enrichment plots** so the offline path emits the same visual artifacts the legacy "Neat RNA-Seq" workflow produced, written into the existing `Enrichment/…` layout. Target helpers:
- `ridgeplot_edited()` / `ridgeplot_edited1()` — GSEA ridgeplots.
- `plot_shared_genes()` — shared-gene plots (gene↔term and term↔term Jaccard heatmaps).
- `plot_expression_heatmap_for_precomputed_clusters` — core-gene / all-gene expression heatmaps (incl. z-score and manual-cluster joins; the heatmap piece is a deferred sub-item).

**Files expected to be modified.**
- `R/core/09_enrichment.R` (or `R/core/06_plots.R`) — add the new plotting helpers (each with a roxygen2 docstring; match the existing `ggplot2`/plotting idiom).
- `R/modules/rnaseq/05_mod_pathway.R` — call the new helpers from `.run_local_enrichment()` (additive; existing routing/return shape preserved).
- `tests/testthat/` — file-existence / dimension / schema smoke tests for the new outputs.
- Possibly `config/templates/rna_config.yaml` — only if a plot toggle is needed (prefer reusing existing display options like `max_terms_in_dotplot`).

**Must remain untouched.**
- The enrichment **results** engine (ORA/GSEA computation, parameters, collections) — Phase 2 is plots only, reading existing results.
- The Phase 1 downstream contract and its consumers (`collect_pipeline_stats`, `extract_enrichment_df`, `get_pathway_highlights`), incl. the A12 heterogeneous-table guarantee.
- The online path, `renv.lock` (unless a new plotting dep is explicitly approved), the HTML report templates (Phase 4), Shiny, proteomics (Phase 3), and the read-only reference branch.
- The output **directory layout** (write into the existing legacy `Enrichment/GSEA/…`, `Enrichment/ORA/…`).

**Dependencies on completed Phase 1 work.**
- Phase 1 is complete + validated (§4.4); the offline enrichment results these plots consume already exist and are stable.
- The **legacy plotting source** lives only in the external legacy Rmd (`Functions_for_clustering_and_enrichment_1.0.R` in `Projects/Uri_Gat/`) — **absent** from the reference branch and current main; that script is the authoritative source and must be available before porting.
- Reusable infra already present: `save_gsea_per_pathway_artifacts()` (reference branch), `generate_pathway_plots()`, `generate_clustered_dotplots()` (core/09), heatmap helpers in `R/core/06_plots.R`.
- **Ordering note:** the §12 parallel-execution task should land first (see §12 rationale).

**Recommended implementation order.**
1. Fetch and read the legacy plotting helpers from the external Rmd; catalog exact inputs/outputs each expects.
2. Port GSEA **ridgeplots** first (self-contained; read GSEA result objects; lowest risk).
3. Port **shared-gene** plots (Jaccard gene↔term / term↔term).
4. Port **expression heatmaps** last (most inputs: expression matrix + cluster joins + z-score) — the deferred sub-item.
5. Wire each into `.run_local_enrichment()` additively, guarded so a missing/failed plot never breaks the pipeline (fail-soft, per A6).
6. Add smoke/dimension tests + do a visual review.

**Expected validation strategy.**
- Unit/smoke tests assert the expected plot files exist with sane dimensions and the right output paths.
- Run the smoke config end-to-end; confirm the new plots are produced and the pipeline stays green (no new errored targets, `tar_validate()` passes).
- **Visual review** of each new plot type against the legacy exemplars.
- Confirm additive-only: enrichment result tables and all Phase 1 downstream consumers are unchanged.

**Phase 2 kick-off checklist (for tomorrow).**
- [ ] Confirm on branch `feature/enrichment-migration-v2`; Phase 1 changes present/uncommitted (or committed per the §11 decision).
- [ ] Resolve the §12 vs §13 ordering with the user (recommended: §12 parallel execution first).
- [ ] Obtain the legacy plotting source (`Functions_for_clustering_and_enrichment_1.0.R`) and confirm it's readable.
- [ ] Catalog each legacy plot helper's inputs/outputs and map them to current result objects.
- [ ] Get explicit Phase 2 go/no-go (Open Decision §9, item 2).
- [ ] Draft the additive integration points in `.run_local_enrichment()` (no engine changes).
- [ ] Decide plot config toggles (prefer reusing existing options; avoid `renv.lock` changes unless a new dep is approved).
- [ ] Implement in the recommended order (ridgeplots → shared-gene → heatmaps), fail-soft.
- [ ] Add tests + run the smoke config; visual-review outputs; confirm green + unchanged downstream.
