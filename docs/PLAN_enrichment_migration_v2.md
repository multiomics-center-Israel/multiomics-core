# Current Status

Status: Planning

Current Branch:
feature/enrichment-migration-v2

Current Approved Phase:
None

Next Planned Phase:
Phase 1 – RNAseq Local Enrichment Core

Implementation Started:
No

Merged to Main:
No

---

This document is the authoritative source of truth for the enrichment migration project. Before making any implementation decisions, review this document first.

> _Maintenance note: keep the **Current Status** block above up to date as the project evolves — update `Status`, `Current Approved Phase`, `Next Planned Phase`, `Implementation Started`, and `Merged to Main` whenever a phase is approved, implementation begins, or anything is merged._

---

# Enrichment Migration V2 — Authoritative Project Plan

> Location: **`docs/PLAN_enrichment_migration_v2.md`** on branch `feature/enrichment-migration-v2`.
> Single source of truth for this migration across sessions. Planning only — no code until each phase is explicitly approved.

**Branch:** `feature/enrichment-migration-v2` — work here only. **No merge to main, no PR, no rebase/cherry-pick** until explicit approval.
**Reference (read-only):** `claude/enrichment-migration-plan-cuduc` — never merged; consulted only to re-apply/validate code.
**Companion docs (already on branch):** `docs/SPEC_enrichment_migration.md`, `docs/HANDOVER_enrichment_migration.md`, `docs/audit-enrichment-layer.md`.

---

## 1. Project goal

- **Behavioral preservation** of the legacy enrichment workflow (cluster-based ORA + multi-method GSEA), not line-by-line code.
- **Integration into the current architecture** (the `skip_outputs` core/single-omics pipeline split; the de-duplicated single `R/core/09_enrichment.R`).
- **Offline / local-table-driven** enrichment: KEGG + GO `.tab` (TERM2GENE / TERM2NAME) only.
- **No runtime online resources**: no KEGG REST, no biomaRt, no OrgDb-driven pathway loading, no GMT as the primary source. The existing online path remains only as an untouched fallback when `annotation_dir` is unset.

## 2. Architecture decisions

- **Single `rna_pathway_res` target** (no extra pathway targets).
- **`clustering_res = NULL` in multiomics** (`skip_outputs = TRUE`) → GSEA only, ORA skipped with a clear warning.
- **`clustering_res = rna_clustering_obj` in single-omics** (`skip_outputs = FALSE`) → GSEA + cluster-based ORA.
- **`pre = rna_pre`** (not `rna_batch_corr`).
- **Preserve downstream compatibility**: return `list(annotation, pathway_results, plot_files)` where `pathway_results` is always a (possibly empty) named list of data.frames carrying a `padj` column — verified against `collect_pipeline_stats()`, `extract_enrichment_df()` (multiomics), `get_pathway_highlights()`.
- Conditional target construction lives in the `pipe_rnaseq()` factory (already the established pattern there).

## 3. Confirmed implementation decisions

- **Phased** implementation (not all at once); each phase separately approved.
- **Sequential execution**; **no `future.apply`**; the parallel branch stays present but guarded off; **no `renv.lock` change** now.
- **Preserve the legacy output directory structure**; do not redesign layout without a clear compatibility reason.
- **No merge to main, no PR**; work only on `feature/enrichment-migration-v2`.
- **No HTML report changes now** (Phase 4, planning only).
- **Shiny payload/export out of scope** unless required for compatibility.

## 4. Error-handling decisions

- **Missing annotation database/file** → skip only that DB, emit a clear warning, continue with the rest (e.g. KEGG + GO_BP + GO_CC run, GO_MF skipped).
- **Low/empty gene overlap** between local tables and pipeline FeatureIDs → warn and return a **valid empty** result structure (no crash).
- **No pipeline crash** from enrichment: all failure modes degrade gracefully to warnings + valid (possibly empty) results.

## 5. GO simplify / OrgDb decisions (with findings)

**Findings (verified from the reference branch `run_cluster_ora` + clusterProfiler behavior):**
1. The OrgDb requirement comes **only** from `clusterProfiler::simplify()`. The ORA itself (`enricher()` over local TERM2GENE/TERM2NAME) needs **no** OrgDb — fully offline.
2. `simplify()` reduces redundant GO terms via **GOSemSim** semantic similarity, which needs `godata(OrgDb, ont)` (organism-specific OrgDb + GO.db). An `enricher()` result lacks the `@organism/@ont/@keytype` slots, so an explicit `semData = GOSemSim::godata(<OrgDb>, ont)` must be built/passed or `simplify()` fails (the reference branch wraps it in `tryCatch`).
3. No local-table-only alternative exists **within** `simplify()` (it needs the GO DAG/IC, absent from `.tab`). The repo's `cluster_enrichment_terms()` (Jaccard/overlap, offline) is a *different*, non-equivalent method — usable only as an explicitly-labeled alternative.
4. Legacy almost certainly relied on a **locally pre-installed OrgDb** loaded by the Rmd, undocumented as a hard dependency.

**Decisions:**
- **Unsimplified GO ORA is always produced** (offline, local tables).
- **`simplify()` is optional** and runs only if explicitly enabled **and** a pre-installed `enrichment.orgdb` is available (`requireNamespace`); build `semData` per ontology and pass it explicitly.
- **No online OrgDb loading. No organism auto-detection. No online fallback.** OrgDb is strictly a **local pre-installed optional dependency**, named explicitly in config.
- If disabled / OrgDb missing / simplify errors → skip simplify, keep unsimplified GO ORA, warn clearly.

## 6. Proteomics investigation (Phase 3 — investigation-first)

**Current findings (from `R/domain/proteomics/07b_pathway.R`):**
- Proteomics enrichment runs on **gene symbols mapped from proteins**, not protein IDs. `map_proteins_to_gene_symbols()` takes the `Genes` column (semicolon-separated → first symbol; strips `" | desc"`), remaps `FeatureID`→symbol, dedups keeping lowest `padj`, stores original in `ProteinID`.
- `extract_de_table_for_pathway()` converts wide `padj.imputs.<cn>` / `pvalue.imputs.<cn>` / `linearFC.imputs.<cn>` → standard `FeatureID, log2FoldChange, pvalue, padj, stat`; carries `Protein.Names / Genes / First.Protein.Description`.
- **Protein→gene mapping is required**; source = the `Genes` column in `de_res$summary_df` (DIA-NN/row_data). (Multiomics separately has `harmonization$gene_protein_mapping`; single-omics proteomics uses the `Genes` column.)
- Local KEGG/GO `.tab` are **gene-based** → compatible with proteomics **after** the existing mapping. No protein-based tables needed.
- `run_proteomics_pathway()` already returns the same nested `pathway_results` (df + `padj`) → reuse unchanged for downstream compatibility.

**Remaining questions (resolve before implementing Phase 3):**
1. Does proteomics have a clustering target (`prot_clustering_obj`)? Does it cluster on protein IDs or mapped gene symbols? (Cluster-ORA needs cluster gene sets in the **gene-symbol** space of the local tables.)
2. `Genes` population rate across real DIA-NN inputs (NA rate) → decide the overlap warn-threshold.
3. Reuse the RNA local-enrichment core (factor functions to be omics-agnostic) vs a proteomics-specific `.run_local_enrichment` wrapper.
4. `fc` ranking parity: confirm `linearFC.imputs` sign/scale vs RNA's `log2FoldChange`.
5. Wiring location: `R/modules/proteomics/05_mod_pathway.R` (+ `R/domain/proteomics/07b_pathway.R`) and the proteomics pathway target (mirror the single-target / conditional-clustering decision).

**Required validations:** protein→gene mapping correctness (dedup, NA handling); gene-space overlap with local tables; downstream `extract_enrichment_df()` parity; identical output schema.

## 7. Phase breakdown

### Phase 1 — RNA-seq local enrichment core *(implement after approval)*
Local KEGG/GO table loading (per-DB skip + overlap-guard) · 4 ranked gene lists (`pval_wo_direction`, `pval_with_direction`, `fc` with `2^log2FC` recovery + `signif(4)`, `min_pval_any_contrast`) · GSEA (`minGSSize=4`, `maxGSSize=length(unique(TERM2GENE[,2]))`, sequential) · **cluster-based ORA where clustering exists** (`enricher` `minGSSize=0/maxGSSize=10000/qvalueCutoff=1`, fold-enrichment, GO simplify optional/gated) · graceful `clustering_res=NULL` (GSEA-only + warning) · missing-DB handling · valid empty results · downstream compatibility · single-target pipeline wiring · config block · tests. Excludes plots/proteomics/report.

### Phase 2 — Legacy enrichment plots *(gated; NOT implemented now)*
Port + adapt the legacy plotting helpers. **These live only in the external legacy enrichment Rmd — confirmed absent from both the reference branch and current main.** Porting requires that external script as the source.

**Legacy helper functions to eventually port/adapt:**
- `ridgeplot_edited()` and `ridgeplot_edited1()` — GSEA ridgeplots.
- `plot_shared_genes()` — shared-gene plots across pathways/contrasts.
- Legacy **heatmap helpers** — core-gene and all-gene heatmaps (with z-score / manual-cluster joins, a deferred sub-item).

**Already present (reuse, not re-port):** `save_gsea_per_pathway_artifacts()` (reference branch) for per-pathway gseaplot2 PNGs + core-gene CSVs; `generate_pathway_plots()`, `generate_clustered_dotplots()` (core/09); heatmap helpers in `R/core/06_plots.R`.

**Integration targets:** new plotting helpers in `R/core/09_enrichment.R` (or `R/core/06_plots.R`), called from `.run_local_enrichment()` in `R/modules/rnaseq/05_mod_pathway.R`. **Outputs:** write into the existing legacy `Enrichment/GSEA/...` and `Enrichment/ORA/...` directories (preserve layout). **Testing:** file-existence smoke tests + dimension/schema checks + visual review.

### Phase 3 — Proteomics adaptation *(gated; investigation-first — see §6)*

### Phase 4 — HTML report integration *(planning only; deferred — no changes now)*
Generated **HTML report** (not Shiny). Likely files: `R/domain/rnaseq/report_template*.Rmd` and `render_rnaseq_report()` (consumed by the `rna_report` target). Surface: GSEA tables (method × contrast × DB), cluster-ORA tables (incl. `Cluster`, `Fold_enrichment`), GO simplified vs unsimplified, dotplots, and (later) ridgeplots/heatmaps/shared-gene plots. Keep additive/guarded so absent outputs never break existing reports. **No template edits until Phase 4 is approved.**

## 8. File-by-file strategy (Phase 1)

| File | Action | Detail |
|---|---|---|
| `R/core/09_enrichment.R` | APPEND | 12 local-enrichment functions (re-applied + re-validated); parallel guarded/off; per-DB skip + overlap + empty-result guards. No existing fn modified. |
| `R/modules/rnaseq/05_mod_pathway.R` | EXTEND (surgical) | `clustering_res = NULL`; local-vs-online routing on `annotation_dir`; `.run_local_enrichment()`, `.store_ora_result()`; preserve online fallback + return shape. |
| `R/pipeline/rnaseq/00_pipe_rnaseq.R` | REWIRE | Single `rna_pathway_res`; `pre = rna_pre`; conditional `clustering_res` by `skip_outputs`. |
| `config/templates/rna_config.yaml` | ADD | `modes.rna.enrichment:` block (`annotation_dir`, optional `databases`, `go_simplify`, `orgdb`, `workers=1`). |
| `tests/testthat/test-enrichment-local.R` | NEW | Loader (missing-DB + empty overlap), rankers (`fc`/`any_contrast`), GSEA params, cluster-ORA columns, `clustering_res=NULL` path, downstream `padj` contract. Synthetic fixtures only. |
| `R/domain/rnaseq/07_pathway.R` | LEAVE | Only `build_pathway_volcano_data()`; untouched. |
| `renv.lock` | NO CHANGE | Sequential; no `future.apply`. |

## 9. Open decisions (require approval before implementation)

1. **GO simplify config** — add `enrichment.go_simplify` (default `false`) + `enrichment.orgdb` (explicit pre-installed pkg, no auto-detect)? Default = unsimplified GO ORA only.
2. **Approval scope** — confirm approving this plan authorizes **Phase 1 only**; Phases 2/3/4 separately gated.
3. **`databases` default** — when unset, default to all four (KEGG, GO_BP, GO_MF, GO_CC), skipping any with missing files?
4. **Overlap warn-threshold** — warn + treat as empty below what fraction (SPEC suggested <5%)?
5. **Config coexistence** — keep `enrichment:` (local activation switch via `annotation_dir`) as a sibling of the existing online `pathway:` config? (Recommended.)
6. **Test fixtures** — commit synthetic `.tab` + DE-table fixtures under `tests/testthat/` (no real data)?

## 10. Verification (Phase 1, when implemented)
`targets::tar_validate()`; `testthat::test_file("tests/testthat/test-enrichment-local.R")`; targeted `tar_make(rna_pathway_res)` in both modes (single-omics → GSEA+ORA; multiomics/`skip_outputs=TRUE` → GSEA-only + warning); regression with `annotation_dir` unset → identical to current online behavior; `rna_exec_summary` / `rna_pipeline_summary` complete. Use `str()`/`dim()`/`colnames()` only (data-safety).
