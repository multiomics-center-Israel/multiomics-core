# Enrichment Migration — Handover Document

**Branch**: `claude/enrichment-migration-plan-cuduc`
**Status**: Not merged. 96 commits behind main. Merge conflicts exist.
**Date**: 2026-06-18
**Prepared for**: Colleague continuing development

---

## What this branch does

This branch implements a **behavior-preserving enrichment migration** for the RNA-seq pipeline. The goal is to replace the current online/dynamic gene set loading (OrgDb, KEGG REST, biomaRt) with a local, offline, table-driven enrichment path that reproduces the legacy enrichment workflow's behavior.

The branch adds two parallel enrichment approaches under a single config switch:

- **Local path** (new): activated when `modes.rna.enrichment.annotation_dir` is set in config. Loads precomputed KEGG/GO TERM2GENE/TERM2NAME tab files. Runs cluster-based ORA and multi-method GSEA. No online resources.
- **Online fallback** (unchanged): when `annotation_dir` is not set, the existing online gene set loading path runs exactly as before.

---

## What was implemented (11 commits)

### Design spec
- `docs/SPEC_enrichment_migration.md` — frozen design note with behavioral requirements, acceptance criteria, and implementation plan

### Phase 1: Local GSEA
- `load_local_pathway_tables()` — reads KEGG + GO tab files, validates structure, checks gene ID overlap
- Four ranked gene list builders matching legacy behavior:
  - `rank_by_pval_wo_direction()` — -log10(pvalue), unsigned
  - `rank_by_pval_with_direction()` — sign(FC) * -log10(pvalue)
  - `rank_by_fc()` — log2 of signed linear FC with signif(digits=4)
  - `rank_by_min_pval_any_contrast()` — cross-contrast minimum pvalue
- `build_ranked_gene_lists()` — convenience wrapper
- `run_gsea_local()` — clusterProfiler::GSEA() with legacy parameters (minGSSize=4, maxGSSize=n_unique_genes)
- `run_gsea_all()` — orchestrator with future.apply parallelization (Windows-safe)
- GSEA dotplots via enrichplot::dotplot() with ggplot2 fallback
- Per-pathway artifacts: gseaplot2 PNGs + core gene CSVs

### Phase 2: Cluster-based ORA
- `run_cluster_ora()` — per-cluster enricher() + merge_result() + GO simplify(), matching legacy Clusters_Enrichment_Test()
- `process_enrichment_table()` — fold enrichment computation (GeneRatio/BgRatio expansion)
- `build_gene_lists()` — constructs gene_lists[[method]][[round]] structure including:
  - `contrasts`: per-contrast up/down DE gene assignments
  - `contrasts_wo_direction`: per-contrast all DE genes
  - `partition`: clustering-derived assignments
  - `binary_patterns`: binary pattern assignments
- Input validation, edge case guards, legacy fidelity hardening
- Pipeline DAG wired: `clustering_res = rna_clustering_obj` added to `rna_pathway_res` target

### Phase 3: Per-pathway artifacts
- `save_gsea_per_pathway_artifacts()` — GSEA_plots/ and Excels_core_genes/ per significant pathway

---

## Files modified on this branch

| File | What changed |
|------|-------------|
| `R/core/09_enrichment.R` | +1016 lines appended: all new enrichment functions. No existing functions removed. |
| `R/modules/rnaseq/05_mod_pathway.R` | Rewritten: adds `clustering_res` parameter, local enrichment routing, ORA + GSEA orchestration |
| `R/pipeline/rnaseq/00_pipe_rnaseq.R` | Added `clustering_res = rna_clustering_obj` to `rna_pathway_res` target |
| `config/templates/rna_config.yaml` | Added `enrichment:` section with `annotation_dir` and `workers` keys |
| `docs/SPEC_enrichment_migration.md` | Design spec (new file) |

---

## CRITICAL: Merge conflicts with main

The branch is **96 commits behind main**. A test merge shows conflicts in:

### 1. `R/pipeline/rnaseq/00_pipe_rnaseq.R` — STRUCTURAL CONFLICT

**Main has restructured `pipe_rnaseq()`** to accept `skip_outputs = FALSE`. Targets are now split:
- **Core targets** (always run, even in multiomics): includes `rna_pathway_res`
- **Single-omics outputs** (skipped in multiomics): includes `rna_clustering_obj`

**Our branch** adds `clustering_res = rna_clustering_obj` to `rna_pathway_res`. But on main, `rna_pathway_res` is in core targets while `rna_clustering_obj` is in single-omics-only targets. This means **in multiomics mode, clustering wouldn't exist when pathway analysis runs**.

**Resolution needed**: The `rna_pathway_res` target on main passes `pre = rna_pre` (not `rna_batch_corr`). Our branch uses `rna_batch_corr`. During merge:
1. Use main's pipeline structure (`skip_outputs` pattern)
2. Accept `clustering_res = NULL` as default (our module handles this gracefully — runs GSEA only, skips ORA with a warning)
3. In the single-omics-only section, either:
   - (a) Move `rna_pathway_res` to single-omics and create a lightweight `rna_pathway_res_core` for multiomics, OR
   - (b) Keep `rna_pathway_res` in core (no clustering), and add a second target `rna_pathway_res_full` in single-omics that passes clustering. The module's return structure is the same either way.
   - (c) Simply don't pass clustering in core mode. ORA only runs in single-omics mode. This is the simplest option.

### 2. `R/core/09_enrichment.R` — MODERATE CONFLICT

Main has modified `load_gene_sets()` (removed "all" shorthand, simplified GO ontology loading) and added two new functions (`cluster_enrichment_terms`, `generate_clustered_dotplots`). Our branch appended new functions at the end.

**Resolution**: Accept main's changes to existing functions (our branch doesn't modify them). Our appended functions (everything after line ~660) should apply cleanly after main's additions since they're pure additions at the end.

### 3. `R/modules/rnaseq/05_mod_pathway.R` — MODERATE CONFLICT

Main made a small change: `pw_cfg` null guard changed from `pw_cfg %||% list()` to `pw_cfg` with an `is.null(pw_cfg) ||` check. Our branch rewrote the whole file.

**Resolution**: Take our branch's version. It's a full rewrite that includes the online fallback path (which matches main's behavior when `annotation_dir` is not set). Just verify the null-guard style matches main's convention.

### 4. `renv.lock` — CONFLICT

Main has many package updates. Our branch added `future` and `future.apply`.

**Resolution**: Take main's renv.lock, then add future/future.apply entries.

---

## How to merge

Recommended approach:

```bash
git checkout claude/enrichment-migration-plan-cuduc
git fetch origin main
git rebase origin/main
# Resolve conflicts file by file using the guidance above
```

Or if rebase is too complex (96 commits behind):

```bash
git checkout -b claude/enrichment-migration-v2 origin/main
# Cherry-pick or manually re-apply the enrichment changes
```

The safest approach may be to **start a fresh branch from main** and re-apply the enrichment additions, because:
1. All our new code is appended to `09_enrichment.R` (no modifications to existing functions)
2. The module is a full rewrite (take ours wholesale)
3. The pipeline DAG needs manual adaptation to main's `skip_outputs` structure
4. The config addition is trivial to re-apply

---

## What is NOT done yet

### Deferred from Phase 3
- **Heatmaps** (core genes, all genes) — directories planned but not implemented
- **Ridgeplots** — `ridgeplot_edited()` / `ridgeplot_edited1()` not ported
- **Shared gene plots** — `plot_shared_genes()` not ported
- **Z-score joins** and manual clustering joins for excel outputs

### Not started
- **Proteomics pathway module adaptation** — uses the same core functions but has its own protein-to-gene mapping
- **Output directory restructuring** — current outputs use a functional structure; legacy used a more hierarchical one
- **Report rendering integration** — current report templates may need updating to find enrichment outputs

### Known gaps vs legacy
1. **Hierarchical clustering as a separate ORA method**: the current clustering module overwrites `objects$clusters` when partition runs. Hierarchical cuts are lost. To fix: modify the clustering module to store them separately.
2. **Multiple partition k values**: legacy may have iterated multiple rounds. Current pipeline produces one.
3. **GO simplify requires GOSemSim + OrgDb**: fails silently if not installed. Unsimplified results still produced.

---

## Key design decisions (frozen in spec)

These were decided during design review and should not be reopened without good reason:

1. **compareClusterResult** is the preferred implementation vehicle for ORA (enables simplify + cluster-aware dotplot), not a strict requirement
2. **GSEA parity** is expected pending verification of three items: linearFC conversion, any_contrast reconstruction, signif(digits=4) rounding
3. **When clustering is disabled**: warn and skip ORA, run GSEA only. Do NOT fall back to per-contrast ORA as a substitute — it's a different analysis.
4. **Enrichment behavior** is the requirement, not the specific implementation. The spec preserves legacy analysis intent, not line-by-line code.

---

## Config reference

```yaml
modes:
  rna:
    enrichment:
      annotation_dir: "path/to/functional_annotation"  # activates local path
      # databases: ["KEGG", "GO_BP", "GO_MF", "GO_CC"]  # optional subset
      # workers: 4  # parallel GSEA workers (needs future + future.apply)
```

Expected files in `annotation_dir`:
```
KEGG_pathway2gene.tab    KEGG_pathway2name.tab
GO2gene_BP.tab           GO2name_BP.tab
GO2gene_MF.tab           GO2name_MF.tab
GO2gene_CC.tab           GO2name_CC.tab
```

Two-column tab-separated. Files may have headers (auto-detected).

---

## Testing checklist

When testing after merge:

1. **Online path unchanged**: run with no `enrichment.annotation_dir` set. Should produce identical results to main.
2. **Local path with GSEA only**: set `annotation_dir`, disable clustering. Should produce GSEA results + warning about ORA skipped.
3. **Local path with ORA + GSEA**: set `annotation_dir`, enable clustering. Should produce:
   - `Enrichment/ORA/{db}/{db}_{method}_{round}.csv` files
   - `Enrichment/ORA/{db}/{db}_{method}_{round}.pdf` dotplots
   - `Enrichment/ORA/{db}/Simplify_*.csv` for GO databases
   - `Enrichment/GSEA/{db}/ranking_by_{method}/{contrast}/GSEA_results_{contrast}.csv`
   - `Enrichment/GSEA/{db}/ranking_by_{method}/{contrast}/GSEA_plots/*.png`
   - `Enrichment/GSEA/{db}/ranking_by_{method}/{contrast}/Excels_core_genes/*.csv`
4. **Downstream consumers**: verify `rna_exec_summary`, `rna_pipeline_summary`, `rna_report` targets complete without error.
5. **Multiomics passthrough**: verify multiomics `extract_enrichment_df()` can read the pathway results.

---

## Function reference (new functions on this branch)

All in `R/core/09_enrichment.R`, appended after existing code:

| Function | Purpose |
|----------|---------|
| `load_local_pathway_tables()` | Read KEGG/GO tab files into TERM2GENE/TERM2NAME pairs |
| `rank_by_pval_wo_direction()` | -log10(pvalue) ranking |
| `rank_by_pval_with_direction()` | Signed pvalue ranking |
| `rank_by_fc()` | Log2 signed FC ranking (with linearFC recovery) |
| `rank_by_min_pval_any_contrast()` | Cross-contrast min pvalue ranking |
| `build_ranked_gene_lists()` | Builds all 4 ranking variants |
| `run_gsea_local()` | Single GSEA run via clusterProfiler::GSEA() |
| `run_gsea_all()` | Parallelized orchestrator for all GSEA jobs |
| `save_gsea_per_pathway_artifacts()` | Per-pathway gseaplot2 PNGs + core gene CSVs |
| `run_cluster_ora()` | Legacy-faithful cluster-based ORA (enricher + merge_result + simplify) |
| `process_enrichment_table()` | Fold enrichment computation from GeneRatio/BgRatio |
| `build_gene_lists()` | Constructs gene_lists[[method]][[round]] from DE + clustering |

Module helpers in `R/modules/rnaseq/05_mod_pathway.R`:

| Function | Purpose |
|----------|---------|
| `.run_local_enrichment()` | Internal orchestrator for the local enrichment path |
| `.store_ora_result()` | Adds padj/pathway columns to ORA results for downstream compatibility |
