# Enrichment Layer Migration — Implementation Spec

Status: **frozen design note**
Scope: RNA-seq enrichment (proteomics adaptation deferred)
Reference: legacy enrichment Rmd + helper functions (provided externally)

---

## 1. Non-Negotiable Behavioral Requirements

### ORA

- The unit of analysis is the **cluster**, not the contrast. Each cluster's gene set is tested independently against the pathway database.
- ORA uses `clusterProfiler::enricher()` with: `minGSSize=0`, `maxGSSize=10000`, `qvalueCutoff=1`. Background universe is the implicit TERM2GENE universe (enricher default when no `universe` argument is passed).
- Only clusters with at least one term passing `p.adjust < pvalueCutoff` are retained before merging.
- Per-cluster results are merged into a single combined result with a `Cluster` column.
- For GO types, semantic similarity-based simplify is applied (`cutoff=0.7, by="p.adjust", select_fun=min`), producing a separate result set.
- Output tables include fold enrichment: `(in_cluster_in_term / in_cluster) / (in_term / in_genome)`, with `GeneRatio` and `BgRatio` expanded into numeric columns.

### GSEA

- Three per-contrast ranking methods, each producing a separate named vector:
  - `pval_wo_direction`: `-log10(pvalue)`, always positive
  - `pval_with_direction`: `sign(linearFC) * -log10(pvalue)`
  - `fc`: `log2(ifelse(fc > 0, fc, -1/fc))`, rounded via `signif(digits=4)`
- One cross-contrast method: `any_contrast` — row-wise `min(pvalue)` across all contrasts, then `-log10`
- GSEA uses `minGSSize=4` and `maxGSSize=length(unique(TERM2GENE[,2]))` (total unique genes in the pathway database, not a fixed cap).
- All ranking methods use raw `pvalue`, not `padj`.
- The `fc` ranking method assumes **linear fold change** input. Current pipeline provides `log2FoldChange` — conversion is required.

### Data loading

- Enrichment consumes only local precomputed tab-separated tables:
  - KEGG: `KEGG_pathway2gene.tab`, `KEGG_pathway2name.tab`
  - GO: `GO2gene_{BP,MF,CC}.tab`, `GO2name_{BP,MF,CC}.tab`
- Each file is two columns (term ID, gene ID or term name). No header required.
- No online resources: no KEGG REST, no biomaRt, no OrgDb-based pathway loading, no GMT as primary format.

### Result organization

- ORA results stored per database type, keyed by clustering method/round.
- GSEA results stored per database type, keyed by ranking method and contrast.
- Both ORA and GSEA produce per-analysis CSV files with the enrichment table.

---

## 2. Implementation Choices

### compareClusterResult

`compareClusterResult` is the **preferred implementation vehicle**, not a strict requirement. The actual requirements are: (a) GO simplify works, (b) cluster-aware dotplot renders, (c) merged result table has a `Cluster` column. Using `merge_result()` to produce a `compareClusterResult` is the lowest-risk path to all three. Do not treat the S4 class as sacred — if the three behaviors are satisfied another way, that is acceptable.

### GSEA parity

GSEA behavioral parity is **expected, pending implementation-time verification** of three specific items:

1. The `fc` ranking method requires converting `log2FoldChange` back to linear FC before applying the legacy transform. The math must be validated.
2. The `any_contrast` ranking must be reconstructed from per-contrast DE tables (current pipeline format), not from a wide `stats_df`. Must produce identical gene rankings.
3. The `signif(digits=4)` rounding in the `fc` method must not be dropped.

### Clustering disabled policy

**Warn and skip ORA.** When `enrichment.annotation_dir` is configured but clustering results are NULL (disabled or empty): run GSEA, skip ORA, emit message: `"Cluster-based ORA requires clustering results. Enable clustering in config to run ORA. GSEA will proceed."`

Do not fall back to per-contrast ORA as a substitute. Per-contrast ORA is a different analysis with different gene sets, backgrounds, and parameters. It does not approximate the legacy cluster-based behavior.

---

## 3. Acceptance Criteria

The migration is behavior-preserving when ALL of the following are true:

1. **Offline only.** Pipeline completes enrichment with no network access. No KEGG REST, biomaRt, or OrgDb-based pathway loading.
2. **Cluster-based ORA.** Output CSVs contain `Cluster`, `Fold_enrichment`, `in_cluster_in_term`, `in_cluster`, `in_term`, `in_genome` columns. ORA runs per cluster with `minGSSize=0, maxGSSize=10000, qvalueCutoff=1`.
3. **GO simplify.** Both unsimplified and simplified CSVs exist for GO_BP, GO_MF, GO_CC.
4. **GSEA coverage.** Output CSVs exist for each combination of {pval_wo_direction, pval_with_direction, fc} x {contrast} x {database}, plus any_contrast x {database}.
5. **GSEA parameters.** `minGSSize=4`, `maxGSSize=length(unique(TERM2GENE[,2]))`.
6. **Downstream intact.** `rna_exec_summary`, `rna_pipeline_summary`, `rna_report`, and multiomics enrichment targets complete without error. `collect_pipeline_stats()` finds data.frames with `padj` columns in `pathway_res$pathway_results`.
7. **Clustering disabled graceful.** With `clustering.enabled: false`, GSEA produces results, ORA is skipped with explicit warning, pipeline completes.

If criteria 1, 4, 5, 6, 7 are met but criteria 2, 3 are not (because clustering is disabled), the state is **"partial — GSEA complete, ORA pending clustering enablement"**.

---

## 4. Known Implementation Risks

| Risk | Impact | Mitigation |
|------|--------|------------|
| `linearFC` vs `log2FoldChange` mismatch in `fc` ranking | Silently wrong gene rankings | Explicit conversion: `2^log2FC` to recover linear FC before applying legacy transform. Validate with a test case. |
| `clusterProfiler::enricher()` or `merge_result()` not available in environment | Enrichment fails at runtime | Check availability at module entry. Emit actionable install message. |
| `simplify()` requires GOSemSim + OrgDb for semantic similarity | GO simplify fails for organisms without OrgDb installed | Make simplify conditional. Warn if OrgDb not available. OrgDb is a local package, not an online resource — document that it must be pre-installed. |
| Gene IDs in local TERM2GENE tables don't match pipeline FeatureIDs | Zero overlap, empty enrichment results | Validation in loader: compute overlap between TERM2GENE gene column and pipeline FeatureIDs. Warn if <5%. |
| `any_contrast` ranking reconstructed from per-contrast tables may differ from legacy wide-format `apply(..., min)` due to gene set differences across tables | Slightly different gene rankings | Align genes across all contrasts before taking row-wise min. Verify identical output on test data. |
| Downstream consumers expect `pathway_res$pathway_results` with specific nesting | Report/summary modules break | Return structure must include `$pathway_results` containing nested data.frames with `padj` columns. Confirmed compatible with `collect_pipeline_stats()` and `extract_enrichment_df()`. |
| Dotplot generation via `enrichplot::dotplot()` on `compareClusterResult` may fail if enrichplot version mismatches | Plot generation fails, not analysis | Wrap in tryCatch. Fall back to basic ggplot2 dotplot if enrichplot fails. |

---

## 5. Minimal Implementation Order

### Phase 1 — Local loading + GSEA (no DAG change)

| Step | File | Change | Confidence |
|------|------|--------|------------|
| 1 | `R/core/09_enrichment.R` | ADD `load_local_pathway_tables()` — reads KEGG + GO TSVs into TERM2GENE/TERM2NAME pairs, validates structure and overlap | High |
| 2 | `R/core/09_enrichment.R` | ADD four ranked gene list builder functions, adapted to accept per-contrast DE tables. Handle linearFC conversion for `fc` method. Handle cross-table alignment for `any_contrast`. | High |
| 3 | `R/core/09_enrichment.R` | ADD `run_gsea_local()` — calls `clusterProfiler::GSEA()` with local TERM2GENE/TERM2NAME, `minGSSize=4`, `maxGSSize=n_unique_genes` | Medium |
| 4 | `R/core/09_enrichment.R` | ADD `process_enrichment_table()` — port of legacy fold enrichment computation from `process_clusterprofiler_results_table()` | High |
| 5 | `config/templates/rna_config.yaml` | ADD `enrichment: annotation_dir: ""` under `modes.rna` | High |
| 6 | `R/modules/rnaseq/05_mod_pathway.R` | MODIFY `mod_rnaseq_pathway()` — if `annotation_dir` is set: load local tables, build ranked lists, run GSEA per method x contrast x db, save results. If not set: existing online path unchanged. Emit clustering-disabled warning for ORA. | Medium |

**Phase 1 delivers:** offline data loading, GSEA with all four ranking methods, fold enrichment tables, backward compatibility.
**Phase 1 does not deliver:** cluster-based ORA, GO simplify, cluster-aware dotplots.

### Phase 2 — Cluster-based ORA (one-line DAG change)

| Step | File | Change | Confidence |
|------|------|--------|------------|
| 7 | `R/core/09_enrichment.R` | ADD `run_cluster_ora()` — per-cluster `enricher()`, merge, simplify, process, dotplot. Returns 4-element list matching legacy structure. | Medium |
| 8 | `R/pipeline/rnaseq/00_pipe_rnaseq.R` | MODIFY `rna_pathway_res` target: add `clustering_res = rna_clustering_obj` argument | High |
| 9 | `R/modules/rnaseq/05_mod_pathway.R` | MODIFY signature: add `clustering_res = NULL`. When clusters available, run cluster-based ORA per db. | Medium |

**Phase 2 delivers:** cluster-based ORA, GO simplify, cluster-aware dotplots, complete behavioral parity.

### Phase 3 — Deferred

- Proteomics pathway module adaptation
- Ridgeplot porting
- Shared-genes heatmap porting
- Output directory restructuring (pending report rendering analysis)
