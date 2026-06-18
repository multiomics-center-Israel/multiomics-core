# Enrichment Layer Audit Report

**Date**: 2026-04-13
**Scope**: All enrichment-related code in multiomics-core
**Type**: Code audit + cleanup plan (no implementation)

---

## 1. Enrichment Audit Map

### 1.1 Files

| Layer | File | Lines | Role |
|-------|------|-------|------|
| **Core** | `R/core/09_enrichment.R` | 651 | Shared enrichment: GMT reading, gene set loading, ORA, fGSEA, plots |
| **Core** | `R/core/09b_enrichment.R` | 651 | **EXACT DUPLICATE** of `09_enrichment.R` (MD5 identical) |
| **Core** | `R/core/13_gmt_utils.R` | ~451 | GMT writing, validation, generation, merging |
| **Core** | `R/core/11_annotation.R` | — | `get_organism_info()` used by enrichment |
| **Domain/rnaseq** | `R/domain/rnaseq/07_pathway.R` | 437 | RNA-seq pathway: fGSEA, ORA, plots, volcano data |
| **Domain/proteomics** | `R/domain/proteomics/07b_pathway.R` | 240 | Proteomics-to-pathway bridge; calls shared functions |
| **Domain/metabolomics** | `R/domain/metabolomics/06_enrichment.R` | ~1300 | Independent metabolomics enrichment (QEA, ssGSEA, ORA, GSEA) |
| **Domain/multiomics** | `R/domain/multiomics/07_enrichment.R` | ~1500 | Cross-omics KEGG enrichment, consensus analysis |
| **Domain/multiomics** | `R/domain/multiomics/07b_multigsea_plots.R` | — | MultiGSEA correlation plots |
| **Domain/multiomics** | `R/domain/multiomics/07c_kegg_pathview.R` | — | KEGG pathview overlay visualization |
| **Module/rnaseq** | `R/modules/rnaseq/05_mod_pathway.R` | 140 | RNA-seq pathway module orchestrator |
| **Module/proteomics** | `R/modules/proteomics/05_mod_pathway.R` | 17 | Proteomics pathway module (thin wrapper) |
| **Module/proteomics** | `R/modules/proteomics/05b_mod_pathway.R` | 17 | **EXACT DUPLICATE** of `05_mod_pathway.R` |
| **Module/metabolomics** | `R/modules/metabolomics/05_mod_enrichment.R` | ~200 | Metabolomics enrichment module orchestrator |
| **Module/multiomics** | `R/modules/multiomics/04_mod_enrichment.R` | 88 | Cross-omics enrichment module |

### 1.2 Key Functions and Ownership

| Function | Defined In | Called From | Runtime Winner |
|----------|-----------|-------------|----------------|
| `read_gmt()` | `core/09_enrichment.R:15`, `core/09b_enrichment.R:15` | `core/09_enrichment.R:63`, `core/13_gmt_utils.R:67,319` | rnaseq domain (last loaded, but identical) |
| `load_gene_sets()` | `core/09_enrichment.R:52`, `core/09b_enrichment.R:52` | `modules/rnaseq/05_mod_pathway.R:89`, `domain/proteomics/07b_pathway.R:194` | Core version (no domain override) |
| `run_ora()` | `core/09_enrichment.R:247`, `core/09b_enrichment.R:247` | Internal to `run_pathway_analysis()` | Whichever `run_pathway_analysis()` calls it |
| `run_pathway_analysis()` | `core/09_enrichment.R:428`, `core/09b_enrichment.R:428`, **`domain/rnaseq/07_pathway.R:109`** | `modules/rnaseq/05_mod_pathway.R:110`, `domain/proteomics/07b_pathway.R:215` | **Domain/rnaseq version** (see Section 2) |
| `save_pathway_results()` | `core/09_enrichment.R:549`, `core/09b_enrichment.R:549`, **`domain/rnaseq/07_pathway.R:243`** | `modules/rnaseq/05_mod_pathway.R:125`, `domain/proteomics/07b_pathway.R:226` | **Domain/rnaseq version** |
| `generate_pathway_plots()` | `core/09_enrichment.R:581`, `core/09b_enrichment.R:581`, **`domain/rnaseq/07_pathway.R:279`** | `modules/rnaseq/05_mod_pathway.R:126`, `domain/proteomics/07b_pathway.R:229` | **Domain/rnaseq version** |
| `lookup_go_term_names()` | `core/09_enrichment.R:337`, `core/09b_enrichment.R:337`, `domain/rnaseq/07_pathway.R:14` | Internal to `add_pathway_names()` | Domain/rnaseq version |
| `add_pathway_names()` | `core/09_enrichment.R:374`, `core/09b_enrichment.R:374`, `domain/rnaseq/07_pathway.R:51` | Internal to `run_pathway_analysis()` | Domain/rnaseq version |
| `build_pathway_volcano_data()` | `domain/rnaseq/07_pathway.R:363` | **Nowhere** (defined but never called) | N/A — dead code |
| `read_gmt_list()` | `domain/metabolomics/06_enrichment.R:22` | `domain/metabolomics/06_enrichment.R:343,447,631,789` | Metabolomics-only (parallel to core `read_gmt()`) |
| `run_metabolomics_*()` | `domain/metabolomics/06_enrichment.R` | `modules/metabolomics/05_mod_enrichment.R` | Metabolomics-only (independent system) |
| `build_per_omics_enrichment()` | `domain/multiomics/07_enrichment.R:22` | `modules/multiomics/04_mod_enrichment.R:31` | Multiomics-only |

### 1.3 Sourcing Order (from `_targets.R`)

```
core_files -> service_files -> domain_files -> module_files -> pipeline_files
```

Files are `sort()`-ed alphabetically within each category, then sourced via `tar_source()`.

**Within core**: `09_enrichment.R` loads first, then `09b_enrichment.R` overwrites it (identical, so no effect).

**Core -> Domain**: Domain files load AFTER core. `R/domain/rnaseq/07_pathway.R` re-defines `run_pathway_analysis`, `save_pathway_results`, `generate_pathway_plots`, `lookup_go_term_names`, and `add_pathway_names`, **overriding the core versions for ALL downstream callers** — including proteomics.

---

## 2. Duplicates and Ambiguity

### 2.1 Byte-Identical Duplicate Files

| Original | Duplicate | Verdict |
|----------|-----------|---------|
| `R/core/09_enrichment.R` | `R/core/09b_enrichment.R` | **Accidental**. Files are MD5-identical. `09b` always overwrites `09` at load time (alphabetical order). One is dead weight. |
| `R/modules/proteomics/05_mod_pathway.R` | `R/modules/proteomics/05b_mod_pathway.R` | **Accidental**. Files are byte-identical. `05b` always overwrites `05` at load time. One is dead weight. |

### 2.2 Diverged Duplicate Functions

Five functions are defined in BOTH `R/core/09_enrichment.R` AND `R/domain/rnaseq/07_pathway.R`:

#### `run_pathway_analysis()` — DIVERGED (correctness issue)

- **Core** (`09_enrichment.R:428-541`): Original version. Does not guard against empty fGSEA results. Does not guard against empty ORA results before appending metadata columns.
- **Domain/rnaseq** (`07_pathway.R:109-231`): Improved version. Adds:
  - Empty fGSEA result guard (lines 164-166): `if (nrow(fgsea_df) == 0) { message(...) } else { ... }`
  - Empty ORA result guard (lines 200, 212): `if (nrow(ora_up) > 0) { ... }`
- **Divergence type**: The domain version has **bug fixes** the core version lacks.
- **Runtime winner**: Domain/rnaseq version (loads later).
- **Impact**: The core version's bugs are masked at runtime but would resurface if anyone removes the domain version or if a future refactor changes load order.

#### `save_pathway_results()` — IDENTICAL

- Core (`09_enrichment.R:549-573`) and Domain (`07_pathway.R:243-267`): Same logic.
- **Divergence type**: None. Pure copy.

#### `generate_pathway_plots()` — IDENTICAL

- Core (`09_enrichment.R:581-651`) and Domain (`07_pathway.R:279-349`): Same logic.
- **Divergence type**: None. Pure copy.

#### `lookup_go_term_names()` — IDENTICAL

- Core (`09_enrichment.R:337-367`) and Domain (`07_pathway.R:14-43`): Same logic.

#### `add_pathway_names()` — IDENTICAL

- Core (`09_enrichment.R:374-414`) and Domain (`07_pathway.R:51-91`): Same logic.

### 2.3 Shadowing / Override Behavior

The core file `09_enrichment.R` contains a comment at line 328-330:
```r
# Moved from R/domain/rnaseq/07_pathway.R — generic DE table interface
```

This confirms the **history**: these functions were copied from the rnaseq domain file to core to share them with proteomics. However, **the rnaseq domain file was never cleaned up** — the original definitions were left in place. The domain version then continued to evolve independently (adding empty-result guards), while the core copy stagnated.

### 2.4 Runtime Ambiguity

**Critical**: `R/domain/proteomics/07b_pathway.R` line 4 documents:
```r
#' (run_pathway_analysis from R/core/09_enrichment.R).
```

This comment is **misleading**. At runtime, proteomics gets the **domain/rnaseq** version of `run_pathway_analysis()`, NOT the core version, because the domain/rnaseq file loads after core and overrides the definition globally.

This is currently benign (the domain version is more correct), but:
- The documented ownership is wrong
- Any future change to the rnaseq domain version could unexpectedly change proteomics behavior
- If someone "fixes" this by removing the domain copy (trusting the comment), proteomics would regress to the buggier core version

### 2.5 Parallel GMT Parsers

| Function | File | Signature |
|----------|------|-----------|
| `read_gmt()` | `core/09_enrichment.R:15` | `read_gmt(gmt_file)` — returns named list, descriptions as attribute |
| `read_gmt_list()` | `domain/metabolomics/06_enrichment.R:22` | `read_gmt_list(gmt_file, include_descriptions=FALSE)` — returns named list or `list(sets, descriptions)` |

These are **not duplicates** — they have different interfaces and serve different callers. But they do the same fundamental job (parse GMT files) with different return contracts. The metabolomics version handles missing files gracefully (returns empty); the core version `stop()`s.

---

## 3. Dead or Suspicious Code

### 3.1 Confirmed Dead Code

| Item | Location | Evidence |
|------|----------|----------|
| `build_pathway_volcano_data()` | `domain/rnaseq/07_pathway.R:363-437` | Defined but **never called** anywhere in the codebase (.R or .Rmd). No callers found via grep. |
| `R/core/09b_enrichment.R` (entire file) | `R/core/09b_enrichment.R` | Byte-identical duplicate of `09_enrichment.R`. All its definitions are immediately overwritten by domain files at load time. Provides zero unique value. |
| `R/modules/proteomics/05b_mod_pathway.R` (entire file) | `R/modules/proteomics/05b_mod_pathway.R` | Byte-identical duplicate of `05_mod_pathway.R`. |

### 3.2 Likely Legacy / Stale Code

| Item | Location | Evidence |
|------|----------|----------|
| Core versions of `run_pathway_analysis`, `save_pathway_results`, `generate_pathway_plots`, `lookup_go_term_names`, `add_pathway_names` | `core/09_enrichment.R:337-651` | These are **never the runtime winner** — always overridden by domain/rnaseq. The core `run_pathway_analysis` has known missing guards compared to the domain version. These ~315 lines of core code are dead at runtime. |

### 3.3 Assessment

The `"b"` suffix files (`09b_enrichment.R`, `05b_mod_pathway.R`, `07b_pathway.R`) appear to follow a naming pattern. However:
- `09b` and `05b_mod` are pure accidental duplicates (byte-identical)
- `07b_pathway.R` (proteomics) is a legitimate separate file (different content — proteomics bridge logic)

The `"b"` suffix is being used inconsistently: sometimes for "variant/extension" (proteomics `07b`), sometimes as an accidental copy.

---

## 4. Interface Inconsistencies

### 4.1 DE Table Schema

The shared enrichment functions expect this standard schema:

```
FeatureID | log2FoldChange | pvalue | padj | stat
```

- **RNA-seq**: DE tables already use this schema natively (from DESeq2).
- **Proteomics**: `extract_de_table_for_pathway()` (in `07b_pathway.R:21-64`) converts from the proteomics-native wide format (`padj.imputs.<cn>`, `linearFC.imputs.<cn>`) to the standard schema. The `stat` column is synthesized as `sign(log2FC) * -log10(pvalue)`.
- **Metabolomics**: Does NOT use the shared enrichment at all — has its own parallel system.
- **Multiomics**: Does NOT use the shared enrichment — has its own KEGG-focused system.

**No inconsistency here** — the adapter pattern in proteomics is correct.

### 4.2 Gene Set Format

- `load_gene_sets()` returns: `list(GO = list(...), KEGG = list(...), custom = list(...))`
- Each value is a named list of character vectors (gene IDs).
- `read_gmt()` returns: named list with `attr(, "descriptions")`.
- `read_gmt_list()` returns: either a plain named list or `list(sets=, descriptions=)`.

**Inconsistency**: The two GMT parsers return descriptions differently (attribute vs list element). This doesn't cause bugs today because metabolomics never calls core `read_gmt()`, but it means any future unification would need to harmonize the return format.

### 4.3 Output Directory Layout

Both RNA-seq and proteomics modules follow the same convention:
```
{out_dir}/Enrichment/           — CSV result files
{out_dir}/Enrichment/plots/     — PNG plot files
{out_dir}/Enrichment/gene_annotation.csv
```

Metabolomics uses a different convention via `create_legacy_output_dirs()`:
```
{out_dir}/enrichment/           — plots
{out_dir}/datasets/             — result TSV files
```

**Inconsistency**: Different capitalization (`Enrichment` vs `enrichment`) and different file format (CSV vs TSV) and different directory layout between the genomics and metabolomics paths.

### 4.4 `run_ora()` Hard-Coded Thresholds

In `run_pathway_analysis()` (both versions), the ORA branch hard-codes DE significance thresholds:
```r
de_cfg_padj <- 0.05
de_cfg_lfc  <- log2(1.5)
```
These are NOT pulled from config and cannot be overridden by the user. This could cause inconsistency if the pipeline's main DE thresholds differ from these values.

---

## 5. Cleanup Priorities

### 5.1 Must Clean Now (Before Any New Integration Work)

| # | Issue | Risk If Deferred | Effort |
|---|-------|-------------------|--------|
| 1 | **Delete `R/core/09b_enrichment.R`** | Byte-identical duplicate. Zero risk to delete. Leaves a confusing ghost file for anyone reading the codebase. Could mislead someone into editing the wrong file. | 1 min |
| 2 | **Delete `R/modules/proteomics/05b_mod_pathway.R`** | Same as above. Byte-identical duplicate of `05_mod_pathway.R`. | 1 min |
| 3 | **Consolidate the 5 duplicated functions**: Move the improved domain/rnaseq versions into core, remove them from `07_pathway.R` | The core version of `run_pathway_analysis` has missing empty-result guards. If load order changes or someone removes the domain copy trusting the "moved to core" comment, proteomics would get the buggier core version. Also blocks any clean integration of new enrichment methods. | 30 min |
| 4 | **Fix the misleading comment** in `07b_pathway.R:4` ("run_pathway_analysis from R/core/09_enrichment.R") | Documents wrong ownership. Could mislead someone into breaking changes. | 1 min |

### 5.2 Safe to Defer

| # | Issue | Why Deferrable |
|---|-------|----------------|
| 5 | Remove or repurpose `build_pathway_volcano_data()` | Dead code, but harmless. Not blocking anything. Can be removed in a future cleanup pass or connected to a UI feature. |
| 6 | Harmonize `read_gmt()` vs `read_gmt_list()` return format | Metabolomics enrichment is independent and works. Only matters if someone tries to unify the GMT parsing path. |
| 7 | Make ORA thresholds configurable | Cosmetic/flexibility issue, not a correctness issue. |
| 8 | Harmonize output directory naming (`Enrichment/` vs `enrichment/`) | Only a problem if cross-omics tools try to find results from single-omics runs by path. |

### 5.3 Leave As-Is

| # | Item | Reason |
|---|------|--------|
| 9 | Metabolomics having its own independent enrichment system | Metabolomics enrichment (QEA, ssGSEA, compound-based ORA/GSEA) is fundamentally different from gene-based enrichment. The parallel system is justified. |
| 10 | Multiomics having its own KEGG enrichment system | Cross-omics enrichment has different inputs (multiple DE tables, compound mappings, consensus logic). The separate system is justified. |
| 11 | Proteomics bridge pattern (`extract_de_table_for_pathway` + `map_proteins_to_gene_symbols`) | Clean adapter pattern. Works correctly. |

---

## 6. Smallest Safe Cleanup Plan

### Step 1: Delete byte-identical duplicate files

- Delete `R/core/09b_enrichment.R`
- Delete `R/modules/proteomics/05b_mod_pathway.R`
- **Verification**: Run `tar_make(callr_function = NULL)` or `tar_manifest()` to confirm no target references these files by name. Since they are `tar_source()`-ed by directory glob, removing them simply removes the redundant load.
- **Risk**: Zero. The remaining copies define the same functions.

### Step 2: Promote domain/rnaseq function improvements to core

For these 5 functions currently duplicated between `core/09_enrichment.R` and `domain/rnaseq/07_pathway.R`:

- `lookup_go_term_names()`
- `add_pathway_names()`
- `run_pathway_analysis()`
- `save_pathway_results()`
- `generate_pathway_plots()`

**Action**:
1. Replace the core versions in `R/core/09_enrichment.R` with the domain/rnaseq versions (which have the improved empty-result guards).
2. Remove these 5 function definitions from `R/domain/rnaseq/07_pathway.R`.
3. Keep `build_pathway_volcano_data()` in `07_pathway.R` for now (it's rnaseq-specific, even if unused).

**Verification**:
- Diff the core and domain versions before replacing. Confirm the only differences are the empty-result guards in `run_pathway_analysis()`.
- Run the RNA-seq and proteomics test pipelines to confirm no behavior change.

**Risk**: Low. The domain version is already the runtime winner. This just makes the code match the runtime reality.

### Step 3: Fix the misleading comment in proteomics

**File**: `R/domain/proteomics/07b_pathway.R`
**Line 4**: Change from:
```r
#' (run_pathway_analysis from R/core/09_enrichment.R).
```
To:
```r
#' Bridges proteomics DE results to the shared pathway analysis
#' functions defined in R/core/09_enrichment.R (run_pathway_analysis,
#' save_pathway_results, generate_pathway_plots).
```

### What NOT to Change Yet

- Do NOT merge metabolomics enrichment into the core pathway system
- Do NOT merge multiomics enrichment into the core pathway system
- Do NOT add new enrichment methods until steps 1-3 are done
- Do NOT rename files (e.g., `07b_pathway.R` -> `07_pathway.R`) in this cleanup — the proteomics `07b` naming is intentional (slot `07` may be reserved or used differently)
- Do NOT refactor the output directory structure

---

## 7. Risks

### If cleanup is done carelessly

| Risk | How to Mitigate |
|------|-----------------|
| Deleting `09b_enrichment.R` could break if any code references it by name | Grep for `09b_enrichment` in the entire repo. Confirmed: no explicit references to the filename. Only loaded via directory glob. |
| Moving functions from `07_pathway.R` to core could change behavior if the versions have silently diverged further | Do a line-by-line diff of each function pair immediately before editing. The audit shows `save_pathway_results` and `generate_pathway_plots` are identical; `run_pathway_analysis` diverges only in empty-result guards. |
| Proteomics pathway behavior could change if `run_pathway_analysis` is modified during consolidation | The domain/rnaseq version is already what proteomics gets at runtime. Promoting it to core preserves exact current behavior. |
| Removing functions from `07_pathway.R` could break if any rnaseq code does `source("R/domain/rnaseq/07_pathway.R")` directly | Grep for direct `source()` calls to this file. Confirmed: none found. All loading is via `tar_source()` directory glob. |

### What must be verified before implementation

1. Run `md5sum` on `09_enrichment.R` vs `09b_enrichment.R` immediately before deletion (done: confirmed identical).
2. Run `md5sum` on `05_mod_pathway.R` vs `05b_mod_pathway.R` immediately before deletion (done: confirmed identical).
3. Diff `run_pathway_analysis()` between core and domain/rnaseq to confirm only the empty-result guards differ (done: confirmed).
4. After consolidation, run at minimum: `tar_make()` or the RNA-seq and proteomics pathway targets in isolation to verify no regression.
5. Verify that `build_pathway_volcano_data()` has no callers before deciding to leave or remove it.
