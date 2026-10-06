# Claude Code Fix Prompts — Multi-omics Core

**Project:** `multiomics-core` pipeline maintenance and bug fixes.
**Updated:** 2026-03-02

> **Note:** Prompts 1 (Report Template) and 2 (Pipeline Tests) from the original version of this
> document have been **removed** — those tasks are already complete:
> - `R/domain/multiomics/report_template_multiomics.Rmd` exists (823 lines, fully implemented)
> - `tests/testthat/test-pipeline-orchestration.R` (16 tests) and
>   `tests/testthat/test-multiomics-harmonization.R` (11 tests) are in place

---

## Prompt: Proteomics Contrast Column Name Mismatch (FIXED)

### Root Cause (Identified & Patched)

When `tar_make()` builds the proteomics DE summary, contrast names with spaces
(e.g., `"1.56ppm vs. 0ppm"`) are stripped of ALL spaces by
`gsub(" ", "", cn)` in two producer functions inside
`R/domain/proteomics/05_de_summary.R`:

| Function | Line | Code |
|---|---|---|
| `summarize_limma_mult_imputation()` | 59 | `contrast_print <- gsub(" ", "", cn)` |
| `load_precomputed_proteomics_de()` | 447 | `contrast_print <- gsub(" ", "", ctr)` |

This creates columns like `linearFC.imputs.1.56ppmvs.0ppm` (no spaces).

However, downstream **consumers** reconstruct column names from the original
config contrast name *without* stripping spaces, producing lookups like
`linearFC.imputs.1.56ppm vs. 0ppm` — which don't exist. The mismatch causes
`stopifnot()` and `stop()` failures during Excel export and pathway analysis.

### Affected Consumer Locations

| File | Function | Issue |
|---|---|---|
| `R/core/05_export_excel.R` | `get_contrast_cols()` | Took raw contrast name, built columns with spaces |
| `R/core/05_export_excel.R` | `fill_manual_cutoffs_formulas_legacy()` | Did not pass `mode` to `get_contrast_cols()` |
| `R/domain/proteomics/07b_pathway.R` | `extract_de_table_for_pathway()` | Built column names with spaces |

### Fix Applied

1. **Created `normalize_contrast_name()`** in `R/core/01_io.R` — a shared utility
   that centralizes `gsub(" ", "", x)` so producers and consumers use the same logic.

2. **Updated producers** in `R/domain/proteomics/05_de_summary.R` (lines 59, 447)
   to call `normalize_contrast_name()` instead of inline `gsub()`.

3. **Updated consumers**:
   - `get_contrast_cols()` now normalizes the contrast name before building column names.
   - `fill_manual_cutoffs_formulas_legacy()` now passes `mode` to `get_contrast_cols()`.
   - `extract_de_table_for_pathway()` now normalizes `contrast_name` before building column names.

### Validation Prompt (use after the fix to confirm)

```text
Run tar_make() for the GT15 config and verify that the proteomics export
step completes without "Missing required columns" errors. Check that:
1. summary_df column names match what get_contrast_cols() constructs
2. The Excel export (build_final_results_generic) finds all expected columns
3. Pathway analysis (extract_de_table_for_pathway) finds padj/FC columns
```
