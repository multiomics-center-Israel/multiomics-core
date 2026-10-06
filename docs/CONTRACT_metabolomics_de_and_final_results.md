# Contract: Metabolomics DE Column Naming & Final Results Schema

> Status: **IMPLEMENTED, with deviations** — the rename below is what the
> pipeline writes today; §0 lists where the code differs from this design.
> Designed: 2026-03-11. Checked against the code: 2026-10-06.

---

## 0. Implementation status

The column rename (§2), the `get_contrast_cols()` branch (§3), the code
changes in §4 and the final-results builder (§5) are in the code. Lipidomics
uses the same `build_de_summary()` and the same `get_contrast_cols()` branch.

Where the code **differs** from this design (recorded, not resolved — the code
has not been changed to match the document, nor the document to hide the gap):

- **`pass_any_contrast` and `n_pass_contrasts` depend on the DE path.** The
  pre-computed loader calls `add_pass_any_contrast()`, so it writes
  `pass_any_contrast` as `1`/`NA` and adds `n_pass_contrasts`, as designed.
  The computed path, `build_de_summary()`, computes `pass_any_contrast` inline
  as `1`/`0` and writes no `n_pass_contrasts`. The final results copy the
  value through unchanged, so their `pass_any_contrast` is `1`/`NA` or `1`/`0`
  depending on which path ran.
- **The log2 fold change is round-tripped through a rounded value.**
  `linearFC.<cn>` is stored as `signif(..., 3)`, and `extract_contrast_table()`
  back-computes `logFC` from it for volcano/MA plots and the report. The
  "no round-trip needed" row in §7 does not hold for those consumers.
- **There is no `config_padj_cutoff` key.** The metabolomics Shiny payload
  carries `padj_cutoff`, `log_fc_cutoff` and `fc_cutoff`
  (see `SHINY_PAYLOAD_HANDOFF.md`). `padj_cutoff` there is read as
  `de$padj_cutoff`, falling back to `de$p_cutoff`; the DE step reads
  `de$p_cutoff` only.
- **Additions not in the design:** `original_id` (when present in
  `row_data`) is placed directly after `feature_id`; the builder takes
  `cv_contrasts_df` and `config` for the CV columns; `row_data` defaults to
  `pre$row_data`.
- **No per-group mean columns.** Unlike RNA-seq and proteomics, the
  metabolomics final results carry no `Mean.<group>` columns. The design never
  promised them; this is noted so nobody assumes parity.

---

## 1. Problem Statement

*Historical: this section describes the naming before the rename. Since then
RNA-seq's pass column has become `<cn>_pass` and its summary has gained
`log2FC.<cn>` columns; see `CONTRACT_rnaseq_final_results.md`.*

Before this contract, the metabolomics DE pipeline (`build_de_summary()` in
`03_differential.R`) used a column naming convention that differed from both
proteomics and RNA-seq:

| Aspect | Proteomics | RNA-seq | Metabolomics (before) |
|--------|-----------|---------|----------------------|
| Separator | `.` | `.` | `_` |
| Infix | `.imputs.` | (none) | (none) |
| FC column | `linearFC.imputs.<cn>` | `linearFC.<cn>` | `logFC_<cn>` |
| P-value | `pvalue.imputs.<cn>` | `pvalue.<cn>` | `P.Value_<cn>` |
| Adj p-value | `padj.imputs.<cn>` | `padj.<cn>` | `adj.P.Val_<cn>` |
| Pass flag | `pass.imputs.<cn>` | `pass.<cn>` | `pass_<cn>` |
| ID column | `FeatureID` (configurable) | `FeatureID` | `feature_id` |

This forces every downstream consumer (Excel export, Shiny payload, `get_contrast_cols()`)
to branch on mode. The `get_contrast_cols()` function already has `mode` routing for
proteomics vs RNA-seq; metabolomics is not yet wired in.

---

## 2. DE Column Mapping Table

### 2a. Summary: current metabolomics → proposed aligned naming

| # | Current column | Proposed column | Change type | Notes |
|---|---------------|----------------|-------------|-------|
| 1 | `feature_id` | `feature_id` | **Keep** | Already canonical. Proteomics uses configurable `de_table.id_col`; metabolomics keeps `feature_id` as default. |
| 2 | `logFC_<cn>` | `linearFC.<cn>` | **Rename + transform** | `sign(logFC) * 2^abs(logFC)`. Matches proteomics signed linear FC. |
| 3 | `AveExpr_<cn>` | `AveExpr.<cn>` | **Rename separator** | `_` → `.`. Metabolomics-specific extra; proteomics doesn't store this in summary_df. |
| 4 | `P.Value_<cn>` | `pvalue.<cn>` | **Rename** | Drop limma naming, adopt pipeline convention. |
| 5 | `adj.P.Val_<cn>` | `padj.<cn>` | **Rename** | Drop limma naming, adopt pipeline convention. |
| 6 | `pass_<cn>` | `pass.<cn>` | **Rename separator** | `_` → `.` for consistency. Value semantics unchanged: `1` = passed, `0` = did not pass. |
| 7 | *(missing)* | `upDown.<cn>` | **Add (computed)** | `"up"` / `"down"` / `""` based on `pass == 1` + sign of `linearFC`. Required by Excel `fill_manual_cutoffs_formulas_legacy()`. |
| 8 | *(missing)* | `manual_cutoffs.<cn>` | **Add (placeholder)** | `NA` — Excel formula column placeholder. |
| 9 | `pass_any_contrast` | `pass_any_contrast` | **Keep** | Already shared across all modes. |
| 10 | *(missing)* | `n_pass_contrasts` | **Add (computed)** | Integer count. Use shared `add_pass_any_contrast()` from proteomics (already parameterized with `pass_prefix`). **Deviation:** only the pre-computed loader does this; see §0. |

### 2b. Detailed classification

#### Columns that can be renamed directly (no data transformation)

| Current | Aligned | Regex migration |
|---------|---------|-----------------|
| `P.Value_<cn>` | `pvalue.<cn>` | `s/^P\\.Value_/pvalue./` |
| `adj.P.Val_<cn>` | `padj.<cn>` | `s/^adj\\.P\\.Val_/padj./` |
| `pass_<cn>` | `pass.<cn>` | `s/^pass_/pass./` (but careful: must not match `pass_any_contrast`) |
| `AveExpr_<cn>` | `AveExpr.<cn>` | `s/^AveExpr_/AveExpr./` |
| `pass_any_contrast` | `pass_any_contrast` | No change |

#### Columns that require reshaping / computation

| Current | Aligned | Transformation |
|---------|---------|---------------|
| `logFC_<cn>` | `linearFC.<cn>` | `lfc <- logFC; ifelse(lfc >= 0, 2^lfc, -(2^abs(lfc)))` — produces signed linear FC identical to proteomics `linearFC.imputs` semantics |
| *(missing)* | `upDown.<cn>` | `ifelse(pass == 1, ifelse(linearFC > 0, "up", "down"), "")` |
| *(missing)* | `manual_cutoffs.<cn>` | Static `NA` — populated by Excel formula at write time |
| *(missing)* | `n_pass_contrasts` | `rowSums(pass_mat == 1, na.rm = TRUE)` via `add_pass_any_contrast()` |

#### Columns that are metabolomics-specific (remain as additional fields)

| Column | Rationale |
|--------|-----------|
| `AveExpr.<cn>` | Needed for MA plots. Proteomics computes this inside limma but doesn't carry it in summary_df. Keep as bonus column. |
| `feature_id` (vs `FeatureID`) | Metabolomics convention. The `feature_id_col` parameter in `build_final_results_generic()` already handles this — no hardcoding needed. |

---

## 3. `get_contrast_cols()` extension

The function in `R/core/05_export_excel.R` currently has two branches: `"proteomics"` (with `.imputs.`) and `"rna"` (without). Metabolomics uses the **same pattern as RNA-seq** (no imputation infix):

```r
get_contrast_cols <- function(contrast, mode = "proteomics") {
    # ... existing normalization ...

    if (mode == "metabolomics") {
        # Identical structure to RNA — no .imputs. infix
        list(
            fc     = paste0("linearFC.", contrast),
            p      = paste0("pvalue.", contrast),
            padj   = paste0("padj.", contrast),
            pass   = paste0("pass.", contrast),
            updown = paste0("upDown.", contrast),
            manual = paste0("manual_cutoffs.", contrast)
        )
    } else if (mode == "rna") {
        # ... existing RNA branch (unchanged) ...
    } else {
        # ... existing proteomics branch (unchanged) ...
    }
}
```

This is the **only** change needed in `core/05_export_excel.R`.

---

## 4. Files that need updating (implementation scope)

### 4a. Primary: `build_de_summary()` rewrite

**File:** `R/domain/metabolomics/03_differential.R`, function `build_de_summary()`

This is the **single origin** of metabolomics summary_df column names. Change from:

```r
# CURRENT
summary_df[[paste0("logFC_", ctr)]]     <- tbl$logFC
summary_df[[paste0("AveExpr_", ctr)]]   <- tbl$AveExpr
summary_df[[paste0("P.Value_", ctr)]]   <- tbl$P.Value
summary_df[[paste0("adj.P.Val_", ctr)]] <- tbl$adj.P.Val
summary_df[[paste0("pass_", ctr)]]      <- pass
```

To:

```r
# ALIGNED
lfc <- tbl$logFC
linear_fc_signed <- ifelse(lfc >= 0, 2^lfc, -(2^abs(lfc)))

summary_df[[paste0("linearFC.", ctr)]] <- signif(linear_fc_signed, 3)
summary_df[[paste0("AveExpr.", ctr)]]  <- tbl$AveExpr
summary_df[[paste0("pvalue.", ctr)]]   <- tbl$P.Value
summary_df[[paste0("padj.", ctr)]]     <- tbl$adj.P.Val
summary_df[[paste0("pass.", ctr)]]     <- pass
```

Then replace the `pass_any_contrast` block at the end with a call to the shared
`add_pass_any_contrast(summary_df, pass_prefix = "^pass\\.")`.

### 4b. Primary: `load_precomputed_metabolomics_de()` — same file

Same renaming in the precomputed loader loop in `03_differential.R`.

### 4c. Secondary: `extract_contrast_table()` — same file

Update column lookups from `logFC_`, `P.Value_`, `adj.P.Val_` to `linearFC.`, `pvalue.`, `padj.`.
Note: returns `logFC` to callers — must back-compute: `sign(linearFC) * log2(abs(linearFC))`.

### 4d. Secondary: `build_de_summary_counts_metabolomics()` — `07_shiny_export.R`

Currently greps `^pass_` and tries `logFC_` FC patterns. Update to `^pass\\.` and `linearFC.`.

### 4e. Core: `get_contrast_cols()` — `R/core/05_export_excel.R`

Add `"metabolomics"` branch (Section 3 above).

### 4f. No change needed

| File | Why no change |
|------|---------------|
| `build_final_results_generic()` | Already mode-aware via `get_contrast_cols(cn, mode)` |
| `fill_manual_cutoffs_formulas_legacy()` | Already uses `get_contrast_cols()` |
| `write_final_results_excels_legacy_generic()` | Already generic |
| `add_pass_any_contrast()` | Already parameterized with `pass_prefix` |

---

## 5. Proposed `Final_results_metabolomics.tsv` Schema

### 5a. Column layout (in order)

```
┌─────────────┬──────────────────────┬───────────────┬────────────┬─────────────────────┬──────────┐
│  Group 1    │      Group 2         │   Group 3     │  Group 3b  │      Group 4        │ Group 5  │
│  Feature ID │   Annotations        │  Expression   │  Group CV  │  Per-contrast DE    │ Aggregate│
├─────────────┼──────────────────────┼───────────────┼────────────┼─────────────────────┼──────────┤
│ feature_id  │ <annot_col_1>        │ <sample_1>    │ CV.<grp1>  │ linearFC.<cn1>      │ pass_any │
│             │ <annot_col_2>        │ <sample_2>    │ CV.<grp2>  │ pvalue.<cn1>        │ _contrast│
│             │ ...                  │ ...           │ ...        │ padj.<cn1>          │          │
│             │                      │ <sample_N>    │            │ upDown.<cn1>        │          │
│             │                      │               │            │ manual_cutoffs.<cn1>│          │
│             │                      │               │            │ linearFC.<cn2>      │          │
│             │                      │               │            │ ...                 │          │
└─────────────┴──────────────────────┴───────────────┴────────────┴─────────────────────┴──────────┘
```

> **Group CV (`CV.<group>`)** — optional, gated by `modes.<omics>.excel.group_cv`
> (default `true`). One column per group appearing in any contrast
> (`union(Numerator, Denominator)`), holding the percent coefficient of
> variation (`100 * sd/mean`, sample SD) across that group's biological
> replicates, computed on a **linear** scale. Single-sample groups, zero means,
> and (proteomics) groups with <2 observed values yield `NA`.
>
> **Per-omics linear scale:** RNA-seq uses linear CPM (`compute_cpm`);
> proteomics uses `2^expr_filt` (observed/unimputed values only); metabolomics
> uses the post-normalization linear matrix reconstructed as
> `2^expr_work − pseudocount`.
>
> ⚠️ **Metabolomics comparability caveat:** the reconstruction is exact for
> tss/pqn/eigenms normalization and approximate for median/eigenms_forced. If
> feature scaling is enabled (`preprocessing.scaling != "none"`), the
> back-transform is invalid and CV falls back to the **pre-normalization**
> matrix (`expr_filt`); in that case metabolomics CV folds in technical
> variation (drift, total-intensity) that normalization removes and is **not**
> directly comparable to the RNA-seq/proteomics CV columns.
>
> *As implemented* (`build_group_cv_metabolomics()`): with scaling off, the
> `2^x − pseudocount` inverse is used only when the workspace is known to be
> log2 — `chosen_norm` of `tss`, `pqn`, `eigenms`, `eigenms_forced` or
> `bio_factor`, or `median`/`none` with `transform: "log2"`. Otherwise the CV
> columns are skipped with a warning. The exact/approximate statement above is
> from the design and was not re-checked.

### 5b. Column specification

| Group | Column(s) | Type | Source | Required |
|-------|-----------|------|--------|----------|
| 1. ID | `feature_id` | character | `summary_df$feature_id` | Yes |
| 1b. Original ID | `original_id` | character | `row_data`, placed directly after `feature_id` (*added after the design*) | No (only when `row_data` has it) |
| 2. Annotations | Configurable from `pre$row_data` (e.g., compound name, HMDB ID, m/z, RT, molecular formula) | character | `row_data` columns via `annot_cols` parameter | No (gracefully empty if no row_data) |
| 3. Expression | One column per sample | numeric | `pre$expr_work` (normalized log2 matrix) | Yes |
| 4a. FC | `linearFC.<cn>` | numeric (signed) | Computed from logFC in summary_df | Yes, per contrast |
| 4b. P-value | `pvalue.<cn>` | numeric [0,1] | summary_df | Yes, per contrast |
| 4c. Adj p-value | `padj.<cn>` | numeric [0,1] | summary_df | Yes, per contrast |
| 4d. Direction | `upDown.<cn>` | character: `"up"` / `"down"` / `""` | Computed from pass + FC sign | Yes, per contrast |
| 4e. Manual cutoffs | `manual_cutoffs.<cn>` | NA (Excel formula placeholder) | Static | Yes, per contrast |
| 5. Aggregate | `pass_any_contrast` | integer: `1` / `NA` (**deviation:** `1` / `0` when DE was computed rather than pre-computed; see §0) | summary_df | Yes |

### 5c. Proteomics comparison (side-by-side)

| Column group | Proteomics `Final_results.tsv` | Metabolomics `Final_results_metabolomics.tsv` | Structural match? |
|-------------|-------------------------------|----------------------------------------------|-------------------|
| ID | `FeatureID` | `feature_id` | Yes — configurable via `feature_id_col` |
| Annotations | `Protein.Names`, `Genes`, `First.Protein.Description` | Compound name, HMDB, m/z, RT (from `row_data`) | Yes — different columns but same mechanism (`annot_cols` parameter) |
| Expression | Sample columns from `pre$expr_filt` | Sample columns from `pre$expr_work` | Yes — both numeric matrices. Proteomics uses pre-imputation (may have NAs); metabolomics uses normalized, which can also hold NAs (missing values are filtered, not imputed). |
| FC | `linearFC.imputs.<cn>` | `linearFC.<cn>` | Yes — same semantics, different infix. Both signed linear FC. |
| P-value | `pvalue.imputs.<cn>` | `pvalue.<cn>` | Yes |
| Adj p-value | `padj.imputs.<cn>` | `padj.<cn>` | Yes |
| Direction | `upDown.imputs.<cn>` | `upDown.<cn>` | Yes |
| Manual | `manual_cutoffs.<cn>` | `manual_cutoffs.<cn>` | Exact match |
| Aggregate | `pass_any_contrast` | `pass_any_contrast` | Exact match |

### 5d. Builder function signature

*The design sketch below is kept for reference. The implemented function in
`R/domain/metabolomics/05_outputs_legacy.R` also takes `cv_contrasts_df` and
`config` (for the CV columns), defaults `row_data` to `pre$row_data`, moves
`original_id` directly after `feature_id`, and passes `cv_cols` to
`build_final_results_generic()`.*

```r
build_final_results_metabolomics <- function(
    pre,
    summary_df,
    contrast_labels,
    row_data      = NULL,
    feature_id_col = "feature_id"
) {
    # Wrap contrast_labels into contrasts_df for generic API
    contrasts_df <- data.frame(
        Contrast_name = contrast_labels,
        stringsAsFactors = FALSE
    )

    # Build annotation column mapping from row_data
    annot_cols <- NULL
    if (!is.null(row_data)) {
        available <- setdiff(colnames(row_data), feature_id_col)
        if (length(available) > 0) {
            annot_cols <- setNames(available, available)
        }
    }

    build_final_results_generic(
        summary_df     = summary_df,
        expr_df        = pre$expr_work,
        contrasts_df   = contrasts_df,
        feature_id_col = feature_id_col,
        annot_cols     = annot_cols,
        row_data       = row_data,
        fc_is_signed   = TRUE,
        mode           = "metabolomics"
    )
}
```

This reuses `build_final_results_generic()` (in `core/05_export_excel.R`) which
already handles:
- Column ordering (ID → annotations → expression → DE stats → aggregate)
- `upDown` computation from signed FC + pass flag
- `manual_cutoffs` placeholder
- `pass_any_contrast` propagation

The only mode-dependent behavior is `get_contrast_cols(cn, mode = "metabolomics")`
which returns column names without the `.imputs.` infix.

---

## 6. Shiny Impact

*Pipeline side implemented: `build_de_summary_counts_metabolomics()` reads
`pass.<cn>` and `linearFC.<cn>`, and `SHINY_PAYLOAD_HANDOFF.md` documents those
names. Whether the app-side audit described below took place is not recorded
in this repository. `config_padj_cutoff` does not exist; the payload key is
`padj_cutoff` (see §0). The table below is the design-time plan.*

The Shiny payload builder (`build_shiny_payload_metabolomics()` in `07_shiny_export.R`)
passes `de_res$summary_df` directly as `payload$de_stats`. After the column rename:

| Shiny consumer | Current pattern | New pattern | Impact |
|----------------|----------------|-------------|--------|
| `build_de_summary_counts_metabolomics()` | `grep("^pass_", ...)` | `grep("^pass\\.", ...)` | Must update regex |
| Same function | `paste0("logFC_", cn)` | `paste0("linearFC.", cn)` | Must update FC lookup |
| `payload$config_padj_cutoff` | Key name only | Key name only | **No change to key name** — but owner must verify the Shiny app doesn't hardcode `pass_<cn>` or `logFC_<cn>` patterns |

**Design-time note (was BLOCKED):** Owner will audit the Shiny app for hardcoded references to
`config_padj_cutoff`, `pass_<cn>`, `logFC_<cn>` patterns before implementation proceeds.

---

## 7. Risks and mitigations

| Risk | Severity | Mitigation |
|------|----------|------------|
| `extract_contrast_table()` callers expect `logFC` in output | Medium | Back-compute: `sign(linearFC) * log2(abs(linearFC))`. Keep return column named `logFC` for plot compatibility. |
| `linearFC` values of exactly 0 | Low | `2^0 = 1` → linearFC = 1 (no change), not 0. This is correct semantics. |
| Precomputed DE tables still use old format | Medium | `load_precomputed_metabolomics_de()` gets the same column rename. Old cached `.tsv` files will be re-read and re-mapped. |
| `pass.` prefix collides with `pass_any_contrast` in regex | Low | Use anchored regex: `"^pass\\.[^_]"` or exclude `pass_any_contrast` explicitly (as `build_de_summary_counts_metabolomics` already does). |
| Log2FC precision loss in round-trip | None (design); **does not hold as implemented** | Per-contrast `de_tables` (narrow format) still store raw `logFC`. Only `summary_df` stores `linearFC`. But `extract_contrast_table()` back-computes `logFC` from the rounded `linearFC` (`signif(..., 3)`) for plots and the report, so those consumers do see the round-trip. See §0. |
