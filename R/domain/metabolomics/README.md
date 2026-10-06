# Metabolomics domain

What the metabolomics mode does and where its code lives. For installing, picking
a config and running `{targets}`, see [onboarding](../../../docs/onboarding.md); for
the layer architecture, see [PROJECT_STRUCTURE](../../../PROJECT_STRUCTURE.md).

## What the mode does

It takes a feature × sample intensity table and its sample metadata. It
filters features on missingness, normalizes and transforms the matrix,
optionally corrects signal drift, and then runs differential abundance per
contrast, feature selection, pathway enrichment, QC, clustering and reports.
It also writes final-results tables (TSV + Excel) and a payload for the Shiny
app.

**There is no pipeline-level imputation.** Missing values are filtered, not
filled. A few steps that cannot take `NA` fill gaps locally for their own use
only: row medians for the PCA display, column medians inside feature
selection, and a display-only fill in the Shiny payload's `expr_norm`. DE runs
on the matrix with its missing values.

Inputs are set by `modes.metabolomics.input.format`: `cd_raw` (default),
`processed_wide`, `long` or `multi_level`. Alternatively, DE can be loaded from
pre-computed per-contrast tables (`files.de_table`) instead of being computed.

## Preprocessing flow

```
metab_inputs ─► met_raw ─► met_filtered ─► met_log ──────────► met_norm_median
                (linear)    (missingness    (log transform)
                             filter)    └─► met_norm_tss, met_norm_pqn,
                                            met_norm_eigenms, met_norm_eigenms_forced,
                                            met_norm_bio_factor            (from the linear matrix)
                                                      │
                          preprocessing.chosen_norm ──┴─► met_corrected ─► metab_pre
                                                          (+ optional LOESS drift correction)
```

`met_missingness_stats` records missingness on the filtered matrix, and
`met_norm_comparison` writes a side-by-side of the normalizations.
`metab_pre` is the adapter that every downstream step reads; its working
matrix is `expr_work`.

## Normalization: `preprocessing` and `chosen_norm`

`modes.metabolomics.preprocessing` is the **only** place normalization is set.
There is no separate `normalization:` block.

| Key | Meaning |
|---|---|
| `chosen_norm` | `none`, `tss`, `median`, `pqn`, `eigenms`, `eigenms_forced`, `bio_factor`, or `null` |
| `transform`, `pseudocount` | `transform` sets the log transform of `met_log`, so it affects only `median` and `none`. `tss`, `pqn`, `eigenms`, `eigenms_forced` and `bio_factor` always use `log2(x + pseudocount)` |
| `scaling` | Optional feature scaling |
| `drift_correction` | Optional LOESS signal-drift correction (QC samples + injection order) |

-   **`chosen_norm: null` is QC-review mode.** The pipeline builds every normalization,
    a QC suite per normalization (`met_*_qc`), `met_qc_comparison` and a QC summary
    report, and stops there: no DE and no downstream analysis. Use it to choose a
    method, then set `chosen_norm` and rerun.
-   **`chosen_norm: none`** skips sample normalization, for tables normalized upstream.
    `transform` and `scaling` still apply, so set `transform: "none"` as well if the
    table is also already log-scaled.

## Main targets after `metab_pre`

| Target | Built by | Domain code |
|---|---|---|
| `metab_de_res` | `mod_metabolomics_de()` | `03_differential.R`: `run_metabolomics_de()` (`de.method`: `limma`, `t_test` (Welch), `t_test_equal`, `wilcoxon`), `load_precomputed_metabolomics_de()`, `build_de_summary()` |
| `metab_feature_sel_res` | `mod_metabolomics_feature_selection()` | `04_feature_selection.R` (random forest, PLS-DA) |
| `metab_enrichment_res` | `mod_metabolomics_enrichment()` | `06_enrichment.R`: QEA, ssGSEA, ORA, GSEA |
| `metab_mummichog_pinned_files` | `mod_mummichog_pinned()` (only when mummichog is enabled) | `06c`–`06e_mummichog_*.R` |
| `metab_final_results` | `write_metabolomics_final_results()` | `05_outputs_legacy.R`: `build_final_results_metabolomics()`, `build_group_cv_metabolomics()` |
| `metab_shiny_payload` | `save_shiny_payload_metabolomics()` → `build_shiny_payload_metabolomics()` | `07_shiny_export.R` |
| `metab_report` | `mod_metabolomics_report()` | `R/pipeline/metabolomics/templates/report_metabolomics.Rmd` |

Other domain files: `08_metabolite_network.R`, `10_drift_correction.R`,
`10_pipeline_summary.R`, `11_powerpoint.R`. The preprocessing modules
(`mod_met_*`) are in `R/modules/metabolomics/00_mod_preprocessing.R`.

User-facing text says "features", "metabolites" and "differential abundance".
Internal names keep `de` (`de:` config key, `metab_de_res`, `de_stats`); renaming
them is a separate, deliberate task (see
[DE_to_DA_rename_map.md](../../../docs/DE_to_DA_rename_map.md)).

## Related documents

-   [CONTRACT_metabolomics_de_and_final_results.md](../../../docs/CONTRACT_metabolomics_de_and_final_results.md):
    DE column names and the final-results schema, with the known deviations listed in its §0.
-   [mummichog.md](../../../docs/mummichog.md): the optional mummichog pathway analysis.
-   [SHINY_PAYLOAD_HANDOFF.md](../../../docs/SHINY_PAYLOAD_HANDOFF.md): what the Shiny
    payload contains, including the metabolomics-specific keys.
