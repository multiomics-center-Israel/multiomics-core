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
`processed_wide` or `multi_level`. Alternatively, DE can be loaded from
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
`metab_pre` is the adapter that the downstream steps read, directly or through `metab_de_res`; its working
matrix is `expr_work`, which carries the optional feature `scaling`. Differential
abundance and the QEA/ssGSEA enrichment test `expr_unscaled` instead: the same
normalized, transformed matrix without feature scaling, chosen by
`resolve_metab_de_matrix()`. Changing `scaling` therefore changes the
multivariate views (PCA, sample distance and correlation heatmaps, Shiny
`expr_norm`) but not the DE statistics.
`expr_log` is the transform alone, with no sample normalization.

## Normalization: `preprocessing` and `chosen_norm`

`modes.metabolomics.preprocessing` is the **only** place normalization is set.
There is no separate `normalization:` block.

| Key | Meaning |
|---|---|
| `chosen_norm` | `none`, `tss`, `median`, `pqn`, `eigenms`, `eigenms_forced`, `bio_factor`, or `null` |
| `transform`, `pseudocount` | `transform` sets the log transform of `met_log`, so it affects only `median` and `none`. `tss`, `pqn`, `eigenms`, `eigenms_forced` and `bio_factor` always use `log2(x + pseudocount)` |
| `pqn_reference` | PQN reference: `all_samples` (default) or `pooled_qc` (the samples flagged in `drift_correction.qc_flag_col`; at least 2, otherwise the run stops). Recorded in `info$normalization` and `normalization_applied.tsv`. Until this was made explicit, PQN used the flagged QC samples automatically whenever at least 2 were present, so such projects can give different results when rerun; set `pooled_qc` to keep that behaviour |
| `scaling` | Optional feature scaling |
| `drift_correction` | Optional LOESS signal-drift correction (QC samples + injection order). Opt-in: runs only with `enabled: true` |

-   **`chosen_norm: null` is QC-review mode.** The pipeline builds the normalizations,
    a QC suite for the log matrix and for `tss`, `median`, `pqn`, `eigenms` and
    `eigenms_forced` (`met_*_qc`), `met_qc_comparison` and a QC summary report. It
    stops there: no DE and no downstream analysis. **`bio_factor` is not part of the
    review-mode QC or the comparison**, even when it is configured. Use review mode
    to choose a method, then set `chosen_norm` and rerun.
-   **Review mode supplies no layer to multi-omics integration.** Metabolomics takes
    part in integration only with a concrete `chosen_norm`.
-   **`chosen_norm: none`** skips sample normalization, for tables normalized upstream.
    `transform` and `scaling` still apply, so set `transform: "none"` as well if the
    table is also already log-scaled.
-   **LOESS drift correction is opt-in.** It runs only when
    `drift_correction.enabled` is explicitly `true` and the metadata has the
    injection-order and QC columns with at least 4 QC samples. A missing or empty
    `drift_correction` block means off. Until this change a missing `enabled` meant
    on, so **rerunning a project that never set it can give different results**
    when its metadata met those conditions. To keep the old behaviour, set
    `enabled: true` explicitly. **LOESS cannot be combined with feature scaling:**
    with `scaling` other than `"none"` the run stops, because the correction is
    not valid on centered or scaled features (see #362).

## Main targets after `metab_pre`

| Target | Built by | Domain code |
|---|---|---|
| `metab_de_res` | `mod_metabolomics_de()` | `03_differential.R`: `run_metabolomics_de()` (`de.method`: `limma`, `t_test` (Welch), `t_test_equal`, `wilcoxon`), `load_precomputed_metabolomics_de()`, `build_de_summary()` |
| `metab_feature_sel_res` | `mod_metabolomics_feature_selection()` | `04_feature_selection.R` (random forest, PLS-DA) |
| `metab_enrichment_res` | `mod_metabolomics_enrichment()` | `06_enrichment.R`: QEA, ssGSEA, ORA, GSEA |
| `metab_mummichog_pinned_files` | `mod_mummichog_pinned()` (only when mummichog is enabled), or `mod_mummichog_metaboanalyst()` with `mummichog.engine: metaboanalyst` | `06c`–`06e_mummichog_*.R`; MetaboAnalyst-compatible engine `06j`–`06n_mummichog_ma_*.R` |
| `metab_mummichog_sections` | `mod_mummichog_sections()` | `06e_mummichog_plots.R`: `build_mummichog_report_sections()` (ORA plot + table, GSEA, evidence, Summary) |
| `metab_mummichog_exports` | `mod_mummichog_exports()` | `06e_mummichog_plots.R`: `save_mummichog_exports()`, `clear_orphan_mummichog_gsea_exports()` |
| `metab_final_results` | `write_metabolomics_final_results()` | `05_outputs_legacy.R`: `build_final_results_metabolomics()`, `build_group_cv_metabolomics()` |
| `metab_shiny_payload` | `save_shiny_payload_metabolomics()` → `build_shiny_payload_metabolomics()` | `07_shiny_export.R` |
| `metab_report` | `mod_metabolomics_report()` | `R/pipeline/metabolomics/templates/report_metabolomics.Rmd` |

## Enrichment reference files

Two files outside the data folder decide enrichment results, and both are
tracked by `{targets}`:

| Config key | Read by | Target that carries its fingerprint |
|---|---|---|
| `enrichment.gmt_file` (one or a list) | `metab_enrichment_res`, `metab_mummichog_sections` | `metab_enrichment_gmt_ref` |
| `enrichment.mapping_file` (HMDB → KEGG) | `met_raw`, `metab_enrichment_res`, `metab_mummichog_sections` | `metab_kegg_mapping_ref` |

The fingerprint is the md5 of each existing file, taken each time the plan is
defined and written into the target's command
(`06q_enrichment_reference_files.R`). So a library edited in place, or a
configured library that was missing and now exists, reruns enrichment; a changed
mapping table reruns `met_raw` and enrichment; and when nothing changed, nothing
is outdated. A changed GMT does not rerun `met_raw`.

When mummichog is enabled, its report sections (`metab_mummichog_sections`)
read the same two files themselves: the GMT libraries to decide which mummichog
pathways can be compared with the identified-metabolite enrichment, and the
mapping table for the supporting evidence. The sections therefore depend on
both fingerprints directly, and rerun when either file changes even if the
enrichment result and the preprocessed data come back unchanged (for example
with `run_enrichment: false`). The mummichog run itself reads neither file, so
a changed GMT does not rerun it.

This is separate from **mummichog's own reference files** — the pathway model
it runs on (the built-in one, or `model_ref` / `model_json`). Those are not
tracked by these two targets and are outside this mechanism.

A configured file that does not exist is not tracked and is not an error, as
before: the enrichment code skips it. For a relative `mapping_file`, `met_raw`
reads the path as written and the enrichment stage resolves it under the raw
data folder; both locations are tracked.

## ORA and GSEA provenance

ORA and GSEA test **compounds**, and several measured features can map to one
compound. Each result therefore carries a `provenance` element beside its
unchanged `table` (and `ranks`), built in `06r_enrichment_provenance.R` from
the objects the test computed. It never changes a statistic, and no column is
added to any table or `results.tsv`.

| Element | What it records |
|---|---|
| `status`, `diagnostics`, `error` | `"ok"`, `"inconsistent"` (a check against the result failed) or `"failed"` (could not be built) |
| `reference_files` | the GMT and mapping files the test was given, whether each was read, and its md5 taken before the read; a library is hashed again after its read, and a file that changed while the test ran makes the provenance inconsistent |
| `feature_compounds` | per measured feature: the compound id used, which step chose it (`kegg_column` > `hmdb_mapping` > `name_column`), every KEGG id in the cell, and whether the cell was ambiguous |
| `pathway_sources` | every GMT line and what became of it: tested, replaced by a later line with the same id, replaced by an equal label earlier in the same file or in a later library, too few matches, or not returned by fgsea |
| `pathway_members` | every listed member of every set given to the test, duplicates kept |
| ORA: `thresholds`, `compounds`, `features` | universe and foreground per compound; significance per feature, separately |
| GSEA: `ranked`, `features` | the final ranking with the representative feature of each compound; every feature's role |

Three distinctions the tables cannot make:

-   **Feature roles in ORA.** A feature can be non-significant and still belong
    to a foreground compound, because another feature of that compound is
    significant. Counts are per compound; `features` keeps significance,
    universe and foreground as separate flags.
-   **The GSEA representative.** When several features share a compound, the
    ranking keeps the one with the largest |logFC|, and on a tie the one that
    comes first in the DE table. `ranked$representative_feature` is the
    feature that was kept, recorded during the deduplication, not worked out
    again from the logFC. Ties are flagged, including ties of opposite sign,
    where the row order decides the sign the compound is ranked with.
-   **Given to the test versus used by it.** A set can list one compound
    twice, typically after HMDB ids are translated to KEGG. `in_input_set`
    marks every occurrence in the vector the test received;
    `in_effective_statistic` marks what counted. ORA counts every occurrence in
    the pathway size, so the two are equal. fgsea keeps each member once, and
    drops a set that is then too small. Both behaviours are as they were.

`metab_pathway_feature_membership(provenance)` lists every measured feature of
every pathway compound for display, with those roles per row. It shows more
than the test used: in GSEA a feature that did not represent its compound, or
a compound that was not ranked, is listed and never marked as contributing to
the score or the leading edge.

Code that needs provenance calls `require_metab_enrichment_provenance(res)`,
which stops unless the status is `"ok"`. A failed or inconsistent provenance
leaves the statistical result untouched and raises a warning. The warning can
never end the analysis: where warnings are turned into errors
(`options(warn = 2)`) it is dropped, and the status stays on the record.

`upstream_origin = "not_tracked"` on an id taken from the `KEGG` column means
the enrichment stage read the column as it found it: it cannot tell whether
the value came with the input annotation or was filled earlier from a mapping
file.

## ORA and GSEA in the Shiny payload

`build_metab_enrichment_payload()` (`07b_shiny_enrichment_payload.R`) turns a
run's ORA and GSEA results into the two keys the Shiny payload carries side by
side. It is a re-encoding of results and provenance: nothing is recomputed, and
a result whose provenance is not `"ok"` stops it.

There is **no option to set**. `build_shiny_payload_metabolomics()` adds the two
keys whenever the run has enrichment results (`run_enrichment: true`):

-   No existing key changes: the legacy `enrichment_ora`, `enrichment_gsea`,
    `enrichment_qea` and `enrichment_ssgsea` tables stay, and `mummichog` is
    built on its own, with or without enrichment results.
-   A run without enrichment adds nothing.
-   A run whose ORA and GSEA produced no result gets a block with
    `available = FALSE` and every contrast listed as `no_result`.
-   Results that cannot be encoded correctly (for example results cached
    before provenance was recorded) are not added: the payload is still
    written with every other key, a warning names the reasons, and
    `enrichment_provenance` holds `available = FALSE`, `reason = "blocked"`
    and the preflight table.

-   **`enrichment`** has exactly the nine components the app reads for RNA-seq
    (`available`, `gene_index`, `config`, `manifest`, `ora`, `gsea`,
    `gsea_leading_edge`, `gsea_rankings`, `pathway_membership`), under database
    `combined` and ranking `fc`. Ids are feature ids.
-   **`enrichment_provenance`** holds what that shape cannot: every tested
    pathway with its full statistics (`all_tested`), compounds, per-feature
    roles, pathway members with multiplicity, and the state of each analysis.

The two are joined by analysis, contrast and pathway id. Four things are kept
apart throughout:

| What | Where |
|---|---|
| Significant measured features behind the ORA foreground | `geneID` and `n_features` of an ORA row; `features$significant` |
| Compounds the ORA test counted | `Count`, `GeneRatio`, `BgRatio`, `FoldEnrichment` |
| Pathways significant in a result | the tables under `enrichment$ora` / `enrichment$gsea` |
| Exploratory membership, significant or not, ranked or not | `enrichment$pathway_membership` |

-   The display tables hold the significant pathways only, at the threshold
    the enrichment stage itself reports with (`.METAB_ENRICHMENT_SIG_FDR`).
    Non-significant tested pathways are in `all_tested`.
-   GSEA rankings and leading edges are in the representative feature of each
    compound; no other feature of a compound appears in them.
-   `gene_index` starts with the feature ids of the exported matrix, in its row
    order. Every other feature the block or its provenance names follows, and
    `enrichment_provenance$feature_index$in_expr_norm` tells them apart.
-   An analysis that ran and found nothing has a manifest row with status
    `empty` and a zero-row table. One with no result has no manifest row and
    is `no_result` in `enrichment_provenance$status`.
-   A pathway that was given to fgsea and came back without a result (its
    members collapse below fgsea's size limit once duplicates are dropped) is
    listed in `enrichment_provenance$gsea[[contrast]]$no_result`, apart from
    `all_tested`, and the detail reader returns it with state `no_result`.
-   The manifest rows are the layout's own; `storage_key` here names the unit's
    place in the block, not a directory.

`preflight_metab_enrichment_payload()` reports, without stopping, what would
block the build on a given run: results without usable provenance, a feature id
containing `/`, duplicated ids in the exported matrix, or two tested pathways
sharing an id. `metab_enrichment_pathway_detail(payload, analysis, contrast,
pathway_id)` returns the result row, the listed compounds and the measured
features behind any tested pathway, from the payload alone.

## Enrichment output layout

`mod_metabolomics_enrichment()` writes each result unit to its own directory
under `<mode>/Enrichment/`. The layout follows the RNA-seq one: the analysis
context is in the path, and file names inside a unit are short and fixed. All
paths are built in `06p_enrichment_layout.R`.

```
Enrichment/
├── enrichment_manifest.csv            the units this run wrote (same columns as RNA-seq)
├── ORA/combined/contrasts/<contrast>/de/       results.tsv, dotplot.png, lollipop.png, barplot.png
├── GSEA/combined/ranking_by_fc/<contrast>/     results.tsv, nes_barplot.png, lollipop.png, barplot.png
├── QEA/combined/all_samples/                   results.tsv, barplot.png, lollipop.png, library_<library>.png
├── QEA/combined/contrasts/<contrast>/          results.tsv, lollipop.png, barplot.png
├── ssGSEA/combined/all_samples/                results.tsv, barplot.png, boxplots.png, lollipop.png
└── ssGSEA/combined/contrasts/<contrast>/       results.tsv, lollipop.png, barplot.png
```

-   **`combined`** is the database level. Each method adjusts FDR across all
    configured GMT libraries together, so there is one result family per unit;
    the library of a pathway is the `library` column of the table where the
    method records it (QEA, ORA, GSEA).
-   **`<contrast>`** is the contrast label made file-safe. Labels that collide
    after that get a numeric suffix, with a warning, and a label that is a
    Windows device name (`CON`, `NUL`, `COM1`, ...) is prefixed with
    `contrast_`.
-   **`enrichment_manifest.csv`** is removed when an enrichment run starts (the
    run stops, before writing anything, if it cannot be removed) and written
    once all results are on disk, through a temporary file that is
    renamed into place, so it is either complete or absent: a run that stopped
    part-way has none. Only units that produced a table are listed. Its columns
    are the RNA-seq ones, **except that `storage_key` means something else
    here**: the unit directory relative to `Enrichment/`, not the key of a table
    inside a results object as in RNA-seq.
-   **Files from an earlier run.** After the manifest is written, the run
    removes the files of this stage that it did not write: a contrast that is
    gone or gave no result, a plot that was not drawn this time, the flat
    `enrichment_*` files that preceded this layout. A file counts as this
    stage's only when both its path and its name are ones the stage produces:
    -   in a unit directory of the tree above: `results.tsv`, the plots listed
        for that method, and `library_*.png` in `QEA/combined/all_samples/`;
    -   directly in `Enrichment/`: `enrichment_<method>_<plot>[_<contrast>].png`
        for those same plot names, and `enrichment_qea_*.png`;
    -   directly in `Datasets/`: `enrichment_<method>[_<contrast>]_results.tsv`.

    Any other file is left alone, including another PNG inside a method folder.
    Directories are never removed. Nothing is removed by a run that fails
    (including a failure to write the manifest), a run that produces no result
    at all, or a run with enrichment disabled. A file that could not be deleted
    is named in a warning, not reported as removed.
-   **Still written for compatibility:** the first contrast's ORA and GSEA table
    in `Datasets/enrichment_{ora,gsea}_results.tsv` and its plots as
    `Enrichment/enrichment_{ora,gsea}_*.png`.
-   The mummichog stage is separate and keeps its own `mummichog_pinned/` tree.

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
