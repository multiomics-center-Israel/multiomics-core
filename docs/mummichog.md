# Mummichog pathway analysis (pinned v2, isolated venv)

> Moved here from the top-level README; the stochasticity caveat has since been corrected to match the code. Back to the [docs index](README.md).

> **Engine.** `engine: metaboanalyst` (the MetaboAnalyst-compatible engine, ORA + GSEA with joint mixed polarity; [below](#metaboanalyst-compatible-engine-joint-mixed-polarity)) is the **default**: an enabled config that names no engine runs it, and the stage says so in its log. `engine: python` selects the pinned python mummichog engine described in the next sections, and must be named to be used. The default is not a fallback: if the MetaboAnalyst-compatible engine's prerequisites are missing (a model JSON, a supported `ionization_mode`, `fitdistrplus`), the stage stops; it never switches to python. A disabled config selects no engine. The engine selected for the last enabled run is recorded in `mummichog_pinned/mummichog_engine.json`. Each stage also checks the selection itself: the python stage stops unless the config selects `python`. **Re-running an older config that has no `engine` may not reproduce it** — see [Backward compatibility](#backward-compatibility-configs-without-engine).

The metabolomics mode runs [mummichog](http://mummichog.org) for m/z-based pathway/network enrichment via a **version-pinned, isolated engine** (`R/domain/metabolomics/06c_mummichog_pinned.R`): `mummichog==2.7.0` invoked as a `{processx}` subprocess in a dedicated venv, depending only on light R packages (`readr`, `processx`, `jsonlite`) — no Bioconductor. It runs on mummichog's built-in `human_mfn` model by default.

## One-time setup

The pinned engine calls Python in a dedicated venv, kept out of git (`envs/` is `.gitignore`d). Once per machine (or checkout):

``` bash
make setup
```

That builds the venv (`envs/mummichog`) and **prints the exact `MUMMICHOG_PYTHON=<path>` line for this checkout**. Add just that line to your `.Renviron` (create the file in the project root if you don't have one — it's `.gitignore`d):

``` bash
# append the line make setup printed, e.g.:
echo 'MUMMICHOG_PYTHON=/abs/path/to/envs/mummichog/bin/python' >> .Renviron
```

> If you'd rather start from the tracked template with `cp .Renviron.example .Renviron`, also **set or remove its `MULTIOMICS_CONFIG=/path/to/your/config.yaml` placeholder** — an active dummy value there overrides the `config.yaml` default and makes `tar_make()` fail before it reaches mummichog.

After that, R reads `.Renviron` on start and `targets::tar_make()` just works — no manual `export` each session. **`.Renviron` is machine-specific and `.gitignore`d — never commit it** (that's why `make setup` prints the line for you to add rather than writing the file itself). A relative path works if you always start R from the project root, but the absolute path `make setup` prints is more robust.

<details>
<summary>Manual / advanced use</summary>

``` bash
make mummichog-venv                 # creates envs/mummichog, writes requirements-mummichog.lock
# or, to reproduce the exact committed tree:
make mummichog-lock                 # installs from requirements-mummichog.lock (USE_LOCK=1)

# instead of .Renviron, you can export the interpreter path per shell:
export MUMMICHOG_PYTHON="$(pwd)/envs/mummichog/bin/python"
```

On Windows the venv interpreter is at `envs\mummichog\Scripts\python.exe` instead (both `make setup` and the pipeline pick the right path per platform). Both `requirements-mummichog.txt` (the top-level pin) and `requirements-mummichog.lock` (the fully-resolved tree) are committed.

If the venv/interpreter is missing when the stage runs, the pipeline fails loudly and names the fix (`run make setup`) — it never silently builds a venv mid-run.
</details>

## How to run

It's wired into the metabolomics DAG (as `metab_mummichog_pinned_*` targets) and is **opt-in via config** — set `enabled: true` under `modes.metabolomics.enrichment.mummichog`. When disabled or omitted, the targets aren't added to the graph and the Python venv is never needed.

``` yaml
modes:
  metabolomics:
    enrichment:
      mummichog:
        enabled: true
        p_cutoff: 0.05
        n_permutations: 100
        tolerance_ppm: 10
        ionization_mode: pos_default   # pos_default | positive | negative | mixed (see below)
        force_primary_ion: true        # require a primary ion; false allows non-primary adducts
```

`force_primary_ion` maps to mummichog's `-z`. mummichog 2.7.0 **requires a primary ion** (`M+H[+]` for positive, `M-H[-]` for negative) to be present before accepting a metabolite prediction — this filters out noise from irrelevant adducts and is the engine's **default**. Set `force_primary_ion: false` to relax that (emits `-z False`, keeping adduct-only predictions); omit the key to keep the default. It maps to MetaboAnalyst's `force_primary_ion` option.

Then run as usual:

``` r
library(targets)
tar_make(names = tidyselect::starts_with("met"))
```

> **Per-contrast:** mummichog runs **independently for each differential-abundance contrast** (each contrast's own p-values define its significant set against all features, sharing one model + params). Every contrast with a result renders as its own tab in the HTML report's mummichog section and gets its own files on disk (see below).
>
> **Organism:** the built-in model is **human only**. A non-human `modes.metabolomics.organism` with no custom model is rejected with a clear error rather than silently run against the human network — supply an organism-specific model (see below).

## Mixed-polarity data

mummichog runs **one ionization mode per run** and matches every submitted m/z against that mode's ion forms; it never reads a feature's own ion label. A single mode over a table that holds both positive- and negative-ion features therefore matches the other polarity's features against the wrong ion forms. For such data, set:

``` yaml
        ionization_mode: mixed
        polarity_column: "Reference Ion"   # required
        positive_mode: positive            # required: pos_default | positive
```

Each contrast is then run **once per polarity**: the positive-ion features in `positive_mode`, the negative-ion features in `negative`, each with **only its own polarity's features** as background and significant list. Polarity is read per feature from `polarity_column` — an ion label ending in its charge sign (`[M+H]+1`, `[M-H]-1`, `[M-2H]2-`, `[M]+`) or the word `positive`/`negative`. A missing column, or any feature whose polarity cannot be read, **stops the stage**; no polarity is assumed. A polarity with no features is skipped with a message and the other still runs.

Each run is its own unit, labelled `<contrast>, positive ions` / `<contrast>, negative ions`: its own directory (`mummichog_pinned/<contrast>_positive_ions/`), its own report tab and its own export suffix. The two runs' p-values are not combined. The identified-metabolite support shown in each tab is contrast-level (full dataset), not polarity-specific. Single-polarity modes (`pos_default`, `positive`, `negative`) behave exactly as before.

This per-polarity behaviour is the **python engine's**. For one joint analysis of both polarities, as MetaboAnalyst runs it, use the MetaboAnalyst-compatible engine below.

## MetaboAnalyst-compatible engine (joint mixed polarity)

`engine: metaboanalyst` replaces the python subprocess with a native R reimplementation of the mummichog algorithm in MetaboAnalyst's MS Peaks-to-Pathways module. The reference is MetaboAnalystR at commit `398476ae2a0c996e390925fa62adc46b5ecde334`. Each step was ported from that commit's source and checked against that commit's own output on synthetic fixtures (`tests/testthat/fixtures/mummichog_ma_parity/`):

- the compound library built from our model JSON;
- peak preparation and per-polarity m/z matching;
- empirical compounds;
- permutations and pathway statistics.

This engine is the default: it runs when `engine` is `metaboanalyst` or absent. The python engine is untouched and selectable with `engine: python`.

What only the python engine provides (lost if it is ever removed): mummichog's modular/network analysis outputs and `result.html`, the `pos_default` ion list, the built-in `human_mfn` model, and its own GSEA variant on python mummichog's m/z-merged ECs. Both engines produce the Pathway Evidence section; see [Pathway Evidence with this engine](#pathway-evidence-with-this-engine).

``` yaml
        engine: metaboanalyst
        ionization_mode: mixed             # or positive / negative
        polarity_column: "Reference Ion"   # required for mixed
        model_ref: { url: ..., sha256: ... }   # or model_json; required
        # optional (MetaboAnalyst defaults): rt_frac: 0.02, min_pathway_size: 3
```

What differs from the python engine:

| | python engine | `engine: metaboanalyst` |
|---|---|---|
| `mixed` | one run per polarity, each with its own background | **one joint run**: one peak list, one EC grouping, one background, one permutation null. Each feature is still matched only against its own polarity's ion forms, and a compound's positive- and negative-ion matches at similar retention time form one EC |
| ion lists | `pos_default` (8) or `positive` (25); `negative` | MetaboAnalyst's: 25 positive, 16 negative (`pos_default` refused) |
| significant features | `p < p_cutoff` | `p <= p_cutoff` (MetaboAnalyst's rule) |
| headline pathway p | mummichog's empirical p-value | **FET** (Fisher's exact test), MetaboAnalyst's plotted p-value |
| other pathway statistics | — | EASE, Gamma-adjusted, empirical (from permutations) |
| adjusted p-values | none | **Holm** (`adj.p Holm *` columns), MetaboAnalyst's `p.adjust()` default — **not BH/FDR** |
| permutation seeds | `params.seed` | `123 + i` for permutation `i` (MetaboAnalyst's own) |
| model | built-in `human_mfn` or a custom JSON | a custom JSON (`model_ref` / `model_json`) is required |
| GSEA | standard fgsea on python mummichog's ECs (m/z-merged scores) | standard fgsea on the **same joint-polarity ECs and pathway sets as the ORA**; EC score = max of its features' signed DE statistic; `padj` is **BH**. Not MetaboAnalyst's GSEA unless `gsea_method: metaboanalyst_legacy` (see [GSEA method](#gsea-method)) |
| Pathway Evidence | from python mummichog's EC tables; each feature's annotation compared with all of its EC's pathway candidates | from the engine's matched-feature table; each feature's annotation compared only with the pathway candidates it matched itself ([below](#pathway-evidence-with-this-engine)) |

Reproduced on purpose, because MetaboAnalyst does it:

- **Duplicate m/z across polarities:** a positive- and a negative-ion feature with the *same* m/z keep only the one with the smaller p-value.
- **Gamma parameterization:** the Gamma p-value is `1 - pgamma(1 - FET, shape, rate = <scale estimate>)`, so the scale estimate is passed as the rate.
- **Row order:** rows are ordered by Gamma as text.
- **Rounding:** values are rounded to 5 significant digits before the Holm adjustment.
- **fitdistrplus warnings:** these may concern the standard errors of the Gamma fit; the engine uses only the estimates.

Outputs, per contrast, in `mummichog_pinned/<contrast>/metaboanalyst/`:

- `ma_pathway_results.tsv`: the pathway table.
- `mummichog_pathway_enrichment_mummichog.csv`: MetaboAnalyst's own result table, verbatim.
- `ma_matched_features.tsv`: every match with its EC, feature id and polarity.
- `ma_input.tsv`: the peak list the analysis was given.
- `ma_gsea_results.tsv`: the GSEA table ordered by p-value (raw p, BH `padj`, NES, leading edge, each pathway's ECs). The default method writes NES, p and `padj` unrounded; the legacy method writes them at MetaboAnalyst's 4 significant digits.
- `ma_run_info.json`: parameters, model, Gamma fit, GSEA settings and provenance, and counts. Its `gsea` block names the `method` that produced the table, how its p-values were estimated (`reference.method`, and `params.nperm_simple` or `params.n_perm`) and, for the default method, the tie summary (`ties`), the number of tested pathways without a p-value and a summary of fgsea's error estimate (`log2err`). A record with no `method` was written before the method was recorded and is read as `metaboanalyst_legacy`. A record is always described by what it stores, never by the current default (see [Earlier results](#earlier-results)).

GSEA is optional: if it fails, the stage warns, the contrast's ORA is still written, and `ma_run_info.json` records the GSEA error. If a whole contrast fails, the stage warns and its earlier result files are replaced by `ma_failed.json`, which holds the error. That contrast gets no report section, and all of its earlier exports are cleared. If every contrast fails, the stage first clears every contrast's earlier exports and then stops.

The report, its exports and the Shiny payload read these through the same targets as before. The report's ORA table adds the EASE, Gamma, empirical and Holm columns; its GSEA tab, scatter and GSEA table export show this engine's GSEA. Python GSEA never runs with this engine, so the two engines' results are never mixed in one report.

By default the GSEA is standard fgsea and is **not** MetaboAnalyst's (see [GSEA method](#gsea-method)); the ORA above remains MetaboAnalyst-compatible. With `gsea_method: metaboanalyst_legacy` the GSEA reproduces `.compute.mummichog.RT.fgsea()` at the pinned commit, with two documented exceptions that change no statistic: MetaboAnalyst's result CSV attaches EC hits in library order to rows sorted by p-value (here each row carries its own pathway's ECs), and a single tested pathway stays a row. That legacy path is validated end to end against that commit on the synthetic fixtures (`generate_gsea_reference.R`).

Limits:

- Retention time is required. MetaboAnalyst's no-RT, compound-level path is not reproduced.
- In mixed mode, both polarities' retention times must be on one scale.
- Pathway names must be unique in the model, because MetaboAnalyst identifies pathways by name. A model that repeats a name stops the stage.
- MetaboAnalyst's custom-adduct and currency-customisation options are not offered. The default currency list is applied.

### Pathway Evidence with this engine

The report's Pathway Evidence tab, the Summary's *Supporting evidence* and *Annotation agreement* columns, and the `mummichog_evidence_{pathways,empirical_compounds,features}_<contrast>` exports are built from this engine's own results (`build_mummichog_ma_pathway_evidence()`, `R/domain/metabolomics/06s_mummichog_ma_evidence.R`). Nothing is recomputed:

- a pathway's supporting ECs are its significant ECs from `ma_pathway_results.tsv` (`ec_hits`), so `Supporting ECs` equals `Sig.EC`; the evidence is refused if they ever disagree;
- each EC's compounds, features, ion forms and polarities come from `ma_matched_features.tsv`; a feature's p-value and input statistic from `ma_input.tsv`; its `Feature log2FC (DE)` from the DE results; candidate names and pathway membership from the metabolic model;
- a **pathway-matching candidate** is a compound matched in the EC that belongs to the pathway in the model.

What differs from the python evidence:

- **Agreement** compares a feature's annotation only with the pathway candidates **that feature itself matched**. MetaboAnalyst merges the ECs of features that match several compounds, so an EC can hold a feature matched only to compounds outside the pathway; that feature reads `Not assessed`, and `Matched via pathway candidate` is `FALSE`. EC membership and the support counts are unchanged.
- The feature table adds `Polarity`, `Matched compound(s)`, `Pathway candidate(s) matched by feature`, `Matched via pathway candidate`, `Feature statistic (MetaboAnalyst input)` (the signed statistic the engine was given: the t statistic when the DE method reports one, else the log2FC) and `Feature log2FC (DE)`. It has no `Feature log2FC (mummichog input)` column. `Feature log2FC (DE)` is joined from the DE results by feature id; when no DE table can be matched to the contrast, the table lacks `feature_id` / `logFC`, or some evidence features have no logFC, its cells are empty, the run warns, and the Pathway Evidence tab (and the feature-level table under Advanced) says which of these happened. The EC table adds `Polarities`.

Feature ids come from the engine's own m/z-key join, written into `ma_matched_features.tsv`. The evidence refuses a contrast whose matches cannot be traced to exactly one input feature (a missing id, an id absent from `ma_input.tsv`, a repeated input id, an m/z mapped to several ids or the reverse, or a polarity that differs between the two files).

Limits: a feature dropped as a duplicate m/z, one outside the m/z range, one matched to no library compound and one whose only matches were removed (currency compounds, the primary-ion rule) were never part of an EC, so they appear in no evidence row; for the first two the engine records how many, not which. Evidence re-reads the model, so it is built only when the configured model's sha256 equals the one in the contrast's `ma_run_info.json`. Otherwise, or when anything above is refused, the stage's ORA and GSEA are untouched, the run warns, and the contrast's Pathway Evidence tab shows the reason instead of disappearing.

## Backward compatibility: configs without `engine`

What an enabled config with no `engine` key does has changed more than once, so **the absence of the key does not tell you which engine an older analysis used**:

| Code | Enabled config with no `engine` |
|---|---|
| Before the `engine` key existed | ran the python engine (the only one) |
| Development commits of the MetaboAnalyst-compatible engine, before it was merged to `main` | first ran **python** (the key was opt-in), then ran **metaboanalyst** (stopping only if the config had `ionization_mode: pos_default` or `positive_mode`), then stopped |
| `main`, from the merge of the MetaboAnalyst-compatible engine until this change | stopped and asked for the key; nothing ran |
| Now | runs **metaboanalyst**, and records `engine_source: "default"` |

In addition, the real-data validation scripts in `tests/testthat/fixtures/mummichog_ma_parity/` (`validate_real_run.R`, `build_ma_report_run.R`) set `engine: metaboanalyst` in memory on top of a stored run, so their outputs are MetaboAnalyst-compatible while the run's config file names no engine.

So an older analysis whose config has no `engine` may have used either engine. **To reproduce an older analysis, set `engine` explicitly to the engine that analysis actually used.** Leaving the key out now always means `metaboanalyst`: an analysis that ran python will give different results (different statistics, headline p-value, ion lists and mixed-polarity handling — see the table above), or will stop on a MetaboAnalyst prerequisite such as the missing model JSON.

To find out which engine an older run used, read the run, not the config:

- **Its report.** The mummichog methodology text names the MetaboAnalyst-compatible engine (joint-polarity ORA and GSEA, FET and Holm columns) when that engine produced the results; a python run, for mixed data, has one tab per polarity. A Pathway Evidence section does not tell them apart: python reports can have one, MetaboAnalyst-compatible reports built before this change have none, and later ones can have one too.
- **Its output tree.** `mummichog_pinned/<contrast>/metaboanalyst/ma_run_info.json` exists only for the MetaboAnalyst-compatible engine; `mummichog_pinned/<contrast>/v2/` (with `mummichog_manifest.tsv`) only for python.
- **Its code version.** `execution_info/git_commit.txt` gives the commit the run was made from.

From this change on, `mummichog_pinned/mummichog_engine.json` answers the question directly:

``` json
{
  "engine": "metaboanalyst",
  "engine_source": "default",
  "configured_engine": null
}
```

`engine` is the engine selected for the run; `engine_source` is `"config"` when the config named it and `"default"` when it did not; `configured_engine` is the config value as written. The file records the engine selected for the **last run in which mummichog was enabled**. It does not say that the analysis completed — the engine's own result files show that — and it does not say whether mummichog is enabled in the current config: if mummichog is later disabled and the pipeline is rerun into the same output directory, the record stays where it is, like the earlier result trees and exports beside it, and still names the engine that produced them. Nothing under `mummichog_pinned/` is removed when mummichog is disabled. It is written by its own target (`metab_mummichog_engine_info`), only when mummichog is enabled, and before the engine stage runs, so a stage that stops still leaves the current selection on disk.

## Choosing a metabolic model

The `-n` model is selected from the `mummichog` config block with this precedence:

1.  **`model_ref`** — a published model fetched by URL and verified against its `sha256`, then cached under `envs/mummichog-models/<sha256>.json` (a gitignored dir). This is the preferred way to run organism-specific models without committing large JSON into the repo: the file is downloaded once, checked, and reused on later runs as long as its content still matches the digest. A sha256 mismatch is a hard error — an unverified model is never used.
2.  **`model_json`** — a path to a local model JSON on the machine running the pipeline.
3.  **built-in `human_mfn`** — mummichog's bundled human model (the default).

``` yaml
modes:
  metabolomics:
    organism: "Caenorhabditis elegans"     # non-human -> a custom model is required
    enrichment:
      mummichog:
        enabled: true
        model_ref:
          url: https://github.com/multiomics-center-Israel/multiomics-annotation-prep/releases/download/cre_kegg_20260711/cre_kegg_20260711.json
          sha256: c403c96fbec8df9ae34b828fec01270c8ea3940acc36e4e5ff770868dc8b912b
```

Supplying any custom model (`model_ref` or `model_json`) also satisfies the human-only guard, so a non-human organism runs against its own network.

## Where outputs land

Under `<metab_out_dir>/mummichog_pinned/`, one subdirectory per contrast (`<contrast>/`, the contrast name sanitised to `A-Za-z0-9_`):

-   `<contrast>/input.tsv` and `<contrast>/input.tsv.idmap.tsv` — the exact table sent to mummichog for that contrast (m/z, retention time, p-value, statistic, **feature\_id as the 5th column**) plus a provenance id-map.
-   `<contrast>/v2/<timestamp>.<project>/` — the mummichog result tree: `result.html`, `tables/` (`mcg_pathwayanalysis_*.tsv`/`.xlsx`, `mcg_modularanalysis_*.tsv`/`.xlsx`, `ListOfEmpiricalCompounds.tsv`, `userInputData.txt`, `userInput_to_EmpiricalCompounds.tsv`), `figures/` and `js/`. Result tables are **`.tsv`/`.xlsx`, never `.csv`**.
-   `<contrast>/v2/mummichog_manifest.tsv` and `<contrast>/v2/runner.log`.

Directly under `mummichog_pinned/`, `mummichog_engine.json` records the engine selected for the last enabled run (see [Backward compatibility](#backward-compatibility-configs-without-engine)).

Plus, directly under `mummichog_pinned/`, the presentation exports per contrast: `mummichog_pathway_bubble_<contrast>.{png,pdf}` (the bubble plot) and `mummichog_pathway_table_<contrast>.{tsv,csv}` (the sorted pathway table), the GSEA exports described below, the `mummichog_evidence_*_<contrast>.{tsv,csv}` supporting-evidence tables, `mummichog_summary_<contrast>.{tsv,csv}` (the Summary table, described below), `mummichog_results_all_contrasts.xlsx` (every contrast in one workbook, described below), and `contrasts.tsv`, which maps each sanitised subdirectory name back to its original DE contrast label (so the report can show real contrast names).

These exports are written by their own target (`metab_mummichog_exports`), not by the report. They are produced wherever the stage runs — **single-omics and multiomics alike** — and they do not need `rmarkdown`: a run without it loses the HTML report and nothing else. Before this was split out (issue #306) the supporting evidence and the GSEA were computed inside the report module, so both vanished in multiomics mode and whenever `rmarkdown` was missing.

To map pathways back to your feature ids, `join_features_to_results()` uses the feature id mummichog echoes into its own tables (via the 5th input column) — not the fragile post-de-duplication row numbers.

## The report section

The section opens with **Methodology**: a short description, the caveats a reader must not miss (what enters the engine, raw versus adjusted p-values, counts are not additive), a per-contrast table of what entered the engine (from each contrast's run record), and blue **ⓘ** panels — native `<details>` elements, no JavaScript — for the secondary detail: input and feature selection, ORA and significance, GSEA, empirical compounds and annotation agreement, parameters. With more than one contrast, an **Across contrasts** tab follows (see below). Then, per contrast, the tabs are:

1.  **Summary** — the top 20 pathways (all of them if fewer), ranked by the smaller of the raw ORA and GSEA p-values (a missing value counts as `Inf`), with a one-line count of what entered the engine. Rows need not be significant; when none is nominal, the tab says so. Columns: Pathway; **ORA** (status, raw p — FET for the MetaboAnalyst-compatible engine — and its Holm-adjusted p; a pathway with no significant EC reads *No significant EC overlap*, never *Not tested*); **GSEA** (status, raw p and BH FDR); **ORA overlap ECs** (`Sig.EC`) and **Features in overlap ECs** (distinct features, including non-significant ones and ones without their own pathway-candidate match), both only for pathways nominal in ORA — never for a GSEA-only pathway; **Annotation agreement (ECs)**; **Identified-metabolite enrichment**; **Interpretation** (*Also supported by identified-metabolite enrichment*, *Mummichog-only hypothesis*, *Not comparable*, *No nominal Mummichog signal*). It is the first rows of the contrast's overview (`build_mummichog_overview()`, `06t_mummichog_presentation.R`), which is built on `build_mummichog_summary()` and relabels it.
2.  **ORA plot** — the bubble plot with a one-line count of nominal pathways.
3.  **GSEA plot** — the scatter with a one-line count of nominal pathways and how many survive BH correction.
4.  **Pathway Evidence** — Pathway → EmpiricalCompound → measured features for the Summary pathways that are nominal in ORA and have evidence (a GSEA-only pathway's ORA evidence row stays under Advanced, so it is never read as ORA support). Every EC's feature table, whatever the EC's agreement (Match, Conflict, Mixed or Not assessed), shows each feature's original annotation, its annotation confidence (when the dataset declares any) and its agreement.
5.  **Advanced / Downloads** — every pathway (the overview, with numeric p-values), the full ORA table with every engine column named for what it holds (including `Hits.*`, putative compound matches), the full GSEA table, the full evidence tables at every grain, and a preview of the feature-level table (the exported file is complete).

**Shiny contract.** The report's Summary and full ORA table are presentation versions. The section objects the Shiny payload carries — `summary` (the top 20, with the column names and labels documented in `docs/SHINY_PAYLOAD_HANDOFF.md`), `table` and `gsea` — and the `mummichog_summary_<contrast>` and `mummichog_pathway_table_<contrast>` exports are unchanged, so their wording (*Not tested*, *Confirms identified signal*, `n ECs / m features`) differs from the report's until the app adopts the new presentation.

### Across contrasts

A tab before the contrast tabs, and the first data sheet of the workbook (`build_mummichog_across_contrasts()`, `06u_mummichog_across.R`). It is a lookup of each contrast's existing results, aligned on the model's pathway id (the name only where a pathway has no id). Per contrast, a cell reads *ORA Holm* or *GSEA FDR* (adjusted p below 0.05: Holm for ORA, BH for GSEA), *ORA nominal* or *GSEA nominal* (raw p below `p_cutoff` only), *—* (a result without nominal evidence) or *No result*. The report lists pathways with nominal evidence in at least two contrasts or adjusted significance in at least one (`mmc_across_report_rows()`); the workbook lists every pathway, with each contrast's four p-values. No p-value is combined across contrasts and no NES direction is shown: a pathway that recurs is not shown to change in the same direction, nor to be supported by the same ECs or features — compare them in the evidence tables.

### The mummichog workbook

`mummichog_pinned/mummichog_results_all_contrasts.xlsx` (`06v_mummichog_workbook.R`), written by the exports target beside the CSV/TSV files, which are all still written. Sheets: **README** (guide, definitions, caveats, and per contrast what entered the engine and whether its evidence was built), **Across Contrasts**, **Pathway Overview** (every pathway of every contrast, not only the report's top 20), **ORA Results** (including the Holm columns and `Hits.*`), **GSEA Results** (raw p and BH FDR), **EC Evidence** and **Feature Evidence**. In Feature Evidence, each feature's own pathway candidates (`Pathway-matched metabolite candidate ID(s)`, the evidence's `Pathway candidate(s) matched by feature`) and their model names in the same order (`Candidate metabolite name(s) (putative)`; `(no name in model)` where the model has none) sit just before `Original annotation`, `Annotation ID`, `Annotation confidence` and `Agreement` — names of putative matches, not identifications. Every data sheet has a `Contrast` column, a frozen header and a filter; p-values are numbers and identifiers are text. Every sheet is built from the same section objects as the report. It never reads or writes the differential-abundance workbook (`Final_results_*.xlsx`). It follows the same retention rule as every other output under `mummichog_pinned/`: a later run with mummichog disabled neither rewrites nor removes it (see [Backward compatibility](#backward-compatibility-configs-without-engine)). The HTML report names the workbook only when the file is present beside it, and says it comes from the exports step: the report does not depend on that step, so a file on disk does not show that the current run wrote it.

The Summary compares a Mummichog pathway with enrichment on identified metabolites (ORA, QEA and GSEA; not ssGSEA) through the **exact pathway id** both sides share — the model's pathway `id` and the GMT id, e.g. `mng00010`. Never by name. Three questions are kept apart:

| | means |
| --- | --- |
| joinable | the model gives the pathway a single id, and that id is in the identified-enrichment GMT |
| tested | the id appears in at least one of this contrast's identified ORA / QEA / GSEA result tables |
| supported | at least one tested method has FDR < 0.05 |

Identified-metabolite support reads `Not comparable` (not joinable), `Not tested` (joinable, not tested), the supporting methods such as `ORA + QEA`, or `None` (tested, not supported). Mummichog contribution reads `No nominal Mummichog signal` (neither ORA nor GSEA below `p_cutoff` — a ranked row, never a candidate), then for nominal rows `Not comparable`, `Confirms identified signal` (supported) or `Mummichog-only candidate` (otherwise) — a hypothesis, not a confirmed finding. Nothing here computes a statistic: it reads the existing results (`build_mummichog_summary()`, `R/domain/metabolomics/06h_mummichog_summary.R`).

## Complementary GSEA (MS peaks-to-pathways)

Alongside the ORA above, the pipeline runs a GSEA over the same mapped EmpiricalCompound universe, an analysis MetaboAnalyst introduced for MS peaks. It is **additive computationally**: the mummichog engine, its ORA statistics and its outputs are untouched, and the ORA bubble plot and table are still produced and still exported.

### GSEA method

`gsea_method` under the `mummichog:` block chooses how the EC scores and pathway EC sets are tested. Both engines use it the same way, and both methods are given the same EC scores, the same ranking statistic and the same pathway sets.

| `gsea_method` | What runs |
| --- | --- |
| unset, or `fgsea` (default) | `fgsea::fgsea()`, fgsea's adaptive multilevel method, on the named, signed EC-score vector (`mmc_gsea_fgsea()`, through the pipeline's shared `run_fgsea_core()`: the same call, with the same `nPermSimple = 10000`, as the GSEA of the other omics). The enrichment score, p-value, NES, BH adjustment and leading edge are fgsea's. |
| `metaboanalyst_legacy` | The reproduction of MetaboAnalystR's vendored peaks GSEA (`mmc_gsea_metaboanalyst()`). Kept only to reproduce results reported with it. It warns when used, and the report says so. |

Any other value stops the run, including when `gsea_enabled` is `false`: the key is validated whenever mummichog runs. So does `gsea_permutations` together with the default method (see [P-value estimation](#p-value-estimation)). The method that ran is recorded in the result, the report, the workbook guide and, for the MetaboAnalyst-compatible engine, `ma_run_info.json`.

**Why the default is not MetaboAnalyst's GSEA.** MetaboAnalystR's implementation has two defects, reported upstream as [xia-lab/MetaboAnalystR#381](https://github.com/xia-lab/MetaboAnalystR/issues/381):

1. It takes the absolute value of the ranked scores and then re-sorts them, while pathway positions still refer to the signed order. A hit is therefore weighted by a score that is not its own, and the leading edge is read back through the wrong names.
2. ECs with equal scores are collapsed into one ranked position. A pathway holding two of them cannot be scored and gets ES = 0 and p = 1.

**What changes with the default method, beyond those two points.**

- **The permutation universe.** Every EC has its own ranking position, including ECs with equal scores, so permutations are drawn over ECs. The legacy method ranks distinct scores, so its universe is smaller whenever scores are tied. The two results are not comparable pathway for pathway.
- **Tied scores.** ECs with equal scores are ordered by EC id, in a locale-independent order, and no score is altered. The result reports how many ECs are tied and their proportion, and names the pathways whose enrichment score can depend on that order: those for which a group of equal-scored ECs holds both members and non-members. A warning is given when there is at least one such pathway. There is no threshold.
- **Missing p-values.** fgsea gives no p-value when fewer than 10 of its 10,000 initial permutations produced an enrichment score of the observed sign. Such a pathway is kept in the table with empty NES, p and `padj`, counted in the result, and shown as `Tested; p-value unavailable`, never as `Not tested`. It has no `padj`. `padj` is fgsea's own BH adjustment and is passed through unchanged: fgsea adjusts over the pathways that have a p-value, except when it settles every p-value in its initial permutation step, where the adjustment was made before those p-values were withheld and so still counts them (slightly conservative for the others).
- **A pathway holding every ranked EC** is not tested (fgsea caps the size at one less than the list).
- **NES sign** is the ordinary GSEA one: which end of the ranked list the pathway's ECs are concentrated at. An EC takes the largest signed score among its features, so it is still not a per-metabolite direction of change.
- **The table schema is unchanged.** `ES defaulted by reference implementation` and `ES fallback reason` stay, as `FALSE` and empty, because the default method has no fallback.

#### P-value estimation

The default method estimates p-values with fgsea's adaptive multilevel scheme rather than from a fixed number of permutations, so a p-value is not bounded below by one over a permutation count. `nPermSimple` — the permutations of its initial step — is fixed at 10,000, the value every other GSEA in the pipeline uses, and is not a config key.

- **`gsea_permutations` is the legacy method's setting only.** With the default method the key has nothing to set. A config that carries it stops the run with a message, rather than having the value ignored or read as something else; remove the key, or set `gsea_method: metaboanalyst_legacy` to run the analysis the key was written for.
- **Reproducibility.** The call runs serially under `gsea_seed`, so the same inputs and seed give the same result.
- **`log2err`.** fgsea reports the expected error of each log2 p-value. The per-pathway values stay in the engine result; the report (Advanced / Downloads) and `ma_run_info.json` carry their maximum and the number of pathways that have a p-value but no error estimate. fgsea withholds the estimate for a p-value it flags as likely overestimated, and for one below its floor of 1e-50, which it reports as 1e-50. The GSEA table itself has the same columns as before.

#### Earlier results

`gsea_method: fgsea` has named two calculations. For a short period the default ran `fgsea::fgseaSimple()` with a fixed number of permutations (100 unless `gsea_permutations` said otherwise); it now runs `fgsea::fgsea()` as above. The earlier calculation can no longer be run. Results it produced are not relabeled: a result or run record states its own method label, `reference.method` (`fgsea::fgseaSimple`) and `params.n_perm`, and the report and workbook describe it from those stored fields. A multilevel result stores `reference.method` `fgsea::fgsea (fgseaMultilevel)` and `params.nperm_simple` instead, and no `n_perm`.

In the HTML report each contrast has five tabs: **Summary → ORA plot → GSEA plot → Pathway Evidence → Advanced / Downloads** (see [The report section](#the-report-section)). The ORA bubble and the GSEA scatter each keep their own tab under their own name; the ORA bubble is never shown under a GSEA heading. A contrast without GSEA gets the explanation in its GSEA plot tab instead of a plot — see [When GSEA has no results](#when-gsea-has-no-results).

Where ORA tests the features that clear the significance cutoff, GSEA uses the **complete ranked EmpiricalCompound list with no cutoff**, so it can reach pathways whose ECs shift coherently without any single feature reaching significance.

### The ranking statistic is chosen by DE method

The metabolomics DE stage writes different quantities into one generic `statistic` column, so the column name alone cannot say what is in it:

| `de_res$method` | what `statistic` holds | GSEA ranks on |
| --- | --- | --- |
| `limma` | `topTable()$t` — signed moderated t | `statistic` (or a literal `t`) |
| `t_test`, `t_test_equal` | `t.test()$statistic` — signed t | `statistic` (or a literal `t`) |
| `wilcoxon` | `wilcox.test()$statistic` — **unsigned rank-sum W** | `logFC` |
| `precomputed` | dropped by the loader | `logFC` |
| unrecognised or absent | unknown provenance | `logFC` |

A signed ranking is the premise of GSEA, and W is non-negative, so ranking on it would order EmpiricalCompounds by magnitude with no direction while the result claimed a t statistic. The report names the column it used and why.

### Requirements and switches

GSEA needs `fgsea` and a **readable** metabolic model — the built-in `human_mfn` model ships inside the Python package, so its pathway membership cannot be read and both GSEA and the supporting-evidence section are skipped with a message. Use `model_ref` or `model_json`.

`gsea_enabled: true` under the `mummichog:` block is the switch; it defaults to on wherever mummichog is enabled, so configs written before GSEA existed keep working. Three optional keys (`gsea_seed`, `gsea_min_size`, `gsea_max_size`) carry MetaboAnalystR's defaults, for both methods; a fourth, `gsea_permutations` (default 100), belongs to the legacy method only and stops the run if set with the default one. With `gsea_method: metaboanalyst_legacy`, overriding any of them is a deviation from the reference and the result records it as one; with the default method the three are plain settings, reported with the result. Note `gsea_seed` is **not** the project `seed:`: 123 is upstream's own value. `gsea_method` is described [above](#gsea-method).

### When GSEA has no results

The GSEA plot tab is rendered for **every** contrast, including ones with nothing to plot, so it always answers what GSEA found. What it says depends on which of three things happened:

| State | The GSEA plot tab shows |
| --- | --- |
| GSEA ran and the scatter was built | the scatter and its NES / -log10(p) legend |
| GSEA ran but the scatter could not be drawn | a sentence saying so, pointing at the full GSEA table under Advanced / Downloads, which is unaffected |
| GSEA produced nothing | the explanation below. No plot. |

That last message deliberately **does not assert a cause**. A missing GSEA result means either that the analysis could not run (no `fgsea`, no readable model, no DE table for the contrast, or `gsea_enabled: false`) *or* that it ran and no pathway passed the size filter — a legitimate negative result. The report cannot tell the two apart, so it names both and points at the run's console log, which can. Reporting a real negative as a misconfiguration would be worse than saying less.

In none of these states is the ORA bubble shown under a GSEA heading. ORA always lives in its own tab, under its own name, with its own legend.

### Outputs

Per contrast, under `mummichog_pinned/`: `mummichog_gsea_scatter_<contrast>.{png,pdf}` and `mummichog_gsea_table_<contrast>.{tsv,csv}` (the complete table; the HTML report caps its view). These four files are cleared by name before each write, so a rerun that legitimately comes back without GSEA cannot leave a stale table reading as current.

### What is verified, and what is not

**Default method (`fgsea`).** The statistic is fgsea's own: nothing in it is reimplemented here. The tests check what this pipeline adds around the call: that the result equals a direct `fgsea::fgsea(nPermSimple = 10000)` call on the same ordered scores under the same seed, that the enrichment score and leading edge equal `fgsea::calcGseaStat()` on the documented ordering of tied ECs, that a hand-calculated example with negative scores gives the expected values, that results are reproducible for a seed and independent of input order, that both engines call the same function, and that a result recorded by the earlier fixed-permutation default is read back as such. How the EC scores are built is unchanged and is covered by the point about ranked inputs below.

**Legacy method (`metaboanalyst_legacy`).** It reproduces the pinned `xia-lab/MetaboAnalystR` Peaks-to-Pathways implementation at **two different strengths**, and the result says which is which:

-   the **GSEA calculation given ranked EC inputs** is checked against a reference generated by running upstream's own code (`tests/testthat/fixtures/mummichog_gsea_parity/`);
-   the **construction of those ranked inputs** — notably how features sharing an m/z are merged — is matched by reading the pinned source and is covered by no fixture. Independent validation is tracked in issue #302.

For the legacy method, NES sign is retained for parity only. The reference transforms and reorders the score vector *after* pathway positions are built, so the sign is not a biological direction. For the default method the sign is the ordinary GSEA one (see [GSEA method](#gsea-method)). The report states the note that applies to the method that ran wherever NES appears.

## Stochasticity caveat

mummichog v2 estimates null distributions by random permutation, drawn from Python's `random` module. mummichog itself has no seed option, so the pipeline's launcher seeds `random` from `params.seed` (1 when unset) before mummichog starts. A rerun on the same input with the same seed therefore reproduces every pathway's p-value and overlap.

What can still differ between runs is **order**: mummichog writes some sets in memory-address order, so the order of rows with tied p-values, and of the IDs within an overlap cell, can change. Compare those columns as sets, not as ordered lists, in tests. `{targets}` only re-runs the stage when its inputs change, so this does not cause spurious rebuilds.
