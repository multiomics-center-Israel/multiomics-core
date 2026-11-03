# Mummichog pathway analysis (pinned v2, isolated venv)

> Moved here from the top-level README; the stochasticity caveat has since been corrected to match the code. Back to the [docs index](README.md).

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

The python engine is untouched and stays the default.

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
| GSEA | on python mummichog's ECs (m/z-merged scores) | MetaboAnalyst's peaks GSEA on the **same joint-polarity ECs and pathway sets as the ORA**; EC score = max of its features' signed DE statistic; `padj` is **BH** |
| Pathway Evidence | yes | not produced (it reads the python engine's EC tables) |

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
- `ma_gsea_results.tsv`: the GSEA table in MetaboAnalyst's order (raw p, BH `padj`, NES, leading edge, each pathway's ECs).
- `ma_run_info.json`: parameters, model, Gamma fit, GSEA settings and provenance, and counts.

The report, its exports and the Shiny payload read these through the same targets as before. The report's ORA table adds the EASE, Gamma, empirical and Holm columns; its GSEA tab, scatter and GSEA table export show this engine's GSEA. Python GSEA never runs with this engine, so the two engines' results are never mixed in one report.

GSEA reproduces `.compute.mummichog.RT.fgsea()` at the pinned commit, with two documented exceptions that change no statistic: MetaboAnalyst's result CSV attaches EC hits in library order to rows sorted by p-value (here each row carries its own pathway's ECs), and a single tested pathway stays a row. Validated end to end against that commit on the synthetic fixtures (`generate_gsea_reference.R`).

Limits:

- Retention time is required. MetaboAnalyst's no-RT, compound-level path is not reproduced.
- In mixed mode, both polarities' retention times must be on one scale.
- MetaboAnalyst's custom-adduct and currency-customisation options are not offered. The default currency list is applied.

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

Plus, directly under `mummichog_pinned/`, the presentation exports per contrast: `mummichog_pathway_bubble_<contrast>.{png,pdf}` (the bubble plot) and `mummichog_pathway_table_<contrast>.{tsv,csv}` (the sorted pathway table), the GSEA exports described below, the `mummichog_evidence_*_<contrast>.{tsv,csv}` supporting-evidence tables, `mummichog_summary_<contrast>.{tsv,csv}` (the Summary table, described below), and `contrasts.tsv`, which maps each sanitised subdirectory name back to its original DE contrast label (so the report can show real contrast names).

These exports are written by their own target (`metab_mummichog_exports`), not by the report. They are produced wherever the stage runs — **single-omics and multiomics alike** — and they do not need `rmarkdown`: a run without it loses the HTML report and nothing else. Before this was split out (issue #306) the supporting evidence and the GSEA were computed inside the report module, so both vanished in multiomics mode and whenever `rmarkdown` was missing.

To map pathways back to your feature ids, `join_features_to_results()` uses the feature id mummichog echoes into its own tables (via the 5th input column) — not the fragile post-de-duplication row numbers.

## The report section

Per contrast, the tabs are:

1.  **Summary** — the top 20 pathways (all of them if fewer), ranked by the smaller of the raw Mummichog ORA and GSEA p-values (a missing value counts as `Inf`). Rows need not be significant: a pathway is *nominal* when ORA or GSEA is below `p_cutoff`, and when none is, the tab says so above the table. Columns: Pathway, Identified-metabolite support, Mummichog ORA (`p = …`; empirical raw p, no adjusted value exists), Mummichog GSEA (`p = …; FDR = …`), Mummichog contribution, Supporting evidence (`n ECs / m features`, distinct measured features, shown only for pathways nominal in ORA; `Not available — GSEA-only signal` for a GSEA-only pathway; `—` for a pathway nominal in neither), Annotation agreement (EC counts, zeros omitted). The `mummichog_summary_<contrast>` export holds the same rows; the complete ORA and GSEA results are under Advanced / Downloads.
2.  **ORA plot** — the bubble plot with a one-line count of nominal pathways.
3.  **GSEA plot** — the scatter with a one-line count of nominal pathways and how many survive BH correction.
4.  **Pathway Evidence** — Pathway → EmpiricalCompound → measured features for the Summary pathways that are nominal in ORA and have evidence (a GSEA-only pathway's ORA evidence row stays under Advanced, so it is never read as ORA support). A Mixed EC also shows each feature's annotation and agreement.
5.  **Advanced / Downloads** — the full ORA and GSEA tables (ES, NES, leading edge, ES-fallback flags, deviations), the full evidence tables at every grain, and a preview of the feature-level table (the exported file is complete).

The Summary compares a Mummichog pathway with enrichment on identified metabolites (ORA, QEA and GSEA; not ssGSEA) through the **exact pathway id** both sides share — the model's pathway `id` and the GMT id, e.g. `mng00010`. Never by name. Three questions are kept apart:

| | means |
| --- | --- |
| joinable | the model gives the pathway a single id, and that id is in the identified-enrichment GMT |
| tested | the id appears in at least one of this contrast's identified ORA / QEA / GSEA result tables |
| supported | at least one tested method has FDR < 0.05 |

Identified-metabolite support reads `Not comparable` (not joinable), `Not tested` (joinable, not tested), the supporting methods such as `ORA + QEA`, or `None` (tested, not supported). Mummichog contribution reads `No nominal Mummichog signal` (neither ORA nor GSEA below `p_cutoff` — a ranked row, never a candidate), then for nominal rows `Not comparable`, `Confirms identified signal` (supported) or `Mummichog-only candidate` (otherwise) — a hypothesis, not a confirmed finding. Nothing here computes a statistic: it reads the existing results (`build_mummichog_summary()`, `R/domain/metabolomics/06h_mummichog_summary.R`).

## Complementary GSEA (MS peaks-to-pathways)

Alongside the ORA above, the pipeline runs a MetaboAnalyst-style GSEA over the same mapped EmpiricalCompound universe. It is **additive computationally**: the pinned mummichog engine, its statistics and its outputs are untouched, and the ORA bubble plot and table are still produced and still exported.

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

`gsea_enabled: true` under the `mummichog:` block is the switch; it defaults to on wherever mummichog is enabled, so configs written before GSEA existed keep working. Four optional keys (`gsea_permutations`, `gsea_seed`, `gsea_min_size`, `gsea_max_size`) carry the pinned MetaboAnalystR defaults. Overriding any of them is a methodological deviation, and the result records it as one. Note `gsea_seed` is **not** the project `seed:`: 123 is upstream's own value.

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

The engine reproduces the pinned `xia-lab/MetaboAnalystR` Peaks-to-Pathways implementation at **two different strengths**, and the result says which is which:

-   the **GSEA calculation given ranked EC inputs** is checked against a reference generated by running upstream's own code (`tests/testthat/fixtures/mummichog_gsea_parity/`);
-   the **construction of those ranked inputs** — notably how features sharing an m/z are merged — is matched by reading the pinned source and is covered by no fixture. Independent validation is tracked in issue #302.

NES sign is retained for parity only. The reference transforms and reorders the score vector *after* pathway positions are built, so the sign is not a biological direction; the report states this wherever NES appears.

## Stochasticity caveat

mummichog v2 estimates null distributions by random permutation, drawn from Python's `random` module. mummichog itself has no seed option, so the pipeline's launcher seeds `random` from `params.seed` (1 when unset) before mummichog starts. A rerun on the same input with the same seed therefore reproduces every pathway's p-value and overlap.

What can still differ between runs is **order**: mummichog writes some sets in memory-address order, so the order of rows with tied p-values, and of the IDs within an overlap cell, can change. Compare those columns as sets, not as ordered lists, in tests. `{targets}` only re-runs the stage when its inputs change, so this does not cause spurious rebuilds.
