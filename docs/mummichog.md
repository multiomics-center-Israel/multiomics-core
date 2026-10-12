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
        ionization_mode: pos_default   # pos_default | positive | negative
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

Plus, directly under `mummichog_pinned/`, the report's presentation exports per contrast: `mummichog_pathway_bubble_<contrast>.{png,pdf}` (the bubble plot) and `mummichog_pathway_table_<contrast>.{tsv,csv}` (the sorted pathway table), the GSEA exports described below, the `mummichog_evidence_*_<contrast>.{tsv,csv}` supporting-evidence tables, and `contrasts.tsv`, which maps each sanitised subdirectory name back to its original DE contrast label (so the report can show real contrast names).

To map pathways back to your feature ids, `join_features_to_results()` uses the feature id mummichog echoes into its own tables (via the 5th input column) — not the fragile post-de-duplication row numbers.

## Complementary GSEA (MS peaks-to-pathways)

Alongside the ORA above, the pipeline runs a MetaboAnalyst-style GSEA over the same mapped EmpiricalCompound universe. It is **additive**: the pinned mummichog engine, its statistics and its outputs are untouched, and the bubble plot remains the section's primary view.

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
