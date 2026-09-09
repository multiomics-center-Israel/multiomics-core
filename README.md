# multiomics-core

A modular and reproducible R framework for single-omics and multi-omics analyses (RNA-seq, proteomics, metabolomics, multi-omics integration).

### The project emphasizes

-   clear separation of concerns (I/O, preprocessing, DE, QC, plotting)
-   configuration-driven workflows with reproducibility via `renv`
-   scalable orchestration via `{targets}`

#### This project relies heavily on {targets} for reproducible pipeline orchestration.

For in-depth documentation and tutorials, see the official targets book: <https://books.ropensci.org/targets/>

------------------------------------------------------------------------

## What this repository provides

-   Standardized data loading and validation
-   Omics-specific preprocessing (filtering, normalization, imputation)
-   Proteomics differential expression via **limma** with multiple imputations and stability filtering
-   RNA-seq differential expression via **DESeq2**
-   Metabolomics preprocessing (missingness classification, MNAR/MAR imputation, TSS/Median/PQN normalization, LOESS drift correction) and DE
-   Pathway enrichment analysis (fGSEA, ORA, QEA, ssGSEA)
-   Unified QC utilities (PCA, heatmaps, sample distance)
-   A central YAML configuration file controlling all parameters
-   A `{targets}` pipeline for reproducible, dependency-aware execution

------------------------------------------------------------------------

## Getting started (new users)

If you are new to **multiomics-core**, start here:

-   📘 **Onboarding guide:** `docs/onboarding.md`
-   📘 **Developer guide:** `docs/developer_guide.md`

The onboarding guide explains:

-   How to open the project in RStudio
-   How to restore the R environment with `renv`
-   How to run analyses interactively or via `{targets}`
-   How to reproduce previous runs

------------------------------------------------------------------------

## Repository structure

```
R/
├── core/         # Generic utilities (I/O, validation, QC, clustering, enrichment, plotting)
├── domain/       # Omics-specific logic (rnaseq, proteomics, metabolomics, multiomics)
├── modules/      # Pipeline steps (wrappers for domain logic)
├── pipeline/     # {targets} pipeline orchestration
├── services/     # External integrations (AI commentary)
config/
├── templates/   # Analysis config templates (rna, proteins, metabolomics, multiomics)
data/            # Example datasets and reference files
docs/            # Onboarding, developer guide, ADRs, migration notes
tests/           # testthat tests
_targets.R       # {targets} pipeline definition
run.R            # CLI entrypoint / wizard launcher
renv.lock        # Locked dependency versions
```

------------------------------------------------------------------------

## Requirements

-   **R ≥ 4.3** (tested with R 4.5.x)
-   **RStudio** (recommended)
-   **`renv`** (for reproducible environments)

### Windows users (IMPORTANT)

On Windows, some CRAN / Bioconductor packages may need to be **compiled from source** (e.g. `SparseArray`, `IRanges`, `Biobase`). Therefore, **Rtools is required**.

**Install Rtools (matching your R version):** 👉 <https://cran.r-project.org/bin/windows/Rtools/>

After installation, **restart R / RStudio**, then verify:

``` r
Sys.which("gcc")
Sys.which("make")
```

Both commands should return a valid path. If they return `""`, Rtools is not correctly installed or not on `PATH`.

Missing Rtools may cause `renv::restore()` to fail.

### Linux / macOS users

-   A standard compiler toolchain is required
-   System libraries commonly needed by Bioconductor (e.g. `libxml2`, `curl`, `openssl`)

------------------------------------------------------------------------

## Package repositories (CRAN + Bioconductor)

This project relies on **CRAN** and **Bioconductor** packages.

We recommend using **Posit Package Manager (PPM)** for CRAN together with standard Bioconductor repositories.

Recommended setup:

``` r
options(repos = c(
  CRAN = "https://packagemanager.posit.co/cran/latest"
))

BiocManager::repositories()
```

The Bioconductor version must match the one recorded in `renv.lock` (e.g. `Bioconductor 3.22`).

To check:

``` r
BiocManager::version()
```

If needed:

``` r
BiocManager::install(version = "3.22")
```

------------------------------------------------------------------------

## Setup

Clone the repository and restore the R environment:

``` r
install.packages("renv")
renv::restore()
```

------------------------------------------------------------------------

## Configuration

The pipeline reads its config path from the **`MULTIOMICS_CONFIG`** environment variable.
If the variable is not set, it defaults to `config.yaml` in the project root.

### Setting up your config path

1.  Copy the example environment file:

``` bash
cp .Renviron.example .Renviron
```

2.  Edit `.Renviron` and set the path to your YAML config:

```
MULTIOMICS_CONFIG=/path/to/your/config.yaml
```

3.  Restart your R session (`.Renviron` is loaded on startup).

> **Note:** `.Renviron` is git-ignored so each collaborator can point to their own config without modifying tracked files.

### Creating a config file

Start from an existing template:

``` bash
cp config/templates/proteins_config.yaml config/<PROJECT>_<ROUND>.yaml
```

The configuration controls:

-   input and output file paths
-   omics-specific parameters
-   filtering, normalization, and imputation settings
-   differential expression thresholds
-   QC aesthetics (color, shape, sample ID columns)

------------------------------------------------------------------------

## Running the pipeline (via `{targets}`)

The recommended way to run analyses is via `{targets}`.

From an R session in the project root:

``` r
library(targets)
tar_make()
```

`tar_make()` runs whichever modes are enabled in your config (proteomics, RNA-seq, metabolomics). Each mode executes its own DAG covering configuration validation, input loading, preprocessing, differential expression, QC, and output generation.

To run a single mode:

``` r
tar_make(names = starts_with("prot_"))  # proteomics only
tar_make(names = starts_with("rna_"))   # RNA-seq only
tar_make(names = starts_with("met"))    # metabolomics only (matches met_* and metab_*)
```

`{targets}` ensures that only steps affected by changes are recomputed.

### Learning more about `{targets}`

This project uses `{targets}` for reproducible, dependency-aware pipeline orchestration.

For a detailed introduction, tutorials, and best practices, see the official **targets** book: <https://books.ropensci.org/targets/>

------------------------------------------------------------------------

## Mummichog pathway analysis (pinned v2, isolated venv)

The metabolomics mode runs [mummichog](http://mummichog.org) for m/z-based pathway/network enrichment via a **version-pinned, isolated engine** (`R/domain/metabolomics/06c_mummichog_pinned.R`): `mummichog==2.7.0` invoked as a `{processx}` subprocess in a dedicated venv, depending only on light R packages (`readr`, `processx`, `jsonlite`) — no Bioconductor. It runs on mummichog's built-in `human_mfn` model by default.

### One-time setup

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

### How to run

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

### Choosing a metabolic model

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

### Where outputs land

Under `<metab_out_dir>/mummichog_pinned/`, one subdirectory per contrast (`<contrast>/`, the contrast name sanitised to `A-Za-z0-9_`):

-   `<contrast>/input.tsv` and `<contrast>/input.tsv.idmap.tsv` — the exact table sent to mummichog for that contrast (m/z, retention time, p-value, statistic, **feature\_id as the 5th column**) plus a provenance id-map.
-   `<contrast>/v2/<timestamp>.<project>/` — the mummichog result tree: `result.html`, `tables/` (`mcg_pathwayanalysis_*.tsv`/`.xlsx`, `mcg_modularanalysis_*.tsv`/`.xlsx`, `ListOfEmpiricalCompounds.tsv`, `userInputData.txt`, `userInput_to_EmpiricalCompounds.tsv`), `figures/` and `js/`. Result tables are **`.tsv`/`.xlsx`, never `.csv`**.
-   `<contrast>/v2/mummichog_manifest.tsv` and `<contrast>/v2/runner.log`.

Plus, directly under `mummichog_pinned/`, the report's presentation exports per contrast:

-   `mummichog_pathway_bubble_<contrast>.{png,pdf}` and `mummichog_pathway_table_<contrast>.{tsv,csv}` — the ORA bubble plot and the sorted ORA pathway table.
-   `mummichog_gsea_scatter_<contrast>.{png,pdf}` and `mummichog_gsea_table_<contrast>.{tsv,csv}` — the complementary GSEA summary scatter and its result table (see below).
-   `mummichog_evidence_pathways_<contrast>.{tsv,csv}`, `mummichog_evidence_empirical_compounds_<contrast>.{tsv,csv}` and `mummichog_evidence_features_<contrast>.{tsv,csv}` — the pathway supporting-evidence tables (see below).
-   `contrasts.tsv`, which maps each sanitised subdirectory name back to its original DE contrast label (so the report can show real contrast names).

The GSEA and evidence exports are written only when they could be produced; the ORA exports are always written.

To map pathways back to your feature ids, `join_features_to_results()` uses the feature id mummichog echoes into its own tables (via the 5th input column) — not the fragile post-de-duplication row numbers.

### Complementary GSEA (MS peaks-to-pathways)

Alongside the pinned 2.7.0 **ORA**, the pipeline runs a MetaboAnalyst-style **GSEA** on the same mapped feature / EmpiricalCompound / pathway universe (`R/domain/metabolomics/06g_mummichog_gsea.R`). The mummichog engine and its statistics are untouched — this is an additional analysis, not a replacement — and the two answer different questions:

| | Mummichog ORA | GSEA |
|---|---|---|
| Features used | only those passing `p_cutoff` | the **complete ranked list**, no cutoff |
| Unit | significant EmpiricalCompounds | all detected EmpiricalCompounds |
| Reports | overlap / detected pathway size, empirical permutation p-value | ES, NES, raw p-value, BH-adjusted p-value, leading-edge ECs |
| Role in the report | complementary evidence table | the **primary** pathway plot |

#### Behavioural reference

The implementation reproduces a **pinned** upstream version, not a paraphrase of it:

| | |
|---|---|
| Repo | `xia-lab/MetaboAnalystR` |
| Commit | `398476ae2a0c996e390925fa62adc46b5ecde334` |
| Files | `R/util_fgsea.R`, `R/peaks_to_function.R` |
| Path | `PerformPSEA` → `.init.RT.Permutations` → `.compute.mummichog.RT.fgsea` → `fgsea2` → `my.fgsea` → `.run_fgsea_inner` |

Both files are byte-identical between that commit and upstream's current default-branch head, so the pin is the live Peaks-to-Pathways code path. A **parity test** (`tests/testthat/test-mummichog-gsea-parity.R`) compares our engine against a fixture produced by running the pinned `.run_fgsea_inner()` itself; see `tests/testthat/fixtures/mummichog_gsea_parity/REFERENCE.md`.

Because the engine calls fgsea **internals**, the fgsea build is pinned too.
The generator refuses to run unless the installed fgsea matches `renv.lock`
exactly (a patch-level difference needs an explicit opt-in that is then recorded
in the fixture, and a different `x.y` series is refused outright), the parity
test **skips rather than compares** across series, and a further test fails if
the suite runs somewhere the locked build is installed while the fixture was
generated on something else. The committed fixture has now been regenerated
under the exact locked build, **fgsea 1.36.2**, with
`fgsea_exact_match = TRUE`; its reference results are identical to the previous
1.36.0 fixture. See
`tests/testthat/fixtures/mummichog_gsea_parity/REFERENCE.md` for provenance and
validation details.

#### Semantics

-   **Ranking statistic** — the moderated `t` statistic (our per-contrast limma tables carry it as `statistic`), which is what MetaboAnalyst ranks on; `logFC` is used only when a DE table has no usable statistic. **The ORA input contract is unchanged** — it still sends `logFC` as mummichog's `statistic` column. Whichever metric GSEA used is printed in the report and stored in the result.
-   **EmpiricalCompound score** — features sharing one m/z are merged first (mean, MetaboAnalyst's default), then an EC takes the **signed maximum** of its member features' scores (`ec.exp.vec <- unlist(lapply(ec_exp_dict, max))` in the retention-time/v2 branch). Not the mean, not `max(abs())`, not the most significant feature.
-   **Gene sets** — pathway → its *detected* EmpiricalCompounds.
-   **Ranked inputs** — ECs sharing a score are collapsed into one ranked position named by their `"; "`-joined ids, and every EC maps to that position's index. This is why the ranked vector can be shorter than the EC universe, and why `Tested size` can differ from `Detected ECs`.
-   **Engine** — `fgsea::calcGseaStat` for the enrichment score and `fgsea:::calcGseaStatCumulativeBatch` for the permutation tallies, orchestrated exactly as upstream does. `fgsea::fgseaSimple()` **cannot** stand in: it takes one per-item vector rather than the stats/ranks pair, de-duplicates pathway positions, and applies `abs()` *after* sorting rather than before. `mmc_gsea_metaboanalyst()` is the smallest local helper that reproduces the upstream orchestration.
-   **Defaults** (all MetaboAnalyst's) — `permNum = 100`, `gseaParam = 1`, `minSize = 1`, `maxSize = Inf`, `set.seed(123)` before the per-batch permutation seeds. Overriding any of them is recorded as a methodological deviation in the result and printed in the report.
-   **Multiple testing** — `p.adjust(..., "fdr")` over exactly the pathways that passed the size filter, i.e. the tested universe. GSEA p-values are never mixed with the ORA's empirical p-values.

#### Reading NES

`ES`/`NES` keep their sign — but **only for parity**, and the sign is not an interpretable direction. Two things in the reference implementation independently break that reading:

1.  the ranked vector is transformed to |score| (`stats <- abs(stats)^gseaParam`), so it carries no direction; and
2.  pathway positions are constructed from the **signed-score ordering** *before* that transform and the subsequent magnitude re-sort, so a pathway's scored indices need not correspond to where its ECs actually sit in the ranking that gets scored.

So the sanctioned wording is deliberately weak:

> NES sign is retained to reproduce the pinned MetaboAnalystR result. Because the reference implementation transforms and reorders the score vector after pathway positions are constructed, NES sign should not be interpreted as biological up/down direction or as a direct statement about the magnitude of the pathway members' original scores.

The report, tables, plot and exports carry that wording and nothing stronger; a regression test keeps up/down, numerator/denominator and high/low-|score|-end phrasing from coming back.

#### Two upstream behaviours we reproduce rather than "fix"

1.  **A failed enrichment score becomes `ES = 0`.** Upstream's `tryCatch` swallows *any* `calcGseaStat()` failure into `ES = 0` with an empty leading edge. The common trigger is tied EC scores: pathway members become tie-group *indices* and are never de-duplicated, so a pathway holding two ECs with the same score passes a non-strictly-increasing `selectedStats`, which `calcGseaStat()` rejects. We reproduce the numbers exactly **and record why** — the results table carries a generic `ES defaulted by reference implementation` flag plus an `ES fallback reason` (classified from the actual condition: duplicate ranked positions, whole-ranking selection, empty/missing members, or the verbatim error prefixed *unexpected*), an unexpected cause also raises an R warning, and the report lists the observed reasons. A placeholder `0` is never read as a measured score, and a non-tie failure is never described as a tie.
2.  **Leading edges are often empty.** Because `abs()` is applied before the sort, pathway indices taken from the signed ordering address a differently-ordered vector; the leading edge is then intersected with the pathway's own ECs and frequently comes back empty.

#### Config

Nothing needs to be set — the defaults above are the reference values. All four keys are optional, and setting any of them is reported as a deviation:

``` yaml
modes:
  metabolomics:
    enrichment:
      mummichog:
        # gsea_permutations: 100    # PerformPSEA(permNum = 100)
        # gsea_seed: 123            # .run_fgsea_inner: set.seed(123)
        # gsea_min_size: 1
        # gsea_max_size: .inf
        # gsea_param: 1
```

GSEA needs `fgsea`, a **readable metabolic model** (`model_ref` or `model_json` — see the caveat below) and a signed DE statistic. When any of those is missing it is skipped with a message and the report falls back to showing the ORA bubble plot.

### Pathway supporting evidence

For every enriched pathway the report can trace the chain

    Pathway -> EmpiricalCompound -> pathway-matching candidate -> measured feature -> original annotation -> agreement

so a biologist can see *which measured features supported the pathway* and *whether the identity mummichog used to place them there agrees with the dataset's own annotation* (`R/domain/metabolomics/06f_mummichog_evidence.R`).

-   **Pathway-matching candidate(s)** — all candidate compounds of the EmpiricalCompound intersected with the compounds of *that* pathway. Every surviving candidate is kept; none is arbitrarily picked. This deliberately does **not** use mummichog's `face_compound` / "Best guess": verified in the 2.7.0 sources, `designate_face_cpd()` picks `chosen_compounds[-1]` ("arbitrarily designated" per its own docstring) and `collect_hit_Trios()` fills `chosen_compounds` from the union of *all* significant pathways — which is also what the pathway table's `overlap_features (id)` column contains. Neither is a per-pathway, evidence-ranked identification. A pathway-matching candidate is a **putative** identity: the identity *through which* this EC maps to this pathway.
-   **Measured features** — every input signal forming the EC is listed with its feature id, m/z, retention time, adduct/ion, p-value, statistic and significance flag. Nothing is summed or collapsed into a representative feature.
-   **Original annotation** — normalised into a schema-agnostic contract (`original_annotation_name`, `original_annotation_id`, `original_annotation_confidence`), so datasets with MSI identification levels, datasets with annotations but no level system, and unannotated features are all handled. Nothing is hard-coded to "Level 1", and missing information stays `NA`.
-   **Agreement** — kept strictly separate from annotation *confidence*, and compared on stable compound ids (KEGG preferred) when both sides have one, otherwise by a conservative normalised-name/synonym comparison against the model's own `";"`-separated synonym lists. **Never** matched on m/z, molecular formula or mass. Conflicts are reported, never removed — this is an evidence layer, not a re-analysis of the enrichment.

    Each measured feature gets `Match`, `Conflict` or `Not assessed`. An EmpiricalCompound rolls its features up into **four** states, so a disagreement is never hidden behind an agreement:

    | features | EC verdict |
    |---|---|
    | some Match, no Conflict | `Match` |
    | some Conflict, no Match | `Conflict` |
    | both Match and Conflict | `Mixed` |
    | neither | `Not assessed` |

    Features with no usable annotation never override assessed evidence (`Match + Not assessed` → `Match`, `Conflict + Not assessed` → `Conflict`).

-   **Counts come at two grains and are labelled as such.** In the pathway summary, `ECs Match/Conflict/Mixed/Not assessed` count EmpiricalCompounds by their roll-up state; `features Match/Conflict/Not assessed` count measured features by their own verdict. On an EC row, `n_match`/`n_conflict`/`n_not_assessed` are that EC's feature-level evidence and sum to its `# Features`. None of these is the pathway overlap, the detected pathway size, or the number of candidate compounds.

> **Caveat — pathway membership needs a readable model.** mummichog does not export pathway → compound membership (in `reporting.py` the `'all_compounds': P.cpds` line is commented out), so it is read from the metabolic model the stage ran on, resolved through the same `mmc_select_model()`. The built-in **`human_mfn` model lives inside the Python package** and is not a JSON file R can read, so with that model both the GSEA and the evidence layers report themselves unavailable and their report sections are omitted rather than guessed at. Configure `model_ref` or `model_json` (Azimuth or mummichog-2 native JSON) to enable them.

In the HTML report each contrast's mummichog tab becomes **Pathway Plot | ORA Results Table | GSEA Results Table | Supporting Evidence**.

The **Pathway Plot** is the MetaboAnalyst-style GSEA scatter: x = NES, y = −log10(raw GSEA p), diverging colour on NES centred at 0, point size = √−log10 p (MetaboAnalyst's own `radi.vec` mapping, which re-encodes the significance already on the y axis rather than adding a quantity — the legend says so and never exposes the upstream variable name), a NES = 0 vertical reference, a p = 0.05 horizontal reference, and `ggrepel` labels on the most significant pathways. It is the **primary and only** pathway plot whenever GSEA ran; the ORA bubble plot is built solely as the fallback for when it could not, so the two never compete as rival primary visualisations. The ORA analysis stays fully available as its own results table plus the supporting-evidence drill-down.

The **Supporting Evidence** tab shows a pathway summary table plus a collapsible per-pathway drill-down (capped at the 25 most significant pathways in the HTML; the full tables are always exported as `mummichog_evidence_*`).

### Stochasticity caveat

mummichog v2 estimates null distributions by **random permutation with no seed control**, so p-values and rankings vary slightly between runs on identical input. `{targets}` only re-runs the stage when its inputs change, so this doesn't cause spurious rebuilds — but do **not** expect bit-identical reruns, and don't assert exact equality in tests.

------------------------------------------------------------------------

## Running preprocessing interactively (example)

For exploratory work or debugging:

``` r
# Load functions in dependency order
# 1. Core utilities
invisible(lapply(list.files("R/core", full.names = TRUE, recursive = TRUE), source))
# 2. Services
invisible(lapply(list.files("R/services", full.names = TRUE, recursive = TRUE), source))
# 3. Domain logic
invisible(lapply(list.files("R/domain", full.names = TRUE, recursive = TRUE), source))
# 4. Modules
invisible(lapply(list.files("R/modules", full.names = TRUE, recursive = TRUE), source))

# Load config
config <- load_config("config/config.yaml")

# --- Proteomics ---
inputs <- load_proteomics_inputs(config)
res    <- preprocess_proteomics(inputs, config)

# --- RNA-seq ---
inputs <- load_rna_inputs(config)
res    <- preprocess_rna(inputs, config)

# --- Metabolomics ---
inputs <- load_metabolomics_inputs(config)
res    <- preprocess_metabolomics(inputs, config)

# Example QC: PCA
qc_pca_scatter(
  expr_mat = res$expr_work,
  meta     = res$meta,
  cfg      = config$modes$proteomics,
  out_file = "outputs/proteomics/qc/pca_pc1_pc2.png"
)
```

------------------------------------------------------------------------

## Reproducibility

-   All package versions are locked in `renv.lock`
-   Outputs and caches are excluded from git
-   `{targets}` provides deterministic, restartable pipelines
-   Each run records execution metadata (config snapshot, git commit, session info)

------------------------------------------------------------------------

## Outputs

All analysis outputs are written under `<project.dir>/<paths.out>/Results_<project.name>_<analysis_round>/<mode>/`, where path components come from your config YAML (defaults: `paths.out: "outputs"`).

-   The project directory lives outside the repository (set via `project.dir` in your config)
-   Results should be shared by zipping the relevant output folder
-   Each run is isolated by its configuration parameters

------------------------------------------------------------------------

## Developer notes

If you want to extend, modify, or maintain **multiomics-core**, see:

-   📘 **Developer guide:** `docs/developer_guide.md`

------------------------------------------------------------------------

## Acknowledgments

-   AI-powered figure commentary uses scientific domain knowledge templates informed by [K-Dense AI claude-scientific-skills](https://github.com/K-Dense-AI/claude-scientific-skills) (MIT License). To enable AI commentary, clone their repository into the project root and configure credentials for your chosen backend (see below).
-   Commentary generation is powered by [Claude](https://www.anthropic.com/) (Anthropic) or [GPT-4o](https://openai.com/) (OpenAI).

### AI commentary setup (optional)

Each backend has its own prerequisite. No credentials are stored in or shared via this repository.

-   `claude-code` (default): the [`claude` CLI](https://claude.ai/claude-code) must be installed, on your `PATH`, and authenticated (e.g. via `claude login`). The pipeline only checks that the CLI is on `PATH` — authentication is handled by the CLI itself.
-   `claude`: set `ANTHROPIC_API_KEY="your-key-here"` in your environment (or `.Renviron`).
-   `openai`: set `OPENAI_API_KEY="your-key-here"` in your environment (or `.Renviron`).

Enable in your config YAML:

```yaml
commentary:
  enabled: true
  backend: "claude-code"   # "claude-code" (default) | "claude" | "openai"
  claude_code_model: "sonnet"
  max_tokens: 1500
  max_retries: 2
```

If the configured backend's prerequisite is missing at runtime (`claude` CLI not on `PATH`, `ANTHROPIC_API_KEY` unset, or `OPENAI_API_KEY` unset), the pipeline emits a warning and silently falls back to data-driven commentary (no AI, no cost). Check the run log for a `Falling back to data-driven commentary` message to verify the backend you configured actually engaged.

------------------------------------------------------------------------

## Status

**Current version:** v0.2.1

### Implemented

-   **Proteomics**: Preprocessing, Multi-imputation DE (Limma), Clustering (Hierarchical, k-means/PAM, Binary patterns), Pathway enrichment, PPI networks, Advanced statistics
-   **RNA-seq**: Full pipeline (DESeq2), Pathway enrichment (fGSEA/ORA)
-   **Metabolomics**: Missingness classification (MNAR/MAR), Imputation (KNN + min/2), Normalization (TSS/Median/PQN with comparison), DE (limma/t-test/Wilcoxon), Feature selection (Random Forest, PLS-DA), Pathway enrichment (QEA, ssGSEA, ORA, GSEA), LOESS drift correction, QC suite, Report generation
-   **Multi-omics**: Integration (DIABLO, MOFA, SNF), Concordance analysis, RNA-protein correlation, Cross-omics enrichment (multiGSEA, multi-ORA), Loadings-based enrichment, Foundational analysis (correlations, WGCNA), Mechanistic inference (COSMOS, TF activity, mediation), Consensus across methods, Stability analysis (bootstrap, k-fold, cluster stability), Integrated reporting, AI commentary
-   **QC**: PCA (2D/3D, multi-resolution), UMAP, Sample distance/correlation, Density plots, Outlier detection
-   **Plots**: Volcano, MA, Heatmaps, Profile plots (3-color Up/Down/NS scheme)
-   **Reporting**: Interactive HTML reports, Executive summaries, Pipeline summaries, AI figure commentary
-   **Infrastructure**: Docker support, CLI wizard (`run.R`), Environment-variable config, Organism auto-detection, Multi-organism annotation
-   **Architecture**: Strict dependency loading, `{targets}` orchestration, Unified config validation

