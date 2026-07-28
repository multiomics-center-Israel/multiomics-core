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

Plus, directly under `mummichog_pinned/`, the report's presentation exports per contrast: `mummichog_pathway_bubble_<contrast>.{png,pdf}` (the bubble plot) and `mummichog_pathway_table_<contrast>.{tsv,csv}` (the sorted pathway table), and `contrasts.tsv`, which maps each sanitised subdirectory name back to its original DE contrast label (so the report can show real contrast names).

To map pathways back to your feature ids, `join_features_to_results()` uses the feature id mummichog echoes into its own tables (via the 5th input column) — not the fragile post-de-duplication row numbers.

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

