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
-   RNA-seq differential expression via **DESeq2** with batch correction (ComBat-Seq/SVA/RUV) and cell-type deconvolution
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

The metabolomics mode can run [mummichog](http://mummichog.org) for m/z-based pathway/network enrichment. There are **two independent engines** in the repo:

-   **`06b_mummichog.R`** — runs mummichog via `{reticulate}` and builds custom KEGG organism models (honey bee, *E. coli*, human…) with `{KEGGREST}`. This is the current default path.
-   **`06c_mummichog_pinned.R`** — runs the **version-pinned** `mummichog==2.7.0` as an isolated subprocess via `{processx}`, depending only on light R packages (`readr`, `processx`, `jsonlite`) — no Bioconductor. This is the new, reproducible engine documented here.

### Setup

The pinned engine calls Python in a dedicated venv, kept out of git (`envs/` is `.gitignore`d):

``` bash
make mummichog-venv                 # creates envs/mummichog, writes requirements-mummichog.lock
# or, to reproduce the exact committed tree:
make mummichog-lock                 # installs from requirements-mummichog.lock (USE_LOCK=1)
```

Then point the R wrapper at that interpreter (the default is `envs/mummichog/bin/python`):

``` bash
export MUMMICHOG_PYTHON="$(pwd)/envs/mummichog/bin/python"
```

Both `requirements-mummichog.txt` (the top-level pin) and `requirements-mummichog.lock` (the fully-resolved tree) are committed.

### How to run

The pinned stage is **not wired into the metabolomics DAG by default** (it coexists with `06b`). To enable it, add these targets to `R/pipeline/metabolomics/00_pipe_metabolomics.R` (they slot in beside `metab_enrichment_res`):

``` r
tar_target(
  metab_mummichog_pinned_files,
  mod_mummichog_pinned(
    pre     = metab_pre,
    de_res  = metab_de_res,
    config  = config,
    out_dir = metab_out_dir
  ),
  format = "file"
),
tar_target(mummichog_pinned_pathways, read_mummichog_pathways(metab_mummichog_pinned_files)),
tar_target(mummichog_pinned_modules,  read_mummichog_modules(metab_mummichog_pinned_files))
```

It reads the same knobs `06b` uses, under `modes.metabolomics.enrichment.mummichog` in your config:

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
        # model_json: /path/to/ame_metabolic_model.json  # optional: reuse a 06b organism model
```

Then run as usual:

``` r
library(targets)
tar_make(names = starts_with("met"))
```

### Where outputs land

Under `<metab_out_dir>/mummichog_pinned/`:

-   `input.tsv` and `input.tsv.idmap.tsv` — the exact table sent to mummichog (m/z, retention time, p-value, statistic, **feature\_id as the 5th column**) plus a provenance id-map.
-   `v2/<timestamp>.<project>/` — the mummichog result tree: `result.html`, `tables/` (`mcg_pathwayanalysis_*.tsv`/`.xlsx`, `mcg_modularanalysis_*.tsv`/`.xlsx`, `ListOfEmpiricalCompounds.tsv`, `userInputData.txt`, `userInput_to_EmpiricalCompounds.tsv`), `figures/` and `js/`. Result tables are **`.tsv`/`.xlsx`, never `.csv`**.
-   `v2/mummichog_manifest.tsv` and `v2/runner.log`.

To map pathways back to your feature ids, `join_features_to_results()` uses the feature id mummichog echoes into its own tables (via the 5th input column) — not the fragile post-de-duplication row numbers.

### Stochasticity caveat

mummichog v2 estimates null distributions by **random permutation with no seed control**, so p-values and rankings vary slightly between runs on identical input. `{targets}` only re-runs the stage when its inputs change, so this doesn't cause spurious rebuilds — but do **not** expect bit-identical reruns, and don't assert exact equality in tests.

### Cross-checking against `06b`

For a like-for-like comparison on the same organism, build the KEGG model with `06b` and pass it to the pinned engine via `model_json` (the `-n` argument), so both engines use the same metabolic network on the same input.

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
-   **RNA-seq**: Full pipeline (DESeq2), Batch correction (ComBat-Seq/SVA/RUV), Cell-type deconvolution (xCell2), Pathway enrichment (fGSEA/ORA)
-   **Metabolomics**: Missingness classification (MNAR/MAR), Imputation (KNN + min/2), Normalization (TSS/Median/PQN with comparison), DE (limma/t-test/Wilcoxon), Feature selection (Random Forest, PLS-DA), Pathway enrichment (QEA, ssGSEA, ORA, GSEA), LOESS drift correction, QC suite, Report generation
-   **Multi-omics**: Integration (DIABLO, MOFA, SNF), Concordance analysis, RNA-protein correlation, Cross-omics enrichment (multiGSEA, multi-ORA), Loadings-based enrichment, Foundational analysis (correlations, WGCNA), Mechanistic inference (COSMOS, TF activity, mediation), Consensus across methods, Stability analysis (bootstrap, k-fold, cluster stability), Integrated reporting, AI commentary
-   **QC**: PCA (2D/3D, multi-resolution), UMAP, Sample distance/correlation, Density plots, Outlier detection
-   **Plots**: Volcano, MA, Heatmaps, Profile plots (3-color Up/Down/NS scheme)
-   **Reporting**: Interactive HTML reports, Executive summaries, Pipeline summaries, AI figure commentary
-   **Infrastructure**: Docker support, CLI wizard (`run.R`), Environment-variable config, Organism auto-detection, Multi-organism annotation
-   **Architecture**: Strict dependency loading, `{targets}` orchestration, Unified config validation

