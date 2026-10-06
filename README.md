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
-   Omics-specific preprocessing (filtering and normalization; imputation in proteomics only)
-   Proteomics differential expression via **limma** with multiple imputations and stability filtering
-   RNA-seq differential expression via **DESeq2**
-   Metabolomics preprocessing (missingness filtering — no imputation — then sample normalization: none/TSS/median/PQN/EigenMS/bio-factor; optional LOESS drift correction) and differential abundance
-   Multi-omics integration, and integration of finished DE tables (`de_integration`)
-   Lipidomics — **in development**: the code exists under `R/*/lipidomics/`, but it is not currently wired into the main `_targets.R` plan, so `tar_make()` does not run it
-   Pathway enrichment analysis (fGSEA, ORA, QEA, ssGSEA)
-   Unified QC utilities (PCA, heatmaps, sample distance)
-   A central YAML configuration file controlling all parameters
-   A `{targets}` pipeline for reproducible, dependency-aware execution

------------------------------------------------------------------------

## Getting started (new users)

Read these in order:

1.  **This README** — what the pipeline does and how to install it.
2.  📘 **[Onboarding](docs/onboarding.md)** — from a fresh clone to a first run on the bundled example data.
3.  📘 **[Project structure](PROJECT_STRUCTURE.md)** — the architecture: layers, pipelines, config layout.
4.  📘 **[Contributing](CONTRIBUTING.md)** — only if you are going to change code.

Everything else (contracts, the Shiny payload, mummichog, Docker) is listed in the **[docs index](docs/README.md)**.

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
├── templates/   # Config templates (rna, proteins, metabolomics, lipidomics, multiomics, de_integration)
data/            # Synthetic example datasets (data/example_*) and reference files
docs/            # Documentation index, onboarding, contracts; archive/ holds historical plans
tests/           # testthat tests
_targets.R       # {targets} pipeline definition
run.R            # CLI entrypoint and browser wizard (Rscript run.R --help)
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
-   filtering and normalization settings (and imputation, for proteomics)
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

Targets are named by mode prefix: `prot_`, `rna_`, `met_`/`metab_`, `multiomics_`, `dei_` (DE-table integration) and `lipid_`. The `lipid_` targets exist in the code but are not currently added to the plan by `_targets.R`.

Alternatively, `Rscript run.R --config path/to/config.yaml` runs the same pipeline with a pre-flight check for missing RNA-seq and proteomics input files, and a separate `{targets}` store per project (see `Rscript run.R --help`).

`{targets}` ensures that only steps affected by changes are recomputed.

### Learning more about `{targets}`

This project uses `{targets}` for reproducible, dependency-aware pipeline orchestration.

For a detailed introduction, tutorials, and best practices, see the official **targets** book: <https://books.ropensci.org/targets/>

------------------------------------------------------------------------

## Mummichog pathway analysis (metabolomics, optional)

The metabolomics mode can run [mummichog](http://mummichog.org) 2.7.0 for m/z-based pathway enrichment, in a pinned Python venv that is kept out of git. It is **off by default**; turn it on with `enabled: true` under `modes.metabolomics.enrichment.mummichog`. When it is off, the venv is never needed.

One-time setup is `make setup`, which builds the venv and prints the `MUMMICHOG_PYTHON=...` line to add to your `.Renviron`.

Setup, model selection (`model_ref` / `model_json` / built-in human model), output layout and the stochasticity caveat are in **[docs/mummichog.md](docs/mummichog.md)**.

------------------------------------------------------------------------

## Running preprocessing interactively (example)

For exploratory work or debugging:

``` r
# Load functions in the same order as _targets.R (core -> services -> domain -> modules).
# The .R pattern matters: R/domain also holds README.md and .Rmd report templates.
for (layer in c("core", "services", "domain", "modules")) {
  files <- sort(list.files(file.path("R", layer), pattern = "\\.R$",
                           full.names = TRUE, recursive = TRUE))
  invisible(lapply(files, source))
}

# Load config (the same file the pipeline uses)
config <- validate_config(load_config(Sys.getenv("MULTIOMICS_CONFIG", "config.yaml")))

# --- Proteomics ---
inputs <- load_proteomics_inputs(config)
res    <- preprocess_proteomics(inputs, config)

# --- RNA-seq ---
inputs <- load_rna_inputs(config)
res    <- preprocess_rna(inputs, config)

# --- Metabolomics ---
# Metabolomics preprocessing is a chain of targets rather than one function;
# after a tar_make(), read its result from the store instead:
res <- targets::tar_read(metab_pre)

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
-   Each run writes `execution_info/` into its results folder: the config used, its path, a timestamp, `sessionInfo()`, a copy of `_targets.R`, and the git commit when git is available

------------------------------------------------------------------------

## Outputs

All analysis outputs are written under `<project.dir>/<paths.out>/Results_<project.name>_<analysis_round>/<mode>/`, where path components come from your config YAML (defaults: `paths.out: "outputs"`).

-   The project directory lives outside the repository (set via `project.dir` in your config)
-   Results should be shared by zipping the relevant output folder
-   Each run is isolated by its configuration parameters

------------------------------------------------------------------------

## Developer notes

If you want to extend, modify, or maintain **multiomics-core**, read [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) and then [CONTRIBUTING.md](CONTRIBUTING.md).

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
-   **Metabolomics**: Missingness filtering (no imputation), Normalization (none/TSS/Median/PQN/EigenMS/bio-factor, with comparison), DE (limma/t-test/Wilcoxon), Feature selection (Random Forest, PLS-DA), Pathway enrichment (QEA, ssGSEA, ORA, GSEA), LOESS drift correction, QC suite, Report generation
-   **Multi-omics**: Integration (DIABLO, MOFA, SNF), Concordance analysis, RNA-protein correlation, Cross-omics enrichment (multiGSEA, multi-ORA), Loadings-based enrichment, Foundational analysis (correlations, WGCNA), Mechanistic inference (COSMOS, TF activity, mediation), Consensus across methods, Stability analysis (bootstrap, k-fold, cluster stability), Integrated reporting, AI commentary
-   **QC**: PCA (2D/3D, multi-resolution), UMAP, Sample distance/correlation, Density plots, Outlier detection
-   **Plots**: Volcano, MA, Heatmaps, Profile plots (3-color Up/Down/NS scheme)
-   **Reporting**: Interactive HTML reports, Executive summaries, Pipeline summaries, AI figure commentary
-   **Infrastructure**: Docker support, CLI wizard (`run.R`), Environment-variable config, Organism auto-detection, Multi-organism annotation
-   **Architecture**: Strict dependency loading, `{targets}` orchestration, Unified config validation

