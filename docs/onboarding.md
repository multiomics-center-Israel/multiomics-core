# Onboarding — multiomics-core

This guide takes you from a fresh clone to a finished pipeline run on the
synthetic example data that ships with the repo, and then to a run of your own.
It assumes you will mostly **run analyses**; if you will also change code, carry
on to [PROJECT_STRUCTURE.md](../PROJECT_STRUCTURE.md) and
[CONTRIBUTING.md](../CONTRIBUTING.md) afterwards.

Read the [README](../README.md) first: its *Requirements* and *Package
repositories* sections cover R, Rtools (Windows) and Bioconductor, and this
guide does not repeat them.

------------------------------------------------------------------------

## 0. One rule before anything else: no real data in git

This repository is **public**. Input data, per-project configs and results
never go into it:

-   `data/*` is git-ignored, except the synthetic `data/example_*` sets and a
    few reference files.
-   `config/*.yaml` is git-ignored; only `config/templates/` is tracked.
-   `outputs/` and `_targets*` stores are git-ignored.

Keep project data and results **outside the repo**, in the folder you set as
`project.dir` (section 3). Never paste sample IDs, values or result counts from
a real project into commit messages, issues or pull requests.

------------------------------------------------------------------------

## 1. Get the code and the R environment

``` bash
git clone <REPO_URL>
```

Open the project in RStudio with *File → Open Project…* and pick
`omics-core.Rproj` in the clone. Always work inside the project, so that paths
resolve relative to the repo root.

Then, in the RStudio **Console**:

``` r
install.packages("renv")   # once per machine
renv::restore()            # installs the exact versions in renv.lock
renv::status()             # should report no problems
```

`renv::restore()` takes a while the first time. If it fails on Windows, check
Rtools (README → *Windows users*).

------------------------------------------------------------------------

## 2. First run: the bundled example data

The quickest way to see the whole pipeline work is the browser wizard, which can
load a synthetic example dataset for you.

From a terminal in the repo root:

``` bash
Rscript run.R --wizard
```

Your browser opens on the wizard (at `http://localhost:8080`, or the next free
port up to 8099). Then:

1.  Click **Load Example (Proteomics)** or **Load Example (Metabolomics)**.
2.  Scroll down and click **Run Pipeline**.
3.  When the run finishes, the HTML report opens.

The wizard saves the config it built under `config/` and writes results under
`outputs/` in the repo; both are git-ignored.

> **Don't use "Load Example (Lipidomics)" yet.** Lipidomics is in development:
> its code exists, but it is not currently wired into the main `_targets.R`
> plan, so that run will not build the lipidomics analysis.

Other entry points of `run.R`: `Rscript run.R --help`.

------------------------------------------------------------------------

## 3. Your own project: the config file

Everything about an analysis — where the data is, which modes run, every
parameter — lives in **one YAML config file**. Starting a new project should
never require changing R code.

### 3.1 Start from a template

| Mode | Template |
|---|---|
| RNA-seq | `config/templates/rna_config.yaml` |
| Proteomics | `config/templates/proteins_config.yaml` |
| Metabolomics | `config/templates/metabolomics_template.yaml` |
| Multi-omics integration | `config/templates/multiomics_config.yaml` |
| Integration of finished DE tables | `config/templates/de_integration_config.yaml` |
| Lipidomics (in development, see above) | `config/templates/lipidomics_config.yaml` |

``` bash
cp config/templates/proteins_config.yaml config/<PROJECT>_<ROUND>.yaml
```

The RNA-seq, proteomics and metabolomics templates mark the fields you must
fill in with `REQUIRED`. The ones every project needs:

-   `project.dir` — absolute path to the project folder, **outside the repo**
-   `project.name` and `project.analysis_round` (e.g. `A01`)
-   `paths.raw` and `paths.out` — folders *relative to* `project.dir` for the
    input data and the results
-   `params.seed` — the seed for every random step, so a rerun gives the same
    result
-   `modes.<mode>.files.*` — the input files, relative to `paths.raw`

### 3.2 Point the pipeline at your config

The pipeline reads the config path from the `MULTIOMICS_CONFIG` environment
variable. If it is not set, it uses `config.yaml` in the repo root.

``` bash
cp .Renviron.example .Renviron
```

Edit `.Renviron` so `MULTIOMICS_CONFIG` holds the **absolute** path of your
config, then restart R (`.Renviron` is read at start-up). `.Renviron` is
git-ignored, so everyone can point at their own config.

> `.Renviron.example` ships with a placeholder path. If you copy it, set the
> real path or delete that line, or `tar_make()` will fail looking for a file
> that does not exist.

------------------------------------------------------------------------

## 4. Running the pipeline

From R, in the repo root:

``` r
library(targets)
tar_make()
```

`tar_make()` builds every mode that has a block under `modes:` in your config,
and on later runs rebuilds only what your changes affect. To build one mode,
filter by target prefix:

``` r
tar_make(names = starts_with("prot_"))  # proteomics
tar_make(names = starts_with("rna_"))   # RNA-seq
tar_make(names = starts_with("met"))    # metabolomics (met_* and metab_*)
```

Useful while it runs, or after:

``` r
tar_visnetwork()        # the dependency graph
tar_progress()          # what has run, what is running
tar_read(prot_de_res)   # load one target's value into your session
```

The alternative is `Rscript run.R --config config/<PROJECT>_<ROUND>.yaml`.
It runs the same pipeline, checks that the RNA-seq and proteomics input files
exist first, and keeps a **separate `{targets}` store per project**, so two
projects never share a cache.

> **Careful:** `tar_destroy()` and `Rscript run.R --fresh` delete the cache, and
> everything is recomputed on the next run, which can take a long time. Use them
> only when you mean it.

------------------------------------------------------------------------

## 5. Where the results are

Results are written to:

```
<project.dir>/<paths.out>/Results_<project.name>_<analysis_round>/<mode>/
```

Each run also writes `execution_info/` next to the mode folders, holding the
config as used, its path, a timestamp, `sessionInfo()`, a copy of `_targets.R`
and, when git is available, the commit hash. Keep it with the results: it is
what makes the run reproducible.

To share results, zip the relevant `Results_*` folder.

------------------------------------------------------------------------

## 6. Running steps by hand (debugging)

To step through part of the pipeline in your session, load the functions in the
same order `_targets.R` does:

``` r
for (layer in c("core", "services", "domain", "modules")) {
  files <- sort(list.files(file.path("R", layer), pattern = "\\.R$",
                           full.names = TRUE, recursive = TRUE))
  invisible(lapply(files, source))
}

config <- validate_config(load_config(Sys.getenv("MULTIOMICS_CONFIG", "config.yaml")))

prot_inputs <- load_proteomics_inputs(config)
prot_pre    <- preprocess_proteomics(prot_inputs, config)
```

For anything that is a chain of targets rather than a single function, run
`tar_make()` once and read the intermediate result with `tar_read()`, e.g.
`tar_read(metab_pre)` for metabolomics preprocessing.

------------------------------------------------------------------------

## 7. Reproducing a previous run

1.  Check out the commit recorded in that run's `execution_info/`.
2.  Run `renv::restore()` with that commit's `renv.lock`.
3.  Use the config from `execution_info/config_used.yaml`, including its
    `params.seed`.

Some steps can still differ slightly between reruns. For example, mummichog
uses unseeded permutations (see [docs/mummichog.md](mummichog.md)).

------------------------------------------------------------------------

## Where next

-   **[PROJECT_STRUCTURE.md](../PROJECT_STRUCTURE.md)** — how the code is
    organised, if you want to understand what you just ran.
-   **[CONTRIBUTING.md](../CONTRIBUTING.md)** — before you change any code.
-   **[docs index](README.md)** — contracts, the Shiny payload, mummichog,
    Docker.
