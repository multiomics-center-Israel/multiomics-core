# Project structure — multiomics-core

A practical map of the code: what each layer holds, how a config turns into
results, where each omics lives, and what not to edit by hand. It is the third
stop on the reading path (README → [onboarding](docs/onboarding.md) → this file →
[CONTRIBUTING](CONTRIBUTING.md)).

-   How to install and run: [README](README.md) and [onboarding](docs/onboarding.md).
-   How to branch, commit and open PRs: [CONTRIBUTING](CONTRIBUTING.md).
-   Rules for AI coding agents: `CLAUDE.md`. Agents read it before this file.

If this map and the code disagree, the code is right. Fix the map.

---

## 1. How a run flows

```
config YAML ──► _targets.R ──► pipe_<mode>() ──► tar_target(...) ──► mod_*() ──► domain functions ──► core helpers
 (modes: …)     reads config     R/pipeline/        one per step        R/modules/    R/domain/<omics>/     R/core/
                path, adds a                                                                │
                pipeline per                                                                ▼
                modes: block                                writer targets (format = "file") ──► Results_<name>_<round>/<mode>/
```

1.  **Config.** One YAML file decides what runs and with which parameters. Its path comes
    from `MULTIOMICS_CONFIG`, falling back to `config.yaml` in the repo root.
2.  **`_targets.R`** defines four foundation targets (`config_file`, `config`, `run_dir`,
    `execution_info_files`). It then reads the raw YAML and, for every `modes:` block it
    finds, appends that mode's pipeline factory. It holds target definitions only, no
    helpers.
3.  **Pipeline factories** (`pipe_<mode>()` in `R/pipeline/<mode>/`) return the list of
    `tar_target()`s for one mode: the DAG.
4.  **Modules** (`mod_*()` in `R/modules/<mode>/`) are what a target calls. Each one
    combines domain functions into a single cacheable step.
5.  **Domain functions** (`R/domain/<mode>/`) do the omics-specific work. They use the
    generic helpers in `R/core/`.
6.  **Writers** (`write_*()`, called from targets with `format = "file"`) put tables,
    plots and reports on disk, and return the file paths so `{targets}` can track them.

A real chain, in proteomics:

```
prot_inputs   load_proteomics_inputs(config)                          domain
prot_pre_raw  preprocess_proteomics(prot_inputs, config)              domain
prot_pre      mod_proteomics_batch_correction(prot_pre_raw, …)        module
prot_de_res   mod_proteomics_de(prot_pre, prot_inputs, config)        module → make_imputations_proteomics(), run_limma_multimp()
prot_exports  mod_proteomics_exports(prot_pre, prot_de_res, …)        module → writers (files)
```

Metabolomics is a longer chain: `metab_inputs` → `met_raw` → `met_filtered` →
`met_norm_<method>` → `metab_pre` → `metab_de_res` → `metab_final_results`, … .
Run `targets::tar_visnetwork()` to see any mode's full graph.

When `modes.multiomics` is present and at least two single-omics modes are
configured, the single-omics factories run with `skip_outputs = TRUE`: they
build only what the integration needs (inputs, preprocessing, DE), with no QC,
reports or exports.

---

## 2. The five layers of `R/`

`_targets.R` sources the layers in this order. Within a layer, files are
sourced in `sort()` order, which is why most file names carry a number.

| Layer | Holds | May use |
|---|---|---|
| `R/core/` | Generic utilities with no omics knowledge | nothing above it |
| `R/services/` | External integrations (AI figure commentary) | `core` |
| `R/domain/<omics>/` | Omics-specific logic | `core`, `services` |
| `R/modules/<omics>/` | `mod_*()` wrappers that a target calls | the layers below |
| `R/pipeline/<omics>/` | `pipe_<omics>()` factories that build the DAG (sourced last) | everything |

**Dependencies point downward only.** A helper shared by several omics belongs
in the lowest layer that fits its meaning; don't move omics-specific logic into
`core` just to make it reachable.

### `R/core/`, grouped by role

| Files | Role |
|---|---|
| `00_paths.R`, `01_io.R` | Path resolution (`project.dir`, `paths.raw`, `paths.out`), loading tables, shared contrast-name helpers |
| `02_validation.R`, `03_alignment.R`, `04_config.R` | Input validation, sample/metadata alignment, `load_config()` / `validate_config()`, execution metadata |
| `05_export_excel.R` | Final-results tables and Excel export, shared by all omics (`build_final_results_generic()`, `get_contrast_cols()`) |
| `06_plots.R`, `08_qc.R`, `17_cutoff_panel.R` | Generic plots (volcano, MA, heatmaps), QC (PCA, distances, outliers) |
| `07_shiny_contract.R` | The Shiny payload contract (see §7) |
| `09_clustering.R`, `09_enrichment.R`, `13_gmt_utils.R` | Clustering, generic enrichment (fGSEA/ORA), GMT files |
| `10_organism_detection.R`, `11_annotation.R` | Organism detection and gene annotation |
| `15_user_summary.R`, `16_pptx_helpers.R`, `de_summary_counts.R`, `de_shrinkage_check.R`, `results_fact_sheet.R` | Reporting helpers and result sanity checks |
| `18_targets_store.R` | The per-project `{targets}` store that `run.R` uses |

### `R/domain/<omics>/`

Numbers follow the analysis flow: `00_inputs.R` loads the inputs, then come
preprocessing, DE, outputs, enrichment and reports, with Shiny export,
PowerPoint and commentary near the end. `90_config_validate.R`, where a mode
has one, holds that mode's config checks and is called by `validate_config()`.
RMarkdown report templates live next to the code (`R/domain/<omics>/report_template*.Rmd`)
or under `R/pipeline/<omics>/templates/` (metabolomics, lipidomics).

`R/domain/metabolomics/README.md` and `R/domain/multiomics/README.md` describe those
two modes; read them before working on either.

---

## 3. Where each omics lives

| Mode | Config block | Template | Code | Target prefix | Status |
|---|---|---|---|---|---|
| RNA-seq | `modes.rna` | `rna_config.yaml` | `R/{domain,modules,pipeline}/rnaseq/` | `rna_` | runs via `tar_make()` |
| Proteomics | `modes.proteomics` | `proteins_config.yaml` | `…/proteomics/` | `prot_` | runs via `tar_make()` |
| Metabolomics | `modes.metabolomics` | `metabolomics_template.yaml` | `…/metabolomics/` | `met_`, `metab_` | runs via `tar_make()` |
| Multi-omics integration | `modes.multiomics` | `multiomics_config.yaml` | `…/multiomics/` | `multiomics_` | runs when ≥2 single-omics modes are configured too |
| DE-table integration | `modes.de_integration` | `de_integration_config.yaml` | `…/de_integration/` | `dei_` | runs on finished DE tables; no single-omics mode needed. In progress: so far it loads and standardises the layers and resolves the contrasts |
| Lipidomics | `modes.lipidomics` | `lipidomics_config.yaml` | `…/lipidomics/` | `lipid_` | **in development: `pipe_lipidomics()` exists but `_targets.R` does not call it** |

Metabolomics specifics that are easy to get wrong:
`modes.metabolomics.preprocessing` is the single place normalization is set
(`chosen_norm`, `transform`, `scaling`, `pseudocount`); there is no
`normalization:` block and no imputation stage, only missingness filtering.

### Pipeline factories

| Factory | Signature |
|---|---|
| `pipe_rnaseq()` | `(skip_outputs = FALSE)` |
| `pipe_proteomics()` | `(skip_outputs = FALSE)` |
| `pipe_metabolomics()` | `(chosen_norm = NULL, skip_outputs = FALSE, mummichog_enabled = FALSE)` |
| `pipe_multiomics()` | `()` |
| `pipe_de_integration()` | `()` |
| `pipe_lipidomics()` | `()`, not called from `_targets.R` |

---

## 4. Entry points

| Entry point | What it does |
|---|---|
| `targets::tar_make()` | Runs the plan in `_targets.R` for the config in `MULTIOMICS_CONFIG`, in the default store `_targets/` (or whatever a local `_targets.yaml` names) |
| `Rscript run.R --config <file>` | Same pipeline, with a pre-flight check of RNA-seq and proteomics input files and a separate store per project (`_targets_<run>_<hash>/`) |
| `Rscript run.R --wizard` | Browser wizard: build a config, load an example dataset, run |
| `Rscript run.R --new` / `--fresh` / `--help` | Terminal wizard / full rerun of the last config (clears that project's cache) / usage |
| `docker-run.bat`, `Dockerfile` | The wizard in Docker; see `DOCKER.md` |
| `install.bat`, `launch_wizard.*`, `configure_wizard.*`, `create_shortcut.vbs` | Windows install and launcher helpers |
| `make setup` | Builds the pinned mummichog venv (see `docs/mummichog.md`) |
| `tests/testthat/` | Unit tests. The end-to-end `test-e2e-*.R` files are opt-in (`RUN_E2E_OMICS=1`). CI: `.github/workflows/ci.yml` |

---

## 5. Configuration

-   `config/templates/*.yaml` are tracked. Copy one to start a project.
-   `config/*.yaml` (per-project configs) are git-ignored. Real configs usually live
    in the project folder, outside the repo.
-   Top-level sections: `project` (`dir`, `name`, `analysis_round`, …), `paths`
    (`raw`, `out`, relative to `project.dir`), `params` (`seed`), `modes`
    (one block per mode), and the optional `commentary`.

Each mode's `90_config_validate.R` (where present), together with
`validate_config()`, is the contract a config must satisfy. The templates show
every key with its default.

---

## 6. Data, outputs and caches

| Path | What it is | In git? |
|---|---|---|
| `data/example_{proteomics,metabolomics,lipidomics}/` | Small synthetic datasets, each with a `generate_example.R` | yes |
| `data/HMDB2kegg_cpd.*.txt` | HMDB → KEGG compound mapping used by the pipeline | yes |
| other files directly under `data/` | Older tracked files; whether they belong in the repo is under review | yes, for now |
| `<project.dir>/<paths.out>/Results_<name>_<round>/<mode>/` | Results, plus `execution_info/` next to the mode folders | no |
| `outputs/` | Where the wizard writes results when `project.dir` is the repo | no |
| `_targets/`, `_targets_*/` | `{targets}` caches | no, never edit |

---

## 7. The Shiny app

The Shiny app is **not in this repository**; it reads payloads this pipeline
writes. The contract is defined in `R/core/07_shiny_contract.R`
(`init_shiny_payload()`, `assert_shiny_payload_contract()`,
`build_feature_annot()`) and filled by the per-omics `*_shiny_export.R` files.
The keys and per-omics differences are documented in
[docs/SHINY_PAYLOAD_HANDOFF.md](docs/SHINY_PAYLOAD_HANDOFF.md). Changing a payload
key breaks the app, so agree such changes with the app owner first.

---

## 8. Design rules

These keep the pipeline reproducible and the layers clean.

1.  **Config holds analysis decisions; code is infrastructure.** Anything that can vary
    between runs is a config key with a default, not a constant in R. Starting a new
    project should never need an R change.
2.  **Computation and disk I/O are separate.** Domain and module functions take objects
    and return objects. Writing files happens in `write_*()` functions called from
    targets with `format = "file"`, which return the paths they wrote. Plot functions
    return plots; a writer saves them.
3.  **Function or target?** Write a *function* for reusable, testable logic that you also
    want to call interactively. Add a *target* when a result is expensive, feeds several
    later steps, or is a file. In short, the function is *what* is computed and the
    target is *when*, with which dependencies.
4.  **Validate early, fail loudly.** Check inputs right after loading, after sample
    mapping and after building matrices (dimensions, numeric type, alignment with the
    metadata, unique IDs). An error should say what was expected, what was found and
    how to fix it.
5.  **Stable return contracts.** Downstream code relies on the fields a stage returns
    (e.g. `expr_work`, `meta`, `row_data` from preprocessing; `summary_df` from DE).
    Add fields rather than rename them, and give new config keys defaults. The
    documented contracts are in `docs/CONTRACT_*.md`.
6.  **One seed.** Random steps derive their seed from `params.seed`; no ad-hoc
    `set.seed()`.
7.  **Debug with the pipeline's own functions.** Don't keep debug-only copies of logic.
8.  **Think across modes.** Every mode has matching `domain/`, `modules/` and
    `pipeline/` folders. When a feature lands in one mode, check whether the others
    need it, and whether it belongs in `core`.

### Naming, as the code does it

-   **Targets:** a mode prefix plus what the target holds: `prot_inputs`, `prot_pre`,
    `prot_de_res`, `metab_final_results`, `rna_pathway_res`. (`CLAUDE.md` currently
    describes target names as `<verb>_<noun>`; the code uses this mode-prefix
    pattern.)
-   **Functions:** `snake_case`. Common shapes: `load_<mode>_inputs()`,
    `preprocess_<mode>()`, `mod_<mode>_<step>()`, `write_<mode>_<artifact>()`,
    `pipe_<mode>()`, `validate_*()`.
-   **Objects:** `meta` for sample metadata; `cfg` for `config$modes$<mode>`; `config`
    for the whole config.
-   **Docs:** every function in `R/` has a roxygen2 docstring.

---

## 9. Don't edit these directly

| What | Why / instead |
|---|---|
| `_targets/`, `_targets_*/` | Pipeline caches. To rebuild, change an input or ask before clearing a cache |
| `renv.lock` | Dependency changes are a deliberate team decision |
| `config/templates/*.yaml` | Shared starting points. Copy one for a project; change a template only for a project-wide default |
| `R/core/07_shiny_contract.R`, `*_shiny_export.R` | Contract with the separate Shiny app (§7) |
| `.github/` | CI workflow and the PR template |
| `data/example_*/` CSVs | Regenerate with that folder's `generate_example.R` |
| `.Renviron`, any key files | Local, may hold secrets; never commit |

---

## 10. Other top-level folders

-   `docs/`: documentation index, contracts, `archive/` for historical plans.
-   `scripts/`, `utils/`, `tools/`: standalone helper scripts (KEGG/GMT generation,
    deck builders, the MOFA bridge, a pandoc installer). Some are called by the
    pipeline (for example `scripts/run_mofa.py`, `utils/build_feature_ko_map.R`,
    `tools/install_pandoc.R`). Note that `.gitignore` ignores new `scripts/*.R`
    files, although some `.R` files there are tracked.
