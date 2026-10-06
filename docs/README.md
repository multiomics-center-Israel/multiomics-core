# Documentation index

## Reading path for new team members

Read these in order. Stop after step 2 if you only run analyses; continue to
step 4 before you change code.

1. [README](../README.md) — what the pipeline does, and how to install it.
2. [Onboarding](onboarding.md) — from a fresh clone to your first pipeline run.
3. [Project structure](../PROJECT_STRUCTURE.md) — the architecture: the five
   layers under `R/`, how `_targets.R` assembles the plan, and the config layout.
4. [Contributing](../CONTRIBUTING.md) — branches, commits, pull requests and
   review. You only need this if you are changing the code.

`CLAUDE.md` in the repository root holds the rules for AI coding agents working
in this repo. It is worth skimming, but it is not part of the reading path.

## Reference

Open these when you need them; you don't have to read them up front.

| Document | What it covers |
|---|---|
| [DOCKER.md](../DOCKER.md) | Running the pipeline and the config wizard in Docker |
| [mummichog.md](mummichog.md) | Optional mummichog pathway analysis for metabolomics: setup, models, outputs |
| [R/domain/multiomics/README.md](../R/domain/multiomics/README.md) | Multi-omics integration: config, data flow, troubleshooting (some output paths and file names are out of date) |
| [R/domain/metabolomics/README.md](../R/domain/metabolomics/README.md) | Metabolomics mode: preprocessing flow, `chosen_norm`, main targets, where the code lives |
| [CONTRACT_rnaseq_final_results.md](CONTRACT_rnaseq_final_results.md) | Columns of the RNA-seq final results and the run artefacts |
| [CONTRACT_metabolomics_de_and_final_results.md](CONTRACT_metabolomics_de_and_final_results.md) | Metabolomics DE column naming and the final results schema |
| [SHINY_PAYLOAD_HANDOFF.md](SHINY_PAYLOAD_HANDOFF.md) | What the pipeline emits for the Shiny app |
| [DE_to_DA_rename_map.md](DE_to_DA_rename_map.md) | Impact survey for a possible DE → DA rename in metabolomics (nothing renamed yet) |

## Archive

Plans, handovers, audits and prompt logs that no longer describe the current
code live in [archive/](archive/README.md).
