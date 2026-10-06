# Archived documents

These files are kept for history. **None of them describes how the pipeline
works today** — plans, handovers and audits record what someone intended or
found at a point in time, and the code has moved on since. When a document here
disagrees with the code, the code is right.

For current documentation, start at the [docs index](../README.md).

Paths quoted inside these files (e.g. `docs/SPEC_enrichment_migration.md`)
refer to where the files lived when they were written, before they were moved
here.

| File | Written | What it was | Why archived |
|---|---|---|---|
| [CLAUDE_CODE_FIX_PROMPTS.md](CLAUDE_CODE_FIX_PROMPTS.md) | 2026-03 | Prompts for Claude Code bug-fix sessions | One-off working notes; nothing references them |
| [CLAUDE_CODE_MIGRATION_PROMPTS.md](CLAUDE_CODE_MIGRATION_PROMPTS.md) | 2026-02 | Prompts for the agents that ran the `multiomics_pipeline` migration | The migration is done |
| [MIGRATION_PLAN_multiomics_pipeline.md](MIGRATION_PLAN_multiomics_pipeline.md) | 2026-02 | Plan for moving `multiomics_pipeline` into this repo | The migration is done; the plan quotes paths on a developer machine |
| [PLAN_enrichment_migration_v2.md](PLAN_enrichment_migration_v2.md) | — | Phase plan and status log for the enrichment migration branch | Describes a branch's status at the time, not `main` |
| [HANDOVER_enrichment_migration.md](HANDOVER_enrichment_migration.md) | 2026-06 | Handover for the enrichment migration branch | Its own header records the branch as unmerged |
| [SPEC_enrichment_migration.md](SPEC_enrichment_migration.md) | — | Design spec for the enrichment migration | Marked "frozen design note"; companion to the two files above |
| [audit-enrichment-layer.md](audit-enrichment-layer.md) | 2026-04 | Audit of the enrichment code at that date | A point-in-time snapshot |
| [ROADMAP_metabolomics_pipeline.md](ROADMAP_metabolomics_pipeline.md) | 2026-03 | Metabolomics roadmap | Last updated in March; describes a pipeline that has since changed (see the metabolomics section of `CLAUDE.md`) |
| [plan-tss-loess-normalization.md](plan-tss-loess-normalization.md) | — | Implementation plan for TSS / CLR / LOESS / QC-RLSC normalization | A plan, not a description; normalization now lives in `R/domain/metabolomics/01_normalization.R` and `10_drift_correction.R` |
| [lipidomics_enhancements_plan.md](lipidomics_enhancements_plan.md) | 2026-05 | Plan for lipidomics enhancements on the `Lipidomics_Oz` branch | A branch plan, not a description of `main` |
| [developer_guide.md](developer_guide.md) | — | Earlier developer guide: principles, naming, adding a step | Its current content now lives in `PROJECT_STRUCTURE.md` (§8, design rules and naming) and `CONTRIBUTING.md`; the rest named directories that no longer exist (`R/plots/`, `R/qc/`) or target patterns the code does not use |
| [MULTIOMICS_QUICKSTART.md](MULTIOMICS_QUICKSTART.md) | — | Multi-omics quick start (formerly in the repo root) | Duplicated `R/domain/multiomics/README.md` (config, flow, `tar_read`, outputs, troubleshooting) plus general setup now in onboarding; it also named a config that is not in the repo, and its flow diagram used abbreviated target names that do not match the code |
