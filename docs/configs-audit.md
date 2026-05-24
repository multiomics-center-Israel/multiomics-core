# Run-config audit — `config/templates/`

**Scope:** every YAML in `config/templates/` cross-referenced against
every R / R-Markdown file under `R/`, plus `_targets.R`, `run.R`,
`configure_wizard.R`, `validate_multiomics.R`.

**This audit is keys-only.** No values were changed, no YAML or R file
was edited. All findings cite `file:line`.

**What I parsed**

| YAML                                          | Drives           | Run type                                   |
|-----------------------------------------------|------------------|--------------------------------------------|
| `config/templates/proteins_config.yaml`       | proteomics       | Single-omic template (DIA-NN / MaxQuant)   |
| `config/templates/rna_config.yaml`            | RNA-seq          | Single-omic template (counts → DESeq2/limma) |
| `config/templates/metabolomics_template.yaml` | metabolomics     | Single-omic template (CD raw or processed wide) |
| `config/templates/multiomics_config.yaml`     | rna + proteomics + multiomics integration | Multi-omics integration template (≥2 omics + DIABLO/SNF) |

R's `validate_config()` (`R/core/04_config.R:13–29`) dispatches to four
per-omic validators in `R/domain/*/90_config_validate.R` (metabolomics
lives in `R/domain/metabolomics/00_inputs.R:127`). Only
`validate_multiomics_config()` mutates the config to apply defaults.

A short note on `prior.count`: YAML 1.1 keys are plain strings, so
`prior.count: 1` is read as a single key whose name is the literal
seven-character string `prior.count`. In R, `cfg$normalization$prior.count`
is equivalent to `cfg$normalization[["prior.count"]]` (R's `$` does not
treat `.` specially — it just takes the rest of the name as one symbol).
So the code does work today, but anyone reading the YAML naturally
expects `prior` to be a sub-object and `count` to be a leaf. This is
called out in section (i) for renaming.

---

## (a) Inventory — every YAML key and its consumers

Status legend:
- **U** = used (at least one R/Rmd read).
- **DU** = defined but unused (no R/Rmd read found).
- **UMC** = used-but-missing-in-some-configs (R reads it with a `%||%`
  default, but the key is not declared in every YAML that exercises that
  code path).
- **C** = consumed via container alias (e.g. `report_cfg <- prot_cfg$report`).

### Global (top-level)

| key_path                       | configs_defining_it                 | R_files_using_it (representative lines) | status |
|--------------------------------|--------------------------------------|-----------------------------------------|--------|
| `project.dir`                  | prot, rna, met, multi (`project.dir`) | `R/core/00_paths.R`, `R/modules/metabolomics/06_mod_report.R:289`, `R/domain/multiomics/04b_integration_mofa.R:1012` | U |
| `project.name`                 | all four                             | `R/core/00_paths.R`, `R/modules/metabolomics/06_mod_report.R:289` | U |
| `project.analysis_round`       | all four                             | `R/core/00_paths.R`, `R/modules/metabolomics/06_mod_report.R:43` | U |
| `project.analyst`              | prot, rna, multi                     | `R/modules/metabolomics/06_mod_report.R:290`, `R/domain/proteomics/10_pipeline_summary.R` | U |
| `project.user_notes`           | prot, rna                            | `R/core/15_user_summary.R:15`, `R/domain/proteomics/report_template_proteomics.Rmd:267` | U |
| `project.technical_report_file`| prot, rna                            | `R/core/15_user_summary.R:16`           | U |
| `paths.raw`                    | all four                             | `R/core/00_paths.R`, `R/pipeline/proteomics/00_pipe_proteomics.R:23` | U |
| `paths.out`                    | all four                             | `R/core/00_paths.R`, `R/domain/rnaseq/03_preprocess.R:189` | U |
| `params.seed`                  | all four                             | `R/domain/proteomics/03_imputation.R:227` | U |
| `design.condition_column`      | multi only                           | `R/domain/multiomics/12_mechanistic_inference.R:845` | U |
| `design.design_formula`        | multi only                           | `R/domain/multiomics/13_commentary.R:460` | U |
| `design.contrasts`             | multi only                           | (multiple multiomics integration sites)  | U |
| `design.reference_level`       | multi only                           | `R/domain/multiomics/06_concordance.R:231` | U |
| `design.batch_column`          | multi only                           | **none** — no grep hit                  | DU |
| `global.omics_present`         | multi only                           | `R/domain/multiomics/07_enrichment.R:26` | U |
| `global.organism`              | multi only                           | `R/domain/multiomics/12_mechanistic_inference.R:674` | U |
| `global.metadata`              | multi only                           | **none** — no grep hit                  | DU |
| `output.output_dir`            | multi only                           | only in docstring comments (`R/domain/proteomics/07c_ppi_networks.R:7`, `07d_advanced_stats.R:8`) — code uses `paths.out` instead | DU |
| `output.report_format`         | multi only                           | **none** — no grep hit                  | DU |

### Proteomics — `modes.proteomics.*`

| key_path | configs_defining_it | R_files_using_it | status |
|----------|----------------------|------------------|--------|
| `engine` | prot, multi | `R/domain/proteomics/10_pipeline_summary.R:207` | U |
| `scale_in` | prot | `R/domain/proteomics/01_expression.R:68,112` | U |
| `files.protein` | prot, multi | `R/pipeline/proteomics/00_pipe_proteomics.R:30` | U |
| `files.sample_map` | prot | `R/domain/proteomics/00_inputs.R` | U |
| `files.metadata` | prot, multi | `R/pipeline/proteomics/00_pipe_proteomics.R:23` | U |
| `files.contrasts` | prot, multi | `R/pipeline/proteomics/00_pipe_proteomics.R:24` | U |
| `files.peptides` | prot | `R/domain/proteomics/report_template_proteomics.Rmd:2284` | U |
| `id_columns.protein_id` | prot, multi | `R/domain/proteomics/01_expression.R:17`, `05_de_summary.R:212` | U |
| `id_columns.sample_col` | prot, multi | `R/core/02_validation.R:33`, `R/domain/proteomics/04_preprocess.R:44` (always behind `effects$samples %||% …`) | U (secondary) |
| `id_columns.protein_annot` | prot, multi | `R/domain/proteomics/01_expression.R:152,201`, `05_de_summary.R:214` | U |
| `id_columns.map_from` | prot | `R/domain/proteomics/00_inputs.R:53`, `R/core/03_alignment.R:42` | U |
| `id_columns.map_to` | prot | `R/domain/proteomics/90_config_validate.R:20`, `R/core/03_alignment.R:42` | U |
| `filtering.remove_contaminants` | prot | `R/domain/proteomics/90_config_validate.R:58` (validated; consumer follows in `02_filtering.R`) | U |
| `filtering.contaminant_prefix` | prot | `R/domain/proteomics/90_config_validate.R:62` | U |
| `filtering.min_count.default` | prot | `R/domain/proteomics/02_filtering.R:33` | U |
| `filtering.min_groups` | prot | `R/domain/proteomics/02_filtering.R:35` | U |
| `normalization.method` | prot, multi | `R/domain/proteomics/04_preprocess.R:72` | U |
| `imputation.multi_imputation` | prot, multi | `R/domain/proteomics/90_config_validate.R:31`, `03_imputation.R` | U |
| `imputation.method` | prot, multi | `R/domain/proteomics/09_executive_summary.R:50`, `03_imputation.R` | U |
| `imputation.no_repetitions` | prot | `R/domain/proteomics/90_config_validate.R:32` | U |
| `imputation.min_no_passed` | prot | `R/domain/proteomics/05_de_summary.R:8` (read via `imp_cfg$min_no_passed`) | U |
| `imputation.width` | prot, multi | `R/modules/proteomics/01_mod_qc_pre.R:114`, `R/domain/proteomics/03_imputation.R` | U |
| `imputation.downshift` | prot, multi | `R/modules/proteomics/01_mod_qc_pre.R:115`, `R/domain/proteomics/03_imputation.R` | U |
| `imputation.dep2_method` | prot | `R/domain/proteomics/03_imputation.R:64` | U |
| `imputation.dep2_random_seed` | prot | `R/domain/proteomics/03_imputation.R:65` | U |
| `imputation.qrilc_random_seed` | prot | `R/domain/proteomics/03_imputation.R:119` | U |
| `de.method` | prot, multi | `R/modules/proteomics/02_mod_de.R:55`, `R/domain/proteomics/05_de_summary.R:…` | U |
| `de.use_adj_for_pass1` | prot | `R/domain/proteomics/05_de_summary.R:16` | U |
| `de.p_cutoff` | prot, multi | `R/domain/proteomics/report_template_proteomics.Rmd:450`, `09b_commentary.R` | U |
| `de.linear_fc_cutoff` | prot, multi | `R/modules/proteomics/04_mod_clustering.R:55`, `report_template_proteomics.Rmd:450`, `09b_commentary.R:262` | U |
| `de.p_adjust_method` | prot | `R/domain/proteomics/05b_de_anova.R:19` | U |
| `de.fdrtool_correction` | prot | `R/domain/proteomics/05_de_summary.R:243` | U |
| `de.paired` | prot | `R/domain/proteomics/05c_de_ttest.R:21` | U |
| `de.pairing_col` | prot | `R/domain/proteomics/05c_de_ttest.R:22` | U |
| `de_table.id_col` | prot, rna | `R/modules/proteomics/03_mod_qc_post.R:81` | U |
| `de_table.pass_any_col` | prot, rna | **none** — only grep'd in `.Rmd` as a string-literal column name, never read from cfg | DU |
| `qc_post.enabled` | prot, rna | `R/modules/proteomics/03_mod_qc_post.R:20` | U |
| `qc_post.de_source` | prot, rna | `R/pipeline/proteomics/00_pipe_proteomics.R:98` | U |
| `qc_post.plots.volcano` | prot, rna | `R/modules/proteomics/03_mod_qc_post.R:33` | U |
| `qc_post.plots.ma` | prot, rna | `R/modules/proteomics/03_mod_qc_post.R:34` | U |
| `qc_post.outputs.write_de_tables` | prot, rna | `R/modules/proteomics/03_mod_qc_post.R:35` | U |
| `clustering.enabled` | prot, rna, met | `R/pipeline/metabolomics/templates/report_metabolomics.Rmd:448`, `R/modules/proteomics/04_mod_clustering.R` | U |
| `clustering.group_col` | prot, rna, met | `R/core/09_clustering.R:528` (`get_clustering_group_col`) | U |
| `clustering.de_source` | prot, rna, met | `R/domain/proteomics/90_config_validate.R:196` (validated) — passed through to feature gating | U (validator only) |
| `clustering.min_groups` | prot, rna, met | `R/core/09_clustering.R:619,637` | U |
| `clustering.steps.hierarchical.enabled` | prot, rna, met | `R/core/09_clustering.R` (hier_enabled) | U |
| `clustering.steps.hierarchical.distance` | prot, rna, met | `R/core/09_clustering.R`, validator | U |
| `clustering.steps.hierarchical.linkage` | prot, rna, met | `R/core/09_clustering.R`, validator | U |
| `clustering.steps.hierarchical.k` | prot, rna | validator (`R/domain/proteomics/90_config_validate.R:208`); accepted | U |
| `clustering.steps.partition.enabled` | prot, rna, met | `R/core/09_clustering.R` (`part_enabled`) | U |
| `clustering.steps.partition.algorithm` | prot, rna, met | `R/core/09_clustering.R`, validator | U |
| `clustering.steps.partition.distance` | prot, rna, met | `R/core/09_clustering.R`, validator | U |
| `clustering.steps.partition.k` | prot, rna | validator | U |
| `clustering.steps.partition.k_max` | prot, rna, met | validator + `R/core/09_clustering.R` | U |
| `clustering.steps.partition.nstart` | prot, rna, met | `R/core/09_clustering.R`, validator | U |
| `clustering.steps.partition.x_axis_col` | prot, rna, met | `R/core/09_clustering.R:1078` (`%||% group_col`) | UMC (always null in YAML) |
| `clustering.steps.partition.color_col` | prot, rna, met | `R/core/09_clustering.R:1077` | UMC (always null in YAML) |
| `clustering.steps.binary_patterns.enabled` | prot, rna, met | `R/core/09_clustering.R:631` | U |
| `clustering.steps.binary_patterns.corr_cutoff` | prot, rna, met | `R/core/09_clustering.R:1507` | U |
| `clustering.steps.binary_patterns.counts_cutoff` | prot | `R/core/09_clustering.R:1508` (proteomics alias for `_high`) | U |
| `clustering.steps.binary_patterns.counts_cutoff_high` | rna only | `R/core/09_clustering.R:1508` | U |
| `clustering.steps.binary_patterns.counts_cutoff_low` | rna only | `R/core/09_clustering.R:1509` | U |
| `clustering.heatmap.max_rows` | prot, rna, met | validator + `R/core/09_clustering.R` | U |
| `clustering.heatmap.cluster_cols` | prot, rna, met | validator + `R/core/09_clustering.R` | U |
| `effects.color` | prot, rna, met, multi | ≈60 sites — see section (b/4). Representative: `R/domain/proteomics/05_de_summary.R:209`, `R/modules/rnaseq/01_mod_qc_pre.R:45` | U |
| `effects.shape` | prot, rna, met | `R/modules/proteomics/01_mod_qc_pre.R:183`, `R/modules/rnaseq/01_mod_qc_pre.R:299` | U |
| `effects.samples` | prot, rna, met, multi | ≈50 sites — see section (b/4). Representative: `R/core/02_validation.R:33`, `R/domain/proteomics/04_preprocess.R:44` | U |
| `effects.heatmap_annotations` | met (only one that sets it; commented in rna) | `R/modules/rnaseq/01_mod_qc_pre.R:181`, `R/modules/rnaseq/03_mod_qc_post.R:164`, `R/modules/metabolomics/01_mod_qc_pre.R:299`, `R/modules/metabolomics/02_mod_qc_suite.R:590` | U |
| `report.project_id` | prot | **none** — no grep hit | DU |
| `report.introduction_text` | prot | `R/domain/proteomics/report_template_proteomics.Rmd:264` (via `report_cfg`) | C |
| `report.introduction_background` | prot | `report_template_proteomics.Rmd:265` | C |
| `report.introduction_aim` | prot | `report_template_proteomics.Rmd:266` | C |
| `report.experimental_design_file` | prot | `report_template_proteomics.Rmd:339` | C |
| `report.raw_files_location` | prot | `report_template_proteomics.Rmd:638` | C |
| `report.analysis_files_location` | prot | `report_template_proteomics.Rmd:639` | C |
| `report.full_results_table_link` | prot | `report_template_proteomics.Rmd:1178` | C |
| `report.expected_proteins_file` | prot | `report_template_proteomics.Rmd:692` | C |
| `report.show_pca` | prot | `report_template_proteomics.Rmd:216` | C |
| `report.show_density` | prot | `report_template_proteomics.Rmd:221` | C |
| `report.show_distance` | prot | `report_template_proteomics.Rmd:222` | C |
| `report.show_correlation` | prot | `report_template_proteomics.Rmd:225` | C |
| `report.show_expression_heatmap` | prot | `report_template_proteomics.Rmd:226` | C |
| `report.show_imputation_qc` | prot | `report_template_proteomics.Rmd:227` | C |
| `report.show_identified_proteins_venn` | prot | `report_template_proteomics.Rmd:1399` | C |
| `report.force_pca` | prot | `report_template_proteomics.Rmd:237` | C |
| `report.stats_methods_text` | prot | `report_template_proteomics.Rmd:1160` | C |
| `report.stats_section_text` | prot | `report_template_proteomics.Rmd:3749` | C |
| `report.disclaimer_text` | prot | `report_template_proteomics.Rmd:3807` | C |
| `excel.annotation_rows` | prot, rna, met | `R/modules/proteomics/05_mod_exports.R:114`, `R/domain/rnaseq/05_outputs_legacy.R:93`, `R/domain/metabolomics/05_outputs_legacy.R:155` | U |
| `excel.sample_label_cols` | prot, rna, met | `R/modules/proteomics/05_mod_exports.R:115`, `R/domain/rnaseq/05_outputs_legacy.R:94` | U |
| `annotation.organism` | prot, multi | `R/modules/rnaseq/08_mod_deconvolution.R:22`, `R/domain/proteomics/…` | U |
| `annotation.id_type` | prot, multi | `R/modules/rnaseq/05_mod_pathway.R:65` (via `ann_cfg$id_type`); proteomics consumer in `R/domain/proteomics/11_annotation.R` chain | U |
| `annotation.skip_annotation` | prot | `R/modules/rnaseq/05_mod_pathway.R:68,72`; not used in proteomics-specific code path | UMC (defined for prot, never read on prot path) |
| `annotation.custom_mapping_file` | prot | `R/modules/rnaseq/05_mod_pathway.R:73`; not used in proteomics-specific code path | UMC (defined for prot, never read on prot path) |
| `qc.run_umap` | prot | `R/modules/proteomics/01_mod_qc_pre.R:187` | U |
| `qc.umap_n_neighbors` | prot | `R/domain/proteomics/07_qc_advanced.R:22` | U |
| `qc.umap_min_dist` | prot | `R/domain/proteomics/07_qc_advanced.R:23` | U |
| `qc.outlier_sd_threshold` | prot | `R/domain/proteomics/07_qc_advanced.R:89` | U |
| `qc.min_sample_correlation` | prot | `R/domain/proteomics/07_qc_advanced.R:90` | U |
| `qc.top_de_heatmap_sizes` | prot | `R/modules/proteomics/03_mod_qc_post.R:101` | U |
| `qc.pca_subsets` | prot | `R/modules/proteomics/01_mod_qc_pre.R:55` | U |
| `pathway.enabled` | prot, rna, multi | `R/modules/proteomics/05_mod_pathway.R:12`, `R/modules/rnaseq/05_mod_pathway.R` | U |
| `pathway.method` | prot, rna, multi | `R/domain/rnaseq/10_pipeline_summary.R:188`, `R/modules/rnaseq/05_mod_pathway.R` | U |
| `pathway.databases` | prot, rna, multi | `R/domain/proteomics/report_template_proteomics.Rmd:451`, `R/domain/rnaseq/10_pipeline_summary.R:189`, validator | U |
| `pathway.gmt_file` | prot | `R/modules/rnaseq/05_mod_pathway.R:87` (`pw_cfg$gmt_file`); same pattern in `R/domain/proteomics/07b_pathway.R` | U |
| `pathway.min_size` | prot, rna | `R/modules/rnaseq/05_mod_pathway.R:107`, `R/domain/proteomics/07b_pathway.R:231` | U |
| `pathway.max_size` | prot, rna | `R/modules/rnaseq/05_mod_pathway.R:108`, `R/domain/proteomics/07b_pathway.R:232` | U |
| `pathway.gsea_ranking` | prot | `R/domain/proteomics/07b_pathway.R:47` | U |
| `pathway.pathway_volcano` | prot, rna | `R/modules/rnaseq/05_mod_pathway.R:139`, `R/domain/proteomics/07b_pathway.R:310` | U |
| `pathway.cluster_threshold` | prot | `R/domain/proteomics/07b_pathway.R:264` | U |
| `ppi.enabled` | prot, multi | `R/modules/proteomics/06_mod_ppi.R:14` | U |
| `ppi.species` | prot (commented), multi | `R/modules/proteomics/06_mod_ppi.R:81` | U |
| `ppi.significance_threshold` | prot, multi | `R/domain/proteomics/07c_ppi_networks.R:93` | U |
| `ppi.lfc_threshold` | prot | `R/domain/proteomics/07c_ppi_networks.R:94` | U |
| `ppi.string_score_threshold` | prot, multi | `R/domain/proteomics/07c_ppi_networks.R:147` | U |
| `ppi.active_subnetwork` | prot | `R/domain/proteomics/07c_ppi_networks.R:178` | U |
| `ppi.complex_analysis` | prot | `R/domain/proteomics/07c_ppi_networks.R:187` | U |
| `domain_analysis.enabled` | prot | `R/modules/proteomics/06b_mod_domains.R:7` | U |
| `advanced_stats.enabled` | prot | `R/modules/proteomics/07_mod_advanced_stats.R:14` | U |
| `advanced_stats.compute_effect_size_ci` | prot | `R/domain/proteomics/07d_advanced_stats.R:63` | U |
| `advanced_stats.bootstrap_n` | prot | validator + `R/domain/multiomics/11_stability_analysis.R:51,88` (multiomics, not prot) | UMC (prot value never reaches a consumer) |
| `advanced_stats.ci_level` | prot | `R/domain/multiomics/11_stability_analysis.R:730` (multiomics, not prot) | UMC (prot value never reaches a consumer) |
| `advanced_stats.run_robust_regression` | prot | `R/domain/proteomics/07d_advanced_stats.R:79` | U |
| `advanced_stats.run_power_analysis` | prot | `R/domain/proteomics/07d_advanced_stats.R:87` | U |
| `advanced_stats.random_effects` | prot | `R/domain/proteomics/07d_advanced_stats.R:277,301` | U |
| `technical_report.*` (facility, sample_prep, ms_acquisition, search_engine, search_parameters, acknowledgment) | prot | `R/domain/proteomics/10_pipeline_summary.R:23` (block read) | U (as a whole list) |
| `commentary.enabled` | prot, rna, multi | `R/modules/rnaseq/06_mod_commentary.R:17,21`, `R/modules/multiomics/07_mod_consensus.R:145,147` | U |
| `commentary.backend` | prot, rna, multi | passed through to `generate_all_commentary()` in `R/services/12_commentary.R` | U |
| `commentary.claude_code_model` | prot | `R/services/12_commentary.R:678` | U |
| `commentary.max_tokens` | prot | `R/services/12_commentary.R:778` | U |
| `commentary.max_retries` | prot | `R/services/12_commentary.R:785` | U |

### RNA — `modes.rna.*` (delta from proteomics; shared keys above)

| key_path | configs_defining_it | R_files_using_it | status |
|----------|----------------------|------------------|--------|
| `files.counts` | rna, multi | `R/pipeline/rnaseq/00_pipe_rnaseq.R`, `R/domain/rnaseq/00_inputs.R` | U |
| `files.sample_map` | rna | `R/domain/rnaseq/03_preprocess.R:94` (`map_from`/`map_to` lookup) | U |
| `files.metadata` | rna, multi | `R/domain/rnaseq/90_config_validate.R:36` | U |
| `files.contrasts` | rna, multi | (loaded via `R/core/01_io.R`) | U |
| `files.annotation` | rna | `R/modules/rnaseq/05_mod_pathway.R` (annotation chain) | U |
| `files.trinotate` | rna | `R/domain/rnaseq/01b_annotation.R:228` (Trinotate path) | U |
| `annotation.source` | rna | `R/domain/rnaseq/90_config_validate.R:73` (validated); consumer in `R/domain/rnaseq/01b_annotation.R` | U |
| `id_columns.gene_id` | rna | `R/domain/rnaseq/00_inputs.R:54`, `R/domain/rnaseq/03_preprocess.R:19,69` | U |
| `id_columns.sample_col` | rna | `R/domain/rnaseq/00_inputs.R:46,55`, `R/domain/rnaseq/03_preprocess.R:5`, `R/domain/rnaseq/04_de_summary.R:40,311` | U |
| `id_columns.map_from` | rna | `R/domain/rnaseq/03_preprocess.R:94` | U |
| `id_columns.map_to` | rna | `R/domain/rnaseq/03_preprocess.R:95` | U |
| `id_columns.gene_id_col` | multi only (rna mode in multiomics YAML) | **none** — R reads `gene_id`, not `gene_id_col` | DU (multiomics-only typo) |
| `normalization.method` | rna, multi | `R/domain/rnaseq/01_expression.R:18`, validator | U |
| `normalization.prior.count` | rna | `R/domain/rnaseq/03_preprocess.R:222`, validator, `R/domain/rnaseq/10_pipeline_summary.R:182` | U (rename — see (i)) |
| `filtering.group_col` | rna, multi | `R/domain/rnaseq/03_preprocess.R:162`, `R/modules/rnaseq/02_mod_de.R:23`, `R/services/12_commentary.R:185`, `R/core/04_config.R:99` | U |
| `filtering.auto_filter.min` | rna | only in `R/domain/rnaseq/report_template.Rmd:294,419,3005` — **not passed into `run_auto_filter_pipeline()`** at `03_preprocess.R:190` → consumer never sees the value (see (f)) | UMC (report-only) |
| `filtering.auto_filter.max` | rna | same as above; report-only | UMC (report-only) |
| `filtering.auto_filter.fallback` | rna | same as above; report-only | UMC (report-only) |
| `filtering.mode` | multi (rna mode) | `R/domain/rnaseq/03_preprocess.R:181,184` (`filter_mode == "deseq2_only"` / `"fixed"`) | U |
| `filtering.fixed_threshold` | (none — code only) | `R/domain/rnaseq/03_preprocess.R:185` | UMC (code-side only) |
| `de.method` | rna, multi | `R/domain/rnaseq/04_de_summary.R`, `report_template.Rmd:296` | U |
| `de.deseq_mode` | rna | `R/domain/rnaseq/report_template.Rmd:297,320`, `R/modules/rnaseq/06_mod_commentary.R:121` | U |
| `de.p_cutoff` | rna, multi | `R/modules/rnaseq/04_mod_clustering.R:85`, `report_template.Rmd:421` | U |
| `de.linear_fc_cutoff` | rna, multi | `R/domain/rnaseq/04_de_summary.R:152` | U |
| `de.use_adj` | rna | `R/domain/rnaseq/report_template.Rmd:300` | U |
| `qc_post.max_top_de_features` | rna | `R/modules/rnaseq/03_mod_qc_post.R:128` | U |
| `qc.adaptive_plots` | rna | `R/modules/rnaseq/01_mod_qc_pre.R:29` | U |
| `qc.thresholds.min_samples_for_pca3d` | rna, met | `R/modules/rnaseq/01_mod_qc_pre.R:38` | U |
| `qc.thresholds.max_samples_for_heatmaps` | rna, met | `R/modules/rnaseq/01_mod_qc_pre.R:39`, `R/modules/metabolomics/01_mod_qc_pre.R:206` | U |
| `qc.thresholds.max_samples_for_expr_heatmap` | rna | `R/modules/rnaseq/01_mod_qc_pre.R:40` | U |
| `qc.heatmap_viz.show_sample_labels` | rna | `R/modules/rnaseq/01_mod_qc_pre.R:189` | U |
| `qc.heatmap_viz.cluster_samples` | rna | `R/modules/rnaseq/01_mod_qc_pre.R:190` | U |
| `qc.heatmap_viz.adjust_correlation_scale` | rna | `R/modules/rnaseq/01_mod_qc_pre.R:191` | U |

### Metabolomics — `modes.metabolomics.*` (delta)

| key_path | configs_defining_it | R_files_using_it | status |
|----------|----------------------|------------------|--------|
| `input.format` | met | `R/domain/metabolomics/00_inputs.R:128`, `R/pipeline/metabolomics/00_pipe_metabolomics.R` | U |
| `input.level_format` | met | `R/modules/metabolomics/00_mod_preprocessing.R:74` | U |
| `input.level_pattern` | met | `R/pipeline/metabolomics/00_pipe_metabolomics.R:59` | U |
| `input.sheet` | met | `R/domain/metabolomics/00_inputs.R` (passed to `read_metab_file`) | U |
| `files.data` | met | `R/domain/metabolomics/00_inputs.R:132` | U |
| `files.data_dir` | met (commented) | `R/domain/metabolomics/00_inputs.R:133` | U |
| `files.metadata` | met | `R/domain/metabolomics/00_inputs.R` | U |
| `files.contrasts` | met | (loaded via `R/core/01_io.R`) | U |
| `parsing.cd_area_prefix` | met | `R/domain/metabolomics/00_inputs.R:215` | U |
| `parsing.cd_sample_regex` | met | `R/domain/metabolomics/00_inputs.R:232` | U |
| `id_columns.name_col` | met | `R/domain/metabolomics/00_inputs.R:654` | U |
| `id_columns.mz_col` | met | `R/domain/metabolomics/00_inputs.R:655` | U |
| `id_columns.rt_col` | met | `R/domain/metabolomics/00_inputs.R:656` | U |
| `id_columns.map_from` | met | `R/modules/metabolomics/00_mod_preprocessing.R:68`, validator | U |
| `id_columns.map_to` | met | `R/modules/metabolomics/00_mod_preprocessing.R:69`, validator | U |
| `preprocessing.feat_missing_threshold` | met | `R/modules/metabolomics/00_mod_preprocessing.R:254` | U |
| `preprocessing.sample_missing_threshold` | met | `R/modules/metabolomics/00_mod_preprocessing.R:255` | U |
| `preprocessing.chosen_norm` | met | `R/modules/metabolomics/00_mod_preprocessing.R:377`, `_targets.R:130` | U |
| `preprocessing.drift_correction.enabled` | met | `R/modules/metabolomics/00_mod_preprocessing.R:481` (`dc_cfg <- pre_cfg$drift_correction`) | U |
| `preprocessing.drift_correction.injection_order_col` | met | `R/domain/metabolomics/10_drift_correction.R` | U |
| `preprocessing.drift_correction.qc_flag_col` | met | `R/domain/metabolomics/10_drift_correction.R` | U |
| `preprocessing.drift_correction.loess_span` | met | `R/domain/metabolomics/10_drift_correction.R` | U |
| `preprocessing.drift_correction.epsilon` | met | `R/domain/metabolomics/10_drift_correction.R` | U |
| `normalization.input_already_normalized` | met | **none** — no grep hit | DU |
| `normalization.sample_norm` | met | validator only (`R/domain/metabolomics/00_inputs.R:163`); actual normalization is driven by `preprocessing.chosen_norm`. Legacy reference per comment in YAML lines 95-101. | DU (legacy reference per comment) |
| `normalization.transform` | met | validator only; actual transform applied via `transform_metab(..., method = "log2")` calls (hard-coded `log2`) — see `R/modules/metabolomics/00_mod_preprocessing.R:290,478,546,583,627` | DU (legacy reference per comment) |
| `normalization.scaling` | met | validator only; YAML comment says "retained for backward compatibility only" | DU (legacy reference per comment) |
| `normalization.pseudocount` | met | `R/modules/metabolomics/00_mod_preprocessing.R:287,478,570,608` | U |
| `normalization.na_policy` | met | validator only (`R/domain/metabolomics/00_inputs.R:172`); never read elsewhere | DU |
| `normalization.is_ref_col` | met | **none** — only `normalize_samples()` (`R/domain/metabolomics/01_normalization.R:25`) takes a `ref_col` parameter, but `normalize_samples()` is never called anywhere | DU |
| `normalization.biological_factor_col` | met | same — `normalize_samples()` is dead code | DU |
| `sample_filter.enabled` | met | `R/domain/metabolomics/00_inputs.R:730` (`get_sample_filter_rules_metab`) | U |
| `sample_filter.rules.exclude_blanks` | met | `R/core/03_alignment.R:120` (`apply_sample_filter`) | U |
| `sample_filter.rules.exclude_qc` | met | same | U |
| `sample_filter.rules.exclude_samples` | met | same | U |
| `sample_filter.exclude_after_norm` | (none — code-side only) | `R/modules/metabolomics/00_mod_preprocessing.R:398` | UMC (code-side only) |
| `qc.enabled` | met | `R/modules/metabolomics/02_mod_qc_suite.R` (gate) | U |
| `qc.qc_flag_column` | met | `R/modules/metabolomics/02_mod_qc_suite.R:315` | U |
| `qc.missingness_heatmap` | met | **none** — no grep hit | DU |
| `qc.thresholds.min_samples_for_pca3d` | met | `R/modules/metabolomics/01_mod_qc_pre.R` | U |
| `qc.thresholds.max_samples_for_heatmaps` | met | `R/modules/metabolomics/02_mod_qc_suite.R:316` | U |
| `de.method` | met | `R/modules/metabolomics/07_mod_commentary.R:251`, `R/pipeline/metabolomics/templates/report_metabolomics.Rmd:474` | U |
| `de.condition_column` | met | `R/modules/metabolomics/05_mod_enrichment.R:143`, `R/domain/metabolomics/03_differential.R:195` | U |
| `de.p_cutoff` | met | `R/modules/metabolomics/04b_mod_network.R:47` | U |
| `de.pval_cutoff` | met | **none** — no grep hit (duplicate of p_cutoff in the same section) | DU |
| `de.linear_fc_cutoff` | met | `R/modules/metabolomics/04_mod_clustering.R:73`, `report_metabolomics.Rmd:392,476,533,836` | U |
| `de.robust_ebayes` | met | **none** — no grep hit | DU |
| `clustering.*` | met | identical to proteomics/rna paths in `R/core/09_clustering.R` | U |
| `enrichment.run_enrichment` | met | (gate exists upstream; mod selector in `R/modules/metabolomics/05_mod_enrichment.R`) | U |
| `enrichment.libraries` | met | `R/modules/metabolomics/05_mod_enrichment.R` (lib list passthrough) | U |
| `enrichment.gmt_file` | met | `R/modules/metabolomics/05_mod_enrichment.R` (gmt path passthrough) | U |
| `enrichment.mapping_file` | met | `R/modules/metabolomics/05_mod_enrichment.R` (compound→KEGG mapping passthrough) | U |
| `enrichment.pvalue_column` | met | `R/domain/metabolomics/06_enrichment.R` | U |

### Multiomics — `modes.multiomics.*` (delta)

| key_path | configs_defining_it | R_files_using_it | status |
|----------|----------------------|------------------|--------|
| `gene_protein_mapping_file` | multi | `R/domain/multiomics/01_id_mapping.R:20` | U |
| `require_one_to_one_mapping` | multi | `R/modules/multiomics/01_mod_harmonization.R:42` | U |
| `feature_selection.method` | multi | validator default + `R/domain/multiomics/report_template_multiomics.Rmd:224` | U |
| `feature_selection.top_n` | multi | validator default + `R/domain/multiomics/03_feature_selection.R`, `report_template_multiomics.Rmd:225` | U |
| `feature_selection.min_var_percentile` | multi | `R/domain/multiomics/03_feature_selection.R:19,49,130` | U |
| `integration.methods` | multi | `R/modules/multiomics/02_mod_integration.R:32`, validator | U |
| `integration.diablo.ncomp` | multi | validator (default 2) + `R/domain/multiomics/04_integration_diablo.R` | U |
| `integration.diablo.design_matrix` | multi | `R/domain/multiomics/04_integration_diablo.R:22` | U |
| `integration.diablo.cv_folds` | multi | `R/domain/multiomics/04_integration_diablo.R:23` | U |
| `integration.snf.K` | multi | validator + `R/domain/multiomics/05_integration_snf.R` | U |
| `integration.snf.alpha` | multi | `R/domain/multiomics/05_integration_snf.R:22` | U |
| `integration.snf.T` | multi | `R/domain/multiomics/05_integration_snf.R:23` | U |
| `integration.snf.n_clusters` | multi | validator + `R/domain/multiomics/05_integration_snf.R`, `R/domain/multiomics/10_integration_consensus.R:948` | U |
| `condition_column` | multi | `R/domain/multiomics/04_integration_diablo.R:35`, `R/domain/multiomics/90_config_validate.R:172` | U |
| `enrichment.run_enrichment` | multi | `R/modules/multiomics/04_mod_enrichment.R:25` | U |

---

## (b) Verdicts on the four prior suspicions

### 1. `clustering.min_groups` — does it gate anything real?

**It gates a real branch.** Read at `R/core/09_clustering.R:619`:

```r
min_groups <- cl$min_groups %||% 3L
```

Used at lines 637–642:

```r
can_multi_group <- isTRUE(n_groups >= as.integer(min_groups))

list(
  hierarchical    = hier_enabled,
  partition       = part_enabled,
  binary_patterns = isTRUE(bin_enabled && can_multi_group)
)
```

So when the design has fewer than `min_groups` groups, `binary_patterns`
is force-disabled (because the binary-pattern algorithm requires multiple
groups to compare). It is **not** a no-op `if()`.

A separate `filtering.min_groups` (proteomics-only, default 1) gates a
totally different thing — feature-level filtering at
`R/domain/proteomics/02_filtering.R:35,37` (passes the value to
`pass_filter()`). That one is also live.

**Verdict: KEEP** both. They are different parameters in different
sections that happen to share a name. The audit recommends in section
(h)/(i) that they get an inline comment so this is obvious to a
first-time reader.

### 2. `effects` section — extract to shared defaults, or eliminate?

The suspicion is that `effects` duplicates values that already live
elsewhere. The data say otherwise.

**`effects.samples` is the de-facto canonical sample-column source.**
Approximately 50 R sites read it; only ~14 read `id_columns.sample_col`,
and almost always behind a `%||%` chain where `effects.samples` is
checked first. Examples:

- `R/core/02_validation.R:33` — `sample_col <- cfg$effects$samples %||% cfg$id_columns$sample_col`
- `R/domain/proteomics/04_preprocess.R:44` — same pattern
- `R/domain/proteomics/03_imputation.R:289` — `p$effects$samples %||% p$id_columns$sample_col %||% "SampleID"`
- `R/modules/proteomics/05_mod_exports.R:100-101` — same

**Metabolomics has no `id_columns.sample_col` at all.** The metabolomics
parser at `R/domain/metabolomics/00_inputs.R:311` uses
`cfg$effects$samples %||% "sample_id"` as the only source. So
`effects.samples` is not removable: it is the only sample column source
for metabolomics today.

**`effects.color` is more entangled than just a duplicate of `group_col`.**
It is the *primary* grouping value (`get_color_config()` at
`R/core/04_config.R:48–54`); `resolve_group_col()` at lines 93–108 walks
`filtering.group_col → de_table.group_col → effects.color[1] → "Condition"`.
About 60 R sites read it, including all PCA/QC/heatmap colouring,
report templates, and commentary. It carries semantics that
`filtering.group_col` does not (it can be a vector to generate multiple
PCAs).

**`effects.shape` and `effects.heatmap_annotations`** have no
counterpart in `id_columns` / `filtering` / `clustering`. They are not
duplicates of anything.

**Option (a) — extract to shared defaults** is not strictly possible:
`effects.color` is `"Condition"` in proteomics but `"line"` in RNA and
`"treatment"` in metabolomics. They are NOT identical across omics.
`effects.samples` *is* `"SampleID"` in proteomics and RNA, and
`"SampleID"` in metabolomics — but metabolomics uses lowercase
fallbacks (`"sample_id"`) all through its code. Pulling them to a
shared file would require resolving that mismatch first.

**Option (b) — eliminate `effects` entirely** is not possible: the
fields are *not* covered by `id_columns` / `filtering` / `clustering`.
`effects.shape`, `effects.heatmap_annotations`, vector-valued
`effects.color` have no other home. `effects.samples` cannot move to
`id_columns.sample_col` without writing the metabolomics consumer
from scratch.

**Verdict: KEEP `effects` as a section. SAFE AFTER R REFACTOR** (a
minor one) — collapse the proteomics three-way `effects.samples` /
`id_columns.sample_col` / `id_columns.map_to` repetition discussed in
suspicion 4 below. Plus: drop the proteomics validator's hard rule
that `id_columns.sample_col` must equal `id_columns.map_to`
(`R/domain/proteomics/90_config_validate.R:19–24`) once that
collapse is done.

### 3. `partition.x_axis_col` and `partition.color_col` — non-null anywhere?

**Always null.** Verified in all four YAMLs:

- `config/templates/proteins_config.yaml:161–162`
- `config/templates/rna_config.yaml:128–129`
- `config/templates/metabolomics_template.yaml:206–207`

And the only R consumer is `R/core/09_clustering.R:1077–1078`:

```r
color_col  <- cfg$clustering$steps$partition$color_col  # NULL if not set
x_axis_col <- cfg$clustering$steps$partition$x_axis_col %||% group_col
```

`group_col` here is from `get_clustering_group_col(cfg, meta)`
(`R/core/09_clustering.R:527`), which already reads
`cfg$clustering$group_col` — so the override mechanism duplicates what
the user has already declared one level up.

To remove safely, replace those two lines at `R/core/09_clustering.R:1077–1078` with:

```r
color_col  <- NULL
x_axis_col <- group_col
```

And remove the two warning blocks at lines 1080–1089 that re-validate
the override. The `build_cluster_long_df()` signature
(`R/core/09_clustering.R:557–559`) keeps `color_col = NULL` as a
function parameter — that's fine; it's just no longer a config-driven
parameter.

**Verdict: SAFE TO REMOVE.** Only `R/core/09_clustering.R` changes (8
lines).

### 4. `SampleID` repetition (proteomics: `id_columns.sample_col` + `id_columns.map_to` + `effects.samples`)

Each of the three is read by different parts of the code:

**`id_columns.sample_col`** (`"SampleID"`):

- `R/core/02_validation.R:33` — `sample_col <- cfg$effects$samples %||% cfg$id_columns$sample_col`
- `R/domain/proteomics/90_config_validate.R:16` — required scalar check
- `R/domain/proteomics/90_config_validate.R:21` — enforces `sample_col == map_to`
- `R/domain/proteomics/04_preprocess.R:44` — same `%||%` fallback as core
- `R/domain/proteomics/03_imputation.R:289` — `p$effects$samples %||% p$id_columns$sample_col %||% "SampleID"`
- `R/modules/proteomics/05_mod_exports.R:101` — `prot_cfg$effects$samples %||% prot_cfg$id_columns$sample_col`

**`id_columns.map_to`** (`"SampleID"`):

- `R/domain/proteomics/90_config_validate.R:20` — scalar check
- `R/domain/proteomics/90_config_validate.R:21` — enforces equality with `sample_col`
- `R/domain/proteomics/00_inputs.R:53,60` — used by `check_has_cols` / `check_all_in` to validate the `sample_map` CSV
- `R/domain/proteomics/01_expression.R:240` — passed to `apply_sample_map_to_colnames(..., id_cols$map_to)` to rename matrix columns
- `R/core/03_alignment.R:42,47` — generic mapper

**`effects.samples`** (`"SampleID"`):

- ≈50 sites, all reading `cfg$effects$samples %||% "SampleID"`.

So `effects.samples` is the dominant, near-universal read; `id_columns.sample_col` is only consulted as a fallback, and `id_columns.map_to` is the *renaming target* for the `sample_map` join, semantically distinct from the others.

**Recommendation: collapse to two keys, not one.**

- Keep `id_columns.map_from` and `id_columns.map_to` only when a
  `sample_map` file is provided. Their job is "rename raw column names
  to unified IDs" — they're about *the sample_map file*, not about the
  sample column in metadata.
- Keep `effects.samples` (or rename it `id_columns.sample_col` — see
  (j) for the open question), and use it everywhere.
- Drop the parallel `id_columns.sample_col` field.

**R changes required (call them M1…M6):**

- M1 `R/core/02_validation.R:33–34` — drop the `%||% cfg$id_columns$sample_col` tail.
- M2 `R/domain/proteomics/90_config_validate.R:16` — remove the `id_columns$sample_col` assertion.
- M3 `R/domain/proteomics/90_config_validate.R:19–24` — remove the "must equal `map_to`" rule.
- M4 `R/domain/proteomics/04_preprocess.R:44` — drop the `%||% cfg$id_columns$sample_col` tail.
- M5 `R/domain/proteomics/03_imputation.R:289` — drop the `%||% p$id_columns$sample_col` tail.
- M6 `R/modules/proteomics/05_mod_exports.R:100–101` — drop the `%||% prot_cfg$id_columns$sample_col` tail.

That's six edits across five files. The metabolomics path is already
single-source (`effects.samples`); the RNA `id_columns.sample_col` is
read at six sites (`R/domain/rnaseq/00_inputs.R:46,55`,
`03_preprocess.R:5`, `04_batch_correction.R:35`,
`04_de_summary.R:40,311`) and would need parallel edits — about six
more.

**Verdict: SAFE AFTER R REFACTOR.** Total edits ≈ 12 across ~7 files.
The rule-of-3 cap in CLAUDE.md is exceeded, so this should be a
dedicated PR with an explicit plan; that's question (j-3) below.

---

## (c) Cross-omic naming conflicts — same concept, different keys

### C1. `binary_patterns.counts_cutoff` (prot) vs `counts_cutoff_high` + `counts_cutoff_low` (rna)

**Real conflict.** Proteomics YAML uses bare `counts_cutoff`. RNA splits
into `counts_cutoff_high` and `counts_cutoff_low`. The R code at
`R/core/09_clustering.R:1507–1509` handles both via a fallback:

```r
corr_cutoff        = bcfg$corr_cutoff %||% 0.8,
counts_cutoff_high = bcfg$counts_cutoff_high %||% bcfg$counts_cutoff %||% 0,
counts_cutoff_low  = bcfg$counts_cutoff_low %||% NULL,
```

The validator (`R/domain/proteomics/90_config_validate.R:230–232`)
accepts all three names, but no validator rejects setting both `counts_cutoff` and `counts_cutoff_high`.

**Recommended canonical name:** `counts_cutoff_high` (RNA's pair).
Rationale: the algorithm has always been *dual-threshold* (a high cutoff
for "on" positions, a low cutoff for "off" positions); the bare
`counts_cutoff` was just an early single-threshold alias. RNA's two-key
form is the more accurate description.

**Files to change:**

- `config/templates/proteins_config.yaml:166` — rename `counts_cutoff: 0` → `counts_cutoff_high: 0`.
- `R/core/09_clustering.R:1508` — drop the `%||% bcfg$counts_cutoff` legacy tail (just `counts_cutoff_high = bcfg$counts_cutoff_high %||% 0`).
- `R/domain/proteomics/90_config_validate.R:230` — remove `assert_scalar_num(b$counts_cutoff, ...)`.
- `run.R:1095–1097` — template emitter, no `counts_cutoff` written today; verify nothing else writes the legacy name (a fresh `grep -rn '\bcounts_cutoff\b'` shows no other writers).

### C2. `binary_patterns.group_col` — written by template emitter, never read

**Stale name.** Neither the proteomics YAML nor the RNA YAML actually
declares `binary_patterns.group_col` or `binary_patterns.group`. But:

- `run.R:676` emits `group_col: "%s"` *inside* the `binary_patterns:`
  block when it scaffolds a new config — this writes a key R never reads.
- `R/core/09_clustering.R:7` (docstring): "Groups are defined by
  `cfg$clustering$steps$binary_patterns$group_col`" — stale docstring;
  the actual read is at line 88 (`get_clustering_group_col(cfg, meta)`),
  which goes to `cfg$clustering$group_col` (top-level), not the
  binary_patterns sub-section.

So the user's suspicion phrasing ("prot vs rna" for this key) doesn't
match the YAML — neither omic actually has it. But the docstring and
the template emitter are misleading.

**Files to change:**

- `R/core/09_clustering.R:7` — fix docstring to say
  "Groups are defined by `cfg$clustering$group_col`".
- `run.R:676` — remove the `group_col: "%s"` line from the templated
  `binary_patterns:` block; the YAML emitter should write it at the
  `clustering:` level if at all.

### C3. `normalization.prior.count` (rna) — dotted YAML key

**Not strictly a conflict** (only one omic has this concept), but the
dotted YAML key is non-standard and confusing. `yaml::read_yaml`
preserves it as a single string-named field, R reads it correctly via
`cfg$normalization$prior.count` (which becomes `[["prior.count"]]`).

This compiles and works today (`R/domain/rnaseq/03_preprocess.R:222`,
`R/domain/rnaseq/90_config_validate.R:43`,
`R/domain/rnaseq/10_pipeline_summary.R:182`). But it reads like a
nested object and isn't — anyone editing the YAML expects
`normalization.prior` to be a section, with `count` as a leaf.

**Recommended canonical name:** `prior_count` (snake_case, no dot).

**Files to change:**

- `config/templates/rna_config.yaml:57` — rename.
- `R/domain/rnaseq/03_preprocess.R:222` — `cfg$normalization$prior_count`.
- `R/domain/rnaseq/90_config_validate.R:43` — `n$prior_count`, message string.
- `R/domain/rnaseq/10_pipeline_summary.R:182` — `rna_cfg$normalization$prior_count`.

### C4. `id_columns.gene_id_col` (multi/rna) vs `id_columns.gene_id` (rna)

Found while writing this audit. `config/templates/multiomics_config.yaml:34`
declares the RNA mode's gene-id field as `gene_id_col: "gene_id"`, but
every R consumer reads `cfg$id_columns$gene_id` (e.g.
`R/domain/rnaseq/00_inputs.R:54`, `03_preprocess.R:19,69`,
`R/domain/rnaseq/90_config_validate.R:21,23`,
`R/domain/multiomics/01_id_mapping.R:103`).

**Real bug: this key is currently dead.** When a user copies the
multiomics template and runs it, the RNA gene_id falls through to the
validator's `allow_null = TRUE` path and downstream code uses the
hard-coded `"gene_id"` fallback. Today that happens to be the correct
value, masking the bug.

**Recommended canonical name:** `gene_id` (matches single-omic RNA template).

**Files to change:**

- `config/templates/multiomics_config.yaml:34` — rename `gene_id_col` → `gene_id`.

(No R changes needed — they already read `gene_id`.)

---

## (d) Dead parameters — declared in YAML, never read

These are YAML keys with no matching read anywhere in `R/`, `_targets.R`,
`run.R`, `configure_wizard.R`, or `validate_multiomics.R`. Each was
confirmed by `grep -rn '<leaf_name>' R/`.

| YAML key | File:line | Notes |
|----------|-----------|-------|
| `report.project_id` | `config/templates/proteins_config.yaml:184` | Comment says "Shown in report title", but no R/Rmd reads `report_cfg$project_id` |
| `de_table.pass_any_col` | `config/templates/proteins_config.yaml:119`, `rna_config.yaml:94` | Only appears in `.Rmd` as a hard-coded grep pattern for `"pass_any"` column names in DE outputs — never read from `cfg$de_table$pass_any_col` |
| `id_columns.gene_id_col` | `config/templates/multiomics_config.yaml:34` | R reads `gene_id`, not `gene_id_col` — see C4 above |
| `design.batch_column` | `config/templates/multiomics_config.yaml:198` | No grep hit |
| `global.metadata` | `config/templates/multiomics_config.yaml:210` | Each omic loads its own metadata path |
| `output.output_dir` | `config/templates/multiomics_config.yaml:213` | Only appears in docstring comments (`R/domain/proteomics/07c_ppi_networks.R:7`, `07d_advanced_stats.R:8`) — code uses `paths.out` |
| `output.report_format` | `config/templates/multiomics_config.yaml:214` | No grep hit |
| `normalization.input_already_normalized` | `metabolomics_template.yaml:91` | YAML comment admits "skip sample_norm but still allow transform/scaling" — but no code reads it |
| `normalization.sample_norm` | `metabolomics_template.yaml:102` | YAML self-documents as "legacy / reference only"; only the validator at `R/domain/metabolomics/00_inputs.R:163` reads it (rejects invalid values). The actual normalization branch is selected by `preprocessing.chosen_norm` |
| `normalization.transform` | `metabolomics_template.yaml:106` | Validator-only; pipeline calls `transform_metab(..., method = "log2", ...)` with a hardcoded method |
| `normalization.scaling` | `metabolomics_template.yaml:113` | YAML self-documents as "retained for backward compatibility only" |
| `normalization.na_policy` | `metabolomics_template.yaml:120` | Validator-only |
| `normalization.is_ref_col` | `metabolomics_template.yaml:123` | `normalize_samples()` (`R/domain/metabolomics/01_normalization.R:25`) takes a `ref_col` parameter, but `normalize_samples()` is never called |
| `normalization.biological_factor_col` | `metabolomics_template.yaml:128` | Same — `normalize_samples()` is dead code |
| `de.pval_cutoff` | `metabolomics_template.yaml:183` | Adjacent to `de.p_cutoff` (line 182) which IS read; this one is a copy-paste duplicate |
| `de.robust_ebayes` | `metabolomics_template.yaml:185` | No grep hit |
| `qc.missingness_heatmap` | `metabolomics_template.yaml:167` | No grep hit |

---

## (e) Missing parameters — keys the R code reads (often with a `%||%`
default) that are absent from at least one config

The risk class here is **silent NULL** — the consumer falls through to a
hard-coded default, the user never sees a warning, and the YAML they
edit is misleading about what knobs exist.

| key_path | Defined in | Missing in | R consumer (representative) | Default |
|----------|-----------|------------|-----------------------------|---------|
| `filtering.group_col` | rna | **prot, met** | `R/domain/proteomics/02_filtering.R:24` (proteomics chain reads it before falling back to `de_table.group_col`) | `effects$color %||% "Condition"` |
| `de_table.group_col` | (none — code-side only) | all | `R/domain/proteomics/02_filtering.R:24`, `R/domain/rnaseq/03_preprocess.R:162`, `R/domain/rnaseq/05_deconvolution.R:70` | `effects$color %||% "Condition"` |
| `effects.heatmap_annotations` | met (commented in rna) | **prot** | `R/modules/rnaseq/01_mod_qc_pre.R:181`, `R/modules/rnaseq/03_mod_qc_post.R:164` (not currently invoked from proteomics — falls back to colour var) | `NULL` |
| `annotation.skip_annotation` | prot | rna, met, multi | `R/modules/rnaseq/05_mod_pathway.R:68,72` (RNA path only — the proteomics declaration is read by no proteomics consumer) | `FALSE` |
| `annotation.custom_mapping_file` | prot | rna, met, multi | `R/modules/rnaseq/05_mod_pathway.R:73` (RNA only) | `NULL` |
| `advanced_stats.bootstrap_n` | prot | (n/a — multiomics-only consumer) | `R/domain/multiomics/11_stability_analysis.R:51,88` | 100 |
| `advanced_stats.ci_level` | prot | (n/a — multiomics-only consumer) | `R/domain/multiomics/11_stability_analysis.R:730` | 0.95 |
| `clustering.steps.partition.k_method` | rna only | prot, met | code path exists (`R/core/09_clustering.R`, gap-statistic branch) but only the RNA YAML documents it | `"silhouette"` |
| `clustering.steps.partition.gap_B` | rna only | prot, met | same | 10 |
| `qc_post.max_top_de_features` | rna | prot, met | `R/modules/rnaseq/03_mod_qc_post.R:128` — RNA-only consumer | 100 |
| `qc.adaptive_plots` | rna | prot, met | `R/modules/rnaseq/01_mod_qc_pre.R:29` — RNA-only consumer | varies |
| `qc.heatmap_viz.*` | rna | prot, met | `R/modules/rnaseq/01_mod_qc_pre.R:189–191` — RNA-only consumer | varies |
| `qc.thresholds.max_samples_for_expr_heatmap` | rna | met | `R/modules/rnaseq/01_mod_qc_pre.R:40` | 60 |
| `report.show_pca_3d` | (none in YAML) | prot | `report_template_proteomics.Rmd:218` | TRUE |
| `report.show_pca_subsets` | (none in YAML) | prot | `report_template_proteomics.Rmd:219` | TRUE |
| `report.show_pca_topvar` | (none in YAML) | prot | `report_template_proteomics.Rmd:220` | FALSE |
| `report.show_umap` | (none in YAML) | prot | `report_template_proteomics.Rmd:217` | FALSE |
| `filtering.mode` | (none in YAML) | rna | `R/domain/rnaseq/03_preprocess.R:181,184` (`"deseq2_only"` / `"fixed"`) | not gated |
| `filtering.fixed_threshold` | (none in YAML) | rna | `R/domain/rnaseq/03_preprocess.R:185` | 1.0 |
| `sample_filter.exclude_after_norm` | (none in YAML) | met | `R/modules/metabolomics/00_mod_preprocessing.R:398` | FALSE |
| `pathway.method`, `pathway.databases` | prot, rna, multi | (consistent) | RNA reads them — proteomics YAML defines them, validator accepts; no `prot` consumer reads `pathway.method` directly (proteomics uses `fgsea+ora` hard-wired in `R/domain/proteomics/07b_pathway.R`) | n/a |

---

## (f) Suspicious values

### F1. RNA `binary_patterns.enabled: true` but no group_col under binary_patterns
**Confirmed but not a bug.** `R/core/09_clustering.R:630–634`:

```r
bin_cfg <- steps$binary_patterns %||% list()
bin_enabled <- isTRUE(bin_cfg$enabled %||% FALSE)
bin_group_col <- cl$group_col
# Only enable if both enabled flag is TRUE AND group_col is provided (non-NULL)
bin_enabled <- isTRUE(bin_enabled && !is.null(bin_group_col))
```

The `group_col` is taken from `clustering.group_col` (the parent
section), NOT from `binary_patterns.group_col`. So `enabled: true` works
because the RNA YAML provides `clustering.group_col: "line"`
(`rna_config.yaml:110`). The reason the suspicion looks valid is that
the docstring at `R/core/09_clustering.R:7` says the wrong thing
(see C2 above).

**Action:** fix the docstring, not the YAML.

### F2. Proteomics `sample_filter: {enabled: false, rules: {}}` — never wired up
**The premise is wrong.** Proteomics YAML doesn't define a
`sample_filter` block at all. The `{enabled: false, rules: {}}` lives
in `metabolomics_template.yaml:133–138`. And it *is* wired up:

- Metabolomics: `R/modules/metabolomics/00_mod_preprocessing.R:124–134` and `398–410`.
- Proteomics (the omic that doesn't declare it): `R/domain/proteomics/04_preprocess.R:49–63` reads it via `get_sample_filter_rules(config, mode = "proteomics")` — when missing, returns `NULL`, the block is a no-op.
- RNA: `R/domain/rnaseq/03_preprocess.R:133–150`.

So `sample_filter` is alive in all three pipelines via
`R/core/03_alignment.R:120` (`apply_sample_filter`). The proteomics
YAML could declare an empty `sample_filter:` block for symmetry —
that's an editorial choice (open question j-7).

### F3. RNA `effects.heatmap_annotations` indentation
**Confirmed.** `config/templates/rna_config.yaml:159–165`:

```yaml
    effects:
      # Color variable for plots
      ...
      samples: "SampleID"

      # Heatmap annotation columns (optional)
      # If not specified, uses color variable only
      # heatmap_annotations:
      #   - "line"
      #   - "time_point"
```

The commented block is correctly indented under `effects:`, so when
uncommented it lands at `effects.heatmap_annotations`. R reads exactly
that path:

- `R/modules/rnaseq/01_mod_qc_pre.R:181` — `annot_cols <- cfg$effects$heatmap_annotations %||% NULL`
- `R/modules/rnaseq/03_mod_qc_post.R:164` — same

**It IS wired correctly.** The suspicion phrasing ("indentation looks
unintentional") doesn't reflect the file — the YAML is fine; it's just
commented out.

### F4. NEW: `filtering.auto_filter.{min,max,fallback}` (RNA) is plumbed only to the report — not to the function that actually filters
The RNA YAML declares (`rna_config.yaml:64–67`):

```yaml
auto_filter:
  min: 0.5
  max: 2.0
  fallback: 1.0
```

The actual call site at `R/domain/rnaseq/03_preprocess.R:190` is:

```r
fr <- run_auto_filter_pipeline(norm_for_filter, meta2, sample_col, group_col, output_plot = plot_path)
```

— **no `min_limit`, `max_limit`, or `fallback` are passed.**
`run_auto_filter_pipeline()` (`R/domain/rnaseq/02_filtering.R:126–130`)
uses its hard-coded function defaults (0.5, 2.0, 1.0).

Today this happens to match the YAML, so the bug is invisible. Any
user who edits `auto_filter.min` will see the report claim it was
applied (`report_template.Rmd:294,419,459,3005`) while the actual
threshold search ignores it.

**Action:** either delete the three YAML keys, or fix
`03_preprocess.R:190` to pass them through. The latter is correct;
proposed change:

```r
fr <- run_auto_filter_pipeline(
  norm_for_filter, meta2, sample_col, group_col,
  output_plot = plot_path,
  min_limit   = cfg$filtering$auto_filter$min      %||% 0.5,
  max_limit   = cfg$filtering$auto_filter$max      %||% 2.0,
  fallback    = cfg$filtering$auto_filter$fallback %||% 1.0
)
```

### F5. NEW: Proteomics validator enforces `id_columns.sample_col == id_columns.map_to`
`R/domain/proteomics/90_config_validate.R:21–23`:

```r
if (!identical(cfg$id_columns$sample_col, cfg$id_columns$map_to)) {
    stop("id_columns$sample_col and id_columns$map_to must be identical if map_to is provided")
}
```

This is a tautology under the current single-source design (any user
who follows the template will set both to `"SampleID"`). If we collapse
those keys (suspicion 4 / verdict M-series), this assertion goes away
naturally. Until then, it's a noisy guard for a self-imposed
duplicate.

### F6. NEW: Metabolomics `normalization.*` is largely a frozen artifact
Seven keys (`input_already_normalized`, `sample_norm`, `transform`,
`scaling`, `na_policy`, `is_ref_col`, `biological_factor_col`) are
either dead or read only by the validator. The active normalization
choice is `preprocessing.chosen_norm`. The block self-documents this
with comments at lines 95–101, 110–113, 119–120 — but the keys are
still present in fresh templates, which is misleading.

**Action:** remove all seven from the template, OR move them into a
commented-out "deprecated" block at the bottom of the metabolomics
config.

### F7. NEW: Proteomics validator declares `clustering.de_source` choices = `c("any_contrast")` — single-value enum
`R/domain/proteomics/90_config_validate.R:196`:

```r
assert_one_of(clust_cfg$de_source, "modes$proteomics$clustering$de_source", c("any_contrast"), allow_null = TRUE)
```

A one-element allowed set is effectively a constant. Either widen the
enum (the comment in the YAML at line 142 promises that future values
will arrive) or drop the key and hard-code the behaviour.

---

## (g) Structural drift — section X exists in omic Y only

| Section | prot | rna | met | multi | Recommendation |
|---------|:----:|:---:|:---:|:-----:|----------------|
| `qc` | ✓ (umap/heatmap-sizes/pca_subsets) | ✓ (adaptive/thresholds/heatmap_viz) | ✓ (qc_flag_column/missingness/thresholds) | — | **Different sub-schemas, same name.** Unify the THRESHOLD sub-schema (`qc.thresholds.*` appears in rna+met but not prot — proteomics would benefit). Leave the omic-specific keys (`run_umap`, `pca_subsets` for prot; `qc_flag_column` for met). |
| `de_table` | ✓ (`id_col`, `pass_any_col`) | ✓ (`id_col`, `pass_any_col`) | — | — | Keep as-is. Both keys are dead in the YAML (the consumer hard-codes them) — see (d). The section can probably be removed entirely. |
| `sample_filter` | — | (commented) | ✓ (`enabled`, `rules.*`) | — | The R code supports it for all three single-omic modes. Make the section present-but-disabled in every template — the consumer pattern (`%||% NULL`) handles missing fine, but explicit declarations are kinder to users. **Unify.** |
| `qc_post` | ✓ | ✓ | — | — | Metabolomics has no QC-post step (DE output structure differs). **Keep omic-specific.** |
| `imputation` | ✓ (perseus / dep2 / qrilc) | — | — | — | Proteomics-specific. **Keep.** |
| `domain_analysis`, `ppi`, `advanced_stats` | ✓ | — | — | — | Proteomics-specific feature areas. **Keep.** |
| `pathway` | ✓ (fgsea + ORA + databases + clustering) | ✓ (fgsea + ORA + databases) | — (uses `enrichment` instead) | — | RNA + prot share the section; metabolomics has its own `enrichment` block with different structure. **Naming drift:** unify name as `pathway` everywhere (metabolomics is the outlier), or as `enrichment` everywhere. See (j-1). |
| `enrichment` | — | — | ✓ (libraries / gmt / mapping / pvalue_column) | ✓ (run_enrichment) | Two different `enrichment` shapes between metabolomics and multiomics — also drift. See (j-1). |
| `commentary` | ✓ | ✓ | — (lives inside `7_mod_commentary.R` but no YAML block) | ✓ | The metabolomics commentary module reads from `metab_cfg$commentary` (via `%||%`) but the YAML doesn't declare it. **Add commentary block to metabolomics template.** |
| `excel` | ✓ | ✓ | ✓ | — | All three are present and consistent. **Keep.** |
| `annotation` | ✓ (`organism`, `id_type`, `skip_annotation`, `custom_mapping_file`) | ✓ (`source` only) | — | — (`global.organism` instead) | **Drift in schema.** RNA's `annotation.source` is a parser hint; prot's `annotation.*` is for ID conversion + pathway organism. They are not the same concept. Suggest renaming RNA's to `annotation.parser_source` to disambiguate. |
| `report` | ✓ (proteomics — 19 keys) | — | — | — | The proteomics report is the most mature; the RNA & metabolomics reports read configuration from elsewhere (`technical_report.*`, `effects.*`). **Pull a smaller shared subset** (`introduction_text`, `experimental_design_file`, `raw_files_location`, `analysis_files_location`, `full_results_table_link`, `disclaimer_text`) into a unified shared block. |
| `technical_report` | ✓ (facility / sample_prep / ms_acquisition / search_engine / search_parameters / acknowledgment) | ✓ (facility / library_prep / sequencing / acknowledgment) | — | — | Same idea, omic-specific fields. **Standardise the parent key (`technical_report`) and let omic-specific subkeys vary.** |
| `params` | ✓ (`seed`) | ✓ | ✓ | ✓ | Consistent. **Keep.** |
| `design`, `global`, `output` | — | — | — | ✓ | Multiomics-only top-level blocks. Some of their keys are dead (see (d)). **Trim and keep.** |
| `sample_alignment` | — | — | — | — | The user's task description mentions this — I found no such section in any of the four YAMLs. May have been removed already; or the user is thinking of `sample_filter`/`sample_map`. |

---

## (h) Documentation gaps in the YAMLs

- **`config/templates/rna_config.yaml:39–41`**: section says "Column
  mappings are currently hardcoded per source type. Future: make these
  config-driven." This is a TODO inline. Should be moved to an issue
  or to docs/, not left as a config comment.
- **`config/templates/rna_config.yaml:69–73`**: a commented-out
  `sample_filter:` block uses a completely different schema
  (`tissue: ["leaf", "root"]`) than the metabolomics template's
  schema (`enabled: false; rules: {exclude_blanks: …}`). The two
  formats are incompatible — one of them does not match what
  `apply_sample_filter()` expects (`R/core/03_alignment.R:120–145`).
  Pick one (the metabolomics schema is the one R actually reads).
- **`config/templates/rna_config.yaml:81`**: trailing whitespace
  after `de:` (no comment text on the next line — visible only as a
  blank line after the method field). Cosmetic.
- **`config/templates/rna_config.yaml:161–165`**: 4-line commented
  example for `heatmap_annotations` should be promoted to a real
  documented key with a sensible default (or removed if not promoted),
  not left as a commented block.
- **`config/templates/proteins_config.yaml:67–71`**: commented-out
  per-condition `min_count` override has no equivalent in
  `R/domain/proteomics/02_filtering.R` — verify the implementation
  supports the documented form (per-condition map under `min_count`)
  before leaving the example.
- **`config/templates/proteins_config.yaml:243`**: trailing line
  `# - name: "Custom subset"` with no following key — looks abandoned.
- **`config/templates/metabolomics_template.yaml:34–38`**:
  `# data_dir: "metabolics/dir_name"` and `# sample_map:` use the
  typo `metabolics/` (should be `metabolomics/`). Same typo on
  line 36, 38, 40. Consistent within metabolomics paths.
- **`config/templates/metabolomics_template.yaml:86–129`**: the entire
  `normalization:` block has long-form explanatory comments which read
  more like a design doc than an operational config. Either move the
  prose to docs/ or shorten the comments to one-liners.
- **`config/templates/multiomics_config.yaml:124–138`**: a 14-line
  commented-out metabolomics block uses an entirely different schema
  than the canonical metabolomics template (uses `input.has_group_row`,
  `normalization.sample_norm: "sum"`, `de.method: "t_test"`,
  `enrichment.run_enrichment` — none of which appear in
  `metabolomics_template.yaml`). This is stale documentation.
- **Comment style is inconsistent across the four YAMLs.**
  - `proteins_config.yaml` uses `============` and `------------`
    banners. `rna_config.yaml` mostly uses `------------`.
    `metabolomics_template.yaml` uses `==============================================================================`.
    `multiomics_config.yaml` uses a mix.
  - Trailing-comment style varies: `# REQUIRED:` vs `# REQUIRED` vs
    `# (REQUIRED)`. Same for inline guidance: `# "limma" | "ttest" | "welch" | "anova"` vs
    `# "deseq2" | "limma_voom" (example; depends what you implement)`.
- **No top-of-file metadata schema version.** None of the four
  templates declares a `schema_version:` field, so a future migration
  has nothing to gate on.

---

## (i) Proposed unified schema

### Canonical section ordering (single-omic templates)

```
project:           # identity & I/O roots
paths:             # raw / out
params:            # seed, etc.
modes:
  <omic>:
    engine:        # data-origin / engine (proteomics, metabolomics input)
    input:         # format selection (metabolomics; could be empty for prot/rna)
    files:         # all input files
    id_columns:    # column-name conventions (REQUIRED columns only)
    parsing:       # parser hints (currently metabolomics-only)
    annotation:    # organism / id_type / skip / custom mapping
    sample_filter: # sample-level filter rules (uniform schema)
    filtering:     # feature-level filtering (per omic)
    normalization: # method choice + parameters
    imputation:    # proteomics-only; absent in others
    preprocessing: # metabolomics-only (chosen_norm + drift)
    de:            # differential expression
    de_table:      # (optional — drop if (d/de_table) confirmed dead)
    qc_post:       # post-DE QC (prot/rna only)
    clustering:    # downstream clustering (3 steps + heatmap)
    effects:       # plotting aesthetics (color / shape / samples / heatmap_annotations)
    excel:         # excel export
    qc:            # QC-specific advanced settings + thresholds + heatmap_viz
    pathway:       # enrichment (unify name across omics — see (j-1))
    ppi:           # proteomics-only
    domain_analysis: # proteomics-only
    advanced_stats: # proteomics-only
    technical_report: # facility + omic-specific subkeys
    commentary:    # AI commentary
    report:        # shared subset (intro + locations + links + show flags)
```

### Naming convention — recommend snake_case everywhere

All keys are already snake_case **except** `normalization.prior.count`
in RNA. Rename to `prior_count`. The dotted form looks nested but isn't;
fresh contributors expect `prior` to be a section and `count` to be a
leaf.

**Rule:** "applies to all omics".

### Comment style

- **Section headers:** one banner line, e.g. `# --- DE settings ---`,
  or fenced `# ===== Section =====` — pick one (recommend the simple
  triple-dash `# --- ... ---` to match the dominant pattern in
  `rna_config.yaml`).
- **Inline comments on a value:** `key: value  # one-line guidance` —
  exactly two spaces before `#`, one space after `#`.
- **Multi-line guidance:** put it above the section header, not above
  each key.
- **Banned:** trailing bare `#` with no text; commented-out blocks
  whose syntax doesn't match a current key; multi-paragraph
  explanations inside the YAML (move to docs/).

**Rule:** "applies to all omics".

### Where defaults should live

Currently mixed:

- Single-file defaults: validators in
  `R/domain/multiomics/90_config_validate.R` apply ~17 defaults to the
  multiomics block; the other three validators apply none.
- Call-site defaults: ~200 `%||%` sites apply leaf defaults at the
  function that reads the key.

This is the single biggest source of "silent NULL" risk: a key can be
missing from the YAML, validate cleanly, and silently use a
function-default that the YAML never documented.

**Recommendation:**

1. **Move all simple scalar defaults into the validators.** Each
   `validate_<omic>_config()` should `config <- modifyList(defaults, config)` before validating, exactly like
   `validate_multiomics_config()` does today (`R/domain/multiomics/90_config_validate.R:99–151`).
2. **Keep call-site `%||%` for derived/conditional defaults** (e.g.
   `effects$color %||% "Condition"` is a real fallback chain that
   would not collapse to a single literal).
3. **Single shared defaults file vs per-omic:** start per-omic
   (kept next to the validators that read them). Cross-omic defaults
   that genuinely are shared (e.g. `params.seed`, clustering structure)
   can be hoisted to `R/core/04_config.R` later if duplication appears.

**Rule:** "applies to all omics" — but per-omic exceptions are allowed
for any key only present in that omic.

---

## (j) Open questions for Michal

These all need decisions before any code touches the YAMLs or R sources.
Numbered for inline answers.

1. **Canonical name for the "pathway / enrichment" section.** Right now
   prot + rna use `pathway:`; metabolomics + multiomics use
   `enrichment:`. Which name wins? My recommendation: `pathway` in
   prot/rna (gene-centric, KEGG/GO/Reactome) and `enrichment` in
   metabolomics (compound-centric, mummichog/KEGG-compounds) — they
   genuinely are different concepts. Multiomics should use
   `enrichment` (it merges across).

2. **Canonical name for "sample column".** Two viable options:
   (a) keep `effects.samples` everywhere (de facto winner today, ~50
   reads vs ~14); or (b) rename `effects.samples` →
   `id_columns.sample_col` and drop `effects.samples`. (b) is more
   semantically correct (a sample-id column is an identifier, not an
   aesthetic). (a) is the smaller refactor.

3. **Authorise the SampleID collapse (suspicion 4 verdict).** Verdict
   was SAFE AFTER R REFACTOR with ≈12 edits across ~7 files. Do you
   want me to plan that PR separately, or fold it into the rename in
   (j-2)?

4. **For the four suspicion verdicts:**
   - 1 (`clustering.min_groups`): KEEP. Confirm?
   - 2 (`effects` section): KEEP, but pick option (a) or (b) from (j-2). Confirm?
   - 3 (`partition.x_axis_col` / `color_col`): SAFE TO REMOVE. Confirm?
   - 4 (SampleID triple): SAFE AFTER R REFACTOR — see (j-3).

5. **Cross-omic conflict resolutions:**
   - C1 — canonical = `counts_cutoff_high` / `counts_cutoff_low`. Confirm?
   - C2 — fix the stale docstring at `R/core/09_clustering.R:7` + drop
     the templated `binary_patterns.group_col` write in `run.R:676`. Confirm?
   - C3 — rename `prior.count` → `prior_count`. Confirm?
   - C4 — rename multiomics-template `gene_id_col` → `gene_id`. Confirm?

6. **Dead-parameter purge — which of these should be removed from the
   YAML templates?** All are confirmed unread by R. Tick the ones you
   want gone:
   - [ ] `report.project_id` (prot)
   - [ ] `de_table.pass_any_col` (prot, rna)
   - [ ] `de_table.id_col` (prot, rna) — wait, this one IS read. Keep.
   - [ ] `design.batch_column` (multi)
   - [ ] `global.metadata` (multi)
   - [ ] `output.output_dir` (multi)
   - [ ] `output.report_format` (multi)
   - [ ] `normalization.input_already_normalized` (met)
   - [ ] `normalization.sample_norm` (met — legacy)
   - [ ] `normalization.transform` (met — legacy)
   - [ ] `normalization.scaling` (met — legacy)
   - [ ] `normalization.na_policy` (met)
   - [ ] `normalization.is_ref_col` (met)
   - [ ] `normalization.biological_factor_col` (met)
   - [ ] `de.pval_cutoff` (met — duplicate of `p_cutoff`)
   - [ ] `de.robust_ebayes` (met)
   - [ ] `qc.missingness_heatmap` (met)

7. **Missing-section symmetry:** add empty (disabled) `sample_filter:`
   to prot + rna templates? Add `commentary:` to metabolomics template?
   Add `qc.thresholds.*` to proteomics template?

8. **`filtering.auto_filter.{min,max,fallback}` (RNA):** fix the
   plumbing at `R/domain/rnaseq/03_preprocess.R:190` (preferred) or
   remove the YAML keys?

9. **`clustering.de_source` enum is `c("any_contrast")` only.** Drop
   the key and hard-code, or widen the enum (to what — `"top_n"`,
   `"manual"`, `"pass_strict"`?)?

10. **Defaults — move scalars into validators?** This is the most
    impactful change in (i). Confirm direction before I plan it.

11. **Schema version:** add `schema_version: 1` to every template?

12. **YAML structural conventions:** I'm proposing a single canonical
    section order (see (i)). The current order is omic-specific. Is it
    acceptable to reorder the existing templates to match the canonical
    order, even though git-blame will show a large reformatting commit?

---

## Method notes (for reviewers)

- All key-by-key reads were verified by `grep -rn '<leaf_name>' R/`
  with at least two grep variants per key (e.g. `cfg$X$leaf`,
  `[["leaf"]]`, alias forms like `report_cfg$leaf`).
- `*.Rmd` files (`R/domain/proteomics/report_template_proteomics.Rmd`,
  `R/domain/rnaseq/report_template.Rmd`,
  `R/pipeline/metabolomics/templates/report_metabolomics.Rmd`,
  `R/domain/multiomics/report_template_multiomics.Rmd`) were
  searched in the same pass — many proteomics `report.*` keys live
  only there, behind a `report_cfg <- prot_cfg$report` alias.
- Validators were read in full
  (`R/core/04_config.R`,
  `R/domain/proteomics/90_config_validate.R`,
  `R/domain/rnaseq/90_config_validate.R`,
  `R/domain/metabolomics/00_inputs.R:127–177`,
  `R/domain/multiomics/90_config_validate.R`).
- No YAML or R file was edited. The only file changed by this audit
  is this one (`docs/configs-audit.md`).
