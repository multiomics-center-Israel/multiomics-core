# DE → DA Rename Impact Map

*Date: 2026-06-16*

> **Read-only impact map for a *potential* DE (Differential Expression) → DA (Differential Abundance) rename in the metabolomics mode. Nothing has been renamed — this is a survey only.**
>
> **The high-risk items are the Shiny export / contract files** (`R/core/07_shiny_contract.R`, `R/domain/metabolomics/07_shiny_export.R`, and the parallel RNA/proteomics exporters): the payload keys are a shared cross-omics contract that the Shiny app depends on by exact name. See §5.

---

COMPREHENSIVE IMPACT MAP: DE → DA RENAME
Based on my thorough read-only investigation of the repository, here is the complete map of DE (Differential Expression) terminology that would be affected by a rename to DA (Differential Abundance) for metabolomics:
 
1. CONFIG KEY de: AND NESTED KEYS
Location and reading:
•	Config files (multiple instances):
o	/home/user/multiomics-core/config/config_GT4.yaml:65-73
o	/home/user/multiomics-core/config/config_GT15.yaml:65-73
o	/home/user/multiomics-core/config/templates/metabolomics_template.yaml:186-194
o	/home/user/multiomics-core/config/templates/rna_config.yaml:79-105 (RNA, not metabolomics)
o	/home/user/multiomics-core/config/templates/proteins_config.yaml:104-131 (Proteomics, not metabolomics)
o	/home/user/multiomics-core/config/templates/multiomics_config.yaml:44-51 (RNA section) and 94-101 (Proteomics section)
o	/home/user/multiomics-core/run.R:626, 1066, 1677 (config generation in wizard)
Nested keys under de: (metabolomics context):
•	method (e.g., "limma", "t_test", "wilcoxon")
•	condition_column (e.g., "condition")
•	contrasts (list of contrast specifications)
•	p_cutoff (e.g., 0.05)
•	pval_cutoff (duplicate of p_cutoff)
•	linear_fc_cutoff (e.g., 1.5)
•	skip / enabled (boolean to disable DE)
Where config de: is READ in code:
•	/home/user/multiomics-core/R/modules/metabolomics/02_mod_de.R:22 → de_cfg <- cfg$de %||% list()
•	/home/user/multiomics-core/R/domain/metabolomics/03_differential.R:190 → de_cfg <- cfg$de %||% list()
•	/home/user/multiomics-core/R/modules/metabolomics/04_mod_clustering.R:24 → (cfg$de %||% list())$condition_column
•	/home/user/multiomics-core/R/modules/metabolomics/04b_mod_network.R:47 → p_cutoff <- metab_cfg$de$p_cutoff
•	/home/user/multiomics-core/R/modules/metabolomics/07_mod_commentary.R:220, 224, 251-254 (multiple references to metab_cfg$de)
•	/home/user/multiomics-core/R/domain/metabolomics/07_shiny_export.R:60 → de_cfg <- metab_cfg$de %||% list()
Risk: The config key name de: is hard-coded in all these locations. A rename to da: would require changes in every file that accesses cfg$de or config$modes$metabolomics$de.
 
2. TARGET NAMES AND FUNCTION NAMES (Metabolomics DE sense)
Target names (in pipeline files):
•	/home/user/multiomics-core/R/pipeline/metabolomics/00_pipe_metabolomics.R:315 → metab_de_res (target name)
•	/home/user/multiomics-core/R/pipeline/proteomics/00_pipe_proteomics.R:56 → prot_de_res (target name)
•	/home/user/multiomics-core/R/pipeline/rnaseq/00_pipe_rnaseq.R:10 → rna_de_res (target name)
Module functions:
•	/home/user/multiomics-core/R/modules/metabolomics/02_mod_de.R:17 → mod_metabolomics_de() function
•	/home/user/multiomics-core/R/modules/proteomics/02_mod_de.R:32 → mod_proteomics_de() function
•	/home/user/multiomics-core/R/modules/rnaseq/02_mod_de.R:3 → mod_rnaseq_de() function
Domain-specific DE computation functions:
•	/home/user/multiomics-core/R/domain/metabolomics/03_differential.R:188 → run_metabolomics_de() function
•	/home/user/multiomics-core/R/domain/rnaseq/04_de_summary.R:16 → run_deseq2_de() function
•	/home/user/multiomics-core/R/domain/rnaseq/04_de_summary.R:291 → run_limma_rna_de() function
•	/home/user/multiomics-core/R/domain/proteomics/05_de_summary.R:306 → run_ttest_de() function
DE output/loader functions:
•	/home/user/multiomics-core/R/domain/metabolomics/03_differential.R:80 → load_precomputed_metabolomics_de() function
•	/home/user/multiomics-core/R/domain/rnaseq/04_de_summary.R:375 → Pre-computed loader
Internal DE helpers and plotting:
•	/home/user/multiomics-core/R/modules/metabolomics/02_mod_de.R:166 → plot_pvalue_histogram() (for DE p-value distribution)
•	/home/user/multiomics-core/R/core/06_plots.R (contains generic plot_volcano, plot_ma for DE plots)
Contract validation:
•	/home/user/multiomics-core/R/core/02_validation.R:162 → assert_de_contract() function (validates DE result structure)
•	/home/user/multiomics-core/R/modules/metabolomics/02_mod_de.R:46 → calls assert_de_contract(de_res, stage = stage)
Risk: All these function names are called throughout the pipeline and referenced in reports, Shiny code, and downstream analyses. They are NOT specific to metabolomics (e.g., mod_rnaseq_de exists for RNA), so renaming metabolomics-specific ones would create inconsistency across omics.
 
3. DATA STRUCTURE FIELD NAMES AND COLUMN NAMES
Internal DE result objects:
•	de_res$summary_df — wide format table with one row per feature, columns per contrast (called in all downstream modules)
•	de_res$de_tables — list of per-contrast data.frames (e.g., metabolomics uses this in /home/user/multiomics-core/R/modules/metabolomics/04_mod_clustering.R:69)
•	de_res$method — stores the DE method used (e.g., "limma")
•	de_res$de_model — fitted model object (limma fit)
Column names in DE output tables (metabolomics):
•	/home/user/multiomics-core/R/domain/metabolomics/03_differential.R:277-287:
o	feature_id
o	Contrast
o	logFC (log fold-change)
o	AveExpr (average expression)
o	t (t-statistic)
o	P.Value (raw p-value)
o	adj.P.Val (adjusted p-value)
o	B (B-statistic from eBayes)
Output file names (Excel, TSV):
•	/home/user/multiomics-core/R/core/05_export_excel.R:79-81:
o	Final_results_ALL_P_*.xlsx (all features)
o	Final_results_DE_P_*.xlsx (DE-significant features only)
•	/home/user/multiomics-core/R/modules/metabolomics/02_mod_de.R:68 → de_summary.tsv
•	/home/user/multiomics-core/R/modules/metabolomics/02_mod_de.R:86 → de_<contrast_label>.tsv
Shiny payload keys (standardized across omics):
•	/home/user/multiomics-core/R/core/07_shiny_contract.R:108-147:
o	de_model
o	de_stats (full DE statistics for all features/contrasts)
o	de_sig_stats (subset where significant)
o	de_expr_norm (expression matrix of DE features only)
o	de_summary (per-contrast summary counts)
o	de_final_table
o	all_final_xlsx
o	de_final_xlsx
Risk: These column names and output file names are hardcoded and depended upon by:
•	Excel export logic (users may reference these in downstream analysis)
•	Shiny app (payload keys must match hardcoded contract keys)
•	Reports (column headers visible to users)
•	Downstream integration analyses (multiomics concordance, enrichment use de_tables)
 
4. USER-FACING STRINGS (Reports, Plots, Messages)
In metabolomics report template:
•	/home/user/multiomics-core/R/pipeline/metabolomics/templates/report_metabolomics.Rmd:717 → "differential expression testing"
•	/home/user/multiomics-core/R/pipeline/metabolomics/templates/report_metabolomics.Rmd:754 → "differential expression testing"
•	/home/user/multiomics-core/R/pipeline/metabolomics/templates/report_metabolomics.Rmd:1586 → "differentially expressed metabolites"
•	/home/user/multiomics-core/R/pipeline/metabolomics/templates/report_metabolomics.Rmd:1685, 1731, 1785, 1794 → "differentially expressed metabolites" (multiple occurrences)
•	/home/user/multiomics-core/R/pipeline/metabolomics/templates/report_metabolomics.Rmd:2237 → "differentially expressed metabolites"
•	/home/user/multiomics-core/R/pipeline/metabolomics/templates/report_metabolomics.Rmd:2374 → "Full differential expression summary table"
•	/home/user/multiomics-core/R/pipeline/metabolomics/templates/report_metabolomics.Rmd:2477-2521 → Multiple references to "differentially expressed" in results narrative
Plot titles and legends (metabolomics):
•	/home/user/multiomics-core/R/modules/metabolomics/02_mod_de.R:94-95 → "Volcano: <contrast> (adj.P.Val)" / "Volcano: <contrast> (P.Value)"
•	/home/user/multiomics-core/R/modules/metabolomics/02_mod_de.R:124 → "MA: <contrast>"
•	/home/user/multiomics-core/R/modules/metabolomics/02_mod_de.R:139 → "P-value: <contrast>"
Diagnostic plot directories:
•	/home/user/multiomics-core/R/modules/metabolomics/02_mod_de.R:31 → Output dir DE_plots (hardcoded)
Console messages:
•	/home/user/multiomics-core/R/modules/metabolomics/02_mod_de.R:26 → "metabolomics DE: disabled in config — skipping"
•	/home/user/multiomics-core/R/domain/metabolomics/03_differential.R:275 → "metabolomics DE [limma]: <contrast>"
•	/home/user/multiomics-core/R/domain/metabolomics/03_differential.R:298 → "metabolomics DE [<method>]: <contrast>"
•	/home/user/multiomics-core/R/domain/metabolomics/03_differential.R:313-314 → Summary message on DE completion
Risk: All user-visible text would change. Reports generated with "differential expression" language would be inconsistent with "differential abundance" terminology. This is acceptable for metabolomics specifically but creates terminology mismatch with RNA-seq/proteomics outputs.
 
5. SHINY APP & PAYLOAD DEPENDENCIES
Shiny contract (canonical v2.0):
•	/home/user/multiomics-core/R/core/07_shiny_contract.R:108-147 defines canonical payload keys with exact names:
o	de_model, de_stats, de_sig_stats, de_expr_norm, de_summary, de_final_table, all_final_xlsx, de_final_xlsx
Shiny payload builder (metabolomics):
•	/home/user/multiomics-core/R/domain/metabolomics/07_shiny_export.R:35-49 → build_shiny_payload_metabolomics() function receives de_res parameter
•	/home/user/multiomics-core/R/domain/metabolomics/07_shiny_export.R:224 → Comment references "Final_results_DE_P_*.xlsx"
•	Lines 60-90: Reads metab_cfg$de, passes it to payload building logic
Payload validation (all omics):
•	/home/user/multiomics-core/R/core/07_shiny_contract.R:446 → assert_de_stats() validates presence of required columns in de_stats
Shiny export for RNA and proteomics (identical contract):
•	/home/user/multiomics-core/R/domain/rnaseq/06_shiny_export.R:21 → receives de_res parameter
•	/home/user/multiomics-core/R/domain/proteomics/07_shiny_export.R (similar structure)
Risk (CRITICAL FOR SHINY):
•	The canonical payload keys are SHARED across all omics (RNA, proteomics, metabolomics)
•	The Shiny app code (unknown location but referenced in CLAUDE.md) depends on these exact key names
•	Renaming metabolomics-specific keys to da_* would break the contract for metabolomics payloads
•	Alternative: only rename metabolomics column/file naming while keeping payload keys consistent (not fully compliant with DA terminology but safer)
