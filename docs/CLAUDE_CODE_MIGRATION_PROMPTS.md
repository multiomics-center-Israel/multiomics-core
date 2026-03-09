# Claude Code Migration Prompts — 3 Agents

**Project:** Migrate multiomics_pipeline into multiomics-core  
**Date:** 2026-02-28

Split across 3 Claude Code agents. Each agent works independently on non-overlapping files.

---

## How to Run

**⚠️ IMPORTANT: Agents must run SEQUENTIALLY, not in parallel.**

Agent 3 depends on functions created by Agents 1 and 2. Running them simultaneously will cause syntax check failures.

**Execution Order:**
1. **Agent 1 FIRST** — Creates MOFA2 integration, enrichment plots, config updates
2. **Agent 2 SECOND** — Creates foundational/mechanistic analysis (after Agent 1 completes)
3. **Agent 3 LAST** — Wires pipeline targets that call functions from Agents 1 & 2

Open Claude Code in `/home/ozsol/multiomics-core/` and run each agent's prompt in sequence.

---

## AGENT 1 — MOFA2 Integration + Enrichment/Pathview + Config Updates

```
You are migrating code from a standalone multi-omics pipeline into the multiomics-core R project.

SOURCE REPO: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/
DESTINATION REPO: /home/ozsol/multiomics-core/

=== CONTEXT: multiomics-core Architecture ===

multiomics-core uses a 3-layer architecture:
  R/core/       → shared utilities (config, io, validation) — DO NOT MODIFY THESE
  R/domain/     → pure analytical functions (stateless, per-omics subdirs)
  R/modules/    → orchestration wrappers (combine domain functions + I/O)
  R/pipeline/   → {targets} DAG definitions
  R/services/   → external services (AI commentary)

Config is accessed via config$modes$multiomics$* (NOT config$integration$* directly).
The top-level _targets.R auto-sources all R/ files via tar_source(), so NO source() calls needed.
Logging uses message() (NOT log_message() from the source).
The %||% null-coalescing operator is already available from R/core/04_config.R.
Output directories are passed as explicit out_dir parameters (NOT config$output$output_dir).
Data is stored as formal MultiAssayExperiment objects (NOT custom lists).

=== EXISTING multiomics-core FILES (for reference, read but don't rewrite unless told) ===

- R/domain/multiomics/04_integration_diablo.R  — DIABLO integration (already exists, pattern to follow)
- R/domain/multiomics/05_integration_snf.R     — SNF integration (already exists, pattern to follow)
- R/domain/multiomics/07_enrichment.R          — Cross-omics enrichment (Fisher meta-analysis)
- R/modules/multiomics/02_mod_integration.R    — Integration module (ADD MOFA2 block here)
- R/pipeline/multiomics/00_pipe_multiomics.R   — Pipeline DAG (ADD new targets here)
- R/domain/multiomics/90_config_validate.R     — Config validation (ADD MOFA2 validation here)
- config/multiomics_GT15_test.yaml             — Test config (ADD new sections here)

=== YOUR TASKS ===

TASK 1: Create R/domain/multiomics/04b_integration_mofa.R
  Source: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/R/06_mofa.R (~419 lines)
         /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/R/mofa_wrapper.R (~316 lines)

  Read these source files carefully. Port them into a SINGLE file following multiomics-core conventions:

  Minimum changes required:
  a) Replace ALL config$integration$mofa2$* with config$modes$multiomics$integration$mofa2$*
  b) Replace log_message(...) with message(...)
  c) Replace save_plot(plot, filename, config, ...) with:
     ggsave(file.path(out_dir, filename), plot, width=w, height=h, dpi=300)
  d) Replace save_table(df, filename, config, ...) with:
     write.csv(df, file.path(out_dir, filename), row.names = FALSE)
  e) The main function signature should be: run_mofa2_integration(mae, config, out_dir = NULL)
     - Input `mae` is a formal MultiAssayExperiment. Extract matrices via:
       lapply(names(experiments(mae)), function(nm) as.matrix(assay(mae[[nm]]))) |> setNames(names(experiments(mae)))
     - This replaces the source's extract_matrices_for_integration(feature_data$filtered_mae)
  f) For the mofa_wrapper.R portion, keep the run_mofa() function that calls scripts/run_mofa.py
     but update the script path to: "scripts/run_mofa.py" (relative from project root).
     Change the default python_exec to: Sys.getenv("RETICULATE_PYTHON", "python3")
  g) Include write_mofa_results(mofa_results, out_dir) function to save CSVs.
     The function must write these outputs (following write_diablo_results() pattern):
       - mofa_factors.csv — sample factor values (samples × factors matrix)
       - mofa_weights_<view>.csv — feature weights per view (features × factors)
       - mofa_variance_explained.csv — variance explained per factor per view
       - mofa_top_features.csv — top contributing features per factor
  h) The run_mofa2_integration() function must return a named list with:
       - model: the MOFA model object (or NULL if using Python wrapper)
       - factors: data.frame of sample factor values
       - weights: named list of weight matrices per view
       - variance_explained: data.frame of variance per factor/view
       - top_features: data.frame of top features per factor
       - plots: list of generated ggplot objects
       - config: the mofa2 config used
  i) Follow the pattern of 04_integration_diablo.R for function documentation and structure.

TASK 2: Copy scripts/run_mofa.py
  Source: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/scripts/run_mofa.py
  Destination: /home/ozsol/multiomics-core/scripts/run_mofa.py
  Copy as-is with no changes.

TASK 3: Create R/domain/multiomics/07b_multigsea_plots.R
  Source: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/R/13_multigsea_plots.R (~401 lines)

  Port with these changes:
  a) Replace config$enrichment$multigsea$* with config$modes$multiomics$enrichment$multigsea$*
  b) Replace config$output$output_dir with explicit out_dir parameter
  c) Replace log_message() with message()
  d) Replace save_plot()/save_table() with ggsave()/write.csv()
  e) Main function: run_multigsea_plots(enrichment_results, config, out_dir = NULL)
  f) Also port run_multigsea_pathview() from same source file.

TASK 4: Create R/domain/multiomics/07c_kegg_pathview.R
  Source: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/R/15_kegg_pathview.R (~770 lines)

  Port with these changes:
  a) Remove install_load_packages() function — packages managed externally in multiomics-core
  b) Replace hard-coded data/HMDB2kegg_cpd.Jan2026.v2.txt with a parameter that defaults to
     file.path("data", "HMDB2kegg_cpd.Jan2026.v2.txt")
  c) Replace log_message() with message()
  d) Keep all analytical functions: map_symbols_to_entrez(), get_hmdb_kegg_mapping(),
     map_hmdb_to_kegg(), load_gene_data(), load_metabolite_data(), plot_kegg_overlay(),
     plot_multiple_pathways()
  e) Also port run_consensus_pathview() from 10_enrichment.R source if it references pathview functions.

TASK 5: Update R/modules/multiomics/02_mod_integration.R
  Add MOFA2 integration block. Read the existing file first, then add after the SNF block:

  ```r
  # MOFA2
  if ("MOFA2" %in% methods) {
      message("\nRunning MOFA2 integration...")
      mofa_dir <- file.path(out_dir, "mofa")
      dir.create(mofa_dir, showWarnings = FALSE)

      results$mofa_results <- tryCatch({
          run_mofa2_integration(
              mae = mae_subset,
              config = config,
              out_dir = mofa_dir
          )
      }, error = function(e) {
          warning("MOFA2 integration failed: ", e$message)
          NULL
      })

      if (!is.null(results$mofa_results)) {
          write_mofa_results(results$mofa_results, mofa_dir)
      }
  }
  ```

  Update the return list to include mofa_results (use consistent `_results` suffix):
  ```r
  list(
      diablo_results = results$diablo,
      snf_results = results$snf,
      mofa_results = results$mofa_results,  # <-- ADD THIS LINE
      feature_selection = feature_selection,
      mae_subset = mae_subset,
      methods_run = methods
  )
  ```

TASK 6: Update R/domain/multiomics/90_config_validate.R

  NOTE: The file already has `valid_methods <- c("DIABLO", "SNF", "MOFA2")` at line 34.
  MOFA2 is listed but not validated. You must ADD the validation block.

  Add MOFA2 config validation block (after the existing SNF block, around line 71):
  ```r
  # Validate MOFA2 config if specified
  if ("MOFA2" %in% methods) {
      mofa2_cfg <- multiomics_cfg$integration$mofa2 %||% list()

      if (is.null(mofa2_cfg$num_factors)) {
          message("  MOFA2 num_factors not specified, defaulting to 10")
          multiomics_cfg$integration$mofa2$num_factors <- 10
      } else if (mofa2_cfg$num_factors < 1 || mofa2_cfg$num_factors > 50) {
          warning("MOFA2 num_factors should be between 1 and 50; got ", mofa2_cfg$num_factors)
      }

      valid_conv_modes <- c("fast", "medium", "slow")
      if (!is.null(mofa2_cfg$convergence_mode) &&
          !(mofa2_cfg$convergence_mode %in% valid_conv_modes)) {
          stop("MOFA2 convergence_mode must be one of: ", paste(valid_conv_modes, collapse=", "))
      }

      if (is.null(mofa2_cfg$seed)) {
          multiomics_cfg$integration$mofa2$seed <- 42
      }
  }
  ```

  Also add validation for new config sections (foundational, consensus, stability, enrichment$multigsea)
  with sensible defaults following the existing pattern.

TASK 7: Update config/multiomics_GT15_test.yaml
  Read the existing file first. Add these sections under modes.multiomics.

  The current file has this structure under integration:
    integration:
      methods: [...]
      snf: {...}

  You need to ADD diablo and mofa2 sections. Full integration block should be:

  ```yaml
    integration:
      methods:
        - "DIABLO"
        - "SNF"
        - "MOFA2"

      diablo:
        ncomp: 2
        design: "full"      # or "null" for no design matrix

      snf:
        K: 15
        alpha: 0.5
        T: 20
        n_clusters: 2

      mofa2:
        num_factors: 10
        convergence_mode: "fast"
        seed: 42
  ```

  Add these NEW sections at the same level as integration (under modes.multiomics):

  ```yaml
    foundational:
      run_foundational: true
      correlation_method: "spearman"
      top_variable_features: 500

    consensus:
      compare_methods: true

    stability:
      run_stability: false
      n_bootstrap: 100

    mechanistic:
      run_mechanistic: false

    # Extend existing enrichment section (currently just run_enrichment: true)
    enrichment:
      run_enrichment: true
      multigsea:
        run_multigsea: true
        pvalue_threshold: 0.05
      pathview:
        run_pathview: true
        top_n: 5

    commentary:
      enabled: false
      backend: "none"
  ```

=== IMPORTANT RULES ===
- READ source files before porting. Don't invent code — adapt what exists.
- Keep ALL analytical logic intact. Only change infrastructure/plumbing.
- Use message() not log_message(). Use %||% freely (already available).
- Functions in R/domain/ must be STATELESS — no side effects except writing to out_dir.
- All new files must have roxygen-style (#') documentation headers.
- Do NOT modify R/core/ files.
- Do NOT modify R/domain/multiomics/04_integration_diablo.R or 05_integration_snf.R.
- After creating all files, run: Rscript -e 'source("R/domain/multiomics/04b_integration_mofa.R")' to check for syntax errors.

=== WHEN COMPLETE ===
After finishing all tasks and passing syntax checks, output this status update:

✅ AGENT 1 COMPLETE:
  - Created: R/domain/multiomics/04b_integration_mofa.R
  - Created: R/domain/multiomics/07b_multigsea_plots.R
  - Created: R/domain/multiomics/07c_kegg_pathview.R
  - Copied: scripts/run_mofa.py
  - Updated: R/modules/multiomics/02_mod_integration.R (added MOFA2)
  - Updated: R/domain/multiomics/90_config_validate.R
  - Updated: config/multiomics_GT15_test.yaml
  - Syntax checks: PASSED
→ Agent 2 can now start.
```

---

## AGENT 2 — Foundational Analysis + RNA-Protein Correlation + Mechanistic Inference

```
You are migrating code from a standalone multi-omics pipeline into the multiomics-core R project.

SOURCE REPO: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/
DESTINATION REPO: /home/ozsol/multiomics-core/

=== CONTEXT: multiomics-core Architecture ===

multiomics-core uses a 3-layer architecture:
  R/core/       → shared utilities (config, io, validation) — DO NOT MODIFY
  R/domain/     → pure analytical functions (stateless, per-omics subdirs)
  R/modules/    → orchestration wrappers (combine domain functions + I/O)
  R/pipeline/   → {targets} DAG definitions — DO NOT MODIFY (Agent 3 handles this)
  R/services/   → external services

Key conventions:
- Config accessed via config$modes$multiomics$* (NOT config$integration$* directly)
- The %||% operator is available from R/core/04_config.R
- Logging: use message() (NOT log_message())
- Output: pass out_dir as parameter (NOT config$output$output_dir)
- Data: formal MultiAssayExperiment (NOT custom mae_data list)
- All R/ files are auto-sourced by _targets.R — NO source() calls in files
- Functions must NOT have side effects beyond writing to out_dir

=== DATA STRUCTURE BRIDGE ===

The source pipeline passes mae_data as a custom list:
  mae_data$harmonized_omics$transcriptomics$normalized_matrix  → a matrix
  mae_data$harmonized_omics$transcriptomics$de_table           → a data.frame
  mae_data$metadata                                             → sample metadata
  mae_data$common_samples                                       → character vector
  mae_data$gene_mapping                                         → gene-protein mapping df

In multiomics-core, the MAE is a formal MultiAssayExperiment.

**OPTION A (Recommended for large files):** Use an adapter to minimize code changes.
Define this adapter at the top of EACH new domain file:

**OPTION B (Cleaner but more work):** Refactor source code to use direct MAE access:
  - names(mae@ExperimentList) instead of names(mae_data$harmonized_omics)
  - SummarizedExperiment::assay(mae[[om]]) instead of mae_data$harmonized_omics[[om]]$normalized_matrix
  - as.data.frame(colData(mae)) instead of mae_data$metadata
  See R/domain/multiomics/04_integration_diablo.R for examples.

For this migration, use OPTION A (adapter) to minimize code changes:

```r
#' Convert formal MAE + extras to legacy mae_data list
#' Used internally to bridge source pipeline data structures
.mae_to_legacy <- function(mae, de_results = NULL, gene_protein_mapping = NULL) {
    harmonized_omics <- lapply(names(experiments(mae)), function(nm) {
        exp_data <- experiments(mae)[[nm]]
        de <- if (!is.null(de_results[[nm]])) de_results[[nm]] else NULL
        list(
            normalized_matrix = as.matrix(assay(exp_data)),
            de_table = de,
            da_table = de,
            feature_annotation = as.data.frame(rowData(exp_data))
        )
    })
    names(harmonized_omics) <- names(experiments(mae))

    list(
        mae = mae,
        harmonized_omics = harmonized_omics,
        metadata = as.data.frame(colData(mae)),
        common_samples = colnames(mae),
        gene_mapping = gene_protein_mapping
    )
}
```

Then each function can call .mae_to_legacy(mae) internally and work with the familiar
mae_data$harmonized_omics$*$normalized_matrix pattern — MINIMAL code changes.

=== YOUR TASKS ===

TASK 1: Create R/domain/multiomics/06b_rna_protein_correlation.R
  Source: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/R/14_rna_protein_correlation.R (~343 lines)

  Port the rna_protein_correlation() function with these changes:
  a) New signature: run_rna_protein_correlation(mae, de_results = NULL, gene_protein_mapping = NULL, config, out_dir = NULL)
  b) Use .mae_to_legacy() internally to convert MAE to legacy format
  c) Replace config$output$output_dir with out_dir parameter
  d) Replace log_message() with message()
  e) Keep ALL analytical logic: expression correlation, log2FC concordance, translation efficiency
  f) Keep ALL plotting code (histogram, scatter, etc.) — just change save path

TASK 2: Create R/domain/multiomics/09_foundational_correlations.R
  Source: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/R/05b_foundational_correlations.R (~1,880 lines)

  This is the LARGEST file. Port it with these changes:
  a) New main function signature: run_foundational_analysis(mae, config, gene_protein_mapping = NULL, out_dir = NULL)
  b) Use .mae_to_legacy() to convert MAE, then pass mae_data internally as before
  c) Replace ALL config$foundational$* with config$modes$multiomics$foundational$*
  d) Replace log_message() with message()
  e) Replace save_plot(plot, filename, config, ...) with ggsave(file.path(out_dir, "plots", filename), plot, ...)
  f) Replace save_table(df, filename, config, ...) with write.csv(df, file.path(out_dir, "tables", filename), row.names=FALSE)
  g) Keep ALL functions exactly as they are:
     - get_foundational_config()
     - compute_crossomics_feature_correlations()
     - map_features_to_common_ids()
     - compute_pairwise_correlations()
     - compute_partial_correlations()
     - build_correlation_network()
     - compute_sample_concordance()
     - compute_sample_rank_correlations()
     - compute_clustering_consistency()
     - compute_ari(), compute_nmi()
     - compute_condition_similarity()
     - compute_consensus_clustering()
     - identify_discordant_samples()
     - analyze_pathway_overlap()
     - run_basic_enrichment()
     - compute_pathway_overlap_matrix()
     - identify_shared_specific_pathways()
     - compute_pathway_direction_concordance()
     - find_crossomics_modules()
     - All plot_*() and save_*() and summarize_*() functions
  h) The key change in get_foundational_config() is:
     Replace: fc <- config$foundational %||% list()
     With:    fc <- config$modes$multiomics$foundational %||% list()

TASK 3: Create R/domain/multiomics/12_mechanistic_inference.R
  Source: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/R/05c_mechanistic_inference.R (~1,403 lines)

  Port with these changes:
  a) New main function: run_mechanistic_analysis(mae, foundational_results = NULL, config, gene_protein_mapping = NULL, out_dir = NULL)
  b) Use .mae_to_legacy() internally
  c) Replace config$mechanistic$* with config$modes$multiomics$mechanistic$*
  d) Replace log_message() with message()
  e) Replace save_plot()/save_table() with ggsave()/write.csv() using out_dir
  f) Keep ALL functions intact:
     - get_mechanistic_config()
     - analyze_rna_protein_regulation()
     - map_rna_protein_features()
     - compute_rna_protein_correlations()
     - compute_translation_efficiency()
     - detect_post_transcriptional_regulation()
     - infer_protein_stability()
     - infer_tf_activity()
     - get_tf_regulons()
     - infer_regulatory_network()
     - run_mediation_analysis()
     - All plot_*() and save_*() functions

TASK 4: Create R/modules/multiomics/06_mod_foundational.R
  New orchestration wrapper:

  ```r
  #' Module: Foundational cross-omics analysis
  #'
  #' Runs pre-integration cross-omics correlations, sample concordance,
  #' pathway overlap, and co-expression module detection.
  #'
  #' @param harmonization_res Output from mod_multiomics_harmonization()
  #' @param config Full config object
  #' @param out_dir Output directory
  #' @return List with foundational analysis results, or NULL if disabled/failed
  mod_multiomics_foundational <- function(harmonization_res, config, out_dir) {
      message("\n=== Foundational Cross-Omics Analysis ===\n")

      fc <- config$modes$multiomics$foundational %||% list()
      if (!(fc$run_foundational %||% TRUE)) {
          message("Foundational analysis disabled in config. Skipping.")
          return(NULL)
      }

      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      tryCatch({
          run_foundational_analysis(
              mae = harmonization_res$mae,
              config = config,
              gene_protein_mapping = harmonization_res$gene_protein_mapping,
              out_dir = out_dir
          )
      }, error = function(e) {
          warning("Foundational analysis failed: ", e$message)
          NULL
      })
  }
  ```

TASK 5: Create R/modules/multiomics/08_mod_mechanistic.R
  Similar wrapper for mechanistic analysis:

  ```r
  #' Module: Mechanistic inference analysis
  #'
  #' @param harmonization_res Output from mod_multiomics_harmonization()
  #' @param foundational_results Output from mod_multiomics_foundational()
  #' @param config Full config object
  #' @param out_dir Output directory
  #' @return List with mechanistic analysis results, or NULL if disabled/failed
  mod_multiomics_mechanistic <- function(harmonization_res, foundational_results = NULL,
                                          config, out_dir) {
      message("\n=== Mechanistic Inference Analysis ===\n")

      mc <- config$modes$multiomics$mechanistic %||% list()
      if (!(mc$run_mechanistic %||% FALSE)) {
          message("Mechanistic analysis disabled in config. Skipping.")
          return(NULL)
      }

      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      tryCatch({
          run_mechanistic_analysis(
              mae = harmonization_res$mae,
              foundational_results = foundational_results,
              config = config,
              gene_protein_mapping = harmonization_res$gene_protein_mapping,
              out_dir = out_dir
          )
      }, error = function(e) {
          warning("Mechanistic analysis failed: ", e$message)
          NULL
      })
  }
  ```

=== IMPORTANT RULES ===
- READ each source file FULLY before porting. Copy analytical logic verbatim.
- Only change infrastructure: config paths, logging, save functions, data access.
- Use .mae_to_legacy() adapter to minimize changes to data access patterns.
- Every function must have roxygen (#') documentation.
- Do NOT modify any existing files in R/core/, R/pipeline/, or other R/domain/ files.
- Do NOT add source() calls — all files are auto-sourced.
- After creating each file, verify syntax: Rscript -e 'source("R/domain/multiomics/<file>.R")'
- The .mae_to_legacy() helper should be defined in each file that needs it (as a dot-prefixed private function) to avoid cross-file dependencies.

=== WHEN COMPLETE ===
After finishing all tasks and passing syntax checks, output this status update:

✅ AGENT 2 COMPLETE:
  - Created: R/domain/multiomics/06b_rna_protein_correlation.R
  - Created: R/domain/multiomics/09_foundational_correlations.R
  - Created: R/domain/multiomics/12_mechanistic_inference.R
  - Created: R/modules/multiomics/06_mod_foundational.R
  - Created: R/modules/multiomics/08_mod_mechanistic.R
  - Syntax checks: PASSED
→ Agent 3 can now start.
```

---

## AGENT 3 — Consensus/Stability Analysis + Pipeline Wiring + Commentary

```
You are migrating code from a standalone multi-omics pipeline into the multiomics-core R project.

SOURCE REPO: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/
DESTINATION REPO: /home/ozsol/multiomics-core/

=== CONTEXT: multiomics-core Architecture ===

multiomics-core uses a 3-layer architecture:
  R/core/       → shared utilities — DO NOT MODIFY
  R/domain/     → pure analytical functions (stateless, per-omics subdirs)
  R/modules/    → orchestration wrappers
  R/pipeline/   → {targets} DAG definitions (this agent OWNS 00_pipe_multiomics.R)
  R/services/   → external services (AI commentary)

Key conventions:
- Config: config$modes$multiomics$* (NOT config$integration$* directly)
- %||% available from R/core/04_config.R
- Logging: message() (NOT log_message())
- Output: out_dir parameter (NOT config$output$output_dir)
- Data: formal MultiAssayExperiment
- All R/ files auto-sourced — NO source() calls
- design$condition_column, design$contrasts, design$reference_level are at top level config$design$*

=== DATA STRUCTURE BRIDGE ===

Define this adapter at top of each new domain file (dot-prefixed for privacy):

```r
.mae_to_legacy <- function(mae, de_results = NULL, gene_protein_mapping = NULL) {
    harmonized_omics <- lapply(names(experiments(mae)), function(nm) {
        exp_data <- experiments(mae)[[nm]]
        de <- if (!is.null(de_results[[nm]])) de_results[[nm]] else NULL
        list(
            normalized_matrix = as.matrix(assay(exp_data)),
            de_table = de, da_table = de,
            feature_annotation = as.data.frame(rowData(exp_data))
        )
    })
    names(harmonized_omics) <- names(experiments(mae))
    list(mae = mae, harmonized_omics = harmonized_omics,
         metadata = as.data.frame(colData(mae)),
         common_samples = colnames(mae), gene_mapping = gene_protein_mapping)
}
```

=== EXISTING PIPELINE FILE (you will MODIFY this) ===

Read R/pipeline/multiomics/00_pipe_multiomics.R. It currently defines pipe_multiomics() which returns
a list of targets: multiomics_out_dir, multiomics_harmonization, multiomics_de_results,
multiomics_enrichment_results, multiomics_integration, multiomics_concordance,
multiomics_cross_enrichment, multiomics_report.

You will ADD new targets between integration and report.

=== EXISTING MODULE FILES (for reference) ===

- R/modules/multiomics/01_mod_harmonization.R — mod_multiomics_harmonization()
- R/modules/multiomics/02_mod_integration.R   — mod_multiomics_integration()
- R/modules/multiomics/03_mod_concordance.R   — mod_multiomics_concordance()
- R/modules/multiomics/04_mod_enrichment.R    — mod_multiomics_enrichment()
- R/modules/multiomics/05_mod_report.R        — mod_multiomics_report()

These return results used by the pipeline. Read them to understand the data flow.

=== YOUR TASKS ===

TASK 1: Create R/domain/multiomics/10_integration_consensus.R
  Source: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/R/09b_integration_consensus.R (~1,041 lines)

  Port with these changes:
  a) New main function: run_integration_consensus(integration_results, mae, config, gene_protein_mapping = NULL, out_dir = NULL)
     The integration_results parameter format from multiomics-core is:
       list(diablo_results = ..., snf_results = ..., mofa_results = ..., feature_selection = ..., mae_subset = ...)
     The source expects: list(mofa = ..., diablo = ..., snf = ...)
     So adapt at the start: normalize the keys.
  b) Use .mae_to_legacy() where mae_data is needed
  c) Replace config$consensus$* with config$modes$multiomics$consensus$*
  d) Replace log_message() with message()
  e) Replace save_plot()/save_table() with ggsave()/write.csv() using out_dir
  f) Keep ALL functions:
     - compare_sample_clustering()
     - get_mofa_clusters(), get_diablo_clusters(), get_snf_clusters()
     - choose_optimal_k(), choose_optimal_k_spectral()
     - compute_ari(), compute_nmi(), compute_consensus_clustering()
     - compare_feature_importance()
     - get_mofa_feature_importance(), get_diablo_feature_importance()
     - identify_robust_patterns(), identify_method_specific_patterns()
     - run_meta_integration()
     - save_consensus_outputs(), all plot_*() functions

TASK 2: Create R/domain/multiomics/11_stability_analysis.R
  Source: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/R/09c_stability_analysis.R (~1,096 lines)

  Port with these changes:
  a) New main function: run_stability_analysis(mae, feature_selection, integration_results, config, out_dir = NULL)
     The feature_selection object is from mod_multiomics_integration()$feature_selection.
     The integration_results keys need the same normalization as Task 1.
  b) Use .mae_to_legacy() adapter
  c) Replace config$stability$* with config$modes$multiomics$stability$*
  d) Replace log_message() with message()
  e) Replace save_plot()/save_table() with ggsave()/write.csv() using out_dir
  f) Keep ALL functions:
     - run_bootstrap_feature_stability()
     - bootstrap_mofa_features(), bootstrap_diablo_features()
     - stratified_bootstrap_indices()
     - create_mofa_for_bootstrap(), run_diablo_for_bootstrap()
     - run_leave_one_out_analysis()
     - run_kfold_stability(), create_balanced_folds()
     - compute_loading_confidence_intervals()
     - assess_cluster_stability()
     - generate_stability_plots()
     - All plot_*() functions
     - create_stability_summary(), save_stability_outputs()

TASK 3: Create R/domain/multiomics/13_commentary.R
  Source: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/R/11_commentary.R (~670 lines)
  Source: /home/ozsol/multiomics/src/pipelines/multiomics_pipeline/R/11b_commentary_fallbacks.R (~512 lines)

  ⚠️ IMPORTANT: multiomics-core already has R/services/12_commentary.R for single-omics commentary.
  READ THAT FILE FIRST to check for existing function names and avoid collisions.

  Port the multiomics-specific commentary as an EXTENSION:
  a) Create R/domain/multiomics/13_commentary.R with multiomics-specific functions.
     PREFIX ALL FUNCTIONS with "multiomics_" to avoid name collisions:
     - multiomics_build_figures_table() — catalog all multi-omics figures
     - multiomics_generate_commentary() — orchestrate commentary for multi-omics figures
     - multiomics_build_figure_context() — build context for multi-omics figure types
  b) Create R/domain/multiomics/13b_commentary_fallbacks.R with:
     - multiomics_generate_fallback_commentary() — router
     - multiomics_fallback_mofa_*() functions for MOFA figure types
     - multiomics_fallback_diablo_*() functions for DIABLO figure types
     - multiomics_fallback_concordance_*() for concordance figures
     - multiomics_fallback_consensus_*() for consensus figures
  c) Replace config$commentary$* with config$modes$multiomics$commentary$*
  d) Replace source("R/11b_commentary_fallbacks.R") — NOT NEEDED (auto-sourced)
  e) For Claude/OpenAI API calls, delegate to the existing R/services/12_commentary.R functions
     if they exist, or use system() calls to scripts/ as fallback.
  f) Replace log_message() with message()

TASK 4: Create R/modules/multiomics/07_mod_consensus.R
  New orchestration wrapper:

  ```r
  #' Module: Integration consensus and stability analysis
  #'
  #' @param integration_res Output from mod_multiomics_integration()
  #' @param harmonization_res Output from mod_multiomics_harmonization()
  #' @param config Full config object
  #' @param out_dir Output directory
  #' @return List with consensus and stability results
  mod_multiomics_consensus <- function(integration_res, harmonization_res, config, out_dir) {
      message("\n=== Integration Consensus & Stability ===\n")
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      results <- list(consensus = NULL, stability = NULL)

      # Normalize integration results to source format
      integration_results <- list(
          mofa = integration_res$mofa_results,
          diablo = integration_res$diablo_results,
          snf = integration_res$snf_results
      )

      # Consensus analysis
      cc <- config$modes$multiomics$consensus %||% list()
      if (cc$compare_methods %||% TRUE) {
          results$consensus <- tryCatch({
              run_integration_consensus(
                  integration_results = integration_results,
                  mae = harmonization_res$mae,
                  config = config,
                  gene_protein_mapping = harmonization_res$gene_protein_mapping,
                  out_dir = file.path(out_dir, "consensus")
              )
          }, error = function(e) {
              warning("Consensus analysis failed: ", e$message)
              NULL
          })
      }

      # Stability analysis
      sc <- config$modes$multiomics$stability %||% list()
      if (sc$run_stability %||% FALSE) {
          results$stability <- tryCatch({
              run_stability_analysis(
                  mae = harmonization_res$mae,
                  feature_selection = integration_res$feature_selection,
                  integration_results = integration_results,
                  config = config,
                  out_dir = file.path(out_dir, "stability")
              )
          }, error = function(e) {
              warning("Stability analysis failed: ", e$message)
              NULL
          })
      }

      results
  }
  ```

TASK 5: Update R/pipeline/multiomics/00_pipe_multiomics.R
  Read the existing file. Add NEW targets while keeping all existing ones.

  Add these targets IN ORDER, after the existing multiomics_integration target
  and BEFORE multiomics_report:

  a) multiomics_foundational — runs foundational analysis (depends on multiomics_harmonization)
     ```r
     tar_target(
         multiomics_foundational,
         {
             if (is.null(multiomics_harmonization)) return(NULL)
             mod_multiomics_foundational(
                 harmonization_res = multiomics_harmonization,
                 config = config,
                 out_dir = file.path(multiomics_out_dir, "foundational")
             )
         }
     ),
     ```

  b) multiomics_rna_protein_corr — extended RNA-protein correlation (depends on harmonization + DE)
     ```r
     tar_target(
         multiomics_rna_protein_corr,
         {
             if (is.null(multiomics_harmonization)) return(NULL)
             tryCatch({
                 run_rna_protein_correlation(
                     mae = multiomics_harmonization$mae,
                     de_results = multiomics_de_results,
                     gene_protein_mapping = multiomics_harmonization$gene_protein_mapping,
                     config = config,
                     out_dir = file.path(multiomics_out_dir, "rna_protein_correlation")
                 )
             }, error = function(e) {
                 warning("RNA-protein correlation failed: ", e$message)
                 NULL
             })
         }
     ),
     ```

  c) multiomics_consensus — consensus + stability (depends on integration + harmonization)
     ```r
     tar_target(
         multiomics_consensus,
         {
             if (is.null(multiomics_integration)) return(NULL)
             mod_multiomics_consensus(
                 integration_res = multiomics_integration,
                 harmonization_res = multiomics_harmonization,
                 config = config,
                 out_dir = file.path(multiomics_out_dir, "consensus")
             )
         }
     ),
     ```

  d) multiomics_mechanistic — mechanistic inference (depends on harmonization + foundational)
     ```r
     tar_target(
         multiomics_mechanistic,
         {
             if (is.null(multiomics_harmonization)) return(NULL)
             mc <- config$modes$multiomics$mechanistic %||% list()
             if (!(mc$run_mechanistic %||% FALSE)) return(NULL)
             mod_multiomics_mechanistic(
                 harmonization_res = multiomics_harmonization,
                 foundational_results = multiomics_foundational,
                 config = config,
                 out_dir = file.path(multiomics_out_dir, "mechanistic")
             )
         }
     ),
     ```

  e) multiomics_multigsea — MultiGSEA plots (depends on enrichment)
     ```r
     tar_target(
         multiomics_multigsea,
         {
             if (is.null(multiomics_cross_enrichment)) return(NULL)
             tryCatch({
                 run_multigsea_plots(
                     enrichment_results = multiomics_cross_enrichment,
                     config = config,
                     out_dir = file.path(multiomics_out_dir, "multigsea")
                 )
             }, error = function(e) {
                 warning("MultiGSEA plots failed: ", e$message)
                 NULL
             })
         }
     ),
     ```

  f) Update the multiomics_report target to force() all new targets:
     ```r
     force(multiomics_foundational)
     force(multiomics_rna_protein_corr)
     force(multiomics_consensus)
     force(multiomics_mechanistic)
     force(multiomics_multigsea)
     ```

=== IMPORTANT RULES ===
- READ every source file FULLY before porting. Copy analytical logic verbatim.
- Only change infrastructure: config paths, logging, save functions, data access.
- Use .mae_to_legacy() adapter in domain files to minimize data access changes.
- Every function must have roxygen (#') documentation.
- Do NOT modify R/core/ files or files created by other agents.
- The pipeline file MUST preserve all existing targets — only ADD new ones.
- After all files are created, verify syntax errors:
  for f in R/domain/multiomics/10_integration_consensus.R R/domain/multiomics/11_stability_analysis.R R/domain/multiomics/13_commentary.R R/domain/multiomics/13b_commentary_fallbacks.R R/modules/multiomics/07_mod_consensus.R; do echo "=== $f ===" && Rscript -e "source('$f')" 2>&1 | head -5; done
- Test the pipeline definition parses: Rscript -e 'library(targets); source("R/pipeline/multiomics/00_pipe_multiomics.R"); cat("OK\n")'

=== WHEN COMPLETE ===
After finishing all tasks and passing syntax/pipeline checks, output this status update:

✅ AGENT 3 COMPLETE:
  - Created: R/domain/multiomics/10_integration_consensus.R
  - Created: R/domain/multiomics/11_stability_analysis.R
  - Created: R/domain/multiomics/13_commentary.R
  - Created: R/domain/multiomics/13b_commentary_fallbacks.R
  - Created: R/modules/multiomics/07_mod_consensus.R
  - Updated: R/pipeline/multiomics/00_pipe_multiomics.R
  - Pipeline parse check: PASSED
→ Ready for full pipeline test: tar_make()
```

---

## Coordination Notes

### Execution Order (CRITICAL)

```
┌─────────────────────────────────────────────────────────────────┐
│  Agent 1 (FIRST)                                                │
│  Creates: run_mofa2_integration(), write_mofa_results(),        │
│           run_multigsea_plots(), MOFA2 config validation        │
└─────────────────────┬───────────────────────────────────────────┘
                      │ Must complete before Agent 2
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│  Agent 2 (SECOND)                                               │
│  Creates: run_foundational_analysis(), run_mechanistic_analysis │
│           run_rna_protein_correlation(), module wrappers        │
└─────────────────────┬───────────────────────────────────────────┘
                      │ Must complete before Agent 3
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│  Agent 3 (LAST)                                                 │
│  Creates: consensus/stability analysis, commentary              │
│  Updates: 00_pipe_multiomics.R (calls functions from 1 & 2)     │
└─────────────────────────────────────────────────────────────────┘
```

### File Ownership

| What | Agent 1 | Agent 2 | Agent 3 |
|------|---------|---------|---------|
| **Domain files** | 04b_integration_mofa, 07b_multigsea, 07c_pathview | 06b_rna_protein_corr, 09_foundational, 12_mechanistic | 10_consensus, 11_stability, 13_commentary, 13b_fallbacks |
| **Module files** | — | 06_mod_foundational, 08_mod_mechanistic | 07_mod_consensus |
| **Pipeline file** | — | — | 00_pipe_multiomics.R (sole owner) |
| **Config updates** | 90_config_validate.R, multiomics_GT15_test.yaml | — | — |
| **Integration module** | 02_mod_integration.R (add MOFA2) | — | — |
| **Scripts** | scripts/run_mofa.py | — | — |

**No file conflicts between agents.** Each agent owns distinct files.

### Inter-Agent Status Updates

When completing major steps, agents should post a brief status update to a shared Slack channel
(or terminal log) so the next agent knows when to start:

**Agent 1 should announce when done:**
```
✅ AGENT 1 COMPLETE:
  - Created: R/domain/multiomics/04b_integration_mofa.R
  - Created: R/domain/multiomics/07b_multigsea_plots.R
  - Created: R/domain/multiomics/07c_kegg_pathview.R
  - Copied: scripts/run_mofa.py
  - Updated: R/modules/multiomics/02_mod_integration.R (added MOFA2)
  - Updated: R/domain/multiomics/90_config_validate.R
  - Updated: config/multiomics_GT15_test.yaml
  - Syntax checks: PASSED
→ Agent 2 can now start.
```

**Agent 2 should announce when done:**
```
✅ AGENT 2 COMPLETE:
  - Created: R/domain/multiomics/06b_rna_protein_correlation.R
  - Created: R/domain/multiomics/09_foundational_correlations.R
  - Created: R/domain/multiomics/12_mechanistic_inference.R
  - Created: R/modules/multiomics/06_mod_foundational.R
  - Created: R/modules/multiomics/08_mod_mechanistic.R
  - Syntax checks: PASSED
→ Agent 3 can now start.
```

**Agent 3 should announce when done:**
```
✅ AGENT 3 COMPLETE:
  - Created: R/domain/multiomics/10_integration_consensus.R
  - Created: R/domain/multiomics/11_stability_analysis.R
  - Created: R/domain/multiomics/13_commentary.R
  - Created: R/domain/multiomics/13b_commentary_fallbacks.R
  - Created: R/modules/multiomics/07_mod_consensus.R
  - Updated: R/pipeline/multiomics/00_pipe_multiomics.R
  - Pipeline parse check: PASSED
→ Ready for full pipeline test: tar_make()
```

### Cross-Agent Dependencies

| Agent 3 Pipeline Target | Calls Function | Created By |
|-------------------------|----------------|------------|
| multiomics_multigsea | `run_multigsea_plots()` | Agent 1 |
| multiomics_foundational | `mod_multiomics_foundational()` | Agent 2 |
| multiomics_rna_protein_corr | `run_rna_protein_correlation()` | Agent 2 |
| multiomics_mechanistic | `mod_multiomics_mechanistic()` | Agent 2 |
| multiomics_consensus | `mod_multiomics_consensus()` → calls `run_integration_consensus()` | Agent 3 |

### Post-Completion Validation

After **Agent 1** completes:
```bash
Rscript -e 'source("R/domain/multiomics/04b_integration_mofa.R")'
Rscript -e 'source("R/domain/multiomics/07b_multigsea_plots.R")'
Rscript -e 'source("R/domain/multiomics/07c_kegg_pathview.R")'
```

After **Agent 2** completes:
```bash
Rscript -e 'source("R/domain/multiomics/06b_rna_protein_correlation.R")'
Rscript -e 'source("R/domain/multiomics/09_foundational_correlations.R")'
Rscript -e 'source("R/domain/multiomics/12_mechanistic_inference.R")'
```

After **Agent 3** completes:
```bash
Rscript -e 'library(targets); source("R/pipeline/multiomics/00_pipe_multiomics.R"); cat("Pipeline OK\n")'
```

After **ALL agents** complete:
```bash
cd /home/ozsol/multiomics-core
Rscript -e 'targets::tar_visnetwork()'  # Verify DAG structure
Rscript -e 'targets::tar_destroy(); targets::tar_make()'  # Full pipeline run
```
