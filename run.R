#!/usr/bin/env Rscript
# =============================================================================
# multiomics-core Pipeline Runner
#
# Usage:
#   Rscript run.R --new              # Interactive wizard for new project
#   Rscript run.R --config path.yaml # Run with existing config
#   Rscript run.R                    # Re-run last config (with caching)
#   Rscript run.R --fresh            # Full re-run (invalidate cache)
#
# In RStudio: source("run.R") runs interactively
# =============================================================================

# --- Helpers ------------------------------------------------------------------

cli_header <- function() {
  cat("\n")
  cat("===========================================================\n")
  cat("   multiomics-core Pipeline\n")
  cat("===========================================================\n\n")
}

## Persistent stdin connection for non-interactive (Rscript) mode
.stdin_con <- NULL
get_stdin_con <- function() {
  if (is.null(.stdin_con)) {
    .stdin_con <<- file("stdin", "r")
  }
  .stdin_con
}

prompt_read <- function(prompt_text) {
  cat(prompt_text)
  if (interactive()) {
    ans <- readline("")
  } else {
    ans <- readLines(get_stdin_con(), n = 1, warn = FALSE)
    if (length(ans) == 0) {
      cat("\n  [End of input — exiting]\n")
      quit(save = "no", status = 1)
    }
  }
  trimws(ans)
}

ask <- function(prompt, default = NULL) {
  if (!is.null(default) && nzchar(default)) {
    prompt_text <- paste0(prompt, " [", default, "]: ")
  } else {
    prompt_text <- paste0(prompt, ": ")
  }
  ans <- prompt_read(prompt_text)
  if (ans == "" && !is.null(default)) return(default)
  if (ans == "") {
    cat("  (input required)\n")
    return(ask(prompt, default))
  }
  ans
}

ask_choice <- function(prompt, choices, default = 1) {
  cat(prompt, "\n")
  for (i in seq_along(choices)) {
    marker <- if (i == default) " (default)" else ""
    cat(sprintf("  %d) %s%s\n", i, choices[i], marker))
  }
  ans <- prompt_read(sprintf("Choose [%d]: ", default))
  if (ans == "") return(default)
  idx <- suppressWarnings(as.integer(ans))
  if (is.na(idx) || idx < 1 || idx > length(choices)) {
    cat("Invalid choice, using default.\n")
    return(default)
  }
  idx
}

ask_yn <- function(prompt, default = TRUE) {
  hint <- if (default) "Y/n" else "y/N"
  ans <- prompt_read(paste0(prompt, " [", hint, "]: "))
  if (ans == "") return(default)
  tolower(ans) %in% c("y", "yes")
}

validate_file <- function(path, label) {
  path <- trimws(path)
  if (!file.exists(path)) {
    cat(sprintf("  WARNING: %s not found at: %s\n", label, path))
    return(NULL)
  }
  normalizePath(path)
}

detect_separator <- function(path) {
  first_line <- readLines(path, n = 1)
  if (grepl("\t", first_line)) return("\t")
  if (grepl(",", first_line)) return(",")
  "\t"
}

detect_columns <- function(path) {
  sep <- detect_separator(path)
  header <- read.table(path, sep = sep, header = FALSE, nrows = 1,
                        stringsAsFactors = FALSE, check.names = FALSE)
  as.character(header[1, ])
}

# --- Interactive Wizard -------------------------------------------------------

wizard_new_project <- function(project_dir) {
  cli_header()
  cat("--- New Project Setup ---\n\n")

  # Project info
  project_name <- ask("Project name (no spaces)")
  analyst <- ask("Analyst name", "Bioinformatics Core")
  round <- ask("Analysis round", "A01")

  # Mode choice
  cat("\n--- Analysis Type ---\n")
  mode_idx <- ask_choice("What type of analysis?",
    c("RNA-seq",
      "Proteomics"),
    default = 1)
  analysis_mode <- c("rna", "proteomics")[mode_idx]

  if (analysis_mode == "rna") {
    wizard_rna(project_dir, project_name, analyst, round)
  } else {
    wizard_proteomics(project_dir, project_name, analyst, round)
  }
}

# --- RNA-seq Wizard -----------------------------------------------------------

wizard_rna <- function(project_dir, project_name, analyst, round) {

  # Data directory
  cat("\n--- RNA-seq Data Files ---\n")
  cat("Provide ABSOLUTE paths to your input files.\n\n")

  counts_path <- ask("Path to counts file (genes x samples, CSV/TSV)")
  counts_path <- validate_file(counts_path, "Counts file")
  if (is.null(counts_path)) {
    counts_path <- ask("Please re-enter the counts file path")
    counts_path <- validate_file(counts_path, "Counts file")
  }

  metadata_path <- ask("Path to metadata file (sample info, CSV/TSV)")
  metadata_path <- validate_file(metadata_path, "Metadata file")

  contrasts_path <- ask("Path to contrasts file (CSV/TSV, or 'none' to skip)", "none")
  if (tolower(contrasts_path) == "none") {
    contrasts_path <- NULL
  } else {
    contrasts_path <- validate_file(contrasts_path, "Contrasts file")
  }

  sample_map_path <- ask("Path to sample map file (optional, or 'none')", "none")
  if (tolower(sample_map_path) == "none") {
    sample_map_path <- NULL
  } else {
    sample_map_path <- validate_file(sample_map_path, "Sample map file")
  }

  # Auto-detect columns from counts and metadata
  cat("\n--- Column Detection ---\n")
  if (!is.null(counts_path)) {
    counts_cols <- detect_columns(counts_path)
    cat("Counts file columns: ", paste(head(counts_cols, 5), collapse = ", "),
        if (length(counts_cols) > 5) "..." else "", "\n")
    gene_id_col <- ask("Gene ID column name in counts", counts_cols[1])
  } else {
    gene_id_col <- ask("Gene ID column name in counts", "gene_id")
  }

  if (!is.null(metadata_path)) {
    meta_cols <- detect_columns(metadata_path)
    cat("Metadata columns: ", paste(meta_cols, collapse = ", "), "\n")
    sample_col <- ask("Sample ID column in metadata", meta_cols[1])

    # Group column for filtering and effects
    cat("\nWhich column defines your biological groups/conditions?\n")
    cat("(Used for filtering, PCA coloring, and DE contrasts)\n")
    group_col <- ask("Group column", meta_cols[min(2, length(meta_cols))])
  } else {
    sample_col <- ask("Sample ID column in metadata", "SampleID")
    group_col <- ask("Group/condition column in metadata")
  }

  # Methods
  cat("\n--- Analysis Methods ---\n")

  norm_idx <- ask_choice("Normalization method:",
    c("TMMlogCPM — TMM scaling + log2 CPM (recommended for most RNA-seq)",
      "VST — Variance Stabilizing Transformation (DESeq2)"),
    default = 1)
  norm_method <- c("TMMlogCPM", "VST")[norm_idx]

  filter_idx <- ask_choice("Feature filtering strategy:",
    c("Adaptive KDE — auto-detect noise/signal threshold (recommended)",
      "DESeq2 built-in only — no pre-filter, matches standard facility pipelines",
      "Fixed CPM threshold — manual cutoff (classic approach)"),
    default = 1)
  filter_mode <- c("adaptive", "deseq2_only", "fixed")[filter_idx]
  filter_cpm_threshold <- ""
  if (filter_mode == "fixed") {
    filter_cpm_threshold <- ask("CPM threshold (genes must exceed this in ≥50% reps per group)", "1.0")
  }

  de_idx <- ask_choice("Differential expression method:",
    c("DESeq2 - default (recommended)",
      "DESeq2 - legacy (betaPrior=TRUE)",
      "limma-voom"),
    default = 1)
  de_method <- c("deseq2", "deseq2", "limma_voom")[de_idx]
  deseq_mode <- c("default", "legacy", "default")[de_idx]

  p_cutoff <- as.numeric(ask("Adjusted p-value cutoff", "0.05"))
  fc_cutoff <- as.numeric(ask("Linear fold-change cutoff", "1.5"))

  pathway_idx <- ask_choice("Pathway analysis:",
    c("fGSEA (recommended)",
      "ORA (Over-Representation Analysis)",
      "Both",
      "Skip"),
    default = 1)
  pathway_method <- c("fgsea", "ora", "both", "none")[pathway_idx]
  pathway_enabled <- pathway_idx != 4

  clustering_on <- ask_yn("Enable clustering?", FALSE)
  batch_corr_on <- ask_yn("Enable batch correction? (detects and corrects batch effects)", FALSE)
  deconv_on <- ask_yn("Enable cell type deconvolution? (human/mouse only, uses xCell2)", FALSE)

  # Organism & Annotation
  cat("\n--- Organism & Annotation ---\n")
  org_idx <- ask_choice("Which organism does your data come from?",
    c("Auto-detect (recommended)",
      "Human (Homo sapiens)",
      "Mouse (Mus musculus)",
      "C. elegans (Caenorhabditis elegans)",
      "Giardia lamblia",
      "Other (type name)"),
    default = 1)

  organism_names <- c("auto", "Homo sapiens", "Mus musculus",
                      "Caenorhabditis elegans", "Giardia lamblia", "other")
  selected_organism <- organism_names[org_idx]

  if (selected_organism == "other") {
    selected_organism <- ask("Organism name (e.g. 'Rattus norvegicus')")
  }

  custom_mapping_file <- "null"
  custom_gmt_file <- "null"

  # For non-model organisms, offer custom annotation/GMT files
  is_non_model <- selected_organism %in% c("Giardia lamblia") ||
    (selected_organism != "auto" &&
     !selected_organism %in% c("Homo sapiens", "Mus musculus",
                                "Caenorhabditis elegans", "Rattus norvegicus",
                                "Danio rerio", "Drosophila melanogaster",
                                "Saccharomyces cerevisiae", "Arabidopsis thaliana"))

  if (is_non_model) {
    cat("\nNote: This is a non-model organism. Annotation via biomaRt/OrgDb may not\n")
    cat("be available. You can provide custom files for gene names and pathways.\n\n")
    mapping_ans <- ask("Path to custom annotation file (CSV/TSV with gene_id + gene_name columns, or 'none')", "none")
    if (tolower(mapping_ans) != "none") {
      mapping_path <- validate_file(mapping_ans, "Custom annotation file")
      if (!is.null(mapping_path)) {
        custom_mapping_file <- paste0('"', mapping_path, '"')
      }
    }
    gmt_ans <- ask("Path to custom GMT file for pathway analysis (or 'none')", "none")
    if (tolower(gmt_ans) != "none") {
      gmt_path <- validate_file(gmt_ans, "Custom GMT file")
      if (!is.null(gmt_path)) {
        custom_gmt_file <- paste0('"', gmt_path, '"')
      }
    }
  }

  # --- AI Commentary ---
  cat("\n--- AI Commentary ---\n")
  cat("Generate AI-powered interpretation for each figure in the report?\n")
  cat("Requires an API key (ANTHROPIC_API_KEY or OPENAI_API_KEY).\n")
  commentary_idx <- ask_choice("Commentary backend",
                                c("None (data-driven fallback only)",
                                  "Claude Code (uses your subscription, no API key needed)",
                                  "Claude API (requires ANTHROPIC_API_KEY)",
                                  "OpenAI API (requires OPENAI_API_KEY)"),
                                default = 1)
  commentary_backends <- c("none", "claude-code", "claude", "openai")
  commentary_backend <- commentary_backends[commentary_idx]
  commentary_enabled <- commentary_backend != "none"

  # --- Technical Report (Optional) ---
  nna <- function(x, default = "") if (is.null(x) || !nzchar(x)) default else x

  cat("\n--- Technical Report (Optional) ---\n")
  cat("Provide a facility technical report (DOCX/PDF) to auto-extract\n")
  cat("library prep, sequencing, and bioinformatics details.\n")
  tech_report_path <- ask("Path to technical report file (or 'none')", "none")
  tech_report_fields <- NULL

  if (tolower(tech_report_path) != "none") {
    tech_report_path <- validate_file(tech_report_path, "Technical report")

    if (!is.null(tech_report_path)) {
      # Check if claude CLI is available
      claude_available <- nzchar(Sys.which("claude"))

      if (!claude_available) {
        cat("  Claude Code CLI not found — cannot auto-extract fields.\n")
        cat("  You can manually add a 'technical_report:' block to the config later.\n")
      } else {
        # Extract text via pandoc (guaranteed available via rmarkdown)
        cat("  Extracting text from report...\n")
        doc_text <- tryCatch(
          paste(system2("pandoc", c("-t", "plain", shQuote(tech_report_path)),
                        stdout = TRUE, stderr = FALSE), collapse = "\n"),
          error = function(e) NULL
        )

        if (is.null(doc_text) || nchar(doc_text) < 50) {
          cat("  WARNING: Could not extract sufficient text from the file.\n")
          cat("  You can manually add a 'technical_report:' block to the config later.\n")
        } else {
          cat(sprintf("  Extracted %d characters. Sending to Claude for parsing...\n",
                      nchar(doc_text)))

          json_schema <- paste0(
            '{"type":"object","properties":{',
            '"facility":{"type":"string","description":"Name of the sequencing facility"},',
            '"library_prep":{"type":"string","description":"Library preparation method, kit, input amount, QC"},',
            '"sequencing":{"type":"string","description":"Sequencing platform, kit, read configuration"},',
            '"upstream_bioinformatics":{"type":"string","description":"QC, trimming, alignment, quantification tools and parameters"},',
            '"upstream_de_note":{"type":"string","description":"Initial DE analysis method, normalization, filtering, and DEG count if mentioned"},',
            '"acknowledgment":{"type":"string","description":"Facility acknowledgment sentence"}',
            '},"required":["facility","library_prep","sequencing","upstream_bioinformatics"]}'
          )

          prompt <- sprintf(
            paste0(
              "Extract structured information from this sequencing facility technical report. ",
              "Return only the fields described in the schema. ",
              "For each field, use the exact details from the report text. ",
              "If a field is not mentioned, use an empty string.\n\n",
              "REPORT TEXT:\n%s"
            ),
            doc_text
          )

          cmd <- sprintf(
            "claude --print --output-format json --model sonnet --json-schema '%s' --no-session-persistence %s",
            json_schema,
            shQuote(prompt)
          )

          parsed <- tryCatch({
            raw <- system(cmd, intern = TRUE, timeout = 120)
            raw_text <- paste(raw, collapse = "\n")
            p <- jsonlite::fromJSON(raw_text, simplifyVector = FALSE)
            if (!is.null(p$structured_output)) p$structured_output
            else if (!is.null(p$result) && nzchar(p$result)) {
              jsonlite::fromJSON(p$result, simplifyVector = FALSE)
            } else NULL
          }, error = function(e) {
            cat(sprintf("  Claude extraction failed: %s\n", e$message))
            NULL
          })

          if (!is.null(parsed) && !is.null(parsed$facility)) {
            cat("\n  Extracted fields:\n")
            cat(sprintf("    Facility:      %s\n", parsed$facility))
            cat(sprintf("    Library prep:  %s\n",
                        substr(nna(parsed$library_prep), 1, 80)))
            cat(sprintf("    Sequencing:    %s\n",
                        substr(nna(parsed$sequencing), 1, 80)))
            cat(sprintf("    Upstream bio:  %s\n",
                        substr(nna(parsed$upstream_bioinformatics), 1, 80)))
            if (nzchar(nna(parsed$upstream_de_note)))
              cat(sprintf("    Upstream DE:   %s\n",
                          substr(parsed$upstream_de_note, 1, 80)))
            if (nzchar(nna(parsed$acknowledgment)))
              cat(sprintf("    Acknowledgment: %s\n",
                          substr(parsed$acknowledgment, 1, 80)))
            cat("  (Truncated for display — full text will be saved in config)\n")

            if (ask_yn("\n  Use these extracted fields?", TRUE)) {
              tech_report_fields <- parsed
            } else {
              cat("  Skipped. You can manually add a 'technical_report:' block later.\n")
            }
          } else {
            cat("  Could not extract structured fields from the report.\n")
            cat("  You can manually add a 'technical_report:' block to the config later.\n")
          }
        }
      }
    }
  }

  # Copy data into project structure
  cat("\n--- Setting Up Project ---\n")

  # Create data directory
  data_dir <- file.path(project_dir, "data", tolower(project_name), "rna")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  # Copy input files
  counts_dest <- file.path("rna", basename(counts_path))
  file.copy(counts_path, file.path(project_dir, "data", tolower(project_name), counts_dest),
            overwrite = TRUE)
  cat(sprintf("  Copied counts -> data/%s/%s\n", tolower(project_name), counts_dest))

  meta_dest <- file.path("rna", basename(metadata_path))
  file.copy(metadata_path, file.path(project_dir, "data", tolower(project_name), meta_dest),
            overwrite = TRUE)
  cat(sprintf("  Copied metadata -> data/%s/%s\n", tolower(project_name), meta_dest))

  contrasts_dest <- ""
  if (!is.null(contrasts_path)) {
    contrasts_dest <- file.path("rna", basename(contrasts_path))
    file.copy(contrasts_path, file.path(project_dir, "data", tolower(project_name), contrasts_dest),
              overwrite = TRUE)
    cat(sprintf("  Copied contrasts -> data/%s/%s\n", tolower(project_name), contrasts_dest))
  }

  sample_map_dest <- ""
  if (!is.null(sample_map_path)) {
    sample_map_dest <- file.path("rna", basename(sample_map_path))
    file.copy(sample_map_path, file.path(project_dir, "data", tolower(project_name), sample_map_dest),
              overwrite = TRUE)
    cat(sprintf("  Copied sample map -> data/%s/%s\n", tolower(project_name), sample_map_dest))
  }

  # Map from/to columns for sample mapping
  map_from <- "sample_original_name"
  map_to <- sample_col
  if (!is.null(sample_map_path)) {
    map_cols <- detect_columns(sample_map_path)
    cat("Sample map columns: ", paste(map_cols, collapse = ", "), "\n")
    map_from <- ask("Map FROM column (original sample names)", map_cols[1])
    map_to <- ask("Map TO column (unified sample IDs)", map_cols[min(2, length(map_cols))])
  }

  # Generate config YAML
  config_name <- paste0("rna_", tolower(gsub("[^a-zA-Z0-9]", "_", project_name)), ".yaml")
  config_path <- file.path(project_dir, "config", config_name)
  dir.create(file.path(project_dir, "config"), showWarnings = FALSE)

  # Shape column - try to detect
  shape_col <- "null"
  if (!is.null(metadata_path)) {
    meta_cols_all <- detect_columns(metadata_path)
    other_cols <- setdiff(meta_cols_all, c(sample_col, group_col))
    if (length(other_cols) > 0) {
      cat("\nOptional: choose a column for shape in PCA plots (or 'none'):\n")
      cat("  Available: ", paste(other_cols, collapse = ", "), "\n")
      shape_ans <- ask("Shape column", "none")
      if (tolower(shape_ans) != "none" && shape_ans %in% other_cols) {
        shape_col <- paste0('"', shape_ans, '"')
      }
    }
  }

  config_yaml <- sprintf('# Auto-generated by run.R wizard
# Project: %s
# Date: %s

project:
  dir: "%s"
  name: "%s"
  analysis_round: "%s"
  analyst: "%s"

paths:
  raw: "data/%s"
  out: "outputs/%s"

modes:
  rna:
    files:
      counts: "%s"
      sample_map: "%s"
      metadata: "%s"
      contrasts: "%s"

    id_columns:
      gene_id_col: "%s"
      sample_col: "%s"
      map_from: "%s"
      map_to: "%s"

    normalization:
      method: "%s"
      prior.count: 1

    filtering:
      group_col: "%s"
      mode: "%s"
      auto_filter:
        min: 0.5
        max: 2.0
        fallback: 1.0
      fixed_threshold: %s

    de:
      method: "%s"
      deseq_mode: "%s"
      p_cutoff: %s
      linear_fc_cutoff: %s
      use_adj: true

    de_table:
      id_col: "FeatureID"
      pass_any_col: "pass_any_contrast"

    qc_post:
      enabled: true
      de_source: "summary"
      max_top_de_features: 100
      plots:
        volcano: true
        ma: true
      outputs:
        write_de_tables: true

    annotation:
      organism: "%s"
      id_type: "auto"
      skip_annotation: false
      custom_mapping_file: %s

    pathway:
      enabled: %s
      method: "%s"
      gmt_file: %s
      databases:
        - "GO"
        - "KEGG"
      min_size: 10
      max_size: 500

    commentary:
      enabled: %s
      backend: "%s"

    clustering:
      enabled: %s
      de_source: "any_contrast"
      min_groups: 2
      steps:
        hierarchical:
          enabled: true
          distance: "euclidean"
          linkage: "complete"
        partition:
          enabled: true
          algorithm: "hclust"
          k_max: 20
        binary_patterns:
          enabled: true
          group_col: "%s"
          corr_cutoff: 0.8

    batch_correction:
      enabled: %s

    deconvolution:
      enabled: %s

    effects:
      color: "%s"
      shape: %s
      samples: "%s"

    qc:
      adaptive_plots: true
      thresholds:
        min_samples_for_pca3d: 10
        max_samples_for_heatmaps: 120
        max_samples_for_expr_heatmap: 60

params:
  seed: 1
',
    project_name, Sys.Date(), project_dir, project_name, round, analyst,
    tolower(project_name), tolower(project_name),
    counts_dest, sample_map_dest, meta_dest, contrasts_dest,
    gene_id_col, sample_col, map_from, map_to,
    norm_method, group_col, filter_mode,
    if (nzchar(filter_cpm_threshold)) filter_cpm_threshold else "null",
    de_method, deseq_mode, p_cutoff, fc_cutoff,
    selected_organism, custom_mapping_file,
    tolower(pathway_enabled), pathway_method, custom_gmt_file,
    tolower(commentary_enabled), commentary_backend,
    tolower(clustering_on), group_col,
    tolower(batch_corr_on),
    tolower(deconv_on),
    group_col, shape_col, sample_col
  )

  # Inject technical_report block if extracted
  if (!is.null(tech_report_fields)) {
    escape_yaml <- function(s) gsub('"', '\\"', nna(s), fixed = TRUE)
    tech_yaml <- sprintf(
      '    technical_report:\n      facility: "%s"\n      library_prep: "%s"\n      sequencing: "%s"\n      upstream_bioinformatics: "%s"\n      upstream_de_note: "%s"\n      acknowledgment: "%s"\n',
      escape_yaml(tech_report_fields$facility),
      escape_yaml(tech_report_fields$library_prep),
      escape_yaml(tech_report_fields$sequencing),
      escape_yaml(tech_report_fields$upstream_bioinformatics),
      escape_yaml(nna(tech_report_fields$upstream_de_note)),
      escape_yaml(nna(tech_report_fields$acknowledgment))
    )
    config_yaml <- sub("    commentary:", paste0(tech_yaml, "\n    commentary:"), config_yaml)
  }

  writeLines(config_yaml, config_path)
  cat(sprintf("\n  Config saved: %s\n", config_path))

  config_path
}

# --- Proteomics Wizard --------------------------------------------------------

wizard_proteomics <- function(project_dir, project_name, analyst, round) {

  # Data files
  cat("\n--- Proteomics Data Files ---\n")
  cat("Provide ABSOLUTE paths to your input files.\n\n")

  protein_path <- ask("Path to protein matrix file (proteins x samples, CSV/TSV)")
  protein_path <- validate_file(protein_path, "Protein matrix file")
  if (is.null(protein_path)) {
    protein_path <- ask("Please re-enter the protein matrix file path")
    protein_path <- validate_file(protein_path, "Protein matrix file")
  }

  metadata_path <- ask("Path to metadata file (sample info, CSV/TSV)")
  metadata_path <- validate_file(metadata_path, "Metadata file")

  contrasts_path <- ask("Path to contrasts file (CSV/TSV, or 'none' to skip)", "none")
  if (tolower(contrasts_path) == "none") {
    contrasts_path <- NULL
  } else {
    contrasts_path <- validate_file(contrasts_path, "Contrasts file")
  }

  sample_map_path <- ask("Path to sample map file (optional, or 'none')", "none")
  if (tolower(sample_map_path) == "none") {
    sample_map_path <- NULL
  } else {
    sample_map_path <- validate_file(sample_map_path, "Sample map file")
  }

  # Column detection
  cat("\n--- Column Detection ---\n")
  protein_id_col <- "Protein.Group"
  if (!is.null(protein_path)) {
    prot_cols <- detect_columns(protein_path)
    cat("Protein matrix columns: ", paste(head(prot_cols, 5), collapse = ", "),
        if (length(prot_cols) > 5) "..." else "", "\n")
    protein_id_col <- ask("Protein ID column", prot_cols[1])
  }

  sample_col <- "SampleID"
  group_col <- "Condition"
  if (!is.null(metadata_path)) {
    meta_cols <- detect_columns(metadata_path)
    cat("Metadata columns: ", paste(meta_cols, collapse = ", "), "\n")
    sample_col <- ask("Sample ID column in metadata", meta_cols[1])
    cat("\nWhich column defines your biological groups/conditions?\n")
    cat("(Used for PCA coloring and DE contrasts)\n")
    group_col <- ask("Group column", meta_cols[min(2, length(meta_cols))])
  }

  # Engine
  cat("\n--- Proteomics Settings ---\n")
  engine_idx <- ask_choice("Search engine used:",
    c("DIA-NN (default)", "MaxQuant", "FragPipe", "Other"),
    default = 1)
  engine <- c("DIANN", "MaxQuant", "FragPipe", "Other")[engine_idx]

  # Normalization
  norm_idx <- ask_choice("Normalization method:",
    c("None (data already normalized, e.g. DIA-NN output)",
      "Median centering"),
    default = 1)
  norm_method <- c("none", "median")[norm_idx]

  # Imputation
  imp_idx <- ask_choice("Imputation method:",
    c("Perseus-style — downshifted normal (recommended for MNAR data)",
      "DEP2 / MinDet — deterministic minimum (experimental)",
      "None (complete cases only)"),
    default = 1)
  imp_method <- c("perseus", "dep2", "none")[imp_idx]

  imp_width <- "0.2"
  imp_downshift <- "1.6"
  dep2_method <- "MinDet"
  if (imp_method == "perseus") {
    cat("  Perseus parameters (press Enter for defaults):\n")
    imp_width <- ask("    Width (SD scaling for imputed distribution)", "0.2")
    imp_downshift <- ask("    Downshift (SDs below mean)", "1.6")
  } else if (imp_method == "dep2") {
    dep2_method <- ask("  DEP2 method", "MinDet")
  }

  # DE cutoffs
  cat("\n--- Differential Expression ---\n")
  p_cutoff <- as.numeric(ask("Adjusted p-value cutoff", "0.05"))
  fc_cutoff <- as.numeric(ask("Linear fold-change cutoff", "1.5"))

  clustering_on <- ask_yn("Enable clustering?", FALSE)

  # Organism & Annotation
  cat("\n--- Organism & Annotation ---\n")
  org_idx <- ask_choice("Which organism does your data come from?",
    c("Human (Homo sapiens)",
      "Mouse (Mus musculus)",
      "Rat (Rattus norvegicus)",
      "Zebrafish (Danio rerio)",
      "Other (type name)"),
    default = 1)
  organism_names <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "other")
  selected_organism <- organism_names[org_idx]
  if (selected_organism == "other") {
    selected_organism <- ask("Organism name (Latin binomial, e.g. 'Giardia lamblia')")
  }

  # Detect non-model organisms and offer custom files
  model_organisms <- c("Homo sapiens", "Mus musculus",
                        "Rattus norvegicus", "Danio rerio",
                        "Drosophila melanogaster", "Saccharomyces cerevisiae",
                        "Caenorhabditis elegans", "Arabidopsis thaliana",
                        "Gallus gallus", "Sus scrofa", "Bos taurus")
  is_non_model <- !(selected_organism %in% model_organisms)

  custom_mapping_file <- "null"
  custom_gmt_file <- "null"

  if (is_non_model) {
    cat("\n")
    cat("  *** Non-model organism detected: ", selected_organism, " ***\n")
    cat("  Standard databases (GO via OrgDb, KEGG via clusterProfiler) are NOT\n")
    cat("  available for this organism. To get meaningful results you should provide:\n")
    cat("    - A custom annotation file (CSV/TSV with protein_id + gene_name columns)\n")
    cat("    - A custom GMT file for pathway enrichment analysis\n")
    cat("  Without these files, pathway analysis will be very limited or empty.\n\n")
    mapping_ans <- ask("Path to custom annotation file (or 'none' to skip)", "none")
    if (tolower(mapping_ans) != "none") {
      mapping_path <- validate_file(mapping_ans, "Custom annotation file")
      if (!is.null(mapping_path)) {
        custom_mapping_file <- paste0('"', mapping_path, '"')
      }
    }
    gmt_ans <- ask("Path to custom GMT file for pathway analysis (or 'none' to skip)", "none")
    if (tolower(gmt_ans) != "none") {
      gmt_path <- validate_file(gmt_ans, "Custom GMT file")
      if (!is.null(gmt_path)) {
        custom_gmt_file <- paste0('"', gmt_path, '"')
      }
    }
    if (custom_mapping_file == "null" && custom_gmt_file == "null") {
      cat("\n  Warning: No custom files provided. Annotation will use KEGG REST API\n")
      cat("  as fallback (limited coverage). Pathway analysis may produce no results.\n")
    }
  }

  # Pathway enrichment
  cat("\n--- Pathway Enrichment ---\n")
  pathway_idx <- ask_choice("Pathway analysis:",
    c("Both fGSEA + ORA (recommended)",
      "fGSEA only",
      "ORA only",
      "Skip"),
    default = 1)
  pathway_method <- c("both", "fgsea", "ora", "none")[pathway_idx]
  pathway_enabled <- pathway_idx != 4

  # PPI network analysis
  cat("\n--- Network & Advanced Analysis ---\n")
  ppi_on <- ask_yn("Enable PPI network analysis (STRING)?", TRUE)
  adv_stats_on <- ask_yn("Enable advanced statistics (effect sizes, power analysis)?", TRUE)

  # AI Commentary
  cat("\n--- AI Commentary ---\n")
  cat("Generate AI-powered interpretation for each figure in the report?\n")
  commentary_idx <- ask_choice("Commentary backend:",
    c("None (data-driven fallback only)",
      "Claude Code (uses your subscription, no API key needed)",
      "Claude API (requires ANTHROPIC_API_KEY)",
      "OpenAI API (requires OPENAI_API_KEY)"),
    default = 1)
  commentary_backends <- c("none", "claude-code", "claude", "openai")
  commentary_backend <- commentary_backends[commentary_idx]
  commentary_enabled <- commentary_backend != "none"

  # --- Technical / Instrument Report (Optional) ---
  nna <- function(x, default = "") if (is.null(x) || !nzchar(x)) default else x

  cat("\n--- Technical / Instrument Report (Optional) ---\n")
  cat("Provide a search engine log file (e.g. DIA-NN .log) or facility report\n")
  cat("(DOCX/PDF) to auto-extract instrument, software, and parameter details.\n")
  tech_report_path <- ask("Path to log file or technical report (or 'none')", "none")
  tech_report_fields <- NULL

  if (tolower(tech_report_path) != "none") {
    tech_report_path <- validate_file(tech_report_path, "Technical report / log file")

    if (!is.null(tech_report_path)) {
      claude_available <- nzchar(Sys.which("claude"))

      if (!claude_available) {
        cat("  Claude Code CLI not found — cannot auto-extract fields.\n")
        cat("  You can manually add a 'technical_report:' block to the config later.\n")
      } else {
        # Extract text: for log/txt files read directly, for DOCX/PDF use pandoc
        cat("  Extracting text from file...\n")
        ext <- tolower(tools::file_ext(tech_report_path))

        doc_text <- tryCatch({
          if (ext %in% c("log", "txt", "tsv", "csv")) {
            paste(readLines(tech_report_path, warn = FALSE), collapse = "\n")
          } else {
            paste(system2("pandoc", c("-t", "plain", shQuote(tech_report_path)),
                          stdout = TRUE, stderr = FALSE), collapse = "\n")
          }
        }, error = function(e) NULL)

        if (is.null(doc_text) || nchar(doc_text) < 50) {
          cat("  WARNING: Could not extract sufficient text from the file.\n")
          cat("  You can manually add a 'technical_report:' block to the config later.\n")
        } else {
          cat(sprintf("  Extracted %d characters. Sending to Claude for parsing...\n",
                      nchar(doc_text)))

          json_schema <- paste0(
            '{"type":"object","properties":{',
            '"facility":{"type":"string","description":"Name of proteomics facility if mentioned"},',
            '"sample_prep":{"type":"string","description":"Sample preparation, digestion, labeling if mentioned"},',
            '"ms_acquisition":{"type":"string","description":"Mass spectrometer model (e.g. Exploris 480), acquisition mode (DIA/DDA)"},',
            '"search_engine":{"type":"string","description":"Search engine name and version (e.g. DIA-NN 2.2.0)"},',
            '"search_parameters":{"type":"string","description":"Key search parameters: FDR, peptide length, mods, missed cleavages, FASTA"},',
            '"acknowledgment":{"type":"string","description":"Facility acknowledgment if present"}',
            '},"required":["facility","sample_prep","ms_acquisition","search_engine","search_parameters"]}'
          )

          prompt <- sprintf(
            paste0(
              "Extract structured information from this proteomics search engine log file or facility report. ",
              "Return only the fields described in the schema. ",
              "For each field, use the exact details from the text. ",
              "If a field is not mentioned, use an empty string. ",
              "For ms_acquisition, look for instrument model names in raw file paths (e.g. '480Ex' = Exploris 480). ",
              "For search_engine, extract the software name and version from headers or version lines.\n\n",
              "LOG/REPORT TEXT:\n%s"
            ),
            doc_text
          )

          cmd <- sprintf(
            "claude --print --output-format json --model sonnet --json-schema '%s' --no-session-persistence %s",
            json_schema,
            shQuote(prompt)
          )

          parsed <- tryCatch({
            raw <- system(cmd, intern = TRUE, timeout = 120)
            raw_text <- paste(raw, collapse = "\n")
            p <- jsonlite::fromJSON(raw_text, simplifyVector = FALSE)
            if (!is.null(p$structured_output)) p$structured_output
            else if (!is.null(p$result) && nzchar(p$result)) {
              jsonlite::fromJSON(p$result, simplifyVector = FALSE)
            } else NULL
          }, error = function(e) {
            cat(sprintf("  Claude extraction failed: %s\n", e$message))
            NULL
          })

          if (!is.null(parsed) && (!is.null(parsed$search_engine) || !is.null(parsed$ms_acquisition))) {
            cat("\n  Extracted fields:\n")
            if (nzchar(nna(parsed$facility)))
              cat(sprintf("    Facility:       %s\n", parsed$facility))
            if (nzchar(nna(parsed$sample_prep)))
              cat(sprintf("    Sample prep:    %s\n",
                          substr(nna(parsed$sample_prep), 1, 80)))
            if (nzchar(nna(parsed$ms_acquisition)))
              cat(sprintf("    MS acquisition: %s\n",
                          substr(nna(parsed$ms_acquisition), 1, 80)))
            if (nzchar(nna(parsed$search_engine)))
              cat(sprintf("    Search engine:  %s\n",
                          substr(nna(parsed$search_engine), 1, 80)))
            if (nzchar(nna(parsed$search_parameters)))
              cat(sprintf("    Search params:  %s\n",
                          substr(nna(parsed$search_parameters), 1, 80)))
            if (nzchar(nna(parsed$acknowledgment)))
              cat(sprintf("    Acknowledgment: %s\n",
                          substr(parsed$acknowledgment, 1, 80)))
            cat("  (Truncated for display — full text will be saved in config)\n")

            if (ask_yn("\n  Use these extracted fields?", TRUE)) {
              tech_report_fields <- parsed
            } else {
              cat("  Skipped. You can manually add a 'technical_report:' block later.\n")
            }
          } else {
            cat("  Could not extract structured fields from the file.\n")
            cat("  You can manually add a 'technical_report:' block to the config later.\n")
          }
        }
      }
    }
  }

  # Copy data into project structure
  cat("\n--- Setting Up Project ---\n")

  data_dir <- file.path(project_dir, "data", tolower(project_name), "proteomics")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  protein_dest <- file.path("proteomics", basename(protein_path))
  file.copy(protein_path, file.path(project_dir, "data", tolower(project_name), protein_dest),
            overwrite = TRUE)
  cat(sprintf("  Copied protein matrix -> data/%s/%s\n", tolower(project_name), protein_dest))

  meta_dest <- file.path("proteomics", basename(metadata_path))
  file.copy(metadata_path, file.path(project_dir, "data", tolower(project_name), meta_dest),
            overwrite = TRUE)
  cat(sprintf("  Copied metadata -> data/%s/%s\n", tolower(project_name), meta_dest))

  contrasts_dest <- ""
  if (!is.null(contrasts_path)) {
    contrasts_dest <- file.path("proteomics", basename(contrasts_path))
    file.copy(contrasts_path, file.path(project_dir, "data", tolower(project_name), contrasts_dest),
              overwrite = TRUE)
    cat(sprintf("  Copied contrasts -> data/%s/%s\n", tolower(project_name), contrasts_dest))
  }

  sample_map_dest <- ""
  if (!is.null(sample_map_path)) {
    sample_map_dest <- file.path("proteomics", basename(sample_map_path))
    file.copy(sample_map_path, file.path(project_dir, "data", tolower(project_name), sample_map_dest),
              overwrite = TRUE)
    cat(sprintf("  Copied sample map -> data/%s/%s\n", tolower(project_name), sample_map_dest))
  }

  # Sample map columns
  map_from <- "sample_orig_name"
  map_to <- sample_col
  if (!is.null(sample_map_path)) {
    map_cols <- detect_columns(sample_map_path)
    cat("Sample map columns: ", paste(map_cols, collapse = ", "), "\n")
    map_from <- ask("Map FROM column (original sample names)", map_cols[1])
    map_to <- ask("Map TO column (unified sample IDs)", map_cols[min(2, length(map_cols))])
  }

  # Generate config YAML
  config_name <- paste0("prot_", tolower(gsub("[^a-zA-Z0-9]", "_", project_name)), ".yaml")
  config_path <- file.path(project_dir, "config", config_name)
  dir.create(file.path(project_dir, "config"), showWarnings = FALSE)

  config_yaml <- sprintf('# Auto-generated by run.R wizard
# Project: %s
# Date: %s

project:
  dir: "%s"
  name: "%s"
  analysis_round: "%s"
  analyst: "%s"

paths:
  raw: "data/%s"
  out: "outputs/%s"

modes:
  proteomics:
    engine: "%s"

    files:
      protein: "%s"
      sample_map: "%s"
      metadata: "%s"
      contrasts: "%s"

    id_columns:
      protein_id: "%s"
      sample_col: "%s"
      protein_annot:
        - "Protein.Names"
        - "Genes"
        - "First.Protein.Description"
      map_from: "%s"
      map_to: "%s"

    filtering:
      min_count:
        default: 3

    normalization:
      method: "%s"

    imputation:
      method: "%s"
      no_repetitions: 10
      min_no_passed: 8
      width: %s
      downshift: %s
      dep2_method: "%s"
      dep2_random_seed: 1

    de:
      method: "limma"
      use_adj_for_pass1: true
      p_cutoff: %s
      linear_fc_cutoff: %s

    de_table:
      id_col: "FeatureID"
      pass_any_col: "pass_any_contrast"

    qc_post:
      enabled: true
      de_source: "summary"
      plots:
        volcano: true
        ma: true
      outputs:
        write_de_tables: true

    clustering:
      enabled: %s
      de_source: "any_contrast"
      min_groups: 2
      steps:
        hierarchical:
          enabled: true
          distance: "euclidean"
          linkage: "complete"
        partition:
          enabled: true
          algorithm: "hclust"
          k_max: 20
        binary_patterns:
          enabled: true
          corr_cutoff: 0.8

    effects:
      color: "%s"
      shape: null
      samples: "%s"

    qc:
      run_umap: true
      umap_n_neighbors: 15
      umap_min_dist: 0.1
      outlier_sd_threshold: 3
      min_sample_correlation: 0.8
      top_de_heatmap_sizes: [25, 50, 100]

    annotation:
      organism: "%s"
      id_type: "auto"
      skip_annotation: false
      custom_mapping_file: %s

    pathway:
      enabled: %s
      method: "%s"
      databases: ["GO", "KEGG"]
      gmt_file: %s
      min_size: 10
      max_size: 500

    ppi:
      enabled: %s
      significance_threshold: 0.05
      lfc_threshold: 0
      string_score_threshold: 400
      active_subnetwork: true
      complex_analysis: true

    advanced_stats:
      enabled: %s
      compute_effect_size_ci: true
      bootstrap_n: 1000
      ci_level: 0.95
      run_robust_regression: false
      run_power_analysis: true
      random_effects: null

    commentary:
      enabled: %s
      backend: "%s"
      claude_code_model: "sonnet"
      max_tokens: 1500
      max_retries: 2

params:
  seed: 1
',
    project_name, Sys.Date(), project_dir, project_name, round, analyst,
    tolower(project_name), tolower(project_name),
    engine,
    protein_dest, sample_map_dest, meta_dest, contrasts_dest,
    protein_id_col, sample_col,
    map_from, map_to,
    norm_method,
    imp_method, imp_width, imp_downshift, dep2_method,
    p_cutoff, fc_cutoff,
    tolower(clustering_on),
    group_col, sample_col,
    selected_organism, custom_mapping_file,
    tolower(pathway_enabled), pathway_method, custom_gmt_file,
    tolower(ppi_on),
    tolower(adv_stats_on),
    tolower(commentary_enabled), commentary_backend
  )

  # Inject technical_report block if extracted
  if (!is.null(tech_report_fields)) {
    escape_yaml <- function(s) gsub('"', '\\"', nna(s), fixed = TRUE)
    tech_yaml <- sprintf(
      '    technical_report:\n      facility: "%s"\n      sample_prep: "%s"\n      ms_acquisition: "%s"\n      search_engine: "%s"\n      search_parameters: "%s"\n      acknowledgment: "%s"\n',
      escape_yaml(nna(tech_report_fields$facility)),
      escape_yaml(nna(tech_report_fields$sample_prep)),
      escape_yaml(nna(tech_report_fields$ms_acquisition)),
      escape_yaml(nna(tech_report_fields$search_engine)),
      escape_yaml(nna(tech_report_fields$search_parameters)),
      escape_yaml(nna(tech_report_fields$acknowledgment))
    )
    config_yaml <- sub("    commentary:", paste0(tech_yaml, "\n    commentary:"), config_yaml)
  }

  writeLines(config_yaml, config_path)
  cat(sprintf("\n  Config saved: %s\n", config_path))

  config_path
}

# --- Pipeline Runner ----------------------------------------------------------

run_pipeline <- function(config_path, fresh = FALSE) {
  config_path <- normalizePath(config_path, mustWork = TRUE)
  cat(sprintf("\n  Config: %s\n", config_path))

  # Set environment variable for _targets.R
  Sys.setenv(MULTIOMICS_CONFIG = config_path)

  if (fresh) {
    cat("  Clearing targets cache (fresh run)...\n")
    targets::tar_invalidate(everything())
  }

  # Quick pre-flight: verify input files exist before starting pipeline
  cfg <- yaml::read_yaml(config_path)
  if (!is.null(cfg$modes$rna)) {
    rna <- cfg$modes$rna
    base_dir <- file.path(cfg$project$dir, cfg$paths$raw)
    missing_files <- character(0)
    for (fname in c("counts", "metadata", "contrasts")) {
      fpath <- rna$files[[fname]]
      if (!is.null(fpath) && nzchar(fpath) && !file.exists(file.path(base_dir, fpath))) {
        missing_files <- c(missing_files, sprintf("  [x] %s: %s", fname, fpath))
      }
    }
    if (length(missing_files) > 0) {
      cat("\n  PRE-FLIGHT: Missing input files:\n")
      cat(paste(missing_files, collapse = "\n"), "\n")
      if (interactive()) {
        proceed <- ask_yn("Continue anyway?", FALSE)
        if (!proceed) return(invisible(NULL))
      } else {
        stop("Pre-flight check failed: missing input files")
      }
    }
  }

  # Pre-flight: proteomics input files
  if (!is.null(cfg$modes$proteomics)) {
    prot <- cfg$modes$proteomics
    base_dir <- file.path(cfg$project$dir, cfg$paths$raw)
    missing_prot <- character(0)
    for (fname in c("protein", "metadata", "contrasts")) {
      fpath <- prot$files[[fname]]
      if (!is.null(fpath) && nzchar(fpath) && !file.exists(file.path(base_dir, fpath))) {
        missing_prot <- c(missing_prot, sprintf("  [x] proteomics/%s: %s", fname, fpath))
      }
    }
    if (length(missing_prot) > 0) {
      cat("\n  PRE-FLIGHT: Missing proteomics input files:\n")
      cat(paste(missing_prot, collapse = "\n"), "\n")
      if (interactive()) {
        proceed <- ask_yn("Continue anyway?", FALSE)
        if (!proceed) return(invisible(NULL))
      } else {
        stop("Pre-flight check failed: missing proteomics input files")
      }
    }
  }

  cat("  Starting pipeline...\n\n")
  targets::tar_make()

  # Find the reports
  cfg <- yaml::read_yaml(config_path)
  run_dir_name <- sprintf("Results_%s_%s",
                           cfg$project$name, cfg$project$analysis_round)
  run_dir_path <- file.path(cfg$project$dir, cfg$paths$out, run_dir_name)

  report_path <- file.path(run_dir_path, "report_rnaseq.html")
  if (file.exists(report_path)) {
    cat(sprintf("\n  RNA-seq Report: %s\n", report_path))
    cat("  Open in browser to view results.\n")
  }

  summary_path <- file.path(run_dir_path, "pipeline_summary.html")
  if (file.exists(summary_path)) {
    cat(sprintf("  Pipeline Summary: %s\n", summary_path))
  }

  prot_report_path <- file.path(run_dir_path, "report_proteomics.html")
  if (file.exists(prot_report_path)) {
    cat(sprintf("\n  Proteomics Report: %s\n", prot_report_path))
    cat("  Open in browser to view results.\n")
  }

  cat("\n  Pipeline complete.\n")
}

# --- Main Entry Point ---------------------------------------------------------

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  project_dir <- getwd()

  if ("--new" %in% args) {
    # Interactive wizard
    config_path <- wizard_new_project(project_dir)
    cat("\n")
    run_now <- ask_yn("Run the pipeline now?", TRUE)
    if (run_now) {
      run_pipeline(config_path, fresh = TRUE)
    } else {
      cat(sprintf("\nTo run later:\n  Rscript run.R --config %s\n\n", config_path))
    }

  } else if ("--config" %in% args) {
    idx <- which(args == "--config")
    if (idx < length(args)) {
      config_path <- args[idx + 1]
    } else {
      stop("--config requires a path argument")
    }
    fresh <- "--fresh" %in% args
    run_pipeline(config_path, fresh = fresh)

  } else if ("--fresh" %in% args) {
    # Re-run with default/last config, fresh
    config_path <- Sys.getenv("MULTIOMICS_CONFIG", unset = "")
    if (config_path == "") {
      # Find the most recent config
      configs <- list.files(file.path(project_dir, "config"),
                             pattern = "\\.yaml$", full.names = TRUE)
      configs <- configs[!grepl("templates", configs)]
      if (length(configs) == 0) stop("No config files found. Run with --new first.")
      config_path <- configs[which.max(file.mtime(configs))]
      cat(sprintf("  Using most recent config: %s\n", basename(config_path)))
    }
    run_pipeline(config_path, fresh = TRUE)

  } else if (length(args) == 0) {
    # No args: if running interactively (RStudio), show menu
    if (interactive()) {
      cli_header()
      cat("What would you like to do?\n")
      choice <- ask_choice("",
        c("Set up a new project (interactive wizard)",
          "Run with an existing config file",
          "Re-run last project (with caching)",
          "Re-run last project (fresh, no cache)"),
        default = 1)

      if (choice == 1) {
        config_path <- wizard_new_project(project_dir)
        run_now <- ask_yn("\nRun the pipeline now?", TRUE)
        if (run_now) run_pipeline(config_path, fresh = TRUE)

      } else if (choice == 2) {
        configs <- list.files(file.path(project_dir, "config"),
                               pattern = "\\.yaml$", full.names = TRUE)
        configs <- configs[!grepl("templates", configs)]
        if (length(configs) == 0) {
          cat("No config files found. Run with option 1 first.\n")
          return(invisible())
        }
        cat("\nAvailable configs:\n")
        for (i in seq_along(configs)) {
          cat(sprintf("  %d) %s\n", i, basename(configs[i])))
        }
        idx <- as.integer(ask("Choose config", "1"))
        run_pipeline(configs[idx])

      } else if (choice == 3) {
        configs <- list.files(file.path(project_dir, "config"),
                               pattern = "\\.yaml$", full.names = TRUE)
        configs <- configs[!grepl("templates", configs)]
        config_path <- configs[which.max(file.mtime(configs))]
        cat(sprintf("  Using: %s\n", basename(config_path)))
        run_pipeline(config_path, fresh = FALSE)

      } else if (choice == 4) {
        configs <- list.files(file.path(project_dir, "config"),
                               pattern = "\\.yaml$", full.names = TRUE)
        configs <- configs[!grepl("templates", configs)]
        config_path <- configs[which.max(file.mtime(configs))]
        cat(sprintf("  Using: %s\n", basename(config_path)))
        run_pipeline(config_path, fresh = TRUE)
      }

    } else {
      # Non-interactive, no args: re-run with caching
      config_path <- Sys.getenv("MULTIOMICS_CONFIG", unset = "")
      if (config_path == "") {
        configs <- list.files(file.path(project_dir, "config"),
                               pattern = "\\.yaml$", full.names = TRUE)
        configs <- configs[!grepl("templates", configs)]
        if (length(configs) == 0) stop("No config files found. Run with --new first.")
        config_path <- configs[which.max(file.mtime(configs))]
      }
      run_pipeline(config_path, fresh = FALSE)
    }
  }
}

# Run
main()
