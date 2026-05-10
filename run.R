#!/usr/bin/env Rscript
# =============================================================================
# multiomics-core Pipeline Runner
#
# Usage:
#   Rscript run.R --wizard           # Open HTML wizard in browser (GUI)
#   Rscript run.R --new              # Interactive wizard for new project (CLI)
#   Rscript run.R --config path.yaml # Run with existing config
#   Rscript run.R                    # Re-run last config (with caching)
#   Rscript run.R --fresh            # Full re-run (invalidate cache)
#   Rscript run.R --help             # Show this help message
#
# In RStudio: source("run.R") runs interactively
# =============================================================================

# --- Fast-start: if launched with --vanilla, add renv library path manually ---
# Also propagate the library path + disable renv's autoloader in any child
# R process the pipeline spawns (callr, parallel::makeCluster, etc.). Without
# this, each subprocess re-runs renv::load() and spends ~27 s scanning the
# project for dependencies — on Windows this caused pathway analysis to
# balloon from ~5 min to ~12 min.
if (!"renv" %in% loadedNamespaces()) {
  lib_dirs <- list.dirs(file.path(getwd(), "renv", "library"), recursive = TRUE)
  for (d in lib_dirs) {
    if (file.exists(file.path(d, "jsonlite"))) {
      .libPaths(c(d, .libPaths()))
      break
    }
  }
  # Export current library paths so grandchild R processes find packages
  # without having to activate renv themselves.
  Sys.setenv(
    R_LIBS                           = paste(.libPaths(), collapse = .Platform$path.sep),
    RENV_CONFIG_AUTOLOADER_ENABLED   = "FALSE"
  )
}

# --- Dependency check --------------------------------------------------------

check_dependencies <- function(mode = "cli") {
  core_pkgs <- c("yaml", "jsonlite", "targets")
  wizard_pkgs <- c("httpuv", "later")
  pkgs <- if (mode == "wizard") c(core_pkgs, wizard_pkgs) else core_pkgs

  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    cat("\n  Missing required packages: ", paste(missing, collapse = ", "), "\n")
    cat("  Install with:  renv::restore()  or  install.packages(c(",
        paste0('"', missing, '"', collapse = ", "), "))\n\n")
    stop("Missing dependencies. See above.", call. = FALSE)
  }
}

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

#' Read unique levels from a column in a metadata file
detect_group_levels <- function(metadata_path, group_col) {
  sep <- detect_separator(metadata_path)
  meta <- read.table(metadata_path, sep = sep, header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)
  if (!group_col %in% colnames(meta)) {
    cat(sprintf("  WARNING: Column '%s' not found in metadata.\n", group_col))
    return(NULL)
  }
  levels <- unique(as.character(meta[[group_col]]))
  levels <- levels[!is.na(levels) & nzchar(levels)]
  sort(levels)
}

#' Generate all pairwise contrasts from a set of group levels
make_pairwise_contrasts <- function(levels, factor_col) {
  pairs <- combn(levels, 2)
  data.frame(
    Contrast_name = apply(pairs, 2, function(p) paste0(p[1], "_vs_", p[2])),
    Factor        = factor_col,
    Numerator     = pairs[1, ],
    Denominator   = pairs[2, ],
    stringsAsFactors = FALSE
  )
}

#' Interactive contrast builder — called when no contrasts file is provided
#' Returns a contrasts data.frame, or NULL if skipped
wizard_build_contrasts <- function(metadata_path, group_col) {
  cat("\n--- Contrast Builder ---\n")
  cat(sprintf("No contrasts file provided. Detecting groups from '%s' column...\n", group_col))

  levels <- detect_group_levels(metadata_path, group_col)
  if (is.null(levels) || length(levels) < 2) {
    cat("  Could not detect enough groups (need at least 2). Skipping contrast generation.\n")
    return(NULL)
  }

  cat(sprintf("  Detected %d groups: %s\n", length(levels), paste(levels, collapse = ", ")))

  all_pairs <- make_pairwise_contrasts(levels, group_col)
  n_pairs <- nrow(all_pairs)

  action_idx <- ask_choice("\nHow would you like to define contrasts?",
    c(sprintf("All pairwise comparisons (%d contrasts) (recommended)", n_pairs),
      "Pick specific comparisons from the list",
      "Build custom contrasts (choose numerator & denominator)",
      "Skip — no contrasts (DE will auto-generate all pairwise at runtime)",
      "Load from a contrasts CSV file"),
    default = 1)

  if (action_idx == 4) {
    cat("  Skipping. Pipeline will auto-generate all pairwise contrasts at runtime.\n")
    return(NULL)
  }

  if (action_idx == 5) {
    contrasts_path <- ask("Path to contrasts CSV file")
    contrasts_path <- validate_file(contrasts_path, "Contrasts file")
    if (is.null(contrasts_path)) {
      cat("  File not found. Falling back to all pairwise contrasts.\n")
      contrasts_df <- all_pairs
    } else {
      sep <- detect_separator(contrasts_path)
      contrasts_df <- read.table(contrasts_path, sep = sep, header = TRUE,
                                  stringsAsFactors = FALSE, check.names = FALSE)
      # Normalize column names: expect Numerator/Denominator or similar
      cn <- tolower(colnames(contrasts_df))
      if (!"Numerator" %in% colnames(contrasts_df)) {
        # Try common alternatives
        if ("numerator" %in% cn) colnames(contrasts_df)[cn == "numerator"] <- "Numerator"
      }
      if (!"Denominator" %in% colnames(contrasts_df)) {
        if ("denominator" %in% cn) colnames(contrasts_df)[cn == "denominator"] <- "Denominator"
      }
      if (!"Factor" %in% colnames(contrasts_df)) {
        if ("factor" %in% cn) {
          colnames(contrasts_df)[cn == "factor"] <- "Factor"
        } else {
          contrasts_df$Factor <- group_col
        }
      }
      cat(sprintf("  Loaded %d contrast(s) from file:\n", nrow(contrasts_df)))
      for (i in seq_len(nrow(contrasts_df))) {
        cat(sprintf("    %d) %s vs %s\n", i, contrasts_df$Numerator[i], contrasts_df$Denominator[i]))
      }
    }
    cat(sprintf("  %d contrast(s) defined.\n", nrow(contrasts_df)))
    return(contrasts_df)
  }

  contrasts_df <- NULL

  if (action_idx == 1) {
    # Auto-generate all pairwise
    contrasts_df <- all_pairs
    cat("  Generated contrasts:\n")
    for (i in seq_len(nrow(contrasts_df))) {
      cat(sprintf("    %d) %s vs %s\n", i, contrasts_df$Numerator[i], contrasts_df$Denominator[i]))
    }

  } else if (action_idx == 2) {
    # Let user pick from all pairwise options
    cat("\nAvailable pairwise comparisons:\n")
    for (i in seq_len(n_pairs)) {
      cat(sprintf("  %d) %s vs %s\n", i, all_pairs$Numerator[i], all_pairs$Denominator[i]))
    }
    selection <- ask("Enter contrast numbers to include (comma-separated, e.g. 1,3)", "all")
    if (tolower(selection) == "all") {
      contrasts_df <- all_pairs
    } else {
      idx <- suppressWarnings(as.integer(trimws(unlist(strsplit(selection, ",")))))
      idx <- idx[!is.na(idx) & idx >= 1 & idx <= n_pairs]
      if (length(idx) == 0) {
        cat("  No valid selections. Using all contrasts.\n")
        contrasts_df <- all_pairs
      } else {
        contrasts_df <- all_pairs[idx, , drop = FALSE]
      }
    }
    cat(sprintf("  Selected %d contrast(s).\n", nrow(contrasts_df)))

  } else if (action_idx == 3) {
    # Build one-at-a-time by selecting numerator & denominator from menus
    cat("\nBuild contrasts by choosing numerator and denominator.\n")
    cat("Groups:\n")
    for (i in seq_along(levels)) cat(sprintf("  %d) %s\n", i, levels[i]))
    rows <- list()
    repeat {
      cat(sprintf("\n  Contrast %d:\n", length(rows) + 1))
      num_idx <- as.integer(ask("    Numerator  (number, or 'done' to finish)"))
      if (is.na(num_idx)) break
      den_idx <- as.integer(ask("    Denominator (number)"))
      if (is.na(den_idx)) break
      if (num_idx < 1 || num_idx > length(levels)) { cat("    Invalid numerator.\n"); next }
      if (den_idx < 1 || den_idx > length(levels)) { cat("    Invalid denominator.\n"); next }
      if (num_idx == den_idx) { cat("    Numerator and denominator must differ.\n"); next }
      num <- levels[num_idx]; den <- levels[den_idx]
      cat(sprintf("    -> %s vs %s\n", num, den))
      rows[[length(rows) + 1]] <- data.frame(
        Contrast_name = paste0(num, "_vs_", den), Factor = group_col,
        Numerator = num, Denominator = den, stringsAsFactors = FALSE
      )
      if (!ask_yn("    Add another contrast?", TRUE)) break
    }
    if (length(rows) == 0) {
      cat("  No contrasts defined. Skipping.\n")
      return(NULL)
    }
    contrasts_df <- do.call(rbind, rows)
  }

  cat(sprintf("  %d contrast(s) defined.\n", nrow(contrasts_df)))
  contrasts_df
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
      "Proteomics",
      "Metabolomics",
      "Lipidomics",
      "Multiomics (integration pipeline)",
      "Combination (choose layers)"),
    default = 1)
  analysis_mode <- c("rna", "proteomics", "metabolomics", "lipidomics",
                      "multiomics", "combination")[mode_idx]

  # For combination mode, let user pick which layers to include
  if (analysis_mode == "combination") {
    cat("\nSelect layers to include:\n")
    layers <- c("rna", "proteomics", "metabolomics", "lipidomics")
    layer_labels <- c("RNA-seq", "Proteomics", "Metabolomics", "Lipidomics")
    active_layers <- character(0)
    for (i in seq_along(layers)) {
      ans <- ask(sprintf("Include %s? (y/n)", layer_labels[i]), "n")
      if (tolower(substr(ans, 1, 1)) == "y") {
        active_layers <- c(active_layers, layers[i])
      }
    }
    if (length(active_layers) < 2) {
      cat("Warning: fewer than 2 layers selected. Multi-omics integration requires at least 2.\n")
      if (length(active_layers) == 1) {
        analysis_mode <- active_layers[1]
      } else {
        stop("No layers selected. Aborting.")
      }
    }
  }

  if (analysis_mode == "multiomics") {
    wizard_multiomics(project_dir, project_name, analyst, round)
  } else if (analysis_mode == "rna") {
    wizard_rna(project_dir, project_name, analyst, round)
  } else if (analysis_mode == "proteomics") {
    wizard_proteomics(project_dir, project_name, analyst, round)
  } else if (analysis_mode == "metabolomics") {
    wizard_metabolomics(project_dir, project_name, analyst, round)
  } else if (analysis_mode == "lipidomics") {
    wizard_lipidomics(project_dir, project_name, analyst, round)
  } else {
    # Multi-mode: run each wizard, merge configs
    configs <- list()

    if ("rna" %in% active_layers) {
      cat("\n=== RNA-seq Setup ===\n")
      configs$rna <- wizard_rna(project_dir, project_name, analyst, round)
    }

    if ("proteomics" %in% active_layers) {
      cat("\n=== Proteomics Setup ===\n")
      configs$prot <- wizard_proteomics(project_dir, project_name, analyst, round)
    }

    if ("metabolomics" %in% active_layers) {
      cat("\n=== Metabolomics Setup ===\n")
      configs$metab <- wizard_metabolomics(project_dir, project_name, analyst, round)
    }

    if ("lipidomics" %in% active_layers) {
      cat("\n=== Lipidomics Setup ===\n")
      configs$lipid <- wizard_lipidomics(project_dir, project_name, analyst, round)
    }

    # Merge configs: start with the first one, add modes from others
    first_path <- configs[[1]]
    merged <- yaml::read_yaml(first_path)

    for (cfg_path in configs[-1]) {
      other_cfg <- yaml::read_yaml(cfg_path)
      for (mode_name in names(other_cfg$modes)) {
        merged$modes[[mode_name]] <- other_cfg$modes[[mode_name]]
      }
    }

    merged_name <- paste0("multi_", tolower(gsub("[^a-zA-Z0-9]", "_", project_name)), ".yaml")
    merged_path <- file.path(project_dir, "config", merged_name)
    yaml::write_yaml(merged, merged_path)

    cat(sprintf("\n  Merged config saved: %s\n", merged_path))
    merged_path
  }
}

# --- RNA-seq Wizard -----------------------------------------------------------

wizard_rna <- function(project_dir, project_name, analyst, round) {

  # Data directory
  cat("\n--- [RNA] Data Files ---\n")
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
  cat("\n--- [RNA] Column Detection ---\n")
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

  # --- Sample Exclusion ---
  excluded_samples <- character(0)
  if (!is.null(metadata_path)) {
    sep <- detect_separator(metadata_path)
    meta_preview <- read.table(metadata_path, sep = sep, header = TRUE,
                               stringsAsFactors = FALSE, check.names = FALSE)
    all_samples <- as.character(meta_preview[[sample_col]])

    cat("\n--- [RNA] Sample Exclusion (optional) ---\n")
    cat("Samples found in metadata:\n")
    for (i in seq_along(all_samples)) {
      cat(sprintf("  %d) %s\n", i, all_samples[i]))
    }
    repeat {
      remaining <- setdiff(all_samples, excluded_samples)
      if (length(remaining) == 0) break
      cat("\nCurrently excluded: ",
          if (length(excluded_samples) == 0) "(none)" else paste(excluded_samples, collapse = ", "), "\n")
      ans <- prompt_read("Enter sample number to exclude (or press Enter to continue): ")
      if (ans == "") break
      idx <- suppressWarnings(as.integer(ans))
      if (!is.na(idx) && idx >= 1 && idx <= length(all_samples)) {
        s <- all_samples[idx]
        if (s %in% excluded_samples) {
          cat(sprintf("  '%s' is already excluded.\n", s))
        } else {
          excluded_samples <- c(excluded_samples, s)
          cat(sprintf("  Excluded: %s (%d remaining)\n", s, length(all_samples) - length(excluded_samples)))
        }
      } else {
        cat("  Invalid selection.\n")
      }
    }
    if (length(excluded_samples) > 0) {
      cat(sprintf("Final excluded samples: %s\n", paste(excluded_samples, collapse = ", ")))
    }
  }

  # Build contrasts if none provided
  generated_contrasts_df <- NULL
  if (is.null(contrasts_path) && !is.null(metadata_path)) {
    generated_contrasts_df <- wizard_build_contrasts(metadata_path, group_col)
  }

  # Methods
  cat("\n--- [RNA] Analysis Methods ---\n")

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
  cat("\n--- [RNA] Organism & Annotation ---\n")
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
  cat("\n--- [RNA] AI Commentary ---\n")
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

  # --- PowerPoint ---
  generate_pptx <- ask_yn("Generate a PowerPoint summary presentation?", TRUE)

  # --- Technical Report (Optional) ---
  nna <- function(x, default = "") if (is.null(x) || !nzchar(x)) default else x

  cat("\n--- [RNA] Technical Report (Optional) ---\n")
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
        # Extract text: PDF via pdftotext, everything else via pandoc
        cat("  Extracting text from report...\n")
        ext <- tolower(tools::file_ext(tech_report_path))
        doc_text <- tryCatch({
          if (ext == "pdf") {
            paste(system2("pdftotext", c("-layout", shQuote(tech_report_path), "-"),
                          stdout = TRUE, stderr = FALSE), collapse = "\n")
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
            "claude --print --output-format json --model sonnet --json-schema %s --no-session-persistence %s",
            shQuote(json_schema),
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
  cat("\n--- [RNA] Setting Up Project ---\n")

  # Create data directory
  data_dir <- file.path(project_dir, "data", tolower(project_name), "rna")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  # Copy input files
  if (is.null(counts_path)) {
    cat("  ERROR: No valid counts file provided. Cannot set up project.\n")
    return(invisible(NULL))
  }
  counts_dest <- file.path("rna", basename(counts_path))
  file.copy(counts_path, file.path(project_dir, "data", tolower(project_name), counts_dest),
            overwrite = TRUE)
  cat(sprintf("  Copied counts -> data/%s/%s\n", tolower(project_name), counts_dest))

  if (is.null(metadata_path)) {
    cat("  ERROR: No valid metadata file provided. Cannot set up project.\n")
    return(invisible(NULL))
  }
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
  } else if (!is.null(generated_contrasts_df)) {
    contrasts_dest <- file.path("rna", "contrasts_auto.csv")
    write.csv(generated_contrasts_df,
              file.path(project_dir, "data", tolower(project_name), contrasts_dest),
              row.names = FALSE)
    cat(sprintf("  Generated contrasts -> data/%s/%s\n", tolower(project_name), contrasts_dest))
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

  # Build exclude_samples YAML fragment
  exclude_yaml <- ""
  if (length(excluded_samples) > 0) {
    exclude_yaml <- paste0("\n    exclude_samples:\n",
        paste0("      - \"", excluded_samples, "\"", collapse = "\n"), "\n")
  } else {
    exclude_yaml <- "\n    exclude_samples: null\n"
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
%s
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

    outputs:
      generate_pptx: %s

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
    exclude_yaml,
    de_method, deseq_mode, p_cutoff, fc_cutoff,
    selected_organism, custom_mapping_file,
    tolower(pathway_enabled), pathway_method, custom_gmt_file,
    tolower(commentary_enabled), commentary_backend,
    tolower(generate_pptx),
    tolower(clustering_on), group_col,
    tolower(batch_corr_on),
    tolower(deconv_on),
    group_col, shape_col, sample_col
  )

  # Inject technical_report block if extracted
  if (!is.null(tech_report_fields)) {
    escape_yaml <- function(s) gsub('"', '\\"', nna(s), fixed = TRUE)
    tech_yaml <- paste0(
      '    technical_report:\n',
      '      facility: "', escape_yaml(tech_report_fields$facility), '"\n',
      '      library_prep: "', escape_yaml(tech_report_fields$library_prep), '"\n',
      '      sequencing: "', escape_yaml(tech_report_fields$sequencing), '"\n',
      '      upstream_bioinformatics: "', escape_yaml(tech_report_fields$upstream_bioinformatics), '"\n',
      '      upstream_de_note: "', escape_yaml(nna(tech_report_fields$upstream_de_note)), '"\n',
      '      acknowledgment: "', escape_yaml(nna(tech_report_fields$acknowledgment)), '"\n'
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
  cat("\n--- [PROT] Data Files ---\n")
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
  cat("\n--- [PROT] Column Detection ---\n")
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

  # --- Sample Exclusion ---
  excluded_samples <- character(0)
  if (!is.null(metadata_path)) {
    sep <- detect_separator(metadata_path)
    meta_preview <- read.table(metadata_path, sep = sep, header = TRUE,
                               stringsAsFactors = FALSE, check.names = FALSE)
    all_samples <- as.character(meta_preview[[sample_col]])

    cat("\n--- [PROT] Sample Exclusion (optional) ---\n")
    cat("Samples found in metadata:\n")
    for (i in seq_along(all_samples)) {
      cat(sprintf("  %d) %s\n", i, all_samples[i]))
    }
    repeat {
      remaining <- setdiff(all_samples, excluded_samples)
      if (length(remaining) == 0) break
      cat("\nCurrently excluded: ",
          if (length(excluded_samples) == 0) "(none)" else paste(excluded_samples, collapse = ", "), "\n")
      ans <- prompt_read("Enter sample number to exclude (or press Enter to continue): ")
      if (ans == "") break
      idx <- suppressWarnings(as.integer(ans))
      if (!is.na(idx) && idx >= 1 && idx <= length(all_samples)) {
        s <- all_samples[idx]
        if (s %in% excluded_samples) {
          cat(sprintf("  '%s' is already excluded.\n", s))
        } else {
          excluded_samples <- c(excluded_samples, s)
          cat(sprintf("  Excluded: %s (%d remaining)\n", s, length(all_samples) - length(excluded_samples)))
        }
      } else {
        cat("  Invalid selection.\n")
      }
    }
    if (length(excluded_samples) > 0) {
      cat(sprintf("Final excluded samples: %s\n", paste(excluded_samples, collapse = ", ")))
    }
  }

  # Build contrasts if none provided
  generated_contrasts_df <- NULL
  if (is.null(contrasts_path) && !is.null(metadata_path)) {
    generated_contrasts_df <- wizard_build_contrasts(metadata_path, group_col)
  }

  # Engine
  cat("\n--- [PROT] Settings ---\n")
  engine_idx <- ask_choice("Search engine used:",
    c("DIA-NN (default)", "MaxQuant", "FragPipe", "Other"),
    default = 1)
  engine <- c("DIANN", "MaxQuant", "FragPipe", "Other")[engine_idx]

  # Input data scale
  scale_idx <- ask_choice("Input data scale (are intensities already log2-transformed?):",
    c("Linear (raw intensities — will apply log2)",
      "Log2 (already log2-transformed, e.g. DIA-NN default output)"),
    default = 1)
  scale_in <- c("linear", "log2")[scale_idx]

  # Normalization
  norm_idx <- ask_choice("Normalization method:",
    c("None (data already normalized, e.g. DIA-NN output)",
      "Median centering",
      "VSN (variance stabilizing — requires Bioconductor 'vsn' package)"),
    default = 1)
  norm_method <- c("none", "median", "vsn")[norm_idx]

  # Imputation
  imp_idx <- ask_choice("Imputation method:",
    c("Perseus-style — downshifted normal (recommended for MNAR data)",
      "DEP2 / MinDet — deterministic minimum (experimental)",
      "QRILC — truncated normal, single imputation (DEP workflow)",
      "MinVal — floor of minimum observed value (deterministic LOD)",
      "None (complete cases only)"),
    default = 1)
  imp_method <- c("perseus", "dep2", "qrilc", "minval", "none")[imp_idx]

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

  multi_imp <- TRUE
  if (imp_method != "none") {
    multi_imp <- ask_yn("Use multiple imputation with consensus filter (10 rounds, 8/10 must pass)?", TRUE)
  }

  # Filtering options
  cat("\n--- [PROT] Filtering ---\n")
  remove_contaminants <- ask_yn("Remove contaminant proteins (e.g. cRAP)?", TRUE)
  min_valid <- ask("Minimum valid values per group", "3")
  min_groups <- ask("Require passing in at least N groups", "1")

  # Peptide file (stub)
  cat("\n--- [PROT] Optional Files ---\n")
  peptide_path <- ask("Path to peptide-level data file (optional, or 'none')", "none")
  peptide_dest <- ""
  if (tolower(peptide_path) != "none") {
    peptide_path <- validate_file(peptide_path, "Peptide data file")
    if (!is.null(peptide_path)) peptide_dest <- file.path("proteomics", basename(peptide_path))
  }

  # DE method & cutoffs
  cat("\n--- [PROT] Differential Expression ---\n")
  de_method_idx <- ask_choice("DE method:",
    c("Limma (recommended)",
      "t-test (equal variance)",
      "Welch test (unequal variance)",
      "ANOVA + Tukey HSD (multi-group)"),
    default = 1)
  de_method <- c("limma", "ttest", "welch", "anova")[de_method_idx]

  paired <- FALSE
  pairing_col <- "null"
  if (de_method %in% c("ttest", "welch")) {
    paired <- ask_yn("Use paired test?", FALSE)
    if (paired && !is.null(metadata_path)) {
      meta_cols_all <- detect_columns(metadata_path)
      cat("  Metadata columns: ", paste(meta_cols_all, collapse = ", "), "\n")
      pairing_col_ans <- ask("Pairing column (e.g. Subject)", "null")
      if (tolower(pairing_col_ans) != "null") {
        pairing_col <- paste0('"', pairing_col_ans, '"')
      }
    }
  }

  p_cutoff <- as.numeric(ask("P-value cutoff", "0.05"))
  use_adj <- ask_yn("Use adjusted p-value (BH) for significance?", TRUE)

  p_adjust_idx <- ask_choice("Multiple testing correction method:",
    c("BH / Benjamini-Hochberg (recommended)",
      "Bonferroni",
      "Holm",
      "BY / Benjamini-Yekutieli"),
    default = 1)
  p_adjust_method <- c("BH", "bonferroni", "holm", "BY")[p_adjust_idx]

  fdrtool_corr <- ask_yn("Apply fdrtool empirical null correction (DEP-style, more conservative)?", FALSE)

  fc_cutoff <- as.numeric(ask("Linear fold-change cutoff (1 = no FC filter)", "1.5"))

  clustering_on <- ask_yn("Enable clustering?", FALSE)

  # Organism & Annotation
  cat("\n--- [PROT] Organism & Annotation ---\n")
  org_idx <- ask_choice("Which organism does your data come from?",
    c("Human (Homo sapiens)",
      "Mouse (Mus musculus)",
      "Rat (Rattus norvegicus)",
      "Zebrafish (Danio rerio)",
      "C. elegans (Caenorhabditis elegans)",
      "Other (type name)"),
    default = 1)
  organism_names <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "Caenorhabditis elegans", "other")
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
  cat("\n--- [PROT] Pathway Enrichment ---\n")
  pathway_idx <- ask_choice("Pathway analysis:",
    c("Both fGSEA + ORA (recommended)",
      "fGSEA only",
      "ORA only",
      "Skip"),
    default = 1)
  pathway_method <- c("both", "fgsea", "ora", "none")[pathway_idx]
  pathway_enabled <- pathway_idx != 4

  gsea_ranking <- "stat"
  pathway_volcano <- TRUE
  if (pathway_enabled) {
    ranking_idx <- ask_choice("GSEA ranking method:",
      c("Signed statistic — sign(log2FC) * -log10(p) (recommended)",
        "Absolute fold-change — |log2FC|",
        "Log2 fold-change — signed, directional"),
      default = 1)
    gsea_ranking <- c("stat", "abs_lfc", "lfc")[ranking_idx]

    pathway_volcano <- ask_yn("Generate pathway-colored volcano data?", TRUE)
  }

  # PPI network analysis
  cat("\n--- [PROT] Network & Advanced Analysis ---\n")
  ppi_on <- ask_yn("Enable PPI network analysis (STRING)?", TRUE)
  adv_stats_on <- ask_yn("Enable advanced statistics (effect sizes, power analysis)?", TRUE)

  # AI Commentary
  cat("\n--- [PROT] AI Commentary ---\n")
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

  # --- PowerPoint ---
  generate_pptx <- ask_yn("Generate a PowerPoint summary presentation?", TRUE)

  # --- Report Customization (Optional) ---
  cat("\n--- [PROT] Report Customization (Optional) ---\n")
  cat("These settings control the HTML report appearance. All are optional.\n")
  cat("You can also edit them later in the config YAML under modes.proteomics.report.\n\n")

  report_project_id <- ask("Project ID for report title (or leave empty)", "")
  report_intro <- ask("Introduction text for the report (or leave empty)", "")
  report_disclaimer <- ask("Legal disclaimer text (or leave empty)", "")

  show_distance_heatmap <- ask_yn("Show distance heatmap in QC section?", FALSE)
  show_venn <- ask_yn("Show identified-proteins Venn diagram?", FALSE)

  # --- Technical / Instrument Report (Optional) ---
  nna <- function(x, default = "") if (is.null(x) || !nzchar(x)) default else x

  cat("\n--- [PROT] Technical / Instrument Report (Optional) ---\n")
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
        # Extract text: plain text via readLines, PDF via pdftotext, DOCX via pandoc
        cat("  Extracting text from file...\n")
        ext <- tolower(tools::file_ext(tech_report_path))

        doc_text <- tryCatch({
          if (ext %in% c("log", "txt", "tsv", "csv")) {
            paste(readLines(tech_report_path, warn = FALSE), collapse = "\n")
          } else if (ext == "pdf") {
            paste(system2("pdftotext", c("-layout", shQuote(tech_report_path), "-"),
                          stdout = TRUE, stderr = FALSE), collapse = "\n")
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
            "claude --print --output-format json --model sonnet --json-schema %s --no-session-persistence %s",
            shQuote(json_schema),
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
  cat("\n--- [PROT] Setting Up Project ---\n")

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
  } else if (!is.null(generated_contrasts_df)) {
    contrasts_dest <- file.path("proteomics", "contrasts_auto.csv")
    write.csv(generated_contrasts_df,
              file.path(project_dir, "data", tolower(project_name), contrasts_dest),
              row.names = FALSE)
    cat(sprintf("  Generated contrasts -> data/%s/%s\n", tolower(project_name), contrasts_dest))
  }

  sample_map_dest <- ""
  if (!is.null(sample_map_path)) {
    sample_map_dest <- file.path("proteomics", basename(sample_map_path))
    file.copy(sample_map_path, file.path(project_dir, "data", tolower(project_name), sample_map_dest),
              overwrite = TRUE)
    cat(sprintf("  Copied sample map -> data/%s/%s\n", tolower(project_name), sample_map_dest))
  }

  # Copy peptide file if provided
  if (nzchar(peptide_dest) && !is.null(peptide_path)) {
    file.copy(peptide_path, file.path(project_dir, "data", tolower(project_name), peptide_dest),
              overwrite = TRUE)
    cat(sprintf("  Copied peptides -> data/%s/%s\n", tolower(project_name), peptide_dest))
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

  # Build exclude_samples YAML fragment
  exclude_yaml <- ""
  if (length(excluded_samples) > 0) {
    exclude_yaml <- paste0("\n    exclude_samples:\n",
        paste0("      - \"", excluded_samples, "\"", collapse = "\n"), "\n")
  } else {
    exclude_yaml <- "\n    exclude_samples: null\n"
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
  proteomics:
    engine: "%s"
    scale_in: "%s"

    files:
      protein: "%s"
      sample_map: "%s"
      metadata: "%s"
      contrasts: "%s"
      peptides: "%s"

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
      remove_contaminants: %s
      contaminant_prefix: "cRAP-"
      min_count:
        default: %s
      min_groups: %s
%s
    normalization:
      method: "%s"

    imputation:
      multi_imputation: %s
      method: "%s"
      no_repetitions: 10
      min_no_passed: 8
      width: %s
      downshift: %s
      dep2_method: "%s"
      dep2_random_seed: 1
      qrilc_random_seed: 1

    de:
      method: "%s"
      use_adj_for_pass1: %s
      p_cutoff: %s
      linear_fc_cutoff: %s
      p_adjust_method: "%s"
      fdrtool_correction: %s
      paired: %s
      pairing_col: %s

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

    report:
      project_id: "%s"
      introduction_text: "%s"
      disclaimer_text: "%s"
      show_pca: true
      show_density: true
      show_distance: %s
      show_correlation: true
      show_expression_heatmap: true
      show_imputation_qc: true
      show_identified_proteins_venn: %s
      force_pca: false

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
      gsea_ranking: "%s"
      pathway_volcano: %s

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

    outputs:
      generate_pptx: %s

params:
  seed: 1
',
    project_name, Sys.Date(), project_dir, project_name, round, analyst,
    tolower(project_name), tolower(project_name),
    engine, scale_in,
    protein_dest, sample_map_dest, meta_dest, contrasts_dest, peptide_dest,
    protein_id_col, sample_col,
    map_from, map_to,
    tolower(remove_contaminants), min_valid, min_groups,
    exclude_yaml,
    norm_method,
    tolower(multi_imp), imp_method, imp_width, imp_downshift, dep2_method,
    de_method, tolower(use_adj), p_cutoff, fc_cutoff, p_adjust_method,
    tolower(fdrtool_corr),
    tolower(paired), pairing_col,
    tolower(clustering_on),
    group_col, sample_col,
    report_project_id,
    gsub('"', '\\"', report_intro, fixed = TRUE),
    gsub('"', '\\"', report_disclaimer, fixed = TRUE),
    tolower(show_distance_heatmap), tolower(show_venn),
    selected_organism, custom_mapping_file,
    tolower(pathway_enabled), pathway_method, custom_gmt_file,
    gsea_ranking, tolower(pathway_volcano),
    tolower(ppi_on),
    tolower(adv_stats_on),
    tolower(commentary_enabled), commentary_backend,
    tolower(generate_pptx)
  )

  # Inject technical_report block if extracted
  if (!is.null(tech_report_fields)) {
    escape_yaml <- function(s) gsub('"', '\\"', nna(s), fixed = TRUE)
    tech_yaml <- paste0(
      '    technical_report:\n',
      '      facility: "', escape_yaml(nna(tech_report_fields$facility)), '"\n',
      '      sample_prep: "', escape_yaml(nna(tech_report_fields$sample_prep)), '"\n',
      '      ms_acquisition: "', escape_yaml(nna(tech_report_fields$ms_acquisition)), '"\n',
      '      search_engine: "', escape_yaml(nna(tech_report_fields$search_engine)), '"\n',
      '      search_parameters: "', escape_yaml(nna(tech_report_fields$search_parameters)), '"\n',
      '      acknowledgment: "', escape_yaml(nna(tech_report_fields$acknowledgment)), '"\n'
    )
    config_yaml <- sub("    commentary:", paste0(tech_yaml, "\n    commentary:"), config_yaml)
  }

  writeLines(config_yaml, config_path)
  cat(sprintf("\n  Config saved: %s\n", config_path))

  config_path
}

# --- Metabolomics Wizard ------------------------------------------------------

wizard_metabolomics <- function(project_dir, project_name, analyst, round) {

  # Data input mode: single table vs multi-level directory
  cat("\n--- [METAB] Input Mode ---\n")
  cat("Compound Discoverer exports separate tables per identification level.\n")
  cat("You can either provide a single (possibly pre-merged) table,\n")
  cat("or point to a folder containing one file per level.\n\n")
  input_mode_idx <- ask_choice("How is your data organized?",
    c("Single table (one file, possibly pre-merged levels)",
      "Multi-level folder (one file per level, pipeline merges them)"),
    default = 1)
  is_multi_level <- (input_mode_idx == 2)

  # Data files
  cat("\n--- [METAB] Data Files ---\n")
  cat("Provide ABSOLUTE paths to your input files.\n\n")

  data_path <- NULL
  data_dir_path <- NULL
  level_pattern <- "\\.xlsx$"

  if (is_multi_level) {
    data_dir_path <- ask("Path to folder containing level files (one file per level)")
    if (!dir.exists(data_dir_path)) {
      cat(sprintf("  WARNING: Directory not found: %s\n", data_dir_path))
      data_dir_path <- ask("Please re-enter the folder path")
    }
    if (dir.exists(data_dir_path)) {
      data_dir_path <- normalizePath(data_dir_path)
      level_files <- list.files(data_dir_path, pattern = level_pattern)
      # Exclude Excel temp/lock files (~$...)
      level_files <- level_files[!grepl("^~\\$", level_files)]
      if (length(level_files) == 0) {
        cat("  No .xlsx files found. Trying CSV/TSV...\n")
        level_pattern <- "\\.(csv|tsv|txt)$"
        level_files <- list.files(data_dir_path, pattern = level_pattern)
        level_files <- level_files[!grepl("^~\\$", level_files)]
      }
      cat(sprintf("  Found %d level file(s): %s\n", length(level_files),
                  paste(level_files, collapse = ", ")))
    }
  } else {
    data_path <- ask("Path to metabolomics data file (CD export or processed table, CSV/TSV/XLSX)")
    data_path <- validate_file(data_path, "Metabolomics data file")
    if (is.null(data_path)) {
      data_path <- ask("Please re-enter the metabolomics data file path")
      data_path <- validate_file(data_path, "Metabolomics data file")
    }
  }

  metadata_path <- ask("Path to metadata file (sample info, CSV/TSV — or 'none')", "none")
  if (tolower(metadata_path) == "none") {
    metadata_path <- NULL
  } else {
    metadata_path <- validate_file(metadata_path, "Metadata file")
  }

  # Input format (for single file, or per-level file format in multi-level)
  cat("\n--- [METAB] Input Format ---\n")
  if (is_multi_level) {
    cat("Each level file in the folder will be parsed with the selected format.\n")
  }
  format_idx <- ask_choice("What format is your data file?",
    c("Compound Discoverer raw export (Area: columns)",
      "Already-processed wide table (clean sample columns)"),
    default = 1)
  input_format <- if (is_multi_level) "multi_level" else c("cd_raw", "processed_wide")[format_idx]
  level_file_format <- c("cd_raw", "processed_wide")[format_idx]

  # Column detection & mapping
  cat("\n--- [METAB] Column Detection ---\n")
  feature_id_col <- "null"
  name_col <- "Name"
  mz_col <- "m/z"
  rt_col <- "RT [min]"

  # For multi-level, use the first file in the directory for column detection
  detect_path <- data_path
  if (is_multi_level && !is.null(data_dir_path) && dir.exists(data_dir_path)) {
    lf <- list.files(data_dir_path, pattern = level_pattern, full.names = TRUE)
    if (length(lf) > 0) {
      detect_path <- lf[1]
      cat(sprintf("  (Using first level file for column detection: %s)\n", basename(detect_path)))
    }
  }

  if (!is.null(detect_path)) {
    data_cols <- detect_columns(detect_path)
    cat("Data file columns: ", paste(head(data_cols, 8), collapse = ", "),
        if (length(data_cols) > 8) "..." else "", "\n")

    if (level_file_format == "processed_wide") {
      cat("\nFor processed_wide format, specify the feature ID column.\n")
      feature_id_col <- ask("Feature ID column (or 'null' to construct from Name/mz/RT)", data_cols[1])
      if (tolower(feature_id_col) == "null") feature_id_col <- "null"
    } else {
      cat("\nFor Compound Discoverer format, specify the annotation columns.\n")
      name_col <- ask("Name column", "Name")
      mz_col <- ask("m/z column", "m/z")
      rt_col <- ask("RT column", "RT [min]")
    }
  }

  sample_col <- "sample_id"
  group_col <- "sample_type"
  if (!is.null(metadata_path)) {
    meta_cols <- detect_columns(metadata_path)
    cat("Metadata columns: ", paste(meta_cols, collapse = ", "), "\n")
    sample_col <- ask("Sample ID column in metadata", meta_cols[1])
    cat("\nWhich column defines your biological groups/conditions?\n")
    cat("(Used for PCA coloring and DE contrasts)\n")
    group_col <- ask("Group/condition column", meta_cols[min(2, length(meta_cols))])
  }

  # Build contrasts
  generated_contrasts <- NULL
  if (!is.null(metadata_path)) {
    generated_contrasts <- wizard_build_contrasts(metadata_path, group_col)
  }

  # --- QC Review vs Full Analysis (ask early to skip unnecessary questions) ---
  cat("\n--- [METAB] Run Mode ---\n")
  cat("QC Review mode runs all normalizations (TSS, Median, PQN, EigenMS),\n")
  cat("compares them, and produces a QC summary report. No DE or downstream.\n")
  cat("Use this first to pick the best normalization, then re-run with Full Analysis.\n\n")
  run_mode_idx <- ask_choice("Run mode:",
    c("QC Review — compare normalizations, no DE (recommended first run)",
      "Full Analysis — use selected normalization, run DE + everything"),
    default = 1)
  qc_only <- (run_mode_idx == 1)

  # Normalization — always needed (QC compares methods, Full uses the chosen one)
  cat("\n--- [METAB] Normalization ---\n")
  if (qc_only) {
    cat("QC mode will compare multiple normalization methods automatically.\n")
    cat("Pick a default normalization (used if you later switch to Full mode).\n\n")
  }
  norm_idx <- ask_choice("Sample normalization method:",
    c("PQN — Probabilistic Quotient Normalization (recommended)",
      "Median centering",
      "Sum normalization",
      "EigenMS — SVD-based systematic bias removal (NOREVA)",
      "EigenMS-forced — EigenMS retaining a fixed % of eigentrends",
      "Internal standard — divide each sample by an IS feature",
      "Biological factor — divide each sample by a metadata covariate (e.g. mg protein)",
      "None (data already normalized)"),
    default = 1)
  sample_norm <- c("pqn", "median", "sum", "eigenms", "eigenms_forced",
                   "is", "bio_factor", "none")[norm_idx]

  # Method-specific follow-ups
  is_ref_col <- "null"
  bio_factor_col <- "null"
  if (sample_norm == "is") {
    is_ref_col <- ask("Row metadata column with the IS feature flag (e.g. 'is_internal_standard')",
                      "is_internal_standard")
    is_ref_col <- sprintf('"%s"', is_ref_col)
  } else if (sample_norm == "bio_factor") {
    bio_factor_col <- ask("Metadata column with per-sample numeric factor to divide by (e.g. 'mg_protein')",
                          "mg_protein")
    bio_factor_col <- sprintf('"%s"', bio_factor_col)
  }

  transform_idx <- ask_choice("Transformation:",
    c("log2 (recommended)",
      "log10",
      "Generalized log10 (MetaboAnalyst default — handles zeros smoothly)",
      "None"),
    default = 1)
  transform <- c("log2", "log10", "glog10", "none")[transform_idx]

  scaling_idx <- ask_choice("Scaling:",
    c("None (recommended for most workflows)",
      "Auto-scaling (unit variance)",
      "Pareto scaling (square root of SD)",
      "Range scaling (min-max)"),
    default = 1)
  scaling <- c("none", "auto", "pareto", "range")[scaling_idx]

  chosen_norm <- if (qc_only) "null" else sprintf('"%s"', sample_norm)

  # --- Defaults for QC-only (skip downstream questions) ---
  de_method <- "limma"
  p_cutoff <- 0.05
  fc_cutoff <- 1.5
  run_rf <- FALSE
  run_plsda <- FALSE
  selected_organism <- "Homo sapiens"
  kegg_org_code <- "hsa"
  run_enrichment <- FALSE
  gmt_file <- "null"
  run_mummichog <- FALSE
  commentary_enabled <- FALSE
  commentary_backend <- "none"
  generate_pptx <- FALSE

  if (!qc_only) {
    # DE method & cutoffs
    cat("\n--- [METAB] Differential Expression ---\n")
    de_method_idx <- ask_choice("DE method:",
      c("Limma — moderated t-test (recommended)",
        "Welch t-test (unequal variance)",
        "Student t-test (equal variance, MetaboAnalyst default)",
        "Wilcoxon rank-sum test (non-parametric)"),
      default = 1)
    de_method <- c("limma", "t_test", "t_test_equal", "wilcoxon")[de_method_idx]

    p_cutoff <- as.numeric(ask("Adjusted p-value cutoff", "0.05"))
    fc_cutoff <- as.numeric(ask("Linear fold-change cutoff", "1.5"))

    # Feature selection
    cat("\n--- [METAB] Feature Selection ---\n")
    run_rf <- ask_yn("Run Random Forest feature importance?", TRUE)
    run_plsda <- ask_yn("Run PLS-DA multivariate analysis?", TRUE)

    # Organism
    cat("\n--- [METAB] Organism ---\n")
    cat("Used for mummichog pathway analysis (KEGG organism-specific models).\n\n")
    org_idx <- ask_choice("Which organism does your data come from?",
      c("Human (Homo sapiens)",
        "Mouse (Mus musculus)",
        "Rat (Rattus norvegicus)",
        "Zebrafish (Danio rerio)",
        "C. elegans (Caenorhabditis elegans)",
        "Other (type KEGG organism code, e.g. 'ame' for honey bee)"),
      default = 1)

    # Latin names -> KEGG organism codes for mummichog
    kegg_codes <- c("hsa", "mmu", "rno", "dre", "cel", "other")
    organism_names <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus",
                        "Danio rerio", "Caenorhabditis elegans", "other")
    selected_organism <- organism_names[org_idx]
    kegg_org_code <- kegg_codes[org_idx]

    if (kegg_org_code == "other") {
      selected_organism <- ask("Organism name (e.g. 'Apis mellifera')")
      kegg_org_code <- ask("KEGG organism code (e.g. 'ame' — see https://www.kegg.jp/kegg/catalog/org_list.html)")
    }

    # Enrichment
    cat("\n--- [METAB] Pathway Enrichment ---\n")
    run_enrichment <- ask_yn("Enable pathway enrichment analysis?", FALSE)
    gmt_file <- "null"
    run_mummichog <- FALSE
    if (run_enrichment) {
      gmt_ans <- ask("Path to GMT file for pathway analysis (or 'none')", "none")
      if (tolower(gmt_ans) != "none") {
        gmt_path <- validate_file(gmt_ans, "GMT file")
        if (!is.null(gmt_path)) gmt_file <- paste0('"', gmt_path, '"')
      }
      run_mummichog <- ask_yn(sprintf("Enable mummichog pathway analysis? (uses KEGG %s model)", kegg_org_code), FALSE)
    }

    # AI Commentary
    cat("\n--- [METAB] AI Commentary ---\n")
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

    # --- PowerPoint ---
    generate_pptx <- ask_yn("Generate a PowerPoint summary presentation?", TRUE)
  }

  # --- Technical / Instrument Report (Optional) ---
  nna <- function(x, default = "") if (is.null(x) || !nzchar(x)) default else x
  cat("\n--- [METAB] Technical / Instrument Report (Optional) ---\n")
  cat("Provide a Compound Discoverer / MS-DIAL parameter dump or facility report\n")
  cat("(LOG/TXT/PDF/DOCX) to auto-extract instrument, software, and parameter details.\n")
  tech_report_path <- ask("Path to log file or technical report (or 'none')", "none")
  tech_report_fields <- NULL

  if (tolower(tech_report_path) != "none") {
    tech_report_path <- validate_file(tech_report_path, "Technical report / log file")
    if (!is.null(tech_report_path) && nzchar(Sys.which("claude"))) {
      cat("  Extracting text from file...\n")
      ext <- tolower(tools::file_ext(tech_report_path))
      doc_text <- tryCatch({
        if (ext %in% c("log", "txt", "tsv", "csv")) {
          paste(readLines(tech_report_path, warn = FALSE), collapse = "\n")
        } else if (ext == "pdf") {
          paste(system2("pdftotext", c("-layout", shQuote(tech_report_path), "-"),
                        stdout = TRUE, stderr = FALSE), collapse = "\n")
        } else {
          paste(system2("pandoc", c("-t", "plain", shQuote(tech_report_path)),
                        stdout = TRUE, stderr = FALSE), collapse = "\n")
        }
      }, error = function(e) NULL)

      if (!is.null(doc_text) && nchar(doc_text) >= 50) {
        cat(sprintf("  Extracted %d characters. Sending to Claude for parsing...\n",
                    nchar(doc_text)))
        json_schema <- paste0(
          '{"type":"object","properties":{',
          '"facility":{"type":"string","description":"Name of metabolomics facility / core if mentioned"},',
          '"sample_prep":{"type":"string","description":"Sample preparation, extraction, derivatization, internal standards"},',
          '"ms_acquisition":{"type":"string","description":"Mass spectrometer model and chromatography (e.g. Q Exactive HF, HILIC, RP)"},',
          '"search_engine":{"type":"string","description":"Software name and version (e.g. Compound Discoverer 3.3, MS-DIAL 5.0, MZmine 4)"},',
          '"search_parameters":{"type":"string","description":"Mass tolerance, RT alignment, library, identification confidence"},',
          '"acknowledgment":{"type":"string","description":"Facility acknowledgment if present"}',
          '},"required":["facility","sample_prep","ms_acquisition","search_engine","search_parameters"]}'
        )
        prompt <- sprintf(paste0(
          "Extract structured information from this metabolomics processing log or facility report. ",
          "Return only the fields described in the schema. Use exact details from the text. ",
          "If a field is not mentioned, use an empty string.\n\nLOG/REPORT TEXT:\n%s"), doc_text)
        cmd <- sprintf(
          "claude --print --output-format json --model sonnet --json-schema %s --no-session-persistence %s",
          shQuote(json_schema), shQuote(prompt))
        parsed <- tryCatch({
          raw <- system(cmd, intern = TRUE, timeout = 120)
          p <- jsonlite::fromJSON(paste(raw, collapse = "\n"), simplifyVector = FALSE)
          if (!is.null(p$structured_output)) p$structured_output
          else if (!is.null(p$result) && nzchar(p$result))
            jsonlite::fromJSON(p$result, simplifyVector = FALSE)
          else NULL
        }, error = function(e) { cat(sprintf("  Claude extraction failed: %s\n", e$message)); NULL })

        if (!is.null(parsed) && (!is.null(parsed$search_engine) || !is.null(parsed$ms_acquisition))) {
          cat("\n  Extracted fields:\n")
          for (k in c("facility", "sample_prep", "ms_acquisition", "search_engine",
                      "search_parameters", "acknowledgment")) {
            v <- nna(parsed[[k]])
            if (nzchar(v)) cat(sprintf("    %-16s %s\n", paste0(k, ":"), substr(v, 1, 80)))
          }
          if (ask_yn("\n  Use these extracted fields?", TRUE)) tech_report_fields <- parsed
        }
      }
    } else if (!is.null(tech_report_path)) {
      cat("  Claude Code CLI not found — add 'technical_report:' block manually later.\n")
    }
  }

  # --- Report Section Toggles (Optional) ---
  cat("\n--- [METAB] Report Section Toggles (Optional) ---\n")
  cat("Each section can be toggled on/off in the final HTML report.\n")
  cat("Press Enter to keep all sections ON.\n\n")
  show_qc          <- ask_yn("Show QC section (PCA, density, boxplots, heatmaps)?", TRUE)
  show_de          <- ask_yn("Show DE section (volcano + tables + explorer)?", TRUE)
  show_de_heatmaps <- ask_yn("Show DE heatmaps section?", TRUE)
  show_clustering  <- ask_yn("Show Clustering section (hierarchical / partition / binary)?", TRUE)
  show_fs          <- ask_yn("Show Feature Selection section (RF + PLS-DA)?", TRUE)
  show_enrichment  <- ask_yn("Show Enrichment section (QEA / ssGSEA / ORA / GSEA)?", TRUE)
  show_network     <- ask_yn("Show Metabolite Network section?", TRUE)
  show_mummichog   <- ask_yn("Show Mummichog Pathway section?", TRUE)
  show_methods     <- ask_yn("Show Methods section?", TRUE)
  show_session     <- ask_yn("Show Session Info section?", TRUE)

  # Copy data into project structure
  cat("\n--- [METAB] Setting Up Project ---\n")

  data_dir <- file.path(project_dir, "data", tolower(project_name), "metabolomics")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  data_dest <- ""
  data_dir_dest <- ""

  if (is_multi_level) {
    # Copy entire level folder into project
    data_dir_dest <- file.path("metabolomics", "levels")
    dest_levels <- file.path(project_dir, "data", tolower(project_name), data_dir_dest)
    dir.create(dest_levels, recursive = TRUE, showWarnings = FALSE)
    level_files <- list.files(data_dir_path, pattern = level_pattern, full.names = TRUE)
    level_files <- level_files[!grepl("^~\\$", basename(level_files))]  # skip Excel temp files
    for (lf in level_files) {
      file.copy(lf, file.path(dest_levels, basename(lf)), overwrite = TRUE)
    }
    cat(sprintf("  Copied %d level files -> data/%s/%s/\n",
                length(level_files), tolower(project_name), data_dir_dest))
  } else {
    data_dest <- file.path("metabolomics", basename(data_path))
    file.copy(data_path, file.path(project_dir, "data", tolower(project_name), data_dest),
              overwrite = TRUE)
    cat(sprintf("  Copied data file -> data/%s/%s\n", tolower(project_name), data_dest))
  }

  meta_dest <- ""
  if (!is.null(metadata_path)) {
    meta_dest <- file.path("metabolomics", basename(metadata_path))
    file.copy(metadata_path, file.path(project_dir, "data", tolower(project_name), meta_dest),
              overwrite = TRUE)
    cat(sprintf("  Copied metadata -> data/%s/%s\n", tolower(project_name), meta_dest))
  }

  # Write auto-generated contrasts if built
  contrasts_yaml <- ""
  if (!is.null(generated_contrasts) && nrow(generated_contrasts) > 0) {
    contrast_strings <- paste0(generated_contrasts$Numerator, " - ", generated_contrasts$Denominator)
    contrasts_yaml <- paste0("      contrasts:\n",
      paste(sprintf('        - "%s"', contrast_strings), collapse = "\n"), "\n")
  } else {
    contrasts_yaml <- ""  # let pipeline auto-generate from metadata groups
    cat("  No contrasts defined. Pipeline will auto-generate all pairwise contrasts at runtime.\n")
  }

  # Generate config YAML
  config_name <- paste0("metab_", tolower(gsub("[^a-zA-Z0-9]", "_", project_name)), ".yaml")
  config_path <- file.path(project_dir, "config", config_name)
  dir.create(file.path(project_dir, "config"), showWarnings = FALSE)

  feature_id_yaml <- if (feature_id_col == "null") "null" else paste0('"', feature_id_col, '"')

  # Build files section depending on input mode
  if (is_multi_level) {
    files_yaml <- sprintf('    files:\n      data_dir: "%s"\n      metadata: "%s"',
                          data_dir_dest, meta_dest)
    # Double backslashes for YAML: R string "\\.xlsx$" -> YAML needs "\\\\.xlsx$"
    yaml_safe_pattern <- gsub("\\\\", "\\\\\\\\", level_pattern)
    input_yaml <- sprintf('    input:\n      format: "multi_level"\n      level_pattern: "%s"\n      level_format: "%s"\n      sheet: null',
                          yaml_safe_pattern, level_file_format)
  } else {
    files_yaml <- sprintf('    files:\n      data: "%s"\n      metadata: "%s"',
                          data_dest, meta_dest)
    input_yaml <- sprintf('    input:\n      format: "%s"\n      sheet: null',
                          input_format)
  }

  # Build mummichog section
  mummichog_yaml <- ""
  if (run_mummichog) {
    mummichog_yaml <- sprintf('
      mummichog:
        enabled: true
        organisms:
          - %s
        p_cutoff: 0.05
        tolerance_ppm: 10
        n_permutations: 100
        ionization_mode: null', kegg_org_code)
  }

  # Build preprocessing section (with relaxed thresholds for multi-level)
  if (is_multi_level) {
    preproc_yaml <- sprintf(
      '    preprocessing:\n      chosen_norm: %s\n      feat_missing_threshold: 0.80\n      sample_missing_threshold: 0.99',
      chosen_norm)
  } else {
    preproc_yaml <- sprintf('    preprocessing:\n      chosen_norm: %s', chosen_norm)
  }

  # Build report_sections block (only emit non-default toggles to keep YAML clean)
  rs_lines <- character(0)
  if (!show_qc)          rs_lines <- c(rs_lines, "      qc: false")
  if (!show_de)          rs_lines <- c(rs_lines, "      de: false")
  if (!show_de_heatmaps) rs_lines <- c(rs_lines, "      de_heatmaps: false")
  if (!show_clustering)  rs_lines <- c(rs_lines, "      clustering: false")
  if (!show_fs)          rs_lines <- c(rs_lines, "      feature_selection: false")
  if (!show_enrichment)  rs_lines <- c(rs_lines, "      enrichment: false")
  if (!show_network)     rs_lines <- c(rs_lines, "      network: false")
  if (!show_mummichog)   rs_lines <- c(rs_lines, "      mummichog: false")
  if (!show_methods)     rs_lines <- c(rs_lines, "      methods: false")
  if (!show_session)     rs_lines <- c(rs_lines, "      session_info: false")
  report_sections_yaml <- if (length(rs_lines) > 0) {
    paste0("\n    report_sections:\n", paste(rs_lines, collapse = "\n"), "\n")
  } else ""

  # Build technical_report block from extracted Claude fields
  technical_report_yaml <- ""
  if (!is.null(tech_report_fields)) {
    yaml_esc <- function(s) gsub('"', '\\\\"', nna(s))
    technical_report_yaml <- paste0(
      "\n    technical_report:\n",
      sprintf('      facility: "%s"\n',          yaml_esc(tech_report_fields$facility)),
      sprintf('      sample_prep: "%s"\n',       yaml_esc(tech_report_fields$sample_prep)),
      sprintf('      ms_acquisition: "%s"\n',    yaml_esc(tech_report_fields$ms_acquisition)),
      sprintf('      search_engine: "%s"\n',     yaml_esc(tech_report_fields$search_engine)),
      sprintf('      search_parameters: "%s"\n', yaml_esc(tech_report_fields$search_parameters)),
      sprintf('      acknowledgment: "%s"\n',    yaml_esc(tech_report_fields$acknowledgment))
    )
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
  metabolomics:

%s

%s

%s

    parsing:
      cd_area_prefix: "Area:"
      cd_sample_regex: null

    id_columns:
      feature_id_col: %s
      name_col: "%s"
      mz_col: "%s"
      rt_col: "%s"
      annotation_cols: null

    organism: "%s"

    normalization:
      input_already_normalized: %s
      sample_norm: "%s"
      transform: "%s"
      scaling: "%s"
      pseudocount: 1
      na_policy: "keep"
      is_ref_col: %s
      biological_factor_col: %s

    sample_filter:
      enabled: false
      rules:
        exclude_blanks: true
        exclude_qc: false
        exclude_samples: null

    effects:
      samples: "%s"
      color: "%s"
      shape: null
      heatmap_annotations: null

    qc:
      adaptive_plots: true
      thresholds:
        min_samples_for_pca3d: 10
        max_samples_for_heatmaps: 120

    de:
      method: "%s"
      condition_column: "%s"
%s      p_cutoff: %s
      pval_cutoff: %s
      linear_fc_cutoff: %s

    rf:
      run_rf: %s
      n_trees: 500
      importance: "permutation"
      top_n: 20
      seed: 1234

    plsda:
      run_plsda: %s
      vip_top_n: 15
      colors: null

    enrichment:
      run_enrichment: %s
      libraries: []
      gmt_file: %s
      mapping_file: null%s

    commentary:
      enabled: %s
      backend: "%s"
      claude_code_model: "sonnet"
      max_tokens: 1500
      max_retries: 2

    outputs:
      generate_pptx: %s
%s%s
params:
  seed: 1
',
    project_name, Sys.Date(), project_dir, project_name, round, analyst,
    tolower(project_name), tolower(project_name),
    preproc_yaml,
    input_yaml,
    files_yaml,
    feature_id_yaml, name_col, mz_col, rt_col,
    selected_organism,
    tolower(sample_norm == "none"),
    sample_norm, transform, scaling,
    is_ref_col, bio_factor_col,
    sample_col, group_col,
    de_method, group_col,
    contrasts_yaml,
    p_cutoff, p_cutoff, fc_cutoff,
    tolower(run_rf),
    tolower(run_plsda),
    tolower(run_enrichment), gmt_file, mummichog_yaml,
    tolower(commentary_enabled), commentary_backend,
    tolower(generate_pptx),
    report_sections_yaml,
    technical_report_yaml
  )

  writeLines(config_yaml, config_path)
  cat(sprintf("\n  Config saved: %s\n", config_path))

  if (qc_only) {
    cat("\n  ** QC REVIEW MODE **\n")
    cat("  The pipeline will compare TSS, Median, and PQN normalizations.\n")
    cat("  Review the QC report, then re-run the wizard in Full Analysis mode\n")
    cat("  with your chosen normalization to run DE + downstream analysis.\n\n")
  } else {
    cat(sprintf("\n  ** FULL ANALYSIS MODE (normalization: %s) **\n\n", sample_norm))
  }

  config_path
}

# --- Lipidomics Wizard --------------------------------------------------------

wizard_lipidomics <- function(project_dir, project_name, analyst, round) {

  # Data files
  cat("\n--- [LIPID] Data Files ---\n")
  cat("Provide ABSOLUTE paths to your input files.\n\n")

  data_path <- ask("Path to lipidomics data file (LipidSearch export, translated CSV, or processed table)")
  data_path <- validate_file(data_path, "Lipidomics data file")
  if (is.null(data_path)) {
    data_path <- ask("Please re-enter the lipidomics data file path")
    data_path <- validate_file(data_path, "Lipidomics data file")
  }

  metadata_path <- ask("Path to metadata file (sample info, CSV/TSV — or 'none')", "none")
  if (tolower(metadata_path) == "none") {
    metadata_path <- NULL
  } else {
    metadata_path <- validate_file(metadata_path, "Metadata file")
  }

  # Input format
  cat("\n--- [LIPID] Input Format ---\n")
  format_idx <- ask_choice("What format is your data file?",
    c("LipidSearch export (Excel or CSV with annotation columns)",
      "Translated CSV (Metabolite col + sample cols, optional group row)",
      "Already-processed wide table (clean sample columns)"),
    default = 1)
  input_format <- c("lipidsearch", "translated_csv", "processed_wide")[format_idx]

  # Format-specific options
  sheet_name <- "null"
  has_group_row <- FALSE
  name_col <- "Metabolite"

  if (input_format == "lipidsearch") {
    sheet_name <- ask("Excel sheet name (or 'null' for first sheet)", "RAW Data")
    if (tolower(sheet_name) == "null") sheet_name <- "null"
  } else if (input_format == "translated_csv") {
    has_group_row <- ask_yn("Does the first data row contain group labels (e.g. Oxy/Ana)?", TRUE)
  }

  # Column detection
  cat("\n--- [LIPID] Column Detection ---\n")
  sample_col <- "sample_id"
  group_col <- "condition"

  if (!is.null(data_path)) {
    data_cols <- detect_columns(data_path)
    cat("Data file columns: ", paste(head(data_cols, 8), collapse = ", "),
        if (length(data_cols) > 8) "..." else "", "\n")

    if (input_format == "translated_csv" || input_format == "processed_wide") {
      name_col <- ask("Feature name column", "Metabolite")
    }
  }

  if (!is.null(metadata_path)) {
    meta_cols <- detect_columns(metadata_path)
    cat("Metadata columns: ", paste(meta_cols, collapse = ", "), "\n")
    sample_col <- ask("Sample ID column in metadata", meta_cols[1])
    cat("\nWhich column defines your biological groups/conditions?\n")
    group_col <- ask("Group/condition column", meta_cols[min(2, length(meta_cols))])
  }

  # Build contrasts
  generated_contrasts <- NULL
  if (!is.null(metadata_path)) {
    generated_contrasts <- wizard_build_contrasts(metadata_path, group_col)
  }

  # Normalization
  cat("\n--- [LIPID] Normalization ---\n")
  norm_idx <- ask_choice("Sample normalization method:",
    c("Sum normalization (recommended for lipidomics)",
      "PQN — Probabilistic Quotient Normalization",
      "Median centering",
      "EigenMS — SVD-based systematic bias removal",
      "EigenMS-forced — EigenMS retaining a fixed % of eigentrends",
      "Internal standard — divide each sample by an IS feature",
      "Biological factor — divide each sample by a metadata covariate (e.g. mg protein)",
      "None (data already normalized)"),
    default = 1)
  sample_norm <- c("sum", "pqn", "median", "eigenms", "eigenms_forced",
                   "is", "bio_factor", "none")[norm_idx]

  # Method-specific follow-ups
  is_ref_col <- "null"
  bio_factor_col <- "null"
  if (sample_norm == "is") {
    is_ref_col <- ask("Row metadata column with the IS feature flag (e.g. 'is_internal_standard')",
                      "is_internal_standard")
    is_ref_col <- sprintf('"%s"', is_ref_col)
  } else if (sample_norm == "bio_factor") {
    bio_factor_col <- ask("Metadata column with per-sample numeric factor to divide by (e.g. 'mg_protein')",
                          "mg_protein")
    bio_factor_col <- sprintf('"%s"', bio_factor_col)
  }

  transform_idx <- ask_choice("Transformation:",
    c("log2 (recommended)",
      "log10",
      "Generalized log10 (MetaboAnalyst-style)",
      "None"),
    default = 1)
  transform <- c("log2", "log10", "glog10", "none")[transform_idx]

  scaling_idx <- ask_choice("Scaling:",
    c("Auto-scaling (unit variance, recommended for lipidomics)",
      "Pareto scaling (square root of SD)",
      "None",
      "Range scaling (min-max)"),
    default = 1)
  scaling <- c("auto", "pareto", "none", "range")[scaling_idx]

  na_idx <- ask_choice("Missing value handling:",
    c("Replace with half-minimum per feature (recommended)",
      "Keep as NA",
      "Replace with zero"),
    default = 1)
  na_policy <- c("min_half", "keep", "zero")[na_idx]

  # DE method & cutoffs
  cat("\n--- [LIPID] Differential Expression ---\n")
  de_method_idx <- ask_choice("DE method:",
    c("Limma — moderated t-test (recommended)",
      "Welch t-test (unequal variance)",
      "Student t-test (equal variance)",
      "Wilcoxon rank-sum test (non-parametric)"),
    default = 1)
  de_method <- c("limma", "t_test", "t_test_equal", "wilcoxon")[de_method_idx]

  p_cutoff <- as.numeric(ask("p-value cutoff", "0.05"))
  fc_cutoff <- as.numeric(ask("Linear fold-change cutoff", "1.5"))

  # Feature selection
  cat("\n--- [LIPID] Feature Selection ---\n")
  run_rf <- ask_yn("Run Random Forest feature importance?", TRUE)
  run_plsda <- ask_yn("Run PLS-DA multivariate analysis?", TRUE)

  # AI Commentary
  cat("\n--- [LIPID] AI Commentary ---\n")
  commentary_idx <- ask_choice("Commentary backend:",
    c("None (data-driven fallback only)",
      "Claude Code (uses your subscription, no API key needed)",
      "Claude API (requires ANTHROPIC_API_KEY)",
      "OpenAI API (requires OPENAI_API_KEY)"),
    default = 1)
  commentary_backends <- c("none", "claude-code", "claude", "openai")
  commentary_backend <- commentary_backends[commentary_idx]
  commentary_enabled <- commentary_backend != "none"

  # PowerPoint
  generate_pptx <- ask_yn("Generate a PowerPoint summary presentation?", TRUE)

  # --- Technical / Instrument Report (Optional) ---
  nna <- function(x, default = "") if (is.null(x) || !nzchar(x)) default else x
  cat("\n--- [LIPID] Technical / Instrument Report (Optional) ---\n")
  cat("Provide a LipidSearch / MS-DIAL parameter dump or facility report\n")
  cat("(LOG/TXT/PDF/DOCX) to auto-extract instrument and software details.\n")
  tech_report_path <- ask("Path to log file or technical report (or 'none')", "none")
  tech_report_fields <- NULL

  if (tolower(tech_report_path) != "none") {
    tech_report_path <- validate_file(tech_report_path, "Technical report / log file")
    if (!is.null(tech_report_path) && nzchar(Sys.which("claude"))) {
      cat("  Extracting text from file...\n")
      ext <- tolower(tools::file_ext(tech_report_path))
      doc_text <- tryCatch({
        if (ext %in% c("log", "txt", "tsv", "csv")) {
          paste(readLines(tech_report_path, warn = FALSE), collapse = "\n")
        } else if (ext == "pdf") {
          paste(system2("pdftotext", c("-layout", shQuote(tech_report_path), "-"),
                        stdout = TRUE, stderr = FALSE), collapse = "\n")
        } else {
          paste(system2("pandoc", c("-t", "plain", shQuote(tech_report_path)),
                        stdout = TRUE, stderr = FALSE), collapse = "\n")
        }
      }, error = function(e) NULL)

      if (!is.null(doc_text) && nchar(doc_text) >= 50) {
        cat(sprintf("  Extracted %d characters. Sending to Claude for parsing...\n",
                    nchar(doc_text)))
        json_schema <- paste0(
          '{"type":"object","properties":{',
          '"facility":{"type":"string","description":"Lipidomics facility / core if mentioned"},',
          '"sample_prep":{"type":"string","description":"Lipid extraction protocol, internal standards (e.g. SPLASH Lipidomix)"},',
          '"ms_acquisition":{"type":"string","description":"Mass spectrometer model and chromatography (e.g. Q Exactive, HILIC, RP, normal-phase)"},',
          '"search_engine":{"type":"string","description":"Software name and version (e.g. LipidSearch 5.0, MS-DIAL 5.0, LipidMatch)"},',
          '"search_parameters":{"type":"string","description":"Mass tolerance, lipid class library, identification grade, FA scoring"},',
          '"acknowledgment":{"type":"string","description":"Facility acknowledgment if present"}',
          '},"required":["facility","sample_prep","ms_acquisition","search_engine","search_parameters"]}'
        )
        prompt <- sprintf(paste0(
          "Extract structured information from this lipidomics processing log or facility report. ",
          "Return only the fields described in the schema. Use exact details from the text. ",
          "If a field is not mentioned, use an empty string.\n\nLOG/REPORT TEXT:\n%s"), doc_text)
        cmd <- sprintf(
          "claude --print --output-format json --model sonnet --json-schema %s --no-session-persistence %s",
          shQuote(json_schema), shQuote(prompt))
        parsed <- tryCatch({
          raw <- system(cmd, intern = TRUE, timeout = 120)
          p <- jsonlite::fromJSON(paste(raw, collapse = "\n"), simplifyVector = FALSE)
          if (!is.null(p$structured_output)) p$structured_output
          else if (!is.null(p$result) && nzchar(p$result))
            jsonlite::fromJSON(p$result, simplifyVector = FALSE)
          else NULL
        }, error = function(e) { cat(sprintf("  Claude extraction failed: %s\n", e$message)); NULL })

        if (!is.null(parsed) && (!is.null(parsed$search_engine) || !is.null(parsed$ms_acquisition))) {
          cat("\n  Extracted fields:\n")
          for (k in c("facility", "sample_prep", "ms_acquisition", "search_engine",
                      "search_parameters", "acknowledgment")) {
            v <- nna(parsed[[k]])
            if (nzchar(v)) cat(sprintf("    %-16s %s\n", paste0(k, ":"), substr(v, 1, 80)))
          }
          if (ask_yn("\n  Use these extracted fields?", TRUE)) tech_report_fields <- parsed
        }
      }
    } else if (!is.null(tech_report_path)) {
      cat("  Claude Code CLI not found — add 'technical_report:' block manually later.\n")
    }
  }

  # --- Report Section Toggles (Optional) ---
  cat("\n--- [LIPID] Report Section Toggles (Optional) ---\n")
  cat("Each section can be toggled on/off in the final HTML report.\n")
  cat("Press Enter to keep all sections ON.\n\n")
  show_qc          <- ask_yn("Show QC section (PCA, density, boxplots, heatmaps)?", TRUE)
  show_de          <- ask_yn("Show DE section (volcano + tables)?", TRUE)
  show_de_heatmaps <- ask_yn("Show DE heatmaps section?", TRUE)
  show_clustering  <- ask_yn("Show Clustering section (hierarchical / partition / binary)?", TRUE)
  show_class       <- ask_yn("Show Lipid Class Analysis section?", TRUE)
  show_fs          <- ask_yn("Show Feature Selection section (RF + PLS-DA)?", TRUE)
  show_biomarker   <- ask_yn("Show Biomarker Discovery section?", TRUE)
  show_pathway     <- ask_yn("Show Lipid Pathway section?", TRUE)
  show_methods     <- ask_yn("Show Methods section?", TRUE)
  show_session     <- ask_yn("Show Session Info section?", TRUE)
  show_technical   <- ask_yn("Show Technical Summary section?", TRUE)

  # Copy data into project structure
  cat("\n--- [LIPID] Setting Up Project ---\n")

  data_dir <- file.path(project_dir, "data", tolower(project_name), "lipidomics")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  data_dest <- file.path("lipidomics", basename(data_path))
  file.copy(data_path, file.path(project_dir, "data", tolower(project_name), data_dest),
            overwrite = TRUE)
  cat(sprintf("  Copied data file -> data/%s/%s\n", tolower(project_name), data_dest))

  meta_dest <- ""
  if (!is.null(metadata_path)) {
    meta_dest <- file.path("lipidomics", basename(metadata_path))
    file.copy(metadata_path, file.path(project_dir, "data", tolower(project_name), meta_dest),
              overwrite = TRUE)
    cat(sprintf("  Copied metadata -> data/%s/%s\n", tolower(project_name), meta_dest))
  }

  # Write auto-generated contrasts if built
  contrasts_yaml <- ""
  if (!is.null(generated_contrasts) && nrow(generated_contrasts) > 0) {
    contrast_strings <- paste0(generated_contrasts$Numerator, " - ", generated_contrasts$Denominator)
    contrasts_yaml <- paste0("      contrasts:\n",
      paste(sprintf('        - "%s"', contrast_strings), collapse = "\n"), "\n")
  } else {
    contrasts_yaml <- ""  # let pipeline auto-generate from metadata groups
    cat("  No contrasts defined. Pipeline will auto-generate all pairwise contrasts at runtime.\n")
  }

  # Generate config YAML
  config_name <- paste0("lipid_", tolower(gsub("[^a-zA-Z0-9]", "_", project_name)), ".yaml")
  config_path <- file.path(project_dir, "config", config_name)
  dir.create(file.path(project_dir, "config"), showWarnings = FALSE)

  sheet_yaml <- if (sheet_name == "null") "null" else paste0('"', sheet_name, '"')

  # Build report_sections block (only emit non-default toggles to keep YAML clean)
  rs_lines <- character(0)
  if (!show_qc)          rs_lines <- c(rs_lines, "      qc: false")
  if (!show_de)          rs_lines <- c(rs_lines, "      de: false")
  if (!show_de_heatmaps) rs_lines <- c(rs_lines, "      de_heatmaps: false")
  if (!show_clustering)  rs_lines <- c(rs_lines, "      clustering: false")
  if (!show_class)       rs_lines <- c(rs_lines, "      class_analysis: false")
  if (!show_fs)          rs_lines <- c(rs_lines, "      feature_selection: false")
  if (!show_biomarker)   rs_lines <- c(rs_lines, "      biomarker: false")
  if (!show_pathway)     rs_lines <- c(rs_lines, "      pathway: false")
  if (!show_methods)     rs_lines <- c(rs_lines, "      methods: false")
  if (!show_session)     rs_lines <- c(rs_lines, "      session_info: false")
  if (!show_technical)   rs_lines <- c(rs_lines, "      technical_summary: false")
  report_sections_yaml <- if (length(rs_lines) > 0) {
    paste0("\n    report_sections:\n", paste(rs_lines, collapse = "\n"), "\n")
  } else ""

  technical_report_yaml <- ""
  if (!is.null(tech_report_fields)) {
    yaml_esc <- function(s) gsub('"', '\\\\"', nna(s))
    technical_report_yaml <- paste0(
      "\n    technical_report:\n",
      sprintf('      facility: "%s"\n',          yaml_esc(tech_report_fields$facility)),
      sprintf('      sample_prep: "%s"\n',       yaml_esc(tech_report_fields$sample_prep)),
      sprintf('      ms_acquisition: "%s"\n',    yaml_esc(tech_report_fields$ms_acquisition)),
      sprintf('      search_engine: "%s"\n',     yaml_esc(tech_report_fields$search_engine)),
      sprintf('      search_parameters: "%s"\n', yaml_esc(tech_report_fields$search_parameters)),
      sprintf('      acknowledgment: "%s"\n',    yaml_esc(tech_report_fields$acknowledgment))
    )
  }

  config_yaml <- sprintf('# Auto-generated by run.R wizard (lipidomics)
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
  lipidomics:

    input:
      format: "%s"
      sheet: %s
      has_group_row: %s

    files:
      data: "%s"
      metadata: "%s"

    id_columns:
      name_col: "%s"

    normalization:
      sample_norm: "%s"
      transform: "%s"
      scaling: "%s"
      pseudocount: 1
      na_policy: "%s"
      is_ref_col: %s
      biological_factor_col: %s

    sample_filter:
      enabled: false

    effects:
      samples: "%s"
      color: "%s"

    qc:
      adaptive_plots: true
      thresholds:
        min_samples_for_pca3d: 10
        max_samples_for_heatmaps: 120

    de:
      method: "%s"
      condition_column: "%s"
%s      p_cutoff: %s
      linear_fc_cutoff: %s

    rf:
      run_rf: %s
      n_trees: 500
      importance: "permutation"
      top_n: 20
      seed: 1234

    plsda:
      run_plsda: %s
      vip_top_n: 15

    commentary:
      enabled: %s
      backend: "%s"
      claude_code_model: "sonnet"
      max_tokens: 1500
      max_retries: 2

    outputs:
      generate_pptx: %s
%s%s
params:
  seed: 1
',
    project_name, Sys.Date(), project_dir, project_name, round, analyst,
    tolower(project_name), tolower(project_name),
    input_format, sheet_yaml, tolower(has_group_row),
    data_dest, meta_dest,
    name_col,
    sample_norm, transform, scaling, na_policy,
    is_ref_col, bio_factor_col,
    sample_col, group_col,
    de_method, group_col,
    contrasts_yaml,
    p_cutoff, fc_cutoff,
    tolower(run_rf),
    tolower(run_plsda),
    tolower(commentary_enabled), commentary_backend,
    tolower(generate_pptx),
    report_sections_yaml,
    technical_report_yaml
  )

  writeLines(config_yaml, config_path)
  cat(sprintf("\n  Config saved: %s\n", config_path))

  config_path
}

# --- Multiomics Integration Wizard --------------------------------------------

wizard_multiomics <- function(project_dir, project_name, analyst, round) {

  cat("\n--- Multi-Omics Integration Setup ---\n")
  cat("This pipeline integrates results from individual omics layers\n")
  cat("(RNA-seq, Proteomics, Metabolomics, Lipidomics) using DIABLO, SNF, MOFA2.\n\n")

  # Input mode: from configs or from existing results
  input_idx <- ask_choice("How will you provide the omics data?",
    c("From existing configs (run individual pipelines first, then integrate)",
      "From existing results (point to results folders or payload RDS files)"),
    default = 1)
  input_mode <- c("pipeline", "outputs")[input_idx]

  # Collect individual layer configs or payload paths
  layer_configs <- list()

  if (input_mode == "pipeline") {
    cat("\nPoint to existing config YAML files for each omics layer.\n")
    cat("Leave blank to skip a layer (at least 2 required).\n\n")

    rna_cfg <- ask("RNA-seq config YAML path (or blank to skip)", "")
    if (nzchar(rna_cfg)) {
      rna_cfg <- validate_file(rna_cfg, "RNA-seq config")
      if (!is.null(rna_cfg)) layer_configs$rna <- rna_cfg
    }

    prot_cfg <- ask("Proteomics config YAML path (or blank to skip)", "")
    if (nzchar(prot_cfg)) {
      prot_cfg <- validate_file(prot_cfg, "Proteomics config")
      if (!is.null(prot_cfg)) layer_configs$proteomics <- prot_cfg
    }

    metab_cfg <- ask("Metabolomics config YAML path (or blank to skip)", "")
    if (nzchar(metab_cfg)) {
      metab_cfg <- validate_file(metab_cfg, "Metabolomics config")
      if (!is.null(metab_cfg)) layer_configs$metabolomics <- metab_cfg
    }

    lipid_cfg <- ask("Lipidomics config YAML path (or blank to skip)", "")
    if (nzchar(lipid_cfg)) {
      lipid_cfg <- validate_file(lipid_cfg, "Lipidomics config")
      if (!is.null(lipid_cfg)) layer_configs$lipidomics <- lipid_cfg
    }

    if (length(layer_configs) < 2) {
      stop("Multi-omics integration requires at least 2 omics layers. Got: ",
           length(layer_configs), call. = FALSE)
    }

    cat(sprintf("\n  Layers selected: %s\n", paste(names(layer_configs), collapse = ", ")))

  } else {
    # Outputs mode: point to results folders or individual payload RDS files
    cat("\nProvide results folders from previous pipeline runs.\n")
    cat("The wizard will scan for shiny_payload_*.rds files automatically.\n")
    cat("You can provide one folder per layer, or multiple folders.\n")
    cat("Enter paths one at a time. Type 'done' when finished.\n\n")

    result_folders <- character(0)
    repeat {
      folder <- ask("Results folder path (or 'done' to finish)", "done")
      if (tolower(folder) == "done") break
      folder <- normalizePath(folder, mustWork = TRUE)
      result_folders <- c(result_folders, folder)
    }

    # Scan all provided folders for shiny_payload_*.rds files
    payload_map <- c(
      shiny_payload_rnaseq       = "rna",
      shiny_payload_proteomics   = "proteomics",
      shiny_payload_metabolomics = "metabolomics",
      shiny_payload_lipidomics   = "lipidomics"
    )

    for (folder in result_folders) {
      rds_files <- list.files(folder, pattern = "shiny_payload_.*\\.rds$",
                              recursive = TRUE, full.names = TRUE)
      for (rds in rds_files) {
        base <- tools::file_path_sans_ext(basename(rds))
        if (base %in% names(payload_map)) {
          layer_name <- payload_map[[base]]
          if (!layer_name %in% names(layer_configs)) {
            layer_configs[[layer_name]] <- rds
            cat(sprintf("  Found %s: %s\n", layer_name, rds))
          }
        }
      }
    }

    # If no payloads found from folder scan, fall back to manual entry
    if (length(layer_configs) == 0) {
      cat("\nNo payload files found in provided folders.\n")
      cat("Enter payload RDS paths manually (blank to skip):\n\n")

      rna_rds <- ask("RNA-seq payload RDS path", "")
      if (nzchar(rna_rds)) layer_configs$rna <- normalizePath(rna_rds, mustWork = TRUE)

      prot_rds <- ask("Proteomics payload RDS path", "")
      if (nzchar(prot_rds)) layer_configs$proteomics <- normalizePath(prot_rds, mustWork = TRUE)

      metab_rds <- ask("Metabolomics payload RDS path", "")
      if (nzchar(metab_rds)) layer_configs$metabolomics <- normalizePath(metab_rds, mustWork = TRUE)

      lipid_rds <- ask("Lipidomics payload RDS path", "")
      if (nzchar(lipid_rds)) layer_configs$lipidomics <- normalizePath(lipid_rds, mustWork = TRUE)
    }

    if (length(layer_configs) < 2) {
      stop("Multi-omics integration requires at least 2 omics layers. Found: ",
           length(layer_configs), call. = FALSE)
    }

    cat(sprintf("\n  Layers discovered: %s\n", paste(names(layer_configs), collapse = ", ")))
  }

  # Condition column
  condition_col <- ask("Condition column name (in metadata)", "condition")

  # Sample mapping (optional, for matching samples across omics)
  sample_mapping <- ask("Sample mapping CSV path (optional, blank to skip)", "")

  # Gene-protein mapping
  gene_prot_map <- ask("Gene-protein mapping file (optional, blank to skip)", "")

  # Organism
  organism <- ask("Organism", "Homo sapiens")

  # Output subdirectory
  out_subdir <- ask("Output subdirectory (under outputs/)",
                    paste0("outputs/", tolower(gsub(" ", "_", project_name))))

  # Commentary
  cat("\n--- Commentary ---\n")
  commentary_enabled <- tolower(substr(ask("Enable AI commentary? (y/n)", "y"), 1, 1)) == "y"
  commentary_backend <- "none"
  if (commentary_enabled) {
    backend_idx <- ask_choice("Commentary backend:",
      c("claude-code", "none (generate placeholders)"),
      default = 1)
    commentary_backend <- c("claude-code", "none")[backend_idx]
  }

  # Integration methods
  cat("\n--- Integration Methods ---\n")
  use_diablo <- tolower(substr(ask("Enable DIABLO (supervised)? (y/n)", "y"), 1, 1)) == "y"
  use_snf    <- tolower(substr(ask("Enable SNF (unsupervised clustering)? (y/n)", "y"), 1, 1)) == "y"
  use_mofa2  <- tolower(substr(ask("Enable MOFA2 (factor analysis)? (y/n)", "y"), 1, 1)) == "y"

  methods <- character(0)
  if (use_diablo) methods <- c(methods, "DIABLO")
  if (use_snf)    methods <- c(methods, "SNF")
  if (use_mofa2)  methods <- c(methods, "MOFA2")

  if (length(methods) == 0) {
    cat("Warning: no methods selected, defaulting to DIABLO + SNF\n")
    methods <- c("DIABLO", "SNF")
  }

  # DIABLO cv_folds
  cv_folds <- 3
  if (use_diablo) {
    cv_folds <- as.integer(ask("DIABLO cross-validation folds", "3"))
  }

  # Build omics_present list from discovered layers
  omics_key_map <- c(rna = "rna", proteomics = "proteomics",
                     metabolomics = "metabolomics", lipidomics = "lipidomics")
  omics_present <- unname(omics_key_map[names(layer_configs)])

  # Build config YAML
  config_lines <- c(
    "# Auto-generated by multiomics-core wizard (multi-omics integration)",
    sprintf("# Generated: %s", Sys.time()),
    "",
    "project:",
    sprintf("  dir: \"%s\"", project_dir),
    sprintf("  name: \"%s\"", project_name),
    sprintf("  analysis_round: \"%s\"", round),
    sprintf("  analyst: \"%s\"", analyst),
    "",
    "paths:",
    "  raw: \"data\"",
    sprintf("  out: \"%s\"", out_subdir),
    "",
    "params:",
    "  seed: 1234",
    "",
    "design:",
    sprintf("  condition_column: \"%s\"", condition_col),
    "",
    "modes:"
  )

  if (input_mode == "pipeline") {
    # Merge individual configs' modes sections
    for (layer_name in names(layer_configs)) {
      layer_cfg <- yaml::read_yaml(layer_configs[[layer_name]])
      if (!is.null(layer_cfg$modes[[layer_name]])) {
        layer_yaml <- yaml::as.yaml(list(modes = layer_cfg$modes[layer_name]))
        layer_lines <- strsplit(layer_yaml, "\n")[[1]]
        layer_lines <- layer_lines[layer_lines != "modes:"]
        config_lines <- c(config_lines, layer_lines)
      } else {
        first_key <- names(layer_cfg$modes)[1]
        if (!is.null(first_key)) {
          renamed <- layer_cfg$modes[1]
          names(renamed) <- layer_name
          layer_yaml <- yaml::as.yaml(renamed)
          config_lines <- c(config_lines,
            paste0("  ", strsplit(layer_yaml, "\n")[[1]]))
        }
      }
    }
  }

  # Multiomics integration section
  config_lines <- c(config_lines,
    "  multiomics:",
    sprintf("    condition_column: \"%s\"", condition_col),
    sprintf("    input_mode: \"%s\"", input_mode)
  )

  if (nzchar(gene_prot_map)) {
    config_lines <- c(config_lines,
      sprintf("    gene_protein_mapping_file: \"%s\"", gene_prot_map))
  } else {
    config_lines <- c(config_lines,
      "    gene_protein_mapping_file: ~")
  }

  config_lines <- c(config_lines,
    "    require_one_to_one_mapping: no")

  config_lines <- c(config_lines,
    "",
    "    feature_selection:",
    "      method: \"variance\"",
    "      top_n: 500",
    "      min_var_percentile: 50"
  )

  # Integration block
  config_lines <- c(config_lines,
    "",
    "    integration:",
    "      methods:"
  )
  for (m in methods) {
    config_lines <- c(config_lines, sprintf("        - \"%s\"", m))
  }

  if (use_diablo) {
    config_lines <- c(config_lines,
      "      diablo:",
      "        ncomp: 2",
      "        design_matrix: \"full\"",
      sprintf("        cv_folds: %d", cv_folds))
  }

  if (use_snf) {
    config_lines <- c(config_lines,
      "      snf:",
      "        K: 20",
      "        alpha: 0.5",
      "        T: 20",
      "        n_clusters: 2")
  }

  # Commentary
  config_lines <- c(config_lines,
    "",
    "    commentary:",
    sprintf("      enabled: %s", if (commentary_enabled) "yes" else "no"),
    sprintf("      backend: \"%s\"", commentary_backend)
  )

  # Enrichment
  config_lines <- c(config_lines,
    "",
    "    enrichment:",
    "      run_enrichment: yes"
  )

  # Payload paths (outputs mode)
  if (input_mode == "outputs") {
    config_lines <- c(config_lines,
      "",
      "    payload_paths:")
    for (nm in names(layer_configs)) {
      config_lines <- c(config_lines,
        sprintf("      %s: \"%s\"", nm, layer_configs[[nm]]))
    }
  }

  # Sample mapping
  if (nzchar(sample_mapping)) {
    config_lines <- c(config_lines,
      sprintf("    sample_mapping: \"%s\"", sample_mapping))
  }

  # Global section
  config_lines <- c(config_lines,
    "",
    "global:",
    "  omics_present:"
  )
  for (op in omics_present) {
    config_lines <- c(config_lines, sprintf("    - %s", op))
  }
  config_lines <- c(config_lines,
    sprintf("  organism: \"%s\"", organism)
  )

  config_yaml <- paste(config_lines, collapse = "\n")

  # Save
  config_name <- paste0("multiomics_", tolower(gsub("[^a-zA-Z0-9]", "_", project_name)), ".yaml")
  config_path <- file.path(project_dir, "config", config_name)
  dir.create(dirname(config_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(config_yaml, config_path)
  cat(sprintf("\n  Multi-omics config saved: %s\n", config_path))

  config_path
}

# --- Pipeline Runner ----------------------------------------------------------

run_pipeline <- function(config_path, fresh = FALSE) {
  config_path <- normalizePath(config_path, mustWork = TRUE)
  cat(sprintf("\n  Config: %s\n", config_path))

  # Copy config to project root config.yaml so _targets.R always finds it.
  # This avoids callr subprocess env-var propagation issues.
  cfg_tmp <- yaml::read_yaml(config_path)
  proj_dir <- cfg_tmp$project$dir %||% getwd()
  root_config <- file.path(proj_dir, "config.yaml")
  if (normalizePath(config_path) != normalizePath(root_config, mustWork = FALSE)) {
    file.copy(config_path, root_config, overwrite = TRUE)
    cat(sprintf("  Linked config -> %s\n", root_config))
  }
  Sys.setenv(MULTIOMICS_CONFIG = config_path)

  if (fresh) {
    cat("  Clearing targets cache (fresh run)...\n")
    targets::tar_destroy(ask = FALSE)
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

  # Pre-flight: metabolomics input files
  if (!is.null(cfg$modes$metabolomics)) {
    metab <- cfg$modes$metabolomics
    base_dir <- file.path(cfg$project$dir, cfg$paths$raw)
    missing_metab <- character(0)
    for (fname in c("data", "metadata")) {
      fpath <- metab$files[[fname]]
      if (!is.null(fpath) && nzchar(fpath) && !file.exists(file.path(base_dir, fpath))) {
        missing_metab <- c(missing_metab, sprintf("  [x] metabolomics/%s: %s", fname, fpath))
      }
    }
    if (length(missing_metab) > 0) {
      cat("\n  PRE-FLIGHT: Missing metabolomics input files:\n")
      cat(paste(missing_metab, collapse = "\n"), "\n")
      if (interactive()) {
        proceed <- ask_yn("Continue anyway?", FALSE)
        if (!proceed) return(invisible(NULL))
      } else {
        stop("Pre-flight check failed: missing metabolomics input files")
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

  metab_report_path <- file.path(run_dir_path, "metabolomics_report.html")
  if (file.exists(metab_report_path)) {
    cat(sprintf("\n  Metabolomics Report: %s\n", metab_report_path))
    cat("  Open in browser to view results.\n")
  }

  metab_summary_path <- file.path(run_dir_path, "metabolomics_pipeline_summary.html")
  if (file.exists(metab_summary_path)) {
    cat(sprintf("  Metabolomics Summary: %s\n", metab_summary_path))
  }

  # PowerPoint presentations
  pptx_files <- list.files(run_dir_path, pattern = "_summary\\.pptx$", full.names = TRUE)
  if (length(pptx_files) > 0) {
    cat("\n  PowerPoint Presentations:\n")
    for (pf in pptx_files) cat(sprintf("    %s\n", pf))
  }

  cat("\n  Pipeline complete.\n")
}

# --- Main Entry Point ---------------------------------------------------------

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  project_dir <- getwd()

  if ("--help" %in% args || "-h" %in% args) {
    cli_header()
    cat("Usage:\n")
    cat("  Rscript run.R --wizard           # Open HTML wizard in browser (GUI)\n")
    cat("  Rscript run.R --new              # Interactive wizard for new project (CLI)\n")
    cat("  Rscript run.R --config path.yaml # Run with existing config\n")
    cat("  Rscript run.R                    # Re-run last config (with caching)\n")
    cat("  Rscript run.R --fresh            # Full re-run (invalidate cache)\n")
    cat("  Rscript run.R --help             # Show this help message\n\n")
    return(invisible())
  }

  if ("--wizard" %in% args) {
    check_dependencies("wizard")

    # Serve HTML wizard via httpuv with API backend
    wizard_path <- file.path(project_dir, "wizard.html")
    if (!file.exists(wizard_path)) stop("wizard.html not found in ", project_dir)

    cli_header()
    cat("--- Config Wizard (browser) ---\n\n")

    # State: pipeline child process handle + log path for live-tailing
    .wizard_state <- new.env(parent = emptyenv())
    .wizard_state$run_config <- NULL
    .wizard_state$child      <- NULL
    .wizard_state$log_path   <- NULL
    .wizard_state$config_path <- NULL

    # Find available port
    port <- 8080L
    for (p in 8080:8099) {
      tryCatch({
        srv_test <- httpuv::startServer("127.0.0.1", p, list(call = function(req) {
          list(status = 200L, headers = list(), body = "")
        }))
        httpuv::stopServer(srv_test)
        port <- as.integer(p)
        break
      }, error = function(e) NULL)
    }

    # Parse multipart form data (minimal implementation)
    parse_multipart <- function(body_raw, content_type) {
      boundary <- sub(".*boundary=", "", content_type)
      body_str <- rawToChar(body_raw)
      parts <- strsplit(body_str, paste0("--", boundary))[[1]]
      parts <- parts[nchar(trimws(parts)) > 4]  # skip empty/closing
      result <- list()
      for (part in parts) {
        if (!grepl("Content-Disposition", part)) next
        name_match <- regmatches(part, regexpr('name="([^"]+)"', part))
        if (length(name_match) == 0) next
        name <- sub('name="', '', sub('"$', '', name_match))
        fname_match <- regmatches(part, regexpr('filename="([^"]*)"', part))
        # Extract body after double newline
        body_part <- sub("^.*?\r?\n\r?\n", "", part)
        body_part <- sub("\r?\n--$", "", sub("\r?\n$", "", body_part))
        if (length(fname_match) > 0) {
          fname <- sub('filename="', '', sub('"$', '', fname_match))
          result[[name]] <- list(filename = fname, data = body_part)
        } else {
          result[[name]] <- body_part
        }
      }
      result
    }

    # HTTP handler
    app <- list(
      call = function(req) {
        path <- req$PATH_INFO
        method <- req$REQUEST_METHOD

        # CORS headers for all responses
        cors <- list(
          "Access-Control-Allow-Origin" = "*",
          "Access-Control-Allow-Methods" = "GET, POST, OPTIONS",
          "Access-Control-Allow-Headers" = "Content-Type"
        )

        if (method == "OPTIONS") {
          return(list(status = 200L, headers = cors, body = ""))
        }

        # Serve wizard.html
        if (path == "/" && method == "GET") {
          html <- paste(readLines(wizard_path, warn = FALSE), collapse = "\n")
          return(list(
            status = 200L,
            headers = c(cors, list("Content-Type" = "text/html; charset=utf-8")),
            body = html
          ))
        }

        # API: serve wizard defaults
        if (path == "/api/defaults" && method == "GET") {
          defaults_path <- file.path(project_dir, "wizard_defaults.json")
          if (file.exists(defaults_path)) {
            defaults_json <- paste(readLines(defaults_path, warn = FALSE), collapse = "\n")
          } else {
            defaults_json <- "{}"
          }
          return(list(
            status = 200L,
            headers = c(cors, list("Content-Type" = "application/json")),
            body = defaults_json
          ))
        }

        # API: upload file (save to data dir, return path + detected columns)
        if (path == "/api/upload-file" && method == "POST") {
          ct <- req$CONTENT_TYPE %||% ""
          body_raw <- req$rook.input$read()
          parts <- parse_multipart(body_raw, ct)

          proj_name <- trimws(parts[["project_name"]] %||% "project")
          mode <- trimws(parts[["mode"]] %||% "proteomics")
          file_role <- trimws(parts[["file_role"]] %||% "data")

          file_part <- parts[["file"]]
          if (is.null(file_part)) {
            return(list(status = 400L, headers = cors,
                        body = '{"error":"No file uploaded"}'))
          }

          # Save file — sanitize project name to match YAML path generation
          proj_name_safe <- gsub("[^a-zA-Z0-9_-]", "", tolower(proj_name))
          sub_dir <- file.path(project_dir, "data", proj_name_safe, mode)
          dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
          dest <- file.path(sub_dir, file_part$filename)
          # Write in binary mode to preserve exact bytes — writeLines in text
          # mode on Windows doubles \r\n into \r\r\n, corrupting large CSVs.
          con <- file(dest, "wb")
          writeBin(charToRaw(file_part$data), con)
          close(con)

          # Detect columns from first line
          first_line <- strsplit(file_part$data, "\r?\n")[[1]][1]
          sep <- if (grepl("\t", first_line)) "\t" else ","
          cols <- trimws(gsub('["\']', '', strsplit(first_line, sep)[[1]]))

          # Return the relative path (for YAML) and columns
          rel_path <- file.path(mode, file_part$filename)
          resp <- jsonlite::toJSON(list(
            path = rel_path,
            full_path = dest,
            columns = cols
          ), auto_unbox = TRUE)

          return(list(status = 200L,
                      headers = c(cors, list("Content-Type" = "application/json")),
                      body = resp))
        }

        # API: load example proteomics dataset (returns paths + columns for auto-fill)
        if (path == "/api/load-example" && method == "GET") {
          ex_dir <- file.path(project_dir, "data", "example_proteomics", "proteomics")
          if (!dir.exists(ex_dir)) {
            return(list(status = 404L, headers = c(cors,
                          list("Content-Type" = "application/json")),
                        body = '{"error":"Example dataset missing. Run data/example_proteomics/proteomics/generate_example.R"}'))
          }

          read_header <- function(fn) {
            fp <- file.path(ex_dir, fn)
            if (!file.exists(fp)) return(list(path = fn, columns = character(0), text = ""))
            txt  <- paste(readLines(fp, warn = FALSE), collapse = "\n")
            first <- strsplit(txt, "\r?\n")[[1]][1]
            sep <- if (grepl("\t", first)) "\t" else ","
            cols <- trimws(gsub('["\']', '', strsplit(first, sep)[[1]]))
            list(path = fn, columns = cols, text = txt)
          }

          resp_obj <- list(
            project_name  = "example_proteomics",
            analyst       = "Example Run",
            round         = "A01",
            organism      = "human",
            mode          = "proteomics",
            protein       = read_header("protein_matrix.csv"),
            metadata      = read_header("metadata.csv"),
            sample_map    = read_header("sample_map.csv"),
            contrasts     = read_header("contrasts.csv")
          )
          # Strip metadata text only for non-metadata files to keep response small
          resp_obj$protein$text    <- NULL
          resp_obj$sample_map$text <- NULL
          resp_obj$contrasts$text  <- NULL

          return(list(status = 200L,
                      headers = c(cors, list("Content-Type" = "application/json")),
                      body = jsonlite::toJSON(resp_obj, auto_unbox = TRUE)))
        }

        # API: open a local file with the OS default handler (for report HTML)
        if (path == "/api/open-file" && method == "POST") {
          body_raw <- req$rook.input$read()
          body <- jsonlite::fromJSON(rawToChar(body_raw))
          target <- body$path
          if (is.null(target) || !nzchar(target) || !file.exists(target)) {
            return(list(status = 404L, headers = cors,
                        body = '{"error":"File not found"}'))
          }
          ok <- tryCatch({
            if (.Platform$OS.type == "windows") {
              shell.exec(target)
            } else {
              cmd <- if (Sys.info()[["sysname"]] == "Darwin") "open" else "xdg-open"
              system2(cmd, shQuote(target), wait = FALSE)
            }
            TRUE
          }, error = function(e) FALSE)
          resp <- jsonlite::toJSON(list(opened = ok), auto_unbox = TRUE)
          return(list(status = if (ok) 200L else 500L,
                      headers = c(cors, list("Content-Type" = "application/json")),
                      body = resp))
        }

        # API: list saved configs
        if (path == "/api/list-configs" && method == "GET") {
          config_dir <- file.path(project_dir, "config")
          files <- if (dir.exists(config_dir)) {
            list.files(config_dir, pattern = "\\.ya?ml$", full.names = FALSE)
          } else character(0)
          resp <- jsonlite::toJSON(list(configs = files), auto_unbox = TRUE)
          if (length(files) == 0) resp <- '{"configs":[]}'
          return(list(status = 200L,
                      headers = c(cors, list("Content-Type" = "application/json")),
                      body = resp))
        }

        # API: load a saved config (returns raw YAML text)
        if (path == "/api/load-config" && method == "GET") {
          qs <- req$QUERY_STRING %||% ""
          fname <- sub(".*filename=([^&]+).*", "\\1", qs)
          if (!nzchar(fname) || grepl("/|\\\\|\\.\\.", fname)) {
            return(list(status = 400L, headers = cors,
                        body = '{"error":"Invalid filename"}'))
          }
          config_path <- file.path(project_dir, "config", utils::URLdecode(fname))
          if (!file.exists(config_path)) {
            return(list(status = 404L, headers = cors,
                        body = '{"error":"Config not found"}'))
          }
          yaml_text <- paste(readLines(config_path, warn = FALSE), collapse = "\n")
          resp <- jsonlite::toJSON(list(yaml = yaml_text), auto_unbox = TRUE)
          return(list(status = 200L,
                      headers = c(cors, list("Content-Type" = "application/json")),
                      body = resp))
        }

        # API: save config (write YAML to config/)
        if (path == "/api/save-config" && method == "POST") {
          body_raw <- req$rook.input$read()
          body <- jsonlite::fromJSON(rawToChar(body_raw))

          yaml_content <- body$yaml
          filename <- body$filename %||% "config.yaml"

          config_dir <- file.path(project_dir, "config")
          dir.create(config_dir, showWarnings = FALSE)
          config_path <- file.path(config_dir, filename)
          writeLines(yaml_content, config_path)

          cat(sprintf("  Config saved: %s\n", config_path))

          resp <- jsonlite::toJSON(list(
            path = config_path,
            message = paste("Config saved to", config_path)
          ), auto_unbox = TRUE)

          return(list(status = 200L,
                      headers = c(cors, list("Content-Type" = "application/json")),
                      body = resp))
        }

        # API: run pipeline (launch as background child, keep server alive for log streaming)
        if (path == "/api/run-pipeline" && method == "POST") {
          body_raw <- req$rook.input$read()
          body <- jsonlite::fromJSON(rawToChar(body_raw))

          yaml_content <- body$yaml
          filename <- body$filename %||% "config.yaml"
          fresh <- isTRUE(body$fresh)

          config_dir <- file.path(project_dir, "config")
          dir.create(config_dir, showWarnings = FALSE)
          config_path <- file.path(config_dir, filename)
          writeLines(yaml_content, config_path)

          # Reject if a run is already in progress (keeps UI state sane).
          if (!is.null(.wizard_state$child) && .wizard_state$child$is_alive()) {
            return(list(status = 409L,
                        headers = c(cors, list("Content-Type" = "application/json")),
                        body = '{"error":"A pipeline run is already in progress."}'))
          }

          # Fresh log for this run.
          logs_dir <- file.path(project_dir, "logs")
          dir.create(logs_dir, showWarnings = FALSE)
          ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
          log_path <- file.path(logs_dir, sprintf("wizard_run_%s.log", ts))
          file.create(log_path)

          if (!requireNamespace("processx", quietly = TRUE)) {
            return(list(status = 500L,
                        headers = c(cors, list("Content-Type" = "application/json")),
                        body = '{"error":"processx package missing. Run install.bat."}'))
          }

          r_args <- c("--vanilla", "run.R", "--config", config_path)
          if (isTRUE(fresh)) r_args <- c(r_args, "--fresh")

          child <- processx::process$new(
            command = "Rscript",
            args    = r_args,
            stdout  = log_path,
            stderr  = "2>&1",
            wd      = project_dir,
            cleanup = TRUE
          )

          .wizard_state$child       <- child
          .wizard_state$log_path    <- log_path
          .wizard_state$config_path <- config_path

          cat(sprintf("  Pipeline launched (pid=%d), logging to %s\n",
                      child$get_pid(), log_path))

          resp <- jsonlite::toJSON(list(
            message  = "Pipeline started",
            log_url  = "/api/log",
            status_url = "/api/run-status",
            pid      = child$get_pid()
          ), auto_unbox = TRUE)

          return(list(status = 200L,
                      headers = c(cors, list("Content-Type" = "application/json")),
                      body = resp))
        }

        # API: tail the current run's log from a byte offset
        if (path == "/api/log" && method == "GET") {
          lp <- .wizard_state$log_path
          if (is.null(lp) || !file.exists(lp)) {
            return(list(status = 200L,
                        headers = c(cors, list("Content-Type" = "application/json")),
                        body = '{"text":"","offset":0,"done":false,"alive":false}'))
          }
          qs <- req$QUERY_STRING %||% ""
          off <- 0L
          m <- regmatches(qs, regexpr("offset=\\d+", qs))
          if (length(m) > 0) off <- as.integer(sub("offset=", "", m))

          sz <- file.info(lp)$size
          if (is.na(sz)) sz <- 0L
          text <- ""
          if (sz > off) {
            con <- file(lp, "rb")
            on.exit(try(close(con), silent = TRUE), add = TRUE)
            seek(con, where = off, origin = "start")
            raw_bytes <- readBin(con, what = "raw", n = sz - off)
            text <- rawToChar(raw_bytes)
          }
          alive <- !is.null(.wizard_state$child) && .wizard_state$child$is_alive()
          resp <- jsonlite::toJSON(list(
            text   = text,
            offset = sz,
            alive  = alive,
            done   = !alive
          ), auto_unbox = TRUE)
          return(list(status = 200L,
                      headers = c(cors, list("Content-Type" = "application/json")),
                      body = resp))
        }

        # API: run status (alive + exit code + report path if we can find one)
        if (path == "/api/run-status" && method == "GET") {
          ch <- .wizard_state$child
          if (is.null(ch)) {
            return(list(status = 200L,
                        headers = c(cors, list("Content-Type" = "application/json")),
                        body = '{"alive":false,"started":false}'))
          }
          alive <- ch$is_alive()
          exit_code <- NA
          if (!alive) {
            exit_code <- tryCatch(ch$get_exit_status(), error = function(e) NA)
          }

          # Best-effort locate a report HTML in outputs dir if run finished
          report_path <- NA_character_
          if (!alive && !is.null(.wizard_state$config_path)) {
            cfg <- tryCatch(yaml::read_yaml(.wizard_state$config_path),
                            error = function(e) NULL)
            if (!is.null(cfg)) {
              out_base <- file.path(cfg$project$dir %||% project_dir,
                                    cfg$paths$out %||% "outputs")
              hits <- list.files(out_base, pattern = "report_.*\\.html$",
                                 recursive = TRUE, full.names = TRUE)
              if (length(hits) > 0) {
                report_path <- hits[which.max(file.mtime(hits))]
              }
            }
          }

          resp <- jsonlite::toJSON(list(
            started   = TRUE,
            alive     = alive,
            exit_code = if (is.na(exit_code)) NULL else exit_code,
            report    = if (is.na(report_path)) NULL else report_path
          ), auto_unbox = TRUE, null = "null")
          return(list(status = 200L,
                      headers = c(cors, list("Content-Type" = "application/json")),
                      body = resp))
        }

        # API: shutdown the wizard server (called by client when user clicks Close)
        if (path == "/api/shutdown" && method == "POST") {
          later::later(function() { httpuv::stopAllServers() }, 0.3)
          return(list(status = 200L,
                      headers = c(cors, list("Content-Type" = "application/json")),
                      body = '{"message":"shutting down"}'))
        }

        # 404 for anything else
        list(status = 404L, headers = cors, body = "Not found")
      }
    )

    # Start server
    server <- httpuv::startServer("127.0.0.1", port, app)
    url <- sprintf("http://localhost:%d", port)
    cat(sprintf("  Wizard running at: %s\n", url))
    cat("  Press Ctrl+C to stop.\n\n")

    # Try to open browser
    opened <- FALSE
    if (nzchar(Sys.which("wslview"))) {
      tryCatch({ system2("wslview", url, wait = FALSE); opened <- TRUE },
               error = function(e) NULL)
    }
    if (!opened) {
      tryCatch({ browseURL(url); opened <- TRUE },
               error = function(e) NULL)
    }
    if (!opened) {
      cat(sprintf("  Open in your browser: %s\n\n", url))
    }

    # Serve until stopped
    tryCatch(
      httpuv::service(Inf),
      interrupt = function(e) {
        cat("\n  Wizard stopped.\n")
      }
    )
    httpuv::stopAllServers()

    # If an older-flow run_config was queued (pre-background runs), honor it.
    if (!is.null(.wizard_state$run_config)) {
      cat("\n")
      run_pipeline(.wizard_state$run_config,
                   fresh = isTRUE(.wizard_state$run_fresh))
    }

    # If a background child is still running when the server stops, wait for it.
    if (!is.null(.wizard_state$child) && .wizard_state$child$is_alive()) {
      cat("\n  Pipeline still running in background. Waiting for it to finish...\n")
      cat(sprintf("  Live log: %s\n", .wizard_state$log_path))
      .wizard_state$child$wait()
      cat(sprintf("  Pipeline exited with code %s.\n",
                  .wizard_state$child$get_exit_status()))
    }

    return(invisible(NULL))

  } else if ("--new" %in% args) {
    check_dependencies("cli")
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
        if (is.na(idx) || idx < 1 || idx > length(configs)) {
          cat("  Invalid selection. Using config 1.\n")
          idx <- 1
        }
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
