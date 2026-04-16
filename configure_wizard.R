#!/usr/bin/env Rscript
# ============================================================
# Configure wizard defaults for a specific installation.
# Run via: Rscript configure_wizard.R
#       or: double-click configure_wizard.bat
# ============================================================

defaults_file <- file.path(getwd(), "wizard_defaults.json")

# --- Helpers ------------------------------------------------------------------

prompt_read <- function(prompt_text) {
  cat(prompt_text)
  if (interactive()) {
    ans <- readline("")
  } else {
    con <- file("stdin", "r")
    ans <- readLines(con, n = 1, warn = FALSE)
    close(con)
    if (length(ans) == 0) return("")
  }
  trimws(ans)
}

ask <- function(prompt, default = NULL) {
  if (!is.null(default) && nzchar(default)) {
    prompt_text <- paste0("  ", prompt, " [", default, "]: ")
  } else {
    prompt_text <- paste0("  ", prompt, ": ")
  }
  ans <- prompt_read(prompt_text)
  if (ans == "" && !is.null(default)) return(default)
  ans
}

ask_yn <- function(prompt, default = TRUE) {
  hint <- if (default) "Y/n" else "y/N"
  ans <- prompt_read(paste0("  ", prompt, " [", hint, "]: "))
  if (ans == "") return(default)
  tolower(ans) %in% c("y", "yes")
}

ask_choice <- function(prompt, choices, default = 1) {
  cat(paste0("  ", prompt, "\n"))
  for (i in seq_along(choices)) {
    marker <- if (i == default) " *" else ""
    cat(sprintf("    %d) %s%s\n", i, choices[i], marker))
  }
  ans <- prompt_read(sprintf("  Choose [%d]: ", default))
  if (ans == "") return(default)
  idx <- suppressWarnings(as.integer(ans))
  if (is.na(idx) || idx < 1 || idx > length(choices)) return(default)
  idx
}

# --- Load existing defaults ---------------------------------------------------

template_file <- file.path(getwd(), "wizard_defaults.template.json")

if (file.exists(defaults_file)) {
  cfg <- jsonlite::fromJSON(defaults_file, simplifyVector = TRUE)
  cat("  Loaded existing configuration.\n\n")
} else if (file.exists(template_file)) {
  cfg <- jsonlite::fromJSON(template_file, simplifyVector = TRUE)
  cat("  Loaded defaults from template.\n\n")
} else {
  cfg <- list(
    theme = "light",
    available_modes = c("rna", "proteomics", "metabolomics", "lipidomics", "multi", "combination"),
    default_mode = "rna",
    defaults = list(
      analyst = "Bioinformatics Core",
      round = "A01",
      organism = "human",
      generate_pptx = TRUE,
      commentary_enabled = FALSE,
      commentary_backend = "claude"
    ),
    mode_defaults = list(
      rna = list(),
      proteomics = list(),
      metabolomics = list(),
      lipidomics = list()
    )
  )
}

# --- Interactive Configuration ------------------------------------------------

cat("--- Theme ---\n")
theme_idx <- ask_choice("Default theme:",
  c("Light (recommended for lab computers)", "Dark"),
  default = if (cfg$theme == "dark") 2 else 1)
cfg$theme <- c("light", "dark")[theme_idx]

cat("\n--- Available Analysis Modes ---\n")
cat("  Which modes should be available in the wizard?\n")
all_modes <- c("rna", "proteomics", "metabolomics", "lipidomics", "multi", "combination")
mode_labels <- c("RNA-seq", "Proteomics", "Metabolomics", "Lipidomics", "Multiomics", "Combination")
active_modes <- character(0)
for (i in seq_along(all_modes)) {
  is_active <- all_modes[i] %in% cfg$available_modes
  if (ask_yn(sprintf("  Enable %s?", mode_labels[i]), is_active)) {
    active_modes <- c(active_modes, all_modes[i])
  }
}
if (length(active_modes) == 0) {
  cat("  WARNING: No modes selected, defaulting to all.\n")
  active_modes <- all_modes
}
cfg$available_modes <- active_modes

cat("\n--- Default Mode ---\n")
avail_labels <- mode_labels[all_modes %in% active_modes]
avail_modes <- all_modes[all_modes %in% active_modes]
default_idx <- which(avail_modes == cfg$default_mode)
if (length(default_idx) == 0) default_idx <- 1
default_idx <- ask_choice("Which mode opens by default?", avail_labels, default = default_idx)
cfg$default_mode <- avail_modes[default_idx]

cat("\n--- General Defaults ---\n")
cfg$defaults$analyst <- ask("Default analyst name", cfg$defaults$analyst)
cfg$defaults$round <- ask("Default analysis round", cfg$defaults$round)

org_choices <- c("human", "mouse", "rat", "other")
org_idx <- which(org_choices == cfg$defaults$organism)
if (length(org_idx) == 0) org_idx <- 1
org_idx <- ask_choice("Default organism:", c("Human", "Mouse", "Rat", "Other"), default = org_idx)
cfg$defaults$organism <- org_choices[org_idx]

cfg$defaults$generate_pptx <- ask_yn("Generate PowerPoint by default?",
  isTRUE(cfg$defaults$generate_pptx))
cfg$defaults$commentary_enabled <- ask_yn("Enable AI commentary by default?",
  isTRUE(cfg$defaults$commentary_enabled))

# --- Mode-specific defaults ---------------------------------------------------

if ("proteomics" %in% active_modes) {
  cat("\n--- Proteomics Defaults ---\n")
  prot <- cfg$mode_defaults$proteomics
  if (is.null(prot)) prot <- list()

  de_choices <- c("limma", "ttest")
  de_idx <- which(de_choices == prot$prot_de_method)
  if (length(de_idx) == 0) de_idx <- 1
  de_idx <- ask_choice("DE method:", c("limma (recommended)", "t-test"), default = de_idx)
  prot$prot_de_method <- de_choices[de_idx]

  imp_choices <- c("perseus", "dep2", "qrilc", "minval", "none")
  imp_labels <- c("Perseus-style (recommended)", "DEP2 / MinDet", "QRILC (DEP workflow)", "MinVal (floor of minimum)", "None (complete cases)")
  imp_idx <- which(imp_choices == prot$prot_imp_method)
  if (length(imp_idx) == 0) imp_idx <- 1
  imp_idx <- ask_choice("Imputation method:", imp_labels, default = imp_idx)
  prot$prot_imp_method <- imp_choices[imp_idx]

  path_choices <- c("both", "fgsea", "ora", "none")
  path_labels <- c("Both fGSEA + ORA (recommended)", "fGSEA only", "ORA only", "None")
  path_idx <- which(path_choices == prot$prot_pathway)
  if (length(path_idx) == 0) path_idx <- 1
  path_idx <- ask_choice("Pathway enrichment:", path_labels, default = path_idx)
  prot$prot_pathway <- path_choices[path_idx]

  prot$prot_clustering <- ask_yn("Enable clustering?", isTRUE(if (is.null(prot$prot_clustering)) TRUE else prot$prot_clustering))
  prot$prot_ppi <- ask_yn("Enable PPI networks?", isTRUE(if (is.null(prot$prot_ppi)) FALSE else prot$prot_ppi))

  cfg$mode_defaults$proteomics <- prot
}

if ("rna" %in% active_modes) {
  cat("\n--- RNA-seq Defaults ---\n")
  rna <- cfg$mode_defaults$rna
  if (is.null(rna)) rna <- list()

  de_choices <- c("deseq2", "edger", "limma")
  de_idx <- which(de_choices == rna$rna_de_method)
  if (length(de_idx) == 0) de_idx <- 1
  de_idx <- ask_choice("DE method:", c("DESeq2 (recommended)", "edgeR", "limma-voom"), default = de_idx)
  rna$rna_de_method <- de_choices[de_idx]

  cfg$mode_defaults$rna <- rna
}

# --- Save ---------------------------------------------------------------------

cat("\n--- Summary ---\n")
cat(sprintf("  Theme:          %s\n", cfg$theme))
cat(sprintf("  Available modes: %s\n", paste(cfg$available_modes, collapse = ", ")))
cat(sprintf("  Default mode:   %s\n", cfg$default_mode))
cat(sprintf("  Analyst:        %s\n", cfg$defaults$analyst))
cat(sprintf("  Organism:       %s\n", cfg$defaults$organism))
cat(sprintf("  PPTX:           %s\n", if (cfg$defaults$generate_pptx) "yes" else "no"))
cat(sprintf("  AI commentary:  %s\n", if (cfg$defaults$commentary_enabled) "yes" else "no"))

# Remove comment fields before writing
cfg[["_comment"]] <- NULL
cfg[["_comment2"]] <- NULL

json_out <- jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE)
writeLines(json_out, defaults_file)

cat(sprintf("\n  Saved to: %s\n", defaults_file))
cat("  The wizard will use these defaults on next launch.\n")
