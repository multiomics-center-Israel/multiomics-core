library(targets)

# ------------------------------------------------------------------------------
# Source R files in a strict dependency order:
# 1) core     – generic utilities, no domain knowledge
# 2) services – external-facing helpers (AI commentary, etc.)
# 3) domain   – omics-specific logic (proteomics / rnaseq)
# 4) modules  – target-ready wrappers (qc, de, clustering, etc.)
# 5) pipeline – targets orchestration only
# ------------------------------------------------------------------------------

core_files <- sort(list.files(
  "R/core",
  pattern = "\\.R$",
  full.names = TRUE,
  recursive = TRUE
))

service_files <- sort(list.files(
  "R/services",
  pattern = "\\.R$",
  full.names = TRUE,
  recursive = TRUE
))

domain_files <- sort(list.files(
  "R/domain",
  pattern = "\\.R$",
  full.names = TRUE,
  recursive = TRUE
))

module_files <- sort(list.files(
  "R/modules",
  pattern = "\\.R$",
  full.names = TRUE,
  recursive = TRUE
))

# Source core -> services -> domain -> modules
invisible(lapply(c(core_files, service_files, domain_files, module_files), tar_source))

# Source pipelines LAST (they depend on everything above)
pipeline_files <- sort(list.files(
  "R/pipeline",
  pattern = "\\.R$",
  full.names = TRUE,
  recursive = TRUE
))

invisible(lapply(pipeline_files, tar_source))

# ------------------------------------------------------------------------------
# Global targets options
# ------------------------------------------------------------------------------

required_pkgs <- c(
  "limma", "dplyr", "yaml", "pheatmap", "cluster", "ggplot2",
  "openxlsx", "readr", "readxl", "tidyr", "tibble",
  "edgeR", "DESeq2", "SummarizedExperiment"
)
# Only require packages that are actually installed (allows running a
# subset of pipelines when some omics-specific packages are absent).
available_pkgs <- required_pkgs[vapply(required_pkgs, requireNamespace,
                                       logical(1), quietly = TRUE)]

tar_option_set(packages = available_pkgs)

# Resolve config path once at plan-definition time so the literal path is
# baked into the target command.  This ensures {targets} detects a change
# when MULTIOMICS_CONFIG points to a different file between runs.
config_path <- Sys.getenv("MULTIOMICS_CONFIG", "config.yaml")

# ------------------------------------------------------------------------------
# Targets definition
# ------------------------------------------------------------------------------

list(
  # Configuration file (tracked as a file dependency)
  # Set via: Sys.setenv(MULTIOMICS_CONFIG = "/path/to/config.yaml")
  # Or defaults to config.yaml
  tar_target(
    config_file,
    {
      cfg_path <- Sys.getenv("MULTIOMICS_CONFIG", unset = "")
      if (cfg_path == "") {
        cfg_path <- file.path(getwd(), "config.yaml")
      }
      normalizePath(cfg_path, mustWork = TRUE)
    },
    format = "file"
  ),

  # Load configuration
  tar_target(
    config,
    load_config(config_file)
  ),

  # Validate configuration early; downstream targets should depend on this
  tar_target(
    cfg_validated,
    {
      validate_config(config)
      TRUE
    }
  ),

  # Resolve run output directory
  tar_target(
    run_dir,
    get_run_out_dir(config)
  ),

  # Write execution metadata (run info, config snapshot, etc.)
  tar_target(
    execution_info_files,
    write_execution_info(
      config = config,
      run_dir = run_dir,
      config_path = config_file,
      targets_file = "_targets.R"
    ),
    format = "file"
  ),

  # Mode-specific pipelines (only run if mode is present in config)
  {
    cfg_path <- Sys.getenv("MULTIOMICS_CONFIG", unset = "")
    if (cfg_path == "") cfg_path <- file.path(getwd(), "config.yaml")
    cfg_raw <- yaml::read_yaml(cfg_path)
    mode_targets <- list()
    if (!is.null(cfg_raw$modes$rna))           mode_targets <- c(mode_targets, pipe_rnaseq())
    if (!is.null(cfg_raw$modes$proteomics))    mode_targets <- c(mode_targets, pipe_proteomics())
    if (!is.null(cfg_raw$modes$metabolomics))  mode_targets <- c(mode_targets, pipe_metabolomics())
    mode_targets
  }
)
