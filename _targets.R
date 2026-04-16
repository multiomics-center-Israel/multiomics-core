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
  "edgeR", "DESeq2", "SummarizedExperiment",
  "impute"   # Bioconductor: KNN imputation for metabolomics MAR features
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

  # Load and validate configuration. validate_config() applies defaults (e.g.
  # multiomics integration methods) and returns the updated config, so all
  # downstream targets that depend on `config` receive the defaulted values.
  tar_target(
    config,
    validate_config(load_config(config_file))
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
      config       = config,
      run_dir      = run_dir,
      config_path  = config_file,
      targets_file = "_targets.R"
    ),
    format = "file"
  ),

  # Mode-specific pipelines — only included when the mode is present in config.
  # Read config at plan-definition time so {targets} can detect mode changes.
  {
    cfg_raw      <- yaml::read_yaml(config_path)
    mode_targets <- list()

    # Single-omics pipelines
    if (!is.null(cfg_raw$modes$rna))           mode_targets <- c(mode_targets, pipe_rnaseq())
    if (!is.null(cfg_raw$modes$proteomics))    mode_targets <- c(mode_targets, pipe_proteomics())
    if (!is.null(cfg_raw$modes$metabolomics)) {
        metab_chosen_norm <- cfg_raw$modes$metabolomics$preprocessing$chosen_norm
        mode_targets <- c(mode_targets, pipe_metabolomics(chosen_norm = metab_chosen_norm))
    }
    if (!is.null(cfg_raw$modes$lipidomics))    mode_targets <- c(mode_targets, pipe_lipidomics())

    # Multi-omics integration pipeline (runs AFTER single-omics pipelines)
    # Enabled if: (a) multiomics mode is configured, AND
    #             (b) >=2 omics modes are present OR input_mode == "outputs" (payload mode)
    n_omics <- sum(!is.null(cfg_raw$modes$rna),
                   !is.null(cfg_raw$modes$proteomics),
                   !is.null(cfg_raw$modes$metabolomics))
    payload_mode <- identical(cfg_raw$modes$multiomics$input_mode, "outputs")
    n_omics_global <- length(cfg_raw$global$omics_present %||% character(0))

    if (!is.null(cfg_raw$modes$multiomics) &&
        (n_omics >= 2 || payload_mode || n_omics_global >= 2)) {
      mode_targets <- c(mode_targets, pipe_multiomics())
    }

    mode_targets
  }
)
