library(targets)

# ------------------------------------------------------------------------------
# Source R files in a strict dependency order:
# 1) core    – generic utilities, no domain knowledge
# 2) domain  – omics-specific logic (proteomics / rnaseq)
# 3) modules – target-ready wrappers (qc, de, clustering, etc.)
# 4) pipeline– targets orchestration only
# ------------------------------------------------------------------------------

core_files <- sort(list.files(
  "R/core",
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

# Source core → domain → modules
invisible(lapply(c(core_files, domain_files, module_files), tar_source))

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

tar_option_set(
  packages = c(
    "limma", "dplyr", "yaml", "pheatmap", "cluster", "ggplot2",
    "openxlsx", "readr", "tidyr", "tibble",
    "edgeR", "DESeq2", "SummarizedExperiment"
  )
)

# ------------------------------------------------------------------------------
# Targets definition
# ------------------------------------------------------------------------------

list(
  # Configuration file (tracked as a file dependency)
  tar_target(
    config_file,
    "C:/Users/sharabmi/Documents/BGU/MultiOmics/projects/01_Aaron_Fait/config.yaml",
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
      config_path = config_file
    ),
    format = "file"
  ),

  # Proteomics pipeline (returns a list of targets)
  # pipe_proteomics()

  # RNA-seq pipeline (enable when ready)
  pipe_rnaseq()
)
