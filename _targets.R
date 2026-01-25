library(targets)

r_files <- sort(list.files("R", pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE))
invisible(lapply(r_files, source))

tar_option_set(
  packages = c(
    "limma", "dplyr", "yaml", "pheatmap", "cluster", "ggplot2",
    "openxlsx", "readr", "tidyr", "tibble",
    "edgeR", "DESeq2", "SummarizedExperiment"
  )
)

list(
  tar_target(config_file, "C:/Users/sharabmi/Documents/BGU/MultiOmics/projects/02_Ella_Pick/config.yaml", format = "file"),
  tar_target(config, load_config(config_file)),
  tar_target(cfg_validated, { validate_config(config); TRUE }),
  
  tar_target(run_dir, get_run_out_dir(config)),
  tar_target(execution_info_files, write_execution_info(config = config, run_dir = run_dir), format = "file"),
  
  pipe_proteomics()
  # pipe_rnaseq()
  
)
