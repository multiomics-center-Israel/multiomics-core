library(targets)

r_files <- list.files("R", pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE)
invisible(lapply(r_files, source))

tar_option_set(
  packages = c("limma", "dplyr", "yaml", "pheatmap", "cluster", "ggplot2", "openxlsx", "readr", "tidyr", "tibble")
)

list(
  tar_target(config_file, "config/config.yaml", format = "file"),
  tar_target(config, load_config(config_file)),
  tar_target(cfg_validated, { validate_config(config); TRUE }),
  
  tar_target(prot_run_dir, get_run_out_dir(config)),
  tar_target(execution_info_files, write_execution_info(config = config, run_dir = prot_run_dir), format = "file"),
  
  pipe_proteomics()
)
