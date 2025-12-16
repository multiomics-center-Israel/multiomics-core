library(targets)

# Load all project R functions
r_files <- list.files("R", pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE)
invisible(lapply(r_files, source))


tar_option_set(
  packages = c("limma", "dplyr")
)

list(
  tar_target(
    config,
    load_config("config/config.yaml")
  ),
  
  tar_target(
    cfg_validated,
    {
      validate_config(config)
      TRUE
    }
  ),
  
  tar_target(
    prot_inputs,
    load_proteomics_inputs(config)
  ),
  
  tar_target(
    prot_pre,
    preprocess_proteomics(prot_inputs, config)
  ),
  
  tar_target(
    prot_de,
    run_proteomics_de_mult_impute(
      expr_mat = prot_pre$expr_imp,
      inputs   = prot_inputs,
      config   = config
    )
  ),
  tar_target(
    prot_run_dir,
    get_run_out_dir(config)
  ),
  
  tar_target(
    prot_qc_files,
    mod_proteomics_qc(
      pre     = prot_pre,
      config  = config,
      run_dir = prot_run_dir
    ),
    format = "file"
  ),
  
  
  tar_target(
    prot_de_files,
    write_proteomics_multimpute_outputs(
      pre       = prot_pre,
      de_res    = prot_de,
      inputs    = prot_inputs,
      config    = config,
      run_dir   = get_run_out_dir(config),
      write_runs = FALSE
    ),
    format = "file"
  )
  
  
)
