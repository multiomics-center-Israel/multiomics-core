library(targets)

# Load all project R functions
r_files <- list.files("R", pattern = "\\.[Rr]$", full.names = TRUE, recursive = TRUE)
invisible(lapply(r_files, source))

tar_option_set(
  packages = c("limma", "dplyr", "yaml")
)

list(
  
  tar_target(
    config_file,
    "config/config.yaml",
    format = "file"
  ),
  
  tar_target(
    config,
    load_config(config_file)
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
    prot_pre_legacy,
    as_legacy_pre(prot_pre)
  ),
  
  # Compute multiple imputations ONCE (cached)
  tar_target(
    prot_imputations,
    make_imputations_proteomics(
      expr_mat   = prot_pre$expr_filt_mat,
      cfg        = config,
      seed_base  = 1234,
      verbose    = TRUE
    )
  ),
  
  # Run limma on each imputed dataset
  tar_target(
    prot_limma_runs,
    run_limma_multimp(
      imputations  = prot_imputations,
      meta         = prot_pre$meta,
      contrasts_df = prot_inputs$contrasts,
      prot_tbl     = prot_inputs$protein,
      cfg          = config,
      verbose      = TRUE
    )
  ),
  
  # Summarize DE across imputations (legacy-compatible)
  tar_target(
    prot_de_res,
    {
      runs_de_tables <- lapply(prot_limma_runs, function(x) x$de_tables)
      
      summary_df <- summarize_limma_mult_imputation(
        runs_de_tables = runs_de_tables,
        config         = config
      )
      
      list(
        runs_de_tables = runs_de_tables,   #
        runs           = prot_limma_runs,  # nice to keep
        summary_df     = summary_df
      )
    }
  ),
  
  tar_target(
    project_rdata_file,
    write_project_rdata(
      run_dir     = prot_run_dir
    ),
    format = "file"
  ), 
  
  tar_target(
    prot_run_dir,
    get_run_out_dir(config)
  ),
  
  tar_target(
    execution_info_files,
    write_execution_info(config = config, run_dir = prot_run_dir),
    format = "file"
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
      pre        = prot_pre,
      de_res     = prot_de_res,
      inputs     = prot_inputs,
      config     = config,
      run_dir    = prot_run_dir,
      write_runs = FALSE
    ),
    format = "file"
  )
)
