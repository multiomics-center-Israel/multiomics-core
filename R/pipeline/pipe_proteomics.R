# R/pipeline/pipe_proteomics.R
pipe_proteomics <- function() {
  list(
    tar_target(prot_inputs, load_proteomics_inputs(config)),
    tar_target(prot_pre, preprocess_proteomics(prot_inputs, config)),

    tar_target(
      prot_qc_pre_obj,
      mod_proteomics_qc_pre(
        pre     = prot_pre,
        config  = config,
        run_dir = prot_run_dir
      )
    ),
    
    
    tar_target(
      prot_de_res,
      mod_proteomics_de(
        pre     = prot_pre,
        inputs  = prot_inputs,
        config  = config,
        verbose = TRUE
      )
    ),
    
    tar_target(
      prot_qc_post_obj,
      mod_proteomics_qc_post(
        pre          = prot_pre,
        de_res       = prot_de_res,
        config       = config,
        run_dir      = prot_run_dir,
        de_source    = "table1"
      )
    ),

    tar_target(
      prot_clustering_obj,
      mod_proteomics_clustering(
        pre     = prot_pre,
        de_res  = prot_de_res,
        config  = config,
        run_dir = prot_run_dir
      )
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
    ),

    tar_target(
      project_rdata_file,
      write_project_rdata(
        run_dir      = prot_run_dir,
        config       = config,
        inputs       = prot_inputs,
        pre_process  = prot_pre,
        imputations  = prot_de_res$imputations,  #  
        de_results   = prot_de_res,
        qc_results   = prot_qc_pre_obj
      ),
      format = "file"
    )
  )
}
