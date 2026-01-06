# R/pipeline/pipe_proteomics.R
pipe_proteomics <- function() {
  list(
    tar_target(prot_inputs, load_proteomics_inputs(config)),
    tar_target(prot_pre, preprocess_proteomics(prot_inputs, config)),

    tar_target(
      prot_imputations,
      make_imputations_proteomics(
        expr_mat = prot_pre$expr_filt,
        cfg      = config,
        verbose  = TRUE
      )
    ),
    
    tar_target(
      prot_de_runs,
      run_limma_multimp(
        imputations  = prot_imputations,
        meta         = prot_pre$meta,
        contrasts_df = prot_inputs$contrasts,
        prot_tbl     = prot_inputs$protein,
        cfg          = config,
        verbose      = TRUE
      )
    ),
    
    tar_target(
      prot_de_res,
      build_proteomics_de_results(prot_de_runs, config)
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
      prot_qc_obj,
      mod_proteomics_qc(
        pre     = prot_pre,
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
        write_runs = TRUE
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
        imputations  = prot_imputations,
        de_results   = prot_de_res,
        qc_results   = prot_qc_obj,
        clustering_results = prot_clustering_obj
      ),
      format = "file"
    )
  )
}
