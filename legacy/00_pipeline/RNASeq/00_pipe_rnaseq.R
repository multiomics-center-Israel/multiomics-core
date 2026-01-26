# 00_pipe_rnaseq.R
pipe_rnaseq <- function() {
  list(
    tar_target(rna_inputs, load_rna_inputs(config)),
    tar_target(rna_pre, preprocess_rna(rna_inputs, config, gene_lengths = NULL, verbose = TRUE)),
    
    tar_target(rna_out_dir, get_mode_out_dir(run_dir, "rna")),
    
    tar_target(
      rna_qc_pre_obj,
      mod_rnaseq_qc_pre(
        pre     = rna_pre,
        config  = config,
        out_dir = rna_out_dir
      )
    ),
    tar_target(rna_de_res, mod_rnaseq_de(rna_pre, rna_inputs, config, verbose = TRUE)),
    
    tar_target(
      rna_de_files,
      write_rnaseq_outputs_legacy(
        pre = rna_pre,
        de_res = rna_de_res, 
        inputs = rna_inputs,
        config = config,
        out_dir = rna_out_dir),
        
      format = "file"
      ),
    
   
    tar_target(
      rna_qc_post_obj,
      mod_rnaseq_qc_post(
        pre     = rna_pre,
        de_res  = rna_de_res,
        config  = config,
        out_dir = rna_out_dir
      )
    )
  )
  
}


