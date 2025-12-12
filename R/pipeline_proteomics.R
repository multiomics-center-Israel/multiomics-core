#' Full proteomics preprocessing pipeline (DIA-NN)
#'
#' @param config full config list from load_config()
#'
#' @return list with:
#'   - expr_raw   : raw expr_mat from DIA-NN (log2)
#'   - expr_filt  : filtered matrix after min_count
#'   - expr_imp   : imputed matrix
#'   - row_data   : protein annotations
#'   - meta       : metadata table
run_proteomics_preprocessing <- function(config) {
  cfg    <- config$modes$proteomics
  inputs <- load_proteomics_inputs(config)
  
  # validate inputs
  validate_proteomics_inputs(inputs, cfg)
  
  # DIA-NN → expression matrix
  expr_mat <- get_measurements_per_sample_diann(
    protein    = inputs$protein,
    sample_map = inputs$sample_map,
    meta       = inputs$metadata,
    cfg        = cfg
  )
  
  # protein annotation
  row_data <- inputs$protein[
    ,
    c(cfg$id_columns$protein_id,
      unlist(cfg$id_columns$protein_annot)),
    drop = FALSE
  ]
  
  # filtering
  filt <- filter_proteomics_by_min_count(
    expr_mat = expr_mat,
    row_data = row_data,
    meta     = inputs$metadata,
    cfg      = cfg
  )
  
  expr_filt  <- filt$expr_mat
  row_data_f <- filt$row_data
  
  # imputation (currently: Perseus)
  expr_imp <- impute_proteomics_perseus(
    expr_mat = expr_filt,
    cfg      = cfg
  )
  
  list(
    expr_raw  = expr_mat,
    expr_filt = expr_filt,
    expr_imp  = expr_imp,
    row_data  = row_data_f,
    meta      = inputs$metadata
  )
}
