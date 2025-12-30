# R/core/proteomics_adapter.R

as_legacy_pre <- function(pre) {
  pre$expr_filt <- pre$expr_filt_mat %||% pre$expr_filt
  pre$expr_imp  <- pre$expr_imp_single %||% pre$expr_imp
  
  if (is.null(pre$imputation)) pre$imputation <- list()
  if (is.null(pre$imputation$imputed_flag) && !is.null(pre$imputation_qc$imputed_flag)) {
    pre$imputation$imputed_flag <- pre$imputation_qc$imputed_flag
  }
  
  pre
}
