# R/pipeline/mod_proteomics_de.R

#' Proteomics DE module (orchestration only)
#'
#' Runs DE according to cfg$de$method (currently supports "limma").
#' No file I/O here.
#'
#' @param pre     output of preprocess_proteomics()
#' @param inputs  output of load_proteomics_inputs()
#' @param config  full config
#' @param verbose logical
#' @return list(method, runs, runs_de_tables, summary_df)
mod_proteomics_de <- function(pre, inputs, config, verbose = FALSE) {
  assert_pre_contract(pre, stage = "proteomics")
  
  cfg <- config$modes$proteomics
  
  # ---- choose DE method ----
  method <- cfg$de$method %||% "limma"
  
  if (identical(method, "limma")) {
    
    # 1) imputations
    imputations <- make_imputations_proteomics(
      expr_mat = pre$expr_filt,
      cfg      = config,
      verbose  = verbose
    )
    
    # 2) run limma on imputations
    runs <- run_limma_multimp(
      imputations  = imputations,
      meta         = pre$meta,
      contrasts_df = inputs$contrasts,
      prot_tbl     = inputs$protein,
      cfg          = config,
      verbose      = verbose
    )
    
    # 3) collect DE tables + summarize
    runs_de_tables <- lapply(runs, function(x) x$de_tables)
    
    summary_df <- summarize_limma_mult_imputation(
      runs_de_tables = runs_de_tables,
      config         = config
    )
    
    return(list(
      method        = "limma",
      imputations   = imputations,     # אם לא בא לך להחזיר, אפשר למחוק
      runs          = runs,
      runs_de_tables = runs_de_tables,
      summary_df    = summary_df
    ))
  }
  
  stop("Unsupported proteomics DE method: ", method)
}
