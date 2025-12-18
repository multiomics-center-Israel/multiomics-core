#' Run proteomics DE with limma across multiple imputations (legacy-compatible wrapper)
#'
#' This is a legacy-compatible convenience wrapper that:
#' 1) Generates N imputed datasets (Perseus-style, controlled by config)
#' 2) Runs limma on each imputed dataset
#' 3) Returns a list of per-imputation limma results
#'
#' The output is designed to be consumed by downstream summarization, e.g.
#' `summarize_limma_mult_imputation_proteomics()` (or similar).
#'
#' @param expr_mat Numeric matrix to impute (typically filtered + normalized + log2).
#'        Rows are features (proteins), columns are samples.
#' @param meta Sample metadata (must match the samples/columns in `expr_mat`).
#' @param contrasts_df Contrasts table (e.g. Contrast_name / Factor / Numerator / Denominator).
#' @param prot_tbl Protein annotation table (must include Protein.Group and any annotation columns used downstream).
#' @param cfg Full config list. Must contain proteomics settings under:
#'        - `cfg$modes$proteomics$imputation` (including `no_repetitions` and Perseus parameters)
#'        - any other proteomics fields required by `impute_proteomics_perseus()`
#' @param seed_base Optional integer base seed for reproducibility.
#'        If provided, each imputation uses `set.seed(seed_base + i)`.
#' @param verbose Logical; if TRUE, prints progress messages.
#'
#' @return A list of length N (= `cfg$modes$proteomics$imputation$no_repetions`).
#'         Each element is the output of `run_limma_proteomics()`
#'         (typically containing DE tables per contrast).
#'
#' @examples
#' \dontrun{
#' res_runs <- run_limma_mult_imputation_proteomics(
#'   expr_mat     = expr_log2_filtered,
#'   meta         = meta,
#'   contrasts_df = contrasts_df,
#'   prot_tbl     = prot_tbl,
#'   cfg          = config,
#'   seed_base    = 1234,
#'   verbose      = TRUE
#' )
#' }
run_limma_mult_imputation_proteomics <- function(expr_mat,
                                                 meta,
                                                 contrasts_df,
                                                 prot_tbl,
                                                 cfg,
                                                 seed_base = NULL,
                                                 verbose = FALSE) {
  
  imputations <- make_imputations_proteomics(
    expr_mat   = expr_mat,
    cfg        = cfg,
    seed_base  = seed_base,
    verbose    = verbose
  )
  
  run_limma_multimp(
    imputations  = imputations,
    meta         = meta,
    contrasts_df = contrasts_df,
    prot_tbl     = prot_tbl,
    cfg          = cfg,
    verbose      = verbose
  )
}




#' Run limma proteomics DE on a precomputed list of imputed datasets
#'
#' Use this function when imputations were already generated and you want to:
#' - avoid recomputing imputations (e.g., when trying a different DE method), or
#' - decouple imputation from DE in a `{targets}` pipeline.
#'
#' @param imputations List of imputed expression matrices.
#'        Each matrix must have the same dimensions and row/column names as the original expression matrix.
#' @param meta Sample metadata (must match the samples/columns in each imputed matrix).
#' @param contrasts_df Contrasts table (e.g. Contrast_name / Factor / Numerator / Denominator).
#' @param prot_tbl Protein annotation table (must include Protein.Group and any annotation columns used downstream).
#' @param verbose Logical; if TRUE, prints progress messages.
#'
#' @return A list of length `length(imputations)`. Each element is the output of
#'         `run_limma_proteomics()` for the corresponding imputed dataset.
#'
#' @examples
#' \dontrun{
#' imputations <- make_imputations_proteomics(expr_mat, cfg = config, seed_base = 1234)
#' limma_runs  <- run_limma_multimp(
#'   imputations  = imputations,
#'   meta         = meta,
#'   contrasts_df = contrasts_df,
#'   prot_tbl     = prot_tbl
#' )
#' }
run_limma_multimp <- function(imputations,
                                                meta,
                                                contrasts_df,
                                                prot_tbl,
                                                cfg,
                                                verbose = FALSE) {
  
  validate_proteomics_imputations(imputations = imputations, meta = meta, cfg = cfg)
  
  runs <- vector("list", length(imputations))
  
  for (i in seq_along(imputations)) {
    if (isTRUE(verbose)) {
      message(sprintf("Limma on imputations: %d / %d", i, length(imputations)))
    }
    
    runs[[i]] <- run_limma_proteomics(
      expr_imp     = imputations[[i]],
      meta         = meta,
      contrasts_df = contrasts_df,
      prot_tbl     = prot_tbl,
      cfg          = cfg
    )
  }
  
  runs
}



