#' Run proteomics imputation N times and run limma on each imputed dataset
#'
#' Legacy-compatible workflow:
#'   - Perform multiple imputations (NO_REPETITIONS)
#'   - Run limma per imputation
#'   - Return a list of per-imputation DE tables (one list per run; one table per contrast)
#'
#' This output is designed to be consumed by summarize_mult_imputation().
#'
#' @param expr_mat Matrix to impute (typically filtered + normalized + log2).
#' @param meta Sample metadata with SampleID and Condition.
#' @param contrasts_df Contrasts table (Contrast_name / Factor / Numerator / Denominator).
#' @param prot_tbl Protein annotation table (must include Protein.Group + annotation columns).
#' @param cfg Full config list OR proteomics sub-config; must contain:
#'        cfg$modes$proteomics$imputation$width/downshift and cfg$modes$proteomics$limma$no_repetitions
#'        If you prefer passing sub-configs, see notes below.
#' @param seed_base Optional integer base seed for reproducibility (seed = seed_base + i).
#'
#' @return A list of length NO_REPETITIONS. Each element is run_limma_proteomics() output,
#'         and contains $de_tables (named by contrast).
run_limma_mult_imputation_proteomics <- function(expr_mat,
                                                 meta,
                                                 contrasts_df,
                                                 prot_tbl,
                                                 cfg,
                                                 seed_base = NULL,
                                                 verbose = FALSE) {
  
  # Pull limma + imputation sub-configs (assumes your YAML structure)
  limma_cfg <- cfg$modes$proteomics$limma
  imp_cfg   <- cfg$modes$proteomics$imputation
  
  NO_REPETITIONS <- as.integer(limma_cfg$no_repetitions)
  stopifnot(NO_REPETITIONS >= 1)

  
  # Run N imputations + limma
  runs <- vector("list", NO_REPETITIONS)
  
  for (i in seq_len(NO_REPETITIONS)) {
    
    if (isTRUE(verbose)) {
      message(sprintf("Limma multi-imputation: run %d / %d", i, NO_REPETITIONS))
    }
    
    
    # Set a deterministic seed per iteration (recommended for targets reproducibility)
    if (!is.null(seed_base)) {
      set.seed(as.integer(seed_base) + i)
    }
    
    # --- Imputation step (Perseus-style) ---
    # impute_proteomics_perseus() returns a matrix by default (return_flags = FALSE)
    expr_imp_i <- impute_proteomics_perseus(expr_mat, cfg = cfg$modes$proteomics, return_flags = FALSE)
    
    stopifnot(is.matrix(expr_imp_i))
    stopifnot(all(dim(expr_imp_i) == dim(expr_mat)))
    stopifnot(all(colnames(expr_imp_i) == colnames(expr_mat)))
    stopifnot(all(rownames(expr_imp_i) == rownames(expr_mat)))
    
    # --- Limma step (single run) ---
    runs[[i]] <- run_limma_proteomics(
      expr_imp     = expr_imp_i,
      meta         = meta,
      contrasts_df = contrasts_df,
      prot_tbl     = prot_tbl
    )
  }
  
  runs
}
