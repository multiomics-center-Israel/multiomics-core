#' Mark differential expression pass for a single limma result table
#'
#' Applies the legacy "pass1" rule: a feature passes if it meets BOTH
#' the p-value cutoff (raw or adjusted) and the absolute log2FC cutoff.
#'
#' Output encoding matches the old script:
#'   - 1  : pass
#'   - NA : fail
#'
#' @param de_table A data.frame from limma::topTable() containing at least
#'        'logFC', 'P.Value', and 'adj.P.Val'.
#' @param p_cutoff Numeric cutoff for p-value or adjusted p-value.
#' @param lfc_cutoff Numeric cutoff for |logFC| (log2 scale).
#' @param use_adj Logical; TRUE uses 'adj.P.Val', FALSE uses 'P.Value'.
#'
#' @return Numeric vector of length nrow(de_table) with values {1, NA}.
mark_pass1 <- function(de_table, p_cutoff, lfc_cutoff, use_adj = TRUE) {
  pcol <- if (isTRUE(use_adj)) "adj.P.Val" else "P.Value"
  stopifnot(all(c("logFC", pcol) %in% colnames(de_table)))
  
  pass <- (de_table[[pcol]] <= p_cutoff) &
    (abs(de_table[["logFC"]]) >= lfc_cutoff)
  
  ifelse(pass, 1, NA)
}


#' Summarize differential expression across multiple imputations (legacy-compatible)
#'
#' Reproduces the old script logic:
#'   - "stable" if passed in >= min_no_passed imputations
#'   - average linear ratio across imputations and convert to signed linearFC
#'   - summarize p-value and padj using quantile at (min_no_passed / no_repetitions)
#'
#' @param runs_de_tables List of length cfg$no_repetitions.
#'        Each element corresponds to one imputation run and contains a named list
#'        of DE tables (one per contrast), as returned by run_limma_proteomics()$de_tables.
#' @param config 
#'
#' @return A data.frame with one row per feature, including annotations and per-contrast
#'         summary columns (sum.pass.*, pass.imputs.*, linearRatio.imputs.*, linearFC.imputs.*,
#'         pvalue.imputs.*, padj.imputs.*).
summarize_limma_mult_imputation <- function(runs_de_tables, config) {
  limma_cfg <- config$modes$proteomics$limma
  imp_cfg   <- config$modes$proteomics$imputation
  
  NO_REPETITIONS <- as.integer(imp_cfg$no_repetitions)
  MIN_NO_PASSED  <- as.integer(imp_cfg$min_no_passed)
  
  use_adj_for_pass1 <- isTRUE(limma_cfg$use_adj_for_pass1)
  p_cutoff <- as.numeric(limma_cfg$p_cutoff)
  
  linear_fc_cutoff <- as.numeric(limma_cfg$linear_fc_cutoff)
  lfc_cutoff <- log2(linear_fc_cutoff)
  
  stopifnot(length(runs_de_tables) == NO_REPETITIONS)
  stopifnot(MIN_NO_PASSED >= 1, MIN_NO_PASSED <= NO_REPETITIONS)
  
  q <- MIN_NO_PASSED / NO_REPETITIONS
  
  contrasts <- names(runs_de_tables[[1]])
  stopifnot(length(contrasts) > 0)
  
  id_cols <- c("FeatureID", "Protein.Names", "Genes", "First.Protein.Description")
  ref_df  <- runs_de_tables[[1]][[contrasts[1]]]
  stopifnot(all(id_cols %in% colnames(ref_df)))
  
  out <- ref_df[, id_cols, drop = FALSE]
  ref_ids <- ref_df$FeatureID
  
  for (n in seq_len(NO_REPETITIONS)) {
    for (cn in contrasts) {
      cur <- runs_de_tables[[n]][[cn]]
      stopifnot(nrow(cur) == length(ref_ids))
      stopifnot(all(cur$FeatureID == ref_ids))
    }
  }
  
  for (cn in contrasts) {
    contrast_print <- gsub(" ", "", cn)
    
    logfc_mat <- sapply(seq_len(NO_REPETITIONS),
                        function(n) runs_de_tables[[n]][[cn]][["logFC"]])
    p_mat     <- sapply(seq_len(NO_REPETITIONS),
                        function(n) runs_de_tables[[n]][[cn]][["P.Value"]])
    padj_mat  <- sapply(seq_len(NO_REPETITIONS),
                        function(n) runs_de_tables[[n]][[cn]][["adj.P.Val"]])
    
    pass1_mat <- sapply(seq_len(NO_REPETITIONS), function(n)
      mark_pass1(runs_de_tables[[n]][[cn]],
                 p_cutoff    = p_cutoff,
                 lfc_cutoff  = lfc_cutoff,
                 use_adj     = use_adj_for_pass1)
    )
    
    sum_pass    <- rowSums(pass1_mat, na.rm = TRUE)
    pass_imputs <- ifelse(sum_pass >= MIN_NO_PASSED, 1, NA)
    
    linearRatio_imputs <- rowMeans(2^logfc_mat, na.rm = TRUE)
    linearFC_imputs <- ifelse(linearRatio_imputs >= 1,
                              linearRatio_imputs,
                              -1 / linearRatio_imputs)
    
    pvalue_imputs <- apply(p_mat,   1, quantile, probs = q, na.rm = TRUE)
    padj_imputs   <- apply(padj_mat,1, quantile, probs = q, na.rm = TRUE)
    
    out[[paste0("sum.pass.", contrast_print)]]           <- sum_pass
    out[[paste0("pass.imputs.", contrast_print)]]        <- pass_imputs
    out[[paste0("linearRatio.imputs.", contrast_print)]] <- linearRatio_imputs
    out[[paste0("linearFC.imputs.", contrast_print)]]    <- signif(linearFC_imputs, 4)
    out[[paste0("pvalue.imputs.", contrast_print)]]      <- pvalue_imputs
    out[[paste0("padj.imputs.", contrast_print)]]        <- padj_imputs
  }
  
  out
}

