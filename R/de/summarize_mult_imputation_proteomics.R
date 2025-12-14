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
#' @param cfg Limma sub-config for proteomics (e.g., config$modes$proteomics$limma) with:
#'        no_repetitions, min_no_passed, use_adj_for_pass1, p_cutoff, linear_fc_cutoff.
#'
#' @return A data.frame with one row per feature, including annotations and per-contrast
#'         summary columns (sum.pass.*, pass.imputs.*, linearRatio.imputs.*, linearFC.imputs.*,
#'         pvalue.imputs.*, padj.imputs.*).
summarize_mult_imputation <- function(runs_de_tables, cfg) {
  
  # Read parameters from the limma sub-config
  NO_REPETITIONS <- as.integer(cfg$no_repetitions)
  MIN_NO_PASSED  <- as.integer(cfg$min_no_passed)
  
  use_adj_for_pass1 <- isTRUE(cfg$use_adj_for_pass1)
  p_cutoff <- as.numeric(cfg$p_cutoff)
  
  # Old script uses a linear FC cutoff but compares on log2FC scale for pass1
  linear_fc_cutoff <- as.numeric(cfg$linear_fc_cutoff)
  lfc_cutoff <- log2(linear_fc_cutoff)
  
  stopifnot(length(runs_de_tables) == NO_REPETITIONS)
  stopifnot(MIN_NO_PASSED >= 1, MIN_NO_PASSED <= NO_REPETITIONS)
  
  # Quantile probability used in the legacy code (e.g. 0.8 for 8/10)
  q <- MIN_NO_PASSED / NO_REPETITIONS
  
  contrasts <- names(runs_de_tables[[1]])
  stopifnot(length(contrasts) > 0)
  
  # Use the first run as reference for feature order and annotations
  id_cols <- c("FeatureID", "Protein.Names", "Genes", "First.Protein.Description")
  ref_df  <- runs_de_tables[[1]][[contrasts[1]]]
  stopifnot(all(id_cols %in% colnames(ref_df)))
  
  out <- ref_df[, id_cols, drop = FALSE]
  ref_ids <- ref_df$FeatureID
  
  # Ensure identical feature set and order across all runs and contrasts
  for (n in seq_len(NO_REPETITIONS)) {
    for (cn in contrasts) {
      cur <- runs_de_tables[[n]][[cn]]
      stopifnot(nrow(cur) == length(ref_ids))
      stopifnot(all(cur$FeatureID == ref_ids))
    }
  }
  
  # Summarize per contrast
  for (cn in contrasts) {
    contrast_print <- gsub(" ", "", cn)  # match legacy column naming
    
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
    
    # Stability summary (legacy logic)
    sum_pass    <- rowSums(pass1_mat, na.rm = TRUE)
    pass_imputs <- ifelse(sum_pass >= MIN_NO_PASSED, 1, NA)
    
    # Average linear ratio and signed linearFC (legacy convention)
    linearRatio_imputs <- rowMeans(2^logfc_mat, na.rm = TRUE)
    linearFC_imputs <- ifelse(linearRatio_imputs >= 1,
                              linearRatio_imputs,
                              -1 / linearRatio_imputs)
    
    # Quantile-based p-value summaries (legacy logic)
    pvalue_imputs <- apply(p_mat,   1, quantile, probs = q, na.rm = TRUE)
    padj_imputs   <- apply(padj_mat,1, quantile, probs = q, na.rm = TRUE)
    
    out[[paste0("sum.pass.", contrast_print)]]            <- sum_pass
    out[[paste0("pass.imputs.", contrast_print)]]         <- pass_imputs
    out[[paste0("linearRatio.imputs.", contrast_print)]]  <- linearRatio_imputs
    out[[paste0("linearFC.imputs.", contrast_print)]]     <- signif(linearFC_imputs, 4)
    out[[paste0("pvalue.imputs.", contrast_print)]]       <- pvalue_imputs
    out[[paste0("padj.imputs.", contrast_print)]]         <- padj_imputs
  }
  
  out
}
