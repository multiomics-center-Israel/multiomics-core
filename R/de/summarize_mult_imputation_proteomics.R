#' Mark differential expression pass for a single result table
#'
#' Applies the legacy "pass1" rule: a feature passes if it meets BOTH
#' the p-value cutoff (raw or adjusted) and the absolute log2FC cutoff.
#'
#' Output encoding matches the old script:
#'   - 1  : pass
#'   - NA : fail
#'
#' @param de_table A data.frame from topTable()
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

#' Helper to add pass_any_contrast column
#' 
#' Scans for columns starting with 'pass.imputs.' and creates a summary flag.
#'
#' @param summary_df Data frame containing pass columns.
#' @param pass_prefix Regex prefix for pass columns.
#' @param out_col Name of the output binary column.
#' @param out_n_col Name of the output count column.
#' 
#' @return The data frame with added columns.
add_pass_any_contrast <- function(summary_df,
                                  pass_prefix = "^pass\\.imputs\\.",
                                  out_col = "pass_any_contrast",
                                  out_n_col = "n_pass_contrasts") {
  
  pass_cols <- grep(pass_prefix, names(summary_df), value = TRUE)
  
  if (length(pass_cols) == 0) {
    # If no pass columns found, initialize with 0/NA to avoid downstream errors
    summary_df[[out_n_col]] <- 0L
    summary_df[[out_col]]   <- NA
    return(summary_df)
  }
  
  pass_mat <- as.matrix(summary_df[, pass_cols, drop = FALSE])
  
  # Count how many contrasts passed (ignoring NAs)
  n_pass <- rowSums(!is.na(pass_mat) & pass_mat == 1, na.rm = TRUE)
  
  summary_df[[out_n_col]] <- as.integer(n_pass)
  summary_df[[out_col]]   <- ifelse(n_pass > 0, 1, NA)
  
  summary_df
}

#' Summarize differential expression across multiple imputations (legacy-compatible)
#'
#' Reproduces the old script logic:
#'   - "stable" if passed in >= min_no_passed imputations
#'   - average linear ratio across imputations and convert to signed linearFC
#'   - summarize p-value and padj using quantile
#'
#' @param runs_de_tables List of length cfg$no_repetitions.
#' @param config Full config list.
#'
#' @return A data.frame with one row per feature.
summarize_limma_mult_imputation <- function(runs_de_tables, config) {
  de_cfg <- config$modes$proteomics$de
  imp_cfg   <- config$modes$proteomics$imputation
  
  NO_REPETITIONS <- as.integer(imp_cfg$no_repetitions)
  MIN_NO_PASSED  <- as.integer(imp_cfg$min_no_passed)
  
  use_adj_for_pass1 <- isTRUE(de_cfg$use_adj_for_pass1)
  p_cutoff <- as.numeric(de_cfg$p_cutoff)
  
  linear_fc_cutoff <- as.numeric(de_cfg$linear_fc_cutoff)
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
  
  # Validate alignment
  for (n in seq_len(NO_REPETITIONS)) {
    for (cn in contrasts) {
      cur <- runs_de_tables[[n]][[cn]]
      stopifnot(nrow(cur) == length(ref_ids))
      stopifnot(all(cur$FeatureID == ref_ids))
    }
  }
  
  # Loop over contrasts to calculate stability
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
  
  # Use the helper to add the global pass flags (DRY principle)
  out <- add_pass_any_contrast(out)
  
  return(out)
}