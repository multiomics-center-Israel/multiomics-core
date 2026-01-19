#' Writes intermediate proteomics matrices
write_proteomics_datasets_legacy <- function(pre, runs = NULL, config, out_dir) {
  dirs <- create_legacy_output_dirs(out_dir)
  cfg  <- config$modes$proteomics
  files <- character(0)
  
  files <- c(files, save_tsv(as.data.frame(pre$expr_filt, check.names = FALSE),
                              dirs$datasets, "protein_log2_filtered_unimputed.tsv"))
  
  fname_imp <- sprintf(
    "protein_log2_filtered_imputed_once_width_%s_shift_%s.tsv",
    cfg$imputation$width, cfg$imputation$downshift
  )
  
  files <- c(files, save_tsv(as.data.frame(pre$expr_imp_single,  check.names = FALSE),
                              dirs$datasets, fname_imp))
 
  
  if (!is.null(runs) && length(runs) > 0) {
    rep_dir <- file.path(dirs$datasets, "imputed_repetitions")
    ensure_dir(rep_dir)
    
    for (i in seq_along(runs)) {
      expr_i <- runs[[i]]$expr_imp
      if (is.null(expr_i)) next
      
      # repetitions
      files <- c(files, save_tsv(as.data.frame(expr_i, check.names = FALSE),
                                  rep_dir, sprintf("protein_log2_filtered_imputed_%02d.tsv", i)))
      
    }
  }
  
  unique(files)
}

#' Write legacy-style  multi-imputation summary
write_limma_multimp_summary_legacy <- function(summary_df, config, out_dir) {
  dirs <- create_legacy_output_dirs(out_dir)
  save_tsv(summary_df, dirs$datasets, sprintf("limma_multimp_summary_p%s.tsv", p_tag(config)))
}

#' Build legacy-style wide limma table across imputations
build_limma_results_multimp_wide <- function(runs_de_tables, contrast_name, 
                                             stats_cols = c("logFC", "P.Value", "adj.P.Val")) {
  stopifnot(length(runs_de_tables) >= 1)
  
  # Anchor with first run
  base <- runs_de_tables[[1]][[contrast_name]]
  
  # FIX: Fail-fast if first run is broken
  if (is.null(base)) stop(sprintf("Contrast '%s' missing in imputation run 1", contrast_name))
  stopifnot("FeatureID" %in% colnames(base))
  
  id_cols <- intersect(c("FeatureID", "Protein.Names", "Genes", "First.Protein.Description", "Contrast"), colnames(base))
  out <- base[, id_cols, drop = FALSE]
  
  # Append stats per imputation
  for (i in seq_along(runs_de_tables)) {
    tab <- runs_de_tables[[i]][[contrast_name]]
    
    # FIX: Strict check - missing imputation data is critical
    if (is.null(tab)) {
      stop(sprintf("Critical: Contrast '%s' is missing in imputation run %d (but existed in run 1)", contrast_name, i))
    }
    
    # Ensure alignment
    tab <- align_de_table_by_feature_id(
      tab = tab,
      ref_ids = out$FeatureID,
      run_i = i,
      contrast_name = contrast_name,
      id_col = "FeatureID"
    )
    
    
    stat_block <- tab[, intersect(stats_cols, colnames(tab)), drop = FALSE]
    colnames(stat_block) <- paste0(colnames(stat_block), ".", i)
    out <- cbind(out, stat_block)
  }
  
  out
}

#' Write legacy-style limma results across multiple imputations
write_limma_results_multimp_legacy <- function(de_res, contrast_name, config, out_dir) {
  dirs <- create_legacy_output_dirs(out_dir)
  
  wide_df <- build_limma_results_multimp_wide(
    runs_de_tables = de_res$runs_de_tables,
    contrast_name  = contrast_name
  )
  
  fname <- sprintf("limma_results_multimp_p%s.tsv", p_tag(config))
  save_tsv(wide_df, dirs$datasets, fname)
}


#' Build legacy-style Final_results table
build_final_results_proteomics <- function(pre, summary_df, contrasts_df, row_data = NULL) {
  
  expr_df <- as.data.frame(pre$expr_filt, check.names = FALSE)
  if (is.null(row_data)) row_data <- pre$row_data
  
  stopifnot(!is.null(row_data), "Protein.Group" %in% colnames(row_data))
  
  # Initialize Base
  base <- data.frame(
    Protein = row_data[["Protein.Group"]],
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  # Match to Summary
  m <- match(base$Protein, summary_df$FeatureID)
  
  # Add Annotations
  for (col in c("Protein.Names", "Genes", "First.Protein.Description")) {
    val <- if (col %in% colnames(summary_df)) summary_df[[col]][m] else NA
    if (col %in% colnames(row_data)) {
      fallback <- row_data[[col]]
      val <- ifelse(is.na(val), fallback, val)
    }
    base[[col]] <- val
  }
  
  # Add Expression
  base <- cbind(base, expr_df)
  
  # Add Contrast Stats
  contrast_names <- contrasts_df$Contrast_name
  
  # FIX: Pre-validate that all required columns exist in summary_df
  # This prevents partial failure inside the loop
  for (cn in contrast_names) {
    cols <- get_contrast_cols(cn)
    needed <- c(cols$fc, cols$p, cols$padj) # pass is optional-ish but usually needed
    missing_cols <- setdiff(needed, colnames(summary_df))
    if (length(missing_cols) > 0) {
      stop(sprintf("Summary DF missing columns for contrast '%s': %s", cn, paste(missing_cols, collapse=", ")))
    }
  }
  
  for (cn in contrast_names) {
    cols <- get_contrast_cols(cn) 
    
    fc_vals   <- summary_df[[cols$fc]][m]
    # FIX: Robust check for 'pass' column existence
    pass_vals <- if (cols$pass %in% colnames(summary_df)) summary_df[[cols$pass]][m] else rep(NA, length(m))
    
    base[[cols$fc]]   <- fc_vals
    base[[cols$p]]    <- summary_df[[cols$p]][m]
    base[[cols$padj]] <- summary_df[[cols$padj]][m]
    
    # Calculate Up/Down logic
    base[[cols$updown]] <- ifelse(!is.na(pass_vals),
                                  ifelse(as.numeric(fc_vals) >= 0, "up", "down"), "")
    base[[cols$manual]] <- NA 
  }
  
  # FIX: Robust pass_any_contrast logic
  pass_cols <- paste0("pass.imputs.", contrast_names)
  existing_pass_cols <- intersect(pass_cols, colnames(summary_df))
  
  if (length(existing_pass_cols) == 0) {
    warning("No 'pass.imputs' columns found in summary_df. pass_any_contrast will be NA.")
    base$pass_any_contrast <- NA
  } else {
    # Create matrix, handle potential NAs in 'm' automatically (returns NA row)
    pass_mat <- summary_df[m, existing_pass_cols, drop = FALSE]
    # rowSums ignoring NAs, check if > 0
    # Note: !is.na(NA) is FALSE, so NA entries don't contribute to the sum
    base$pass_any_contrast <- ifelse(rowSums(!is.na(pass_mat)) > 0, 1, NA)
  }
  
  base
}

#' Orchestrator: write all proteomics multi-imputation outputs (legacy-compatible)
#'
#' @param pre Preprocessed object
#' @param de_res List with $runs_de_tables, $runs, $summary_df
#' @param inputs proteomics inputs (expects $contrasts at least)
#' @param config global config
#' @param out_dir run directory
#' @param write_runs logical; if TRUE writes per-imputation matrices (if available in de_res$runs)
#'
#' @return character vector of written file paths
write_proteomics_multimpute_outputs <- function(pre, de_res,
                                                inputs, config, out_dir,excel_order = NULL,
                                                write_runs = FALSE) {
  
  files <- character(0)
  dirs  <- create_legacy_output_dirs(out_dir)
  
  # 1) datasets
  runs_for_datasets <- if (isTRUE(write_runs)) de_res$runs else NULL
  files <- c(files, write_proteomics_datasets_legacy(pre, runs_for_datasets, config, out_dir))
  
  # 2) summary
  if (!is.null(de_res$summary_df)) {
    files <- c(files, write_limma_multimp_summary_legacy(de_res$summary_df, config, out_dir))
  }
  
  # 3) wide limma per contrast (optional)
  if (!is.null(de_res$runs_de_tables) && length(de_res$runs_de_tables) > 0) {
    contrast_names <- names(de_res$runs_de_tables[[1]])
    for (cn in contrast_names) {
      files <- c(files, write_limma_results_multimp_legacy(de_res = de_res, contrast_name = cn, config = config, out_dir = out_dir))
    }
  }
  
  # 4) final results TSV (optional but useful)
  if (!is.null(inputs$contrasts) && !is.null(de_res$summary_df)) {
    final_results <- build_final_results_proteomics(
      pre          = pre,
      summary_df   = de_res$summary_df,
      contrasts_df = inputs$contrasts,
      row_data     = pre$row_data
    )
    files <- c(files, save_tsv(final_results, dirs$datasets, "final_results.tsv"))
    
    # 5) Excel outputs (delegated to dedicated file)
    # files <- c(files, write_final_results_excels_proteomics(final_results, pre, config, out_dir))
    files <- c(files, write_final_results_excels_proteomics(
      final_results = final_results,
      pre          = pre,
      config       = config,
      out_dir      = out_dir,
      excel_order  = excel_order   # NEW
    ))
    
  }
  
  unique(files)
}

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
