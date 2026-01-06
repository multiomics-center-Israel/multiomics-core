#' Writes intermediate proteomics matrices
write_proteomics_datasets_legacy <- function(pre, runs = NULL, config, run_dir) {
  dirs <- create_legacy_output_dirs(run_dir)
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
write_limma_multimp_summary_legacy <- function(summary_df, config, run_dir) {
  dirs <- create_legacy_output_dirs(run_dir)
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
write_limma_results_multimp_legacy <- function(de_res, contrast_name, config, run_dir) {
  dirs <- create_legacy_output_dirs(run_dir)
  
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
#' @param run_dir run directory
#' @param write_runs logical; if TRUE writes per-imputation matrices (if available in de_res$runs)
#'
#' @return character vector of written file paths
write_proteomics_multimpute_outputs <- function(pre, de_res,
                                                inputs, config, run_dir, write_runs = FALSE) {
  
  files <- character(0)
  dirs  <- create_legacy_output_dirs(run_dir)
  
  # 1) datasets
  runs_for_datasets <- if (isTRUE(write_runs)) de_res$runs else NULL
  files <- c(files, write_proteomics_datasets_legacy(pre, runs_for_datasets, config, run_dir))
  
  # 2) summary
  if (!is.null(de_res$summary_df)) {
    files <- c(files, write_limma_multimp_summary_legacy(de_res$summary_df, config, run_dir))
  }
  
  # 3) wide limma per contrast (optional)
  if (!is.null(de_res$runs_de_tables) && length(de_res$runs_de_tables) > 0) {
    contrast_names <- names(de_res$runs_de_tables[[1]])
    for (cn in contrast_names) {
      files <- c(files, write_limma_results_multimp_legacy(de_res = de_res, contrast_name = cn, config = config, run_dir = run_dir))
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
    files <- c(files, write_final_results_excels_legacy(final_results, pre, config, run_dir))
  }
  
  unique(files)
}
