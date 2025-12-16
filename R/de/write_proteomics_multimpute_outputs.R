#' Write legacy-style proteomics datasets (with improved copies)
#'
#' Writes intermediate proteomics matrices in both:
#' 1) Legacy filenames (*.txt, TSV) for backward compatibility
#' 2) Improved filenames (*.tsv) with clearer naming
#'
#' @param pre Output of preprocess_proteomics()
#' @param runs Optional list of limma/imputation runs
#' @param config Full config list
#' @param run_dir Run output root directory
#' @return Character vector of written file paths
write_proteomics_datasets_legacy <- function(pre, runs = NULL, config, run_dir) {
  
  cfg  <- config$modes$proteomics
  dirs <- create_legacy_output_dirs(run_dir)
  
  files_written <- c()
  
  ## 1) Unimputed (log2) expression matrix
  files_written <- c(
    files_written,
    write_dual_tsv(
      x = pre$expr_filt,
      legacy_path   = file.path(dirs$datasets, "LFQ_data_unimputed.txt"),
      improved_path = file.path(dirs$datasets, "protein_log2_filtered_unimputed.tsv")
      
    )
  )
  
  ## 2) Single imputed matrix
  width     <- cfg$imputation$width
  downshift <- cfg$imputation$downshift
  
  files_written <- c(
    files_written,
    write_dual_tsv(
      x = pre$expr_imp,
      legacy_path = file.path(
        dirs$datasets,
        sprintf("Filtered_imputed_data_once_WIDTH_%s_SHIFT_%s.txt", width, downshift)
      ),
      improved_path = file.path(
        dirs$datasets,
        sprintf(
          "protein_log2_filtered_imputed_once_width_%s_shift_%s.tsv",
          width, downshift
        )
      )
      
    )
  )
  
  ## 3) Multiple imputations (if provided)
  if (!is.null(runs)) {
    rep_dir_legacy   <- file.path(dirs$datasets, "Imputed_repetitions")
    rep_dir_improved <- file.path(dirs$datasets, "imputed_repetitions")
    
    dir.create(rep_dir_legacy,   showWarnings = FALSE, recursive = TRUE)
    dir.create(rep_dir_improved, showWarnings = FALSE, recursive = TRUE)
    
    for (i in seq_along(runs)) {
      expr_i <- runs[[i]]$expr_imp %||% runs[[i]]$expr_imputed %||% NULL
      if (is.null(expr_i)) next
      
      files_written <- c(
        files_written,
        write_dual_tsv(
          x = expr_i,
          legacy_path   = file.path(rep_dir_legacy,   sprintf("Imputed.%d.txt", i)),
          improved_path = file.path(
            rep_dir_improved,
            sprintf("protein_log2_filtered_imputed_%02d.tsv", i)
          )
          
        )
      )
    }
  }
  
  files_written
}
#' Write legacy-style Limma multi-imputation summary (with improved copy)
#'
#' Creates the old-pipeline summary table file under Datasets/:
#'   - Limma_summary_mult_imputs_P_<p>.txt  (legacy)
#' and also a clearer TSV filename (improved).
#'
#' @param summary_df data.frame produced by summarize_mult_imputation()
#' @param config Full config list
#' @param run_dir Run output root directory
#' @return Character vector of written file paths
write_limma_multimp_summary_legacy <- function(summary_df, config, run_dir) {
  dirs <- create_legacy_output_dirs(run_dir)
  
  p_cut <- config$modes$proteomics$limma$p_cutoff
  p_tag <- format(p_cut, trim = TRUE, scientific = FALSE)
  
  legacy_path <- file.path(dirs$datasets, sprintf("Limma_summary_mult_imputs_P_%s.txt", p_tag))
  improved_path <- file.path(dirs$datasets, sprintf("limma_multimp_summary_p%s.tsv", p_tag))
  
  write_dual_tsv(summary_df, legacy_path, improved_path)
}
#' Build legacy-style wide limma table across imputations
#'
#' The old pipeline concatenated per-imputation stats columns (with ".n" suffix).
#' This function reproduces that structure from runs_de_tables.
#'
#' @param runs_de_tables List (length N imputations) of per-contrast DE tables.
#'        Each element is a named list of contrasts -> data.frame.
#' @param contrast_name Which contrast to export (e.g., "AHA_vs_OXY").
#' @param stats_cols Which columns to export per imputation.
#' @return data.frame wide legacy-style results (row-aligned by FeatureID)
build_limma_results_multimp_wide <- function(
    runs_de_tables,
    contrast_name,
    stats_cols = c("logFC", "P.Value", "adj.P.Val")
) {
  stopifnot(length(runs_de_tables) >= 1)
  
  # Take base table from first run to anchor FeatureID + annotations
  base <- runs_de_tables[[1]][[contrast_name]]
  stopifnot(!is.null(base))
  stopifnot("FeatureID" %in% colnames(base))
  
  # Keep an "ID block" (what you want on the left side)
  id_cols <- intersect(
    c("FeatureID", "Protein.Names", "Genes", "First.Protein.Description", "Contrast"),
    colnames(base)
  )
  out <- base[, id_cols, drop = FALSE]
  
  # Append stats per imputation with suffix .n
  for (i in seq_along(runs_de_tables)) {
    tab <- runs_de_tables[[i]][[contrast_name]]
    stopifnot(!is.null(tab))
    
    # Ensure same row order via FeatureID
    tab <- tab[match(out$FeatureID, tab$FeatureID), , drop = FALSE]
    
    stat_block <- tab[, intersect(stats_cols, colnames(tab)), drop = FALSE]
    colnames(stat_block) <- paste0(colnames(stat_block), ".", i)
    
    out <- cbind(out, stat_block)
  }
  
  out
}
#' Write legacy-style limma results across multiple imputations (with improved copy)
#'
#' Produces:
#' - Datasets/Limma_results_mult_imput_P_<p>.txt (legacy)
#' - Datasets/limma_results_multimp_p<p>.tsv (improved)
#'
#' @param de_res Output of run_proteomics_de_mult_impute()
#' @param contrast_name Contrast to export
#' @param config Full config list
#' @param run_dir Run output root directory
#' @return Character vector of written file paths
write_limma_results_multimp_legacy <- function(de_res, contrast_name, config, run_dir) {
  dirs <- create_legacy_output_dirs(run_dir)
  
  p_cut <- config$modes$proteomics$limma$p_cutoff
  p_tag <- format(p_cut, trim = TRUE, scientific = FALSE)
  
  wide_df <- build_limma_results_multimp_wide(
    runs_de_tables = de_res$runs_de_tables,
    contrast_name  = contrast_name,
    stats_cols     = c("logFC", "P.Value", "adj.P.Val")
  )
  
  legacy_path   <- file.path(dirs$datasets, sprintf("Limma_results_mult_imput_P_%s.txt", p_tag))
  improved_path <- file.path(dirs$datasets, sprintf("limma_results_multimp_p%s.tsv", p_tag))
  
  write_dual_tsv(wide_df, legacy_path, improved_path)
}
#' Build legacy-style Final_results table (all proteins)
#'
#' Mirrors the old Neat proteomics "final_results" table:
#' - Protein ID + annotation
#' - expression matrix (log2, filtered, unimputed)
#' - per-contrast summary stats (linearFC/pvalue/padj) + up/down label
#' - pass_any_contrast flag
#' - manual_cutoffs.<contrast> placeholder (legacy column)
#'
#' @param pre Output of preprocess_proteomics() (needs expr_filt and meta)
#' @param summary_df Multi-imputation summary table (prot_de$summary_df)
#' @param contrasts_df Contrasts table (inputs$contrasts)
#' @return data.frame final_results (legacy-compatible columns)
build_final_results_proteomics <- function(pre, summary_df, contrasts_df, row_data = NULL) {
  
  expr_df <- as.data.frame(pre$expr_filt, check.names = FALSE)
  
  # Use row_data to map numeric rownames -> Protein.Group IDs
  if (is.null(row_data)) row_data <- pre$row_data
  stopifnot(!is.null(row_data))
  stopifnot(nrow(row_data) == nrow(expr_df))
  stopifnot("Protein.Group" %in% colnames(row_data))
  
  protein_ids <- row_data[["Protein.Group"]]
  
  base <- data.frame(
    Protein = protein_ids,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Add annotation (prefer summary_df columns if present; fallback to row_data)
  m <- match(base$Protein, summary_df$FeatureID)
  
  # Preferred annot from summary_df (already aligned to FeatureID)
  base$Protein.Names <- summary_df$Protein.Names[m]
  base$Genes <- summary_df$Genes[m]
  base$First.Protein.Description <- summary_df$First.Protein.Description[m]
  
  # If any annotation is missing, fill from row_data when available
  for (col in c("Protein.Names", "Genes", "First.Protein.Description")) {
    if (any(is.na(base[[col]])) && (col %in% colnames(row_data))) {
      base[[col]][is.na(base[[col]])] <- row_data[[col]][is.na(base[[col]])]
    }
  }
  
  # Add expression columns (log2, filtered, unimputed)
  base <- cbind(base, expr_df)
  
  # Add per-contrast stats + up/down + manual_cutoffs placeholder
  contrast_names <- contrasts_df$Contrast_name
  
  for (cn in contrast_names) {
    fc_col   <- paste0("linearFC.imputs.", cn)
    p_col    <- paste0("pvalue.imputs.",  cn)
    padj_col <- paste0("padj.imputs.",    cn)
    pass_col <- paste0("pass.imputs.",    cn)
    
    updown_col <- paste0("upDown.imputs.", cn)
    manual_col <- paste0("manual_cutoffs.", cn)
    
    fc_vals   <- summary_df[[fc_col]][m]
    pass_vals <- summary_df[[pass_col]][m]
    
    base[[fc_col]]     <- fc_vals
    base[[p_col]]      <- summary_df[[p_col]][m]
    base[[padj_col]]   <- summary_df[[padj_col]][m]
    base[[updown_col]] <- ifelse(!is.na(pass_vals),
                                 ifelse(as.numeric(fc_vals) >= 0, "up", "down"),
                                 "")
    base[[manual_col]] <- NA
  }
  
  # pass_any_contrast
  pass_mat <- sapply(contrast_names, function(cn) summary_df[[paste0("pass.imputs.", cn)]][m])
  base$pass_any_contrast <- ifelse(rowSums(!is.na(pass_mat)) > 0, 1, NA)
  
  base
}

#' Write legacy Final_results.txt (with improved copy)
#'
#' Exports the full "final_results" table (annotation + expression + per-contrast stats)
#' to the run root, matching the old pipeline filename.
#'
#' @param final_results data.frame from build_final_results_proteomics()
#' @param run_dir Run output root directory
#' @return Character vector of written file paths
write_final_results_table_legacy <- function(final_results, run_dir) {
  legacy_path   <- file.path(run_dir, "Final_results.txt")
  improved_path <- file.path(run_dir, "final_results.tsv")
  write_dual_tsv(final_results, legacy_path, improved_path)
}
#' Write legacy Excel outputs for proteomics final results
#'
#' Produces two Excel files at the run root:
#'  - Final_results_P_<p>.xlsx    (all proteins)
#'  - Final_results_DE_P_<p>.xlsx (DE proteins only; pass_any_contrast == 1)
#'
#' @param final_results data.frame from build_final_results_proteomics()
#' @param config Full config list
#' @param run_dir Run output root directory
#' @return Character vector of written file paths
write_final_results_excels_legacy <- function(final_results, config, run_dir) {
  stopifnot(requireNamespace("openxlsx", quietly = TRUE))
  
  p_cut <- config$modes$proteomics$limma$p_cutoff
  p_tag <- format(p_cut, trim = TRUE, scientific = FALSE)
  
  f_all <- file.path(run_dir, sprintf("Final_results_P_%s.xlsx", p_tag))
  f_de  <- file.path(run_dir, sprintf("Final_results_DE_P_%s.xlsx", p_tag))
  
  # ALL
  wb1 <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb1, "Results")
  openxlsx::writeData(wb1, "Results", final_results)
  add_cutoffs_sheet_legacy(wb1, config)
  fill_manual_cutoffs_formulas_legacy(wb1, "Results", final_results, config)
  openxlsx::saveWorkbook(wb1, f_all, overwrite = TRUE)
  
  # DE only
  de_df <- final_results[!is.na(final_results$pass_any_contrast) & final_results$pass_any_contrast == 1, , drop = FALSE]
  wb2 <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb2, "Results")
  openxlsx::writeData(wb2, "Results", de_df)
  add_cutoffs_sheet_legacy(wb2, config)
  fill_manual_cutoffs_formulas_legacy(wb2, "Results", de_df, config)
  openxlsx::saveWorkbook(wb2, f_de, overwrite = TRUE)
  
  c(f_all, f_de)
}

#' Add legacy-style "Cutoffs" sheet to an openxlsx workbook
#'
#' Reproduces the old Neat proteomics Cutoffs tab (layout + styles),
#' but reads values from config$modes$proteomics$limma.
#'
#' NOTE:
#' In the legacy script, FDR_ADJ controlled which cutoff cell
#' (p-value vs adjusted p-value) was considered "active".
#' Here we map that behavior to limma$use_adj_for_pass1
#' for backward compatibility.
#'
#' @param wb openxlsx workbook object
#' @param config Full config list
add_cutoffs_sheet_legacy <- function(wb, config) {
  stopifnot(requireNamespace("openxlsx", quietly = TRUE))
  
  limma_cfg <- config$modes$proteomics$limma
  
  # Legacy mapping
  FDR_ADJ <- isTRUE(limma_cfg$use_adj_for_pass1)
  
  P_CUTOFF <- as.numeric(limma_cfg$p_cutoff)
  LINEAR_FC_CUTOFF <- as.numeric(limma_cfg$linear_fc_cutoff)
  
  if (is.na(P_CUTOFF))
    stop("Cutoffs sheet: limma$p_cutoff is NA or not numeric.")
  if (is.na(LINEAR_FC_CUTOFF))
    stop("Cutoffs sheet: limma$linear_fc_cutoff is NA or not numeric.")
  
  openxlsx::addWorksheet(wb, sheetName = "Cutoffs", gridLines = TRUE)
  
  openxlsx::writeData(
    wb, "Cutoffs",
    x = c("p-value", "Adjusted pvalue (FDR)", "linear Fold Change (linearFC)"),
    startCol = 2,
    startRow = 4
  )
  
  openxlsx::writeData(
    wb, "Cutoffs",
    x = c(
      ifelse(FDR_ADJ, "", P_CUTOFF),
      ifelse(FDR_ADJ, P_CUTOFF, ""),
      LINEAR_FC_CUTOFF
    ),
    startCol = 3,
    startRow = 4
  )
  
  openxlsx::createNamedRegion(wb, "Cutoffs", cols = 3, rows = 4, name = "PVAL_CO")
  openxlsx::createNamedRegion(wb, "Cutoffs", cols = 3, rows = 5, name = "FDR_CO")
  openxlsx::createNamedRegion(wb, "Cutoffs", cols = 3, rows = 6, name = "LFC_CO")
  
  style_used <- openxlsx::createStyle(
    border = "TopBottomLeftRight",
    borderStyle = "thick",
    fgFill = "green",
    halign = "center"
  )
  style_not <- openxlsx::createStyle(
    border = "TopBottomLeftRight",
    borderStyle = "thick",
    fgFill = "red",
    halign = "center"
  )
  
  used_rows    <- if (FDR_ADJ) 5:6 else c(4, 6)
  notused_rows <- if (FDR_ADJ) 4    else 5
  
  openxlsx::addStyle(
    wb, "Cutoffs", style_not,
    rows = notused_rows, cols = 3,
    gridExpand = FALSE, stack = FALSE
  )
  openxlsx::addStyle(
    wb, "Cutoffs", style_used,
    rows = used_rows, cols = 3,
    gridExpand = FALSE, stack = FALSE
  )
  
  openxlsx::setColWidths(wb, "Cutoffs", cols = 2, widths = "auto")
}
#' Write legacy-style proteomics multi-imputation outputs (orchestrator)
#'
#' This is a convenience wrapper used by {targets}. It writes the legacy-compatible
#' output set (datasets + summary + wide limma + final results + excel files)
#' and returns a character vector of file paths.
#'
#' @param pre Preprocess output from preprocess_proteomics()
#' @param de_res Output of run_proteomics_de_mult_impute()
#' @param inputs Proteomics inputs list (for contrasts)
#' @param config Full config list
#' @param run_dir Run output root directory (e.g., outputs/Results_<project>_<round>)
#' @param write_runs Logical; if TRUE, also write per-imputation run outputs (optional; not implemented here)
#' @return Character vector of file paths written
write_proteomics_multimpute_outputs <- function(pre, de_res, inputs, config, run_dir, write_runs = FALSE) {
  
  # Ensure legacy directory structure exists
  create_legacy_output_dirs(run_dir)
  
  # 1) Datasets (unimputed + imputed once)
  files_ds <- write_proteomics_datasets_legacy(
    pre     = pre,
    runs    = NULL,
    config  = config,
    run_dir = run_dir
  )
  
  # 2) Limma multi-imputation stability summary
  files_sum <- write_limma_multimp_summary_legacy(
    summary_df = de_res$summary_df,
    config     = config,
    run_dir    = run_dir
  )
  
  # 3) Wide legacy limma results across imputations (per contrast)
  contrast_names <- as.data.frame(inputs$contrasts)$Contrast_name
  files_lr <- unlist(lapply(contrast_names, function(cn) {
    write_limma_results_multimp_legacy(
      de_res        = de_res,
      contrast_name = cn,
      config        = config,
      run_dir       = run_dir
    )
  }), use.names = FALSE)
  
  # 4) Final results table (all proteins)
  final_results <- build_final_results_proteomics(
    pre          = pre,
    summary_df   = de_res$summary_df,
    contrasts_df = as.data.frame(inputs$contrasts),
    row_data     = pre$row_data
  )
  
  files_final <- write_final_results_table_legacy(
    final_results = final_results,
    run_dir       = run_dir
  )
  
  # 5) Excel outputs (all + DE) + Cutoffs tab
  files_xlsx <- write_final_results_excels_legacy(
    final_results = final_results,
    config        = config,
    run_dir       = run_dir
  )
  
  # (Optional) write_runs reserved for future (kept for API compatibility)
  invisible(write_runs)
  
  c(files_ds, files_sum, files_lr, files_final, files_xlsx)
}

#' Fill manual_cutoffs.* columns in an Excel sheet with legacy-style formulas
#'
#' Writes Excel formulas that use named ranges (PVAL_CO, FDR_CO, LFC_CO) from the Cutoffs sheet.
#' The formula is written per-row and per-contrast into the manual_cutoffs column.
#'
#' @param wb openxlsx workbook
#' @param sheet Sheet name where results table was written (e.g. "Results")
#' @param final_results data.frame that was written to the sheet
#' @param config Full config list
fill_manual_cutoffs_formulas_legacy <- function(wb, sheet, final_results, config) {
  stopifnot(requireNamespace("openxlsx", quietly = TRUE))
  
  limma_cfg <- config$modes$proteomics$limma
  use_fdr   <- isTRUE(limma_cfg$use_adj_for_pass1)  # legacy mapping
  
  # Find all contrasts by manual_cutoffs.imputs.<contrast>
  manual_cols <- grep("^manual_cutoffs\\.", names(final_results), value = TRUE)
  if (length(manual_cols) == 0) return(invisible(NULL))
  
  # Header row is 1, data starts at row 2 (because writeData writes colnames)
  start_row <- 2
  n <- nrow(final_results)
  
  # Helper to turn numeric column index into Excel column letters
  num_to_excel_col <- function(num) {
    letters <- c(LETTERS)
    out <- ""
    while (num > 0) {
      r <- (num - 1) %% 26
      out <- paste0(letters[r + 1], out)
      num <- (num - 1) %/% 26
    }
    out
  }
  
  # For each contrast, write a formula in manual_cutoffs column
  for (mcol in manual_cols) {
    contrast <- sub("^manual_cutoffs\\.", "", mcol)
    
    fc_col   <- paste0("linearFC.imputs.", contrast)
    p_col    <- paste0("pvalue.imputs.",  contrast)
    padj_col <- paste0("padj.imputs.",    contrast)
    
    # Sanity: required cols exist
    needed <- c(fc_col, p_col, padj_col)
    if (!all(needed %in% names(final_results))) next
    
    # Column indices within the written sheet (1-based)
    fc_i   <- match(fc_col,   names(final_results))
    p_i    <- match(p_col,    names(final_results))
    padj_i <- match(padj_col, names(final_results))
    m_i    <- match(mcol,     names(final_results))
    
    fc_L   <- num_to_excel_col(fc_i)
    p_L    <- num_to_excel_col(p_i)
    padj_L <- num_to_excel_col(padj_i)
    m_L    <- num_to_excel_col(m_i)
    
    rows <- start_row:(start_row + n - 1)
    
    # Build formulas per row
    formulas <- vapply(rows, function(r) {
      fc_ref   <- paste0(fc_L, r)
      p_ref    <- paste0(p_L, r)
      padj_ref <- paste0(padj_L, r)
      
      if (use_fdr) {
        # Use adjusted p-value cutoff (FDR_CO) + LFC_CO
        paste0(
          'IF(AND(ISNUMBER(', padj_ref, '),',
          padj_ref, '<=FDR_CO,',
          'ABS(', fc_ref, ')>=LFC_CO),',
          'IF(', fc_ref, '>0,"up","down"),"")'
        )
      } else {
        # Use raw p-value cutoff (PVAL_CO) + LFC_CO
        paste0(
          'IF(AND(ISNUMBER(', p_ref, '),',
          p_ref, '<=PVAL_CO,',
          'ABS(', fc_ref, ')>=LFC_CO),',
          'IF(', fc_ref, '>0,"up","down"),"")'
        )
      }
    }, character(1))
    
    # Write formulas into the manual_cutoffs column
    openxlsx::writeFormula(
      wb, sheet = sheet,
      x = formulas,
      startCol = m_i,
      startRow = start_row
    )
  }
  
  invisible(NULL)
}

