
write_final_results_excels_legacy <- function(final_results, pre, config, run_dir) {
    stopifnot(requireNamespace("openxlsx", quietly = TRUE))
    cfg <- config$modes$proteomics
    p_tag <- p_tag(config)
    f_all <- file.path(run_dir, sprintf("Final_results_P_%s.xlsx", p_tag))
    f_de <- file.path(run_dir, sprintf("Final_results_DE_P_%s.xlsx", p_tag))
    
    # Helper to write workbook to avoid duplication
    save_wb_results <- function(df, path, with_cutoffs = FALSE) {
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Results")
      openxlsx::writeData(wb, "Results", df)
      
      if (with_cutoffs) {
        add_cutoffs_sheet_legacy(wb, config)
        fill_manual_cutoffs_formulas_legacy(wb, "Results", df, config)
      }
      openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
    }
    
    # 1. ALL Results (With formulas and cutoffs sheet)
    save_wb_results(final_results, f_all, with_cutoffs = TRUE)
    
    # 2. DE Only (Clean table)
    # Filter for pass_any_contrast == 1 and remove formula columns
    is_de <- !is.na(final_results$pass_any_contrast) & final_results$pass_any_contrast == 1
    de_df <- final_results[is_de, , drop = FALSE]
    
    # Remove manual_cutoffs columns for the clean DE file
    de_df <- de_df[, !startsWith(names(de_df), "manual_cutoffs"), drop = FALSE]
    
    id_col <- "Protein"
    expr_for_de <- as.matrix(pre$expr_imp)
    stopifnot(!is.null(rownames(expr_for_de)))
    de_ids <- intersect(de_df[[id_col]], rownames(expr_for_de))
    mat_de <- expr_for_de[de_ids, , drop = FALSE]
    
    z <- zscore_rows(mat_de)
    colnames(z) <- paste0(colnames(z), ".zscore")
    
    ordered_ids <- get_hclust_row_order(z, row_distance = "correlation", hclust_method = "complete")
    
    order_tbl <- data.frame(
      Protein = ordered_ids,
      order = seq_along(ordered_ids),
      stringsAsFactors = FALSE
    )
    
    z_tbl <- as.data.frame(z, check.names = FALSE)
    z_tbl$Protein <- rownames(z_tbl)
    
    de_df <- dplyr::as_tibble(de_df) |>
      dplyr::left_join(order_tbl, by = "Protein") |>
      dplyr::left_join(z_tbl, by = "Protein") |>
      dplyr::distinct()
    
    
    
    save_wb_results(de_df, f_de, with_cutoffs = FALSE)
    
    c(f_all, f_de)
  }
#' Fill manual_cutoffs.* columns with Excel formulas (legacy)
#'
#' Writes one formula for the entire manual_cutoffs.* column using relative
#' row references and named cutoff cells: PVAL_CO / FDR_CO / LFC_CO.
#'
#' @param wb openxlsx workbook
#' @param sheet sheet name (usually "Results")
#' @param final_results data.frame written to the sheet
#' @param config full config (uses modes$proteomics$de$use_adj_for_pass1)
fill_manual_cutoffs_formulas_legacy <- function(wb, sheet, final_results, config) {
  de_cfg <- config$modes$proteomics$de
  use_fdr   <- isTRUE(de_cfg$use_adj_for_pass1)
  
  manual_cols <- grep("^manual_cutoffs\\.", names(final_results), value = TRUE)
  if (length(manual_cols) == 0) return(invisible(NULL))
  
  start_row <- 2
  n_rows <- nrow(final_results)
  
  if (n_rows < 1) return(invisible(NULL))
  
  # Generate a sequence of row numbers (e.g., 2, 3, 4...)
  rows_seq <- start_row:(start_row + n_rows - 1)
  
  for (mcol in manual_cols) {
    contrast <- sub("^manual_cutoffs\\.", "", mcol)
    cols <- get_contrast_cols(contrast)
    
    # Require the stats columns for this contrast (skip quietly if missing)
    if (!all(c(cols$fc, cols$p, cols$padj) %in% names(final_results))) next
    
    # Excel column letters
    fc_L   <- openxlsx::int2col(match(cols$fc,   names(final_results)))
    p_L    <- openxlsx::int2col(match(cols$p,    names(final_results)))
    padj_L <- openxlsx::int2col(match(cols$padj, names(final_results)))
    m_i    <- match(mcol, names(final_results))
    
    fc_refs   <- paste0(fc_L,   rows_seq)
    p_refs    <- paste0(p_L,    rows_seq)
    padj_refs <- paste0(padj_L, rows_seq)
    
    #  Use vectorized sprintf to create a unique formula for each row
    cond <- if (use_fdr) {
      sprintf('AND(ISNUMBER(%s),%s<=FDR_CO,ABS(%s)>=LFC_CO)', 
              padj_refs, padj_refs, fc_refs)
    } else {
      sprintf('AND(ISNUMBER(%s),%s<=PVAL_CO,ABS(%s)>=LFC_CO)', 
              p_refs, p_refs, fc_refs)
    }
    
    formulas <- sprintf('IF(%s,IF(%s>0,"up","down"),"")', cond, fc_refs)
    
    # Write the vector of formulas
    openxlsx::writeFormula(
      wb, sheet,
      x        = formulas,  # <--- Now passing a vector of different formulas
      startCol = m_i,
      startRow = start_row
    )
  }
  
  invisible(TRUE)
}

