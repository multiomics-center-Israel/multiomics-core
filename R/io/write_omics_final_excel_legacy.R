#' Write legacy-style Final_results excels (generic for any mode)
#'
#' @param final_results data.frame (must include pass_any_contrast + manual_cutoffs.* optionally)
#' @param config config
#' @param out_dir output dir (run_dir/mode)
#' @param mode string, e.g. "proteomics" or "rna"
#' @param id_col column in final_results used as feature id (e.g., "Protein" / "Gene")
#' @param expr_for_de matrix used for zscore/order in DE-only excel (rows are feature IDs)
#' @param with_cutoffs logical; if TRUE add cutoffs sheet + fill manual formulas (if manual_cutoffs.* exists)
write_final_results_excels_legacy_generic <- function(final_results,
                                                      config,
                                                      out_dir,
                                                      mode,
                                                      id_col,
                                                      expr_for_de,
                                                      with_cutoffs = TRUE,
                                                      excel_order = NULL) {
  stopifnot(requireNamespace("openxlsx", quietly = TRUE))
  
  ptag <- p_tag(config)
  
  f_all <- file.path(out_dir, sprintf("Final_results_P_%s.xlsx", ptag))
  f_de  <- file.path(out_dir, sprintf("Final_results_DE_P_%s.xlsx", ptag))
  
  save_wb_results <- function(df, path, with_cutoffs = FALSE) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Results")
    openxlsx::writeData(wb, "Results", df)
    
    if (with_cutoffs) {
      add_cutoffs_sheet_legacy(wb, config, mode = mode)
      fill_manual_cutoffs_formulas_legacy(wb, "Results", df, config, mode = mode)
    }
    
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  }
  
  # 1) ALL
  save_wb_results(final_results, f_all, with_cutoffs = isTRUE(with_cutoffs))
  
  # 2) DE-only
  is_de <- !is.na(final_results$pass_any_contrast) & final_results$pass_any_contrast == 1
  de_df <- final_results[is_de, , drop = FALSE]
  de_df <- de_df[, !startsWith(names(de_df), "manual_cutoffs"), drop = FALSE]
  
  # optional: add order + zscore columns if expr_for_de is provided
  if (!is.null(expr_for_de)) {
    expr_for_de <- as.matrix(expr_for_de)
    stopifnot(!is.null(rownames(expr_for_de)))
    
    de_ids <- intersect(de_df[[id_col]], rownames(expr_for_de))
    mat_de <- expr_for_de[de_ids, , drop = FALSE]
    
    if (nrow(mat_de) == 0) {
      save_wb_results(mat_de, f_de, with_cutoffs = FALSE)
      return(c(f_all, f_de))
    }
    
    if (!is.null(excel_order)) {
      ordered_ids <- excel_order$ordered_ids
      z <- excel_order$zscore_mat
      
      ordered_ids <- intersect(ordered_ids, de_df[[id_col]])
      order_tbl <- data.frame(
        tmp_id = ordered_ids,
        order = seq_along(ordered_ids),
        stringsAsFactors = FALSE
      )
      names(order_tbl)[1] <- id_col
      
      z_tbl <- as.data.frame(z, check.names = FALSE)
      z_tbl[[id_col]] <- rownames(z_tbl)
      
      de_df <- dplyr::as_tibble(de_df) |>
        dplyr::left_join(order_tbl, by = id_col) |>
        dplyr::left_join(z_tbl, by = id_col) |>
        dplyr::distinct()
    } else {
      z <- zscore_rows(mat_de)
      colnames(z) <- paste0(colnames(z), ".zscore")
      
      ordered_ids <- get_hclust_row_order(z, row_distance = "correlation", hclust_method = "complete")
      
      order_tbl <- data.frame(
        feature_id = ordered_ids,
        order = seq_along(ordered_ids),
        stringsAsFactors = FALSE
      )
      names(order_tbl)[1] <- id_col
      
      z_tbl <- as.data.frame(z, check.names = FALSE)
      z_tbl[[id_col]] <- rownames(z_tbl)
      
      de_df <- dplyr::as_tibble(de_df) |>
        dplyr::left_join(order_tbl, by = id_col) |>
        dplyr::left_join(z_tbl, by = id_col) |>
        dplyr::distinct()
    }
    
    
  }
  
  
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

fill_manual_cutoffs_formulas_legacy <- function(wb, sheet, final_results, config,
                                                mode = "proteomics") {
  de_cfg <- config$modes[[mode]]$de
  use_fdr <- isTRUE(de_cfg$use_adj_for_pass1)
 
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
# ---- Proteomics final_file ----
write_final_results_excels_proteomics <- function(final_results, pre, config, out_dir, excel_order = NULL) {
  write_final_results_excels_legacy_generic(
    final_results = final_results,
    config       = config,
    out_dir      = out_dir,
    mode         = "proteomics",
    id_col       = "Protein",
    expr_for_de  = pre$expr_imp_single,  
    with_cutoffs = TRUE,
    excel_order  = excel_order            
  )
}


# ---- RNASeq final_file ----

write_final_results_excels_legacy_rna <- function(final_results, pre, config, out_dir) {
  write_final_results_excels_legacy_generic(
    final_results = final_results,
    config       = config,
    out_dir      = out_dir,
    mode         = "rna",
    id_col       = "Gene",
    expr_for_de  = pre$expr_work,
    with_cutoffs = TRUE 
  )
}
