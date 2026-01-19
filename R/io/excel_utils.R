
#' Add legacy Cutoffs sheet + named regions (PVAL_CO/FDR_CO/LFC_CO)
#'
#' This mirrors the old pipeline exactly:
#' - labels at Cutoffs!B4:B6
#' - values at Cutoffs!C4:C6
#' - named regions point to C4/C5/C6
#' - green/red styles mark which cutoff is used (p-value vs FDR)
#' Add legacy Cutoffs sheet + named regions (PVAL_CO/FDR_CO/LFC_CO)
#'
#' Generic version for any omics mode (proteomics / rna / metabolomics ...)
#'
#' @param wb openxlsx workbook
#' @param config full config
#' @param mode string; e.g. "proteomics", "rna"
#' @param sheet sheet name (default "Cutoffs")
add_cutoffs_sheet_legacy <- function(wb,
                                     config,
                                     mode = "proteomics",
                                     sheet = "Cutoffs") {
  stopifnot(requireNamespace("openxlsx", quietly = TRUE))
  
  de_cfg <- config$modes[[mode]]$de
  if (is.null(de_cfg)) stop("No de config found for mode: ", mode)
  
  FDR_ADJ <- isTRUE(de_cfg$use_adj_for_pass1)
  
  P_CUTOFF         <- de_cfg$p_cutoff %||% 0.05
  LINEAR_FC_CUTOFF <- de_cfg$linear_fc_cutoff %||% 1.5
  
  # recreate sheet to avoid stale named regions
  if (sheet %in% names(wb)) openxlsx::removeWorksheet(wb, sheet)
  openxlsx::addWorksheet(wb, sheetName = sheet, gridLines = TRUE)
  
  # labels (B4:B6)
  openxlsx::writeData(
    wb, sheet,
    x = c("p-value", "Adjusted pvalue (FDR)", "linear Fold Change (linearFC)"),
    startCol = 2, startRow = 4,
    colNames = FALSE, rowNames = FALSE
  )
  
  # values (C4:C6)
  openxlsx::writeData(
    wb, sheet,
    x = c(
      ifelse(FDR_ADJ, "", P_CUTOFF),
      ifelse(FDR_ADJ, P_CUTOFF, ""),
      LINEAR_FC_CUTOFF
    ),
    startCol = 3, startRow = 4,
    colNames = FALSE, rowNames = FALSE
  )
  
  # named regions (MUST match legacy)
  openxlsx::createNamedRegion(wb, sheet = sheet, cols = 3, rows = 4, name = "PVAL_CO")
  openxlsx::createNamedRegion(wb, sheet = sheet, cols = 3, rows = 5, name = "FDR_CO")
  openxlsx::createNamedRegion(wb, sheet = sheet, cols = 3, rows = 6, name = "LFC_CO")
  
  # styles
  style_used <- openxlsx::createStyle(
    border = "TopBottomLeftRight",
    borderStyle = "thick",
    fgFill = "green",
    halign = "center"
  )
  style_not_used <- openxlsx::createStyle(
    border = "TopBottomLeftRight",
    borderStyle = "thick",
    fgFill = "red",
    halign = "center"
  )
  
  if (FDR_ADJ) {
    openxlsx::addStyle(wb, sheet, style_not_used, rows = 4, cols = 3, stack = FALSE)
    openxlsx::addStyle(wb, sheet, style_used,     rows = 5:6, cols = 3, stack = FALSE)
  } else {
    openxlsx::addStyle(wb, sheet, style_not_used, rows = 5, cols = 3, stack = FALSE)
    openxlsx::addStyle(wb, sheet, style_used,     rows = c(4, 6), cols = 3, stack = FALSE)
  }
  
  openxlsx::setColWidths(wb, sheet, cols = 2, widths = "auto")
  
  invisible(TRUE)
}



