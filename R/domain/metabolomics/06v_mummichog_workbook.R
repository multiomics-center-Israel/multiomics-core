# R/domain/metabolomics/06v_mummichog_workbook.R
#
# One Excel workbook for the whole mummichog analysis: every contrast, every
# pathway, every EC and feature of the evidence, with a guide.
#
# The sheets are built ONLY from the report sections (06e) -- the same overview,
# ORA detail, GSEA table, evidence and input counts the HTML report shows -- so
# the workbook and the report cannot disagree on a value or a definition. The
# existing per-contrast CSV/TSV exports are written as before; this workbook is
# an addition. It is not the differential-abundance workbook
# (Final_results_*.xlsx), which this file never reads or writes.
#
# Dependency: openxlsx (already in renv.lock, used by R/core/05_export_excel.R).

# The workbook's file name, under <metabolomics out dir>/mummichog_pinned/.
.MMC_WORKBOOK_FILE <- "mummichog_results_all_contrasts.xlsx"

# Excel refuses cells longer than this; longer text is cut with a marker.
.MMC_XLSX_MAX_CHARS <- 32000L

#' Stack per-contrast tables, with a leading Contrast column
#'
#' Columns are united across contrasts (a column absent from one contrast is
#' NA there), so one engine's tables stack whatever their shape.
#'
#' @param tables Named list, contrast -> data.frame or NULL.
#' @return A data.frame (zero rows when every table is NULL or empty).
#' @noRd
.mmc_stack_contrasts <- function(tables) {
  tables <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, tables)
  if (length(tables) == 0L) {
    return(data.frame(Contrast = character(0), stringsAsFactors = FALSE))
  }
  cols <- unique(unlist(lapply(tables, names), use.names = FALSE))
  parts <- lapply(names(tables), function(ct) {
    x <- flatten_list_columns(as.data.frame(tables[[ct]], check.names = FALSE,
                                            stringsAsFactors = FALSE))
    for (cl in setdiff(cols, names(x))) x[[cl]] <- NA
    data.frame(check.names = FALSE, stringsAsFactors = FALSE,
               Contrast = rep(ct, nrow(x)), x[, cols, drop = FALSE])
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

#' A contrast's evidence status, in words
#' @param s One report section.
#' @return One sentence.
#' @noRd
.mmc_evidence_status <- function(s) {
  if (!is.null(s$evidence)) {
    if (is.null(s$evidence_log2fc_note)) "Evidence built" else
      paste("Evidence built.", s$evidence_log2fc_note)
  } else if (!is.null(s$evidence_note)) {
    paste("Evidence not built:", s$evidence_note)
  } else {
    "Evidence not available for this contrast"
  }
}

#' The guide text of the workbook
#'
#' @param engine   `"metaboanalyst"` or `"python"`.
#' @param p_cutoff Nominal threshold.
#' @param fdr      Adjusted-p threshold.
#' @return A data.frame with `Topic` and `Explanation`.
#' @noRd
.mmc_workbook_guide <- function(engine, p_cutoff, fdr = .MMC_SUMMARY_FDR) {
  ma <- identical(engine, "metaboanalyst")
  rows <- c(
    "What this is" = paste0(
      "Mummichog pathway analysis of LC-MS features, one analysis per differential-",
      "abundance (DE) contrast, all contrasts in one workbook. Mummichog gives ",
      "pathway-level hypotheses from features without confirmed identities; its ",
      "candidate compounds are putative matches, not identified metabolites. ",
      "The DE results themselves are in the separate DE workbook, which this file ",
      "does not replace."),
    "What enters Mummichog" = paste0(
      "For each contrast, every feature that reached DE (the contrast's complete DE ",
      "table), not only the features passing the overall DE selection ",
      "(pass_any_contrast): that selection is not a Mummichog filter. ",
      if (ma) paste0(
        "The engine then leaves out features missing an m/z, retention time, p-value ",
        "or statistic, features with m/z outside 50-2001, and features whose m/z ",
        "exactly equals another's (it keeps the one with the smaller p-value, across ",
        "polarities). The README table below gives each contrast's counts, read from ",
        "the engine's run record.") else
        "The python engine's own input filtering applies; it records no counts."),
    "ORA" = paste0(
      "A feature is significant when its raw (unadjusted) DE p-value is <= ", p_cutoff,
      " -- not adjusted DE significance. An empirical compound (EC) is significant when ",
      "it holds at least one significant feature. Each pathway is tested for ",
      "over-representation of significant ECs",
      if (ma) paste0(
        " with Fisher's exact test (FET, 'ORA raw p'). 'ORA Holm-adjusted p' is ",
        "MetaboAnalyst's Holm adjustment, applied over the pathways in the ORA table ",
        "only -- those with at least one significant EC -- not over every pathway ",
        "tested. A pathway with no significant EC is not in the ORA table ('No ",
        "significant EC overlap').") else
        ". The python engine's ORA p-value is empirical; it reports no adjusted value."),
    "GSEA" = paste0(
      "Uses every EC, ranked by a signed statistic (an EC's score is the maximum of ",
      "its features'), with no significance cutoff. 'GSEA raw p' is the permutation ",
      "p-value and 'GSEA BH FDR' its Benjamini-Hochberg adjustment over the pathways ",
      "GSEA tested. The sign of NES is not a biological direction. A pathway can have ",
      "a GSEA signal without an ORA signal; its ORA evidence is then not shown as ",
      "support."),
    "Evidence levels" = paste0(
      "'Adjusted': ORA Holm-adjusted p or GSEA BH FDR < ", fdr, ". 'Nominal': raw p < ",
      p_cutoff, " without adjusted significance -- a hypothesis, not a finding. ORA ",
      "(Holm, family-wise) and GSEA (BH, false discovery rate) control different error ",
      "rates; their adjusted values are not interchangeable."),
    "Counts" = paste0(
      "'ORA overlap ECs (Sig.EC)' is the number of significant ECs in the pathway -- ",
      "the unit ORA counts. 'Features in overlap ECs' counts the distinct measured ",
      "features in those ECs, including features that are not significant",
      if (ma) paste0(" and features that matched no pathway candidate themselves ",
                     "(Feature Evidence: 'Matched via pathway candidate')") else "",
      ". One EC is not one metabolite. ",
      if (ma) "'Hits.*' columns count putative compound matches, not identified metabolites. " else "",
      "Counts are not additive across pathways or contrasts: the same EC and ",
      "feature can support several pathways."),
    "Annotation agreement" = paste0(
      "Compares a feature's dataset annotation with the pathway candidates",
      if (ma) " that feature itself matched" else " of its EC",
      " (KEGG, then HMDB, then name; never mass). Match, Conflict, Mixed (both, in one ",
      "EC) or Not assessed. It does not change any result."),
    "Across Contrasts" = paste0(
      "A lookup of existing results, aligned on the model's pathway ID. No test is ",
      "done across contrasts. A pathway recurring across contrasts does not mean it ",
      "changed in the same direction, nor that the same ECs or features support it; ",
      "compare them in the evidence sheets (filter on Contrast and Pathway)."),
    "Sheets" = paste0(
      "Across Contrasts (one row per pathway); Pathway Overview (one row per pathway ",
      "and contrast, every pathway, not only the report's top 20); ORA Results; GSEA ",
      "Results; EC Evidence; Feature Evidence. Each has a Contrast column and filters. ",
      "p-values are numbers; identifiers are text."))
  data.frame(Topic = names(rows), Explanation = unname(rows),
             stringsAsFactors = FALSE)
}

#' Build the sheets of the mummichog workbook from the report sections
#'
#' Pure: no file is read or written. Every value comes from the sections.
#'
#' @param sections Report sections (`build_mummichog_report_sections()`).
#' @return NULL when there are no sections; otherwise a named list of
#'   data.frames: `README` (guide), `README_input` (what entered the engine,
#'   per contrast), `Across Contrasts`, `Pathway Overview`, `ORA Results`,
#'   `GSEA Results`, `EC Evidence`, `Feature Evidence`.
#' @examples
#' build_mummichog_workbook_sheets(list())
build_mummichog_workbook_sheets <- function(sections) {
  if (length(sections) == 0L) return(NULL)
  s1 <- sections[[1]]
  engine <- s1$engine %||% "python"
  p_cut  <- s1$p_cutoff %||% 0.05

  status <- vapply(sections, .mmc_evidence_status, character(1))
  input <- .mmc_stack_contrasts(lapply(sections, `[[`, "input_counts"))
  if (nrow(input) > 0L) input[["Evidence status"]] <- unname(status[input$Contrast])

  overview <- .mmc_stack_contrasts(lapply(sections, `[[`, "overview"))
  if (nrow(overview) > 0L) overview[["Evidence status"]] <- unname(status[overview$Contrast])

  gsea <- .mmc_stack_contrasts(lapply(sections, function(s) s$gsea$table))
  names(gsea)[names(gsea) == "P.Value"] <- "GSEA raw p"
  names(gsea)[names(gsea) == "padj"]    <- "GSEA BH FDR"

  across <- build_mummichog_across_contrasts(sections)
  list(
    "README"           = .mmc_workbook_guide(engine, p_cut),
    "README_input"     = input,
    "Across Contrasts" = if (is.null(across)) data.frame() else across$wide,
    "Pathway Overview" = overview,
    "ORA Results"      = .mmc_stack_contrasts(lapply(sections, `[[`, "ora_detail")),
    "GSEA Results"     = gsea,
    "EC Evidence"      = .mmc_stack_contrasts(lapply(sections, function(s) s$evidence$ec_table)),
    "Feature Evidence" = .mmc_stack_contrasts(lapply(sections, function(s) s$evidence$feature_table))
  )
}

#' Write the mummichog workbook
#'
#' Header row bold and frozen, a filter on every data sheet, numbers kept
#' numeric and identifiers written as text (they are character columns). The
#' README sheet holds the guide, then the per-contrast input table.
#'
#' @param sheets Output of `build_mummichog_workbook_sheets()`.
#' @param path   Destination `.xlsx`.
#' @return `path`, invisibly.
#' @examples
#' \dontrun{
#' write_mummichog_workbook(build_mummichog_workbook_sheets(sections), "out.xlsx")
#' }
write_mummichog_workbook <- function(sheets, path) {
  wb <- .mmc_workbook_object(sheets)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}

#' The openxlsx workbook object for the sheets (separate for testing)
#' @param sheets Output of `build_mummichog_workbook_sheets()`.
#' @return An openxlsx Workbook.
#' @noRd
.mmc_workbook_object <- function(sheets) {
  clip <- function(df) {
    for (cl in names(df)) {
      if (is.character(df[[cl]])) {
        long <- !is.na(df[[cl]]) & nchar(df[[cl]]) > .MMC_XLSX_MAX_CHARS
        df[[cl]][long] <- paste0(substr(df[[cl]][long], 1L, .MMC_XLSX_MAX_CHARS),
                                 " [truncated; see the CSV/TSV exports]")
      }
    }
    df
  }
  wb <- openxlsx::createWorkbook()
  head_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#DCE6F1",
                                      border = "bottom", wrapText = TRUE,
                                      valign = "top")
  wrap <- openxlsx::createStyle(wrapText = TRUE, valign = "top")

  openxlsx::addWorksheet(wb, "README")
  guide <- sheets[["README"]]
  openxlsx::writeData(wb, "README", guide, headerStyle = head_style)
  openxlsx::addStyle(wb, "README", wrap, rows = seq_len(nrow(guide)) + 1L,
                     cols = 2, gridExpand = TRUE)
  openxlsx::setColWidths(wb, "README", cols = 1:2, widths = c(24, 110))
  input <- sheets[["README_input"]]
  if (is.data.frame(input) && nrow(input) > 0L) {
    at <- nrow(guide) + 3L
    openxlsx::writeData(wb, "README", "What entered Mummichog, per contrast",
                        startRow = at)
    openxlsx::addStyle(wb, "README", openxlsx::createStyle(textDecoration = "bold"),
                       rows = at, cols = 1)
    openxlsx::writeData(wb, "README", clip(input), startRow = at + 1L,
                        headerStyle = head_style)
  }

  for (nm in setdiff(names(sheets), c("README", "README_input"))) {
    df <- sheets[[nm]]
    openxlsx::addWorksheet(wb, nm)
    if (!is.data.frame(df) || ncol(df) == 0L) {
      openxlsx::writeData(wb, nm, "No results for this sheet.")
      next
    }
    openxlsx::writeData(wb, nm, clip(df), headerStyle = head_style,
                        withFilter = TRUE)
    openxlsx::freezePane(wb, nm, firstRow = TRUE)
    openxlsx::setColWidths(wb, nm, cols = seq_len(ncol(df)),
                           widths = ifelse(names(df) %in% c("Pathway", "Explanation"),
                                           40, 16))
  }
  wb
}
