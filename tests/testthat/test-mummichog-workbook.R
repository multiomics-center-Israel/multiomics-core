# tests/testthat/test-mummichog-workbook.R
#
# The all-contrast mummichog workbook (06v): every sheet is built from the
# report sections, so its rows and values must equal the objects the HTML
# report shows. Built from ma_two_contrast_sections() (helper-mummichog-ma.R).
# All data SYNTHETIC.

wb_sheets <- function() {
  tc <- ma_two_contrast_sections()
  list(tc = tc, sheets = build_mummichog_workbook_sheets(tc$sections))
}

strip_contrast <- function(sheet, ct) {
  x <- sheet[sheet$Contrast == ct, setdiff(names(sheet), c("Contrast", "Evidence status")),
             drop = FALSE]
  rownames(x) <- NULL
  x
}

test_that("the workbook has the agreed sheets, each detailed one with a Contrast column", {
  sh <- wb_sheets()$sheets
  expect_named(sh, c("README", "README_input", "Across Contrasts", "Pathway Overview",
                     "ORA Results", "GSEA Results", "EC Evidence", "Feature Evidence"))
  for (nm in c("Pathway Overview", "ORA Results", "GSEA Results", "EC Evidence",
               "Feature Evidence", "README_input")) {
    expect_identical(names(sh[[nm]])[1], "Contrast", label = nm)
    expect_setequal(unique(sh[[nm]]$Contrast), c("HL_vs_LL", "LL_vs_DD"))
  }
})

test_that("every sheet equals the section objects the report shows -- all rows", {
  w <- wb_sheets(); secs <- w$tc$sections; sh <- w$sheets
  for (ct in names(secs)) {
    s <- secs[[ct]]
    ov <- s$overview; rownames(ov) <- NULL
    expect_identical(strip_contrast(sh[["Pathway Overview"]], ct), ov)
    od <- s$ora_detail; rownames(od) <- NULL
    expect_identical(strip_contrast(sh[["ORA Results"]], ct), od)
    ec <- s$evidence$ec_table; rownames(ec) <- NULL
    expect_identical(strip_contrast(sh[["EC Evidence"]], ct), ec)
    ft <- s$evidence$feature_table; rownames(ft) <- NULL
    expect_identical(strip_contrast(sh[["Feature Evidence"]], ct), ft)
    g <- strip_contrast(sh[["GSEA Results"]], ct)
    expect_equal(g[["GSEA raw p"]], s$gsea$table$P.Value)
    expect_equal(g[["GSEA BH FDR"]], s$gsea$table$padj)
  }
  # Not limited to the HTML's top 20: the overview holds every pathway,
  # including the GSEA-only pathway with no significant EC.
  expect_true("TCA cycle" %in% sh[["Pathway Overview"]]$Pathway)
  expect_identical(sh[["Across Contrasts"]],
                   build_mummichog_across_contrasts(secs)$wide)
})

test_that("GSEA columns are named for raw p and BH FDR; ORA keeps Holm apart", {
  sh <- wb_sheets()$sheets
  expect_true(all(c("GSEA raw p", "GSEA BH FDR") %in% names(sh[["GSEA Results"]])))
  expect_false(any(c("P.Value", "padj") %in% names(sh[["GSEA Results"]])))
  expect_true(all(c("ORA raw p (FET)", "Holm-adjusted FET",
                    "Putative compound matches, overlap ECs (Hits.sig)") %in%
                    names(sh[["ORA Results"]])))
})

test_that("repeated features and repeated candidates are preserved, not collapsed", {
  sh <- wb_sheets()$sheets
  fe <- sh[["Feature Evidence"]]
  hl <- fe[fe$Contrast == "HL_vs_LL", ]
  # F01 sits in EC0001 and EC0003 of Glycolysis: two rows.
  expect_identical(sum(hl$Feature == "F01" & hl$Pathway == "Glycolysis"), 2L)
  # EC0002 supports both pathways in HL_vs_LL, each with its own candidates.
  ec <- sh[["EC Evidence"]]
  e2 <- ec[ec$EmpiricalCompound == "EC0002" & ec$Contrast == "HL_vs_LL", ]
  expect_setequal(e2[["Pathway-matching candidate(s)"]], c("C00022", "C00022; C00186"))
})

test_that("evidence status and input counts reach the README", {
  sh <- wb_sheets()$sheets
  inp <- sh[["README_input"]]
  expect_identical(nrow(inp), 2L)
  expect_true(all(inp[["Evidence status"]] == "Evidence built"))
  guide <- paste(sh$README$Explanation, collapse = " ")
  expect_match(guide, "pass_any_contrast", fixed = TRUE)
  expect_match(guide, "not a Mummichog filter", fixed = TRUE)
  expect_match(guide, "only -- those with at least one significant EC", fixed = TRUE)
  expect_match(guide, "Benjamini-Hochberg", fixed = TRUE)
  expect_match(guide, "not additive", fixed = TRUE)
})

test_that("before writing, statistics are numbers, identifiers text, missing values NA", {
  sh <- wb_sheets()$sheets
  num <- list(
    "Pathway Overview" = c("ORA raw p", "ORA Holm-adjusted p", "GSEA raw p", "GSEA BH FDR",
                           "ORA overlap ECs (Sig.EC)", "Features in overlap ECs"),
    "ORA Results"      = c("ORA raw p (FET)", "Holm-adjusted FET",
                           "Putative compound matches, overlap ECs (Hits.sig)"),
    "GSEA Results"     = c("GSEA raw p", "GSEA BH FDR"),
    "Feature Evidence" = c("Feature p-value", "Feature statistic (MetaboAnalyst input)",
                           "Feature log2FC (DE)", "m/z", "RT"))
  for (nm in names(num)) for (cl in num[[nm]]) {
    expect_true(is.numeric(sh[[nm]][[cl]]), label = paste0(nm, ": ", cl))
  }
  txt <- list(
    "Pathway Overview" = "Pathway ID",
    "ORA Results"      = c("Pathway ID", "Overlap EC IDs"),
    "EC Evidence"      = c("EmpiricalCompound", "Pathway-matching candidate(s)",
                           "Candidate KEGG ID(s)"),
    "Feature Evidence" = c("Feature", "EmpiricalCompound", "Matched compound(s)",
                           "Pathway candidate(s) matched by feature"))
  for (nm in names(txt)) for (cl in txt[[nm]]) {
    expect_true(is.character(sh[[nm]][[cl]]), label = paste0(nm, ": ", cl))
  }
  # The GSEA-only pathway has no ORA p-value: a real NA, not a string.
  ov <- sh[["Pathway Overview"]]
  expect_true(is.na(ov[["ORA raw p"]][ov$Pathway == "TCA cycle"]))
})

test_that("the written workbook keeps numbers numeric, ids as text, headers frozen and filtered", {
  skip_if_not_installed("openxlsx")
  w <- wb_sheets()
  path <- file.path(withr::local_tempdir(), "mm.xlsx")
  write_mummichog_workbook(w$sheets, path)
  expect_true(file.exists(path))
  expect_setequal(openxlsx::getSheetNames(path),
                  c("README", "Across Contrasts", "Pathway Overview", "ORA Results",
                    "GSEA Results", "EC Evidence", "Feature Evidence"))
  # read.xlsx() replaces blanks in column names with "." by default (sep.names)
  # and drops all-empty columns (skipEmptyCols); both would hide the columns
  # this test checks, so the workbook is read back as written.
  rx <- function(s) openxlsx::read.xlsx(path, sheet = s, check.names = FALSE,
                                        sep.names = " ", skipEmptyCols = FALSE)
  fe <- rx("Feature Evidence")
  expect_true("0012" %in% fe$Feature)                         # leading zeros kept
  expect_true(is.numeric(fe[["Feature p-value"]]))
  ov <- rx("Pathway Overview")
  expect_true(is.numeric(ov[["ORA raw p"]]))
  expect_identical(nrow(ov), nrow(w$sheets[["Pathway Overview"]]))

  # Every data sheet round-trips: same columns, numbers back as numbers (to
  # Excel's 15 significant digits), text as text, missing values as missing.
  for (nm in c("Across Contrasts", "Pathway Overview", "ORA Results", "GSEA Results",
               "EC Evidence", "Feature Evidence")) {
    orig <- w$sheets[[nm]]
    back <- rx(nm)
    expect_identical(names(back), names(orig), label = nm)
    expect_identical(nrow(back), nrow(orig), label = nm)
    for (cl in names(orig)) {
      lab <- paste0(nm, ": ", cl)
      if (all(is.na(orig[[cl]]))) {
        expect_true(all(is.na(back[[cl]])), label = lab)
      } else if (is.numeric(orig[[cl]])) {
        expect_true(is.numeric(back[[cl]]), label = lab)
        expect_equal(back[[cl]], as.numeric(orig[[cl]]), tolerance = 1e-12, label = lab)
      } else if (is.logical(orig[[cl]])) {
        expect_equal(as.logical(back[[cl]]), orig[[cl]], label = lab)
      } else {
        expect_identical(back[[cl]], as.character(orig[[cl]]), label = lab)
      }
    }
  }

  xml_dir <- withr::local_tempdir()
  utils::unzip(path, exdir = xml_dir)
  sheets_xml <- list.files(file.path(xml_dir, "xl", "worksheets"), "\\.xml$",
                           full.names = TRUE)
  txt <- vapply(sheets_xml, function(f) paste(readLines(f, warn = FALSE), collapse = ""),
                character(1))
  # Six data sheets, each with a filter and a frozen header row.
  expect_identical(sum(grepl("<autoFilter", txt, fixed = TRUE)), 6L)
  expect_identical(sum(grepl("state=\"frozen\"", txt, fixed = TRUE)), 6L)
})

test_that("the exports target writes the workbook beside the CSVs and touches nothing else", {
  skip_if_not_installed("openxlsx")
  tc <- ma_two_contrast_sections()
  out <- withr::local_tempdir()
  # A stand-in for the DE workbook in the same output directory.
  de_xlsx <- file.path(out, "Final_results_DE_P_0.05.xlsx")
  writeLines("not a real workbook", de_xlsx)
  before <- tools::md5sum(c(de_xlsx, tc$a$files, tc$b$files))
  written <- mod_mummichog_exports(tc$sections, c(tc$a$files, tc$b$files), out)
  wb <- file.path(out, "mummichog_pinned", .MMC_WORKBOOK_FILE)
  expect_true(wb %in% written)
  expect_true(all(dirname(written) == file.path(out, "mummichog_pinned")))
  expect_false(de_xlsx %in% written)
  expect_identical(tools::md5sum(c(de_xlsx, tc$a$files, tc$b$files)), before)
  # The per-contrast CSV/TSV exports are all still written.
  for (slug in c("HL_vs_LL", "LL_vs_DD")) {
    expect_true(all(file.exists(file.path(out, "mummichog_pinned", sprintf(
      c("mummichog_summary_%s.tsv", "mummichog_pathway_table_%s.csv",
        "mummichog_evidence_features_%s.tsv", "mummichog_gsea_table_%s.tsv"), slug)))))
  }
  # No sections: a stale workbook is removed rather than left reading as current.
  mod_mummichog_exports(list(), c(tc$a$files, tc$b$files), out)
  expect_false(file.exists(wb))
})
