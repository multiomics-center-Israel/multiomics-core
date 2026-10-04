# Tests for the Cutoffs sheet and the workbook zip written by
# R/core/05_export_excel.R.
#
#   testthat::test_file("tests/testthat/test-export-excel-cutoffs.R")
#
# Test map:
#   C1 de_uses_adjusted_p resolves all three config spellings + per-mode default
#   C2 a saved workbook has no relationship or [Content_Types] Override
#      pointing at a part the archive does not contain
#   C3 the FDR cutoff cell is populated and NUMERIC for mode "rna"
#   C4 the manual_cutoffs formulas reference FDR_CO (not PVAL_CO) for "rna"
#   C5 the pruning is inert: a workbook holding a real image keeps its drawing
#
# All fixtures are synthetic (g1/g2/g3, S1..S4) — no project data.

`%||%` <- function(a, b) if (is.null(a)) b else a

find_repo_file <- function(rel) {
    dir <- normalizePath(".", mustWork = FALSE)
    for (i in 1:8) {
        cand <- file.path(dir, rel)
        if (file.exists(cand)) return(cand)
        dir <- dirname(dir)
    }
    stop("Could not locate ", rel, " from working dir ", getwd())
}

source(find_repo_file("R/core/05_export_excel.R"))

# ---- helpers ---------------------------------------------------------------

unpack_xlsx <- function(path) {
    dir <- tempfile("unpacked_")
    dir.create(dir)
    utils::unzip(path, exdir = dir)
    dir
}

read_part <- function(dir, ...) paste(readLines(file.path(dir, ...), warn = FALSE), collapse = "")

#' Targets referenced by worksheet rels that are absent from the archive
dangling_rel_targets <- function(dir) {
    rels_dir <- file.path(dir, "xl", "worksheets", "_rels")
    if (!dir.exists(rels_dir)) return(character(0))
    out <- character(0)
    for (rf in list.files(rels_dir, pattern = "\\.rels$", full.names = TRUE)) {
        one <- paste(readLines(rf, warn = FALSE), collapse = "")
        rels <- regmatches(one, gregexpr("<Relationship\\b[^>]*/>", one, perl = TRUE))[[1]]
        if (length(rels) == 0L) next
        targets <- sub('.*Target="([^"]*)".*', "\\1", rels)
        # OPC: a rel target is relative to the part's own directory
        # (xl/worksheets/), not to the _rels/ folder holding the .rels file.
        resolved <- normalizePath(file.path(dirname(rels_dir), targets), mustWork = FALSE)
        out <- c(out, targets[!file.exists(resolved)])
    }
    out
}

#' [Content_Types] Overrides naming parts that are absent from the archive
dangling_overrides <- function(dir) {
    ct <- read_part(dir, "[Content_Types].xml")
    ovs <- regmatches(ct, gregexpr("<Override\\b[^>]*/>", ct, perl = TRUE))[[1]]
    parts <- sub('.*PartName="([^"]*)".*', "\\1", ovs)
    parts[!file.exists(file.path(dir, sub("^/", "", parts)))]
}

sheet_xml_path <- function(dir, sheet_name) {
    wbx <- read_part(dir, "xl", "workbook.xml")
    sheets <- regmatches(wbx, gregexpr("<sheet\\b[^>]*/>", wbx, perl = TRUE))[[1]]
    nms <- sub('.*name="([^"]*)".*', "\\1", sheets)
    rid <- sub('.*r:id="([^"]*)".*', "\\1", sheets)[match(sheet_name, nms)]
    rl <- read_part(dir, "xl", "_rels", "workbook.xml.rels")
    r <- regmatches(rl, gregexpr("<Relationship\\b[^>]*/>", rl, perl = TRUE))[[1]]
    tgt <- sub('.*Target="([^"]*)".*', "\\1", r)[match(rid, sub('.*Id="([^"]*)".*', "\\1", r))]
    file.path(dir, "xl", sub("^/", "", tgt))
}

cell_xml <- function(sheet_path, ref) {
    s <- paste(readLines(sheet_path, warn = FALSE), collapse = "")
    pat <- sprintf('<c r="%s"[^>]*/>|<c r="%s"[^>]*>.*?</c>', ref, ref)
    m <- regmatches(s, gregexpr(pat, s, perl = TRUE))[[1]]
    if (length(m) == 0L) NA_character_ else m[[1]]
}

synthetic_rna_results <- function() {
    data.frame(
        Gene                      = c("g1", "g2", "g3"),
        S1                        = c(100, 200, 300),
        S2                        = c(110, 190, 320),
        S3                        = c(150, 100, 310),
        S4                        = c(160, 105, 300),
        "log2FC.B_vs_A"           = c(0.60, -1.00, 0.01),
        "linearFC.B_vs_A"         = c(1.52, -2.00, 1.01),
        "pvalue.B_vs_A"           = c(0.001, 0.002, 0.900),
        "padj.B_vs_A"             = c(0.010, 0.020, 0.950),
        "upDown.B_vs_A"           = c("up", "down", ""),
        "manual_cutoffs.B_vs_A"   = NA,
        pass_any_contrast         = c(1, 1, NA),
        check.names               = FALSE,
        stringsAsFactors          = FALSE
    )
}

rna_config <- function(...) {
    list(modes = list(rna = list(de = c(
        list(p_cutoff = 0.05, linear_fc_cutoff = 1.5), list(...)
    ))))
}

# ---- C1 --------------------------------------------------------------------

test_that("de_uses_adjusted_p reads every mode's spelling of the switch", {
    expect_true(de_uses_adjusted_p(list(use_adj_for_pass1 = TRUE), "proteomics"))
    expect_false(de_uses_adjusted_p(list(use_adj_for_pass1 = FALSE), "proteomics"))
    expect_true(de_uses_adjusted_p(list(use_adj = TRUE), "rna"))
    expect_true(de_uses_adjusted_p(list(use_adjusted_pval = TRUE), "metabolomics"))
    expect_false(de_uses_adjusted_p(list(use_adjusted_pval = FALSE), "metabolomics"))
})

test_that("de_uses_adjusted_p falls back to what each DE step actually does", {
    # RNA gates on padj unconditionally; metabolomics defaults use_adjusted_pval
    # to TRUE; proteomics treats a missing use_adj_for_pass1 as FALSE.
    expect_true(de_uses_adjusted_p(list(p_cutoff = 0.05), "rna"))
    expect_true(de_uses_adjusted_p(list(p_cutoff = 0.05), "metabolomics"))
    expect_true(de_uses_adjusted_p(list(p_cutoff = 0.05), "lipidomics"))
    expect_false(de_uses_adjusted_p(list(p_cutoff = 0.05), "proteomics"))
    expect_false(de_uses_adjusted_p(NULL, "proteomics"))
})

test_that("de_uses_adjusted_p warns rather than misreport RNA as raw-p gated", {
    expect_warning(res <- de_uses_adjusted_p(list(use_adj = FALSE), "rna"),
                   "gates on padj")
    expect_true(res)
})

# ---- C2-C4: one real export, several assertions ----------------------------

rna_export <- function(config) {
    out_dir <- tempfile("xlsx_out_")
    dir.create(out_dir)
    files <- write_final_results_excels_legacy_generic(
        final_results = synthetic_rna_results(),
        config        = config,
        out_dir       = out_dir,
        mode          = "rna",
        id_col        = "Gene",
        expr_for_de   = NULL,
        with_cutoffs  = TRUE
    )
    files
}

test_that("saved workbooks reference no missing parts", {
    skip_if_not_installed("openxlsx")
    files <- rna_export(rna_config(use_adj = TRUE))
    expect_length(files, 2L)

    for (f in files) {
        dir <- unpack_xlsx(f)
        expect_equal(dangling_rel_targets(dir), character(0),
                     info = paste("dangling worksheet rels in", basename(f)))
        expect_equal(dangling_overrides(dir), character(0),
                     info = paste("dangling Content_Types Overrides in", basename(f)))
    }
})

test_that("the FDR cutoff cell is populated and numeric for mode 'rna'", {
    skip_if_not_installed("openxlsx")
    f_all <- rna_export(rna_config(use_adj = TRUE))[[1]]
    dir <- unpack_xlsx(f_all)
    cutoffs <- sheet_xml_path(dir, "Cutoffs")

    fdr_cell <- cell_xml(cutoffs, "C5")   # FDR_CO
    lfc_cell <- cell_xml(cutoffs, "C6")   # LFC_CO
    expect_false(is.na(fdr_cell))
    expect_true(grepl("<v>", fdr_cell, fixed = TRUE))
    # t="s" would make it a shared string, and Excel sorts every number below
    # every string, so p<=FDR_CO would be true for every row.
    expect_false(grepl('t="s"', fdr_cell, fixed = TRUE))
    expect_false(grepl('t="str"', fdr_cell, fixed = TRUE))
    expect_false(grepl('t="s"', lfc_cell, fixed = TRUE))
})

test_that("manual_cutoffs formulas gate on FDR_CO for mode 'rna'", {
    skip_if_not_installed("openxlsx")
    f_all <- rna_export(rna_config(use_adj = TRUE))[[1]]
    dir <- unpack_xlsx(f_all)
    results <- paste(readLines(sheet_xml_path(dir, "Results"), warn = FALSE), collapse = "")
    expect_true(grepl("FDR_CO", results, fixed = TRUE))
    expect_false(grepl("PVAL_CO", results, fixed = TRUE))
})

test_that("a raw-p mode still gates on PVAL_CO", {
    skip_if_not_installed("openxlsx")
    cfg <- list(modes = list(proteomics = list(de = list(
        p_cutoff = 0.05, linear_fc_cutoff = 1.5, use_adj_for_pass1 = FALSE
    ))))
    res <- synthetic_rna_results()
    # proteomics column naming: the DE stats carry the .imputs. infix
    names(res) <- sub("^(log2FC|linearFC|pvalue|padj|upDown)\\.", "\\1.imputs.", names(res))
    out_dir <- tempfile("xlsx_prot_"); dir.create(out_dir)
    f_all <- write_final_results_excels_legacy_generic(
        final_results = res, config = cfg, out_dir = out_dir, mode = "proteomics",
        id_col = "Gene", expr_for_de = NULL, with_cutoffs = TRUE
    )[[1]]

    dir <- unpack_xlsx(f_all)
    results <- paste(readLines(sheet_xml_path(dir, "Results"), warn = FALSE), collapse = "")
    expect_true(grepl("PVAL_CO", results, fixed = TRUE))
    expect_false(grepl("FDR_CO", results, fixed = TRUE))
    expect_false(grepl('t="s"', cell_xml(sheet_xml_path(dir, "Cutoffs"), "C4"), fixed = TRUE))
    expect_equal(dangling_rel_targets(dir), character(0))
})

# ---- C5 --------------------------------------------------------------------

test_that("pruning is inert when the sheet really holds a drawing", {
    skip_if_not_installed("openxlsx")
    skip_if_not(isTRUE(capabilities("png")), "no png device")

    img <- tempfile(fileext = ".png")
    grDevices::png(img, width = 200, height = 200)
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    grDevices::dev.off()

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "WithImage")
    openxlsx::addWorksheet(wb, "Plain")
    openxlsx::writeData(wb, "Plain", data.frame(a = 1:2))
    openxlsx::insertImage(wb, "WithImage", file = img, width = 1, height = 1)

    drop_dangling_drawing_rels(wb)
    expect_true(any(grepl("/relationships/drawing\"", wb$worksheets_rels[[1]], fixed = TRUE)))
    expect_true(any(grepl("/xl/drawings/drawing1.xml", wb$Content_Types, fixed = TRUE)))
    # the untouched sheet is still pruned
    expect_false(any(grepl("/relationships/drawing\"", wb$worksheets_rels[[2]], fixed = TRUE)))

    out <- tempfile(fileext = ".xlsx")
    openxlsx::saveWorkbook(wb, out, overwrite = TRUE)
    dir <- unpack_xlsx(out)
    expect_true(file.exists(file.path(dir, "xl", "drawings", "drawing1.xml")))
    expect_equal(dangling_rel_targets(dir), character(0))
    expect_equal(dangling_overrides(dir), character(0))
})
