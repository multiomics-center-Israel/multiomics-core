# tests/testthat/test-metabolomics-precomputed-de.R
#
# Reading a metabolomics DE table another run already wrote. The reader used to
# recognise only an unnamed first column, so this pipeline's own export -- whose
# id column is named `feature_id` -- fell through to rownames() and keyed every
# feature on its row number: nothing downstream matched, and nothing said so.
#
# The loader now fails fast on anything it cannot read unambiguously -- no id
# column, missing, blank or duplicated ids, no fold-change or raw p-value
# column, the wide summary shape -- and records where each contrast's adjusted
# p-values came from.
#
# Synthetic tables only.

write_de <- function(df, path = tempfile(fileext = ".tsv"), ...) {
    utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, ...)
    path
}

de_frame <- function(id_name = "feature_id", padj = TRUE, ids = c("F1", "F2", "F3")) {
    df <- data.frame(
        id = ids,
        Name = c("Alpha", "Beta", "Gamma"),
        logFC = c(1.5, -0.2, 0.9),
        AveExpr = c(20, 21, 19),
        P.Value = c(0.001, 0.4, 0.02),
        adj.P.Val = c(0.02, 0.5, 0.09),
        stringsAsFactors = FALSE)
    names(df)[1] <- id_name
    if (!padj) df$adj.P.Val <- NULL
    df
}

de_export <- function(...) write_de(de_frame(...))

de_config <- function(path, project_dir = tempdir()) {
    list(project = list(dir = project_dir), paths = list(raw = "data"),
         modes = list(metabolomics = list(
             files = list(de_table = as.list(path)),
             de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5,
                       use_adjusted_pval = TRUE))))
}

load_quietly <- function(path, ...) {
    suppressMessages(load_precomputed_metabolomics_de(de_config(path, ...)))
}

# ---- feature ids -------------------------------------------------------------

test_that("a named feature id column is read, not the row numbers", {
    res <- load_quietly(de_export())
    tbl <- res$de_tables[[1]]

    expect_identical(tbl$feature_id, c("F1", "F2", "F3"))
    expect_equal(tbl$logFC, c(1.5, -0.2, 0.9))
    expect_identical(res$summary_df$feature_id, c("F1", "F2", "F3"))
})

test_that("an unnamed first column still works", {
    # write.table with an empty name quotes nothing; readr names it "...1".
    res <- load_quietly(de_export(id_name = ""))
    expect_identical(res$de_tables[[1]]$feature_id, c("F1", "F2", "F3"))
})

test_that("a named id column wins over an unnamed row-number column", {
    # row.names = TRUE writes an unnamed first column holding 1, 2, 3 beside
    # the real feature_id column -- taking it would mis-key every feature.
    path <- tempfile(fileext = ".tsv")
    utils::write.table(de_frame(), path, sep = "\t", quote = FALSE,
                       row.names = TRUE, col.names = NA)
    res <- load_quietly(path)
    expect_identical(res$de_tables[[1]]$feature_id, c("F1", "F2", "F3"))
})

test_that("a table with no id column at all is refused by name", {
    path <- write_de(data.frame(logFC = 1, P.Value = 0.01))
    expect_error(load_quietly(path), "no feature id column")
    expect_error(load_quietly(path), basename(path), fixed = TRUE)
})

test_that("a later column named X or V1 is not taken as the feature id", {
    # Only column 1 can be the unnamed/index-like fallback; here it is Name.
    df <- data.frame(Name = c("Alpha", "Beta"), logFC = c(1.5, -0.2),
                     P.Value = c(0.001, 0.4), X = c("F1", "F2"), V1 = c("a", "b"),
                     stringsAsFactors = FALSE)
    expect_error(load_quietly(write_de(df)), "no feature id column")
})

test_that("a missing feature id is refused, not dropped", {
    path <- de_export(ids = c("F1", NA, "F3"))
    expect_error(load_quietly(path), "missing feature id in column 'feature_id'")
})

test_that("a blank feature id is refused, not dropped", {
    path <- tempfile(fileext = ".csv")
    df <- de_frame(ids = c("F1", "  ", "F3"))
    utils::write.csv(df, path, row.names = FALSE)
    # readr may read a whitespace-only field as empty or as NA; either way the
    # row is refused rather than kept under a blank key.
    expect_error(load_quietly(path), "(missing|empty) feature id")
})

test_that("duplicated feature ids are refused, not de-duplicated", {
    path <- de_export(ids = c("F1", "F2", "F1"))
    expect_error(load_quietly(path), "appear more than once \\(e\\.g\\. F1\\)")
})

# ---- statistics columns ------------------------------------------------------

test_that("a table with no fold-change column is refused", {
    df <- de_frame()
    df$logFC <- NULL
    expect_error(load_quietly(write_de(df)), "no log2 fold-change column")
})

test_that("a table with no raw p-value column is refused", {
    df <- de_frame()
    df$P.Value <- NULL
    expect_error(load_quietly(write_de(df)), "no raw p-value column")
})

test_that("the wide de_summary shape is refused, not read as all-NA", {
    wide <- data.frame(feature_id = c("F1", "F2"),
                       `linearFC.A_vs_B` = c(2.8, -1.1),
                       `pvalue.A_vs_B` = c(0.001, 0.4),
                       `padj.A_vs_B` = c(0.01, 0.5),
                       check.names = FALSE)
    expect_error(load_quietly(write_de(wide)), "wide DE summary")
})

# ---- adjusted p-values and their provenance ----------------------------------

test_that("the file's own adjusted p-value is kept, and recorded as input", {
    res <- load_quietly(de_export())
    tbl <- res$de_tables[[1]]

    # BH over these three raw p-values would give 0.0015, 0.4, 0.03 -- the point
    # is that the numbers come from the run that wrote the table.
    expect_equal(tbl$adj.P.Val, c(0.02, 0.5, 0.09))
    expect_identical(res$padj_provenance$padj_source, "input")
    expect_identical(res$padj_provenance$padj_column, "adj.P.Val")
})

test_that("with no adjusted column the reader applies BH, and records it", {
    res <- load_quietly(de_export(padj = FALSE))
    expect_equal(res$de_tables[[1]]$adj.P.Val,
                 stats::p.adjust(c(0.001, 0.4, 0.02), method = "BH"))
    expect_identical(res$padj_provenance$padj_source, "computed_bh")
    expect_true(is.na(res$padj_provenance$padj_column))
})

test_that("provenance is recorded per contrast", {
    dir <- tempfile()
    dir.create(dir)
    a <- write_de(de_frame(), file.path(dir, "de_A.tsv"))
    b <- write_de(de_frame(padj = FALSE), file.path(dir, "de_B.tsv"))
    res <- load_quietly(c(a, b))
    expect_identical(res$padj_provenance$contrast, c("A", "B"))
    expect_identical(res$padj_provenance$padj_source, c("input", "computed_bh"))
})

test_that("the Methods sentence follows the recorded provenance", {
    prov <- function(...) data.frame(padj_source = c(...), stringsAsFactors = FALSE)
    expect_identical(describe_padj_provenance(prov("computed_bh", "computed_bh")),
                     "Adjusted p-values were computed using the Benjamini-Hochberg procedure.")
    expect_match(describe_padj_provenance(prov("input")),
                 "taken from the precomputed input tables and were not recomputed")
    expect_match(describe_padj_provenance(prov("input", "computed_bh")),
                 "preserved; where absent, they were computed using the Benjamini-Hochberg")
    # Nothing recorded, or something unrecognised: say nothing rather than guess.
    expect_null(describe_padj_provenance(NULL))
    expect_null(describe_padj_provenance(prov(character(0))))
    expect_null(describe_padj_provenance(prov("input", "somewhere")))
})

# ---- paths -------------------------------------------------------------------

test_that("an absolute path is read as given", {
    path <- de_export()
    expect_silent(load_quietly(path, project_dir = "/nonexistent"))
})

test_that("a relative path is resolved under paths.raw, not the project root", {
    project <- tempfile()
    dir.create(file.path(project, "data", "prev"), recursive = TRUE)
    write_de(de_frame(), file.path(project, "data", "prev", "de_A.tsv"))
    res <- load_quietly("prev/de_A.tsv", project_dir = project)
    expect_identical(names(res$de_tables), "A")

    # The same relative path is not looked up at the project root.
    dir.create(file.path(project, "elsewhere"))
    write_de(de_frame(), file.path(project, "elsewhere", "de_B.tsv"))
    expect_error(load_quietly("elsewhere/de_B.tsv", project_dir = project),
                 "not found: .*data/elsewhere/de_B\\.tsv")
})
