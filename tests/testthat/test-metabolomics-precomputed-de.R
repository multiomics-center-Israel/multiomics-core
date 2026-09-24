# tests/testthat/test-metabolomics-precomputed-de.R
#
# Reading a metabolomics DE table another run already wrote. The reader used to
# recognise only an unnamed first column, so this pipeline's own export -- whose
# id column is named `feature_id` -- fell through to rownames() and keyed every
# feature on its row number: nothing downstream matched, and nothing said so.
#
# Synthetic tables only.

de_export <- function(id_name = "feature_id", padj = TRUE) {
    df <- data.frame(
        id = c("F1", "F2", "F3"),
        Name = c("Alpha", "Beta", "Gamma"),
        logFC = c(1.5, -0.2, 0.9),
        AveExpr = c(20, 21, 19),
        P.Value = c(0.001, 0.4, 0.02),
        adj.P.Val = c(0.02, 0.5, 0.09),
        stringsAsFactors = FALSE)
    names(df)[1] <- id_name
    if (!padj) df$adj.P.Val <- NULL
    path <- tempfile(fileext = ".tsv")
    utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
    path
}

de_config <- function(path) {
    list(project = list(dir = tempdir()), paths = list(raw = "data"),
         modes = list(metabolomics = list(
             files = list(de_table = list(path)),
             de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5,
                       use_adjusted_pval = TRUE))))
}

test_that("a named feature id column is read, not the row numbers", {
    res <- suppressMessages(load_precomputed_metabolomics_de(de_config(de_export())))
    tbl <- res$de_tables[[1]]

    expect_identical(tbl$feature_id, c("F1", "F2", "F3"))
    expect_equal(tbl$logFC, c(1.5, -0.2, 0.9))
    expect_identical(res$summary_df$feature_id, c("F1", "F2", "F3"))
})

test_that("an unnamed first column still works", {
    path <- de_export(id_name = "")
    # write.table with an empty name quotes nothing; readr names it "...1".
    res <- suppressMessages(load_precomputed_metabolomics_de(de_config(path)))
    expect_identical(res$de_tables[[1]]$feature_id, c("F1", "F2", "F3"))
})

test_that("the file's own adjusted p-value is kept, not recomputed", {
    res <- suppressMessages(load_precomputed_metabolomics_de(de_config(de_export())))
    tbl <- res$de_tables[[1]]

    # BH over these three raw p-values would give 0.0015, 0.4, 0.03 -- the point
    # is that the numbers come from the run that wrote the table.
    expect_equal(tbl$adj.P.Val, c(0.02, 0.5, 0.09))
})

test_that("with no adjusted column the reader adjusts, and says which rule", {
    res <- suppressMessages(load_precomputed_metabolomics_de(
        de_config(de_export(padj = FALSE))))
    expect_equal(res$de_tables[[1]]$adj.P.Val,
                 stats::p.adjust(c(0.001, 0.4, 0.02), method = "BH"))
})

test_that("a table with no id column at all is refused by name", {
    df <- data.frame(logFC = 1, P.Value = 0.01)
    path <- tempfile(fileext = ".tsv")
    utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)

    expect_error(suppressMessages(load_precomputed_metabolomics_de(de_config(path))),
                 "no feature id column")
})

test_that("an absolute path is read as given", {
    path <- de_export()
    cfg <- de_config(path)
    cfg$project$dir <- "/nonexistent"
    expect_silent(suppressMessages(load_precomputed_metabolomics_de(cfg)))
})
