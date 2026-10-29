# tests/testthat/test-rnaseq-report-de-called.R
#
# The RNA report's "Top DE genes" table and "Top DE" heatmaps read
# deseq2_summary_p*.tsv, which holds every tested gene. They ranked that whole
# table by padj, so the "All DE" heatmap reported every filtered gene as DE and
# showed genes below the fold-change cutoff. Both now keep only the rows the
# pipeline called DE (pass_any_contrast) via select_de_called_rows().

test_that("select_de_called_rows keeps only rows called DE", {
    df <- data.frame(Gene = c("g1", "g2", "g3", "g4"),
                     padj.A_vs_B = c(0.001, 0.002, 0.2, NA),
                     pass_any_contrast = c(TRUE, FALSE, FALSE, NA))
    out <- select_de_called_rows(df)
    expect_equal(out$Gene, "g1")
})

test_that("select_de_called_rows accepts 0/1 and text encodings", {
    num <- data.frame(Gene = c("g1", "g2"), pass_any_contrast = c(1L, 0L))
    expect_equal(select_de_called_rows(num)$Gene, "g1")
    txt <- data.frame(Gene = c("g1", "g2", "g3"),
                      pass_any_contrast = c("TRUE", "FALSE", "1"))
    expect_equal(select_de_called_rows(txt)$Gene, c("g1", "g3"))
})

test_that("select_de_called_rows returns zero rows when nothing passes", {
    df <- data.frame(Gene = c("g1", "g2"), pass_any_contrast = c(FALSE, FALSE))
    out <- select_de_called_rows(df)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 0L)
})

test_that("select_de_called_rows returns NULL without a DE call column", {
    expect_null(select_de_called_rows(NULL))
    expect_null(select_de_called_rows(data.frame(Gene = "g1", padj.A_vs_B = 0.01)))
})

test_that("the RNA report's DE table and heatmaps read the DE-called rows", {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "rnaseq", "report_template.Rmd"),
        "R/domain/rnaseq/report_template.Rmd"
    )
    f <- candidates[file.exists(candidates)][1]
    if (is.na(f)) stop("R/domain/rnaseq/report_template.Rmd not found from the test directory")
    src <- readLines(f, warn = FALSE)

    expect_true(any(grepl("^de_summary_pass <- select_de_called_rows\\(de_summary_full\\)", src)))
    expect_true(any(grepl("^heatmap_de_source <- de_summary_pass$", src)))
    expect_false(any(grepl("heatmap_de_source <- de_summary_full", src, fixed = TRUE)))
    expect_true(any(grepl("de_top <- head(de_summary_pass, 50)", src, fixed = TRUE)))
    # Gene Lookup stays on the full table: it is meant to find any tested gene.
    expect_true(any(grepl("de_summary_full[, display_cols]", src, fixed = TRUE)))
})
