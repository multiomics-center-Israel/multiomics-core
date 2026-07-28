# Tests for check_metadata_delimiter_safety().
#
# Regression test for a real failure: a tab-separated sample sheet whose leading
# path column carried a comma parsed fine in the pipeline (read_table_auto
# routes by extension) but killed the report render, because the template read
# it with read.csv(). read.csv() saw one field in the header and two in every
# data row, inferred column 1 as row names, found them all identical, and
# aborted with "duplicate 'row.names' are not allowed" — with no traceback
# pointing back at the sample sheet.
#
# The check must fire at load time, name the offending columns, and stay quiet
# for sheets that are actually fine. All fixtures below are synthetic.

# Synthetic sheet; `source` mimics the shape of a path that carries a comma.
make_sheet <- function(source = c("run-1", "run-2", "run-3")) {
    data.frame(
        SampleName = c("S1", "S2", "S3"),
        Group      = c("ctrl", "trt", "trt"),
        Source     = source,
        stringsAsFactors = FALSE
    )
}

write_sheet <- function(df, ext) {
    path <- withr::local_tempfile(fileext = paste0(".", ext), .local_envir = parent.frame())
    sep <- if (ext == "csv") "," else "\t"
    utils::write.table(df, path, sep = sep, row.names = FALSE, quote = FALSE)
    path
}

test_that("commas in a tab-separated sheet raise an actionable warning", {
    df <- make_sheet(source = c("a,b", "c,d", "e,f"))

    expect_warning(
        check_metadata_delimiter_safety(df, "samples.txt", "proteomics"),
        "Source \\(3/3 rows\\)"
    )
    # The message must name the mode and the downstream symptom, so whoever
    # hits it can connect the warning to the render failure.
    w <- tryCatch(check_metadata_delimiter_safety(df, "samples.txt", "proteomics"),
                  warning = function(e) conditionMessage(e))
    expect_match(w, "\\[proteomics\\]")
    expect_match(w, "duplicate 'row.names' are not allowed")
})

test_that("only the columns that actually contain commas are reported", {
    df <- make_sheet(source = c("a,b", "clean", "clean"))

    w <- tryCatch(check_metadata_delimiter_safety(df, "samples.tsv", "rna"),
                  warning = function(e) conditionMessage(e))
    expect_match(w, "Source \\(1/3 rows\\)")
    expect_false(grepl("SampleName", w))
    expect_false(grepl("Group", w))
})

test_that("cell values never appear in the warning", {
    df <- make_sheet(source = c("subject-042,visit-3", "b,c", "d,e"))

    w <- tryCatch(check_metadata_delimiter_safety(df, "samples.txt", "proteomics"),
                  warning = function(e) conditionMessage(e))
    expect_false(grepl("subject-042", w, fixed = TRUE))
})

test_that("a clean tab-separated sheet is silent", {
    expect_silent(check_metadata_delimiter_safety(make_sheet(), "samples.txt", "proteomics"))
})

test_that("a genuine CSV is silent even with commas, since they are quoted", {
    df <- make_sheet(source = c("a,b", "c,d", "e,f"))
    expect_silent(check_metadata_delimiter_safety(df, "samples.csv", "proteomics"))
})

test_that("empty and non-data-frame inputs are tolerated", {
    expect_silent(check_metadata_delimiter_safety(make_sheet()[0, ], "samples.txt", "x"))
    expect_silent(check_metadata_delimiter_safety(NULL, "samples.txt", "x"))
})

test_that("the check catches the exact sheet shape that broke read.csv()", {
    # Reproduce the original failure end to end. Two details matter, and both
    # come from the real sheet: the comma-bearing column is FIRST, and the text
    # before the comma is identical on every row (a shared directory path).
    # read.csv() then sees header=1 field vs data=2 fields, takes column 1 as
    # row names, and finds them all identical. Move the column or vary the
    # prefix and the same file loads without complaint — which is exactly why
    # this failure mode is easy to miss.
    df <- data.frame(
        FileName   = paste0("reports/batch-1,batch-2/", c("a", "b", "c"), ".raw"),
        SampleName = c("S1", "S2", "S3"),
        Group      = c("ctrl", "trt", "trt"),
        stringsAsFactors = FALSE
    )
    path <- write_sheet(df, "txt")

    expect_error(utils::read.csv(path), "duplicate 'row.names' are not allowed")
    expect_warning(check_metadata_delimiter_safety(df, path, "proteomics"), "FileName")
})
