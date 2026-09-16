# tests/testthat/test-report-samplesheet-read.R
#
# Regression test for the sample-sheet read in the proteomics report template.
# read.csv treats column 1 as row names whenever a data row has MORE fields
# than the header. A single unquoted comma in a free-text column is enough, and
# the result is silent: the frame keeps the right column NAMES but every value
# sits one column to the left. In the report that made SampleName hold the raw
# file path, so no expression column ever matched and the protein explorer
# rendered empty with no error anywhere.

write_ragged_samplesheet <- function() {
    f <- tempfile(fileext = ".csv")
    writeLines(c(
        "SeqNum,SampleName,Group,FileName,Explain",
        '1,S1_ctl,ctl,"D:\\runs\\a.raw",sample one, no treatment',
        '2,S1_trt,trt,"D:\\runs\\b.raw",sample one, treated'
    ), f)
    f
}

test_that("read.csv silently shifts columns on a ragged sample sheet", {
    f <- write_ragged_samplesheet(); on.exit(unlink(f), add = TRUE)
    bad <- read.csv(f, stringsAsFactors = FALSE)
    # The column names survive, which is what makes this so hard to notice.
    expect_true("SampleName" %in% names(bad))
    # But the values are wrong.
    expect_false(identical(as.character(bad$SampleName), c("S1_ctl", "S1_trt")))
})

test_that("row.names=NULL does not rescue read.csv here", {
    f <- write_ragged_samplesheet(); on.exit(unlink(f), add = TRUE)
    alt <- tryCatch(read.csv(f, stringsAsFactors = FALSE, row.names = NULL),
                    error = function(e) NULL)
    skip_if(is.null(alt), "read.csv errored outright")
    expect_false(identical(as.character(alt$SampleName), c("S1_ctl", "S1_trt")))
})

test_that("readr keeps the header's column count and the right values", {
    skip_if_not_installed("readr")
    f <- write_ragged_samplesheet(); on.exit(unlink(f), add = TRUE)
    ok <- suppressWarnings(as.data.frame(readr::read_csv(f, show_col_types = FALSE)))
    expect_equal(as.character(ok$SampleName), c("S1_ctl", "S1_trt"))
    expect_equal(as.character(ok$Group), c("ctl", "trt"))
})

test_that("a well-formed sample sheet is read identically either way", {
    f <- tempfile(fileext = ".csv"); on.exit(unlink(f), add = TRUE)
    writeLines(c("SampleName,Group", "S1,ctl", "S2,trt"), f)
    base_r <- read.csv(f, stringsAsFactors = FALSE)
    skip_if_not_installed("readr")
    rdr <- as.data.frame(readr::read_csv(f, show_col_types = FALSE))
    expect_equal(as.character(base_r$SampleName), as.character(rdr$SampleName))
    expect_equal(as.character(base_r$Group), as.character(rdr$Group))
})

test_that("the template defines read_samplesheet and no longer read.csv's the sample sheet", {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "proteomics",
                            "report_template_proteomics.Rmd"),
        "R/domain/proteomics/report_template_proteomics.Rmd"
    )
    f <- candidates[file.exists(candidates)][1]
    skip_if(is.na(f), "report template not found from the test working directory")

    src <- readLines(f, warn = FALSE)
    # The helper itself lives in R/core/01_io.R, not here: the template must
    # call it, and must no longer read.csv() any of the three metadata sites.
    expect_length(grep("read_samplesheet <- function", src, fixed = TRUE), 0)
    expect_gte(length(grep("read_samplesheet(", src, fixed = TRUE)), 3)
    expect_length(grep("read.csv(meta_file", src, fixed = TRUE), 0)
    expect_length(grep("read.csv(meta_file_exp", src, fixed = TRUE), 0)
    # The Technical Summary reads the same sheet and was the site this PR
    # originally missed; it must go through the helper too.
    expect_length(grep("read.csv(file.path(proj_dir_resolved, cfg$paths$raw",
                       src, fixed = TRUE), 0)
    expect_gte(length(grep("meta_ts <- read_samplesheet(", src, fixed = TRUE)), 1)
})


# --- Second failure mode: the wrong delimiter ------------------------------
#
# A tab-separated sample sheet read as CSV collapses into a single column. This
# is a different bug from the ragged-row shift above, and neither reader handles
# the other: read.delim() with the right separator still applies read.csv's
# row-name heuristic, and readr::read_csv() on a tab file still returns one
# column. read_samplesheet() detects the separator first, then reads with a
# reader that respects the header.

write_tab_samplesheet <- function() {
    f <- tempfile(fileext = ".txt")
    writeLines(c(
        "SeqNum\tSampleName\tGroup\tFileName",
        "1\tS1_ctl\tctl\tD:/runs/a.raw",
        "2\tS1_trt\ttrt\tD:/runs/b.raw"
    ), f)
    f
}

test_that("read.csv collapses a tab-separated sample sheet into one column", {
    f <- write_tab_samplesheet(); on.exit(unlink(f), add = TRUE)
    bad <- read.csv(f, stringsAsFactors = FALSE)
    expect_equal(ncol(bad), 1L)
    expect_false("SampleName" %in% names(bad))
})

test_that("the separator is detected from the header, not assumed", {
    f <- write_tab_samplesheet(); on.exit(unlink(f), add = TRUE)
    ok <- read_samplesheet(f)
    expect_equal(names(ok), c("SeqNum", "SampleName", "Group", "FileName"))
    expect_equal(as.character(ok$SampleName), c("S1_ctl", "S1_trt"))
    expect_equal(as.character(ok$Group), c("ctl", "trt"))
})

test_that("detecting the separator does not regress the ragged-CSV case", {
    # Both fixes have to hold at once: this sheet is comma-separated AND ragged.
    f <- write_ragged_samplesheet(); on.exit(unlink(f), add = TRUE)
    ok <- read_samplesheet(f)
    expect_equal(as.character(ok$SampleName), c("S1_ctl", "S1_trt"))
    expect_equal(as.character(ok$Group), c("ctl", "trt"))
})

test_that("a missing or empty path yields NULL rather than an error", {
    expect_null(read_samplesheet(file.path(tempdir(), "does-not-exist.csv")))
    expect_null(read_samplesheet(""))
    expect_null(read_samplesheet(NULL))
})
