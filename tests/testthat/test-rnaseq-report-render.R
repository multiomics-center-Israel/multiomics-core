# tests/testthat/test-rnaseq-report-render.R
#
# render_rnaseq_report() renders with quiet = TRUE, which drops knitr's
# "Quitting from lines ..." note. A failure surfaced only as "missing value
# where TRUE/FALSE needed" with no hint of where in a 3,000-line template.
# The renderer now records the failing chunk's label and maps it to lines with
# find_rmd_chunk_lines().

write_chunk_rmd <- function() {
    f <- tempfile(fileext = ".Rmd")
    writeLines(c(
        "---",
        "title: test",
        "---",
        "```{r setup, include=FALSE}",
        "x <- 1",
        "```",
        "",
        "```{r}",
        "y <- 2",
        "```",
        "",
        "```{r technical-summary, results='asis', echo=FALSE}",
        "if (NA) cat('never')",
        "```"
    ), f)
    f
}

test_that("find_rmd_chunk_lines returns the fence lines of a labelled chunk", {
    rmd <- write_chunk_rmd()
    expect_equal(find_rmd_chunk_lines(rmd, "setup"), c(4L, 6L))
    expect_equal(find_rmd_chunk_lines(rmd, "technical-summary"), c(12L, 14L))
})

test_that("find_rmd_chunk_lines returns NULL for labels it cannot place", {
    rmd <- write_chunk_rmd()
    expect_null(find_rmd_chunk_lines(rmd, "unnamed-chunk-1"))
    expect_null(find_rmd_chunk_lines(rmd, "no-such-chunk"))
    # A prefix of a real label must not match it.
    expect_null(find_rmd_chunk_lines(rmd, "technical"))
})

test_that("the renderer renders the template once, inside its error handler", {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "rnaseq", "08_report.R"),
        "R/domain/rnaseq/08_report.R"
    )
    src <- readLines(candidates[file.exists(candidates)][1], warn = FALSE)
    body <- src[seq_len(grep("^find_rmd_chunk_lines <- ", src) - 1L)]
    expect_length(grep("rmarkdown::render(", body, fixed = TRUE), 1L)
})
