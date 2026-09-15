# tests/testthat/test-read-gmt-multi.R
#
# Unit tests for read_gmt() accepting multiple GMT paths (R/core/09_enrichment.R).
# A config may point one omic's `gmt_file` at several GMTs (e.g. GO + KEGG); the
# reader must merge them into one collection, keep descriptions, and keep the
# first occurrence when a set name repeats across files.

write_gmt_lines <- function(lines) {
    path <- tempfile(fileext = ".gmt")
    writeLines(lines, path)
    path
}

test_that("read_gmt reads a single file (scalar path) unchanged", {
    gmt <- write_gmt_lines(c(
        paste(c("A", "Set A", paste0("g", 1:3)), collapse = "\t"),
        paste(c("B", "Set B", paste0("g", 4:6)), collapse = "\t")
    ))
    on.exit(unlink(gmt), add = TRUE)

    gs <- read_gmt(gmt)

    expect_setequal(names(gs), c("A", "B"))
    expect_equal(gs$A, paste0("g", 1:3))
    expect_equal(attr(gs, "descriptions")[["B"]], "Set B")
})

test_that("read_gmt merges two files into one collection with combined descriptions", {
    go   <- write_gmt_lines(paste(c("GO:1", "go term", paste0("g", 1:5)), collapse = "\t"))
    kegg <- write_gmt_lines(paste(c("ehi001", "kegg term", paste0("g", 3:8)), collapse = "\t"))
    on.exit(unlink(c(go, kegg)), add = TRUE)

    gs <- read_gmt(c(go, kegg))

    expect_setequal(names(gs), c("GO:1", "ehi001"))
    expect_equal(gs[["GO:1"]], paste0("g", 1:5))
    expect_equal(gs[["ehi001"]], paste0("g", 3:8))
    descr <- attr(gs, "descriptions")
    expect_equal(descr[["GO:1"]], "go term")
    expect_equal(descr[["ehi001"]], "kegg term")
})

test_that("read_gmt keeps the first occurrence when a set name repeats across files", {
    f1 <- write_gmt_lines(paste(c("SHARED", "first", paste0("g", 1:3)), collapse = "\t"))
    f2 <- write_gmt_lines(paste(c("SHARED", "second", paste0("h", 1:9)), collapse = "\t"))
    on.exit(unlink(c(f1, f2)), add = TRUE)

    gs <- read_gmt(c(f1, f2))

    expect_equal(names(gs), "SHARED")
    expect_equal(gs$SHARED, paste0("g", 1:3))                  # first file wins
    expect_equal(attr(gs, "descriptions")[["SHARED"]], "first")
})

test_that("read_gmt accepts a list (YAML list-form) as well as a vector", {
    a <- write_gmt_lines(paste(c("A", "a", "g1", "g2"), collapse = "\t"))
    b <- write_gmt_lines(paste(c("B", "b", "g3", "g4"), collapse = "\t"))
    on.exit(unlink(c(a, b)), add = TRUE)

    gs <- read_gmt(list(a, b))

    expect_setequal(names(gs), c("A", "B"))
})

test_that("read_gmt still errors on a missing single file", {
    expect_error(read_gmt(tempfile(fileext = ".gmt")), "GMT file not found")
})
