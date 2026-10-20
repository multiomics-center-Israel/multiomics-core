# tests/testthat/test-validate-pooled-script.R
#
# topk_overlap() in scripts/validate_pooled_enrichment.R. The script runs on
# load, so only the function definition is evaluated here.

load_topk_overlap <- function() {
    path <- testthat::test_path("..", "..", "scripts", "validate_pooled_enrichment.R")
    exprs <- parse(path)
    def <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                      identical(e[[2]], as.name("topk_overlap")), as.list(exprs))
    env <- new.env()
    eval(def[[1]], env)
    env$topk_overlap
}

test_that("top-k overlap reports the k it compared", {
    topk_overlap <- load_topk_overlap()
    a <- c(5, 4, 3, 2, 1, NA)
    b <- c(5, 4, 1, 2, 3, 0)
    ov <- topk_overlap(a, b, k = 100)
    expect_equal(ov[["k"]], 5)          # five features finite in both, not 100
    expect_equal(ov[["top"]], 1)
    ov2 <- topk_overlap(a, b, k = 2)
    expect_equal(ov2[["k"]], 2)
    expect_equal(ov2[["top"]], 1)
})

test_that("top-k overlap with nothing finite is empty, not a perfect overlap", {
    topk_overlap <- load_topk_overlap()
    ov <- topk_overlap(c(NA, Inf), c(1, 2), k = 100)
    expect_equal(ov[["k"]], 0)
    expect_true(is.na(ov[["top"]]) && is.na(ov[["bottom"]]))
})
