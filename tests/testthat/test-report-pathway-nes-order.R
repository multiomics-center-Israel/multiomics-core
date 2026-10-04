# tests/testthat/test-report-pathway-nes-order.R
#
# Tests for the pathway-volcano dropdown ordering in the proteomics report
# template. The list was sorted alphabetically, which buried the most strongly
# shifted gene sets among hundreds of entries. It now sorts by |NES| descending.
# NES is direction-signed, so the magnitude is what ranks a set's shift; a
# strongly down set must rank above a weakly up one.

order_pathways_by_abs_nes <- function(pathways, pathway_info) {
    abs_nes <- vapply(pathways, function(pw) {
        nes <- pathway_info[[pw]]$NES
        if (is.null(nes) || !is.finite(nes)) NA_real_ else abs(nes)
    }, numeric(1))
    pathways[order(-abs_nes, pathways, na.last = TRUE)]
}

test_that("pathways sort by |NES| descending, not alphabetically", {
    info <- list(Alpha = list(NES = 0.5), Beta = list(NES = -2.4), Gamma = list(NES = 1.8))
    expect_equal(order_pathways_by_abs_nes(c("Alpha", "Beta", "Gamma"), info),
                 c("Beta", "Gamma", "Alpha"))
})

test_that("a strongly down set outranks a weakly up one", {
    info <- list(Down = list(NES = -3.1), Up = list(NES = 1.2))
    expect_equal(order_pathways_by_abs_nes(c("Up", "Down"), info)[1], "Down")
})

test_that("sets without NES sort last but are never dropped", {
    info <- list(WithNES = list(NES = 1.1))
    out <- order_pathways_by_abs_nes(c("NoNES", "WithNES"), info)
    expect_equal(out, c("WithNES", "NoNES"))
    expect_length(out, 2)
})

test_that("non-finite NES is treated as missing, not as infinity", {
    info <- list(Inf1 = list(NES = Inf), Real = list(NES = 2))
    expect_equal(order_pathways_by_abs_nes(c("Inf1", "Real"), info), c("Real", "Inf1"))
})

test_that("ties break alphabetically for a stable order", {
    info <- list(B = list(NES = 2), A = list(NES = -2))
    expect_equal(order_pathways_by_abs_nes(c("B", "A"), info), c("A", "B"))
})

test_that("every input pathway survives the ordering", {
    info <- list(A = list(NES = 1))
    out <- order_pathways_by_abs_nes(c("A", "B", "C"), info)
    expect_setequal(out, c("A", "B", "C"))
})
