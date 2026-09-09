# tests/testthat/test-volcano-labels.R
#
# Tests for the top-N labelling on the static volcano (R/core/06_plots.R).
# Labels are chosen from features passing BOTH cutoffs, ranked per direction by
# |log2FC| * -log10(p). Ranking on p alone crowds labels onto whatever is most
# precisely measured; ranking on fold change alone favours features with little
# evidence. These pin the selection rule, which is the part that decides what a
# reader's eye lands on.

mk_df <- function(logfc, neglog10p, direction, genes = NULL) {
    n <- length(logfc)
    data.frame(
        .logFC = logfc, .neglog10p = neglog10p,
        .direction = factor(direction, levels = c("NS", "Down", "Up")),
        Genes = genes %||% paste0("G", seq_len(n)),
        stringsAsFactors = FALSE
    )
}
`%||%` <- function(a, b) if (is.null(a)) b else a

test_that("labelling is off by default and for n_label <= 0", {
    df <- mk_df(c(2, -2), c(5, 5), c("Up", "Down"))
    expect_null(.volcano_label_data(df, 0))
    expect_null(.volcano_label_data(df, NULL))
    expect_null(.volcano_label_data(df, NA_integer_))
})

test_that("non-significant features are never labelled", {
    df <- mk_df(c(9, -9), c(9, 9), c("NS", "NS"))
    expect_null(.volcano_label_data(df, 5))
})

test_that("at most n_label are taken from each direction", {
    df <- mk_df(c(rep(2, 20), rep(-2, 20)), rep(5, 40),
                c(rep("Up", 20), rep("Down", 20)))
    out <- .volcano_label_data(df, 3)
    expect_equal(nrow(out), 6)
    expect_equal(sum(out$.direction == "Up"), 3)
    expect_equal(sum(out$.direction == "Down"), 3)
})

test_that("ranking uses |log2FC| * -log10(p), not either alone", {
    # A: most significant but a small effect. B: largest effect, weak evidence.
    # C: moderate on both, and the highest product. Only C should be picked.
    df <- mk_df(logfc     = c(0.6,  8.0,  3.0),
                neglog10p = c(20.0, 1.5,  9.0),
                direction = rep("Up", 3),
                genes     = c("A", "B", "C"))
    expect_equal(.volcano_label_data(df, 1)$Genes, "C")
})

test_that("both directions are labelled even when one dominates", {
    df <- mk_df(c(rep(5, 10), -1.2), c(rep(10, 10), 2), c(rep("Up", 10), "Down"))
    out <- .volcano_label_data(df, 15)
    expect_true("Down" %in% as.character(out$.direction))
    expect_equal(sum(out$.direction == "Down"), 1)
})

test_that("fewer features than n_label is fine", {
    df <- mk_df(c(2, -2), c(5, 5), c("Up", "Down"))
    expect_equal(nrow(.volcano_label_data(df, 15)), 2)
})

test_that("only one direction present still returns labels", {
    df <- mk_df(c(2, 3), c(5, 6), c("Up", "Up"))
    out <- .volcano_label_data(df, 5)
    expect_equal(nrow(out), 2)
    expect_true(all(out$.direction == "Up"))
})

test_that("protein groups are trimmed to the first symbol", {
    df <- mk_df(c(2, -2), c(5, 5), c("Up", "Down"),
                genes = c("TSPY1;TSPY10;TSPY2", "ACTB"))
    out <- .volcano_label_data(df, 5)
    expect_true("TSPY1" %in% out$.label)
    expect_false(any(grepl(";", out$.label)))
})

test_that("features with no usable label are skipped, not labelled blank", {
    df <- mk_df(c(2, 3), c(5, 6), c("Up", "Up"), genes = c("", "REAL"))
    out <- .volcano_label_data(df, 5)
    expect_equal(nrow(out), 1)
    expect_equal(out$.label, "REAL")
})

test_that("FeatureID is used when no gene column exists", {
    df <- mk_df(c(2, -2), c(5, 5), c("Up", "Down"))
    df$Genes <- NULL
    df$FeatureID <- c("P1", "P2")
    expect_setequal(.volcano_label_data(df, 5)$.label, c("P1", "P2"))
})
