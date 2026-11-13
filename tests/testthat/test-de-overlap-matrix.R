# Shared / unique DE counts between contrasts.

sets <- list(a = c("g1", "g2", "g3"), b = c("g2", "g4"), c = c("g2", NA, "g2"))

test_that("off-diagonal counts shared features, diagonal unique ones", {
    m <- de_overlap_matrix(sets)
    expect_equal(dimnames(m), list(c("a", "b", "c"), c("a", "b", "c")))
    expect_equal(unname(diag(m)), c(2L, 1L, 0L))
    expect_equal(m["a", "b"], 1L)
    expect_equal(m["b", "c"], 1L)
    expect_true(isSymmetric(m))
})

test_that("NA and duplicates are ignored", {
    expect_equal(de_overlap_matrix(sets)["c", "a"], 1L)
})

test_that("empty sets give zeros", {
    m <- de_overlap_matrix(list(a = character(0), b = "g1"))
    expect_equal(unname(m), matrix(c(0L, 0L, 0L, 1L), 2))
})

test_that("plot labels carry the per-contrast totals and one tile per cell", {
    p <- plot_de_overlap_matrix(sets, title = "x", feature_label = "genes")
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), 9L)
    expect_setequal(levels(p$data$col), c("a (n = 3)", "b (n = 2)", "c (n = 1)"))
    d <- p$data[p$data$is_diag, ]
    expect_equal(d$value[match(c("a (n = 3)", "b (n = 2)", "c (n = 1)"), d$col)], c(2L, 1L, 0L))
})
