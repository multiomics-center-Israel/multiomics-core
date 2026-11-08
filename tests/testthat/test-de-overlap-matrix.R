# Shared / unique DE counts between contrasts.

sets <- list(a = c("g1", "g2", "g3"), b = c("g2", "g4"), c = c("g2", NA, "g2"))

test_that("off-diagonal counts shared features, diagonal the contrast totals", {
    m <- de_overlap_matrix(sets)
    expect_equal(dimnames(m), list(c("a", "b", "c"), c("a", "b", "c")))
    expect_equal(unname(diag(m)), c(3L, 2L, 1L))
    expect_equal(m["a", "b"], 1L)
    expect_equal(m["b", "c"], 1L)
    expect_true(isSymmetric(m))
})

test_that("each row peaks on the diagonal", {
    m <- de_overlap_matrix(sets)
    expect_true(all(diag(m) == apply(m, 1, max)))
})

test_that("NA and duplicates are ignored", {
    expect_equal(de_overlap_matrix(sets)["c", "a"], 1L)
    expect_equal(de_overlap_matrix(sets)["c", "c"], 1L)
})

test_that("empty sets give zeros", {
    m <- de_overlap_matrix(list(a = character(0), b = "g1"))
    expect_equal(unname(m), matrix(c(0L, 0L, 0L, 1L), 2))
})

test_that("unique counts are features DE in one contrast only", {
    expect_equal(de_unique_counts(sets), c(a = 2L, b = 1L, c = 0L))
    expect_equal(de_unique_counts(list(a = character(0), b = "g1")), c(a = 0L, b = 1L))
})

test_that("plot labels carry totals and unique counts, diagonal shows totals", {
    p <- plot_de_overlap_matrix(sets, title = "x", feature_label = "genes")
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), 9L)
    labs <- c("a (n = 3, unique = 2)", "b (n = 2, unique = 1)", "c (n = 1, unique = 0)")
    expect_setequal(levels(p$data$col), labs)
    d <- p$data[p$data$is_diag, ]
    expect_equal(d$value[match(labs, d$col)], c(3L, 2L, 1L))
})
