# Volcano plots of several contrasts share one axis range.

volcano_cfg <- list(de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5))

tbl_a <- data.frame(logFC = c(-1, 2), P.Value = c(0.1, 1e-4), adj.P.Val = c(0.5, 0.01))
tbl_b <- data.frame(logFC = c(-3, 0.5), P.Value = c(1e-5, 0.3), adj.P.Val = c(0.001, 0.6))

test_that("limits span all contrasts, x symmetric, y from 0", {
    lims <- volcano_axis_limits(list(a = tbl_a, b = tbl_b))
    expect_equal(lims$x, c(-3, 3))
    expect_equal(lims$y, c(0, 3))
})

test_that("pval limits read the raw p-value column", {
    lims <- volcano_axis_limits(list(a = tbl_a, b = tbl_b), "pval")
    expect_equal(lims$y, c(0, 5))
})

test_that("fewer than two usable tables gives NULL", {
    expect_null(volcano_axis_limits(list(a = tbl_a)))
    expect_null(volcano_axis_limits(list(a = tbl_a, b = NULL)))
    expect_null(volcano_axis_limits(list(a = tbl_a, b = tbl_a[0, ])))
})

test_that("p of 0 is floored at 1e-300 as in plot_volcano", {
    t0 <- tbl_b
    t0$adj.P.Val[1] <- 0
    expect_equal(volcano_axis_limits(list(a = tbl_a, b = t0))$y, c(0, 300))
})

test_that("plot_volcano applies xlim and ylim", {
    p <- plot_volcano(tbl_a, volcano_cfg, xlim = c(-3, 3), ylim = c(0, 5))
    r <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]
    expect_true(r$x.range[1] <= -3 && r$x.range[2] >= 3)
    expect_true(r$y.range[2] >= 5)
})

test_that("plot_volcano without limits is unchanged in range", {
    p <- plot_volcano(tbl_a, volcano_cfg)
    r <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]
    expect_lt(r$x.range[1], -1)
    expect_lt(r$x.range[2], 3)
})
