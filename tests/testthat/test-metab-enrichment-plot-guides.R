# tests/testthat/test-metab-enrichment-plot-guides.R
#
# Reference lines on the metabolomics enrichment plots: a dashed red line at
# FDR = 0.05 on -log10(FDR) axes, dashed red lines at NES = -2 and 2 on NES
# axes, and a symmetric NES axis of -3..3 that widens to max |NES| + 0.5.
# All data here are synthetic.

# Computed data of every layer, read through layer_data() so the checks do not
# depend on the internals of the installed ggplot2 version.
plot_layer_frames <- function(p) {
    out <- list()
    i <- 1L
    repeat {
        d <- tryCatch(ggplot2::layer_data(p, i), error = function(e) NULL)
        if (is.null(d)) break
        out[[i]] <- d
        i <- i + 1L
    }
    out
}

red_line_positions <- function(p, column) {
    frames <- plot_layer_frames(p)
    vals <- lapply(frames, function(d) {
        if (column %in% names(d) && "colour" %in% names(d) &&
            isTRUE(all(d$colour == "red"))) d[[column]] else NULL
    })
    sort(unlist(vals))
}

synthetic_gsea <- function(nes) {
    data.frame(
        pathway = paste0("map0000", seq_along(nes), " - synthetic set ", seq_along(nes)),
        NES     = nes,
        FDR     = seq(0.01, 0.5, length.out = length(nes)),
        stringsAsFactors = FALSE
    )
}

test_that("nes_axis_limit keeps -3..3 and widens past max |NES| + 0.5", {
    expect_equal(nes_axis_limit(c(1.2, -1.8)), 3)
    expect_equal(nes_axis_limit(c(2.5, -0.3)), 3)
    expect_equal(nes_axis_limit(c(2.7, 1)), 3.2)
    expect_equal(nes_axis_limit(c(-3.4, 1)), 3.9)
    expect_equal(nes_axis_limit(c(NA, -1)), 3)
    expect_equal(nes_axis_limit(numeric(0)), 3)
})

test_that("NES barplot has red lines at -2 and 2 and a symmetric axis", {
    p <- plot_gsea_nes_barplot(synthetic_gsea(c(1.2, -1.8, 0.4)))
    expect_equal(red_line_positions(p, "yintercept"), c(-2, 2))
    expect_equal(ggplot2::layer_scales(p)$y$get_limits(), c(-3, 3))

    p_wide <- plot_gsea_nes_barplot(synthetic_gsea(c(3.1, -0.5)))
    expect_equal(ggplot2::layer_scales(p_wide)$y$get_limits(), c(-3.6, 3.6))
})

test_that("GSEA lollipop has red lines at -2 and 2 and a symmetric axis", {
    p <- plot_gsea_lollipop(synthetic_gsea(c(-2.9, 0.8)))
    expect_equal(red_line_positions(p, "xintercept"), c(-2, 2))
    expect_equal(ggplot2::layer_scales(p)$x$get_limits(), c(-3.4, 3.4))
})

test_that("-log10(FDR) plots mark FDR = 0.05 with a red line", {
    enrich <- data.frame(
        pathway = c("map00001 - synthetic a", "map00002 - synthetic b"),
        FDR     = c(0.01, 0.3),
        raw_p   = c(0.001, 0.1),
        hits    = c(4, 3),
        stringsAsFactors = FALSE
    )
    expect_equal(red_line_positions(plot_enrichment_barplot(enrich), "yintercept"),
                 -log10(0.05))
    expect_equal(red_line_positions(plot_qea_lollipop(enrich), "xintercept"),
                 -log10(0.05))
})
