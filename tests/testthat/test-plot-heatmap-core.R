# Which rows plot_heatmap_core() keeps when max_rows caps the heatmap.
#
# The RNA-seq and proteomics reports caption the QC expression heatmap as based
# on the top variable features, so the cap must keep the most variable rows,
# not a random subset. All fixtures are synthetic.

heatmap_rows_kept <- function(m, ...) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    ph <- suppressMessages(plot_heatmap_core(m, cluster_rows = TRUE, ...))
    ph$tree_row$labels
}

test_that("max_rows keeps the most variable rows", {
    skip_if_not_installed("pheatmap")
    m <- rbind(
        flat1 = c(5, 5.1, 5, 5.1),
        loud1 = c(1, 9, 1, 9),
        flat2 = c(3, 3.1, 3, 3.1),
        loud2 = c(2, 8, 2, 8),
        mid   = c(4, 6, 4, 6)
    )
    colnames(m) <- paste0("S", 1:4)

    expect_setequal(heatmap_rows_kept(m, max_rows = 2), c("loud1", "loud2"))
    expect_setequal(heatmap_rows_kept(m, max_rows = 3), c("loud1", "loud2", "mid"))
})

test_that("max_rows keeps the input order of the rows it keeps", {
    skip_if_not_installed("pheatmap")
    m <- rbind(
        b_loud = c(0, 10, 0, 10),
        a_flat = c(1, 1.1, 1, 1.1),
        c_loud = c(1, 9, 1, 9)
    )
    colnames(m) <- paste0("S", 1:4)

    expect_identical(heatmap_rows_kept(m, max_rows = 2), c("b_loud", "c_loud"))
})

test_that("rows under the cap are all kept", {
    skip_if_not_installed("pheatmap")
    m <- rbind(r1 = c(1, 2, 3, 4), r2 = c(4, 3, 2, 1), r3 = c(1, 3, 2, 4))
    colnames(m) <- paste0("S", 1:4)

    expect_setequal(heatmap_rows_kept(m, max_rows = 10), c("r1", "r2", "r3"))
    expect_setequal(heatmap_rows_kept(m), c("r1", "r2", "r3"))
})
