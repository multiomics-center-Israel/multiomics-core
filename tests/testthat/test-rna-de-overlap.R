# The RNA-seq report's DE overlap section: one DE rule for the All / Up / Down
# tabs, and a Venn diagram for up to three contrasts, an UpSet plot beyond.
# All fixtures are synthetic.

de_fixture <- function() {
    data.frame(
        GeneName = paste0("G", 1:5),
        padj.A_vs_B = c(0.01, 0.01, 0.20, NA, 0.01),
        linearFC.A_vs_B = c(2, -2, 3, 3, 1.1),
        padj.A_vs_B_long = c(0.01, 0.50, 0.01, 0.01, 0.01),
        linearFC.A_vs_B_long = c(-3, 2, 2, -2, 2),
        check.names = FALSE
    )
}

test_that("collect_de_gene_sets applies padj and fold-change cutoffs by direction", {
    d <- de_fixture()
    cn <- c("A_vs_B", "A_vs_B_long")
    lfc <- log2(1.5)

    all <- collect_de_gene_sets(d, cn, padj_cut = 0.05, lfc_cut = lfc, direction = "all")
    expect_identical(names(all), c("A vs B", "A vs B long"))
    expect_identical(all[["A vs B"]], c("G1", "G2"))
    expect_identical(all[["A vs B long"]], c("G1", "G3", "G4", "G5"))

    up <- collect_de_gene_sets(d, cn, padj_cut = 0.05, lfc_cut = lfc, direction = "up")
    expect_identical(up[["A vs B"]], "G1")
    down <- collect_de_gene_sets(d, cn, padj_cut = 0.05, lfc_cut = lfc, direction = "down")
    expect_identical(down[["A vs B long"]], c("G1", "G4"))
})

test_that("collect_de_gene_sets matches whole contrast names, not prefixes", {
    d <- de_fixture()[, c("GeneName", "padj.A_vs_B_long", "linearFC.A_vs_B_long")]
    expect_length(collect_de_gene_sets(d, "A_vs_B", 0.05, log2(1.5)), 0)
})

draw_quietly <- function(...) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    out <- NULL
    utils::capture.output(out <- draw_de_overlap(...))
    out
}

test_that("draw_de_overlap uses an UpSet plot above three sets", {
    skip_if_not_installed("ComplexHeatmap")
    sets <- list(a = c("x", "y"), b = c("y", "z"), c = "z", d = character(0))
    expect_identical(draw_quietly(sets, title = "t"), "upset")
})

test_that("draw_de_overlap uses a Venn diagram for up to three sets", {
    skip_if_not_installed("ggVennDiagram")
    sets <- list(a = c("x", "y"), b = c("y", "z"))
    expect_identical(draw_quietly(sets, title = "t"), "venn")
})

test_that("draw_de_overlap draws nothing when no gene is DE", {
    expect_identical(draw_quietly(list(a = character(0), b = character(0)), title = "t"), "none")
})
