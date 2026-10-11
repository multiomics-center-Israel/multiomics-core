# pathway_dotplot(): GSEA dotplots carry red dashed NES reference lines, ORA
# dotplots do not, and nes_lines = NULL turns them off. Synthetic data.

gsea_top <- function() {
    data.frame(pathway_label = c("P1", "P2", "P3"), NES = c(-2.5, 1.2, 2.8),
               padj = c(0.01, 0.2, 0.001), size = c(20, 35, 50))
}

vline_layers <- function(p) {
    Filter(function(l) inherits(l$geom, "GeomVline"), p$layers)
}

test_that("GSEA dotplots mark NES -2 and 2 with red dashed lines", {
    p <- pathway_dotplot(gsea_top(), "A_vs_B", "GO_BP_fgsea")
    v <- vline_layers(p)
    expect_length(v, 1)
    expect_setequal(v[[1]]$data$xintercept, c(-2, 2))
    expect_identical(v[[1]]$aes_params$colour, "red")
    expect_identical(v[[1]]$aes_params$linetype, "dashed")
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("nes_lines = NULL draws no lines, and ORA dotplots never get them", {
    expect_length(vline_layers(pathway_dotplot(gsea_top(), "A_vs_B", "x", nes_lines = NULL)), 0)
    ora <- data.frame(pathway_label = c("P1", "P2", "P3"), pvalue = c(1e-4, 1e-3, 0.01),
                      padj = c(1e-3, 0.01, 0.05), overlap = c(5, 8, 3))
    p <- pathway_dotplot(ora, "A_vs_B", "GO_BP_ora_up")
    expect_length(vline_layers(p), 0)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})
