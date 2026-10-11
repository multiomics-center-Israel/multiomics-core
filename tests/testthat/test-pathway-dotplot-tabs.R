# list_pathway_dotplots(): how the report groups pathway dotplots into
# contrast tabs and orders the plot sub-tabs. Synthetic file names.

test_that("dotplots are grouped by contrast and ordered within each", {
    d <- withr::local_tempdir()
    files <- c("pathway_B_vs_A_KEGG_ora_down.png", "pathway_B_vs_A_KEGG_fgsea.png",
               "pathway_B_vs_A_GO_BP_ora_up.png", "pathway_B_vs_A_GO_BP_fgsea.png",
               "pathway_B_long_vs_A_GO_BP_fgsea.png", "pathway_mystery_KEGG_fgsea.png")
    for (f in files) file.create(file.path(d, f))

    out <- list_pathway_dotplots(d, c("B_long_vs_A", "B_vs_A"))

    expect_identical(out$contrast, c("B_long_vs_A", rep("B_vs_A", 4), NA))
    expect_identical(out$analysis[2:5], c("GO_BP_fgsea", "GO_BP_ora_up", "KEGG_fgsea", "KEGG_ora_down"))
    expect_identical(out$label[2], "GO BP fgsea")
    expect_identical(basename(out$path[6]), "pathway_mystery_KEGG_fgsea.png")
})

test_that("a contrast whose name prefixes another is not matched to its files", {
    d <- withr::local_tempdir()
    file.create(file.path(d, "pathway_X_4h_vs_base_heme_KEGG_fgsea.png"))
    out <- list_pathway_dotplots(d, c("X_4h_vs_base", "X_4h_vs_base_heme"))
    expect_identical(out$contrast, "X_4h_vs_base_heme")
    expect_identical(out$analysis, "KEGG_fgsea")
})

test_that("no dotplots gives zero rows", {
    expect_equal(nrow(list_pathway_dotplots(withr::local_tempdir(), "A_vs_B")), 0)
})
