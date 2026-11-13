# Cross-contrast enrichment heatmaps: reading, selecting, plotting.

write_res <- function(dir, contrast, db, method, rows, direction = NULL) {
    d <- data.frame(pathway = rows$pathway, pathway_name = rows$pathway,
                    padj = rows$padj, contrast = contrast, database = db,
                    method = method, stringsAsFactors = FALSE)
    if (method == "fgsea") d$NES <- rows$NES
    if (!is.null(direction)) d$direction <- direction
    suffix <- if (method == "fgsea") "fgsea" else paste0("ora_", direction)
    utils::write.csv(d, file.path(dir, sprintf("pathway_%s_%s_%s.csv", contrast, db, suffix)),
                     row.names = FALSE)
}

test_that("collect_enrichment_results reads GSEA and ORA, skips volcano data", {
    dir <- withr::local_tempdir()
    write_res(dir, "c1", "KEGG", "fgsea", data.frame(pathway = c("p1", "p2"), padj = c(0.01, 0.5), NES = c(2, -1)))
    write_res(dir, "c1", "KEGG", "ora", data.frame(pathway = "p1", padj = 0.02), direction = "up")
    utils::write.csv(data.frame(x = 1), file.path(dir, "pathway_volcano_data_c1.csv"))
    res <- collect_enrichment_results(dir)
    expect_setequal(unique(res$analysis), c("GSEA", "ORA up"))
    expect_equal(nrow(res), 3L)
    expect_true(all(is.na(res$NES[res$analysis == "ORA up"])))
})

test_that("missing directory gives zero rows", {
    expect_equal(nrow(collect_enrichment_results(file.path(tempdir(), "nope"))), 0L)
})

test_that("shared, unique and top follow the significance counts", {
    cs <- paste0("c", 1:6)
    res <- rbind(
        data.frame(pathway = "wide", contrast = cs, padj = c(rep(0.01, 4), 0.5, 0.5)),
        data.frame(pathway = "one", contrast = cs, padj = c(0.001, rep(0.9, 5))),
        data.frame(pathway = "three", contrast = cs, padj = c(0.03, 0.03, 0.03, rep(0.9, 3))),
        data.frame(pathway = "none", contrast = cs, padj = 0.9)
    )
    s <- select_enrichment_heatmap_sets(res, cs, 0.05, n = 30)
    expect_equal(s$shared, "wide")
    expect_equal(s$unique, "one")
    expect_equal(s$top, c("one", "wide", "three"))
    expect_equal(select_enrichment_heatmap_sets(res, cs, 0.05, n = 2)$top, c("one", "wide"))
})

test_that("heatmap marks significant cells and leaves untested ones NA", {
    res <- data.frame(pathway = c("p1", "p1", "p2"), pathway_name = c("P one", "P one", "P two"),
                      contrast = c("a", "b", "a"), analysis = "GSEA",
                      NES = c(2, -1.5, 1), padj = c(0.01, 0.2, 0.04))
    p <- plot_enrichment_contrast_heatmap(res, c("p1", "p2"), c("a", "b"), 0.05)
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), 4L)
    expect_equal(sum(p$data$sig), 2L)
    expect_true(is.na(p$data$value[p$data$pathway == "p2" & p$data$contrast == "b"]))
    expect_null(plot_enrichment_contrast_heatmap(res, character(), c("a", "b")))
})

test_that("ORA heatmap uses -log10(padj)", {
    res <- data.frame(pathway = "p1", pathway_name = "P", contrast = "a",
                      analysis = "ORA up", NES = NA_real_, padj = 0.001)
    p <- plot_enrichment_contrast_heatmap(res, "p1", "a", 0.05)
    expect_equal(p$data$value, 3)
})
