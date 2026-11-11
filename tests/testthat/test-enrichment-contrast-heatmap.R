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

# --- build_enrichment_contrast_matrix(): one collection, every contrast -------

test_that("GSEA cells hold NES, untested cells stay NA, rows are sets significant anywhere", {
    res <- data.frame(
        pathway = c("p1", "p1", "p2", "p2", "p3"), pathway_name = c("P one", "P one", "P two", "P two", "P three"),
        contrast = c("a", "b", "a", "b", "a"), analysis = "GSEA",
        NES = c(2.1, -1.4, 1.2, 0.3, 0.9), padj = c(0.001, 0.2, 0.04, 0.9, 0.5),
        direction = NA_character_, stringsAsFactors = FALSE)
    b <- build_enrichment_contrast_matrix(res, c("a", "b", "c"), padj_cut = 0.05, n_max = 40)
    expect_identical(b$kind, "NES")
    expect_equal(rownames(b$value), c("p1", "p2"))           # p3 never < 0.05
    expect_equal(colnames(b$value), c("a", "b", "c"))
    expect_equal(b$value["p1", "b"], -1.4)                    # tested, not significant: shown
    expect_true(all(is.na(b$value[, "c"])))                   # contrast c not tested: blank
    expect_equal(c(b$n_sig, b$n_shown, b$n_left_out), c(2L, 2L, 0L))
    expect_equal(b$labels, c("P one", "P two"))
})

test_that("ORA up and down combine into a signed -log10(padj)", {
    res <- data.frame(
        pathway = c("p1", "p1", "p2"), pathway_name = c("P one", "P one", "P two"),
        contrast = c("a", "a", "b"), analysis = c("ORA up", "ORA down", "ORA down"),
        NES = NA_real_, padj = c(0.01, 0.2, 0.001), direction = c("up", "down", "down"),
        stringsAsFactors = FALSE)
    b <- build_enrichment_contrast_matrix(res, c("a", "b"), 0.05)
    expect_identical(b$kind, "signed_log10")
    expect_equal(b$value["p1", "a"], 2)      # up list has the smaller padj: positive
    expect_equal(b$value["p2", "b"], -3)     # down list: negative
    expect_true(is.na(b$value["p1", "b"]))
    expect_equal(rownames(b$value), c("p2", "p1"))  # ordered by smallest padj
})

test_that("without a direction column ORA is unsigned -log10(padj)", {
    res <- data.frame(pathway = "p1", pathway_name = "P", contrast = "a",
                      analysis = "ORA", NES = NA_real_, padj = 0.001,
                      direction = NA_character_, stringsAsFactors = FALSE)
    b <- build_enrichment_contrast_matrix(res, "a", 0.05)
    expect_identical(b$kind, "log10")
    expect_equal(b$value[1, 1], 3)
})

test_that("rows are capped at n_max and the rest counted as left out", {
    res <- data.frame(pathway = paste0("p", 1:5), pathway_name = paste0("p", 1:5),
                      contrast = "a", analysis = "GSEA", NES = 1:5,
                      padj = c(0.04, 0.01, 0.03, 0.02, 0.001),
                      direction = NA_character_, stringsAsFactors = FALSE)
    b <- build_enrichment_contrast_matrix(res, "a", 0.05, n_max = 3)
    expect_equal(rownames(b$value), c("p5", "p2", "p4"))
    expect_equal(c(b$n_sig, b$n_shown, b$n_left_out), c(5L, 3L, 2L))
    expect_null(build_enrichment_contrast_matrix(res, "a", padj_cut = 0.0001))
})

test_that("collect_enrichment_results reads cameraPR tables with FDR and Direction", {
    dir <- withr::local_tempdir()
    utils::write.csv(data.frame(pathway = c("p1", "p2"), NGenes = c(10, 12),
                                Direction = c("Up", "Down"), PValue = c(0.001, 0.01),
                                FDR = c(0.002, 0.02), contrast = "c1", database = "KEGG",
                                method = "camerapr", stringsAsFactors = FALSE),
                     file.path(dir, "pathway_c1_KEGG_camerapr.csv"), row.names = FALSE)
    res <- collect_enrichment_results(dir)
    expect_equal(unique(res$analysis), "cameraPR")
    expect_equal(res$padj, c(0.002, 0.02))
    expect_equal(res$direction, c("up", "down"))
    b <- build_enrichment_contrast_matrix(res, "c1", 0.05)
    expect_equal(unname(b$value[, 1]), c(-log10(0.002), log10(0.02)))
})

test_that("the matrix heatmap leaves untested cells blank and dots significant ones", {
    res <- data.frame(pathway = c("p1", "p1"), pathway_name = "P one", contrast = c("a", "b"),
                      analysis = "GSEA", NES = c(2, -1), padj = c(0.01, 0.3),
                      direction = NA_character_, stringsAsFactors = FALSE)
    b <- build_enrichment_contrast_matrix(res, c("a", "b", "c"), 0.05)
    p <- plot_enrichment_matrix_heatmap(b, title = "t")
    expect_s3_class(p, "ggplot")
    expect_equal(nrow(p$data), 3L)
    expect_equal(sum(p$data$sig), 1L)
    expect_true(is.na(p$data$value[p$data$contrast == "c"]))
})
