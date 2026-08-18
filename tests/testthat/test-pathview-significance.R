# Thresholding of the multi-omics KEGG (pathview) maps: which pathways get
# drawn, and which features are allowed to colour a node.

make_ora <- function(pathway, pvalue, padj = NULL) {
    d <- data.frame(pathway = pathway, pvalue = pvalue, stringsAsFactors = FALSE)
    if (!is.null(padj)) d$padj <- padj
    d
}

test_that("select_ora_pathways_by_fdr selects on adjusted p, not raw p", {
    f <- withr::local_tempfile(fileext = ".csv")
    write.csv(make_ora(c("map00010", "map00020", "map00030"),
                       pvalue = c(0.001, 0.01, 0.30),
                       padj = c(0.02, 0.40, 0.90)),
              f, row.names = FALSE)

    expect_equal(select_ora_pathways_by_fdr(f, "^map[0-9]+$"), "map00010")
})

test_that("select_ora_pathways_by_fdr unions across layers", {
    dir <- withr::local_tempdir()
    rna <- file.path(dir, "rna_ora_up.csv")
    prot <- file.path(dir, "prot_ora_down.csv")
    write.csv(make_ora(c("map00010", "map00020"), c(0.001, 0.001),
                       padj = c(0.01, 0.90)), rna, row.names = FALSE)
    write.csv(make_ora(c("map00020", "map00030"), c(0.001, 0.001),
                       padj = c(0.90, 0.01)), prot, row.names = FALSE)

    expect_setequal(select_ora_pathways_by_fdr(c(rna, prot), "^map[0-9]+$"),
                    c("map00010", "map00030"))
})

test_that("select_ora_pathways_by_fdr keeps only ids in the requested ID space", {
    f <- withr::local_tempfile(fileext = ".csv")
    write.csv(make_ora(c("map00010", "GO:0006412", "spu00010"),
                       pvalue = rep(0.001, 3), padj = rep(0.001, 3)),
              f, row.names = FALSE)

    expect_equal(select_ora_pathways_by_fdr(f, "^map[0-9]+$"), "map00010")
})

test_that("select_ora_pathways_by_fdr announces a raw-p fallback", {
    f <- withr::local_tempfile(fileext = ".csv")
    write.csv(make_ora(c("map00010", "map00020"), pvalue = c(0.01, 0.60)),
              f, row.names = FALSE)

    expect_message(res <- select_ora_pathways_by_fdr(f, "^map[0-9]+$"),
                   "no usable adjusted p")
    expect_equal(res, "map00010")
})

test_that("select_ora_pathways_by_fdr treats an all-NA adjusted p as absent", {
    f <- withr::local_tempfile(fileext = ".csv")
    write.csv(make_ora(c("map00010", "map00020"), pvalue = c(0.01, 0.60),
                       padj = c(NA_real_, NA_real_)),
              f, row.names = FALSE)

    expect_message(res <- select_ora_pathways_by_fdr(f, "^map[0-9]+$"),
                   "no usable adjusted p")
    expect_equal(res, "map00010")
})

test_that("select_ora_pathways_by_fdr skips a table with no p-value column", {
    f <- withr::local_tempfile(fileext = ".csv")
    write.csv(data.frame(pathway = c("map00010", "map00020"), size = c(5, 7)),
              f, row.names = FALSE)

    expect_message(res <- select_ora_pathways_by_fdr(f, "^map[0-9]+$"),
                   "no p-value column")
    expect_length(res, 0)
})

test_that("filter_changed_features requires BOTH a large and a significant change", {
    # The rule is AND, not OR: a big but unsupported change is usually a noisy
    # low-abundance feature, and a confident but tiny change is not what a
    # pathway map is meant to highlight.
    de <- data.frame(
        feature_id = c("big_sig", "down_big_sig", "big_ns", "small_sig",
                       "small_ns", "borderline_sig"),
        log2fc = c(1.5, -1.5, 1.5, 0.1, 0.1, log2(1.5)),
        pvalue = c(0.01, 0.01, 0.90, 0.01, 0.90, 0.01),
        stringsAsFactors = FALSE
    )

    kept <- suppressMessages(filter_changed_features(de))
    # borderline sits exactly at the cut-off, and the fold-change rule is a
    # strict ">", so being significant does not rescue it.
    expect_setequal(kept$feature_id, c("big_sig", "down_big_sig"))
})

test_that("filter_changed_features drops features with no usable fold change", {
    # NaN is what log2() of a signed linear fold change leaves behind for every
    # down-regulated feature -- such a node must stay uncoloured, not be coloured
    # on its p-value alone.
    de <- data.frame(
        feature_id = c("nan_fc", "na_fc", "clean"),
        log2fc = c(NaN, NA_real_, 2),
        pvalue = c(0.001, 0.001, 0.001),
        stringsAsFactors = FALSE
    )

    kept <- suppressMessages(filter_changed_features(de))
    expect_equal(kept$feature_id, "clean")
})

test_that("filter_changed_features falls back to an adjusted p column", {
    # Fold changes are large enough to pass on their own, so this isolates which
    # p column the filter picked up.
    de <- data.frame(
        feature_id = c("a", "b"),
        log2fc = c(2, 2),
        adj.P.Val = c(0.01, 0.90),
        stringsAsFactors = FALSE
    )

    expect_message(kept <- filter_changed_features(de), "adj\\.P\\.Val")
    expect_equal(kept$feature_id, "a")
})

test_that("filter_changed_features handles empty and malformed tables", {
    empty <- data.frame(feature_id = character(0), log2fc = numeric(0),
                        pvalue = numeric(0), stringsAsFactors = FALSE)
    expect_equal(nrow(suppressMessages(filter_changed_features(empty))), 0)

    no_fc <- data.frame(feature_id = "a", pvalue = 0.001, stringsAsFactors = FALSE)
    expect_equal(nrow(suppressMessages(filter_changed_features(no_fc))), 0)

    expect_null(filter_changed_features(NULL))
})

test_that("pathview_significance_caption states the thresholds it was given", {
    cap <- pathview_significance_caption()
    expect_true(grepl("FDR < 0.05", cap, fixed = TRUE))
    expect_true(grepl("1.5-fold", cap, fixed = TRUE))
    expect_true(grepl("0.58", cap, fixed = TRUE))
    expect_true(grepl("raw p < 0.05", cap, fixed = TRUE))
    # An uncoloured node is ambiguous on the map, and the caption must say so
    # rather than promise a distinction pathview does not draw.
    expect_true(grepl("no measured feature", cap, fixed = TRUE))
    expect_true(grepl("does not separate the two", cap, fixed = TRUE))
})

test_that("pathview_significance_caption follows non-default thresholds", {
    cap <- pathview_significance_caption(
        list(fdr_alpha = 0.1, node_fc = 2, node_p = 0.01))
    expect_true(grepl("FDR < 0.1", cap, fixed = TRUE))
    expect_true(grepl("2-fold", cap, fixed = TRUE))
    expect_true(grepl("raw p < 0.01", cap, fixed = TRUE))
})
