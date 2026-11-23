# Which pathways get a GSEA/pair pathway map, and why.
#
# ORA on the DE list misses pathways that rank-based GSEA calls enriched and
# pathways holding a changed enzyme from the enzyme-metabolite table. These
# selectors decide which of those get a map. All fixtures synthetic; no
# pathview call, no network.

pairs_fixture <- function() {
    data.frame(
        contrast = c("A_vs_B", "A_vs_B", "A_vs_B", "A_vs_B"),
        gene_symbol = c("Enz1", "Enz1", "Enz2", "Enz3"),
        enzyme_hit = c(TRUE, TRUE, TRUE, FALSE),
        enzyme_padj = c(0.01, 0.01, 0.04, 0.50),
        pathway_id = c("00310;01100", "00310;01100", "00230", "00010"),
        stringsAsFactors = FALSE
    )
}

gsea_fixture <- function() {
    list(prot = data.frame(
        pathway = c("xyz03008", "xyz00500", "xyz00900", "GO:0000001", "xyz01240"),
        padj = c(0.001, 0.30, 0.02, 0.001, 0.01),
        method = c("fgsea", "fgsea", "ora", "fgsea", "fgsea"),
        contrast = "A_vs_B",
        stringsAsFactors = FALSE
    ))
}


# ---- enzyme-metabolite pairs ------------------------------------------------

test_that("a changed enzyme's pathways are selected, global maps are not", {
    hits <- pair_pathways_by_contrast(pairs_fixture())

    expect_length(hits, 1)
    h <- hits[[1]]
    expect_equal(h$label, "A_vs_B")
    # 01100 is a global map pathview cannot overlay; 00010 belongs to an enzyme
    # that is not a hit.
    expect_equal(h$pathways, c("00310", "00230"))
    expect_equal(unname(h$genes["00310"]), "Enz1")
})

test_that("a pair table with nothing usable selects nothing", {
    expect_equal(pair_pathways_by_contrast(NULL), list())
    p <- pairs_fixture()
    p$enzyme_hit <- FALSE
    expect_equal(pair_pathways_by_contrast(p), list())
    expect_equal(pair_pathways_by_contrast(p[, c("contrast", "gene_symbol")]), list())
})


# ---- combined selection -----------------------------------------------------

test_that("GSEA rows count only when rank-based, KEGG, significant and not global", {
    sel <- select_gsea_pair_pathways(gsea_fixture(), NULL, kegg_org = "xyz")

    # 00500 fails padj, 00900 is ORA, the GO term is not KEGG, 01240 is global.
    expect_equal(sel$pathway, "03008")
    expect_equal(sel$source, "GSEA")
    expect_equal(sel$gsea_padj, 0.001)
})

test_that("a pathway reached both ways is listed once with both reasons", {
    g <- gsea_fixture()
    g$prot$pathway[1] <- "xyz00310"
    sel <- select_gsea_pair_pathways(g, pairs_fixture(), kegg_org = "xyz")

    expect_equal(sum(sel$pathway == "00310"), 1)
    expect_equal(sel$source[sel$pathway == "00310"],
                 "GSEA + enzyme-metabolite pair")
    expect_equal(sel$source[sel$pathway == "00230"], "enzyme-metabolite pair")
})

test_that("top_n caps GSEA pathways but not pair pathways", {
    g <- gsea_fixture()
    g$prot$method <- "fgsea"
    g$prot$padj <- 0.001
    sel <- select_gsea_pair_pathways(g, pairs_fixture(), kegg_org = "xyz",
                                     top_n = 1)

    expect_equal(sum(grepl("GSEA", sel$source)), 1)
    expect_true(all(c("00310", "00230") %in% sel$pathway))
})

test_that("no input yields an empty, well-formed frame", {
    sel <- select_gsea_pair_pathways(list(), NULL, kegg_org = "xyz")
    expect_equal(nrow(sel), 0)
    expect_true(all(c("pathway", "source", "contrast") %in% names(sel)))
})
