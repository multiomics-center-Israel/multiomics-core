# Which terms the pooled ORA bar plot shows, and which collection each came
# from.
#
# A pooled ORA over several GMTs mixes namespaces. GO contributes thousands of
# sets where KEGG contributes a few hundred, so a plain top-n by significance
# leaves a figure that looks like a GO-only analysis however much evidence the
# other collections hold. Terms are therefore drawn round-robin across
# collections -- which decides membership only; the bars stay in global evidence
# order so a reader is never told the second bar outranks the third when it does
# not.
#
# Two things this must not do. It must not decide what is KEGG on the shape of
# an identifier: that goes through the #201 identity contract, because a rule
# like "two to four letters then five digits" claims a custom set named
# abcd12345. And it must not depend on the order rows happened to be bound in:
# every sort key here is derived from the data.
#
# Selection is display only. The class exclusion has already been applied to
# these tables upstream, no p-value or adjustment is touched, and the collection
# is read off the identifier rather than written onto the table.
#
# All fixtures are synthetic.

ora_fixture <- function(id, padj = NA_real_, pvalue = NA_real_,
                        pathway = paste("Set", seq_along(id))) {
    # length.out rather than data.frame()'s own recycling, so a zero-row table
    # is expressible: recycling a length-one default against no ids is an error.
    n <- length(id)
    data.frame(
        pathway = rep(pathway, length.out = n),
        ID      = id,
        pvalue  = rep(pvalue,  length.out = n),
        padj    = rep(padj,    length.out = n),
        stringsAsFactors = FALSE
    )
}

# The bar plot tests below open a throwaway device the same way the cross-omics
# heatmap tests do -- these are about which rows are chosen, not about pixels.


# ---- classification: KEGG by contract, nothing by shape --------------------

test_that("every KEGG spelling this pipeline produces is recognised", {
    df <- ora_fixture(c("hsa04110", "map00010", "ko00020", "00030",
                        "map00040 Pentose phosphate pathway"))

    expect_identical(classify_pathway_collection(df, kegg_org = "hsa"),
                     rep("KEGG", 5))
})

test_that("an organism prefix is only KEGG for the organism of the run", {
    # The contract's answer, not a gap: hsa04110 is not this run's accession
    # when the run is on mouse, and there is no organism code at all on the
    # GMT fallback path, which is where this selector is actually used.
    df <- ora_fixture("hsa04110")

    expect_identical(classify_pathway_collection(df, kegg_org = "mmu"), "Other")
    expect_identical(classify_pathway_collection(df, kegg_org = NULL), "Other")
    # map/ko are species-neutral, so they survive a NULL organism.
    expect_identical(
        classify_pathway_collection(ora_fixture(c("map00010", "ko00020")),
                                    kegg_org = NULL),
        c("KEGG", "KEGG"))
})

test_that("a custom gene-set name shaped like an accession is not claimed", {
    # Both of these were taken by the regexes in the old stack: ^[a-z]{2,4}[0-9]{5}$
    # claimed abcd12345, and an unanchored ^GO:?[0-9]+ claimed GO12345_signalling.
    df <- ora_fixture(c("abcd12345", "GO12345_signalling", "PFAM00069x",
                        "IPR00071", "my favourite set"))

    expect_identical(classify_pathway_collection(df, kegg_org = "hsa"),
                     rep("Other", 5))
})

test_that("GO, Pfam and InterPro are distinguishable from each other and from KEGG", {
    df <- ora_fixture(c("GO:0006915", "PF00069", "IPR000719", "hsa04110",
                        "something else"))

    expect_identical(classify_pathway_collection(df, kegg_org = "hsa"),
                     c("GO", "Pfam", "InterPro", "KEGG", "Other"))
})

test_that("identity follows the contract's ladder, not the ID column alone", {
    # No ID at all: pathway_join_key() falls through to `pathway`, so a GMT
    # table that names its sets in the readable column still classifies.
    df <- data.frame(pathway = c("GO:0006915", "hsa04110"),
                     pvalue = c(0.01, 0.02), stringsAsFactors = FALSE)

    expect_identical(classify_pathway_collection(df, kegg_org = "hsa"),
                     c("GO", "KEGG"))
})

test_that("a row with no identity is Other, and no row is dropped or rewritten", {
    # Nothing in any identity column, so the contract has nothing to resolve.
    df <- ora_fixture(c(NA, "", "   ", "GO:0006915"),
                      pathway = c(NA, "", "  ", "Apoptosis"))
    before <- df

    expect_identical(classify_pathway_collection(df, kegg_org = "hsa"),
                     c("Other", "Other", "Other", "GO"))
    # Classification reads; it does not edit.
    expect_identical(df, before)
})


# ---- selection: membership by round-robin ----------------------------------

test_that("with one collection the selection is the plain top-n", {
    df <- ora_fixture(c("GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004"),
                      padj = c(0.04, 0.01, 0.03, 0.02))

    got <- select_top_ora_per_collection(df, top_n = 2)

    expect_identical(got$ID, c("GO:0000002", "GO:0000004"))
    # Which is exactly what ranking on the score alone would have given.
    expect_identical(got$ID, df$ID[order(df$padj)][1:2])
})

test_that("a collection that contributes few sets is not crowded out by a large one", {
    # Four GO terms all stronger than the single KEGG map. A plain top-n with
    # two slots shows GO twice and the figure reads as a GO-only analysis.
    df <- ora_fixture(c("GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004",
                        "hsa04110"),
                      padj = c(0.001, 0.002, 0.003, 0.004, 0.010))

    got <- select_top_ora_per_collection(df, top_n = 2, kegg_org = "hsa")

    expect_identical(got$ID, c("GO:0000001", "hsa04110"))
})

test_that("slots a short collection cannot use go to the collections that can", {
    # Five GO terms against one KEGG map, four slots. KEGG takes its one and the
    # other three go back to GO -- no quota is left unspent.
    df <- ora_fixture(c(paste0("GO:000000", 1:5), "hsa04110"),
                      padj = c(0.001, 0.002, 0.003, 0.004, 0.005, 0.010))

    got <- select_top_ora_per_collection(df, top_n = 4, kegg_org = "hsa")

    expect_equal(nrow(got), 4L)
    expect_identical(got$ID, c("GO:0000001", "GO:0000002", "GO:0000003",
                               "hsa04110"))
})

test_that("with fewer slots than collections the strongest collections are shown", {
    df <- ora_fixture(c("GO:0000001", "hsa04110", "PF00069"),
                      padj = c(0.001, 0.002, 0.003))

    got <- select_top_ora_per_collection(df, top_n = 2, kegg_org = "hsa")

    expect_identical(got$ID, c("GO:0000001", "hsa04110"))
})

test_that("round-robin decides membership, and the rows still come back strongest-first", {
    # Picked in the order GO1, KEGG1, GO2, KEGG2 -- and returned in the order
    # their scores put them, which is not the same sequence.
    df <- ora_fixture(c("GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004",
                        "hsa04110", "hsa04115"),
                      padj = c(1e-9, 1e-8, 1e-7, 1e-6, 1e-3, 1e-2))

    got <- select_top_ora_per_collection(df, top_n = 4, kegg_org = "hsa")

    expect_identical(got$ID, c("GO:0000001", "GO:0000002", "hsa04110",
                               "hsa04115"))
    expect_false(is.unsorted(got$padj))
})

test_that("a table that fits keeps every row", {
    df <- ora_fixture(c("GO:0000002", "hsa04110", "GO:0000001"),
                      padj = c(0.02, 0.03, 0.01))

    got <- select_top_ora_per_collection(df, top_n = 10, kegg_org = "hsa")

    expect_equal(nrow(got), 3L)
    expect_identical(got$ID, c("GO:0000001", "GO:0000002", "hsa04110"))
})

test_that("the selection does not depend on the order rows arrived in", {
    # Two GO terms tied on the score exactly. With arrival index as the final
    # key this flips when the table is bound the other way round; with the
    # identity as the key it cannot.
    df <- ora_fixture(c("GO:0000002", "GO:0000001", "hsa04110"),
                      padj = c(1e-5, 1e-5, 1e-3))
    reversed <- df[rev(seq_len(nrow(df))), , drop = FALSE]

    a <- select_top_ora_per_collection(df, top_n = 2, kegg_org = "hsa")
    b <- select_top_ora_per_collection(reversed, top_n = 2, kegg_org = "hsa")

    expect_identical(a$ID, c("GO:0000001", "hsa04110"))
    expect_identical(a$ID, b$ID)

    # And with a single slot, the tie inside GO resolves the same way either way.
    expect_identical(
        select_top_ora_per_collection(df, top_n = 1, kegg_org = "hsa")$ID,
        select_top_ora_per_collection(reversed, top_n = 1, kegg_org = "hsa")$ID)
})

test_that("selection leaves its input alone", {
    df <- ora_fixture(c("GO:0000001", "hsa04110"), padj = c(0.01, 0.02))
    before <- df

    select_top_ora_per_collection(df, top_n = 1, kegg_org = "hsa")

    expect_identical(df, before)
})


# ---- the displayed statistic -----------------------------------------------

test_that("adjusted p is used when the table carries usable values", {
    df <- ora_fixture(c("GO:0000001", "GO:0000002"),
                      padj = c(0.20, 0.02), pvalue = c(1e-4, 2e-4))

    # Ranked on padj, so the row with the larger raw p leads.
    expect_identical(select_top_ora_per_collection(df, top_n = 1)$ID,
                     "GO:0000002")
})

test_that("a padj column present but empty falls back to the raw p-value", {
    # Values, not presence: an all-NA padj column would otherwise rank every
    # row equal and leave the figure ordered on identity alone.
    df <- ora_fixture(c("GO:0000001", "GO:0000002"),
                      padj = c(NA_real_, NA_real_), pvalue = c(0.04, 0.01))

    expect_identical(select_top_ora_per_collection(df, top_n = 1)$ID,
                     "GO:0000002")
})


# ---- the bar plot: opt-in selection, unconditional score -------------------

test_that("the default path ranks on padj, plain top-n, and classifies nothing", {
    # padj and the raw p order these in opposite directions, which is what a
    # table row-bound from several BH families can do. The old behaviour --
    # ranking on raw p -- would have led with hsa04110 and GO:0000003.
    df <- ora_fixture(c("GO:0000001", "GO:0000002", "GO:0000003", "hsa04110"),
                      padj   = c(0.001, 0.002, 0.003, 0.010),
                      pvalue = c(5e-4, 4e-4, 3e-4, 1e-4))
    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    drawn <- plot_multi_ora_barplot(df, "t", top_n = 2, kegg_org = "hsa")

    # Ranked on padj.
    expect_identical(drawn$key, c("GO:0000001", "GO:0000002"))
    expect_identical(drawn$score, c(0.001, 0.002))
    # Plain top-n: the KEGG map is displaced by the second GO term, because no
    # collection-aware selection ran.
    expect_false("04110" %in% drawn$key)
    # And nothing was classified, so no bar could have been tagged.
    expect_true(all(is.na(drawn$collection)))
})

test_that("opting in keeps the same score but changes which terms are shown", {
    df <- ora_fixture(c("GO:0000001", "GO:0000002", "GO:0000003", "hsa04110"),
                      padj   = c(0.001, 0.002, 0.003, 0.010),
                      pvalue = c(5e-4, 4e-4, 3e-4, 1e-4))
    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    drawn <- plot_multi_ora_barplot(df, "t", top_n = 2, by_collection = TRUE,
                                    kegg_org = "hsa")

    # Same statistic as the default path, still strongest-first.
    expect_identical(drawn$score, c(0.001, 0.010))
    # Different membership: the KEGG map takes the slot the second GO term had.
    expect_identical(drawn$key, c("GO:0000001", "04110"))
    expect_identical(drawn$collection, c("GO", "KEGG"))
    # Both collections named on a mixed figure.
    expect_true(all(grepl("^\\[(GO|KEGG)\\] ", drawn$label)))
})

test_that("a single-collection figure is not tagged with its one collection", {
    df <- ora_fixture(c("hsa04110", "hsa04115"), padj = c(0.01, 0.02))
    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    drawn <- plot_multi_ora_barplot(df, "t", by_collection = TRUE,
                                    kegg_org = "hsa")

    expect_identical(drawn$collection, c("KEGG", "KEGG"))
    expect_false(any(grepl("^\\[", drawn$label)))
})

test_that("the bar plot falls back to the raw p-value with no usable padj", {
    df <- ora_fixture(c("GO:0000001", "GO:0000002"),
                      padj = c(NA_real_, NA_real_), pvalue = c(0.04, 0.01))
    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    drawn <- plot_multi_ora_barplot(df, "t", top_n = 2)

    expect_identical(drawn$key, c("GO:0000002", "GO:0000001"))
    expect_identical(drawn$score, c(0.01, 0.04))
})

test_that("a row with no readable name still gets a bar", {
    df <- data.frame(ID = c("GO:0000001", NA), pathway = c(NA, NA),
                     padj = c(0.01, 0.02), stringsAsFactors = FALSE)
    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    drawn <- plot_multi_ora_barplot(df, "t")

    expect_equal(nrow(drawn), 2L)
    expect_identical(drawn$label, c("GO:0000001", "(unnamed pathway)"))
})

test_that("nothing to draw returns NULL rather than an empty figure", {
    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    expect_null(plot_multi_ora_barplot(NULL, "t"))
    expect_null(plot_multi_ora_barplot(ora_fixture(character(0)), "t"))
    # Neither significance column: there is no statistic to put on the axis.
    expect_null(plot_multi_ora_barplot(
        data.frame(ID = "GO:0000001", pathway = "a",
                   stringsAsFactors = FALSE), "t"))
})
