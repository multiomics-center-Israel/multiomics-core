# Which statistic fgsea is handed, and what happens to ranks it cannot use.
#
# .build_fgsea_ranks() picks one ranking source for the whole contrast table:
# the Wald/t statistic where the table carries usable values, otherwise
# sign(log2FoldChange) * -log10(pvalue + 1e-300).
#
# The bug this covers was silent. load_precomputed_rna_de() always emits a
# `stat` column and fills it with NA when the source export carries no
# statistic, so a gate testing only whether the column exists selected a column
# of NAs, dropped every rank, handed fgsea an empty vector, and reported "no
# gene set overlap -- skipping" for every collection. Nothing in the log
# distinguished that from a layer with no enrichment, and the fallback that
# exists precisely for this case never fired.
#
# Two properties hold throughout: the choice is per table, never per row, and
# a rank the chosen source cannot produce is dropped rather than filled from
# the other source. Mixing two ranking scales in one vector is not a ranking.
#
# Fixtures are built so the two sources give visibly different numbers, which
# is how each test says which branch ran. All synthetic; no fgsea, no network.

# stat says one thing, the fallback says another, so the returned values name
# the branch: the statistic branch yields 3.2 / -1.7, the fallback yields
# 2 / -2 from log2FC = +/-1 at p = 0.01.
rank_fixture <- function(stat = NULL, log2FC = c(1, -1), pvalue = c(0.01, 0.01),
                         ids = c("g1", "g2")) {
    df <- data.frame(FeatureID = ids,
                     log2FoldChange = log2FC,
                     pvalue = pvalue,
                     stringsAsFactors = FALSE)
    if (!is.null(stat)) df$stat <- stat
    df
}

FALLBACK_RANKS <- c(g1 = 2, g2 = -2)


test_that("a usable statistic is the ranking source", {
    res <- rank_fixture(stat = c(3.2, -1.7))

    ranks <- .build_fgsea_ranks(res)

    expect_equal(ranks, c(g1 = 3.2, g2 = -1.7))
    # Not the fallback, which this fixture would have made 2 / -2.
    expect_false(isTRUE(all.equal(unname(ranks), unname(FALLBACK_RANKS))))
})

test_that("a stat column of all NA falls back instead of ranking on nothing", {
    # The live failure: the column is present, so a presence-only gate took it.
    res <- rank_fixture(stat = c(NA_real_, NA_real_))

    ranks <- .build_fgsea_ranks(res)

    expect_equal(ranks, FALLBACK_RANKS)
})

test_that("an all-NA stat never yields an empty ranking while the fallback is usable", {
    # Stated separately from the branch choice because this is the symptom that
    # reached the log: an empty vector into fgsea, reported as no overlap.
    res <- rank_fixture(stat = rep(NA_real_, 2))

    ranks <- .build_fgsea_ranks(res)

    expect_gt(length(ranks), 0)
    expect_true(all(is.finite(ranks)))
})

test_that("a partly usable statistic keeps the statistic branch and drops the rest", {
    # Michal's worked example: choose the stat branch because finite values
    # exist, drop the unusable row, do not fill it from log2FC/pvalue.
    res <- rank_fixture(stat = c(3.2, NA, -1.7),
                        log2FC = c(1, 1, -1),
                        pvalue = c(0.01, 0.01, 0.01),
                        ids = c("g1", "g2", "g3"))

    ranks <- .build_fgsea_ranks(res)

    expect_equal(ranks, c(g1 = 3.2, g3 = -1.7))
    # g2 is absent rather than carrying the fallback's value for that row.
    expect_false("g2" %in% names(ranks))
})

test_that("no stat column at all leaves the fallback as it was", {
    res <- rank_fixture(stat = NULL)

    expect_equal(.build_fgsea_ranks(res), FALLBACK_RANKS)
})

test_that("a non-numeric stat is unusable, and is not an error", {
    # A character column would sort lexically and rank nonsensically; a factor
    # is an integer vector underneath, which is why the gate tests is.numeric()
    # rather than relying on is.finite() to reject it.
    res_chr <- rank_fixture(stat = c("3.2", "-1.7"))
    res_fct <- rank_fixture(stat = factor(c("a", "b")))

    expect_no_error(ranks_chr <- .build_fgsea_ranks(res_chr))
    expect_equal(ranks_chr, FALLBACK_RANKS)

    expect_no_error(ranks_fct <- .build_fgsea_ranks(res_fct))
    expect_equal(ranks_fct, FALLBACK_RANKS)
})

test_that("Inf, -Inf and NaN never reach fgsea", {
    # !is.na() let the infinities through; fgsea does not reject them, it ranks
    # on them, seating a feature at an end of the list on no evidence.
    res <- rank_fixture(stat = c(Inf, 1, -Inf, NaN, 2),
                        log2FC = rep(1, 5), pvalue = rep(0.01, 5),
                        ids = paste0("g", 1:5))

    ranks <- .build_fgsea_ranks(res)

    expect_equal(ranks, c(g5 = 2, g2 = 1))
    expect_true(all(is.finite(ranks)))
})

test_that("the fallback drops its own unusable rows too", {
    # Same filter, other branch: a feature with no direction or no p-value has
    # no rank, and is left out rather than defaulted.
    res <- rank_fixture(stat = NULL,
                        log2FC = c(1, NA, -1),
                        pvalue = c(0.01, 0.01, NA),
                        ids = c("g1", "g2", "g3"))

    ranks <- .build_fgsea_ranks(res)

    expect_equal(ranks, c(g1 = 2))
})

test_that("ranks come back sorted strongest-first", {
    res <- rank_fixture(stat = c(-1, 5, 2), log2FC = rep(1, 3),
                        pvalue = rep(0.01, 3), ids = c("g1", "g2", "g3"))

    ranks <- .build_fgsea_ranks(res)

    expect_equal(ranks, c(g2 = 5, g3 = 2, g1 = -1))
    expect_false(is.unsorted(rev(ranks)))
})
