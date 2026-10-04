# The ORA producers say what they produced.
#
# Consumers must not infer over-representation from column shape or from the
# absence of GSEA fields -- build_ora_adjusted_p_matrix() requires an explicit
# method == "ora" and refuses a layer without it. That rule is only useful if
# the functions that genuinely run ORA, and nothing but ORA, stamp their own
# results. This file pins that they do, and that the functions which are not
# ORA-only were left alone.
#
# The stamp is metadata. Nothing about membership, p-values, adjustment,
# filtering or row order changes with it, which is what the compound ORA test
# below checks alongside the stamp itself.
#
# All fixtures are synthetic.


# ---- run_compound_ora() ----------------------------------------------------

test_that("compound ORA stamps its result without disturbing it", {
    # get_kegg_compound_pathways() reads this file when it is there, which is
    # what lets the whole producer run with no network. Ten measured compounds,
    # three of them significant.
    cache <- withr::local_tempdir()
    saveRDS(
        data.frame(
            pathway  = c(rep("map00010", 4), rep("map00020", 4), rep("map00030", 2)),
            compound = c(sprintf("C%05d", 1:4), sprintf("C%05d", 5:8),
                         "C00001", "C00009"),
            name     = c(rep("Glycolysis", 4), rep("Citrate cycle", 4),
                         rep("Pentose phosphate", 2)),
            stringsAsFactors = FALSE
        ),
        file.path(cache, "kegg_compound_pathways.rds")
    )

    universe <- sprintf("C%05d", 1:10)
    de_mapped <- data.frame(
        KEGG_ID = universe,
        pvalue  = c(rep(1e-4, 3), rep(0.5, 7)),
        padj    = c(rep(0.01, 3), rep(0.9, 7)),
        stringsAsFactors = FALSE
    )

    res <- suppressMessages(run_compound_ora(
        de_mapped, cache_dir = cache, min_gs = 2, max_gs = 500,
        pval_cutoff = 1, universe = universe))

    expect_true("method" %in% names(res))
    expect_identical(res$method, rep("ora", nrow(res)))

    # map00020 holds no significant compound and is not tested; the two that are
    # come back ordered by p-value, as they were before the stamp existed.
    expect_identical(res$ID, c("map00010", "map00030"))
    expect_true(all(diff(res$pvalue) >= 0))

    # The adjustment is still BH over what was tested, untouched by the stamp.
    expect_equal(res$padj, stats::p.adjust(res$pvalue, method = "BH"))

    # And the stamp was appended, not put in place of anything.
    expect_true(all(c("pathway", "ID", "pvalue", "padj", "GeneRatio", "BgRatio",
                      "setSize", "compounds") %in% names(res)))
})


# ---- the stamp is what the consumer acts on --------------------------------

test_that("a stamped producer frame reaches the ORA figure and an unstamped one does not", {
    # Same table twice, differing only in the column this change adds. This is
    # the whole point of stamping at the producer: without it the layer is
    # refused, and the refusal is not silent.
    stamped <- data.frame(
        pathway = "Glycolysis", ID = "map00010",
        pvalue = 1e-4, padj = 1e-3, method = "ora",
        stringsAsFactors = FALSE
    )
    unstamped <- stamped
    unstamped$method <- NULL

    m <- build_ora_adjusted_p_matrix(list(metabolomics = stamped), "00010",
                                     "metabolomics")
    expect_equal(unname(m["00010", "metabolomics"]), 1e-3)

    expect_warning(
        m2 <- build_ora_adjusted_p_matrix(list(metabolomics = unstamped), "00010",
                                          "metabolomics"),
        "metabolomics"
    )
    expect_true(is.na(m2["00010", "metabolomics"]))
})


# ---- run_ora_kegg() and run_ora_kegg_fisher() ------------------------------

test_that("return() inside tryCatch() leaves the enclosing function", {
    # The reason the stamp in run_ora_kegg() must sit in the frame construction
    # and not after the tryCatch(): the expression is a promise evaluated in the
    # caller's frame, so a return() inside it exits that frame, not just the
    # tryCatch. A stamp placed after the call is unreachable -- and reads as
    # though it were working, which is exactly how it got there.
    f <- function() {
        got <- tryCatch({
            return("left from inside")
            "never evaluated"
        }, error = function(e) "handler")
        "after the tryCatch"
    }

    expect_identical(f(), "left from inside")
})

test_that("the KEGG gene ORA producers stamp every result they return", {
    # Not reachable behaviourally here: every non-NULL path through these two
    # goes out to KEGG REST or through clusterProfiler, and this suite does not
    # make network calls or mock bindings. The assignment is pinned instead.
    ora_kegg_src <- paste(deparse(body(run_ora_kegg)), collapse = " ")
    fisher_src   <- paste(deparse(body(run_ora_kegg_fisher)), collapse = " ")

    # In the frame the clusterProfiler branch builds, so that both filtered
    # subsets carry it out through their own return().
    expect_true(grepl('method = "ora"', ora_kegg_src, fixed = TRUE))
    expect_true(grepl('df$method <- "ora"', fisher_src, fixed = TRUE))

    # And NOT after the tryCatch, where it cannot run. This is the assertion the
    # previous version of this test lacked: it checked only that a stamp existed
    # somewhere, which a dead one satisfies.
    expect_false(grepl('ora_res$method <- "ora"', ora_kegg_src, fixed = TRUE))

    # run_ora_kegg() has exactly two ways out: its own clusterProfiler branch,
    # stamped above, and run_ora_kegg_fisher(), which stamps its own.
    expect_true(grepl("run_ora_kegg_fisher(", ora_kegg_src, fixed = TRUE))
})

test_that("a producer that is not ORA-only was left unstamped", {
    # run_gsea_kegg() returns GSEA from its own branch and delegates to
    # run_ora_kegg() when that fails, so it is not authoritative about which of
    # the two a caller is holding. It must not claim ORA; the delegate stamps
    # the frames that really are.
    gsea_src <- paste(deparse(body(run_gsea_kegg)), collapse = " ")

    expect_false(grepl('method <- "ora"', gsea_src, fixed = TRUE))
    expect_true(grepl("run_ora_kegg(", gsea_src, fixed = TRUE))
})
