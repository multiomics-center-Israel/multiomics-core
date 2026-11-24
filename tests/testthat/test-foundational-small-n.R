# tests/testthat/test-foundational-small-n.R
#
# Two foundational side analyses used to abort on ordinary runs:
#   - consensus clustering passed maxK = floor(n / 3), which is 2 below nine
#     samples, and ConsensusClusterPlus aborts on maxK = 2;
#   - pathway overlap read each layer's DE result as one table, but the
#     pipeline passes per-contrast lists, and the fallback GO step was fixed to
#     the human OrgDb whatever the organism.
# All data here are synthetic.

small_n_fixture <- function(n) {
    samples <- paste0("S", seq_len(n))
    withr::with_seed(7, list(
        omicsA = list(normalized_matrix = matrix(rnorm(40 * n), nrow = 40,
                      dimnames = list(paste0("a", 1:40), samples))),
        omicsB = list(normalized_matrix = matrix(rnorm(30 * n), nrow = 30,
                      dimnames = list(paste0("b", 1:30), samples)))
    ))
}

test_that("consensus clustering runs on 8 samples", {
    skip_if_not_installed("ConsensusClusterPlus")
    h <- small_n_fixture(8)
    res <- suppressMessages(
        compute_foundational_consensus_clustering(h, paste0("S", 1:8), list()))
    expect_false(is.null(res))
    expect_true(res$k_best >= 2)
    expect_length(res$clusters, 8)
})

test_that("consensus clustering skips with a message below 6 samples", {
    h <- small_n_fixture(5)
    expect_message(
        res <- compute_foundational_consensus_clustering(h, paste0("S", 1:5), list()),
        "at least 6 common samples")
    expect_null(res)
})

per_contrast_harmonized <- function() {
    tab <- data.frame(gene_symbol = paste0("G", 1:20), padj = rep(0.01, 20))
    list(
        omicsA = list(de_table = list(B_vs_A = tab, C_vs_A = tab)),
        omicsB = list(de_table = list(B_vs_A = tab))
    )
}

test_that("fallback enrichment skips an organism without an OrgDb", {
    skip_if_not_installed("clusterProfiler")
    expect_message(
        res <- run_basic_enrichment(per_contrast_harmonized(), NULL, list(),
                                    organism = "Entamoeba histolytica"),
        "No OrgDb")
    expect_identical(res, list())
    expect_message(
        res2 <- run_basic_enrichment(per_contrast_harmonized(), NULL, list(),
                                     organism = NULL),
        "No OrgDb")
    expect_identical(res2, list())
})

test_that("fallback enrichment skips per-contrast DE lists instead of aborting", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("org.Hs.eg.db")
    msgs <- character(0)
    res <- withCallingHandlers(
        run_basic_enrichment(per_contrast_harmonized(), NULL, list(),
                             organism = "Homo sapiens"),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    expect_length(res, 0)
    expect_true(any(grepl("per contrast", msgs)))
})

test_that("pathway overlap returns NULL, not an error, on per-contrast DE lists", {
    skip_if_not_installed("clusterProfiler")
    res <- suppressMessages(
        analyze_pathway_overlap(per_contrast_harmonized(), NULL, list(),
                                out_dir = NULL, organism = "Entamoeba histolytica"))
    expect_null(res)
})
