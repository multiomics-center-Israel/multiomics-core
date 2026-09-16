# tests/testthat/test-enrichment-go-classification.R
#
# cluster_enrichment_terms() (R/core/09_enrichment.R) chooses between rrvgo
# semantic clustering, which is only meaningful for GO terms, and Jaccard
# clustering on gene overlap. Collections are now named after their GMT file, so
# the name is not evidence of what the identifiers are: the choice is made from
# the pathway ids.
#
# Both paths are separated cleanly here by using an organism with no OrgDb,
# whether or not rrvgo is installed: .cluster_go_terms() cannot resolve one and
# returns NULL, while .cluster_by_jaccard() returns the clustered rows. All data
# are synthetic.

clust_df <- function(ids) {
    data.frame(pathway = ids, padj = rep(0.001, length(ids)),
               stringsAsFactors = FALSE)
}

clust_sets <- function(ids) {
    stats::setNames(lapply(seq_along(ids), function(i) paste0("g", i:(i + 2))), ids)
}

cluster_with <- function(ids, database) {
    suppressMessages(cluster_enrichment_terms(
        enrichment_df = clust_df(ids),
        database      = database,
        gene_sets     = clust_sets(ids),
        organism      = "Example organism",
        threshold     = 0.7
    ))
}

test_that("a collection named like GO but carrying non-GO ids takes the Jaccard path", {
    # GOLD_domains is not Gene Ontology. Under the previous name-based test it
    # was sent to semantic clustering regardless of its identifiers.
    ids <- c("GOLD001", "GOLD002", "GOLD003")

    res <- cluster_with(ids, database = "GOLD_domains")

    expect_false(is.null(res))
    expect_true(all(c("cluster", "parentTerm") %in% names(res)))
    expect_setequal(res$pathway, ids)
})

test_that("real GO ids take the GO path whatever the collection is called", {
    # Named as a custom collection, but the ids are GO terms. The GO path cannot
    # resolve an OrgDb for this organism and returns NULL; Jaccard would have
    # returned the rows, as the test above shows on the same shape of input.
    ids <- c("GO:0006915", "GO:0008150", "GO:0016020")

    res <- cluster_with(ids, database = "my_custom_sets")

    expect_null(res)
})

test_that("the database argument no longer decides the path", {
    # Same ids under opposite database names: the outcome follows the ids.
    non_go <- c("IPR001", "IPR002", "IPR003")
    go     <- c("GO:0006915", "GO:0008150", "GO:0016020")

    expect_setequal(cluster_with(non_go, "GO_something")$pathway, non_go)
    expect_setequal(cluster_with(non_go, "InterPro")$pathway, non_go)
    expect_null(cluster_with(go, "GO_something"))
    expect_null(cluster_with(go, "InterPro"))
})
