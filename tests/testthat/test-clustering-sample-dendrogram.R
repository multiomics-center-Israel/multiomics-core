# tests/testthat/test-clustering-sample-dendrogram.R
#
# Tests for clustering$cluster_samples (R/core/09_clustering.R). The
# hierarchical step used to hard-code cluster_cols = FALSE, so the heatmap
# never carried a sample dendrogram and there was no way to ask for one. The
# flag must default to FALSE, because every omic shares this code path and
# existing reports should not change silently.

resolve_cluster_samples <- function(cfg) isTRUE(cfg$clustering$cluster_samples)

test_that("sample clustering is off when the flag is absent", {
    expect_false(resolve_cluster_samples(list()))
    expect_false(resolve_cluster_samples(list(clustering = list())))
    expect_false(resolve_cluster_samples(list(clustering = list(enabled = TRUE))))
})

test_that("sample clustering is off when explicitly FALSE or NULL", {
    expect_false(resolve_cluster_samples(list(clustering = list(cluster_samples = FALSE))))
    expect_false(resolve_cluster_samples(list(clustering = list(cluster_samples = NULL))))
})

test_that("sample clustering is on only for TRUE", {
    expect_true(resolve_cluster_samples(list(clustering = list(cluster_samples = TRUE))))
})

test_that("non-logical values do not enable clustering", {
    # isTRUE() guards against a YAML string sneaking through as "true".
    expect_false(resolve_cluster_samples(list(clustering = list(cluster_samples = "true"))))
    expect_false(resolve_cluster_samples(list(clustering = list(cluster_samples = 1))))
})

test_that("the hierarchical step reads the flag and passes it on", {
    # testthat's working directory varies by runner, so locate the source
    # relative to the test file and skip rather than fail if it moved.
    candidates <- c(
        testthat::test_path("..", "..", "R", "core", "09_clustering.R"),
        "R/core/09_clustering.R"
    )
    f <- candidates[file.exists(candidates)][1]
    skip_if(is.na(f), "09_clustering.R not found from the test working directory")

    src <- readLines(f, warn = FALSE)
    expect_length(grep("cluster_samples <- isTRUE(cfg$clustering$cluster_samples)",
                       src, fixed = TRUE), 1)
    # The value must actually reach pheatmap, not merely be computed.
    expect_gte(length(grep("cluster_cols           = cluster_samples",
                           src, fixed = TRUE)), 1)
})
