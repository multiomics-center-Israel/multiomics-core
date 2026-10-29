# tests/testthat/test-foundational-concordance-seed.R
#
# The sample-concordance step of the foundational analysis draws random
# permutations (Mantel test) and runs k-means. Both take fc$seed, so repeated
# runs give identical p-values and clusters whatever the global RNG state.
# All data here are synthetic.

concordance_fixture <- function() {
    samples <- paste0("S", 1:6)
    withr::with_seed(11, list(
        omicsA = list(normalized_matrix = matrix(rnorm(180), nrow = 30,
                      dimnames = list(paste0("a", 1:30), samples))),
        omicsB = list(normalized_matrix = matrix(rnorm(120), nrow = 20,
                      dimnames = list(paste0("b", 1:20), samples)))
    ))
}

test_that("get_foundational_config resolves the seed from foundational, params, default", {
    expect_equal(get_foundational_config(list())$seed, 42)
    expect_equal(get_foundational_config(list(params = list(seed = 5)))$seed, 5)
    cfg <- list(params = list(seed = 5),
                modes = list(multiomics = list(foundational = list(seed = 9))))
    expect_equal(get_foundational_config(cfg)$seed, 9)
})

test_that("Mantel permutation p-values repeat exactly with the same seed", {
    harmonized <- concordance_fixture()
    fc <- list(top_variable_features = 20, correlation_method = "spearman", seed = 3)
    samples <- paste0("S", 1:6)

    r1 <- suppressMessages(compute_sample_rank_correlations(harmonized, samples, fc))
    invisible(stats::runif(5))   # move the global RNG state between the calls
    r2 <- suppressMessages(compute_sample_rank_correlations(harmonized, samples, fc))

    expect_identical(r1$mantel_pvalues, r2$mantel_pvalues)
})

test_that("clustering consistency repeats exactly with the same seed", {
    harmonized <- concordance_fixture()
    fc <- list(seed = 3)
    samples <- paste0("S", 1:6)

    c1 <- suppressMessages(compute_clustering_consistency(harmonized, samples, NULL, fc, list()))
    invisible(stats::runif(5))
    c2 <- suppressMessages(compute_clustering_consistency(harmonized, samples, NULL, fc, list()))

    expect_identical(c1, c2)
})
