# tests/testthat/test-validate-pooled-script.R
#
# Helpers in scripts/validate_pooled_enrichment.R. The script runs on load, so
# only the named function definition is evaluated here.

load_script_fn <- function(name) {
    path <- testthat::test_path("..", "..", "scripts", "validate_pooled_enrichment.R")
    exprs <- parse(path)
    def <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                      identical(e[[2]], as.name(name)), as.list(exprs))
    env <- new.env()
    eval(def[[1]], env)
    env[[name]]
}

load_topk_overlap <- function() load_script_fn("topk_overlap")

test_that("top-k overlap reports the k it compared", {
    topk_overlap <- load_topk_overlap()
    a <- c(5, 4, 3, 2, 1, NA)
    b <- c(5, 4, 1, 2, 3, 0)
    ov <- topk_overlap(a, b, k = 100)
    expect_equal(ov[["k"]], 5)          # five features finite in both, not 100
    expect_equal(ov[["top"]], 1)
    ov2 <- topk_overlap(a, b, k = 2)
    expect_equal(ov2[["k"]], 2)
    expect_equal(ov2[["top"]], 1)
})

test_that("top-k overlap with nothing finite is empty, not a perfect overlap", {
    topk_overlap <- load_topk_overlap()
    ov <- topk_overlap(c(NA, Inf), c(1, 2), k = 100)
    expect_equal(ov[["k"]], 0)
    expect_true(is.na(ov[["top"]]) && is.na(ov[["bottom"]]))
})

test_that("Spearman is NA, not an error, with fewer than two finite pairs", {
    spearman_finite <- load_script_fn("spearman_finite")
    expect_identical(spearman_finite(numeric(0), numeric(0)), NA_real_)
    expect_identical(spearman_finite(c(NA, 1, Inf), c(1, 2, 3)), NA_real_)
    a <- c(1, 2, 3, 4, NA)
    b <- c(2, 1, 4, 3, 5)
    expect_equal(spearman_finite(a, b),
                 stats::cor(a, b, method = "spearman", use = "complete.obs"))
})

test_that("option B runs on a single imputation draw", {
    # sapply() over one fit returns a features x 1 matrix, so rowMeans() and
    # rubin_pool_estimates() see the same shape as with several draws; with one
    # draw nothing is pooled and the result is the single fit squeezed once.
    option_b <- load_script_fn("option_b")
    ids <- paste0("P", 1:6)
    fit <- list(
        coefficients   = matrix(c(1.2, -0.4, 0.3, 2.1, -1.5, 0.1), ncol = 1,
                                dimnames = list(ids, "A-B")),
        stdev.unscaled = matrix(c(0.5, 0.6, 0.5, 0.7, 0.5, 0.6), ncol = 1,
                                dimnames = list(ids, "A-B")),
        sigma          = c(0.8, 1.1, 0.9, 1.4, 0.7, 1.0),
        df.residual    = rep(4, 6)
    )
    z <- option_b(list(list(fit2 = fit)), "A-B", "residual")

    su2 <- fit$stdev.unscaled[, 1]^2
    sq <- limma::squeezeVar(fit$sigma^2, df = fit$df.residual)
    expected <- limma::zscoreT(fit$coefficients[, 1] / sqrt(sq$var.post * su2),
                               df = fit$df.residual + sq$df.prior)
    expect_length(z, 6L)
    expect_true(all(is.finite(z)))
    expect_equal(unname(z), unname(expected))
})
