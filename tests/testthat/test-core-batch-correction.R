# Tests for the shared batch-correction backends (R/core/19_batch_correction.R)
# test-core-batch-correction.R
#
# All data here is synthetic.

# Eight samples, two batches. `layout = "balanced"` puts two Ctrl and two Treat
# in each batch; `layout = "unbalanced"` puts three Ctrl + one Treat in B1 and
# one Ctrl + three Treat in B2.
make_batch_fixture <- function(n_feat = 20, layout = c("balanced", "unbalanced"),
                               group_effect = 4, batch_shift = 0, noise = 0) {
    layout <- match.arg(layout)
    batch <- factor(rep(c("B1", "B2"), each = 4))
    group <- factor(if (layout == "balanced") {
        c("Ctrl", "Ctrl", "Treat", "Treat", "Ctrl", "Ctrl", "Treat", "Treat")
    } else {
        c("Ctrl", "Ctrl", "Ctrl", "Treat", "Ctrl", "Treat", "Treat", "Treat")
    })
    expr <- withr::with_seed(7, {
        base <- matrix(rnorm(n_feat, mean = 10, sd = 1), nrow = n_feat, ncol = 8)
        base + noise * matrix(rnorm(n_feat * 8), nrow = n_feat)
    })
    expr <- sweep(expr, 2, group_effect * (group == "Treat") +
                           batch_shift * (batch == "B2"), "+")
    dimnames(expr) <- list(paste0("F", seq_len(n_feat)), paste0("S", 1:8))
    list(expr = expr, batch = batch, group = group)
}

group_diff <- function(expr, group) {
    mean(expr[, group == "Treat"]) - mean(expr[, group == "Ctrl"])
}
batch_diff <- function(expr, batch) {
    mean(expr[, batch == "B2"]) - mean(expr[, batch == "B1"])
}

# --- Confounding check -------------------------------------------------------

test_that("assert_batch_not_confounded passes a design where batches share groups", {
    fx <- make_batch_fixture(layout = "unbalanced")
    expect_true(assert_batch_not_confounded(fx$batch, fx$group))
    expect_true(assert_batch_not_confounded(fx$batch))   # nothing to keep
})

test_that("assert_batch_not_confounded stops when each batch holds one group", {
    batch <- factor(rep(c("B1", "B2"), each = 4))
    group <- factor(rep(c("Ctrl", "Treat"), each = 4))
    expect_error(assert_batch_not_confounded(batch, group), "confounded")
})

test_that("assert_batch_not_confounded stops when the second batch repeats the first", {
    fx <- make_batch_fixture()
    expect_error(assert_batch_not_confounded(fx$batch, fx$group, batch2 = fx$batch),
                 "confounded")
})

test_that("unused factor levels are not read as confounding", {
    # A level left over from subsetting samples would be an all-zero column.
    fx <- make_batch_fixture(layout = "unbalanced", batch_shift = 2)
    group <- factor(as.character(fx$group), levels = c("Ctrl", "Treat", "Dropped"))
    batch <- factor(as.character(fx$batch), levels = c("B1", "B2", "B9"))
    covs <- data.frame(Sex = factor(rep(c("F", "M"), 4), levels = c("F", "M", "X")))

    expect_true(assert_batch_not_confounded(batch, group, covariates = covs))
    out <- correct_batch_limma(fx$expr, batch, group_var = group, covariates = covs)
    # The shift of 2 is removed; what is left between the batches is the
    # 3-vs-1 Treat composition, 4 * (3/4 - 1/4) = 2.
    expect_equal(batch_diff(fx$expr, fx$batch), 4, tolerance = 1e-8)
    expect_equal(batch_diff(out, fx$batch), 2, tolerance = 1e-8)
    expect_equal(group_diff(out, fx$group), 4, tolerance = 1e-8)
})

# --- limma backend -----------------------------------------------------------

test_that("correct_batch_limma removes an additive batch shift", {
    fx <- make_batch_fixture(batch_shift = 5, noise = 0.1)
    out <- correct_batch_limma(fx$expr, fx$batch, group_var = fx$group)

    expect_equal(dim(out), dim(fx$expr))
    expect_identical(dimnames(out), dimnames(fx$expr))
    expect_lt(abs(batch_diff(out, fx$batch)), 0.1 * abs(batch_diff(fx$expr, fx$batch)))
    expect_equal(group_diff(out, fx$group), 4, tolerance = 0.1)
})

test_that("correct_batch_limma keeps the group effect when batches are unbalanced", {
    # No batch effect at all. Batch means still differ, because B2 holds more
    # Treat samples. With the group in the design nothing is removed; estimated
    # from batch means alone, a quarter of the group effect would be lost.
    fx <- make_batch_fixture(layout = "unbalanced", batch_shift = 0)

    out <- correct_batch_limma(fx$expr, fx$batch, group_var = fx$group)
    expect_equal(group_diff(out, fx$group), 4, tolerance = 1e-8)
    expect_equal(out, fx$expr, tolerance = 1e-8)

    blind <- limma::removeBatchEffect(fx$expr, batch = fx$batch)
    expect_equal(group_diff(blind, fx$group), 3, tolerance = 1e-8)
})

test_that("correct_batch_limma corrects two batch factors and returns the result", {
    # Guards against a corrected matrix that is computed and then dropped.
    fx <- make_batch_fixture(batch_shift = 3)
    run <- factor(rep(c("R1", "R2"), times = 4))
    shifted <- sweep(fx$expr, 2, 2 * (run == "R2"), "+")

    out <- correct_batch_limma(shifted, fx$batch, batch2 = run, group_var = fx$group)

    expect_false(isTRUE(all.equal(out, shifted)))
    expect_lt(abs(batch_diff(out, fx$batch)), 1e-8)
    expect_lt(abs(mean(out[, run == "R2"]) - mean(out[, run == "R1"])), 1e-8)
    expect_equal(group_diff(out, fx$group), 4, tolerance = 1e-8)
})

test_that("correct_batch_limma stops on a confounded design", {
    fx <- make_batch_fixture()
    group <- factor(rep(c("Ctrl", "Treat"), each = 4))
    expect_error(correct_batch_limma(fx$expr, fx$batch, group_var = group), "confounded")
})

test_that("correct_batch_limma rejects missing and non-finite values", {
    fx <- make_batch_fixture()
    na_expr <- fx$expr
    na_expr[1, 1] <- NA
    expect_error(correct_batch_limma(na_expr, fx$batch, group_var = fx$group), "finite")
    inf_expr <- fx$expr
    inf_expr[2, 3] <- -Inf
    expect_error(correct_batch_limma(inf_expr, fx$batch, group_var = fx$group), "finite")
})

# --- ComBat backend (requires sva) -------------------------------------------

test_that("correct_batch_combat reduces the batch shift and is reproducible", {
    skip_if_not_installed("sva")
    fx <- make_batch_fixture(batch_shift = 5, noise = 0.3)

    out1 <- suppressMessages(correct_batch_combat(fx$expr, fx$batch, fx$group, seed = 1))
    out2 <- suppressMessages(correct_batch_combat(fx$expr, fx$batch, fx$group, seed = 1))

    expect_equal(dim(out1), dim(fx$expr))
    expect_identical(dimnames(out1), dimnames(fx$expr))
    expect_identical(out1, out2)
    expect_lt(abs(batch_diff(out1, fx$batch)), abs(batch_diff(fx$expr, fx$batch)))
})

test_that("correct_batch_combat with mean_only shifts batches without rescaling them", {
    skip_if_not_installed("sva")
    fx <- make_batch_fixture(batch_shift = 5, noise = 0.3)
    # Make B2 twice as spread out as B1, so a scale adjustment would show.
    b2 <- fx$batch == "B2"
    centre <- rowMeans(fx$expr[, b2])
    fx$expr[, b2] <- centre + 2 * (fx$expr[, b2] - centre)

    sd_in_b2 <- function(m) apply(m[, b2], 1, stats::sd)
    mean_only <- suppressMessages(
        correct_batch_combat(fx$expr, fx$batch, fx$group, mean_only = TRUE, seed = 1))
    full <- suppressMessages(
        correct_batch_combat(fx$expr, fx$batch, fx$group, mean_only = FALSE, seed = 1))

    expect_equal(sd_in_b2(mean_only), sd_in_b2(fx$expr), tolerance = 1e-8)
    expect_false(isTRUE(all.equal(sd_in_b2(full), sd_in_b2(fx$expr))))
})

test_that("correct_batch_combat needs mean_only for a single-sample batch", {
    # sva::ComBat switches to mean-only by itself when a batch has one sample;
    # refusing instead keeps the configured flag equal to what ran.
    skip_if_not_installed("sva")
    fx <- make_batch_fixture(noise = 0.3)
    batch <- factor(c("B1", "B1", "B1", "B1", "B2", "B2", "B2", "B3"))

    expect_error(correct_batch_combat(fx$expr, batch, fx$group, mean_only = FALSE),
                 "single sample")
    out <- suppressMessages(
        correct_batch_combat(fx$expr, batch, fx$group, mean_only = TRUE, seed = 1))
    expect_equal(dim(out), dim(fx$expr))
})

test_that("correct_batch_combat stops on a confounded design", {
    skip_if_not_installed("sva")
    fx <- make_batch_fixture()
    group <- factor(rep(c("Ctrl", "Treat"), each = 4))
    expect_error(correct_batch_combat(fx$expr, fx$batch, group), "confounded")
})

test_that("correct_batch_combat rejects missing and non-finite values", {
    skip_if_not_installed("sva")
    fx <- make_batch_fixture()
    fx$expr[1, 1] <- NA
    expect_error(correct_batch_combat(fx$expr, fx$batch, fx$group), "finite")
})
