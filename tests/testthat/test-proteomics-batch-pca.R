# tests/testthat/test-proteomics-batch-pca.R
#
# The before/after batch-correction PCA (issue #354):
# mod_proteomics_batch_pca() and select_shared_pca_features(), both in
# R/modules/proteomics/01b_mod_batch_correction.R, loaded by helper.R.
#
# All data here is synthetic. The "correction" is faked by subtracting the
# injected batch shift, so the tests need neither sva nor proBatch.

make_batch_pca_fixture <- function(n_feat = 30, n_per_batch = 4, shift = 5) {
    withr::with_seed(42, {
        samples <- paste0("S", seq_len(2 * n_per_batch))
        batch   <- rep(c(1L, 2L), each = n_per_batch)   # numeric, like pair_no
        group   <- rep(rep(c("Ctrl", "Treat"), each = n_per_batch / 2), times = 2)
        raw <- matrix(rnorm(n_feat * length(samples), mean = 20, sd = 1),
                      nrow = n_feat,
                      dimnames = list(paste0("PROT", seq_len(n_feat)), samples))
    })
    raw[, batch == 2L] <- raw[, batch == 2L] + shift
    corrected <- raw
    corrected[, batch == 2L] <- corrected[, batch == 2L] - shift

    meta <- data.frame(SampleID = samples, Condition = group, pair_no = batch,
                       stringsAsFactors = FALSE)
    pre_raw <- list(expr_imp_single = raw, meta = meta)
    pre     <- list(expr_imp_single = corrected, meta = meta,
                    batch_corrected = TRUE, batch_method = "combat")
    config  <- list(modes = list(proteomics = list(
        effects = list(color = "Condition", samples = "SampleID"),
        batch_correction = list(method = "combat", batch_col = "pair_no")
    )))
    list(pre_raw = pre_raw, pre = pre, config = config)
}

test_that("select_shared_pca_features keeps proteins usable on both sides", {
    before <- matrix(c(1, 2, 3,
                       5, 5, 5,     # zero variance before
                       1, 2, 3,
                       1, NA, NA),  # one finite value
                     nrow = 4, byrow = TRUE,
                     dimnames = list(c("P1", "P2", "P3", "P4"), c("S1", "S2", "S3")))
    after <- before
    after["P2", ] <- c(1, 2, 3)     # usable after, not before
    after["P3", ] <- c(7, 7, 7)     # usable before, not after

    expect_identical(select_shared_pca_features(before, after), "P1")
})

test_that("select_shared_pca_features drops proteins missing from one side", {
    m <- matrix(c(1, 2, 3, 4, 6, 9), nrow = 2, byrow = TRUE,
                dimnames = list(c("P1", "P2"), c("S1", "S2", "S3")))
    expect_identical(select_shared_pca_features(m, m["P2", , drop = FALSE]), "P2")
})

test_that("no correction ran: nothing is written and a stale panel is removed", {
    fx <- make_batch_pca_fixture()
    out_dir <- withr::local_tempdir()
    stale <- file.path(out_dir, "Diagnostic_plots", "PCA_batch_before_after.png")
    dir.create(dirname(stale), recursive = TRUE)
    file.create(stale)

    # Configured, but the runtime record says it did not run.
    fx$pre$batch_corrected <- NULL
    res <- mod_proteomics_batch_pca(fx$pre_raw, fx$pre, fx$config, out_dir)

    expect_identical(res, character(0))
    expect_false(file.exists(stale))
})

test_that("a correction that ran writes the four-panel PNG", {
    skip_if_not_installed("patchwork")
    fx <- make_batch_pca_fixture()
    out_dir <- withr::local_tempdir()

    res <- mod_proteomics_batch_pca(fx$pre_raw, fx$pre, fx$config, out_dir)

    expect_identical(res, file.path(out_dir, "Diagnostic_plots",
                                    "PCA_batch_before_after.png"))
    expect_true(file.exists(res))
    expect_gt(file.size(res), 0)
})
