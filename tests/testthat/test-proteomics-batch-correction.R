# Tests for proteomics batch correction
# test-proteomics-batch-correction.R
#
# All data here is synthetic.

# --- Synthetic fixtures ---

make_synth_prot <- function(n_feat = 30, n_per_batch = 4, batch_shift = 5) {
    withr::with_seed(123, {
        samples <- c(
            paste0("S", seq_len(2 * n_per_batch))
        )
        batch <- rep(c("B1", "B2"), each = n_per_batch)
        group <- rep(rep(c("Ctrl", "Treat"), each = n_per_batch / 2), times = 2)

        base <- matrix(rnorm(n_feat * length(samples), mean = 20, sd = 1),
                       nrow = n_feat,
                       dimnames = list(paste0("PROT", seq_len(n_feat)), samples))
        # inject an additive batch effect on B2
        base[, batch == "B2"] <- base[, batch == "B2"] + batch_shift
    })

    meta <- data.frame(
        SampleID  = samples,
        Condition = group,
        Batch     = batch,
        stringsAsFactors = FALSE
    )

    list(
        expr_raw        = base,
        expr_filt       = base,
        expr_work       = base,
        expr_imp_single = base,
        meta            = meta,
        row_data        = data.frame(Protein.Group = rownames(base),
                                     stringsAsFactors = FALSE)
    )
}

make_config <- function(method = "none", batch_col = "Batch", ...) {
    extra <- list(...)
    bc <- c(list(method = method, batch_col = batch_col), extra)
    list(modes = list(proteomics = list(
        effects = list(color = "Condition", samples = "SampleID"),
        id_columns = list(sample_col = "SampleID"),
        batch_correction = bc
    )))
}

# --- Config resolution ---

test_that("get_proteomics_batch_config applies defaults", {
    cfg <- list(qc = list(batch_col = "Plate"))
    bc <- get_proteomics_batch_config(cfg)
    expect_equal(bc$method, "none")
    expect_false(bc$enabled)
    expect_equal(bc$batch_col, "Plate")        # falls back to qc$batch_col
    expect_true(bc$preserve_condition)
})

test_that("get_proteomics_batch_config enables implicitly when method set", {
    cfg <- list(batch_correction = list(method = "combat", batch_col = "Batch"))
    bc <- get_proteomics_batch_config(cfg)
    expect_true(bc$enabled)
    expect_equal(bc$method, "combat")
})

# --- Pass-through / dispatch ---

test_that("correct_batch_proteomics is a no-op when disabled", {
    pre <- make_synth_prot()
    config <- make_config(method = "none")
    out <- correct_batch_proteomics(pre, config)
    expect_identical(out$expr_filt, pre$expr_filt)
    expect_null(out$batch_corrected)
})

test_that("correct_batch_proteomics errors on unknown method", {
    pre <- make_synth_prot()
    config <- make_config(method = "bogus")
    expect_error(correct_batch_proteomics(pre, config), "Unknown")
})

test_that("correct_batch_proteomics errors on missing batch column", {
    pre <- make_synth_prot()
    config <- make_config(method = "combat", batch_col = "DoesNotExist")
    expect_error(correct_batch_proteomics(pre, config), "not found")
})

test_that("correct_batch_proteomics errors when batch has <2 levels", {
    pre <- make_synth_prot()
    pre$meta$Batch <- "OnlyOne"
    config <- make_config(method = "combat", batch_col = "Batch")
    expect_error(correct_batch_proteomics(pre, config), "fewer than 2 levels")
})

# --- ComBat backend (requires sva) ---

test_that("ComBat reduces the injected batch effect and is reproducible", {
    skip_if_not_installed("sva")
    pre <- make_synth_prot(batch_shift = 6)
    config <- make_config(method = "combat", batch_col = "Batch",
                          group_col = "Condition")

    b2 <- pre$meta$Batch == "B2"
    diff_before <- mean(pre$expr_filt[, b2]) - mean(pre$expr_filt[, !b2])

    out1 <- correct_batch_proteomics(pre, config, seed = 1)
    out2 <- correct_batch_proteomics(pre, config, seed = 1)

    expect_true(out1$batch_corrected)
    expect_equal(out1$batch_method, "combat")
    expect_equal(dim(out1$expr_filt), dim(pre$expr_filt))
    # work / imputed matrices are propagated identically
    expect_identical(out1$expr_work, out1$expr_filt)
    expect_identical(out1$expr_imp_single, out1$expr_filt)
    # seeded -> deterministic
    expect_identical(out1$expr_filt, out2$expr_filt)

    diff_after <- mean(out1$expr_filt[, b2]) - mean(out1$expr_filt[, !b2])
    expect_lt(abs(diff_after), abs(diff_before))
})

# --- probatch backend (requires proBatch) ---

test_that("probatch dispatch errors cleanly when package is absent", {
    if (requireNamespace("proBatch", quietly = TRUE)) {
        skip("proBatch installed; absence-path not exercised")
    }
    pre <- make_synth_prot()
    config <- make_config(method = "probatch", batch_col = "Batch")
    expect_error(correct_batch_proteomics(pre, config), "proBatch")
})

test_that("probatch centers per-batch medians when available", {
    skip_if_not_installed("proBatch")
    pre <- make_synth_prot(batch_shift = 6)
    config <- make_config(method = "probatch", batch_col = "Batch")

    b2 <- pre$meta$Batch == "B2"
    diff_before <- mean(pre$expr_filt[, b2]) - mean(pre$expr_filt[, !b2])

    out <- correct_batch_proteomics(pre, config)
    expect_true(out$batch_corrected)
    expect_equal(out$batch_method, "probatch")
    expect_equal(dim(out$expr_filt), dim(pre$expr_filt))

    diff_after <- mean(out$expr_filt[, b2]) - mean(out$expr_filt[, !b2])
    expect_lt(abs(diff_after), abs(diff_before))
})

# --- Module wrapper ---

test_that("mod_proteomics_batch_correction passes through when disabled", {
    pre <- make_synth_prot()
    config <- make_config(method = "none")
    out <- mod_proteomics_batch_correction(pre, config, verbose = FALSE)
    expect_identical(out$expr_filt, pre$expr_filt)
})
