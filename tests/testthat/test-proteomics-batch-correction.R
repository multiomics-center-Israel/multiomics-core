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

# --- Condition-preservation column resolution (P1: no silent drop) ---

test_that("resolve_group_variable distinguishes unset from configured-but-missing", {
    pre  <- make_synth_prot()
    expr <- pre$expr_imp_single
    meta <- pre$meta

    # Nothing configured -> NULL (genuinely nothing to preserve)
    expect_null(resolve_group_variable(meta, expr, NULL, "SampleID", preserve = TRUE))
    expect_null(resolve_group_variable(meta, expr, "", "SampleID", preserve = TRUE))

    # Configured but absent + preservation requested -> loud error, not silence
    expect_error(
        resolve_group_variable(meta, expr, "Conditon", "SampleID", preserve = TRUE),
        "not found in metadata"
    )

    # Configured but absent + preservation NOT requested -> NULL (won't be used)
    expect_null(resolve_group_variable(meta, expr, "Conditon", "SampleID", preserve = FALSE))

    # Configured and present -> factor aligned to columns
    gv <- resolve_group_variable(meta, expr, "Condition", "SampleID", preserve = TRUE)
    expect_s3_class(gv, "factor")
    expect_length(gv, ncol(expr))
})

test_that("correct_batch_proteomics errors on a misspelled preservation column", {
    pre <- make_synth_prot()
    # valid batch column, but a typo'd group column with preserve_condition (default TRUE)
    config <- make_config(method = "combat", batch_col = "Batch", group_col = "Conditon")
    expect_error(correct_batch_proteomics(pre, config), "not found in metadata")
})

# --- Non-finite guard (P2: reject Inf, not just NA) ---

test_that("correct_batch_proteomics rejects non-finite intensities before correction", {
    pre <- make_synth_prot()
    pre$expr_imp_single[1, 1] <- -Inf   # e.g. log2(0) leaking through upstream
    config <- make_config(method = "combat", batch_col = "Batch")
    # aborts on the finite guard before any backend/package is touched
    expect_error(correct_batch_proteomics(pre, config), "finite")
})

# --- expr_filt keeps its unimputed contract after correction (P1/P2) ---

test_that("batch correction restores the NA pattern in expr_filt", {
    skip_if_not_installed("sva")
    pre <- make_synth_prot(batch_shift = 6)
    # A missing-value pattern in the (unimputed) filtered matrix; expr_imp_single
    # stays complete, exactly as preprocess_proteomics() guarantees.
    miss <- cbind(c(1L, 5L, 10L), c(2L, 3L, 7L))  # (row, col) positions
    pre$expr_filt[miss] <- NA
    config <- make_config(method = "combat", batch_col = "Batch", group_col = "Condition")

    out <- correct_batch_proteomics(pre, config, seed = 1)

    # expr_filt keeps NAs exactly where they were; work/imputed are complete
    expect_true(anyNA(out$expr_filt))
    expect_equal(which(is.na(out$expr_filt)), which(is.na(pre$expr_filt)))
    expect_false(anyNA(out$expr_work))
    expect_false(anyNA(out$expr_imp_single))
    # observed positions carry the corrected values; expr_filt_pre_imp mirrors expr_filt
    obs <- !is.na(out$expr_filt)
    expect_equal(out$expr_filt[obs], out$expr_work[obs])
    expect_identical(out$expr_filt_pre_imp, out$expr_filt)
})

