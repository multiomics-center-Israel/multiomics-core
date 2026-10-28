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
    # qc$batch_col alone only colours QC plots; it must not put batch in DE.
    expect_false(bc$corrects)
    expect_equal(bc$de_mode, "none")
    expect_false(bc$de_mode_explicit)
})

test_that("get_proteomics_batch_config refuses an unknown de_mode", {
    for (m in c("model", "corrected", "")) {   # "corrected" was a draft name
        cfg <- list(batch_correction = list(method = "combat", de_mode = m))
        expect_error(get_proteomics_batch_config(cfg), "de_mode")
    }
})

test_that("get_proteomics_batch_config enables implicitly when method set", {
    cfg <- list(batch_correction = list(method = "combat", batch_col = "Batch"))
    bc <- get_proteomics_batch_config(cfg)
    expect_true(bc$enabled)
    expect_true(bc$corrects)
    expect_equal(bc$method, "combat")
    # A correction that runs puts batch in the limma model by default.
    expect_equal(bc$de_mode, "covariate")
})

test_that("de_mode defaults follow whether a correction runs, and explicit values win", {
    resolve <- function(...) get_proteomics_batch_config(list(batch_correction = list(...)))
    expect_equal(resolve(method = "combat", enabled = FALSE)$de_mode, "none")
    expect_equal(resolve(method = "none", enabled = TRUE)$de_mode, "none")
    expect_equal(resolve(method = "none", de_mode = "covariate")$de_mode, "covariate")
    lg <- resolve(method = "probatch", de_mode = "legacy_corrected")
    expect_equal(lg$de_mode, "legacy_corrected")
    expect_true(lg$de_mode_explicit)
    # Pure: resolving never warns, whatever the mode (validation owns that).
    expect_silent(resolve(method = "combat", de_mode = "legacy_corrected"))
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
    expect_equal(dim(out1$expr_imp_single), dim(pre$expr_filt))
    # The two QC matrices carry the correction identically ...
    expect_identical(out1$expr_work, out1$expr_imp_single)
    # ... and the measured matrix is untouched.
    expect_identical(out1$expr_filt, pre$expr_filt)
    # seeded -> deterministic
    expect_identical(out1$expr_imp_single, out2$expr_imp_single)

    diff_after <- mean(out1$expr_imp_single[, b2]) - mean(out1$expr_imp_single[, !b2])
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
    expect_equal(dim(out$expr_imp_single), dim(pre$expr_filt))
    expect_identical(out$expr_filt, pre$expr_filt)

    diff_after <- mean(out$expr_imp_single[, b2]) - mean(out$expr_imp_single[, !b2])
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
        "preservation"   # single word: robust to cli line-wrapping at any console width
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
    expect_error(correct_batch_proteomics(pre, config), "preservation")
})

# --- Non-finite guard (P2: reject Inf, not just NA) ---

test_that("correct_batch_proteomics rejects non-finite intensities before correction", {
    pre <- make_synth_prot()
    pre$expr_imp_single[1, 1] <- -Inf   # e.g. log2(0) leaking through upstream
    config <- make_config(method = "combat", batch_col = "Batch")
    # aborts on the finite guard before any backend/package is touched
    expect_error(correct_batch_proteomics(pre, config), "finite")
})

# --- Matrix contract: correction touches the QC matrices only ----------------

test_that("batch correction leaves the measured matrices untouched", {
    skip_if_not_installed("sva")
    pre <- make_synth_prot(batch_shift = 6)
    # A missing-value pattern in the (unimputed) filtered matrix; expr_imp_single
    # stays complete, exactly as preprocess_proteomics() guarantees.
    miss <- cbind(c(1L, 5L, 10L), c(2L, 3L, 7L))  # (row, col) positions
    pre$expr_filt[miss] <- NA
    pre$expr_filt_pre_imp <- pre$expr_filt
    config <- make_config(method = "combat", batch_col = "Batch", group_col = "Condition")

    out <- correct_batch_proteomics(pre, config, seed = 1)

    # Measured values -- the DE input, exports, CV, Shiny expr_raw -- as given.
    expect_identical(out$expr_filt, pre$expr_filt)
    expect_identical(out$expr_filt_pre_imp, pre$expr_filt_pre_imp)
    # QC matrices: complete and corrected.
    expect_false(anyNA(out$expr_imp_single))
    expect_identical(out$expr_work, out$expr_imp_single)
    obs <- !is.na(pre$expr_filt)
    expect_false(isTRUE(all.equal(out$expr_imp_single[obs], pre$expr_filt[obs])))
    # No parallel copy: the old draft's expr_filt_uncorrected is gone.
    expect_null(out$expr_filt_uncorrected)
})

# --- Shiny payload contract (the Shiny app itself is not touched) -------------

shiny_cfg <- function() {
    list(modes = list(proteomics = list(
        effects       = list(color = "Condition", samples = "SampleID"),
        id_columns    = list(protein_id = "Protein.Group", sample_col = "SampleID"),
        de_table      = list(id_col = "FeatureID"),
        normalization = list(method = "none"),
        de            = list()
    )))
}

payload_for <- function(pre) {
    suppressMessages(suppressWarnings(build_shiny_payload_proteomics(
        pre = pre, de_res = list(), inputs = list(contrasts = NULL),
        config = shiny_cfg(), final_results = NULL)))
}

test_that("Shiny expr_raw carries measured values and expr_norm the QC matrix", {
    # Package-free: a corrected QC matrix that differs from the measured one
    # everywhere, as after any correction.
    pre <- make_synth_prot()
    pre$expr_filt[cbind(c(2L, 4L), c(1L, 6L))] <- NA
    pre$expr_imp_single <- pre$expr_imp_single - 3
    pre$expr_work <- pre$expr_imp_single

    p <- payload_for(pre)
    expect_identical(p$expr_raw, pre$expr_filt)          # NAs where unmeasured
    expect_identical(p$expr_norm, pre$expr_imp_single)   # complete QC matrix
})

test_that("after a real correction the Shiny payload splits measured from corrected", {
    skip_if_not_installed("sva")
    pre <- make_synth_prot(batch_shift = 6)
    pre$expr_filt[cbind(c(2L, 4L), c(1L, 6L))] <- NA
    out <- correct_batch_proteomics(
        pre, make_config(method = "combat", batch_col = "Batch", group_col = "Condition"),
        seed = 1)

    p <- payload_for(out)
    # Semantic change from the pre-de_mode contract, where expr_raw held the
    # corrected values: it now holds what was measured.
    expect_identical(p$expr_raw, pre$expr_filt)
    expect_identical(p$expr_norm, out$expr_imp_single)
    obs <- !is.na(pre$expr_filt)
    expect_false(isTRUE(all.equal(p$expr_raw[obs], p$expr_norm[obs])))
})

