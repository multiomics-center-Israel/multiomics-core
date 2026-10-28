# tests/testthat/test-proteomics-de-batch-covariate.R
#
# Batch as a fixed effect in the pooled limma fit (batch_correction$de_mode):
#   - resolve_de_batch_covariate(), build_proteomics_de_design(),
#     assess_batch_design(), run_limma_proteomics()  (05_de_summary.R)
#   - mod_proteomics_de(): DE input per mode, multiple imputation, diagnostics
#   - pool_limma_imputations() on batch-adjusted runs  (05d_rubin_pooled_stat.R)
#   - validate_proteomics_config()                      (90_config_validate.R)
#
# Every check here is deterministic: seeded synthetic data, compared with a
# direct limma fit or with exact structural properties. The statistical
# comparison of old vs new false-positive rates is an experiment, not a test:
# utils/experiment_batch_de_false_positives.R.

# Two conditions crossed with two batches, three replicates per cell. Batch B2
# carries an additive shift; the first 30 proteins carry a treatment effect.
make_batch_fixture <- function(n_prot = 120L, batch_shift = 3, effect = 0.8) {
    samples <- sprintf("S%02d", 1:12)
    group   <- rep(c("Ctrl", "Trt"), times = 6)
    batch   <- rep(c("B1", "B2"), each = 6)

    expr <- withr::with_seed(7, {
        m <- matrix(rnorm(n_prot * 12, mean = 20, sd = 0.3), nrow = n_prot)
        m[, batch == "B2"] <- m[, batch == "B2"] + batch_shift
        m[1:30, group == "Trt"] <- m[1:30, group == "Trt"] + effect
        m
    })
    rownames(expr) <- paste0("P", seq_len(n_prot))
    colnames(expr) <- samples

    meta <- data.frame(SampleName = samples, Group = group, Batch = batch,
                       stringsAsFactors = FALSE)
    prot_tbl <- data.frame(
        Protein.Group = rownames(expr), Protein.Names = rownames(expr),
        Genes = rownames(expr), First.Protein.Description = rownames(expr),
        stringsAsFactors = FALSE
    )
    contrasts_df <- data.frame(Contrast_name = "Trt_vs_Ctrl", Factor = "Group",
                               Numerator = "Trt", Denominator = "Ctrl",
                               stringsAsFactors = FALSE)
    list(expr = expr, meta = meta, prot_tbl = prot_tbl, contrasts_df = contrasts_df)
}

batch_cfg <- function(de_mode = "covariate", method = "none", batch_col = "Batch",
                      block_col = NULL, imputation = NULL) {
    list(
        params = list(seed = 42),
        modes = list(proteomics = list(
            effects    = list(samples = "SampleName", color = "Group"),
            id_columns = list(
                protein_id    = "Protein.Group",
                sample_col    = "SampleName",
                protein_annot = c("Protein.Group", "Protein.Names", "Genes",
                                  "First.Protein.Description")
            ),
            de_table   = list(id_col = "FeatureID"),
            de         = list(method = "limma", p_adjust_method = "BH", p_cutoff = 0.05,
                              linear_fc_cutoff = 1.5, block_col = block_col),
            imputation = imputation %||% list(method = "perseus_like", width = 0.3,
                                              downshift = 1.8, multi_imputation = TRUE,
                                              no_repetitions = 3, min_no_passed = 2),
            batch_correction = list(method = method, batch_col = batch_col,
                                    de_mode = de_mode)
        ))
    )
}

run_fixture <- function(fx, cfg, meta = fx$meta) {
    run_limma_proteomics(expr_imp = fx$expr, meta = meta,
                         contrasts_df = fx$contrasts_df, prot_tbl = fx$prot_tbl,
                         cfg = cfg)
}

# The reference fit, written out by hand rather than through any pipeline helper.
direct_limma <- function(expr, meta) {
    X <- cbind(Ctrl = as.numeric(meta$Group == "Ctrl"),
               Trt  = as.numeric(meta$Group == "Trt"),
               B2   = as.numeric(meta$Batch == "B2"))
    cm <- limma::makeContrasts(Trt_vs_Ctrl = Trt - Ctrl, levels = X)
    limma::eBayes(limma::contrasts.fit(limma::lmFit(expr, X), cm))
}

# Small design helper for assess_batch_design(): a contrast matrix over a design.
design_and_contrasts <- function(group, batch, contrasts) {
    g <- factor(group)
    d <- build_proteomics_de_design(g, factor(batch))
    cm <- limma::makeContrasts(contrasts = unname(contrasts), levels = d)
    colnames(cm) <- names(contrasts)   # as run_limma_proteomics() names them
    list(design = d, group = g, batch = factor(batch), cm = cm)
}

# --- resolver -------------------------------------------------------------

test_that("resolve_de_batch_covariate is NULL unless de_mode resolves to covariate", {
    meta <- make_batch_fixture()$meta
    p <- function(...) batch_cfg(...)$modes$proteomics
    expect_null(resolve_de_batch_covariate(meta, p(de_mode = "none")))
    expect_null(resolve_de_batch_covariate(meta, p(de_mode = "legacy_corrected",
                                                   method = "combat")))
    # qc-only batch column, no correction, no de_mode: not modelled.
    qc_only <- p(de_mode = NULL)
    qc_only$batch_correction <- NULL
    qc_only$qc <- list(batch_col = "Batch")
    expect_null(resolve_de_batch_covariate(meta, qc_only))
    # A configured correction makes covariate the default.
    expect_s3_class(resolve_de_batch_covariate(meta, p(de_mode = NULL, method = "combat")),
                    "factor")
})

test_that("resolve_de_batch_covariate refuses unusable batch columns", {
    meta <- make_batch_fixture()$meta
    p <- function(...) batch_cfg(...)$modes$proteomics
    expect_error(resolve_de_batch_covariate(meta, p(batch_col = NULL)), "no batch column")
    expect_error(resolve_de_batch_covariate(meta, p(batch_col = "Plate")),
                 "not in the sample metadata")
    na_meta <- meta
    na_meta$Batch[3] <- NA
    expect_error(resolve_de_batch_covariate(na_meta, p()), "missing values")
    one <- meta
    one$Batch <- "B1"
    expect_error(resolve_de_batch_covariate(one, p()), "single level")
})

# --- design build and assessment -------------------------------------------

test_that("without batch the design is one column per condition, as before", {
    d <- build_proteomics_de_design(factor(c("Ctrl", "Ctrl", "Trt", "Trt")))
    expect_equal(colnames(d), c("Ctrl", "Trt"))
    expect_equal(unname(d[, "Trt"]), c(0, 0, 1, 1))
})

test_that("a batch column clashing with a condition name never renames the condition", {
    d <- build_proteomics_de_design(factor(c("batch_B2", "batch_B2", "X", "X", "batch_B2", "X")),
                                    factor(c("B1", "B1", "B1", "B2", "B2", "B2")))
    expect_equal(colnames(d)[1:2], c("batch_B2", "X"))
    expect_false(anyDuplicated(colnames(d)) > 0)
})

test_that("a balanced design passes with no notes and no SE inflation", {
    skip_if_not_installed("limma")
    fx <- make_batch_fixture()
    x <- design_and_contrasts(fx$meta$Group, fx$meta$Batch, c(Trt_vs_Ctrl = "Trt - Ctrl"))
    a <- assess_batch_design(x$design, x$cm, x$group, x$batch)
    expect_equal(a$df_residual, 9L)
    expect_length(a$notes, 0)
    expect_equal(a$contrasts$SE_inflation, 1)     # batch orthogonal to condition
    expect_equal(a$contrasts$Shared_batches, "B1, B2")
    expect_equal(a$contrasts$Note, "")
    expect_equal(a$cells$B1, c(3L, 3L))
})

test_that("total confounding fails, naming the rank deficiency and the lost contrast", {
    skip_if_not_installed("limma")
    x <- design_and_contrasts(c("A", "A", "A", "B", "B", "B"),
                              c("B1", "B1", "B1", "B2", "B2", "B2"),
                              c(B_vs_A = "B - A"))
    err <- expect_error(assess_batch_design(x$design, x$cm, x$group, x$batch), "rank-deficient")
    expect_match(conditionMessage(err), "not estimable in this design: B_vs_A", fixed = TRUE)
    expect_match(conditionMessage(err), "B1 = {A}; B2 = {B}", fixed = TRUE)
})

test_that("partial confounding names only the contrasts that are not estimable", {
    skip_if_not_installed("limma")
    # C exists only in B2, so batch_B2 is aliased with C. B - A stays estimable.
    grp <- c("A", "A", "B", "B", "C", "C")
    bat <- c("B1", "B1", "B1", "B1", "B2", "B2")
    x <- design_and_contrasts(grp, bat, c(B_vs_A = "B - A", C_vs_A = "C - A"))
    err <- expect_error(assess_batch_design(x$design, x$cm, x$group, x$batch), "rank-deficient")
    expect_match(conditionMessage(err), "not estimable in this design: C_vs_A.", fixed = TRUE)
    expect_false(grepl("B_vs_A", conditionMessage(err), fixed = TRUE))

    # Rank-deficient although every requested contrast is estimable: still
    # refused, and the message says why rather than calling it confounded.
    y <- design_and_contrasts(grp, bat, c(B_vs_A = "B - A"))
    err2 <- expect_error(assess_batch_design(y$design, y$cm, y$group, y$batch), "rank-deficient")
    expect_match(conditionMessage(err2), "Every requested contrast is estimable in principle",
                 fixed = TRUE)
})

test_that("a design with no residual degrees of freedom is refused", {
    skip_if_not_installed("limma")
    x <- design_and_contrasts(c("A", "B", "A", "B"), c("B1", "B1", "B2", "B3"),
                              c(B_vs_A = "B - A"))
    expect_error(assess_batch_design(x$design, x$cm, x$group, x$batch),
                 "no residual degrees of freedom")
})

test_that("a contrast bridged only through other conditions is noted, not refused", {
    skip_if_not_installed("limma")
    # A only in B1, B only in B2, one C in each: B - A is estimable only through
    # the C pair, which carries the whole batch estimate. By hand: Var(B - A)
    # goes from 1/2 + 1/2 to 1/2 + 1/2 + 2, an SE inflation of sqrt(3).
    x <- design_and_contrasts(c("A", "A", "C", "B", "B", "C"),
                              c("B1", "B1", "B1", "B2", "B2", "B2"),
                              c(B_vs_A = "B - A", C_vs_A = "C - A"))
    a <- assess_batch_design(x$design, x$cm, x$group, x$batch)
    ba <- a$contrasts[a$contrasts$Contrast == "B_vs_A", ]
    ca <- a$contrasts[a$contrasts$Contrast == "C_vs_A", ]
    expect_equal(ba$Shared_batches, "")
    expect_match(ba$Note, "No batch holds both sides", fixed = TRUE)
    expect_equal(ba$SE_inflation, round(sqrt(3), 3))
    expect_gt(ba$SE_inflation, batch_se_inflation_note_at())
    expect_match(ba$Note, "heuristic note", fixed = TRUE)
    expect_equal(ca$Shared_batches, "B1")
    expect_true(any(grepl("cells are empty", a$notes, fixed = TRUE)))
    expect_true(any(grepl("2 residual degree", a$notes, fixed = TRUE)))
})

# --- fixed batch vs random block -------------------------------------------

test_that("a block partition identical to the batch partition is refused by grouping, not name", {
    skip_if_not_installed("limma")
    fx <- make_batch_fixture()
    x <- design_and_contrasts(fx$meta$Group, fx$meta$Batch, c(Trt_vs_Ctrl = "Trt - Ctrl"))
    same_groups <- factor(ifelse(fx$meta$Batch == "B1", "run_a", "run_b"))
    expect_error(assess_batch_design(x$design, x$cm, x$group, x$batch, block = same_groups),
                 "exactly as the batch column")
})

test_that("blocks nested in batches fit, and duplicateCorrelation sees the batch term", {
    skip_if_not_installed("limma")
    fx <- make_batch_fixture()
    meta <- fx$meta
    meta$Donor <- rep(c("D1", "D2", "D3", "D4"), each = 3)   # two donors per batch
    cfg <- batch_cfg(block_col = "Donor")

    X <- cbind(Ctrl = as.numeric(meta$Group == "Ctrl"), Trt = as.numeric(meta$Group == "Trt"),
               batch_B2 = as.numeric(meta$Batch == "B2"))
    expected <- limma::duplicateCorrelation(fx$expr, X, block = meta$Donor)$consensus
    if (!is.finite(expected)) skip("duplicateCorrelation consensus not finite on this fixture")

    expect_message(res <- run_fixture(fx, cfg, meta = meta),
                   sprintf("consensus within-block correlation = %.3f", expected), fixed = TRUE)
    expect_equal(colnames(res$design), c("Ctrl", "Trt", "batch_B2"))
    expect_length(res$batch_design$notes, 0)
})

# --- the single fit --------------------------------------------------------

test_that("run_limma_proteomics with batch equals a direct ~ 0 + group + batch fit", {
    skip_if_not_installed("limma")
    fx <- make_batch_fixture()
    res <- run_fixture(fx, batch_cfg())
    ref <- direct_limma(fx$expr, fx$meta)

    expect_equal(colnames(res$design), c("Ctrl", "Trt", "batch_B2"))
    expect_equal(unname(res$contrast_matrix["batch_B2", 1]), 0)
    for (f in c("coefficients", "t", "p.value")) {
        expect_equal(unname(res$fit2[[f]][, "Trt_vs_Ctrl"]), unname(ref[[f]][, 1]), info = f)
    }
    expect_equal(res$fit2$df.total, ref$df.total)
    tbl <- res$de_tables[["Trt_vs_Ctrl"]]
    expect_equal(tbl$FeatureID, rownames(fx$expr))
    expect_equal(tbl$logFC, unname(ref$coefficients[, 1]))
})

test_that("old and new fits differ exactly by the batch term", {
    skip_if_not_installed("limma")
    fx <- make_batch_fixture()
    with_batch <- run_fixture(fx, batch_cfg("covariate"))
    without    <- run_fixture(fx, batch_cfg("none"))
    expect_equal(colnames(without$design), c("Ctrl", "Trt"))
    expect_null(without$batch_design)
    # One more coefficient, one fewer residual df: limma is charged for batch.
    expect_equal(unique(with_batch$fit2$df.residual), unique(without$fit2$df.residual) - 1)
})

# --- the full path: mod_proteomics_de() with multiple imputation ----------

make_de_pre <- function(fx) {
    filt <- fx$expr
    withr::with_seed(3, filt[sample(length(filt), 40)] <- NA)
    # Stand-in for a corrected QC matrix: different from the measured values
    # everywhere, which is all the contract needs. No sva required.
    qc <- fx$expr - 2
    list(expr_raw = filt, expr_filt = filt, expr_filt_pre_imp = filt,
         expr_work = qc, expr_imp_single = qc, meta = fx$meta,
         row_data = fx$prot_tbl, batch_corrected = TRUE, batch_method = "combat")
}

run_de <- function(fx, pre, cfg) {
    suppressMessages(mod_proteomics_de(pre = pre, inputs = list(contrasts = fx$contrasts_df),
                                       config = cfg))
}

test_that("every imputation draw starts from the measured values and every fit has batch", {
    skip_if_not_installed("limma")
    fx  <- make_batch_fixture()
    pre <- make_de_pre(fx)
    de  <- run_de(fx, pre, batch_cfg(de_mode = NULL, method = "combat"))   # default -> covariate

    obs <- !is.na(pre$expr_filt)
    expect_length(de$imputations, 3L)
    for (i in seq_along(de$imputations)) {
        imp <- de$imputations[[i]]
        expect_false(anyNA(imp))
        expect_identical(imp[obs], pre$expr_filt[obs])          # measured, not corrected
        expect_identical(de$runs[[i]]$design, de$runs[[1]]$design)
        expect_equal(colnames(de$runs[[i]]$design), c("Ctrl", "Trt", "batch_B2"))

        ref <- direct_limma(imp, fx$meta)
        for (f in c("coefficients", "t", "p.value")) {
            expect_equal(unname(de$runs[[i]]$fit2[[f]][, "Trt_vs_Ctrl"]),
                         unname(ref[[f]][, 1]), info = sprintf("run %d %s", i, f))
        }
    }
    # Draws differ (perseus_like is stochastic), so the consensus measures something.
    expect_false(identical(de$imputations[[1]], de$imputations[[2]]))
    expect_identical(de$batch_design, de$runs[[1]]$batch_design)
})

test_that("Rubin pooling over batch-adjusted runs equals pooling the direct fits", {
    skip_if_not_installed("limma")
    fx  <- make_batch_fixture()
    pre <- make_de_pre(fx)
    de  <- run_de(fx, pre, batch_cfg(method = "combat"))
    observed <- !is.na(pre$expr_filt)

    pooled <- pool_limma_imputations(de$runs, "Trt_vs_Ctrl", observed = observed)
    direct_runs <- lapply(de$imputations, function(imp) list(fit2 = direct_limma(imp, fx$meta)))
    ref <- pool_limma_imputations(direct_runs, "Trt_vs_Ctrl")

    for (col in c("pooled_logFC", "W", "B", "T", "t_pooled", "df_com", "df_BR")) {
        expect_equal(pooled[[col]], ref[[col]], info = col)
    }
    # The complete-data df is the fits' df.total, which already carries the
    # degree of freedom the batch term costs.
    expect_equal(pooled$df_com,
                 rowMeans(sapply(direct_runs, function(r) r$fit2$df.total)))
    # Observed counts follow the contrast's conditions; batch columns carry
    # weight 0 and add no samples.
    expect_true(all(pooled$n_samples == 12L))
    expect_equal(pooled$n_observed, unname(rowSums(observed)))
})

test_that("legacy_corrected fits on the corrected values with measured missingness", {
    skip_if_not_installed("limma")
    fx  <- make_batch_fixture()
    pre <- make_de_pre(fx)
    de  <- run_de(fx, pre, batch_cfg(de_mode = "legacy_corrected", method = "combat"))

    obs <- !is.na(pre$expr_filt)
    for (imp in de$imputations) {
        # Exactly what correct_batch_proteomics() used to write into expr_filt:
        # the corrected values, re-imputed where the measurement was missing.
        expect_identical(imp[obs], pre$expr_imp_single[obs])
    }
    expect_equal(colnames(de$runs[[1]]$design), c("Ctrl", "Trt"))
    expect_null(de$batch_design)

    no_corr <- pre
    no_corr$batch_corrected <- NULL
    expect_error(run_de(fx, no_corr, batch_cfg(de_mode = "legacy_corrected", method = "combat")),
                 "no batch correction ran")
})

test_that("the resolved de_mode is logged, including when the YAML says null", {
    skip_if_not_installed("limma")
    fx  <- make_batch_fixture()
    pre <- make_de_pre(fx)
    log_for <- function(cfg) {
        testthat::capture_messages(mod_proteomics_de(
            pre = pre, inputs = list(contrasts = fx$contrasts_df), config = cfg))
    }
    m <- log_for(batch_cfg(de_mode = NULL, method = "combat"))
    expect_true(any(grepl(paste0("de_mode 'covariate' (resolved from de_mode: null) -- ",
                                 "batch column 'Batch' in the limma design, fitted on ",
                                 "uncorrected values"), m, fixed = TRUE)))
    m2 <- log_for(batch_cfg(de_mode = "none", method = "combat"))
    expect_true(any(grepl("de_mode 'none' -- batch is not in the DE model", m2, fixed = TRUE)))
})

test_that("de_mode none fits the measured values without a batch term", {
    skip_if_not_installed("limma")
    fx  <- make_batch_fixture()
    pre <- make_de_pre(fx)
    de  <- run_de(fx, pre, batch_cfg(de_mode = "none", method = "combat"))
    obs <- !is.na(pre$expr_filt)
    expect_identical(de$imputations[[1]][obs], pre$expr_filt[obs])
    expect_equal(colnames(de$runs[[1]]$design), c("Ctrl", "Trt"))
})

test_that("design diagnostics are raised once per DE run, not once per imputation", {
    skip_if_not_installed("limma")
    fx <- make_batch_fixture()
    # Unbalanced: B1 holds 6 Ctrl + 1 Trt, B2 holds 5 Trt. Fitted fine, but
    # it carries an empty cell (design note) and a large SE inflation (contrast
    # note): Var goes from 1/6 + 1/6 to 1 + 1/6, about 1.87-fold in SE.
    fx$meta$Group <- rep(c("Ctrl", "Trt"), each = 6)
    fx$meta$Batch <- c(rep("B1", 7), rep("B2", 5))
    pre <- make_de_pre(fx)
    cfg <- batch_cfg(method = "combat")

    w <- testthat::capture_warnings(de <- run_de(fx, pre, cfg))
    bd_w <- grep("^Batch design: ", w, value = TRUE)
    expected <- with(de$batch_design$contrasts, sum(nzchar(Note))) + length(de$batch_design$notes)
    expect_gt(expected, 0)
    expect_length(bd_w, expected)            # not expected * no_repetitions
    expect_match(bd_w, "Trt_vs_Ctrl", fixed = TRUE, all = FALSE)
})

# --- config validation ------------------------------------------------------

test_that("validation: a configured correction defaults DE to the batch covariate", {
    p <- function(...) batch_cfg(...)$modes$proteomics
    expect_true(validate_proteomics_config(p(de_mode = NULL, method = "combat")))
    expect_true(validate_proteomics_config(p(de_mode = "covariate")))   # model-only
    expect_true(validate_proteomics_config(p(de_mode = NULL)))          # nothing configured
})

test_that("validation: unsupported methods fail unless legacy is chosen explicitly", {
    for (m in c("limma_percontrast", "ttest", "welch", "anova")) {
        cfg <- batch_cfg(de_mode = NULL, method = "combat")$modes$proteomics
        cfg$de$method <- m
        expect_error(validate_proteomics_config(cfg), "legacy_corrected", info = m)
        expect_error(validate_proteomics_config(cfg), "defaults to 'covariate'", info = m)

        cfg$batch_correction$de_mode <- "legacy_corrected"
        expect_warning(validate_proteomics_config(cfg), "anti-conservative", info = m)
    }
})

test_that("validation: legacy needs a correction, precomputed ignores the default", {
    p <- function(...) batch_cfg(...)$modes$proteomics
    expect_error(validate_proteomics_config(p(de_mode = "legacy_corrected", method = "none")),
                 "no batch correction")

    pre_default <- p(de_mode = NULL, method = "combat")
    pre_default$de$method <- "ttest"     # irrelevant once tables are precomputed
    pre_default$files <- list(de_table = list(Trt_vs_Ctrl = "de.tsv"))
    expect_true(validate_proteomics_config(pre_default))

    pre_explicit <- p(de_mode = "covariate")
    pre_explicit$files <- list(de_table = list(Trt_vs_Ctrl = "de.tsv"))
    expect_error(validate_proteomics_config(pre_explicit), "precomputed")
})

test_that("validation: covariate mode needs a batch column distinct from block_col", {
    p <- function(...) batch_cfg(...)$modes$proteomics
    expect_error(validate_proteomics_config(p(batch_col = NULL)), "no batch column")
    expect_error(validate_proteomics_config(p(block_col = "Batch")), "both de\\$block_col")
})

test_that("the shipped template resolves to no batch in DE and stays valid", {
    f <- file.path(root_dir, "config", "templates", "proteins_config.yaml")
    skip_if_not(file.exists(f))
    prot <- yaml::read_yaml(f)$modes$proteomics
    bc <- get_proteomics_batch_config(prot)
    expect_equal(bc$de_mode, "none")
    expect_false(bc$de_mode_explicit)     # `de_mode: null` reads as unset
})
