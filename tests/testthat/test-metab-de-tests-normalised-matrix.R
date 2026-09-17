# Regression tests: on a variance-scaled metabolomics lane, the DE tests must
# run on the matrix that carries the chosen SAMPLE NORMALISATION.
#
# The branch taken when scaling is auto/pareto/range used to reach for
# pre$expr_log. In metabolomics expr_log is transform_metab(met_filtered$mat) --
# the TRANSFORM ALONE. norm_pqn / norm_total_sum / the median shift run on
# parallel targets that feed met_corrected -> expr_work, so on every such lane
# the p-values were computed on data no sample normalisation had ever touched.
# (Lipidomics never had the defect: its expr_log IS the normalisation re-run
# with scaling forced to "none" -- R/domain/lipidomics/02_preprocess.R.)
#
# Measured on the one shipped lane running that branch (Yossi Tam A03, PQN +
# glog10 + auto): its own expr_tested.tsv is bit-identical to the un-normalised
# matrix (max|diff| = 0), while 164 of the 505 metabolites it called
# significant do not survive PQN.
#
# The correct target is NOT expr_work either. Per-feature variance scaling is
# exactly invariant for t_test / t_test_equal / wilcoxon but not for limma,
# whose eBayes borrows variance across features and sees every autoscaled
# residual variance flattened to 1 -- and it leaves the fold change in SD units.
# So the matrix must be normalised AND unscaled: mod_met_corrected() now builds
# it as mat_pre_scale and the pre-contract carries it as expr_pre_scale.

# ---- fixture -----------------------------------------------------------------
# A 4x sample-loading inflation on half the samples of each group. PQN removes
# it; the transform alone does not, so DE on expr_log sees pure loading noise
# and DE on expr_pre_scale does not. The two matrices therefore disagree on the
# p-value of essentially every feature.
make_norm_fixture <- function(n_feat = 60, seed = 7) {
    set.seed(seed)
    samples <- c(sprintf("A%d", 1:6), sprintf("B%d", 1:6))
    meta <- data.frame(sample_id = samples,
                       Group = rep(c("A", "B"), each = 6),
                       stringsAsFactors = FALSE)
    lin <- matrix(stats::rlnorm(n_feat * length(samples), meanlog = 10, sdlog = 0.25),
                  nrow = n_feat,
                  dimnames = list(sprintf("F%03d", seq_len(n_feat)), samples))
    # a real group effect on the first third of the features
    lin[1:20, meta$Group == "B"] <- lin[1:20, meta$Group == "B"] * 2

    # technical sample loading: alternating samples run 4x hot, in BOTH groups,
    # so it is nuisance variance and nothing else.
    raw <- lin
    hot <- seq(1, length(samples), by = 2)
    raw[, hot] <- raw[, hot] * 4

    list(meta = meta, raw = raw, samples = samples)
}

metab_cfg_scaled <- function(method, scaling = "auto") {
    list(modes = list(metabolomics = list(
        effects       = list(samples = "sample_id", color = "Group"),
        preprocessing = list(transform = "log2", pseudocount = 1,
                             chosen_norm = "pqn", scaling = scaling),
        de            = list(method = method, p_cutoff = 0.05,
                             linear_fc_cutoff = 1)
    )))
}

# A pre-contract shaped exactly like the one pipe_metabolomics() builds.
metab_pre_scaled <- function(f, scaling = "auto") {
    expr_log       <- transform_metab(f$raw, method = "log2", pseudocount = 1)
    expr_pre_scale <- transform_metab(norm_pqn(f$raw), method = "log2",
                                      pseudocount = 1)
    list(expr_raw   = f$raw,
         expr_filt  = f$raw,
         expr_log   = expr_log,
         expr_pre_scale = expr_pre_scale,
         expr_work  = scale_metab(expr_pre_scale, method = scaling),
         meta = f$meta, row_data = NULL,
         info = list(normalization = list(sample_norm = "pqn",
                                          transform   = "log2",
                                          scaling     = scaling)))
}

ctr_tbl <- data.frame(Contrast_name = "B_vs_A", Numerator = "B",
                      Denominator = "A", stringsAsFactors = FALSE)

de_on <- function(mat, cond, method, ctr = "B - A") {
    switch(method,
        limma        = de_limma(mat, cond, ctr),
        t_test       = de_t_test(mat, cond, ctr),
        t_test_equal = de_t_test_equal(mat, cond, ctr),
        wilcoxon     = de_wilcoxon(mat, cond, ctr))
}

# ---- 1. the tested matrix is the NORMALISED one ------------------------------
for (m in c("limma", "t_test", "t_test_equal", "wilcoxon")) {
    test_that(paste0("metabolomics DE [", m,
                     "]: an auto-scaled lane is tested on the NORMALISED ",
                     "pre-scaling matrix, not on the transform alone"), {
        f   <- make_norm_fixture()
        pre <- metab_pre_scaled(f)
        cond <- factor(f$meta$Group)

        res <- suppressMessages(
            run_metabolomics_de(pre, metab_cfg_scaled(m), ctr_tbl))
        got <- res$de_tables[["B_vs_A"]]

        want <- de_on(pre$expr_pre_scale, cond, m)   # normalised, unscaled
        bad  <- de_on(pre$expr_log,       cond, m)   # transform only (the bug)

        expect_equal(got$P.Value[match(want$feature_id, got$feature_id)],
                     want$P.Value)
        # ...and the two really are different analyses, so the assertion above
        # is not satisfiable by accident.
        expect_gt(max(abs(want$P.Value - bad$P.Value), na.rm = TRUE), 0.05)
    })
}

test_that("metabolomics DE: expr_work (scaled) is not the tested matrix either", {
    # Variance scaling is p-value-invariant for the two-group tests but NOT for
    # limma, and it puts the fold change in SD units for all four.
    f    <- make_norm_fixture()
    pre  <- metab_pre_scaled(f)
    cond <- factor(f$meta$Group)

    unscaled <- de_on(pre$expr_pre_scale, cond, "limma")
    scaled   <- de_on(pre$expr_work,      cond, "limma")
    expect_gt(max(abs(unscaled$P.Value - scaled$P.Value), na.rm = TRUE), 0.01)

    res <- suppressMessages(
        run_metabolomics_de(pre, metab_cfg_scaled("limma"), ctr_tbl))
    got <- res$de_tables[["B_vs_A"]]
    expect_equal(got$P.Value[match(unscaled$feature_id, got$feature_id)],
                 unscaled$P.Value)
})

test_that("metabolomics DE: a scaling='none' lane still tests expr_work", {
    # The else branch was already correct -- it must not move.
    f   <- make_norm_fixture()
    pre <- metab_pre_scaled(f, scaling = "none")
    cond <- factor(f$meta$Group)
    res <- suppressMessages(
        run_metabolomics_de(pre, metab_cfg_scaled("limma", scaling = "none"),
                            ctr_tbl))
    got  <- res$de_tables[["B_vs_A"]]
    want <- de_on(pre$expr_work, cond, "limma")
    expect_equal(got$P.Value[match(want$feature_id, got$feature_id)],
                 want$P.Value)
})

test_that("metabolomics DE: a pre-contract without expr_pre_scale warns", {
    # Contracts built outside pipe_metabolomics() have only expr_log on offer.
    # Falling back to it silently is how the defect stayed invisible.
    f   <- make_norm_fixture()
    pre <- metab_pre_scaled(f)
    pre$expr_pre_scale <- NULL
    expect_warning(
        suppressMessages(
            run_metabolomics_de(pre, metab_cfg_scaled("t_test"), ctr_tbl)),
        "expr_pre_scale is missing")
})

# ---- 2. mod_met_corrected builds that matrix ---------------------------------
corrected_fixture <- function(f, config) {
    filtered <- list(mat = f$raw, meta = f$meta, row_data = NULL)
    logged   <- mod_met_log(filtered, config)
    list(
        # no is_QC column on most fixtures -> PQN falls back to all samples
        norm_pqn = suppressWarnings(
            mod_met_normalize_linear(filtered, "pqn", config)),
        logged   = logged
    )
}

test_that("mod_met_corrected: mat_pre_scale is the chosen norm WITHOUT scaling", {
    f   <- make_norm_fixture()
    cfg <- metab_cfg_scaled("limma")
    cfg$modes$metabolomics$preprocessing$drift_correction <- list(enabled = FALSE)
    parts <- corrected_fixture(f, cfg)

    out <- suppressWarnings(suppressMessages(mod_met_corrected(
        norm_tss = NULL, norm_median = NULL, norm_pqn = parts$norm_pqn,
        logged = parts$logged, meta = f$meta,
        out_dir = withr::local_tempdir(), config = cfg)))

    expect_false(is.null(out$mat_pre_scale))
    # it IS the chosen normalisation, unscaled...
    expect_equal(out$mat_pre_scale, parts$norm_pqn$mat)
    # ...and expr_work is that same matrix, scaled -- unchanged behaviour.
    expect_equal(out$mat, scale_metab(parts$norm_pqn$mat, method = "auto"))
    # ...and it is NOT the transform-alone matrix the DE used to test.
    expect_gt(max(abs(out$mat_pre_scale - parts$logged$mat)), 0.1)
})

test_that("mod_met_corrected: mat_pre_scale carries the post-normalisation pool filter", {
    # Latent until a lane uses exclude_after_norm: expr_log keeps the pool
    # columns that pre$meta has already dropped, so the metadata match in
    # run_metabolomics_de() would hand the pools an NA condition.
    f <- make_norm_fixture()
    pools <- c("A1", "B1")
    f$meta$is_QC <- f$meta$sample_id %in% pools

    cfg <- metab_cfg_scaled("limma")
    cfg$modes$metabolomics$preprocessing$drift_correction <- list(enabled = FALSE)
    cfg$modes$metabolomics$sample_filter <- list(
        enabled = TRUE, exclude_after_norm = TRUE,
        rules = list(exclude_qc = TRUE))
    parts <- corrected_fixture(f, cfg)

    out <- suppressWarnings(suppressMessages(mod_met_corrected(
        norm_tss = NULL, norm_median = NULL, norm_pqn = parts$norm_pqn,
        logged = parts$logged, meta = f$meta,
        out_dir = withr::local_tempdir(), config = cfg)))

    expect_false(any(pools %in% colnames(out$mat_pre_scale)))
    expect_equal(colnames(out$mat_pre_scale), colnames(out$mat))
    expect_equal(colnames(out$mat_pre_scale), out$meta$sample_id)
})

test_that("mod_met_corrected: mat_pre_scale carries the LOESS drift correction", {
    # drift_applied has been FALSE on every lane measured, so this half of the
    # contract is unexercised in production -- assert it here instead.
    f <- make_norm_fixture()
    # the QC pools must bracket the whole injection run -- LOESS does not
    # extrapolate, and a sample past the last pool is silently skipped.
    f$meta$is_QC <- f$meta$sample_id %in% c("A1", "A3", "A6", "B3", "B6")
    f$meta$injection_order <- seq_len(nrow(f$meta))

    cfg <- metab_cfg_scaled("limma")
    cfg$modes$metabolomics$preprocessing$drift_correction <-
        list(enabled = TRUE, injection_order_col = "injection_order",
             qc_flag_col = "is_QC")
    parts <- corrected_fixture(f, cfg)

    out <- suppressWarnings(suppressMessages(mod_met_corrected(
        norm_tss = NULL, norm_median = NULL, norm_pqn = parts$norm_pqn,
        logged = parts$logged, meta = f$meta,
        out_dir = withr::local_tempdir(), config = cfg)))

    expect_true(out$info$drift_applied)
    expected <- suppressWarnings(suppressMessages(
        apply_drift_correction(parts$norm_pqn$mat, f$meta,
                               cfg$modes$metabolomics)$mat))
    expect_equal(out$mat_pre_scale, expected)
    # the correction actually moved the data -- otherwise the check is vacuous
    expect_gt(max(abs(out$mat_pre_scale - parts$norm_pqn$mat)), 1e-6)
})

# ---- 3. composing with the QC/blank exclusion the fit already does -----------

test_that("metabolomics DE: expr_pre_scale composes with the QC/blank exclusion in the fit", {
    # run_metabolomics_de() drops QC/blank/pool samples before fitting
    # (filter_to_biological). Routing the test onto expr_pre_scale must not
    # disturb that: expr_pre_scale is column-aligned with pre$meta, so the
    # exclusion lands once, on columns whose condition is known, and the fit
    # sees exactly the biological columns.
    f  <- make_norm_fixture()
    qc <- c("A2", "B5")
    f$meta$Group[f$meta$sample_id %in% qc] <- "QC"

    pre <- metab_pre_scaled(f)
    got <- suppressMessages(
        run_metabolomics_de(pre, metab_cfg_scaled("t_test_equal"),
                            ctr_tbl))$de_tables[["B_vs_A"]]

    # the QC columns really are present in the matrix handed to the DE, so the
    # assertion below is not vacuous
    expect_true(all(qc %in% colnames(pre$expr_pre_scale)))
    expect_true(all(qc %in% pre$meta$sample_id))

    bio  <- setdiff(colnames(pre$expr_pre_scale), qc)
    cond <- factor(f$meta$Group[match(bio, f$meta$sample_id)])
    want <- de_on(pre$expr_pre_scale[, bio, drop = FALSE], cond, "t_test_equal")

    expect_equal(got$P.Value[match(want$feature_id, got$feature_id)], want$P.Value)
    expect_false(any(is.na(got$P.Value)))
})

test_that("metabolomics DE: a lane that already ran the post-normalisation pool filter has no column without a condition", {
    # exclude_after_norm drops the pools from mat_pre_scale AND from meta, so
    # the two stay aligned and the in-fit QC exclusion has nothing left to
    # remove -- the filter is not applied twice. expr_log, the matrix this
    # commit moves the test OFF, is the pre-filter matrix and still carries the
    # pool columns: the metadata match in run_metabolomics_de() would hand them
    # an NA condition.
    f     <- make_norm_fixture()
    pools <- c("A1", "B1")
    keep  <- setdiff(f$meta$sample_id, pools)

    pre <- metab_pre_scaled(f)
    pre$expr_pre_scale <- pre$expr_pre_scale[, keep, drop = FALSE]
    pre$expr_work      <- pre$expr_work[, keep, drop = FALSE]
    pre$meta           <- f$meta[f$meta$sample_id %in% keep, , drop = FALSE]
    expect_true(all(pools %in% colnames(pre$expr_log)))   # and not in the rest
    expect_false(any(pools %in% colnames(pre$expr_pre_scale)))

    got <- suppressMessages(
        run_metabolomics_de(pre, metab_cfg_scaled("limma"),
                            ctr_tbl))$de_tables[["B_vs_A"]]

    want <- de_on(pre$expr_pre_scale, factor(pre$meta$Group), "limma")
    expect_equal(got$P.Value[match(want$feature_id, got$feature_id)], want$P.Value)
    expect_false(any(is.na(got$P.Value)))
})
