# tests/testthat/test-proteomics-de-covariates.R
#
# Continuous covariates in the pooled limma fit (modes$proteomics$de$covariates):
#   - get_proteomics_de_covariates(), resolve_de_covariates(),
#     build_proteomics_de_design(), contrast_se_inflation(),
#     assess_covariate_design(), run_limma_proteomics()  (05_de_summary.R)
#   - mod_proteomics_de(): every imputation draw carries the covariate,
#     diagnostics raised once                            (02_mod_de.R)
#   - pool_limma_imputations() on covariate-adjusted runs (05d_rubin_pooled_stat.R)
#   - validate_proteomics_config()                       (90_config_validate.R)
#
# Synthetic, seeded data only. Fits are compared with a direct limma fit
# written out by hand, and the SE inflation with its closed form.

# Two conditions, six replicates each. A per-sample covariate `depth` shifts
# every protein along a slope; the first 20 proteins carry a treatment effect.
make_cov_fixture <- function(n_prot = 100L, slope = 0.6, effect = 0.8) {
    samples <- sprintf("S%02d", 1:12)
    group   <- rep(c("Ctrl", "Trt"), each = 6)
    depth   <- c(0.1, -0.4, 0.3, 0.9, -0.2, 0.5, 0.0, 0.7, -0.6, 0.4, 1.1, -0.3)

    expr <- withr::with_seed(11, {
        m <- matrix(rnorm(n_prot * 12, mean = 20, sd = 0.3), nrow = n_prot)
        m <- sweep(m, 2, slope * depth, "+")
        m[1:20, group == "Trt"] <- m[1:20, group == "Trt"] + effect
        m
    })
    rownames(expr) <- paste0("P", seq_len(n_prot))
    colnames(expr) <- samples

    meta <- data.frame(SampleName = samples, Group = group, depth = depth,
                       Line = rep(c("L1", "L2", "L3"), 4), stringsAsFactors = FALSE)
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

cov_cfg <- function(covariates = list("depth"), method = "limma", block_col = NULL,
                    batch_correction = list(method = "none")) {
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
            de         = list(method = method, p_adjust_method = "BH", p_cutoff = 0.05,
                              linear_fc_cutoff = 1.5, block_col = block_col,
                              covariates = covariates),
            imputation = list(method = "perseus_like", width = 0.3, downshift = 1.8,
                              multi_imputation = TRUE, no_repetitions = 3,
                              min_no_passed = 2),
            batch_correction = batch_correction
        ))
    )
}

p_cov <- function(...) cov_cfg(...)$modes$proteomics

direct_cov_limma <- function(expr, meta) {
    X <- cbind(Ctrl = as.numeric(meta$Group == "Ctrl"),
               Trt  = as.numeric(meta$Group == "Trt"),
               cov_depth = meta$depth - mean(meta$depth))
    cm <- limma::makeContrasts(Trt_vs_Ctrl = Trt - Ctrl, levels = X)
    limma::eBayes(limma::contrasts.fit(limma::lmFit(expr, X), cm))
}

cov_design_and_contrasts <- function(group, cov, contrasts) {
    g <- factor(group)
    if (is.null(dim(cov))) cov <- cbind(V1 = cov)
    covariates <- scale(cov, center = TRUE, scale = FALSE)
    d <- build_proteomics_de_design(g, covariates = covariates)
    cm <- limma::makeContrasts(contrasts = unname(contrasts), levels = d)
    colnames(cm) <- names(contrasts)
    list(design = d, group = g, covariates = covariates, cm = cm)
}

# --- resolvers --------------------------------------------------------------

test_that("get_proteomics_de_covariates reads the one key, NULL when unset or empty", {
    expect_null(get_proteomics_de_covariates(p_cov(covariates = NULL)))
    expect_null(get_proteomics_de_covariates(p_cov(covariates = list())))
    expect_identical(get_proteomics_de_covariates(p_cov(covariates = list("a", "b"))),
                     c("a", "b"))
    expect_identical(get_proteomics_de_covariates(p_cov(covariates = "a")), "a")
})

test_that("resolve_de_covariates centres the values and keeps the means", {
    meta <- make_cov_fixture()$meta
    m <- resolve_de_covariates(meta, p_cov())
    expect_equal(colnames(m), "depth")
    expect_equal(unname(m[, "depth"]), meta$depth - mean(meta$depth))
    expect_equal(unname(attr(m, "scaled:center")), mean(meta$depth))
    expect_null(resolve_de_covariates(meta, p_cov(covariates = NULL)))
})

test_that("resolve_de_covariates refuses unusable columns with an actionable message", {
    meta <- make_cov_fixture()$meta
    expect_error(resolve_de_covariates(meta, p_cov(covariates = list("Depth"))),
                 "not in the sample metadata")
    err <- expect_error(resolve_de_covariates(meta, p_cov(covariates = list("Line"))),
                        "not numeric")
    expect_match(conditionMessage(err), "batch_correction$de_mode: 'covariate'", fixed = TRUE)
    expect_match(conditionMessage(err), "batch_col: 'Line'", fixed = TRUE)
    expect_match(conditionMessage(err), "~ 0 + condition + Line", fixed = TRUE)

    na_meta <- meta
    na_meta$depth[3] <- NA
    expect_error(resolve_de_covariates(na_meta, p_cov()), "missing or non-finite values for 1 of 12")
    inf_meta <- meta
    inf_meta$depth[3] <- Inf
    expect_error(resolve_de_covariates(inf_meta, p_cov()), "non-finite")

    const <- meta
    const$depth <- 2
    expect_error(resolve_de_covariates(const, p_cov()), "same value for every sample")

    expect_error(resolve_de_covariates(meta, p_cov(batch_correction = list(
        method = "none", batch_col = "Line", de_mode = "covariate"))),
        "cannot yet be combined with a batch covariate")
})

# --- design -----------------------------------------------------------------

test_that("covariate columns follow the conditions, prefixed and syntactic", {
    g <- factor(c("Ctrl", "Ctrl", "Trt", "Trt", "Ctrl", "Trt"))
    cv <- scale(cbind(depth = c(1, 2, 3, 4, 5, 7), `log ratio` = c(0, 1, 0, 2, 1, 1)),
                center = TRUE, scale = FALSE)
    d <- build_proteomics_de_design(g, covariates = cv)
    expect_equal(colnames(d), c("Ctrl", "Trt", "cov_depth", "cov_log.ratio"))
    expect_equal(unname(d[, "cov_depth"]), unname(cv[, "depth"]))
    expect_null(attr(d, "scaled:center"))
    # No covariates: the condition-only design, as before.
    expect_equal(colnames(build_proteomics_de_design(g)), c("Ctrl", "Trt"))
})

test_that("a condition called like a covariate column keeps its name", {
    g <- factor(c("cov_depth", "cov_depth", "X", "X", "cov_depth", "X"))
    cv <- scale(cbind(depth = c(1, 2, 3, 4, 6, 5)), center = TRUE, scale = FALSE)
    d <- build_proteomics_de_design(g, covariates = cv)
    expect_equal(colnames(d)[1:2], c("cov_depth", "X"))
    expect_equal(colnames(d)[3], "cov_depth.1")
})

# --- assessment -------------------------------------------------------------

test_that("SE inflation matches the closed form for one covariate", {
    skip_if_not_installed("limma")
    # Var(B - A) = s^2 (1/nA + 1/nB + d^2 / Sxx_within): d = 1 (mean B minus
    # mean A), Sxx_within = 2 + 2, 1/nA + 1/nB = 2/3 -> inflation sqrt(1.375).
    x <- cov_design_and_contrasts(rep(c("A", "B"), each = 3), c(1, 2, 3, 2, 3, 4),
                                  c(B_vs_A = "B - A"))
    a <- assess_covariate_design(x$design, x$cm, x$group, x$covariates)
    expect_equal(a$contrasts$SE_inflation, round(sqrt(1.375), 3))
    expect_equal(a$contrasts$Difference, 1)
    expect_equal(a$contrasts$Covariate, "V1")
    expect_equal(a$df_residual, 3)
    expect_equal(unname(contrast_se_inflation(x$design, x$cm, x$group)), sqrt(1.375))
    # R^2 of x on group: between SS 1.5 over total SS 5.5.
    expect_equal(a$contrasts$R2_on_condition, round(1.5 / 5.5, 3))
    expect_identical(a$contrasts$Note, "")
})

test_that("a balanced covariate costs nothing in SE", {
    skip_if_not_installed("limma")
    x <- cov_design_and_contrasts(rep(c("A", "B"), each = 3), c(1, 2, 3, 1, 2, 3),
                                  c(B_vs_A = "B - A"))
    a <- assess_covariate_design(x$design, x$cm, x$group, x$covariates)
    expect_equal(a$contrasts$SE_inflation, 1)
    expect_equal(a$contrasts$Difference, 0)
})

test_that("one side entirely outside the others' range is noted as extrapolation", {
    skip_if_not_installed("limma")
    grp <- rep(c("A", "B", "C"), each = 3)
    cv  <- c(1, 2, 3, 1.5, 2.5, 3.5, 8, 9, 10)          # C above every other sample
    x <- cov_design_and_contrasts(grp, cv, c(C_vs_A = "C - A", B_vs_A = "B - A"))
    a <- assess_covariate_design(x$design, x$cm, x$group, x$covariates)
    note_c <- a$contrasts$Note[a$contrasts$Contrast == "C_vs_A"]
    expect_match(note_c, "Every numerator sample (C) has V1 above every other sample",
                 fixed = TRUE)
    expect_match(note_c, "inflates the standard error", fixed = TRUE)
    expect_false(grepl("Every", a$contrasts$Note[a$contrasts$Contrast == "B_vs_A"]))

    # Below, on the denominator side.
    x2 <- cov_design_and_contrasts(grp, -cv, c(A_vs_C = "A - C"))
    a2 <- assess_covariate_design(x2$design, x2$cm, x2$group, x2$covariates)
    expect_match(a2$contrasts$Note, "Every denominator sample (C) has V1 below", fixed = TRUE)
})

test_that("a covariate fixed by condition is refused as collinear", {
    skip_if_not_installed("limma")
    x <- cov_design_and_contrasts(rep(c("A", "B"), each = 3), c(5, 5, 5, 9, 9, 9),
                                  c(B_vs_A = "B - A"))
    err <- expect_error(assess_covariate_design(x$design, x$cm, x$group, x$covariates),
                        "rank-deficient")
    expect_match(conditionMessage(err), "V1 = 1.000", fixed = TRUE)
    expect_match(conditionMessage(err), "batch_correction$de_mode: 'covariate'", fixed = TRUE)
})

test_that("a design with no residual df is refused", {
    skip_if_not_installed("limma")
    x <- cov_design_and_contrasts(c("A", "A", "B", "B"),
                                  cbind(u = c(1, 2, 4, 3), v = c(2, 7, 1, 3)),
                                  c(B_vs_A = "B - A"))
    expect_error(assess_covariate_design(x$design, x$cm, x$group, x$covariates),
                 "no residual degrees of freedom")
})

test_that("few residual df are a design-level note", {
    skip_if_not_installed("limma")
    x <- cov_design_and_contrasts(c("A", "A", "A", "B", "B"), c(1, 2, 4, 3, 5),
                                  c(B_vs_A = "B - A"))
    a <- assess_covariate_design(x$design, x$cm, x$group, x$covariates)
    expect_equal(a$df_residual, 2)
    expect_match(a$notes, "2 residual degree(s) of freedom", fixed = TRUE)
})

# --- the single fit ---------------------------------------------------------

test_that("run_limma_proteomics with a covariate equals a direct ~ 0 + group + cov fit", {
    skip_if_not_installed("limma")
    fx <- make_cov_fixture()
    res <- run_limma_proteomics(expr_imp = fx$expr, meta = fx$meta,
                                contrasts_df = fx$contrasts_df, prot_tbl = fx$prot_tbl,
                                cfg = cov_cfg())
    ref <- direct_cov_limma(fx$expr, fx$meta)

    expect_equal(colnames(res$design), c("Ctrl", "Trt", "cov_depth"))
    expect_equal(unname(res$contrast_matrix["cov_depth", 1]), 0)
    for (f in c("coefficients", "t", "p.value")) {
        expect_equal(unname(res$fit2[[f]][, "Trt_vs_Ctrl"]), unname(ref[[f]][, 1]), info = f)
    }
    expect_equal(res$fit2$df.total, ref$df.total)
    expect_null(res$batch_design)
    expect_equal(res$covariate_design$df_residual, 9)
    expect_equal(res$covariate_design$contrasts$Covariate, "depth")

    without <- run_limma_proteomics(expr_imp = fx$expr, meta = fx$meta,
                                    contrasts_df = fx$contrasts_df, prot_tbl = fx$prot_tbl,
                                    cfg = cov_cfg(covariates = NULL))
    expect_equal(colnames(without$design), c("Ctrl", "Trt"))
    expect_null(without$covariate_design)
    expect_equal(unique(res$fit2$df.residual), unique(without$fit2$df.residual) - 1)
})

# --- multiple imputation ----------------------------------------------------

make_cov_pre <- function(fx) {
    filt <- fx$expr
    withr::with_seed(5, filt[sample(length(filt), 40)] <- NA)
    list(expr_raw = filt, expr_filt = filt, expr_filt_pre_imp = filt,
         expr_work = filt, expr_imp_single = fx$expr, meta = fx$meta,
         row_data = fx$prot_tbl)
}

test_that("every imputation run carries the covariate, and Rubin pooling needs no change", {
    skip_if_not_installed("limma")
    fx  <- make_cov_fixture()
    pre <- make_cov_pre(fx)
    de  <- suppressMessages(mod_proteomics_de(
        pre = pre, inputs = list(contrasts = fx$contrasts_df), config = cov_cfg()))

    expect_length(de$imputations, 3L)
    for (i in seq_along(de$imputations)) {
        expect_equal(colnames(de$runs[[i]]$design), c("Ctrl", "Trt", "cov_depth"))
        ref <- direct_cov_limma(de$imputations[[i]], fx$meta)
        expect_equal(unname(de$runs[[i]]$fit2$coefficients[, "Trt_vs_Ctrl"]),
                     unname(ref$coefficients[, 1]), info = i)
    }
    expect_false(identical(de$imputations[[1]], de$imputations[[2]]))
    expect_identical(de$covariate_design, de$runs[[1]]$covariate_design)

    observed <- !is.na(pre$expr_filt)
    pooled <- pool_limma_imputations(de$runs, "Trt_vs_Ctrl", observed = observed)
    direct_runs <- lapply(de$imputations,
                          function(imp) list(fit2 = direct_cov_limma(imp, fx$meta)))
    ref <- pool_limma_imputations(direct_runs, "Trt_vs_Ctrl")
    for (col in c("pooled_logFC", "W", "B", "T", "t_pooled", "df_com", "df_BR")) {
        expect_equal(pooled[[col]], ref[[col]], info = col)
    }
    # The covariate column has contrast weight 0, so the observed counts still
    # come from the two conditions' samples only.
    expect_true(all(pooled$n_samples == 12L))
    expect_equal(pooled$n_observed, unname(rowSums(observed)))
})

test_that("covariate notes are warned once per DE run, not once per imputation", {
    skip_if_not_installed("limma")
    fx <- make_cov_fixture()
    fx$meta$depth <- c(1, 2, 3, 2, 1, 2.5, 6, 7, 8, 6.5, 7.5, 9)   # Trt above all Ctrl
    pre <- make_cov_pre(fx)
    w <- testthat::capture_warnings(de <- suppressMessages(mod_proteomics_de(
        pre = pre, inputs = list(contrasts = fx$contrasts_df), config = cov_cfg())))
    cd_w <- grep("^Covariate design: ", w, value = TRUE)
    expected <- sum(nzchar(de$covariate_design$contrasts$Note)) +
        length(de$covariate_design$notes)
    expect_gt(expected, 0)
    expect_length(cd_w, expected)
    expect_match(cd_w, "Contrast 'Trt_vs_Ctrl', covariate 'depth'", fixed = TRUE, all = FALSE)
})

# --- config validation ------------------------------------------------------

test_that("validation: covariates with pooled limma pass, empty list too", {
    expect_true(validate_proteomics_config(p_cov()))
    expect_true(validate_proteomics_config(p_cov(covariates = list("a", "b"))))
    expect_true(validate_proteomics_config(p_cov(covariates = list())))
    expect_true(validate_proteomics_config(p_cov(covariates = NULL)))
})

test_that("validation: covariates must be column names", {
    expect_error(validate_proteomics_config(p_cov(covariates = list(1))),
                 "list of sample-metadata column names")
    expect_error(validate_proteomics_config(p_cov(covariates = TRUE)),
                 "list of sample-metadata column names")
    expect_error(validate_proteomics_config(p_cov(covariates = list(""))),
                 "list of sample-metadata column names")
    expect_error(validate_proteomics_config(p_cov(covariates = list(c("a", "b")))),
                 "list of sample-metadata column names")
    expect_error(validate_proteomics_config(p_cov(covariates = list("a", "a"))),
                 "more than once")
})

test_that("validation: only the pooled limma fit takes covariates", {
    for (m in c("limma_percontrast", "ttest", "welch", "anova")) {
        expect_error(validate_proteomics_config(p_cov(method = m)),
                     sprintf("de\\$method is '%s'", m), info = m)
    }
    pre <- p_cov()
    pre$files <- list(de_table = list(Trt_vs_Ctrl = "de.tsv"))
    expect_error(validate_proteomics_config(pre), "precomputed")
})

test_that("validation: covariates are refused next to a batch covariate", {
    explicit <- p_cov(batch_correction = list(method = "none", batch_col = "Batch",
                                              de_mode = "covariate"))
    expect_error(validate_proteomics_config(explicit),
                 "batch_correction\\$de_mode is 'covariate'")
    dflt <- p_cov(batch_correction = list(method = "combat", batch_col = "Batch"))
    expect_error(validate_proteomics_config(dflt), "defaults to 'covariate'")
    # Batch kept out of the model: allowed.
    none <- p_cov(batch_correction = list(method = "combat", batch_col = "Batch",
                                          de_mode = "none"))
    expect_true(validate_proteomics_config(none))
})

test_that("validation: a covariate cannot also be block_col or the batch column", {
    expect_error(validate_proteomics_config(p_cov(covariates = list("Donor"),
                                                  block_col = "Donor")),
                 "both de\\$block_col and in de\\$covariates")
    legacy <- p_cov(covariates = list("Batch"),
                    batch_correction = list(method = "combat", batch_col = "Batch",
                                            de_mode = "legacy_corrected"))
    expect_error(suppressWarnings(validate_proteomics_config(legacy)),
                 "both the batch column and in de\\$covariates")
})
