# tests/testthat/test-metab-de-matrix-contract.R
#
# The matrix contract for metabolomics differential abundance (#287):
#   - DE and the QEA/ssGSEA enrichment test the normalized, transformed (and,
#     when it ran, drift-corrected) matrix WITHOUT feature scaling;
#   - changing preprocessing.scaling must not change any DE statistic;
#   - nothing falls back to expr_log (the transform alone, no normalization).
#
# The fixtures run the REAL preprocessing chain: the command of each target in
# pipe_metabolomics() (met_filtered ... met_corrected, met_missingness_stats,
# metab_pre) is evaluated against synthetic inputs, so the metab_pre adapter is
# the pipeline's own, not a hand-built list.
#
# These tests fail against each version of this defect seen on main:
#   - Mar-Apr 2026, DE always on expr_log            -> normalization tests
#   - May-Jun 2026, DE on the feature-scaled matrix  -> invariance, limma logFC
#   - Jun-Oct 2026, expr_log under auto/pareto/range -> invariance, normalization,
#     enrichment parity
# including the June regression, where an older copy of 03_differential.R was
# restored over the fix.
#
# All data are synthetic.

suppressPackageStartupMessages(library(targets))

.metab_plan <- local({
  plan <- pipe_metabolomics(chosen_norm = "pqn", skip_outputs = TRUE)
  stats::setNames(plan, vapply(plan, function(t) t$settings$name, character(1)))
})

SCALINGS <- c("none", "center", "auto", "pareto", "range")
BIO_GROUPS <- c("A", "B", "C")

# --- fixture ------------------------------------------------------------------
# 12 biological samples (A, B, C x 4) + 8 pooled QC samples, linear scale.
# Features f1-f5 truly differ (B = 2x A, C = 1.5x A). Every B sample is also
# diluted to half concentration: a global, technical effect that only sample
# normalization removes. On the log-only matrix the null features therefore
# show a B-vs-A fold change near log2(0.5); after PQN they show none.
QC_INJECTIONS <- c(1, 4, 6, 9, 12, 15, 17, 20)

.contract_fixture <- function(with_na = FALSE, with_drift = FALSE) {
  groups <- c(rep(BIO_GROUPS, each = 4), rep("QC", 8))
  ids    <- c(sprintf("S%02d", 1:12), sprintf("QC%d", 1:8))
  n_samp <- length(ids)
  n_feat <- 40
  rnd <- withr::with_seed(287, list(
    base   = exp(stats::rnorm(n_feat, mean = 10, sd = 1)),
    noise  = matrix(exp(stats::rnorm(n_feat * n_samp, 0, 0.1)), nrow = n_feat),
    jitter = exp(stats::rnorm(n_samp, 0, 0.05)),
    bio_order = sample(setdiff(seq_len(n_samp), QC_INJECTIONS))
  ))
  # 8 QC samples, including the first and last injection: loess() cannot
  # extrapolate, so a sample injected outside the QC range would make every
  # feature skip drift correction, and with the 4-QC minimum the quadratic
  # fit has no residual degrees of freedom left (#362). With no missing QC
  # values and intensities far above epsilon, every feature is corrected.
  injection_order <- integer(n_samp)
  injection_order[groups == "QC"] <- QC_INJECTIONS
  injection_order[groups != "QC"] <- rnd$bio_order
  raw <- rnd$base * rnd$noise
  dimnames(raw) <- list(sprintf("f%d", seq_len(n_feat)), ids)
  raw[1:5, groups == "B"] <- raw[1:5, groups == "B"] * 2
  raw[1:5, groups == "C"] <- raw[1:5, groups == "C"] * 1.5
  dilution <- ifelse(groups == "B", 0.5, 1) * rnd$jitter
  raw <- sweep(raw, 2, dilution, FUN = "*")

  meta <- data.frame(sample_id = ids, sample_type = groups,
                     is_QC = groups == "QC", stringsAsFactors = FALSE)
  if (with_drift) {
    # Signal decays along injection order; the QC samples track it.
    meta$injection_order <- injection_order
    raw <- sweep(raw, 2, seq(1, 0.6, length.out = n_samp)[injection_order], FUN = "*")
  }
  if (with_na) raw[cbind(c(8, 15, 22), c(2, 7, 11))] <- NA
  list(raw = raw, meta = meta)
}

.contract_config <- function(scaling, drift_enabled = FALSE, method = "limma") {
  list(modes = list(metabolomics = list(
    preprocessing = list(chosen_norm = "pqn", scaling = scaling,
                         transform = "log2", pseudocount = 1,
                         drift_correction = list(enabled = drift_enabled)),
    effects = list(samples = "sample_id", color = "sample_type"),
    de      = list(method = method, condition_column = "sample_type",
                   p_cutoff = 0.05, linear_fc_cutoff = 1.5)
  )))
}

.eval_target <- function(name, env) {
  tgt <- .metab_plan[[name]]
  if (is.null(tgt)) stop("target not in the plan: ", name)
  value <- NULL
  for (e in as.list(tgt$command$expr)) value <- eval(e, env)
  assign(name, value, envir = env)
}

# Runs the pipeline's own targets from met_filtered to metab_pre.
.run_preprocessing <- function(fx, config) {
  env <- new.env(parent = environment())
  env$config  <- config
  env$met_raw <- list(
    expr_raw   = fx$raw,
    meta       = fx$meta,
    row_data   = data.frame(feature_id = rownames(fx$raw),
                            row.names  = rownames(fx$raw)),
    sample_col = "sample_id"
  )
  env$metab_out_dir           <- tempdir()
  env$met_norm_eigenms        <- NULL
  env$met_norm_eigenms_forced <- NULL
  env$met_norm_bio_factor     <- NULL
  env$met_norm_comparison     <- "normalization_comparison.tsv"
  for (nm in c("met_filtered", "met_log", "met_norm_tss", "met_norm_median",
               "met_norm_pqn", "met_corrected", "met_missingness_stats",
               "metab_pre")) {
    suppressMessages(.eval_target(nm, env))
  }
  env
}

.contrasts <- data.frame(Numerator = c("B", "C"), Denominator = c("A", "B"),
                         Contrast_name = c("B_vs_A", "C_vs_B"),
                         stringsAsFactors = FALSE)

.run_de <- function(env, method = "limma") {
  cfg <- env$config
  cfg$modes$metabolomics$de$method <- method
  # No warnings are suppressed: the fixture is continuous, so the t-tests and
  # Wilcoxon have no ties to warn about. Any warning here is worth seeing.
  suppressMessages(run_metabolomics_de(env$metab_pre, cfg, .contrasts))
}

.bio_ids <- function(fx) fx$meta$sample_id[!fx$meta$is_QC]


# --- the adapter and the matrices it hands downstream --------------------------

test_that("metab_pre carries expr_unscaled, aligned with every other matrix", {
  fx <- .contract_fixture()
  for (sc in SCALINGS) {
    env <- .run_preprocessing(fx, .contract_config(sc))
    pre <- env$metab_pre
    expect_false(is.null(pre$expr_unscaled), info = sc)
    for (m in c("expr_filt", "expr_log", "expr_work", "expr_unscaled")) {
      expect_identical(colnames(pre[[m]]), as.character(pre$meta$sample_id),
                       info = paste(sc, m))
      expect_identical(rownames(pre[[m]]), rownames(pre$expr_work),
                       info = paste(sc, m))
    }
    # Scaling stays on the working matrix, exactly as before #287.
    expect_equal(pre$expr_work, scale_metab(env$met_norm_pqn$mat, method = sc),
                 info = sc)
  }
})

test_that("contrast-specific subsetting keeps expr_unscaled aligned for QEA/ssGSEA", {
  fx  <- .contract_fixture()
  env <- .run_preprocessing(fx, .contract_config("auto"))
  pre <- env$metab_pre

  sub <- .subset_pre_to_groups(pre, "sample_type", c("C", "B"), "sample_id")
  for (m in c("expr_raw", "expr_filt", "expr_log", "expr_work", "expr_unscaled")) {
    expect_identical(colnames(sub[[m]]), as.character(sub$meta$sample_id), info = m)
  }

  dm <- .metab_de_matrix_condition(sub, env$config)
  expect_identical(colnames(dm$mat), as.character(dm$meta$sample_id))
  expect_setequal(as.character(dm$condition), c("B", "C"))
  expect_identical(dm$mat, pre$expr_unscaled[, colnames(dm$mat), drop = FALSE])
})

test_that("exclude_after_norm subsets expr_unscaled with the metadata", {
  fx  <- .contract_fixture()
  cfg <- .contract_config("auto")
  # Post-normalization filter: drop the QC samples and one biological sample.
  cfg$modes$metabolomics$sample_filter <- list(
    enabled = TRUE, exclude_after_norm = TRUE,
    rules = list(exclude_qc = TRUE, exclude_samples = "S03")
  )
  sample_col <- cfg$modes$metabolomics$effects$samples
  env  <- .run_preprocessing(fx, cfg)
  pre  <- env$metab_pre
  kept <- setdiff(.bio_ids(fx), "S03")
  pqn  <- env$met_norm_pqn$mat

  # Only the retained samples, in the metadata's sample-ID order.
  expect_identical(as.character(pre$meta[[sample_col]]), kept)
  expect_identical(colnames(pre$expr_unscaled), kept)
  expect_identical(pre$expr_unscaled, pqn[, kept, drop = FALSE])
  # Still aligned with the working matrix, whose scaling (applied before the
  # filter, as before #287) is unchanged.
  expect_identical(dimnames(pre$expr_work), dimnames(pre$expr_unscaled))
  expect_equal(pre$expr_work, scale_metab(pqn, method = "auto")[, kept, drop = FALSE])

  # What DE and QEA/ssGSEA receive.
  expect_identical(resolve_metab_de_matrix(pre), pre$expr_unscaled)
  dm <- .metab_de_matrix_condition(pre, env$config)
  expect_identical(as.character(dm$meta[[sample_col]]), kept)
  expect_identical(dm$mat, pqn[, kept, drop = FALSE])

  # limma's AveExpr is the row mean of the matrix it fitted, so it shows the
  # fit used exactly the retained samples of the unscaled matrix.
  tab <- .run_de(env, "limma")$de_tables$B_vs_A
  expect_equal(tab$AveExpr,
               unname(rowMeans(pqn[tab$feature_id, kept], na.rm = TRUE)))
})


# --- scaling invariance ---------------------------------------------------------

test_that("feature scaling does not change any DE statistic", {
  fx <- .contract_fixture(with_na = TRUE)
  envs <- lapply(stats::setNames(SCALINGS, SCALINGS),
                 function(sc) .run_preprocessing(fx, .contract_config(sc)))

  # limma: the model's logFC is part of the contract.
  ref <- .run_de(envs$none, "limma")
  for (sc in SCALINGS[-1]) {
    res <- .run_de(envs[[sc]], "limma")
    for (ctr in names(ref$de_tables)) {
      for (col in c("logFC", "statistic", "P.Value", "adj.P.Val")) {
        expect_identical(res$de_tables[[ctr]][[col]], ref$de_tables[[ctr]][[col]],
                         info = paste(sc, ctr, col))
      }
    }
  }

  # t-tests and Wilcoxon: only what they compute on the DE matrix. Their logFC
  # is taken from expr_filt, which no scaling touches, so equality there would
  # say nothing; whether that fold-change definition is right is a separate
  # methodological question and is deliberately NOT asserted here.
  for (method in c("t_test", "t_test_equal", "wilcoxon")) {
    ref_m <- .run_de(envs$none, method)
    for (sc in SCALINGS[-1]) {
      res <- .run_de(envs[[sc]], method)
      for (ctr in names(ref_m$de_tables)) {
        for (col in c("statistic", "P.Value", "adj.P.Val")) {
          expect_identical(res$de_tables[[ctr]][[col]], ref_m$de_tables[[ctr]][[col]],
                           info = paste(method, sc, ctr, col))
        }
      }
    }
  }
})


# --- normalization sensitivity -----------------------------------------------------

test_that("DE tests the normalized matrix, not the log-only one, under every scaling", {
  fx   <- .contract_fixture()
  bio  <- .bio_ids(fx)
  grp  <- fx$meta$sample_type[match(bio, fx$meta$sample_id)]
  null <- sprintf("f%d", 6:40)

  # Fixture check: without normalization the dilution is visible on the nulls.
  env0 <- .run_preprocessing(fx, .contract_config("none"))
  lg   <- env0$met_log$mat[null, bio]
  lo_fc <- rowMeans(lg[, grp == "B"]) - rowMeans(lg[, grp == "A"])
  expect_lt(stats::median(lo_fc), -0.6)

  for (sc in SCALINGS) {
    env <- .run_preprocessing(fx, .contract_config(sc))
    expect_identical(resolve_metab_de_matrix(env$metab_pre), env$met_norm_pqn$mat,
                     info = sc)
    tab <- .run_de(env, "limma")$de_tables$B_vs_A
    fc  <- tab$logFC[match(null, tab$feature_id)]
    expect_lt(abs(stats::median(fc)), 0.2, label = paste("median null logFC,", sc))
  }
})


# --- limma logFC is the model contrast on the DE input matrix ----------------------

test_that("limma logFC equals the contrast from an independent per-feature lm()", {
  fx  <- .contract_fixture(with_na = TRUE)
  env <- .run_preprocessing(fx, .contract_config("pareto"))
  res <- .run_de(env, "limma")

  bio  <- .bio_ids(fx)
  cond <- factor(fx$meta$sample_type[match(bio, fx$meta$sample_id)],
                 levels = BIO_GROUPS)
  # The expected DE input, taken from the PQN target directly rather than
  # through the resolver under test. QC samples are not part of the fit.
  y <- env$met_norm_pqn$mat[, bio]

  contrast_coef <- function(num, den) {
    apply(y, 1, function(v) {
      cf <- stats::coef(stats::lm(v ~ 0 + cond))
      unname(cf[paste0("cond", num)] - cf[paste0("cond", den)])
    })
  }
  for (i in seq_len(nrow(.contrasts))) {
    ctr <- .contrasts$Contrast_name[i]
    tab <- res$de_tables[[ctr]]
    expected <- contrast_coef(.contrasts$Numerator[i], .contrasts$Denominator[i])
    expect_equal(tab$logFC[match(names(expected), tab$feature_id)],
                 unname(expected), tolerance = 1e-8, info = ctr)
  }
})


# --- drift correction, as far as the DE matrix is concerned -----------------------
# Whether the LOESS correction itself is right is #362; here only which matrix
# DE receives.

test_that("with explicit drift correction and no scaling, DE tests the corrected matrix", {
  fx  <- .contract_fixture(with_drift = TRUE)
  env <- .run_preprocessing(fx, .contract_config("none", drift_enabled = TRUE))
  pre <- env$metab_pre
  expect_true(pre$info$drift$drift_applied)

  # The fixture guarantees LOESS corrects every feature (see .contract_fixture),
  # so a skipped or invalid correction cannot pass for a real one.
  dinfo <- pre$info$drift$drift_info
  expect_identical(dinfo$n_corrected, nrow(env$met_norm_pqn$mat))
  expect_identical(dinfo$n_skipped, 0L)
  expect_true(all(is.finite(pre$expr_work)))
  expect_true(all(rowSums(pre$expr_work != env$met_norm_pqn$mat) > 0))

  expect_identical(resolve_metab_de_matrix(pre), pre$expr_work)
})

test_that("without drift correction, scaled runs test the uncorrected normalized matrix", {
  fx  <- .contract_fixture(with_drift = TRUE)  # metadata meets the preconditions
  env <- .run_preprocessing(fx, .contract_config("auto", drift_enabled = FALSE))
  pre <- env$metab_pre
  expect_false(pre$info$drift$drift_applied)
  expect_identical(resolve_metab_de_matrix(pre), env$met_norm_pqn$mat)
})


# --- DE and enrichment read the same matrix ----------------------------------------

test_that("QEA/ssGSEA start from the matrix DE tests, under every scaling", {
  fx  <- .contract_fixture()
  bio <- .bio_ids(fx)
  ref <- NULL
  for (sc in SCALINGS) {
    env <- .run_preprocessing(fx, .contract_config(sc))
    dm  <- .metab_de_matrix_condition(env$metab_pre, env$config)
    expect_identical(colnames(dm$mat), bio, info = sc)
    expect_identical(dm$mat,
                     resolve_metab_de_matrix(env$metab_pre)[, bio, drop = FALSE],
                     info = sc)
    if (is.null(ref)) ref <- dm$mat else expect_identical(dm$mat, ref, info = sc)
  }
})


# --- workflow: chosen_norm still decides review vs analysis mode ---------------------

test_that("chosen_norm: null is QC-review only; a method runs the analysis", {
  target_names <- function(plan) vapply(plan, function(t) t$settings$name, character(1))

  review <- target_names(pipe_metabolomics(chosen_norm = NULL))
  expect_true("met_norm_comparison" %in% review)
  expect_false(any(c("met_corrected", "metab_pre", "metab_de_res") %in% review))

  analysis <- names(.metab_plan)
  expect_true(all(c("met_corrected", "metab_pre", "metab_de_res") %in% analysis))
})
