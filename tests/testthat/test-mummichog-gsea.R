# tests/testthat/test-mummichog-gsea.R
#
# Unit tests for the MetaboAnalyst-style MS peaks-to-pathways GSEA layer (06g):
# the ranking-metric choice, the EmpiricalCompound score construction (signed
# max over member features, with same-m/z features merged first), the
# pathway -> detected-EC gene sets, the MetaboAnalyst-equivalent engine, the
# NES interpretation, and the summary scatter's encodings.
#
# Parity against the pinned upstream implementation lives in
# test-mummichog-gsea-parity.R; this file covers the surrounding contract.
#
# The EC-score semantics asserted here are MetaboAnalystR's, verified in
# peaks_to_function.R (`ec.exp.vec <- unlist(lapply(ec_exp_dict, max))` in the
# mumRT/v2 branch) — deliberately NOT mean or max(abs()).
#
# All fixtures are synthetic.

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Per-(EC, feature) frame in the shape read_mummichog_ec_features() returns.
ec_features_fixture <- function() {
  data.frame(
    stringsAsFactors = FALSE,
    EID        = c("E1", "E1", "E2", "E2", "E3"),
    input_row  = paste0("row", 1:5),
    feature_id = paste0("f", 1:5),
    mz         = c(100, 200, 300, 300, 400),
    retention_time = c(1, 1, 2, 2, 3),
    p_value    = c(0.01, 0.02, 0.03, 0.04, 0.20),
    statistic  = c(3, -5, 1, 3, -2),
    adduct     = rep("M+H[1+]", 5)
  )
}

feature_scores_fixture <- function() {
  stats::setNames(c(3, -5, 1, 3, -2), paste0("f", 1:5))
}

# A GSEA result table in the shape run_mummichog_gsea()$table returns.
gsea_table_fixture <- function() {
  data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    "Pathway"                         = c("PW_pos", "PW_neg", "PW_mid"),
    "Pathway size (model compounds)"  = c(8, 8, 8),
    "Detected ECs"                    = c(8, 8, 8),
    "Tested size"                     = c(8, 8, 8),
    "ES"                              = c(0.80, -0.75, 0.10),
    "NES"                             = c(2.10, -1.90, 0.00),
    "P.Value"                         = c(0.001, 0.004, 0.700),
    "padj"                            = c(0.003, 0.006, 0.700),
    "Leading-edge EmpiricalCompounds" = c("E001; E002", "E024; E023", NA)
  )
}


# ---------------------------------------------------------------------------
# ranking statistic
# ---------------------------------------------------------------------------

test_that("a signed-statistic DE method ranks on the t statistic", {
  de <- data.frame(feature_id = c("f1", "f2"), logFC = c(1, -1),
                   statistic = c(4, -6), stringsAsFactors = FALSE)

  for (m in c("limma", "t_test", "t_test_equal")) {
    r <- mmc_gsea_rank_metric(de, de_method = m)
    expect_identical(r$metric, "t_statistic")
    expect_equal(unname(r$values[c("f1", "f2")]), c(4, -6))
    expect_match(r$label, m, fixed = TRUE)
  }
})

test_that("a signed-statistic method falls back to logFC without a statistic", {
  de <- data.frame(feature_id = c("f1", "f2"), logFC = c(1, -1),
                   stringsAsFactors = FALSE)
  r <- mmc_gsea_rank_metric(de, de_method = "limma")

  expect_identical(r$metric, "logFC")
  expect_equal(unname(r$values[c("f1", "f2")]), c(1, -1))
  expect_match(r$label, "no usable t statistic")

  # an all-NA statistic column is not usable either
  de2 <- transform(de, statistic = c(NA_real_, NA_real_))
  expect_identical(mmc_gsea_rank_metric(de2, de_method = "limma")$metric,
                   "logFC")
})

test_that("the GSEA ranking metric is NULL when nothing is usable", {
  expect_null(mmc_gsea_rank_metric(NULL))
  expect_null(mmc_gsea_rank_metric(data.frame()))
  expect_null(mmc_gsea_rank_metric(data.frame(feature_id = "f1",
                                              AveExpr = 1)))
})

test_that("a literal 't' column is accepted as the signed statistic", {
  de <- data.frame(feature_id = c("f1", "f2"), logFC = c(1, -1),
                   t = c(2, -3), stringsAsFactors = FALSE)
  r <- mmc_gsea_rank_metric(de, de_method = "limma")
  expect_identical(r$metric, "t_statistic")
  expect_equal(unname(r$values[["f2"]]), -3)
})

test_that("an unusable 'statistic' column falls through to 't', not to logFC", {
  # A DE table can carry an all-NA `statistic` next to a perfectly good `t`.
  # Picking the column by name alone committed to `statistic` before testing
  # it, so the real statistic was skipped and the ranking dropped to logFC.
  de <- data.frame(stringsAsFactors = FALSE,
                   feature_id = c("f1", "f2"),
                   statistic  = c(NA_real_, NA_real_),
                   t          = c(2, -3),
                   logFC      = c(1, -1))
  r <- mmc_gsea_rank_metric(de, de_method = "limma")
  expect_identical(r$metric, "t_statistic")
  expect_match(r$label, "'t'", fixed = TRUE)
  expect_equal(unname(r$values[["f2"]]), -3)

  # With no usable statistic at all, logFC is still the documented fallback.
  de_no_t <- de
  de_no_t$t <- NA_real_
  expect_identical(mmc_gsea_rank_metric(de_no_t, de_method = "limma")$metric,
                   "logFC")
})


# ---- the ranking statistic is chosen by DE method, never by column name ----

test_that("Wilcoxon never ranks on, or is labelled as, a t statistic", {
  # de_wilcoxon() stores wilcox.test()'s W in the same generic `statistic`
  # column limma uses for its moderated t. W is non-negative and scales with
  # group size, so ranking on it orders EmpiricalCompounds by rank-sum
  # magnitude with no direction -- while the report would claim a t statistic.
  de <- data.frame(stringsAsFactors = FALSE,
                   feature_id = c("f1", "f2", "f3"),
                   logFC      = c(2, -1, 0.5),
                   statistic  = c(9, 7, 8))       # W: never negative

  r <- mmc_gsea_rank_metric(de, de_method = "wilcoxon")

  expect_identical(r$metric, "logFC")
  expect_equal(unname(r$values[c("f1", "f2", "f3")]), c(2, -1, 0.5))
  expect_false(any(grepl("t statistic", r$label, fixed = TRUE)))
  expect_match(r$label, "wilcoxon", fixed = TRUE)
  expect_match(r$label, "unsigned")

  # A literal `t` column does not re-open the door either: the method, not the
  # column name, decides.
  expect_identical(
    mmc_gsea_rank_metric(transform(de, t = c(3, -2, 1)),
                         de_method = "wilcoxon")$metric,
    "logFC")
})

test_that("Wilcoxon without logFC is refused, not ranked on W", {
  # There is no second-choice signed metric here, and W is not one. Refusing is
  # the only honest outcome -- silently ranking on W is the defect itself.
  de <- data.frame(feature_id = c("f1", "f2"), statistic = c(9, 7),
                   stringsAsFactors = FALSE)
  expect_null(mmc_gsea_rank_metric(de, de_method = "wilcoxon"))
})

test_that("an unknown or missing DE method does not trust 'statistic'", {
  # `statistic` means different things per method, so without knowing which one
  # produced the table there is no basis for calling it signed.
  de <- data.frame(feature_id = c("f1", "f2"), logFC = c(1, -1),
                   statistic = c(4, -6), stringsAsFactors = FALSE)

  r_na <- mmc_gsea_rank_metric(de)                        # method not recorded
  expect_identical(r_na$metric, "logFC")
  expect_match(r_na$label, "was not recorded")

  r_unknown <- mmc_gsea_rank_metric(de, de_method = "precomputed")
  expect_identical(r_unknown$metric, "logFC")
  expect_match(r_unknown$label, "precomputed", fixed = TRUE)

  # ... and with no logFC to fall back on, the contrast is refused rather than
  # ranked on a column of unknown provenance.
  expect_null(mmc_gsea_rank_metric(de[, c("feature_id", "statistic")],
                                   de_method = "precomputed"))
})

test_that("the chosen metric and the method that chose it are reported back", {
  de <- data.frame(feature_id = c("f1", "f2"), logFC = c(1, -1),
                   statistic = c(4, -6), stringsAsFactors = FALSE)
  expect_identical(mmc_gsea_rank_metric(de, de_method = "limma")$de_method,
                   "limma")
  expect_true(is.na(mmc_gsea_rank_metric(de)$de_method))
})


# ---- gseaParam: the engine runs only where the reference is self-consistent -

test_that("the engine refuses any gsea_param but the reference value", {
  # The pinned reference raises the observed score's weights to the exponent a
  # second time (its stats vector is already abs()^gseaParam) while hardcoding
  # gseaParam = 1 in every permutation path. Off 1, observed and permuted
  # scores are weighted differently, so NES and the permutation p-value are not
  # comparable. Reproducing that silently would publish invalid statistics;
  # "fixing" it would stop matching the reference.
  pw <- list(PW = c("E1", "E2"))
  st <- stats::setNames(c(2, 1, -1), c("E1", "E2", "E3"))
  rk <- stats::setNames(c(1, 2, 3), c("E1", "E2", "E3"))

  expect_error(mmc_gsea_metaboanalyst(pw, st, rk, gsea_param = 2),
               "gsea_param must be 1")
  expect_error(mmc_gsea_metaboanalyst(pw, st, rk, gsea_param = 0.5),
               "gsea_param must be 1")
  # The message has to say why, not just refuse.
  expect_error(mmc_gsea_metaboanalyst(pw, st, rk, gsea_param = 2),
               "permutation")
})

test_that("run_mummichog_gsea rejects gsea_param early, before any work", {
  # Early and loud: inside the run it would be caught by the engine's tryCatch
  # and degrade into a warning plus a NULL result, which reads like "GSEA was
  # not applicable" rather than "you asked for something invalid".
  expect_error(
    run_mummichog_gsea(character(0), data.frame(feature_id = "f", logFC = 1),
                       model = list(pathways = list(P = "C1")),
                       de_method = "limma", gsea_param = 2),
    "gsea_param must be 1")
})

test_that("an ES = 0 fallback propagates into NES, the p-value and the BH set", {
  for (pkg in c("fgsea", "fastmatch", "data.table", "BiocParallel")) {
    skip_if_not_installed(pkg)
  }
  # E2 and E3 share a score, so they collapse into ONE ranked position; a
  # pathway holding both hands calcGseaStat a duplicated index and takes the
  # reference implementation's ES = 0 fallback. What that zero then does to the
  # rest of the row is the point of this test -- it is reproduced, not fixed,
  # so it has to be pinned down rather than left to be rediscovered.
  ec  <- c(E1 = 3, E2 = 1, E3 = 1, E4 = -2, E5 = -4)
  inp <- mmc_gsea_ranked_inputs(ec)
  res <- mmc_gsea_metaboanalyst(
    list(PW_tie = c("E2", "E3"), PW_ok = c("E1", "E4")),
    inp$stats, inp$ranks)

  tie <- res[res$pathway == "PW_tie", ]
  ok  <- res[res$pathway == "PW_ok", ]
  expect_equal(nrow(tie), 1L)
  expect_equal(nrow(ok), 1L)

  expect_true(tie$es_defaulted)
  expect_match(tie$es_reason, "^tied EC scores")
  expect_identical(tie$ES, 0)

  # p = 1 is exact, not approximate: with ES = 0 every permuted score is both
  # <= and >= it, so leEs == leZero and geEs == geZero and both ratios are 1.
  expect_identical(tie$pval, 1)

  # NES is 0 / (mean permuted score of that sign) -- zero whenever that mean is
  # finite and non-zero, NaN in the degenerate case where no permuted score
  # falls on that side. Either way it is not a measured enrichment; what must
  # never happen is a defaulted ES yielding a real-looking NES.
  expect_true(tie$NES == 0 || is.nan(tie$NES))

  # The pathway that did score is unaffected and keeps a genuine NES.
  expect_false(ok$es_defaulted)
  expect_true(is.na(ok$es_reason))
  expect_true(is.finite(ok$ES))

  # The defaulted row is NOT dropped from the multiple-testing set: BH runs over
  # every tested pathway, so a fallback still consumes a test and raises the
  # adjusted p-values of the pathways that did score.
  expect_equal(res$padj, stats::p.adjust(res$pval, "fdr"))
  expect_equal(nrow(res), 2L)
})

# ---- provenance: what ran must be what the result says ran ----------------

test_that("an empty engine result carries the same columns as a full one", {
  for (pkg in c("fgsea", "fastmatch", "data.table", "BiocParallel")) {
    skip_if_not_installed(pkg)
  }
  # Callers rbind runs together and read the documented schema; a valid empty
  # run must not hand back a narrower contract than a populated one.
  st <- stats::setNames(c(2, 1, -1), c("E1", "E2", "E3"))
  rk <- stats::setNames(c(1, 2, 3), c("E1", "E2", "E3"))

  no_pathways <- mmc_gsea_metaboanalyst(list(), st, rk)
  filtered    <- mmc_gsea_metaboanalyst(list(PW = c("E1", "E2")), st, rk,
                                        min_size = 99)
  full        <- mmc_gsea_metaboanalyst(list(PW = c("E1", "E2")), st, rk)

  expect_equal(nrow(no_pathways), 0L)
  expect_equal(nrow(filtered), 0L)
  expect_gt(nrow(full), 0L)

  expect_identical(names(no_pathways), names(full))
  expect_identical(names(filtered), names(full))
  expect_true(all(c("es_defaulted", "es_reason") %in% names(no_pathways)))
  expect_type(no_pathways$es_defaulted, "logical")
  expect_type(no_pathways$es_reason, "character")
})

test_that("min_size is normalised once and reported as the value that ran", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  # The engine clamps min_size to 1. Reporting the caller's 0 would describe a
  # run that never happened, and make the result unreproducible from its own
  # provenance.
  res <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                            de_method = "limma", min_size = 0)

  expect_equal(res$params$min_size, 1)
  # ... and 0 is not a departure from the reference either: 1 is the reference.
  expect_false(any(grepl("min_size", res$deviations)))

  same <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                             de_method = "limma", min_size = 1)
  expect_equal(res$table, same$table)
})

test_that("a non-default rank_metric is recorded as a deviation", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  # Merging features that share an m/z by something other than the reference's
  # "mean" changes every EC score, and so every statistic downstream of it.
  res <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                            de_method = "limma", rank_metric = "median")
  expect_true(any(grepl("rank_metric = 'median'", res$deviations, fixed = TRUE)))

  ref <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                            de_method = "limma")
  expect_false(any(grepl("rank_metric", ref$deviations)))
  expect_identical(ref$params$rank_metric, "mean")
})

test_that("an abbreviated rank_metric resolves before it is judged or reported", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  # match.arg() accepts abbreviations, and mmc_ec_scores() resolves them in its
  # own frame only. Left unresolved at the entry point, "mea" would run the
  # reference "mean" aggregation while being announced as a methodological
  # deviation that never happened, and reported under a name that is not a
  # valid setting.
  abbrev <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                               de_method = "limma", rank_metric = "mea")

  expect_identical(abbrev$params$rank_metric, "mean")
  expect_false(any(grepl("rank_metric", abbrev$deviations)))

  full <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                             de_method = "limma", rank_metric = "mean")
  expect_equal(abbrev$table, full$table)

  # A non-default abbreviation still resolves, and is still a deviation -- the
  # normalisation must not swallow the departure along with the shorthand.
  med <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                            de_method = "limma", rank_metric = "medi")
  expect_identical(med$params$rank_metric, "median")
  expect_true(any(grepl("rank_metric = 'median'", med$deviations, fixed = TRUE)))

  # Both entry points answer to the same list of allowed values.
  expect_identical(.MMC_GSEA_RANK_METRICS[[1]], .MMC_GSEA_DEFAULTS$rank_metric)
  expect_error(run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                                  de_method = "limma", rank_metric = "nope"))
})

test_that("the single-pathway branch is reproducible and declared unvalidated", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)
  # One pathway -> the engine takes its seeded sample.int() loop instead of the
  # batch call. Upstream leaves that loop to a non-deterministic worker RNG, so
  # there is no reference draw to match: we offer reproducibility instead, and
  # the result has to say so rather than imply parity.
  one <- list(pathways = f$model$pathways["PW_pos"])

  a <- run_mummichog_gsea(f$files, f$de_table, one, contrast = "c",
                          de_method = "limma")
  b <- run_mummichog_gsea(f$files, f$de_table, one, contrast = "c",
                          de_method = "limma")

  expect_equal(nrow(a$table), 1L)
  expect_equal(a$table, b$table)          # seeded: identical across runs

  expect_true(any(grepl("one tested pathway", a$deviations)))
  expect_true(any(grepl("not covered by the parity fixture", a$deviations)))

  # More than one tested pathway goes through the batch call and claims parity.
  many <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                             de_method = "limma")
  expect_gt(nrow(many$table), 1L)
  expect_false(any(grepl("one tested pathway", many$deviations)))
})

test_that("the fgsea compatibility policy is decided in one place", {
  skip_if_not_installed("fgsea")
  installed <- as.character(utils::packageVersion("fgsea"))
  series <- function(v) paste(strsplit(v, ".", fixed = TRUE)[[1]][1:2],
                              collapse = ".")

  if (identical(installed, .MMC_GSEA_VALIDATED_FGSEA)) {
    # The validated build: supported, and silent about it.
    expect_null(.mmc_gsea_check_fgsea())
  } else if (identical(series(installed), series(.MMC_GSEA_VALIDATED_FGSEA))) {
    # Same series, different patch: allowed, warned about once, and the exact
    # version recorded so the result says what it actually ran on.
    expect_warning(note <- .mmc_gsea_check_fgsea(), installed, fixed = TRUE)
    expect_match(note, installed, fixed = TRUE)
    expect_match(note, .MMC_GSEA_VALIDATED_FGSEA, fixed = TRUE)
    # The provenance path must not warn a second time for the same run.
    expect_silent(.mmc_gsea_check_fgsea(warn = FALSE))
  } else {
    # A different series is refused outright -- no statistics are computed on a
    # build whose internals the parity suite declines to validate.
    expect_error(.mmc_gsea_check_fgsea(), "validated compatibility range")
    expect_error(mmc_gsea_metaboanalyst(
      list(PW = c("E1", "E2")),
      stats::setNames(c(2, 1), c("E1", "E2")),
      stats::setNames(c(1, 2), c("E1", "E2"))),
      "validated compatibility range")
  }
})

test_that("fractional n_perm and seed are refused, not truncated", {
  pw <- list(PW = c("E1", "E2"))
  st <- stats::setNames(c(2, 1, -1), c("E1", "E2", "E3"))
  rk <- stats::setNames(c(1, 2, 3), c("E1", "E2", "E3"))

  # Truncating would put the engine back where min_size and rank_metric were:
  # running 100 permutations while reporting 100.5, or running the default
  # seed 123 while announcing a departure from it.
  expect_error(mmc_gsea_metaboanalyst(pw, st, rk, n_perm = 100.5),
               "n_perm must be a single whole number")
  expect_error(mmc_gsea_metaboanalyst(pw, st, rk, seed = 123.9),
               "seed must be a single whole number")
  expect_error(mmc_gsea_metaboanalyst(pw, st, rk, n_perm = 0), "n_perm")
  expect_error(mmc_gsea_metaboanalyst(pw, st, rk, n_perm = c(10, 20)), "n_perm")
  expect_error(mmc_gsea_metaboanalyst(pw, st, rk, n_perm = NA), "n_perm")

  # The entry point refuses before any work, so the caller gets the error
  # rather than "GSEA skipped".
  de <- data.frame(feature_id = "f", logFC = 1, stringsAsFactors = FALSE)
  model <- list(pathways = list(P = "C1"))
  expect_error(run_mummichog_gsea(character(0), de, model, de_method = "limma",
                                  n_perm = 100.5),
               "n_perm must be a single whole number")
  expect_error(run_mummichog_gsea(character(0), de, model, de_method = "limma",
                                  seed = 123.9),
               "seed must be a single whole number")

  # Whole numbers in either storage mode are fine, and 0 is a legal seed.
  expect_identical(.mmc_gsea_check_count(100, "n_perm"), 100L)
  expect_identical(.mmc_gsea_check_count(100L, "n_perm"), 100L)
  expect_identical(.mmc_gsea_check_count(0, "seed", min = 0L), 0L)
})

test_that("a pathway covering the whole ranking still returns its row", {
  for (pkg in c("fgsea", "fastmatch", "data.table", "BiocParallel")) {
    skip_if_not_installed(pkg)
  }
  # One pathway, holding every ranked position. calcGseaStat refuses that set,
  # so the observed score correctly takes the ES = 0 fallback -- but with one
  # tested pathway the permutation loop draws the whole ranking every time and
  # hit the same refusal unguarded, aborting the engine. The documented
  # ES-defaulted row was lost and the contrast read as "GSEA not applicable".
  ec  <- c(E1 = 3, E2 = 2, E3 = 1)
  inp <- mmc_gsea_ranked_inputs(ec)
  res <- mmc_gsea_metaboanalyst(list(PW_all = c("E1", "E2", "E3")),
                                inp$stats, inp$ranks)

  expect_equal(nrow(res), 1L)
  expect_true(res$es_defaulted)
  expect_identical(res$es_reason, "pathway selects every ranked position")
  expect_identical(res$ES, 0)
  expect_identical(res$pval, 1)
  # Classified by its actual cause, not filed under ties.
  expect_false(grepl("tie", res$es_reason, ignore.case = TRUE))
})

test_that("the reference gsea_param is accepted in both numeric shapes", {
  for (pkg in c("fgsea", "fastmatch", "data.table", "BiocParallel")) {
    skip_if_not_installed(pkg)
  }
  pw <- list(PW = c("E1", "E2"))
  st <- stats::setNames(c(2, 1, -1), c("E1", "E2", "E3"))
  rk <- stats::setNames(c(1, 2, 3), c("E1", "E2", "E3"))

  # 1L must pass as readily as 1: the gate compares values, not storage modes.
  expect_s3_class(mmc_gsea_metaboanalyst(pw, st, rk, gsea_param = 1),
                  "data.frame")
  expect_s3_class(mmc_gsea_metaboanalyst(pw, st, rk, gsea_param = 1L),
                  "data.frame")
})


# ---------------------------------------------------------------------------
# EmpiricalCompound scores (MetaboAnalyst v2 semantics)
# ---------------------------------------------------------------------------

test_that("an EC score is the SIGNED maximum of its member features", {
  ec <- mmc_ec_scores(ec_features_fixture(), feature_scores_fixture())

  # E1 = max(3, -5) = 3 — not the mean (-1) and not max(abs()) (5 / -5)
  expect_equal(unname(ec[["E1"]]), 3)
  expect_false(isTRUE(all.equal(unname(ec[["E1"]]), -1)))
  expect_false(isTRUE(all.equal(unname(ec[["E1"]]), 5)))
  # E3 keeps its negative score (no absolute value anywhere)
  expect_equal(unname(ec[["E3"]]), -2)
})

test_that("features sharing one m/z are merged before the EC max", {
  # f3 and f4 both sit at m/z 300 with scores 1 and 3 -> merged to mean 2,
  # so E2 scores 2 rather than 3.
  ec <- mmc_ec_scores(ec_features_fixture(), feature_scores_fixture())
  expect_equal(unname(ec[["E2"]]), 2)

  # the merge rule is selectable, matching MetaboAnalyst's rank.metric
  expect_equal(unname(mmc_ec_scores(ec_features_fixture(),
                                    feature_scores_fixture(),
                                    rank_metric = "max")[["E2"]]), 3)
  expect_equal(unname(mmc_ec_scores(ec_features_fixture(),
                                    feature_scores_fixture(),
                                    rank_metric = "min")[["E2"]]), 1)
})

test_that("the EC score vector is sorted decreasing and skips unscored features", {
  ec <- mmc_ec_scores(ec_features_fixture(), feature_scores_fixture())
  expect_identical(names(ec), c("E1", "E2", "E3"))
  expect_true(all(diff(unname(ec)) <= 0))

  # a feature with no DE score simply does not contribute
  partial <- feature_scores_fixture()[c("f1", "f5")]
  ec2 <- mmc_ec_scores(ec_features_fixture(), partial)
  expect_setequal(names(ec2), c("E1", "E3"))
  expect_equal(unname(ec2[["E1"]]), 3)
})

test_that("EC scores are empty when there is nothing to score", {
  expect_length(mmc_ec_scores(NULL, feature_scores_fixture()), 0)
  expect_length(mmc_ec_scores(ec_features_fixture(),
                              stats::setNames(numeric(0), character(0))), 0)
  expect_length(mmc_ec_scores(ec_features_fixture()[0, ],
                              feature_scores_fixture()), 0)
})


# ---------------------------------------------------------------------------
# pathway -> detected EC sets
# ---------------------------------------------------------------------------

test_that("pathway EC sets collect every EC with a candidate in the pathway", {
  cand <- data.frame(
    stringsAsFactors = FALSE,
    EID         = c("E1", "E1", "E2", "E3"),
    compound_id = c("C1", "C2", "C2", "C9"),
    compound_name = c("a", "b", "b", "i")
  )
  sets <- mmc_pathway_ec_sets(list(P1 = c("C1", "C3"), P2 = c("C2"),
                                   P3 = c("C7")),
                              cand, ec_ids = c("E1", "E2", "E3"))

  expect_setequal(names(sets), c("P1", "P2"))       # P3 has no detected EC
  expect_setequal(sets$P1, "E1")
  expect_setequal(sets$P2, c("E1", "E2"))           # both ECs carry C2
})

test_that("pathway EC sets are restricted to the ranked universe", {
  cand <- data.frame(EID = c("E1", "E2"), compound_id = c("C1", "C1"),
                     compound_name = c("a", "a"), stringsAsFactors = FALSE)
  sets <- mmc_pathway_ec_sets(list(P1 = "C1"), cand, ec_ids = "E1")
  expect_setequal(sets$P1, "E1")
  expect_length(mmc_pathway_ec_sets(list(P1 = "C1"), cand, character(0)), 0)
  expect_length(mmc_pathway_ec_sets(list(), cand, "E1"), 0)
})


# ---------------------------------------------------------------------------
# NES interpretation (must NOT be a biological direction)
# ---------------------------------------------------------------------------

test_that("the NES note is conservative: parity only, no direction claim", {
  note <- mmc_gsea_nes_note()

  # what it MUST say
  expect_match(note, "retained to reproduce the pinned MetaboAnalystR result",
               fixed = TRUE)
  expect_match(note, "transforms and reorders the score vector after pathway positions are constructed",
               fixed = TRUE)
  expect_match(note, "should not be interpreted as biological up/down direction",
               fixed = TRUE)
  expect_match(note, "magnitude of the pathway members' original scores",
               fixed = TRUE)

  # what it must NOT say: pathway indices come from the SIGNED ordering and are
  # then scored against a transformed, re-sorted vector, so even a claim about
  # which |score| END a pathway sits at is unsupported.
  expect_false(grepl("high-|score| end", note, fixed = TRUE))
  expect_false(grepl("low-|score| end", note, fixed = TRUE))
  expect_false(grepl("enrichment toward", note, fixed = TRUE))
})

test_that("no biological-direction wording survives anywhere in the GSEA layer", {
  # Regression guard: MetaboAnalyst ranks on |score|, so NES sign must never be
  # presented as treatment direction. This locks out the wording that used to be
  # here (a `Direction` column and a numerator/denominator note) and anything
  # equivalent creeping back in.
  expect_false(exists("mmc_gsea_direction_note"))

  root <- normalizePath(if (dir.exists("R")) "." else "../..", mustWork = FALSE)
  targets <- c(
    file.path(root, "R", "domain", "metabolomics", "06g_mummichog_gsea.R"),
    file.path(root, "R", "domain", "metabolomics", "06e_mummichog_plots.R"),
    file.path(root, "R", "pipeline", "metabolomics", "templates",
              "report_metabolomics.Rmd")
  )
  banned <- c("toward numerator", "toward denominator",
              "up-regulated pathway", "down-regulated pathway",
              "upregulated", "downregulated",
              "direction_note", "mmc_gsea_direction_note",
              # over-strong magnitude claims: pathway positions come from the
              # signed ordering, so even an |score|-end reading is unsupported
              "high-|score| end", "low-|score| end",
              "enrichment toward the high", "enrichment toward the low")
  for (f in targets) {
    skip_if_not(file.exists(f), paste("missing", basename(f)))
    txt <- tolower(paste(readLines(f, warn = FALSE), collapse = "\n"))
    for (b in banned) {
      expect_false(grepl(tolower(b), txt, fixed = TRUE),
                   info = paste0("banned direction wording '", b,
                                 "' found in ", basename(f)))
    }
    # the only sanctioned NES sentence, or none at all
    expect_false(grepl("toward hl|toward ll", txt, fixed = FALSE),
                 info = paste("group-direction wording in", basename(f)))
  }
})

test_that("the GSEA result table carries no Direction column", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)
  res  <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "HL_vs_LL",
                             de_method = "limma")

  expect_false("Direction" %in% names(res$table))
  expect_false(any(grepl("direction", tolower(names(res$table)))))
  expect_null(res$direction_note)
  expect_identical(res$nes_note, mmc_gsea_nes_note())
  # the contrast is carried as provenance only
  expect_identical(res$contrast, "HL_vs_LL")
})


# ---------------------------------------------------------------------------
# end-to-end GSEA run
# ---------------------------------------------------------------------------

test_that("GSEA runs over the full ranked EC list and preserves NES sign", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  res <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "HL_vs_LL",
                             de_method = "limma")
  expect_false(is.null(res))

  # the FULL EC universe is ranked, not only the ORA-significant overlap
  expect_equal(res$n_ec, 24L)
  expect_identical(res$metric, "t_statistic")
  expect_setequal(res$table$Pathway, c("PW_pos", "PW_mid", "PW_neg"))

  expect_true(all(c("Pathway", "Pathway size (model compounds)", "Detected ECs",
                    "Tested size", "ES", "NES", "P.Value", "padj",
                    "Leading-edge EmpiricalCompounds") %in% names(res$table)))

  # signed ES/NES are preserved (not flattened to absolute values)
  expect_true(any(res$table$NES < 0) || any(res$table$ES < 0))
  expect_true(all(is.finite(res$table$NES)))

  # engine settings are the MetaboAnalyst-equivalent defaults
  expect_identical(res$params$n_perm, 100L)
  expect_equal(res$params$gsea_param, 1)
  expect_equal(res$params$min_size, 1)
  expect_equal(res$params$max_size, Inf)
  expect_identical(res$params$seed, 123L)

  # sorted by ascending raw p, and the size columns are distinct concepts
  expect_equal(res$table[["P.Value"]], sort(res$table[["P.Value"]]))
  expect_true(all(res$table[["Pathway size (model compounds)"]] == 8))
  expect_true(all(res$table[["Detected ECs"]] <= 8))

  # BH is taken over the tested pathways only
  expect_equal(res$table$padj,
               stats::p.adjust(res$table[["P.Value"]], method = "fdr"),
               tolerance = 1e-12)
})

test_that("GSEA is reproducible for a fixed seed", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  a <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "HL_vs_LL",
                             de_method = "limma")
  b <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "HL_vs_LL",
                             de_method = "limma")
  expect_equal(a$table[["P.Value"]], b$table[["P.Value"]])
  expect_equal(a$table$NES, b$table$NES)
  expect_equal(a$table$ES, b$table$ES)
  # and the default seed is upstream's, not an invented one
  expect_identical(a$params$seed, 123L)
})

test_that("GSEA skips gracefully when a prerequisite is missing", {
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  # no model
  expect_message(expect_null(
    run_mummichog_gsea(f$files, f$de_table, NULL, contrast = "c")),
    "no metabolic model")
  # no usable ranking statistic
  expect_message(expect_null(
    run_mummichog_gsea(f$files, data.frame(feature_id = "x", AveExpr = 1),
                       f$model, contrast = "c")),
    "no usable ranking statistic")
  # no EC tables among the files
  empty <- withr::local_tempdir()
  expect_message(expect_null(
    run_mummichog_gsea(list_mummichog_files(empty), f$de_table, f$model,
                       contrast = "c")),
    "no EmpiricalCompound tables")
})


# ---------------------------------------------------------------------------
# the summary scatter
# ---------------------------------------------------------------------------

test_that("the GSEA scatter maps NES on x, -log10(p) on y and NES to colour", {
  g <- gsea_table_fixture()
  p <- plot_mummichog_gsea_scatter(g, title = "T", subtitle = "S",
                                   p_cutoff = 0.05)

  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "T")
  expect_identical(p$labels$subtitle, "S")
  expect_match(p$labels$x, "NES")
  expect_match(p$labels$y, "-log10\\(GSEA p-value\\)")
  expect_identical(
    p$scales$get_scales("fill")$name,
    "NES"
  )

  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomPoint" %in% geoms)
  expect_true("GeomHline" %in% geoms)   # the p = 0.05 reference line

  b  <- ggplot2::ggplot_build(p)
  pd <- b$data[[which(geoms == "GeomPoint")[1]]]
  expect_equal(sort(pd$x), sort(g$NES))
  expect_equal(sort(pd$y), sort(-log10(g[["P.Value"]])), tolerance = 1e-6)
})

test_that("the GSEA scatter draws the p = 0.05 reference line", {
  p <- plot_mummichog_gsea_scatter(gsea_table_fixture(), title = "T",
                                   p_cutoff = 0.05)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  b  <- ggplot2::ggplot_build(p)
  hl <- b$data[[which(geoms == "GeomHline")[1]]]

  expect_equal(unique(hl$yintercept), -log10(0.05))
  expect_identical(unique(hl$linetype), "dashed")

  # a different cutoff moves the line
  p2 <- plot_mummichog_gsea_scatter(gsea_table_fixture(), title = "T",
                                    p_cutoff = 0.01)
  geoms2 <- vapply(p2$layers, function(l) class(l$geom)[1], character(1))
  hl2 <- ggplot2::ggplot_build(p2)$data[[which(geoms2 == "GeomHline")[1]]]
  expect_equal(unique(hl2$yintercept), -log10(0.01))
})

test_that("the GSEA scatter's colour scale is diverging and centred at NES = 0", {
  g <- gsea_table_fixture()
  p <- plot_mummichog_gsea_scatter(g, title = "T")
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  pd <- ggplot2::ggplot_build(p)$data[[which(geoms == "GeomPoint")[1]]]

  # NES == 0 lands exactly on the scale's midpoint colour
  mid <- pd$fill[which(abs(pd$x) < 1e-12)]
  expect_equal(toupper(substr(mid, 1, 7)), "#FFFEE0")
  # positive and negative NES sit on opposite arms of the scale
  expect_false(identical(pd$fill[which.max(pd$x)], pd$fill[which.min(pd$x)]))
})

test_that("the GSEA scatter's point size is MetaboAnalyst's sqrt(-log10 p)", {
  g <- gsea_table_fixture()
  p <- plot_mummichog_gsea_scatter(g, title = "T")

  expect_equal(p$data$significance_size, sqrt(-log10(g[["P.Value"]])),
               tolerance = 1e-9)
  # the legend names the quantity honestly and never leaks the source's
  # internal variable name
  size_name <- p$scales$get_scales("size")$name
  expect_match(size_name, "Significance")
  expect_match(size_name, "log10 p")
  expect_false(grepl("radi", size_name, ignore.case = TRUE))
})

test_that("the GSEA scatter returns NULL on empty / unusable input", {
  expect_null(plot_mummichog_gsea_scatter(NULL, title = "T"))
  expect_null(plot_mummichog_gsea_scatter(data.frame(), title = "T"))
  expect_null(plot_mummichog_gsea_scatter(
    data.frame(Pathway = "x", NES = 1, check.names = FALSE), title = "T"))
})
