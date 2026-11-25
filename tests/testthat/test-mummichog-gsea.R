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

test_that("a non-numeric statistic column is not ranked on", {
  # as.numeric() on a factor returns LEVEL CODES, not the labels. Sorted
  # alphabetically, "-6" and "4" become 1 and 2 -- a ranking that inverts the
  # real one here, and would have been reported as the signed t statistic.
  de <- data.frame(stringsAsFactors = FALSE,
                   feature_id = c("f1", "f2"),
                   statistic  = factor(c("-6", "4")),
                   logFC      = c(-2, 1))
  r <- mmc_gsea_rank_metric(de, de_method = "limma")

  expect_identical(r$metric, "logFC")
  expect_match(r$label, "no usable t statistic")
  expect_equal(unname(r$values[["f1"]]), -2)   # the real values, not 1 and 2
  # The same applies to the fallback column, which then has nothing to offer.
  de$logFC <- factor(c("-2", "1"))
  expect_null(mmc_gsea_rank_metric(de, de_method = "limma"))

  # A character column is refused on the same rule -- type is checked, never
  # coerced -- and an integer one is a perfectly good statistic.
  de_chr <- data.frame(stringsAsFactors = FALSE,
                       feature_id = c("f1", "f2"),
                       statistic  = c("-6", "4"),
                       logFC      = c(-2, 1))
  expect_identical(mmc_gsea_rank_metric(de_chr, de_method = "limma")$metric,
                   "logFC")

  de_int <- data.frame(stringsAsFactors = FALSE,
                       feature_id = c("f1", "f2"),
                       statistic  = c(-6L, 4L),
                       logFC      = c(-2, 1))
  r_int <- mmc_gsea_rank_metric(de_int, de_method = "limma")
  expect_identical(r_int$metric, "t_statistic")
  expect_equal(unname(r_int$values[["f1"]]), -6)
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
               "gsea_param must be the number 1")
  expect_error(mmc_gsea_metaboanalyst(pw, st, rk, gsea_param = 0.5),
               "gsea_param must be the number 1")
  # The message has to say why, not just refuse.
  expect_error(mmc_gsea_metaboanalyst(pw, st, rk, gsea_param = 2),
               "permutation")

  # Type before value. "1" would pass a numeric comparison and then fail inside
  # abs(stats)^"1", where run_mummichog_gsea()'s tryCatch turns a caller error
  # into "GSEA not applicable"; TRUE would run as 1 but be stored as TRUE.
  for (bad in list("1", TRUE, NA, c(1, 1), NULL)) {
    expect_error(mmc_gsea_metaboanalyst(pw, st, rk, gsea_param = bad),
                 "gsea_param must be the number 1")
  }
})

test_that("run_mummichog_gsea rejects gsea_param early, before any work", {
  # Early and loud: inside the run it would be caught by the engine's tryCatch
  # and degrade into a warning plus a NULL result, which reads like "GSEA was
  # not applicable" rather than "you asked for something invalid".
  expect_error(
    run_mummichog_gsea(character(0), data.frame(feature_id = "f", logFC = 1),
                       model = list(pathways = list(P = "C1")),
                       de_method = "limma", gsea_param = 2),
    "gsea_param must be the number 1")
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

  # Legacy: 0 is not a departure from the reference either, because 1 is the
  # reference. (Asserted on the legacy method: only it reports deviations, so
  # on the default method this could not fail.)
  old <- with_legacy_gsea(run_mummichog_gsea(
    f$files, f$de_table, f$model, contrast = "c", de_method = "limma",
    min_size = 0, method = "metaboanalyst_legacy"))
  expect_equal(old$params$min_size, 1)
  expect_false(any(grepl("min_size", old$deviations)))
  # ... while a real override IS reported there, so the check above can fail.
  three <- with_legacy_gsea(run_mummichog_gsea(
    f$files, f$de_table, f$model, contrast = "c", de_method = "limma",
    min_size = 3, method = "metaboanalyst_legacy"))
  expect_true(any(grepl("min_size = 3", three$deviations)))

  same <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                             de_method = "limma", min_size = 1)
  expect_equal(res$table, same$table)
})

test_that("legacy: a non-default rank_metric is recorded as a deviation", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  # Merging features that share an m/z by something other than the reference's
  # "mean" changes every EC score, and so every statistic downstream of it.
  res <- with_legacy_gsea(run_mummichog_gsea(
    f$files, f$de_table, f$model, contrast = "c", de_method = "limma",
    rank_metric = "median", method = "metaboanalyst_legacy"))
  expect_true(any(grepl("rank_metric = 'median'", res$deviations, fixed = TRUE)))

  ref <- with_legacy_gsea(run_mummichog_gsea(
    f$files, f$de_table, f$model, contrast = "c", de_method = "limma",
    method = "metaboanalyst_legacy"))
  expect_false(any(grepl("rank_metric", ref$deviations)))
  expect_identical(ref$params$rank_metric, "mean")

  # Under the default method it is a setting, reported as the value that ran.
  std <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                            de_method = "limma", rank_metric = "median")
  expect_length(std$deviations, 0L)
  expect_identical(std$params$rank_metric, "median")
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
  # Run on the legacy method, the only one that reports deviations: on the
  # default method the "is not a deviation" check below could not fail.
  legacy <- function(metric) with_legacy_gsea(run_mummichog_gsea(
    f$files, f$de_table, f$model, contrast = "c", de_method = "limma",
    rank_metric = metric, method = "metaboanalyst_legacy"))
  abbrev <- legacy("mea")

  expect_identical(abbrev$params$rank_metric, "mean")
  expect_false(any(grepl("rank_metric", abbrev$deviations)))

  full <- legacy("mean")
  expect_equal(abbrev$table, full$table)

  # The default method resolves the abbreviation the same way.
  std <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                            de_method = "limma", rank_metric = "mea")
  expect_identical(std$params$rank_metric, "mean")

  # A non-default abbreviation still resolves, and is still a deviation -- the
  # normalisation must not swallow the departure along with the shorthand.
  med <- with_legacy_gsea(run_mummichog_gsea(
    f$files, f$de_table, f$model, contrast = "c", de_method = "limma",
    rank_metric = "medi", method = "metaboanalyst_legacy"))
  expect_identical(med$params$rank_metric, "median")
  expect_true(any(grepl("rank_metric = 'median'", med$deviations, fixed = TRUE)))

  # Both entry points answer to the same list of allowed values.
  expect_identical(.MMC_GSEA_RANK_METRICS[[1]], .MMC_GSEA_DEFAULTS$rank_metric)
  expect_error(run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                                  de_method = "limma", rank_metric = "nope"))
})

test_that("legacy: the single-pathway branch is reproducible and declared unvalidated", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)
  # One pathway -> the engine takes its seeded sample.int() loop instead of the
  # batch call. Upstream leaves that loop to a non-deterministic worker RNG, so
  # there is no reference draw to match: we offer reproducibility instead, and
  # the result has to say so rather than imply parity.
  one <- list(pathways = f$model$pathways["PW_pos"])

  legacy <- function(model) with_legacy_gsea(run_mummichog_gsea(
    f$files, f$de_table, model, contrast = "c", de_method = "limma",
    method = "metaboanalyst_legacy"))
  a <- legacy(one)
  b <- legacy(one)

  expect_equal(nrow(a$table), 1L)
  expect_equal(a$table, b$table)          # seeded: identical across runs

  expect_true(any(grepl("one permutable pathway", a$deviations)))
  expect_true(any(grepl("not covered by the parity fixture", a$deviations)))

  # More than one permutable pathway goes through the batch call and claims
  # parity.
  many <- legacy(f$model)
  expect_gt(nrow(many$table), 1L)
  expect_false(any(grepl("one permutable pathway", many$deviations)))
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

test_that("the m/z merge weights a feature by its EC multiplicity (as built)", {
  # DOCUMENTS CURRENT BEHAVIOUR. It is NOT an upstream-parity assertion: the
  # m/z merge has no counterpart in the single-dataset path this engine's stage
  # runs (upstream drops duplicate m/z rows before scores exist, and uses
  # rank.metric only in its meta-analysis function). Pinned here so the choice
  # is visible and cannot change unnoticed while the question is open.
  #
  # fA is mapped to TWO EmpiricalCompounds, so read_mummichog_ec_features()
  # gives it two rows; fB shares its m/z and has one. With "mean" over rows,
  # the shared m/z resolves to (2*3 + 0)/3 = 2, not the per-feature (3 + 0)/2.
  ec_features <- data.frame(
    stringsAsFactors = FALSE,
    EID        = c("E1", "E2", "E3", "E4"),
    feature_id = c("fA", "fA", "fB", "fC"),
    mz         = c(100, 100, 100, 200))
  scores <- c(fA = 3, fB = 0, fC = -5)

  got <- mmc_ec_scores(ec_features, scores)

  expect_equal(unname(got[["E1"]]), 2)   # row-weighted, not 1.5
  expect_equal(unname(got[["E2"]]), 2)
  expect_equal(unname(got[["E3"]]), 2)
  expect_equal(unname(got[["E4"]]), -5)  # its own m/z, untouched by the merge

  # The per-feature mean a de-duplicating implementation would give.
  expect_false(isTRUE(all.equal(unname(got[["E1"]]), 1.5)))

  # Step 2 is upstream's rule and is unaffected: an EC takes the SIGNED max of
  # its members, so the negative score is kept rather than replaced by a
  # larger magnitude.
  expect_equal(unname(got[["E4"]]), -5)
  expect_true(all(diff(unname(got)) <= 0))   # returned in decreasing order
})

test_that("the engine's whole-number settings are validated, not coerced", {
  pw <- list(PW = c("E1", "E2"))
  st <- stats::setNames(c(2, 1, -1), c("E1", "E2", "E3"))
  rk <- stats::setNames(c(1, 2, 3), c("E1", "E2", "E3"))
  de    <- data.frame(feature_id = "f", logFC = 1, stringsAsFactors = FALSE)
  model <- list(pathways = list(P = "C1"))

  # All four go through .mmc_as_seed(), the pinned stage's existing validator:
  # type and range before coercion, because as.integer() turns 1.5 into 1 and
  # "7" into 7 without a word, and a value past .Machine$integer.max into NA.
  for (bad in list(100.5, "100", TRUE, NA, c(10, 20), Inf,
                   .Machine$integer.max + 1)) {
    expect_error(mmc_gsea_metaboanalyst(pw, st, rk, n_perm = bad), "n_perm")
  }
  expect_error(mmc_gsea_metaboanalyst(pw, st, rk, seed = 123.9), "seed")
  expect_error(mmc_gsea_metaboanalyst(pw, st, rk, seed = .Machine$integer.max + 1),
               "seed")
  expect_error(mmc_gsea_metaboanalyst(pw, st, rk, n_perm = 0),
               "n_perm must be at least 1")

  # Both entry points refuse before any work, so the caller gets the error
  # rather than "GSEA skipped".
  expect_error(run_mummichog_gsea(character(0), de, model, de_method = "limma",
                                  n_perm = 100.5), "n_perm")
  expect_error(run_mummichog_gsea(character(0), de, model, de_method = "limma",
                                  seed = 123.9), "seed")
})

test_that("an invalid size bound errors instead of emptying the result", {
  pw <- list(PW = c("E1", "E2"))
  st <- stats::setNames(c(2, 1, -1), c("E1", "E2", "E3"))
  rk <- stats::setNames(c(1, 2, 3), c("E1", "E2", "E3"))
  de    <- data.frame(feature_id = "f", logFC = 1, stringsAsFactors = FALSE)
  model <- list(pathways = list(P = "C1"))

  # The worst shape this class takes. `min_size <= sizes & sizes <= max_size`
  # is all-NA for an NA bound, which() returns nothing, and the run reports
  # "no pathway passed the size filter" -- a valid experiment disappearing as
  # a plausible empty result that reads as a finding about the data.
  for (bad in list(NA, NA_real_, "3", TRUE, c(1, 2), -Inf)) {
    expect_error(mmc_gsea_metaboanalyst(pw, st, rk, min_size = bad), "min_size")
    expect_error(mmc_gsea_metaboanalyst(pw, st, rk, max_size = bad), "max_size")
  }
  expect_error(run_mummichog_gsea(character(0), de, model, de_method = "limma",
                                  max_size = NA), "max_size")

  # Inf is the reference default for the upper bound and must stay legal; the
  # lower bound has no such case.
  expect_identical(.mmc_gsea_check_max_size(Inf), Inf)
  expect_identical(.mmc_gsea_check_max_size(5), 5L)
  expect_error(.mmc_gsea_check_max_size(-Inf), "max_size")
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
  # The null was substituted, not drawn, so the row is one the reference never
  # produces -- counted as such rather than passed off as a permuted result.
  expect_identical(attr(res, "mmc_unpermutable"), 1L)
  # And no seeded-RNG claim: with nothing to draw, that loop never ran.
  expect_false(attr(res, "mmc_single_pathway_rng"))
})

test_that("a sole pathway whose ties outnumber the ranking keeps its row", {
  for (pkg in c("fgsea", "fastmatch", "data.table", "BiocParallel")) {
    skip_if_not_installed(pkg)
  }
  # The second route to K >= length(stats), and the one the all-position fix
  # did not cover. Two ECs sharing a score collapse into ONE ranked position,
  # so length(stats) == 1, while the pathway's tie-group index vector is
  # c(1, 1) and K == 2. sample.int(1, 2) then raises *cannot take a sample
  # larger than the population* -- above the calcGseaStat catch, so the engine
  # aborted and the documented ES-defaulted row was lost.
  ec  <- c(E1 = 5, E2 = 5)
  inp <- mmc_gsea_ranked_inputs(ec)
  expect_length(inp$stats, 1L)          # the tie really does collapse

  res <- mmc_gsea_metaboanalyst(list(PW_tied = c("E1", "E2")),
                                inp$stats, inp$ranks)

  expect_equal(nrow(res), 1L)
  expect_true(res$es_defaulted)
  expect_match(res$es_reason, "^tied EC scores")
  expect_identical(res$ES, 0)
  expect_identical(res$pval, 1)
})

test_that("a tie-expanded pathway keeps its row beside a permutable one", {
  for (pkg in c("fgsea", "fastmatch", "data.table", "BiocParallel")) {
    skip_if_not_installed(pkg)
  }
  # The reported case: the guard added for a SOLE oversized pathway sat inside
  # the single-pathway branch, so a second tested pathway sent the run down the
  # other path and the engine aborted again on the same shape.
  ec  <- c(E1 = 5, E2 = 5)
  inp <- mmc_gsea_ranked_inputs(ec)
  expect_length(inp$stats, 1L)          # the tie really does collapse

  res <- mmc_gsea_metaboanalyst(
    list(PW_tied = c("E1", "E2"), PW_one = c("E1")), inp$stats, inp$ranks)

  expect_equal(nrow(res), 2L)
  tied <- res[res$pathway == "PW_tied", ]
  expect_true(tied$es_defaulted)
  expect_match(tied$es_reason, "^tied EC scores")
  expect_identical(tied$ES, 0)
  expect_identical(tied$pval, 1)
  # Neither row was permuted: PW_tied has no draw at all, and PW_one is the
  # sole remaining pathway while covering the whole (one-position) ranking, so
  # its null is substituted too. Both are counted, and no seeded-RNG deviation
  # is claimed for a loop that never drew anything.
  expect_identical(attr(res, "mmc_unpermutable"), 2L)
  expect_false(attr(res, "mmc_single_pathway_rng"))
})

test_that("an unpermutable pathway does not abort the permutation batch", {
  for (pkg in c("fgsea", "fastmatch", "data.table", "BiocParallel")) {
    skip_if_not_installed(pkg)
  }
  # Two permutable pathways put the run on the cumulative-batch path, where an
  # oversized set is handed straight to calcGseaStatCumulativeBatch(): it
  # cannot draw 4 distinct positions from a 3-position ranking, so the whole
  # engine threw and every row was lost, not just the oversized one.
  ec  <- c(E1 = 5, E2 = 5, E3 = 3, E4 = 1)
  inp <- mmc_gsea_ranked_inputs(ec)
  expect_length(inp$stats, 3L)          # E1/E2 share one position

  res <- mmc_gsea_metaboanalyst(
    list(PW_big = c("E1", "E2", "E3", "E4"),   # index vector of length 4
         PW_a   = c("E3", "E4"),
         PW_b   = c("E1", "E3")),
    inp$stats, inp$ranks)

  expect_equal(nrow(res), 3L)
  big <- res[res$pathway == "PW_big", ]
  expect_identical(big$size, 4)         # tie-expanded beyond the ranking
  expect_true(big$es_defaulted)
  expect_match(big$es_reason, "^tied EC scores")
  expect_identical(big$ES, 0)
  expect_identical(big$pval, 1)
  # The pathways that CAN be permuted are scored normally: the oversized one is
  # excluded from the draw, not allowed to default its neighbours.
  rest <- res[res$pathway != "PW_big", ]
  expect_false(any(rest$es_defaulted))
  # Two permutable pathways means the batch ran, not the seeded loop.
  expect_false(attr(res, "mmc_single_pathway_rng"))
  # And the engine says how many rows it produced where the reference aborts,
  # so the provenance does not have to be reconstructed from the sizes.
  expect_identical(attr(res, "mmc_unpermutable"), 1L)
})

test_that("an engine failure is not reported as a size-filter miss", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  # Stub by assignment into the function's own environment: these are sourced
  # functions, not a package, so with_mocked_bindings() has nothing to rebind.
  # The default method's engine is the one stubbed: it is what a run uses.
  target <- environment(run_mummichog_gsea)
  orig   <- get("mmc_gsea_fgsea", envir = target)
  assign("mmc_gsea_fgsea",
         function(...) stop("synthetic engine failure"), envir = target)
  withr::defer(assign("mmc_gsea_fgsea", orig, envir = target))

  msgs <- capture_messages(
    expect_warning(
      res <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                                de_method = "limma"),
      "engine failed"))

  expect_null(res)
  # The warning already named the real cause; saying the size filter emptied
  # the result would dress an engine failure up as a finding about the data.
  expect_false(any(grepl("size filter", msgs, fixed = TRUE)))
})

test_that("a genuinely empty result still names the size filter", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  # The other side of the split: when filtering really did run and kept
  # nothing, that is exactly what the message should say.
  msgs <- capture_messages(
    res <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                              de_method = "limma", min_size = 1000))

  expect_null(res)
  expect_true(any(grepl("size filter", msgs, fixed = TRUE)))
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

test_that("the default NES note describes the ranking, not a per-metabolite direction", {
  note <- mmc_gsea_nes_note()
  expect_identical(note, mmc_gsea_nes_note("fgsea"))

  # what it MUST say
  expect_match(note, "concentrated among the highest-scoring ECs", fixed = TRUE)
  expect_match(note, "largest signed score among its features", fixed = TRUE)
  expect_match(note, "not a per-metabolite direction of change", fixed = TRUE)

  # what it must NOT say: it is not the legacy method's parity sentence, and it
  # makes no up/down claim about metabolites or pathways.
  expect_false(grepl("retained to reproduce", note, fixed = TRUE))
  expect_false(grepl("MetaboAnalyst", note, fixed = TRUE))
  for (banned in c("up-regulated", "down-regulated", "upregulated",
                   "downregulated", "numerator", "denominator")) {
    expect_false(grepl(banned, note, fixed = TRUE), info = banned)
  }
  # an unknown or abbreviated method is refused, never guessed
  expect_error(mmc_gsea_nes_note("metaboanalyst"), "gsea_method")
})

test_that("legacy: the NES note is conservative: parity only, no direction claim", {
  note <- mmc_gsea_nes_note("metaboanalyst_legacy")

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

  # engine settings: the default method records fgsea's nPermSimple and has
  # no fixed permutation count; the rest are the MetaboAnalyst-equivalent
  # defaults
  expect_identical(res$params$nperm_simple, 10000L)
  expect_null(res$params$n_perm)
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


# ---------------------------------------------------------------------------
# method selection (gsea_method)
# ---------------------------------------------------------------------------

test_that("gsea_method resolves to standard fgsea unless legacy is named", {
  expect_identical(mmc_gsea_method(list()), "fgsea")
  expect_identical(mmc_gsea_method(list(gsea_enabled = TRUE)), "fgsea")
  expect_identical(mmc_gsea_method(list(gsea_method = "fgsea")), "fgsea")
  expect_identical(mmc_gsea_method(list(gsea_method = "metaboanalyst_legacy")),
                   "metaboanalyst_legacy")
  # A longer key that merely starts with the name must not be read as it.
  expect_identical(mmc_gsea_method(list(gsea_method_note = "metaboanalyst_legacy")),
                   "fgsea")
})

test_that("an unknown gsea_method stops instead of falling back", {
  # Exact matching: an abbreviation, the ORA engine's name, another case, a
  # blank and a non-character value are all refused. The key decides which
  # statistic runs, so none of them may quietly become the default or legacy.
  for (bad in list("metaboanalyst", "legacy", "FGSEA", "", NA_character_, 1,
                   TRUE, c("fgsea", "metaboanalyst_legacy"))) {
    expect_error(mmc_gsea_method(list(gsea_method = bad)), "gsea_method")
  }
  # ... and the functions that take the method validate it the same way.
  expect_error(mmc_gsea_method_label("metaboanalyst"), "gsea_method")
  expect_error(run_mummichog_gsea(character(0), data.frame(feature_id = "f", logFC = 1),
                                  list(pathways = list()), method = "metaboanalyst"),
               "gsea_method")
})

test_that("each method is described once, and the descriptions do not overlap", {
  std <- mmc_gsea_method_label("fgsea", "1.36.2")
  old <- mmc_gsea_method_label("metaboanalyst_legacy", "1.36.2")
  expect_match(std, "fgsea::fgsea, adaptive multilevel, fgsea 1.36.2", fixed = TRUE)
  # "fgsea::fgsea" is a prefix of the earlier default's function name, so the
  # match above would not notice that name coming back.
  expect_false(grepl("fgseaSimple", std, fixed = TRUE))
  expect_match(std, "one ranking position per EC", fixed = TRUE)
  expect_false(grepl("MetaboAnalyst", std, fixed = TRUE))
  expect_match(old, "legacy MetaboAnalyst-compatible GSEA", fixed = TRUE)
  expect_match(old, "known defects", fixed = TRUE)
  expect_match(old, .MMC_GSEA_LEGACY_ISSUE, fixed = TRUE)
  # an unknown version is left out rather than printed as NA
  expect_false(grepl("NA", mmc_gsea_method_label("fgsea", NA_character_), fixed = TRUE))
  expect_false(grepl("fgsea ,", mmc_gsea_method_label("fgsea"), fixed = TRUE))
})


# ---------------------------------------------------------------------------
# the default method: standard fgsea::fgsea() (multilevel) on the signed EC
# scores, through the pipeline's shared run_fgsea_core()
# ---------------------------------------------------------------------------

test_that("standard fgsea weights each EC by its own signed score (hand-calculated)", {
  skip_if_not_installed("fgsea")
  # signed order: A=+3, B=+1, D=+0.5, E=-2, C=-5 (N = 5), two members each
  #   P1 = {A, B}: positions 1, 2; weights 3, 1 -> ES = 1, leading edge A, B
  #   P2 = {A, C}: positions 1, 5; weights 3, 5; NR = 8
  #                tops 3/8, 0; bottoms 0, -5/8 -> ES = -5/8, leading edge C
  ec <- c(A = 3, B = 1, D = 0.5, E = -2, C = -5)
  pw <- list(P1 = c("A", "B"), P2 = c("A", "C"))

  res <- mmc_gsea_fgsea(pw, ec)
  expect_identical(res$pathway, c("P1", "P2"))
  expect_equal(res$ES, c(1, -5 / 8), tolerance = 1e-9)
  expect_identical(res$leadingEdge, c("A; B", "C"))
  expect_identical(res$size, c(2, 2))
  expect_false(any(res$es_defaulted))
  expect_true(all(is.na(res$es_reason)))
  # The legacy engine's columns, in its order, then fgsea's log2err.
  expect_identical(names(res), c(
    names(mmc_gsea_metaboanalyst(list(), numeric(0), numeric(0))), "log2err"))
  expect_true(all(is.na(res$nMoreExtreme)))            # multilevel reports none

  # Historical compatibility, on the same input: the legacy engine weights P2's
  # hits by scores that are not theirs, so its ES is positive and its leading
  # edge names C, the most negative EC.
  ri  <- mmc_gsea_ranked_inputs(ec)
  old <- mmc_gsea_metaboanalyst(pw, ri$stats, ri$ranks, n_perm = 100)
  expect_equal(old$ES, c(1, 10 / 11), tolerance = 1e-9)
  expect_identical(old$leadingEdge, c("A", "C"))
})

# 60 synthetic ECs with both signs and tied scores (rounded to one decimal),
# and 12 random pathway sets. Shared by the tests below.
gsea_tied_fixture <- function() {
  ec <- withr::with_seed(11, stats::setNames(round(stats::rnorm(60), 1),
                                             sprintf("E%02d", sample(60))))
  sets <- withr::with_seed(12, lapply(1:12, function(i)
    sample(names(ec), sample(4:15, 1))))
  names(sets) <- sprintf("S%02d", 1:12)
  list(ec = ec, sets = sets)
}

test_that("standard fgsea equals fgsea::calcGseaStat with tied ECs ordered by id", {
  skip_if_not_installed("fgsea")
  fx <- gsea_tied_fixture()
  expect_gt(anyDuplicated(fx$ec), 0L)                  # tied scores are present
  expect_true(any(fx$ec < 0) && any(fx$ec > 0))        # and both signs

  res <- suppressWarnings(mmc_gsea_fgsea(fx$sets, fx$ec))

  # The documented order: score decreasing, equal scores by EC id (C locale).
  # fgsea's own primitive on that order is the expected value, so this fails
  # if tied ECs were ever scored in any other order.
  ordered <- fx$ec[order(-fx$ec, names(fx$ec), method = "radix")]
  want <- lapply(fx$sets, function(m) fgsea::calcGseaStat(
    ordered, sort(match(m, names(ordered))), returnLeadingEdge = TRUE))

  expect_identical(res$pathway, names(fx$sets))
  expect_equal(res$ES, vapply(want, function(x) x$res, numeric(1),
                              USE.NAMES = FALSE), tolerance = 1e-6)
  expect_identical(res$size, as.numeric(lengths(fx$sets)))
  for (i in seq_along(fx$sets)) {
    got <- strsplit(res$leadingEdge[i], "; ", fixed = TRUE)[[1]]
    # leading-edge ids are the pathway's own ECs, and the ones fgsea reports
    expect_true(all(got %in% fx$sets[[i]]), info = names(fx$sets)[i])
    expect_setequal(got, names(ordered)[want[[i]]$leadingEdge])
  }
  # BH over the tested pathways, as fgsea computes it
  expect_equal(res$padj, stats::p.adjust(res$pval, method = "BH"))
})

test_that("standard fgsea is reproducible and does not depend on input order", {
  skip_if_not_installed("fgsea")
  fx <- gsea_tied_fixture()
  run <- function(sets, ec, seed = 5) {
    suppressWarnings(mmc_gsea_fgsea(sets, ec, seed = seed))
  }
  a <- run(fx$sets, fx$ec)
  expect_identical(a, run(fx$sets, fx$ec))             # same seed, same result

  # Shuffling the score vector and each pathway's members changes nothing:
  # the order of tied ECs comes from their ids, not from how they arrived.
  shuffled <- withr::with_seed(99, fx$ec[sample(length(fx$ec))])
  expect_equal(a, run(lapply(fx$sets, rev), shuffled))

  # A single pathway is reproducible too.
  one <- run(fx$sets[1], fx$ec)
  expect_identical(nrow(one), 1L)
  expect_identical(one, run(fx$sets[1], fx$ec))

  # The scores are passed unchanged: no jitter, so tied ECs stay tied.
  ties <- attr(a, "mmc_ties")
  expect_identical(ties$n_tied_ec,
                   sum(duplicated(fx$ec) | duplicated(fx$ec, fromLast = TRUE)))
})

test_that("tied ECs keep their own positions: scored, counted and flagged exactly", {
  skip_if_not_installed("fgsea")
  ec <- c(E1 = 3, E2 = 2, E3 = 2, E4 = 1, E5 = -1, E6 = -1, E7 = -4)
  pw <- list(within = c("E2", "E3"),   # both ECs of one tie group
             mixed  = c("E1", "E2"),   # E2 is tied with the non-member E3
             clean  = c("E1", "E7"),   # no tied member
             neg    = c("E5", "E7"))   # E5 is tied with the non-member E6

  expect_warning(res <- mmc_gsea_fgsea(pw, ec),
                 "can depend on the order")
  expect_identical(res$pathway, names(pw))

  # "within": positions 2 and 3 of 7, weights 2 and 2 (NR = 4)
  #   tops 0.5 - 1/5, 1 - 1/5 = 0.3, 0.8 -> ES = 0.8, a real score, not a fallback
  expect_equal(res$ES[res$pathway == "within"], 0.8, tolerance = 1e-9)
  expect_identical(res$leadingEdge[res$pathway == "within"], "E2; E3")
  expect_false(any(res$es_defaulted))

  ties <- attr(res, "mmc_ties")
  expect_identical(ties$n_ec, 7L)
  expect_identical(ties$n_tied_ec, 4L)
  expect_equal(ties$prop_tied_ec, 4 / 7)
  expect_identical(ties$n_tie_groups, 2L)
  # Exactly the pathways with a tie group holding members AND non-members.
  # "within" holds its whole tie group, so the order inside it cannot matter.
  expect_identical(ties$tie_sensitive_pathways, c("mixed", "neg"))

  # Historical compatibility: the legacy engine puts E2 and E3 at ONE ranked
  # position, cannot score "within" and falls back to ES = 0.
  ri  <- mmc_gsea_ranked_inputs(ec)
  old <- mmc_gsea_metaboanalyst(pw["within"], ri$stats, ri$ranks, n_perm = 100)
  expect_true(old$es_defaulted)
  expect_equal(old$ES, 0)
  expect_lt(length(ri$stats), length(ec))   # 5 ranked positions for 7 ECs
})

test_that("the tie summary counts ties exactly and needs no threshold", {
  # No fgsea needed: this is pure bookkeeping on the scores and the sets.
  none <- mmc_gsea_tie_summary(c(E1 = 3, E2 = 1, E3 = -2), list(P = c("E1", "E2")))
  expect_identical(none$n_tied_ec, 0L)
  expect_identical(none$n_tie_groups, 0L)
  expect_identical(none$tie_sensitive_pathways, character(0))
  expect_equal(none$prop_tied_ec, 0)

  # Doubles that differ only beyond 15 significant digits are NOT a tie.
  near <- c(E1 = 1, E2 = 1 + 2e-16 * 4, E3 = 0)
  expect_identical(mmc_gsea_tie_summary(near, list())$n_tied_ec, 0L)

  # A pathway with no ranked member, and a zero-score tie, are both handled.
  z <- mmc_gsea_tie_summary(c(E1 = 0, E2 = 0, E3 = 1),
                            list(absent = "EX", zero = "E1", both = c("E1", "E2")))
  expect_identical(z$n_tied_ec, 2L)
  expect_identical(z$tie_sensitive_pathways, "zero")

  # Empty inputs give an empty summary, not an error.
  e <- mmc_gsea_tie_summary(stats::setNames(numeric(0), character(0)), list())
  expect_identical(e$n_ec, 0L)
  expect_identical(e$tie_sensitive_pathways, character(0))

  # The sentence states the same counts, and nothing for the legacy method.
  expect_null(mmc_gsea_tie_sentence(NULL))
  expect_match(mmc_gsea_tie_sentence(none),
               "None of the 3 ranked EmpiricalCompounds share a score", fixed = TRUE)
  expect_match(mmc_gsea_tie_sentence(z),
               "2 of the 3 ranked EmpiricalCompounds (66.7%)", fixed = TRUE)
  expect_match(mmc_gsea_tie_sentence(z), "1 tested pathway(s) can depend",
               fixed = TRUE)
})

test_that("fgsea's own warnings reach the caller, except the tie warning it replaces", {
  skip_if_not_installed("fgsea")
  # All scores positive: fgsea's scoreType advice has nothing to do with ties
  # and must not be swallowed.
  pos <- c(E1 = 4, E2 = 3, E3 = 2, E4 = 1, E5 = 0.5)
  expect_warning(mmc_gsea_fgsea(list(P = c("E1", "E2")), pos),
                 "greater than zero")

  # Tied scores: fgsea would say their order is arbitrary. It is not here, so
  # that one warning is replaced by ours, which says what the order is.
  tied <- c(E1 = 3, E2 = 2, E3 = 2, E4 = -1, E5 = -2)
  w <- testthat::capture_warnings(
    mmc_gsea_fgsea(list(P = c("E1", "E2")), tied))
  expect_false(any(grepl("arbitrary", w, fixed = TRUE)))
  expect_true(any(grepl("can depend on the order", w, fixed = TRUE)))

  # No ties, both signs: nothing to warn about.
  quiet <- c(E1 = 3, E2 = 2, E3 = 1, E4 = -1, E5 = -2)
  expect_warning(mmc_gsea_fgsea(list(P = c("E1", "E2")), quiet), NA)
})

test_that("standard fgsea validates its settings and returns one shape when empty", {
  skip_if_not_installed("fgsea")
  ec <- c(E1 = 3, E2 = 1, E3 = -1, E4 = -2)
  pw <- list(P = c("E1", "E2"))
  # No permutation count to set: nPermSimple is fixed (see the test below).
  expect_false("n_perm" %in% names(formals(mmc_gsea_fgsea)))
  expect_error(mmc_gsea_fgsea(pw, ec, seed = "7"), "seed")
  expect_error(mmc_gsea_fgsea(pw, ec, min_size = NA), "min_size")
  expect_error(mmc_gsea_fgsea(pw, ec, max_size = NA), "max_size")
  expect_error(mmc_gsea_fgsea(pw, c(E1 = 3, E1 = 1, E2 = -1)), "one score per EC id")

  cols <- names(mmc_gsea_fgsea(pw, ec))
  expect_identical(names(mmc_gsea_fgsea(list(), ec)), cols)
  expect_identical(names(mmc_gsea_fgsea(pw, ec, min_size = 3)), cols)   # filtered out
  expect_identical(nrow(mmc_gsea_fgsea(pw, ec, min_size = 3)), 0L)

  # ECs without a finite score are not ranked, as in the legacy method; a
  # pathway holding every ranked EC is not tested (fgsea's size cap).
  with_na <- c(ec, E5 = NA_real_)
  expect_equal(mmc_gsea_fgsea(pw, with_na)$ES, mmc_gsea_fgsea(pw, ec)$ES)
  expect_identical(nrow(mmc_gsea_fgsea(list(all = names(ec)), ec)), 0L)
})

test_that("the default method is the pipeline's shared fgsea::fgsea() multilevel call", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("BiocParallel")
  # One nPermSimple for every GSEA in the pipeline: the mummichog constant is
  # the shared wrapper's own default, so neither can drift from the other.
  expect_identical(.MMC_GSEA_NPERM_SIMPLE, 10000L)
  expect_identical(.MMC_GSEA_NPERM_SIMPLE, formals(run_fgsea_core)$n_perm_simple)

  fx  <- gsea_tied_fixture()
  res <- suppressWarnings(mmc_gsea_fgsea(fx$sets, fx$ec, seed = 7))

  # The same call made directly, on the documented ordering and seed. Every
  # statistic is fgsea's: nothing is recomputed around it.
  ordered <- fx$ec[order(-fx$ec, names(fx$ec), method = "radix")]
  want <- suppressWarnings(withr::with_seed(7, fgsea::fgsea(
    pathways = fx$sets, stats = ordered, minSize = 1L, maxSize = Inf,
    nPermSimple = 10000L, BPPARAM = BiocParallel::SerialParam())))
  want <- as.data.frame(want)
  # Multilevel output: an error estimate, and no permutation tally.
  expect_true("log2err" %in% names(want))
  expect_false("nMoreExtreme" %in% names(want))
  want <- want[match(names(fx$sets), want$pathway), , drop = FALSE]

  expect_identical(res$pathway, as.character(want$pathway))
  expect_equal(res$ES, want$ES)
  expect_equal(res$NES, want$NES)
  expect_equal(res$pval, want$pval)
  expect_equal(res$padj, want$padj)
  expect_equal(res$log2err, want$log2err)
  expect_equal(res$size, as.numeric(want$size))
  expect_identical(res$leadingEdge,
                   vapply(want$leadingEdge, paste, character(1),
                          collapse = "; ", USE.NAMES = FALSE))
})

test_that("the shared fgsea wrapper joins the leading edge with ',' unless asked", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("BiocParallel")
  ranks <- c(g1 = 3, g2 = 2, g3 = 1, g4 = -1, g5 = -2, g6 = -3)
  sets  <- list(top = c("g1", "g2"))
  run <- function(...) run_fgsea_core(ranks, sets, min_size = 1, seed = 1L,
                                      bpparam = BiocParallel::SerialParam(), ...)
  expect_identical(formals(run_fgsea_core)$leading_edge_sep, ",")
  expect_identical(run()$leadingEdge, "g1,g2")
  expect_identical(run(leading_edge_sep = "; ")$leadingEdge, "g1; g2")
})

test_that("gsea_permutations is the legacy method's setting and is refused otherwise", {
  # With the default method the key has nothing to set. It is not ignored and
  # not read as nPermSimple: the run stops and says what to do.
  expect_error(mmc_gsea_method(list(gsea_permutations = 100)),
               "gsea_permutations")
  expect_error(mmc_gsea_method(list(gsea_method = "fgsea",
                                    gsea_permutations = 10000)),
               "metaboanalyst_legacy")
  # The message does not offer the legacy method as a way back to a result of
  # the default method's former fgseaSimple calculation.
  msg <- tryCatch(mmc_gsea_method(list(gsea_permutations = 100)),
                  error = conditionMessage)
  expect_match(msg, "it does not reproduce a result", fixed = TRUE)
  expect_false(grepl("written for", msg, fixed = TRUE))
  # The legacy method keeps the key, with its meaning unchanged.
  expect_identical(mmc_gsea_method(list(gsea_method = "metaboanalyst_legacy",
                                        gsea_permutations = 200)),
                   "metaboanalyst_legacy")
  expect_identical(mmc_gsea_method(list()), "fgsea")
  expect_identical(.MMC_GSEA_DEFAULTS$n_perm, 100L)

  # The same rule for direct callers of the adapter.
  expect_error(run_mummichog_gsea(character(0),
                                  data.frame(feature_id = "f", logFC = 1),
                                  list(pathways = list()), n_perm = 100L),
               "n_perm is the legacy method")
})

test_that("a result's p-value estimation is described from its own settings", {
  # multilevel result
  expect_identical(mmc_gsea_permutation_label(list(nperm_simple = 10000L)),
                   "adaptive multilevel p-values, nPermSimple = 10000")
  # fixed-permutation result: the legacy method, or the earlier default
  expect_identical(mmc_gsea_permutation_label(list(n_perm = 100L)),
                   "100 permutations")
  # A field that merely starts with "n_perm" is not a permutation count.
  expect_identical(mmc_gsea_permutation_label(list(n_perm_other = 5L)),
                   "p-value estimation not recorded")
  # neither recorded: said, not guessed
  expect_identical(mmc_gsea_permutation_label(list(seed = 1L)),
                   "p-value estimation not recorded")
  expect_identical(mmc_gsea_permutation_label(NULL),
                   "p-value estimation not recorded")
})

test_that("fgsea's error estimate is summarised, and absent for the legacy method", {
  s <- mmc_gsea_log2err_summary(data.frame(pval = c(0.01, 0.2, 1e-50, NA),
                                           log2err = c(0.1, 0.3, NA, NA)))
  expect_equal(s$max, 0.3)
  # one pathway has a p-value without an estimate; the one without a p-value
  # is not counted here (it is counted as having no p-value)
  expect_identical(s$n_unavailable, 1L)

  none <- mmc_gsea_log2err_summary(data.frame(pval = c(0.5, NA),
                                              log2err = c(NA_real_, NA_real_)))
  expect_null(none$max)
  expect_identical(none$n_unavailable, 1L)

  # The legacy engine's table has no log2err column.
  expect_null(mmc_gsea_log2err_summary(
    mmc_gsea_metaboanalyst(list(), numeric(0), numeric(0))))
  expect_null(mmc_gsea_log2err_summary(NULL))
})


# ---------------------------------------------------------------------------
# the python adapter: default method, legacy method, missing p-values
# ---------------------------------------------------------------------------

test_that("the python adapter runs standard fgsea by default and records it", {
  skip_if_not_installed("fgsea")
  f   <- build_gsea_fixture(withr::local_tempdir())
  res <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                            de_method = "limma")

  expect_identical(res$method, "fgsea")
  expect_match(res$method_label, "fgsea::fgsea, adaptive multilevel", fixed = TRUE)
  expect_false(grepl("fgseaSimple", res$method_label, fixed = TRUE))
  expect_identical(res$reference$method, "fgsea::fgsea (fgseaMultilevel)")
  expect_identical(res$params$nperm_simple, 10000L)
  expect_null(res$params$n_perm)                 # no fixed permutation count
  expect_null(res$reference$commit)              # no MetaboAnalystR pin claimed
  expect_length(res$deviations, 0L)
  expect_identical(res$nes_note, mmc_gsea_nes_note("fgsea"))
  expect_null(res$stats)
  expect_null(res$tie_ranks)
  expect_identical(res$n_positions, res$n_ec)    # one ranked position per EC
  expect_identical(res$ties$n_ec, 24L)
  expect_identical(res$ties$n_tied_ec, 0L)
  expect_identical(res$n_p_unavailable, 0L)

  # It is mmc_gsea_fgsea() on the adapter's own EC scores and pathway sets.
  sets <- mmc_pathway_ec_sets(f$model$pathways,
                              read_mummichog_ec_candidates(f$files),
                              names(res$ranks))
  direct <- suppressWarnings(mmc_gsea_fgsea(sets, res$ranks))
  k <- match(res$table$Pathway, direct$pathway)
  expect_equal(res$table[["P.Value"]], direct$pval[k])
  expect_equal(res$table$ES, round(direct$ES[k], 4))

  # fgsea's error estimate travels beside the table, one value per tested
  # pathway, and is not a column of it: the table is what the exports and the
  # Shiny payload carry, and its columns are the same for both methods.
  expect_identical(names(res$log2err), direct$pathway)
  expect_equal(unname(res$log2err), direct$log2err)
  expect_identical(res$log2err_summary, mmc_gsea_log2err_summary(direct))
  expect_false("log2err" %in% names(res$table))

  # The fixture's three disjoint pathways hold the 8 highest, the 8 middle and
  # the 8 lowest scores, so the signed ranking gives +1 and -1 at the ends and
  # the leading edge names the pathway's own ECs.
  es <- stats::setNames(res$table$ES, res$table$Pathway)
  expect_equal(es[["PW_pos"]], 1)
  expect_equal(es[["PW_neg"]], -1)
  le <- stats::setNames(res$table[["Leading-edge EmpiricalCompounds"]],
                        res$table$Pathway)
  expect_identical(le[["PW_pos"]], paste(sprintf("E%03d", 1:8), collapse = "; "))
  expect_setequal(strsplit(le[["PW_neg"]], "; ", fixed = TRUE)[[1]],
                  sprintf("E%03d", 17:24))

  # The table schema is the one the legacy method writes.
  expect_true(all(c("ES defaulted by reference implementation",
                    "ES fallback reason") %in% names(res$table)))
  expect_false(any(res$table[["ES defaulted by reference implementation"]]))
})

test_that("the legacy method stays reachable, warns, and is recorded as legacy", {
  skip_if_not_installed("fgsea")
  f <- build_gsea_fixture(withr::local_tempdir())

  expect_warning(
    res <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                              de_method = "limma",
                              method = "metaboanalyst_legacy"),
    "known defects")
  expect_identical(res$method, "metaboanalyst_legacy")
  expect_match(res$method_label, "legacy MetaboAnalyst-compatible", fixed = TRUE)
  expect_identical(res$reference$commit, .MMC_GSEA_REF_COMMIT)
  expect_identical(res$nes_note, mmc_gsea_nes_note("metaboanalyst_legacy"))
  expect_null(res$ties)
  expect_identical(res$n_p_unavailable, 0L)
  # Its own setting and default, and nothing of the default method's.
  expect_identical(res$params$n_perm, 100L)
  expect_null(res$params$nperm_simple)
  expect_null(res$log2err)
  expect_null(res$log2err_summary)

  # It is the engine the parity fixture pins, on the adapter's own inputs.
  ri <- mmc_gsea_ranked_inputs(res$ranks)
  expect_identical(res$stats, ri$stats)
  # (values only: fastmatch leaves a lookup-table attribute on the names of
  # the vector the engine was given)
  expect_equal(as.numeric(res$tie_ranks), as.numeric(ri$ranks))
  expect_identical(res$n_positions, length(ri$stats))

  # The default run gives no such warning (fgsea's own, if any, are not
  # this test's subject), and both test the same pathways on the same EC
  # scores.
  w <- testthat::capture_warnings(
    std <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                              de_method = "limma"))
  expect_false(any(grepl("known defects", w, fixed = TRUE)))
  expect_identical(std$ranks, res$ranks)
  expect_setequal(std$table$Pathway, res$table$Pathway)
})

test_that("a tested pathway without a p-value is counted and kept, with NA", {
  skip_if_not_installed("fgsea")
  f <- build_gsea_fixture(withr::local_tempdir())

  # fgsea gives no p-value when too few permutations produced a score of the
  # observed sign. That cannot be arranged reliably from inputs, so the engine
  # is stubbed to return one such row, in the shape mmc_gsea_fgsea() returns.
  target <- environment(run_mummichog_gsea)
  orig   <- get("mmc_gsea_fgsea", envir = target)
  assign("mmc_gsea_fgsea", function(pathways, ec_scores, ...) {
    out <- data.frame(
      pathway = c("PW_pos", "PW_mid", "PW_neg"), pval = c(0.01, NA, 0.02),
      padj = c(0.03, NA, 0.03), ES = c(1, 0.2, -1), NES = c(2, NA, -2),
      nMoreExtreme = NA_real_, size = c(8, 8, 8),
      leadingEdge = c("E001", NA, "E024"), es_defaulted = FALSE,
      es_reason = NA_character_, log2err = c(0.1, NA, NA),
      stringsAsFactors = FALSE)
    attr(out, "mmc_ties") <- mmc_gsea_tie_summary(ec_scores, pathways)
    out
  }, envir = target)
  withr::defer(assign("mmc_gsea_fgsea", orig, envir = target))

  msgs <- capture_messages(
    res <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "c",
                              de_method = "limma"))
  expect_identical(res$n_p_unavailable, 1L)
  expect_true(any(grepl("have no p-value", msgs, fixed = TRUE)))
  # kept as a tested pathway, sorted after the pathways that have a p-value
  expect_identical(res$table$Pathway, c("PW_pos", "PW_neg", "PW_mid"))
  expect_true(is.na(res$table[["P.Value"]][3]))
  expect_identical(res$n_pathways, 3L)
  # PW_neg has a p-value without an error estimate; PW_mid has neither and is
  # counted once, as a pathway without a p-value.
  expect_equal(res$log2err_summary$max, 0.1)
  expect_identical(res$log2err_summary$n_unavailable, 1L)
})
