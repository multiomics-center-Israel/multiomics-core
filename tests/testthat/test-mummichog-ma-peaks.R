# tests/testthat/test-mummichog-ma-peaks.R
#
# MetaboAnalyst-compatible mummichog, component 2: peak-list preparation,
# significant-peak selection and polarity-specific m/z matching (06k). Unit
# tests on small synthetic inputs, plus parity against the pinned MetaboAnalystR
# (fixtures/mummichog_ma_parity/generate_peaks_reference.R) in all three modes.
# All data synthetic.

PROTON <- 1.00727646677

ma_peaks <- function(feature_id, mz, rt, p_value, statistic = 1,
                     polarity = "positive") {
  data.frame(feature_id = feature_id, mz = mz, rt = rt, p_value = p_value,
             statistic = statistic, polarity = polarity,
             stringsAsFactors = FALSE)
}

# ---- mmc_ma_prepare_peaks ------------------------------------------------------

test_that("peaks are ordered by p, ties kept in input order", {
  pk <- ma_peaks(c("A", "B", "C", "D"), c(200.1, 300.2, 400.3, 500.4),
                 c(1, 2, 3, 4), c(0.5, 0.01, 0.2, 0.01))
  prep <- mmc_ma_prepare_peaks(pk, "positive", min_features = 1)
  expect_identical(prep$table$feature_id, c("B", "D", "C", "A"))
  expect_identical(prep$ref_mzlist, c("300.2", "500.4", "400.3", "200.1"))
  expect_identical(prep$expr_dic, stats::setNames(c(1, 1, 1, 1), prep$ref_mzlist))
})

test_that("m/z outside 50 < m/z < 2001 is trimmed", {
  pk <- ma_peaks(c("lo", "edge_lo", "ok", "edge_hi", "hi"),
                 c(49.9, 50, 1000, 2000.9, 2001), 1:5, rep(0.1, 5))
  prep <- mmc_ma_prepare_peaks(pk, "positive", min_features = 1)
  expect_identical(prep$table$feature_id, c("ok", "edge_hi"))
  expect_identical(prep$removed[["out_of_range"]], 3L)
})

test_that("an m/z shared across polarities keeps only the smaller p", {
  pk <- ma_peaks(c("pos", "neg", "other"), c(400.123456, 400.123456, 500),
                 c(9, 9.4, 3), c(0.03, 0.02, 0.5),
                 polarity = c("positive", "negative", "positive"))
  prep <- mmc_ma_prepare_peaks(pk, "mixed", min_features = 1)
  expect_identical(prep$table$feature_id, c("neg", "other"))
  expect_identical(prep$pos_inx, c(FALSE, TRUE))
  expect_identical(prep$removed[["duplicate_mz"]], 1L)
})

test_that("m/z printing alike at 15 digits are duplicates, the later dropped", {
  pk <- ma_peaks(c("a", "b", "c"), c(300.1234567890121, 300.1234567890123, 400),
                 1:3, c(0.1, 0.2, 0.3))
  prep <- mmc_ma_prepare_peaks(pk, "negative", min_features = 1)
  expect_identical(prep$table$feature_id, c("a", "c"))
  expect_identical(prep$ref_mzlist, c("300.123456789012", "400"))
  expect_identical(prep$removed[["duplicate_mz"]], 1L)
})

test_that("RT ranks are reproducible and leave the caller's RNG untouched", {
  pk <- ma_peaks(paste0("F", 1:6), 100 + 1:6, c(5, 5, 5, 1, 2, 2), rep(0.1, 6))
  set.seed(99)
  before <- .Random.seed
  r1 <- mmc_ma_prepare_peaks(pk, "positive", min_features = 1)$rt_rank
  expect_identical(.Random.seed, before)
  r2 <- mmc_ma_prepare_peaks(pk, "positive", min_features = 1)$rt_rank
  expect_identical(r1, r2)
  expect_setequal(unname(r1), 1:6)
  expect_identical(unname(r1[4]), 1L)      # the one unique smallest RT
})

test_that("RT tolerance defaults to max(rt) * rt_frac and honours a fixed value", {
  pk <- ma_peaks(paste0("F", 1:3), c(110, 120, 130), c(2, 10, 30), rep(0.1, 3))
  expect_identical(mmc_ma_prepare_peaks(pk, "positive", min_features = 1)$rt_tol,
                   30 * 0.02)
  expect_identical(mmc_ma_prepare_peaks(pk, "positive", rt_frac = 0.05,
                                        min_features = 1)$rt_tol, 30 * 0.05)
  expect_identical(mmc_ma_prepare_peaks(pk, "positive", rt_tol = 0.3,
                                        min_features = 1)$rt_tol, 0.3)
})

test_that("single-polarity modes ignore the polarity column", {
  pk <- ma_peaks(c("a", "b"), c(110, 120), 1:2, c(0.1, 0.2),
                 polarity = c("negative", "bogus"))
  expect_identical(mmc_ma_prepare_peaks(pk, "positive", min_features = 1)$pos_inx,
                   c(TRUE, TRUE))
  expect_identical(mmc_ma_prepare_peaks(pk, "negative", min_features = 1)$pos_inx,
                   c(FALSE, FALSE))
})

test_that("the peak list refuses what MetaboAnalyst would silently reinterpret", {
  ok <- ma_peaks(c("a", "b"), c(110, 120), 1:2, c(0.1, 0.2),
                 polarity = c("positive", "negative"))
  bad_pol <- ok
  bad_pol$polarity[2] <- "Negative"
  expect_error(mmc_ma_prepare_peaks(bad_pol, "mixed", min_features = 1),
               "polarity other than")
  expect_error(mmc_ma_prepare_peaks(ok[, names(ok) != "polarity"], "mixed",
                                    min_features = 1), "polarity")
  expect_error(mmc_ma_prepare_peaks(ok[, names(ok) != "rt"], "mixed",
                                    min_features = 1), "Retention time is required")
  na_rt <- ok
  na_rt$rt[1] <- NA
  expect_error(mmc_ma_prepare_peaks(na_rt, "mixed", min_features = 1),
               "no usable 'rt'")
  expect_error(mmc_ma_prepare_peaks(ok, "mixed"), "under 100 features")
})

# ---- mmc_ma_input_mzlist --------------------------------------------------------

test_that("significant peaks use one inclusive cutoff over both polarities", {
  n <- 40
  pk <- ma_peaks(paste0("F", 1:n), 100 + 1:n, 1:n,
                 c(rep(0.05, 6), rep(0.01, 6), rep(0.5, n - 12)),
                 polarity = rep(c("positive", "negative"), length.out = n))
  prep <- mmc_ma_prepare_peaks(pk, "mixed", min_features = 1)
  sig <- mmc_ma_input_mzlist(prep, 0.05)
  expect_length(sig, 12L)
  expect_identical(sig, prep$ref_mzlist[1:12])
  expect_error(mmc_ma_input_mzlist(prep, 0.01), "fewer than 11")
  expect_error(mmc_ma_input_mzlist(prep, 0.001), "no feature")
  expect_error(mmc_ma_input_mzlist(prep, NA_real_), "single number")
})

# ---- mmc_ma_match_peaks ----------------------------------------------------------

ma_match_lib <- function() {
  cpd <- data.frame(id = c("C90001", "C00002"), name = c("syn A", "ATP"),
                    mw = c(180.063388, 506.995745), stringsAsFactors = FALSE)
  mmc_ma_library_from_model(ma_write_model(cpd, list(syn01 = cpd$id)))
}

test_that("each feature is matched against its own polarity's ion forms", {
  lib <- ma_match_lib()
  pk <- ma_peaks(c("pos", "neg"), c(180.063388 + PROTON, 180.063388 - PROTON),
                 c(4, 4.05), c(0.01, 0.02), polarity = c("positive", "negative"))
  mixed <- mmc_ma_match_peaks(mmc_ma_prepare_peaks(pk, "mixed", min_features = 1),
                              lib, ppm = 10)
  expect_identical(names(mixed), c("Query.Mass", "Matched.Compound",
                                   "Matched.Form", "Mass.Diff",
                                   "Retention.Time", "RT.Rank"))
  expect_true(all(vapply(mixed, is.character, logical(1))))
  expect_identical(mixed$Matched.Form, c("M+H[1+]", "M-H[-]"))   # positive first
  expect_identical(mixed$Matched.Compound, c("C90001", "C90001"))

  # Read as all-positive, the [M-H]- feature finds no positive ion form.
  pos_only <- mmc_ma_match_peaks(
    mmc_ma_prepare_peaks(pk, "positive", min_features = 1), lib, ppm = 10)
  expect_identical(pos_only$Matched.Form, "M+H[1+]")
})

test_that("the ppm tolerance is ppm * 1e-6 * m/z, strict", {
  lib <- ma_match_lib()
  mh <- 180.063388 + PROTON
  pk <- ma_peaks(c("in", "out"), c(mh * (1 + 8e-6), mh * (1 + 15e-6)), 1:2,
                 c(0.01, 0.02))
  m <- mmc_ma_match_peaks(mmc_ma_prepare_peaks(pk, "positive", min_features = 1),
                          lib, ppm = 10)
  expect_identical(m$Query.Mass, as.character(mh * (1 + 8e-6)))
  expect_error(mmc_ma_match_peaks(
    mmc_ma_prepare_peaks(pk[2, ], "positive", min_features = 1), lib, ppm = 10),
    "no feature matched")
  expect_error(mmc_ma_match_peaks(
    mmc_ma_prepare_peaks(pk, "positive", min_features = 1), lib, ppm = 0), "ppm")
})

test_that("currency compounds are dropped unless they are every match", {
  lib <- ma_match_lib()
  pk <- ma_peaks(c("syn", "atp"), c(180.063388 + PROTON, 506.995745 - PROTON),
                 1:2, c(0.01, 0.02), polarity = c("positive", "negative"))
  prep <- mmc_ma_prepare_peaks(pk, "mixed", min_features = 1)
  expect_false("C00002" %in% mmc_ma_match_peaks(prep, lib, 10)$Matched.Compound)
  expect_true("C00002" %in% mmc_ma_match_peaks(prep, lib, 10,
              exclude_currency = FALSE)$Matched.Compound)
  only_atp <- mmc_ma_prepare_peaks(pk[2, ], "mixed", min_features = 1)
  expect_identical(mmc_ma_match_peaks(only_atp, lib, 10)$Matched.Compound, "C00002")
})

# ---- parity with the pinned MetaboAnalystR ---------------------------------------

ma_peaks_ref <- function() {
  ref_file <- file.path(ma_fixture_dir, "peaks_reference.rds")
  skip_if_not(file.exists(ref_file), "peaks parity reference not generated yet")
  ref <- readRDS(ref_file)
  expect_identical(ref$metaboanalystr_commit,
                   "398476ae2a0c996e390925fa62adc46b5ecde334")
  ref
}

ma_ref_library <- function(inp) {
  pws <- stats::setNames(inp$pathways, names(inp$pathways))
  mmc_ma_library_from_model(ma_write_model(inp$compounds, pws,
                                           pathway_names = inp$pathway_names))
}

for (md in c("mixed", "positive", "negative")) {
  test_that(sprintf("peak preparation matches MetaboAnalystR (%s mode)", md), {
    ref <- ma_peaks_ref()
    up <- ref$modes[[md]]
    inp <- ref$inputs
    prep <- mmc_ma_prepare_peaks(inp$peaks, md, rt_frac = inp$rt_frac)

    expect_identical(prep$ref_mzlist, up$ref_mzlist)
    expect_identical(prep$pos_inx, up$pos_inx)
    expect_identical(prep$ret_time, up$ret_time)
    expect_identical(prep$expr_dic, up$expr_dic)
    expect_identical(prep$rt_tol, up$rt_tol)
    expect_identical(prep$rt_rank, up$ord_rt)
    proc <- up$mummi.proc
    expect_identical(prep$table$p_value, unname(as.numeric(proc$p.value)))
    expect_identical(prep$table$mz, unname(as.character(proc$m.z)))
    expect_identical(prep$table$statistic, unname(as.numeric(proc$t.score)))
    expect_identical(prep$table$rt, unname(as.numeric(proc$r.t)))
    expect_identical(mmc_ma_input_mzlist(prep, inp$cutoff), up$input_mzlist)
  })

  test_that(sprintf("m/z matching matches MetaboAnalystR (%s mode)", md), {
    ref <- ma_peaks_ref()
    up <- ref$modes[[md]]
    inp <- ref$inputs
    prep <- mmc_ma_prepare_peaks(inp$peaks, md, rt_frac = inp$rt_frac)
    matched <- mmc_ma_match_peaks(prep, ma_ref_library(inp), ppm = inp$ppm)
    expect_identical(matched, up$matched_res)
  })
}
