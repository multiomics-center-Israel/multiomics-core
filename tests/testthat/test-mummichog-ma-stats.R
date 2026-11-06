# tests/testthat/test-mummichog-ma-stats.R
#
# MetaboAnalyst-compatible mummichog, component 4: permutations and pathway
# statistics (06m). Unit tests on small synthetic EC sets, plus parity against
# the pinned MetaboAnalystR (fixtures/mummichog_ma_parity/
# generate_stats_reference.R), including the Gamma parameterization and the
# Holm adjustment as upstream actually outputs them. All data synthetic.

# Ten ECs E1..E10, one compound each; pathways A = E1-E4, B = E5-E7, C = E8-E10;
# significant ECs E1, E2, E5 (q = 3).
st_pathways <- list(
  cpds = list(paste0("C", 1:4), paste0("C", 5:7), paste0("C", 8:10)),
  name = c("A", "B", "C"),
  emp_cpds = list(paste0("E", 1:4), paste0("E", 5:7), paste0("E", 8:10)))
st_ecs <- list(
  input_ecpdlist = c("E1", "E2", "E5"),
  ecpd_cpd_dict = stats::setNames(as.list(paste0("C", 1:10)), paste0("E", 1:10)),
  total_matched_ecpds = paste0("E", 1:10),
  mz2ec_dict = stats::setNames(as.list(paste0("E", 1:10)),
                               as.character(100 + 1:10)))
st_prep <- list(ref_mzlist = as.character(100 + 1:20))
st_sig <- as.character(100 + c(1, 2, 5, 11, 12))     # N = 5

# Permutation records against the MetaboAnalystR fixture. The fixture was
# generated on one platform; phyper() goes through the C maths library
# (lgamma, log1p, exp), whose last bits can differ on another (CI is Linux),
# so the p-values are compared to a relative 1e-12 rather than bit for bit.
# Everything else stays exact: the record's structure (per-permutation list,
# matrix dimensions and attributes) and its NA pattern; the overlap counts
# (perm_hits) and all reported ORA / Holm outputs are compared with
# expect_identical() by the callers. The measured difference is printed so the
# CI log shows it.
expect_perm_record_parity <- function(got, ref, what) {
  shape <- function(rec) lapply(rec, function(x) lapply(x, attributes))
  expect_identical(length(got), length(ref), label = paste(what, "length"))
  expect_identical(shape(got), shape(ref), label = paste(what, "structure"))
  a <- unlist(got); b <- unlist(ref)
  expect_identical(is.na(a), is.na(b), label = paste(what, "NA pattern"))
  ok <- !is.na(a) & !is.na(b)
  rel <- if (any(ok)) max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), .Machine$double.xmin)) else 0
  cat(sprintf("\n[perm_record parity] %s: %d value(s), %d differ, max rel diff %.3g (%s, %s)\n",
              what, sum(ok), sum(a[ok] != b[ok]), rel, R.version.string,
              Sys.info()[["sysname"]]))
  expect_lte(rel, 1e-12, label = paste(what, "max relative difference"))
}

# ---- permutation p-values --------------------------------------------------------

test_that("a permutation's p-value is upstream's unclamped phyper", {
  out <- .mmc_ma_perm_pvals(c("E1", "E2", "E5"), paste0("E", 1:10), st_pathways)
  # A: u 2, s 4, nn = 10 + 2 - 4 - 3 = 5; B: u 1, s 3, nn = 5; C: u 0, s 3, nn = 4
  expect_identical(out[[1]][, 1], c(
    stats::phyper(1, 4, 3 + 5, 3, lower.tail = FALSE),
    stats::phyper(0, 3, 3 + 5, 3, lower.tail = FALSE),
    stats::phyper(-1, 3, 3 + 4, 3, lower.tail = FALSE)))
  expect_identical(out[[2]], as.matrix(c(2L, 1L, 0L)))
})

test_that("permutation i draws N m/z from the whole list under seed 123 + i", {
  set.seed(5)
  before <- .Random.seed
  perms <- mmc_ma_permutations(st_prep, st_ecs, st_pathways, st_sig, n_perm = 3)
  expect_identical(.Random.seed, before)
  expect_length(perms$perm_record, 3L)
  for (i in 1:3) {
    draw <- withr::with_seed(123 + i, sample(st_prep$ref_mzlist, 5))
    hit <- unique(unlist(st_ecs$mz2ec_dict[draw]))
    exp <- .mmc_ma_perm_pvals(hit, st_ecs$total_matched_ecpds, st_pathways)
    expect_identical(perms$perm_record[[i]], exp[1])
    expect_identical(perms$perm_hits[[i]], exp[2])
  }
  expect_identical(perms, mmc_ma_permutations(st_prep, st_ecs, st_pathways,
                                              st_sig, n_perm = 3))
  expect_error(mmc_ma_permutations(st_prep, st_ecs, st_pathways, st_sig,
                                   n_perm = 0), "whole number")
})

# ---- observed statistics ---------------------------------------------------------

test_that("FET and EASE follow upstream's contingency tables", {
  skip_if_not_installed("fitdistrplus")
  perms <- mmc_ma_permutations(st_prep, st_ecs, st_pathways, st_sig, n_perm = 20)
  st <- suppressWarnings(mmc_ma_pathway_stats(st_ecs, st_pathways, perms))
  pv <- st$pvals
  # observed: nn = max(0, T + u - s - q); FET has "- u", EASE uses max(0, u - 1)
  expect_equal(pv["A", "FET"], stats::phyper(1, 4, 3 + 5 - 2, 3, lower.tail = FALSE))
  expect_equal(pv["B", "FET"], stats::phyper(0, 3, 3 + 5 - 1, 3, lower.tail = FALSE))
  expect_equal(pv["A", "EASE"], stats::fisher.test(
    matrix(c(1, 3, 1, 5), nrow = 2), alternative = "greater")$p.value)
  expect_equal(pv["B", "EASE"], stats::fisher.test(
    matrix(c(0, 3, 2, 5), nrow = 2), alternative = "greater")$p.value)
  expect_equal(unname(pv[, "Expected"]), 3 * c(4, 3, 3) / 10)
  expect_identical(unname(pv[, "Sig.EC"]), c(2, 1, 0))
  expect_identical(unname(pv[, "Hits.sig"]), c(2, 1, 0))
})

test_that("the result drops pathways without hits, rounds, and adjusts by Holm", {
  skip_if_not_installed("fitdistrplus")
  perms <- mmc_ma_permutations(st_prep, st_ecs, st_pathways, st_sig, n_perm = 20)
  st <- suppressWarnings(mmc_ma_pathway_stats(st_ecs, st_pathways, perms))
  expect_setequal(rownames(st$resmat), c("A", "B"))
  expect_equal(st$resmat[, "FET"], signif(st$pvals[rownames(st$resmat), "FET"], 5),
               tolerance = 1e-12)
  expect_identical(st$table$adj_fet,
                   unname(stats::p.adjust(st$resmat[, "FET"], method = "holm")))
  expect_identical(st$table$adj_ease,
                   unname(stats::p.adjust(st$resmat[, "EASE"], method = "holm")))
  expect_identical(colnames(st$restable)[c(5, 6, 7, 10, 11, 12, 13, 14)],
                   c("P(Fisher)", "P(EASE)", "P(Gamma)", "AdjP.Fisher",
                     "AdjP.EASE", "AdjP.Gamma", "Pathway Number", "cpd.hits"))

  # empirical p: share of permutation p-values <= the observed FET
  rec <- do.call(cbind, do.call(cbind, perms$perm_record))
  emp_a <- mean(rec[1, ] <= st$pvals["A", "FET"])
  expect_equal(st$resmat["A", "Empirical"], signif(emp_a, 5), tolerance = 1e-12)
  expect_identical(st$resmat["A", "Emp.Hits"], sum(rec[1, ] <= st$pvals["A", "FET"]) + 0)
})

test_that("Gamma passes the scale estimate as the rate, on FET", {
  skip_if_not_installed("fitdistrplus")
  perms <- mmc_ma_permutations(st_prep, st_ecs, st_pathways, st_sig, n_perm = 20)
  st <- suppressWarnings(mmc_ma_pathway_stats(st_ecs, st_pathways, perms))
  skip_if_not(st$gamma$fitted, "Gamma fit did not converge on this tiny null")
  est <- st$gamma$estimate
  fet <- st$pvals[rownames(st$resmat), "FET"]
  expect_equal(unname(st$resmat[, "Gamma"]),
               unname(signif(1 - stats::pgamma(1 - fet, shape = est["shape"],
                                               rate = est["scale"]), 5)),
               tolerance = 1e-12)
})

test_that("one pathway with hits or fewer stops, as upstream does", {
  skip_if_not_installed("fitdistrplus")
  ecs1 <- st_ecs
  ecs1$input_ecpdlist <- c("E1", "E2")
  perms <- mmc_ma_permutations(st_prep, ecs1, st_pathways, st_sig[1:2], n_perm = 5)
  expect_error(suppressWarnings(mmc_ma_pathway_stats(ecs1, st_pathways, perms)),
               "only 1 pathway")
})

# ---- parity with the pinned MetaboAnalystR ---------------------------------------

ma_stats_ref <- function() {
  ref_file <- file.path(ma_fixture_dir, "stats_reference.rds")
  skip_if_not(file.exists(ref_file), "stats parity reference not generated yet")
  ref <- readRDS(ref_file)
  expect_identical(ref$metaboanalystr_commit,
                   "398476ae2a0c996e390925fa62adc46b5ecde334")
  skip_if_not(identical(Sys.getlocale("LC_COLLATE"), ref$locale_collate),
              paste0("stats reference generated under LC_COLLATE '",
                     ref$locale_collate, "'"))
  skip_if_not_installed("fitdistrplus")
  skip_if_not(identical(as.character(utils::packageVersion("fitdistrplus")),
                        ref$fitdistrplus_version),
              paste0("stats reference generated under fitdistrplus ",
                     ref$fitdistrplus_version))
  ref
}

ma_stats_run <- function(ref, up) {
  inp <- ref$inputs
  lib <- mmc_ma_library_from_model(ma_write_model(
    inp$compounds, inp$pathways, pathway_names = inp$pathway_names))
  prep <- mmc_ma_prepare_peaks(inp$peaks, up$mode, rt_frac = inp$rt_frac)
  matched <- mmc_ma_match_peaks(prep, lib, ppm = inp$ppm)
  sig <- mmc_ma_input_mzlist(prep, inp$cutoff)
  ecs <- mmc_ma_build_ecs(matched, prep, sig,
                          force_primary_ion = identical(up$primary_ion, "yes"))
  pw <- mmc_ma_pathway_ecs(lib, ecs, min_lib = ref$min_lib)
  perms <- mmc_ma_permutations(prep, ecs, pw, sig, n_perm = ref$perm_num)
  list(perms = perms, stats = mmc_ma_pathway_stats(ecs, pw, perms))
}

for (case in c("mixed", "positive", "negative", "mixed_no_primary")) {
  test_that(sprintf("pathway statistics match MetaboAnalystR (%s)", case), {
    ref <- ma_stats_ref()
    up <- ref$cases[[case]]
    got <- ma_stats_run(ref, up)
    st <- got$stats

    expect_perm_record_parity(got$perms$perm_record, up$perm_record, case)
    expect_identical(got$perms$perm_hits, up$perm_hits)
    expect_identical(st$pvals, up$pvals)
    expect_identical(st$resmat[, -11, drop = FALSE], up$mummi.resmat)
    expect_identical(rownames(st$resmat), up$path.nms)
    expect_identical(unname(st$resmat[, "FET"]), up$path.pval)
    expect_identical(st$restable, up$restable)
  })

  test_that(sprintf("upstream's Gamma and Holm are what its output shows (%s)", case), {
    ref <- ma_stats_ref()
    up <- ref$cases[[case]]
    st <- ma_stats_run(ref, up)$stats
    res <- up$mummi.resmat
    tab <- up$restable

    # Holm, not BH: the AdjP.* columns are p.adjust(<rounded column>, "holm").
    for (pair in list(c("FET", "AdjP.Fisher"), c("EASE", "AdjP.EASE"),
                      c("Gamma", "AdjP.Gamma"))) {
      expect_equal(suppressWarnings(as.numeric(tab[, pair[2]])),
                   unname(stats::p.adjust(res[, pair[1]], method = "holm")),
                   tolerance = 1e-12, label = pair[2])
    }

    # Gamma: 1 - pgamma(1 - FET, shape, rate = <scale estimate>).
    if (all(is.na(res[, "Gamma"]))) {
      expect_false(st$gamma$fitted)
    } else {
      expect_true(st$gamma$fitted)
      est <- st$gamma$estimate
      fet <- st$pvals[rownames(res), "FET"]
      as_rate  <- signif(1 - stats::pgamma(1 - fet, shape = est["shape"],
                                           rate = est["scale"]), 5)
      expect_equal(unname(res[, "Gamma"]), unname(as_rate), tolerance = 1e-12)
      # This fixture reports FET = 1 for every pathway, where rate and scale
      # (and Holm and BH) coincide; the targeted fixture below tells them apart.
    }
  })
}

test_that("the end-to-end stats fixture fits Gamma at least once", {
  ref <- ma_stats_ref()
  expect_true(any(vapply(ref$cases, function(up)
    any(!is.na(up$mummi.resmat[, "Gamma"])), logical(1))))
})

# ---- targeted fixture: informative p-values ----------------------------------------
# A hand-built EC universe (fixtures/mummichog_ma_parity/
# generate_stats_targeted_reference.R) whose pathway p-values separate Gamma's
# rate from scale and Holm from BH, run through the pinned upstream
# permutations and statistics.

ma_targeted_ref <- function() {
  ref_file <- file.path(ma_fixture_dir, "stats_targeted_reference.rds")
  skip_if_not(file.exists(ref_file), "targeted stats reference not generated yet")
  ref <- readRDS(ref_file)
  expect_identical(ref$metaboanalystr_commit,
                   "398476ae2a0c996e390925fa62adc46b5ecde334")
  skip_if_not(identical(Sys.getlocale("LC_COLLATE"), ref$locale_collate),
              paste0("targeted reference generated under LC_COLLATE '",
                     ref$locale_collate, "'"))
  skip_if_not_installed("fitdistrplus")
  skip_if_not(identical(as.character(utils::packageVersion("fitdistrplus")),
                        ref$fitdistrplus_version),
              paste0("targeted reference generated under fitdistrplus ",
                     ref$fitdistrplus_version))
  ref
}

ma_targeted_run <- function(ref) {
  inp <- ref$inputs
  ecs <- inp[c("input_ecpdlist", "total_matched_ecpds", "mz2ec_dict",
               "ecpd_cpd_dict")]
  perms <- mmc_ma_permutations(list(ref_mzlist = inp$ref_mzlist), ecs,
                               inp$pathways, inp$input_mzlist,
                               n_perm = ref$perm_num)
  # fitdistrplus warns while estimating standard errors; see the test below.
  st <- suppressWarnings(mmc_ma_pathway_stats(ecs, inp$pathways, perms))
  list(perms = perms, stats = st)
}

test_that("targeted: pathway statistics match MetaboAnalystR exactly", {
  ref <- ma_targeted_ref()
  up <- ref$upstream
  got <- ma_targeted_run(ref)
  st <- got$stats
  expect_perm_record_parity(got$perms$perm_record, up$perm_record, "targeted")
  expect_identical(got$perms$perm_hits, up$perm_hits)
  expect_identical(st$pvals, up$pvals)
  expect_identical(st$resmat[, -11, drop = FALSE], up$mummi.resmat)
  expect_identical(rownames(st$resmat), up$path.nms)
  expect_identical(unname(st$resmat[, "FET"]), up$path.pval)
  expect_identical(st$restable, up$restable)
})

test_that("targeted: the fixture is informative", {
  ref <- ma_targeted_ref()
  res <- ref$upstream$mummi.resmat
  expect_gte(nrow(res), 3L)
  expect_true(any(res[, "FET"] > 0 & res[, "FET"] < 0.05))
  expect_true(all(is.finite(res[, "Gamma"])))
  expect_true(any(res[, "Gamma"] > 0 & res[, "Gamma"] < 1))
})

test_that("targeted: upstream passes the Gamma scale estimate as the rate", {
  ref <- ma_targeted_ref()
  st <- ma_targeted_run(ref)$stats
  res <- ref$upstream$mummi.resmat
  est <- st$gamma$estimate
  fet <- st$pvals[rownames(res), "FET"]
  as_rate  <- signif(1 - stats::pgamma(1 - fet, shape = est["shape"],
                                       rate = est["scale"]), 5)
  as_scale <- signif(1 - stats::pgamma(1 - fet, shape = est["shape"],
                                       scale = est["scale"]), 5)
  expect_false(isTRUE(all.equal(unname(as_rate), unname(as_scale))))
  expect_equal(unname(res[, "Gamma"]), unname(as_rate), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(unname(res[, "Gamma"]), unname(as_scale))))
})

test_that("targeted: upstream adjusts by Holm, not BH", {
  ref <- ma_targeted_ref()
  res <- ref$upstream$mummi.resmat
  tab <- ref$upstream$restable
  for (pair in list(c("FET", "AdjP.Fisher"), c("EASE", "AdjP.EASE"),
                    c("Gamma", "AdjP.Gamma"))) {
    adj  <- as.numeric(tab[, pair[2]])
    holm <- unname(stats::p.adjust(res[, pair[1]], method = "holm"))
    bh   <- unname(stats::p.adjust(res[, pair[1]], method = "BH"))
    expect_false(isTRUE(all.equal(holm, bh)), label = paste(pair[1], "Holm vs BH"))
    expect_equal(adj, holm, tolerance = 1e-12, label = pair[2])
    expect_false(isTRUE(all.equal(adj, bh)), label = paste(pair[2], "is not BH"))
  }
})

test_that("targeted: Gamma uses the fitted estimates; the fit warnings concern only their SEs", {
  ref <- ma_targeted_ref()
  st <- ma_targeted_run(ref)$stats
  dg <- ref$diagnostics
  # the optimiser converged, and the estimates are finite and positive ...
  expect_identical(dg$convergence, 0L)
  expect_true(all(is.finite(dg$estimate) & dg$estimate > 0))
  # ... our fit is the same call on the same null, estimate for estimate ...
  expect_true(st$gamma$fitted)
  expect_identical(st$gamma$estimate, dg$estimate)
  # ... and the Gamma column is a function of those estimates alone (already
  # shown to equal upstream's). Any warnings come from fitdistrplus inverting
  # the Hessian for the standard errors (cov2cor / sqrt(diag(varcovar))),
  # which neither upstream nor this engine reads.
  if (length(dg$warnings) > 0) {
    expect_true(all(grepl("cov2cor|sqrt|NaN", dg$warnings)))
  }
})
