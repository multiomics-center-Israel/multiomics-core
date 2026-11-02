# tests/testthat/test-mummichog-ma-ecs.R
#
# MetaboAnalyst-compatible mummichog, component 3: empirical compounds and
# their pathway sets (06l). Unit tests on small synthetic inputs, plus parity
# against the pinned MetaboAnalystR (fixtures/mummichog_ma_parity/
# generate_ecs_reference.R). All data synthetic.

PROTON <- 1.00727646677
MW_A <- 180.063388   # C90001
MW_B <- 146.105528   # C90002

ma_ec_lib <- function() {
  cpd <- data.frame(id = sprintf("C9000%d", 1:5),
                    name = sprintf("syn %d", 1:5),
                    mw = c(MW_A, MW_B, 300.9995, 400.5, 500.25),
                    stringsAsFactors = FALSE)
  mmc_ma_library_from_model(ma_write_model(cpd, list(
    syn01 = c("C90001", "C90002", "C90003"),
    syn02 = c("C90001", "C90002"),                 # 2 compounds: < minLib
    syn03 = c("C90003", "C90004", "C90005"))))     # no detected compound
}

# Run prepare -> match -> ECs on a small peak list.
ma_ecs_for <- function(pk, mode = "mixed", rt_frac = 0.02, rt_tol = NA,
                       force_primary_ion = TRUE, sig = NULL) {
  prep <- mmc_ma_prepare_peaks(pk, mode, rt_frac = rt_frac, rt_tol = rt_tol,
                               min_features = 1)
  matched <- mmc_ma_match_peaks(prep, ma_ec_lib(), ppm = 10)
  mmc_ma_build_ecs(matched, prep, sig %||% prep$ref_mzlist[1],
                   force_primary_ion = force_primary_ion)
}

ec_of <- function(ecs, form) {
  unique(ecs$matched_res$Empirical.Compound[ecs$matched_res$Matched.Form == form])
}

pk_both <- data.frame(
  feature_id = c("A_pos", "A_neg", "B_pos"),
  mz = c(MW_A + PROTON, MW_A - PROTON, MW_B + PROTON),
  rt = c(4, 4.05, 20), p_value = c(0.01, 0.02, 0.5), statistic = c(2, -3, 1),
  polarity = c("positive", "negative", "positive"), stringsAsFactors = FALSE)

# ---- grouping across polarities -------------------------------------------------

test_that("a positive and a negative ion of one compound at one RT form one EC", {
  ecs <- ma_ecs_for(pk_both)
  a <- ecs$matched_res[ecs$matched_res$Matched.Compound == "C90001", ]
  expect_setequal(a$Matched.Form, c("M+H[1+]", "M-H[-]"))
  expect_length(unique(a$Empirical.Compound), 1L)
  expect_length(ecs$total_matched_ecpds, 2L)            # C90001's EC, C90002's
})

test_that("read as one polarity, the other polarity's ion joins no EC", {
  ecs <- ma_ecs_for(pk_both, mode = "positive")
  expect_false("M-H[-]" %in% ecs$matched_res$Matched.Form)
})

pk_far <- data.frame(
  feature_id = c("A_early", "B_mid", "A_late"),
  mz = c(MW_A + PROTON, MW_B + PROTON, MW_A + PROTON + 1e-6),
  rt = c(2, 15, 28), p_value = c(0.01, 0.02, 0.03), statistic = 1,
  polarity = "positive", stringsAsFactors = FALSE)

test_that("one compound at distant retention times gives two ECs", {
  ecs <- ma_ecs_for(pk_far)
  a <- ecs$matched_res[ecs$matched_res$Matched.Compound == "C90001", ]
  expect_identical(nrow(a), 2L)
  expect_length(unique(a$Empirical.Compound), 2L)
})

test_that("an RT cut needs both the RT gap and the RT-rank gap", {
  # RT gap 26 > rt_tol 1, but rank gap 2 <= 3 features * 0.9 = 2.7: no cut.
  ecs <- ma_ecs_for(pk_far, rt_frac = 0.9, rt_tol = 1)
  a <- ecs$matched_res[ecs$matched_res$Matched.Compound == "C90001", ]
  expect_length(unique(a$Empirical.Compound), 1L)
})

# ---- primary-ion rule ----------------------------------------------------------

pk_rescue <- data.frame(
  feature_id = c("A_k", "A_neg", "B_pos"),
  mz = c(MW_A + 37.9555 + PROTON, MW_A - PROTON, MW_B + PROTON),
  rt = c(4, 4.05, 20), p_value = c(0.01, 0.02, 0.5), statistic = 1,
  polarity = c("positive", "negative", "positive"), stringsAsFactors = FALSE)

test_that("a negative primary ion keeps an EC whose positive ion is not primary", {
  ecs <- ma_ecs_for(pk_rescue)
  a <- ecs$matched_res[ecs$matched_res$Matched.Compound == "C90001", ]
  expect_setequal(a$Matched.Form, c("M+K[1+]", "M-H[-]"))
  expect_length(unique(a$Empirical.Compound), 1L)
})

test_that("without a primary ion the EC is dropped, unless the rule is off", {
  k_only <- pk_rescue[pk_rescue$feature_id != "A_neg", ]
  expect_false("C90001" %in% ma_ecs_for(k_only)$matched_res$Matched.Compound)
  off <- ma_ecs_for(k_only, force_primary_ion = FALSE)
  expect_true("M+K[1+]" %in% off$matched_res$Matched.Form)
})

test_that("an EC set with no primary ion at all stops with a clear message", {
  k_only <- pk_rescue[pk_rescue$feature_id == "A_k", ]
  expect_error(ma_ecs_for(k_only), "no empirical compound contains a primary ion")
})

# ---- dictionaries and outputs ---------------------------------------------------

test_that("EC statistic is the max of its members; significant ECs follow the m/z", {
  ecs <- ma_ecs_for(pk_both, sig = as.character(MW_B + PROTON))
  ec_a <- ec_of(ecs, "M+H[1+]")
  expect_identical(unname(ecs$ec_exp[ec_a]), 2)
  ec_b <- unique(ecs$matched_res$Empirical.Compound[
    ecs$matched_res$Matched.Compound == "C90002"])
  expect_identical(ecs$input_ecpdlist, ec_b)
  expect_setequal(ecs$total_matched_ecpds, c(ec_a, ec_b))
  expect_setequal(names(ecs$mz2ec_dict), unique(ecs$matched_res$Query.Mass))
  expect_named(ecs$matched_res, c("Query.Mass", "Matched.Compound",
                                  "Matched.Form", "Retention.Time",
                                  "Mass.Diff", "Empirical.Compound"))
})

test_that("ECs are labelled EC000<n> in first-appearance order", {
  ecs <- ma_ecs_for(pk_both)
  expect_identical(unique(ecs$matched_res$Empirical.Compound),
                   paste0("EC000", seq_along(ecs$total_matched_ecpds)))
})

test_that("pathways keep >= min_lib compounds and at least one EC", {
  ecs <- ma_ecs_for(pk_both)
  pw <- mmc_ma_pathway_ecs(ma_ec_lib(), ecs, min_lib = 3)
  expect_identical(pw$id, "syn01")
  expect_named(pw, c("cpds", "name", "emp_cpds", "id"))
  expect_setequal(pw$emp_cpds[[1]], ecs$total_matched_ecpds)
  expect_identical(mmc_ma_pathway_ecs(ma_ec_lib(), ecs, min_lib = 2)$id,
                   c("syn01", "syn02"))
})

# ---- parity with the pinned MetaboAnalystR ---------------------------------------

ma_ecs_ref <- function() {
  ref_file <- file.path(ma_fixture_dir, "ecs_reference.rds")
  skip_if_not(file.exists(ref_file), "EC parity reference not generated yet")
  ref <- readRDS(ref_file)
  expect_identical(ref$metaboanalystr_commit,
                   "398476ae2a0c996e390925fa62adc46b5ecde334")
  # EC labels follow locale collation (see 06l): compare under the same one.
  skip_if_not(identical(Sys.getlocale("LC_COLLATE"), ref$locale_collate),
              paste0("EC reference generated under LC_COLLATE '",
                     ref$locale_collate, "'"))
  ref
}

for (case in c("mixed", "positive", "negative", "mixed_no_primary")) {
  test_that(sprintf("EC construction matches MetaboAnalystR (%s)", case), {
    ref <- ma_ecs_ref()
    up <- ref$cases[[case]]
    inp <- ref$inputs
    lib <- mmc_ma_library_from_model(ma_write_model(
      inp$compounds, inp$pathways, pathway_names = inp$pathway_names))
    prep <- mmc_ma_prepare_peaks(inp$peaks, up$mode, rt_frac = inp$rt_frac)
    matched <- mmc_ma_match_peaks(prep, lib, ppm = inp$ppm)
    ecs <- mmc_ma_build_ecs(matched, prep, mmc_ma_input_mzlist(prep, inp$cutoff),
                            force_primary_ion = identical(up$primary_ion, "yes"))

    expect_identical(ecs$initial_ecs, up$initial_ecs)
    expect_identical(ecs$matched_res, up$matched_res)
    for (d in c("ec_exp_dict", "cpd_exp_dict", "ec_exp", "mz2cpd_dict",
                "mz2ec_dict", "ec2mz_dict", "ecpd_cpd_dict", "cpd_ecpd_dict",
                "cpd_ecpd_counts", "cpd_form_dict", "input_ecpdlist",
                "total_matched_ecpds")) {
      expect_identical(ecs[[d]], up[[d]], label = d)
    }

    pw <- mmc_ma_pathway_ecs(lib, ecs, min_lib = ref$min_lib)
    expect_identical(pw$cpds, up$pathways$cpds)
    expect_identical(pw$name, up$pathways$name)
    expect_identical(pw$emp_cpds, up$pathways$emp_cpds)
  })
}
