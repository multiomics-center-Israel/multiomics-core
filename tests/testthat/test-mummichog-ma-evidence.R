# tests/testthat/test-mummichog-ma-evidence.R
#
# Pathway Evidence for the MetaboAnalyst-compatible mummichog engine (06s, and
# the shared assembler in 06f): built from the engine's own result files, the
# model and the dataset's annotations, with expected values worked out by hand
# from the synthetic fixture `ma_evidence_fixture()` (helper-mummichog-ma.R,
# which documents every EC and feature). Also checks that building the
# evidence changes no statistic, writes no engine file, and cannot invalidate
# the engine stage's {targets} dependencies.
#
# None of these tests needs fitdistrplus or the parity reference: the engine is
# not run, only its output files are read. All data SYNTHETIC.

ma_ev <- function(fx = ma_evidence_fixture(), ...) {
  build_mummichog_ma_pathway_evidence(fx$pathways, fx$files, fx$model, fx$annot,
                                      de_table = fx$de, ...)
}

feature_row <- function(ev, pathway, ec, feature) {
  ft <- ev$feature_table
  r <- ft[ft$Pathway == pathway & ft$EmpiricalCompound == ec & ft$Feature == feature, ,
          drop = FALSE]
  expect_identical(nrow(r), 1L, label = paste(pathway, ec, feature))
  r
}

ec_row <- function(ev, pathway, ec) {
  et <- ev$ec_table
  r <- et[et$Pathway == pathway & et$EmpiricalCompound == ec, , drop = FALSE]
  expect_identical(nrow(r), 1L, label = paste(pathway, ec))
  r
}

# ---- support counts -------------------------------------------------------------

test_that("Supporting ECs equal each pathway's Sig.EC, and the ECs are ec_hits", {
  fx <- ma_evidence_fixture()
  ev <- ma_ev(fx)
  ps <- ev$pathway_summary
  expect_setequal(ps$Pathway, fx$pathways$pathway)
  sig <- stats::setNames(fx$pathways$Sig.EC, fx$pathways$pathway)
  expect_equal(ps[["Supporting ECs"]], unname(sig[ps$Pathway]))
  for (i in seq_len(nrow(fx$pathways))) {
    pw <- fx$pathways$pathway[i]
    expect_setequal(ev$ec_table$EmpiricalCompound[ev$ec_table$Pathway == pw],
                    strsplit(fx$pathways$ec_hits[i], ";")[[1]])
  }
})

test_that("distinct supporting features and feature-EC links are counted apart", {
  ev <- ma_ev()
  ps <- ev$pathway_summary
  gly <- ps[ps$Pathway == "Glycolysis", ]
  # F01 sits in EC0001 and EC0003: one feature, two links.
  expect_equal(gly[["Supporting features"]], 5)
  expect_equal(gly[["Feature-EC links"]], 6)
  pyr <- ps[ps$Pathway == "Pyruvate metabolism", ]
  expect_equal(pyr[["Supporting ECs"]], 1)
  expect_equal(pyr[["Supporting features"]], 2)
  expect_equal(pyr[["Feature-EC links"]], 2)
  # The ORA statistics are passed through untouched.
  expect_equal(gly$p.value, 0.01)
  expect_equal(gly$Overlap, 3)
  expect_equal(gly[["Detected pathway size"]], 4)
})

# ---- candidates -----------------------------------------------------------------

test_that("pathway-matching candidates are the EC's compounds in THIS pathway", {
  ev <- ma_ev()
  # EC0001 also holds C99999 (in no pathway): never a candidate.
  e1 <- ec_row(ev, "Glycolysis", "EC0001")
  expect_identical(e1[["Pathway-matching candidate(s)"]], "C00031")
  expect_identical(e1[["Candidate name(s)"]], "D-Glucose")      # from the model
  # EC0002 reaches Glycolysis through C00022 only ...
  expect_identical(ec_row(ev, "Glycolysis", "EC0002")[["Pathway-matching candidate(s)"]],
                   "C00022")
  # ... and Pyruvate metabolism through both of its compounds, counted once.
  e2 <- ec_row(ev, "Pyruvate metabolism", "EC0002")
  expect_identical(e2[["Pathway-matching candidate(s)"]], "C00022; C00186")
  expect_identical(e2[["Candidate name(s)"]], "Pyruvate; (S)-Lactate")
  expect_identical(e2[["Candidate KEGG ID(s)"]], "C00022; C00186")
  expect_identical(e2[["# Features"]], 2L)
})

# ---- agreement (option B: only the feature's own pathway candidates) ------------

test_that("each feature is compared only with the pathway candidates it matched", {
  ev <- ma_ev()
  expect_identical(feature_row(ev, "Glycolysis", "EC0001", "F01")$Agreement, "Match")
  expect_identical(feature_row(ev, "Glycolysis", "EC0001", "F02")$Agreement, "Conflict")
  # F04 matched only C99999, outside the pathway: not assessed, though its
  # annotation differs from the EC's candidate -- and it stays in the EC.
  f04 <- feature_row(ev, "Glycolysis", "EC0001", "F04")
  expect_identical(f04$Agreement, "Not assessed")
  expect_false(f04[["Matched via pathway candidate"]])
  expect_true(is.na(f04[["Pathway candidate(s) matched by feature"]]))
  expect_identical(f04[["Matched compound(s)"]], "C99999")
  # F01 again, in EC0003, where its only pathway match is C00267.
  f01b <- feature_row(ev, "Glycolysis", "EC0003", "F01")
  expect_identical(f01b[["Pathway candidate(s) matched by feature"]], "C00267")
  expect_identical(f01b$Agreement, "Conflict")
  # Unannotated, though it matched a candidate.
  f12 <- feature_row(ev, "Glycolysis", "EC0002", "0012")
  expect_true(f12[["Matched via pathway candidate"]])
  expect_identical(f12$Agreement, "Not assessed")
  expect_identical(feature_row(ev, "Glycolysis", "EC0002", "F06")$Agreement, "Match")
})

test_that("ECs roll their features up into the four states, counted at two grains", {
  ev <- ma_ev()
  e1 <- ec_row(ev, "Glycolysis", "EC0001")
  expect_identical(e1$Agreement, "Mixed")
  expect_identical(c(e1$n_match, e1$n_conflict, e1$n_not_assessed), c(1L, 1L, 1L))
  expect_identical(ec_row(ev, "Glycolysis", "EC0002")$Agreement, "Match")
  expect_identical(ec_row(ev, "Glycolysis", "EC0003")$Agreement, "Conflict")
  expect_identical(ec_row(ev, "Pyruvate metabolism", "EC0002")$Agreement, "Match")

  gly <- ev$pathway_summary[ev$pathway_summary$Pathway == "Glycolysis", ]
  expect_equal(unlist(gly[paste("ECs", c("Match", "Conflict", "Mixed", "Not assessed"))],
                      use.names = FALSE), c(1, 1, 1, 0))
  expect_equal(unlist(gly[paste("feature-EC links", c("Match", "Conflict", "Not assessed"))],
                      use.names = FALSE), c(2, 2, 2))
})

test_that("the python agreement rule is unchanged: all EC candidates, no own-match filter", {
  fx <- ma_evidence_fixture()
  inp <- read_mummichog_ma_evidence_inputs(fx$files, fx$de)
  py <- .mmc_assemble_evidence(
    fx$pathways, lapply(fx$pathways$ec_hits, .mmc_split_cell, sep = ";"),
    split(inp$cand, inp$cand$EID), split(inp$feat, inp$feat$EID),
    fx$model, fx$annot, own_matches = FALSE)
  # Under the python rule F04 is compared with EC0001's candidate and conflicts;
  # the MetaboAnalyst rule leaves it unassessed. The switch is the only change.
  ft <- py$feature_table
  expect_identical(ft$Agreement[ft$EmpiricalCompound == "EC0001" & ft$Feature == "F04"],
                   "Conflict")
  expect_named(ft, c("Pathway", "EmpiricalCompound", "Feature", "m/z", "RT",
                     "Adduct/ion", "Feature p-value",
                     "Feature log2FC (mummichog input)", "Original annotation",
                     "Annotation ID", "Annotation confidence", "Agreement"))
  expect_false("Polarities" %in% names(py$ec_table))
  expect_identical(py$pathway_summary[["Supporting ECs"]],
                   ma_ev(fx)$pathway_summary[["Supporting ECs"]])
})

# ---- polarity, statistic, traceability --------------------------------------------

test_that("mixed polarity: one EC holds both polarities, each feature keeps its own", {
  ev <- ma_ev()
  expect_identical(ec_row(ev, "Glycolysis", "EC0001")$Polarities, "negative; positive")
  expect_identical(ec_row(ev, "Glycolysis", "EC0002")$Polarities, "positive")
  f02 <- feature_row(ev, "Glycolysis", "EC0001", "F02")
  expect_identical(f02$Polarity, "negative")
  expect_identical(f02[["Adduct/ion"]], "M-H[-]")
  expect_identical(feature_row(ev, "Glycolysis", "EC0001", "F01")$Polarity, "positive")
})

test_that("the input statistic and the DE log2FC are separate, correctly named columns", {
  fx <- ma_evidence_fixture()
  ev <- ma_ev(fx)
  ft <- ev$feature_table
  expect_false("Feature log2FC (mummichog input)" %in% names(ft))
  f01 <- feature_row(ev, "Glycolysis", "EC0001", "F01")
  expect_equal(f01[["Feature statistic (MetaboAnalyst input)"]], 4.0)
  expect_equal(f01[["Feature log2FC (DE)"]], 1.5)
  expect_equal(f01[["Feature p-value"]], 0.001)
  expect_equal(f01[["m/z"]], 181.0707)
  # Every row: statistic from ma_input.tsv, log2FC from the DE table.
  inp <- utils::read.delim(file.path(fx$ma_dir, "ma_input.tsv"),
                           colClasses = c(feature_id = "character"))
  expect_equal(ft[["Feature statistic (MetaboAnalyst input)"]],
               inp$statistic[match(ft$Feature, inp$feature_id)])
  expect_equal(ft[["Feature log2FC (DE)"]],
               fx$de$logFC[match(ft$Feature, fx$de$feature_id)])
  # Without a DE table the log2FC column stays, empty, and nothing else moves.
  no_de <- build_mummichog_ma_pathway_evidence(fx$pathways, fx$files, fx$model,
                                               fx$annot)
  expect_true(all(is.na(no_de$feature_table[["Feature log2FC (DE)"]])))
  expect_identical(no_de$feature_table[["Feature statistic (MetaboAnalyst input)"]],
                   ft[["Feature statistic (MetaboAnalyst input)"]])
})

test_that("feature ids stay traceable, including the duplicate m/z case", {
  fx <- ma_evidence_fixture()
  ev <- ma_ev(fx)
  ft <- ev$feature_table
  # Read as text: the leading zeros of "0012" survive.
  expect_true("0012" %in% ft$Feature)
  expect_true(all(ft$Feature %in% c("F01", "F02", "F04", "0012", "F06")))
  # F06 and F07 share an m/z; the engine kept F06 (smaller p). Its row is F06's
  # own -- polarity, p-value, annotation -- and F07 is nowhere.
  expect_false("F07" %in% ft$Feature)
  f06 <- feature_row(ev, "Glycolysis", "EC0002", "F06")
  expect_identical(f06$Polarity, "positive")
  expect_equal(f06[["Feature p-value"]], 0.005)
  expect_identical(f06[["Original annotation"]], "Pyruvate")
  expect_false(any(ft[["Original annotation"]] %in% "(S)-Lactate"))
})

test_that("rows are ordered explicitly within an EC (p-value, then feature id)", {
  ev <- ma_ev()
  ft <- ev$feature_table
  e1 <- ft[ft$Pathway == "Glycolysis" & ft$EmpiricalCompound == "EC0001", ]
  expect_identical(e1$Feature, c("F01", "F02", "F04"))
})

# ---- refusals -------------------------------------------------------------------

rewrite <- function(fx, file, edit) {
  f <- file.path(fx$ma_dir, file)
  d <- utils::read.delim(f, colClasses = "character", check.names = FALSE)
  readr::write_tsv(edit(d), f, na = "NA")
  fx
}

test_that("a different model than the one analysed refuses the evidence", {
  fx <- ma_evidence_fixture()
  attr(fx$model, "sha256") <- strrep("0", 64)
  expect_error(ma_ev(fx), "not the model this contrast was analysed with")
  attr(fx$model, "sha256") <- NULL
  expect_error(ma_ev(fx), "currently configured metabolic model is unknown")
  expect_error(ma_ev(ma_evidence_fixture(record_sha = FALSE)), "names no model")
})

test_that("an untraceable or inconsistent match refuses the evidence", {
  no_id <- rewrite(ma_evidence_fixture(), "ma_matched_features.tsv",
                   function(d) { d$feature_id[3] <- NA; d })
  expect_error(ma_ev(no_id), "carry no feature id")

  unknown <- rewrite(ma_evidence_fixture(), "ma_matched_features.tsv",
                     function(d) { d$feature_id[d$feature_id == "F04"] <- "F99"; d })
  expect_error(ma_ev(unknown), "not in 'ma_input.tsv'")

  dup <- rewrite(ma_evidence_fixture(), "ma_input.tsv",
                 function(d) { d$feature_id[6] <- "F06"; d })
  expect_error(ma_ev(dup), "more than once")

  two_ids <- rewrite(ma_evidence_fixture(), "ma_matched_features.tsv",
                     function(d) { d$Query.Mass[d$feature_id == "F02"] <- "181.0707"; d })
  expect_error(ma_ev(two_ids), "not traceable")

  pol <- rewrite(ma_evidence_fixture(), "ma_matched_features.tsv",
                 function(d) { d$polarity[d$feature_id == "F02"] <- "positive"; d })
  expect_error(ma_ev(pol), "different polarity")

  lost_ec <- rewrite(ma_evidence_fixture(), "ma_matched_features.tsv",
                     function(d) d[d$Empirical.Compound != "EC0003", ])
  expect_error(ma_ev(lost_ec), "are not in 'ma_matched_features.tsv'")
})

test_that("support counts that disagree with Sig.EC refuse the evidence", {
  expect_error(ma_ev(ma_evidence_fixture(sig_ec = c(4, 1))),
               "different number of supporting ECs")
})

# ---- report sections: wiring, and nothing statistical changes --------------------

ma_sections <- function(fx, model = fx$model, de_res = fx$de_res) {
  suppressMessages(build_mummichog_report_sections(
    list(HL_vs_LL = fx$pathways), fx$config, files = fx$files,
    row_data = fx$row_data, de_res = de_res, model = model))
}

test_that("MA sections carry the evidence and fill the Summary's evidence columns", {
  fx <- ma_evidence_fixture()
  s <- ma_sections(fx)[["HL_vs_LL"]]
  expect_identical(s$engine, "metaboanalyst")
  expect_null(s$evidence_note)
  expect_null(s$evidence_log2fc_note)          # every feature has a DE logFC
  expect_identical(s$evidence, ma_ev(fx))
  sm <- s$summary
  gly <- sm[sm$Pathway == "Glycolysis", ]
  expect_identical(gly[["Supporting evidence"]], "3 ECs / 5 features")
  expect_identical(gly[["Annotation agreement"]], "1 Match · 1 Conflict · 1 Mixed")
  # Not nominal in ORA (FET 0.2) or GSEA (p 0.6): no evidence claimed.
  expect_identical(sm[["Supporting evidence"]][sm$Pathway == "Pyruvate metabolism"], "—")
})

test_that("the section resolves and checks the configured model itself", {
  fx <- ma_evidence_fixture()
  s <- ma_sections(fx, model = NULL)[["HL_vs_LL"]]
  expect_false(is.null(s$evidence))
  expect_identical(s$evidence, ma_ev(fx))
})

test_that("evidence changes no ORA, GSEA or Summary statistic", {
  fx <- ma_evidence_fixture()
  s <- ma_sections(fx)[["HL_vs_LL"]]
  # ORA table: exactly what the engine's table gives, in its order.
  expect_identical(s$table, build_mummichog_pathway_table(fx$pathways))
  expect_identical(s$table$p.value, fx$pathways$FET)
  # GSEA: the engine's own result, untouched; python GSEA never ran.
  expect_identical(s$gsea, attr(fx$pathways, "ma_gsea"))
  # Summary: everything except the evidence cells equals a Summary built with
  # no evidence at all.
  base <- build_mummichog_summary(
    s$table, gsea_table = s$gsea$table, evidence = NULL,
    pathway_ids = fx$model$pathway_ids,
    identified = mmc_identified_support(NULL, "HL_vs_LL"),
    library_ids = character(0), p_cutoff = 0.05)
  ev_cols <- c("Supporting evidence", "Annotation agreement", "Supporting ECs",
               "Supporting features")
  keep <- setdiff(names(base), ev_cols)
  expect_identical(s$summary[, keep], base[, keep])
})

test_that("building sections and exports leaves every engine file byte-identical", {
  fx <- ma_evidence_fixture()
  before <- tools::md5sum(fx$files)
  s <- ma_sections(fx)
  out <- mod_mummichog_exports(s, fx$files, fx$dir)
  expect_identical(tools::md5sum(fx$files), before)
  expect_setequal(list.files(fx$ma_dir), basename(fx$files))
  # The three evidence exports are written beside the others.
  ev_files <- file.path(fx$dir, "mummichog_pinned", sprintf(
    "mummichog_evidence_%s_HL_vs_LL.%s",
    rep(c("pathways", "empirical_compounds", "features"), each = 2),
    c("tsv", "csv")))
  expect_true(all(file.exists(ev_files)))
  expect_true(all(ev_files %in% out))
  ft <- utils::read.delim(ev_files[5], colClasses = "character", check.names = FALSE)
  expect_true("0012" %in% ft$Feature)
})

test_that("a refused contrast keeps ORA and GSEA, and says why", {
  fx <- ma_evidence_fixture()
  bad <- fx$model
  attr(bad, "sha256") <- strrep("a", 64)
  expect_warning(secs <- ma_sections(fx, model = bad),
                 "no Pathway Evidence built .*not the model this contrast was analysed with")
  s <- secs[["HL_vs_LL"]]
  expect_null(s$evidence)
  expect_match(s$evidence_note, "not the model this contrast was analysed with")
  expect_identical(s$table, build_mummichog_pathway_table(fx$pathways))
  expect_identical(s$gsea, attr(fx$pathways, "ma_gsea"))
  expect_identical(s$summary[["Supporting evidence"]][s$summary$Pathway == "Glycolysis"],
                   "Not available")
})

test_that("an unmatched DE contrast leaves log2FC empty and says why", {
  fx <- ma_evidence_fixture()
  # Two tables, neither named for this contrast: none can be picked safely.
  other <- list(de_tables = list(A_vs_B = fx$de, C_vs_D = fx$de))
  expect_warning(secs <- ma_sections(fx, de_res = other),
                 "No differential-abundance table could be matched to contrast 'HL_vs_LL'")
  s <- secs[["HL_vs_LL"]]
  expect_false(is.null(s$evidence))                 # evidence is still built
  expect_null(s$evidence_note)
  ft <- s$evidence$feature_table
  expect_true(all(is.na(ft[["Feature log2FC (DE)"]])))
  expect_match(s$evidence_log2fc_note, "Feature log2FC (DE) column is empty", fixed = TRUE)
  # Nothing else in the evidence moves.
  full <- ma_ev(fx)$feature_table
  same <- setdiff(names(ft), "Feature log2FC (DE)")
  expect_identical(ft[, same], full[, same])
})

test_that("a DE table without logFC, or missing some features, is explained", {
  fx <- ma_evidence_fixture()
  no_col <- list(de_tables = list(HL_vs_LL = fx$de[, c("feature_id", "P.Value")]))
  expect_warning(s <- ma_sections(fx, de_res = no_col)[["HL_vs_LL"]],
                 "has no logFC column")
  expect_true(all(is.na(s$evidence$feature_table[["Feature log2FC (DE)"]])))

  partial <- list(de_tables = list(HL_vs_LL = fx$de[fx$de$feature_id != "F04", ]))
  expect_warning(s <- ma_sections(fx, de_res = partial)[["HL_vs_LL"]],
                 "1 of 5 evidence feature\\(s\\) have no log2FC")
  ft <- s$evidence$feature_table
  expect_true(all(is.na(ft[["Feature log2FC (DE)"]][ft$Feature == "F04"])))
  expect_false(anyNA(ft[["Feature log2FC (DE)"]][ft$Feature != "F04"]))
})

test_that("a contrast without result files gets a reason, not a silent gap", {
  fx <- ma_evidence_fixture()
  r <- suppressWarnings(mmc_ma_section_evidence(fx$pathways, NULL, fx$model, fx$annot))
  expect_null(r$evidence)
  expect_match(r$note, "no MetaboAnalyst-compatible result files")
})

test_that("python sections never get an evidence note", {
  pw <- data.frame(pathway = c("a", "b"), overlap_size = c(1, 2),
                   pathway_size = c(4, 5), `p-value` = c(0.2, 0.01),
                   check.names = FALSE)
  s <- suppressMessages(build_mummichog_report_sections(
    list(A_vs_B = pw), list(modes = list(metabolomics = list(enrichment = list(
      mummichog = list(enabled = TRUE, engine = "python")))))))[["A_vs_B"]]
  expect_identical(s$engine, "python")
  expect_null(s$evidence)
  expect_null(s$evidence_note)
})

# ---- {targets}: the engine stage's dependencies are untouched ---------------------

test_that("no evidence or report-section function is reachable from the engine stage", {
  skip_if_not_installed("codetools")
  env <- environment(mod_mummichog_metaboanalyst)
  # Every object of ours the stage's command can reach, as {targets} walks it:
  # the globals of each function, recursively.
  seen <- character(0)
  queue <- c("mod_mummichog_metaboanalyst", "write_mummichog_engine_info")
  while (length(queue) > 0) {
    nm <- queue[[1]]
    queue <- queue[-1]
    if (nm %in% seen || !exists(nm, envir = env, inherits = FALSE)) next
    seen <- c(seen, nm)
    obj <- get(nm, envir = env, inherits = FALSE)
    if (is.function(obj)) queue <- c(queue, codetools::findGlobals(obj))
  }
  expect_true("read_mummichog_model_pathways" %in% seen)   # the walk works
  changed <- c("build_mummichog_pathway_evidence", ".mmc_assemble_evidence",
               "build_mummichog_ma_pathway_evidence",
               "read_mummichog_ma_evidence_inputs", "mmc_ma_section_evidence",
               "mmc_load_model_pathways", "build_mummichog_report_sections")
  expect_identical(intersect(changed, seen), character(0))
})
