# tests/testthat/test-mummichog-presentation.R
#
# The reader-facing presentation layer (06t): the overview, the detailed ORA
# table and the input counts, and that adding them leaves the Shiny contract
# (the Summary, the ORA table, the GSEA table) exactly as it was. Expected
# values are worked out by hand from ma_two_contrast_sections() and
# ma_evidence_fixture() (helper-mummichog-ma.R). All data SYNTHETIC.

ov_row <- function(ov, pathway) {
  r <- ov[ov$Pathway == pathway, , drop = FALSE]
  expect_identical(nrow(r), 1L, label = pathway)
  r
}

# ---- the overview ----------------------------------------------------------------

test_that("the overview lists every pathway, and its first rows are the Summary's", {
  tc <- ma_two_contrast_sections()
  s  <- tc$sections[["HL_vs_LL"]]
  expect_setequal(s$overview$Pathway, c("Glycolysis", "Pyruvate metabolism", "TCA cycle"))
  expect_identical(utils::head(s$overview$Pathway, nrow(s$summary)), s$summary$Pathway)
  expect_named(s$overview, .MMC_OVERVIEW_COLS)
  # Ranked by the smaller raw p of ORA and GSEA: TCA (GSEA 0.005) first.
  expect_identical(s$overview$Pathway,
                   c("TCA cycle", "Glycolysis", "Pyruvate metabolism"))
})

test_that("raw and adjusted p-values of both tests are kept apart", {
  ov <- ma_two_contrast_sections()$sections[["HL_vs_LL"]]$overview
  g <- ov_row(ov, "Glycolysis")
  expect_equal(g[["ORA raw p"]], 0.01)
  expect_equal(g[["ORA Holm-adjusted p"]], 0.02)
  expect_identical(g[["ORA status"]], "Holm-adjusted p < 0.05")
  expect_equal(g[["GSEA raw p"]], 0.03)
  expect_equal(g[["GSEA BH FDR"]], 0.06)
  expect_identical(g[["GSEA status"]], "Nominal (raw p < 0.05)")
  expect_identical(g[["ORA raw p type"]], "FET")
  p <- ov_row(ov, "Pyruvate metabolism")
  expect_identical(p[["ORA status"]], "Not nominal")
  expect_identical(p[["GSEA status"]], "Not nominal")
})

test_that("a pathway with no significant EC is never called 'Not tested'", {
  ov <- ma_two_contrast_sections()$sections[["HL_vs_LL"]]$overview
  tca <- ov_row(ov, "TCA cycle")
  expect_identical(tca[["ORA status"]], "No significant EC overlap (not in ORA table)")
  expect_true(is.na(tca[["ORA raw p"]]))
  expect_identical(tca[["GSEA status"]], "BH FDR < 0.05")
  expect_false(any(grepl("Not tested", unlist(ov[c("ORA status")]), fixed = TRUE)))
})

test_that("ORA evidence is shown only for ORA-nominal pathways, never for GSEA-only ones", {
  tc <- ma_two_contrast_sections()
  a <- tc$sections[["HL_vs_LL"]]$overview
  g <- ov_row(a, "Glycolysis")
  expect_equal(g[["ORA overlap ECs (Sig.EC)"]], 3)            # the engine's Sig.EC
  expect_equal(g[["Features in overlap ECs"]], 5)
  expect_identical(g[["ORA evidence"]], "See Pathway Evidence")
  expect_identical(g[["Annotation agreement (ECs)"]], "1 Match · 1 Conflict · 1 Mixed")
  tca <- ov_row(a, "TCA cycle")
  expect_true(is.na(tca[["ORA overlap ECs (Sig.EC)"]]))
  expect_identical(tca[["ORA evidence"]], "GSEA-only signal: no ORA evidence shown")

  # LL_vs_DD: Glycolysis has an ORA table row (Sig.EC 3, FET 0.3) and evidence,
  # but its signal is GSEA's -- the ORA counts must not appear as its support.
  b <- ov_row(tc$sections[["LL_vs_DD"]]$overview, "Glycolysis")
  expect_identical(b[["GSEA status"]], "BH FDR < 0.05")
  expect_identical(b[["ORA status"]], "Not nominal")
  expect_true(is.na(b[["ORA overlap ECs (Sig.EC)"]]))
  expect_true(is.na(b[["Features in overlap ECs"]]))
  expect_identical(b[["Annotation agreement (ECs)"]], "—")
  expect_identical(b[["ORA evidence"]], "GSEA-only signal: no ORA evidence shown")
})

test_that("identified-metabolite wording is cautious and complete", {
  fx <- ma_evidence_fixture()
  tbl <- build_mummichog_pathway_table(fx$pathways)
  ident <- data.frame(pathway_id = c("syn00010", "syn00620"),
                      method = c("ORA", "QEA"), supported = c(TRUE, FALSE),
                      stringsAsFactors = FALSE)
  ov <- build_mummichog_overview(tbl, pathway_ids = fx$model$pathway_ids,
                                 identified = ident,
                                 library_ids = c("syn00010", "syn00620"),
                                 engine = "metaboanalyst")
  g <- ov_row(ov, "Glycolysis")
  expect_identical(g[["Identified-metabolite enrichment"]], "ORA (FDR < 0.05)")
  expect_identical(g$Interpretation, "Also supported by identified-metabolite enrichment")
  p <- ov_row(ov, "Pyruvate metabolism")
  expect_identical(p[["Identified-metabolite enrichment"]], "Tested; none at FDR < 0.05")
  expect_identical(p$Interpretation, "No nominal Mummichog signal")
  expect_false(any(grepl("Confirms", unlist(ov), fixed = TRUE)))
  # Every label the Summary can give has reader-facing wording.
  expect_setequal(names(.MMC_INTERPRETATION_LABELS),
                  c("No nominal Mummichog signal", "Not comparable",
                    "Confirms identified signal", "Mummichog-only candidate"))
  expect_error(.mmc_relabel("Something new", .MMC_INTERPRETATION_LABELS),
               "no reader-facing wording")
})

test_that("the python engine keeps its own ORA semantics in the overview", {
  pw <- data.frame(pathway = c("a", "b"), overlap_size = c(2, 1),
                   pathway_size = c(5, 4), `p-value` = c(0.01, 0.3),
                   check.names = FALSE)
  tbl <- build_mummichog_pathway_table(pw)
  gsea <- data.frame(Pathway = c("a", "c"), P.Value = c(0.2, 0.01),
                     padj = c(0.3, 0.02), stringsAsFactors = FALSE)
  ov <- build_mummichog_overview(tbl, gsea_table = gsea, engine = "python")
  a <- ov_row(ov, "a")
  expect_true(is.na(a[["ORA Holm-adjusted p"]]))
  expect_identical(a[["ORA status"]], "Nominal (raw p < 0.05)")   # never "Holm"
  expect_identical(a[["ORA raw p type"]], "empirical")
  expect_identical(ov_row(ov, "c")[["ORA status"]], "Not in ORA results")
  expect_identical(ov_row(ov, "b")[["GSEA status"]], "Not tested by GSEA")
})

test_that("the evidence level rule is the one shared rule", {
  expect_identical(.mmc_evidence_level(c(0.01, 0.01, 0.2, NA), c(0.02, 0.2, 0.01, NA),
                                       p_cutoff = 0.05),
                   c("adjusted", "nominal", "none", "absent"))
  # An adjusted value that is NA can never give "adjusted".
  expect_identical(.mmc_evidence_level(0.001, NA, 0.05), "nominal")
  expect_identical(.mmc_evidence_level(0.06, NA, 0.05), "none")
  # With a stricter p_cutoff, "adjusted" still requires the raw p to be nominal,
  # so the level never contradicts the Summary's own nominal rule.
  expect_identical(.mmc_evidence_level(0.02, 0.04, p_cutoff = 0.01), "none")
})

test_that("with p_cutoff below 0.05 the overview never contradicts the Summary", {
  pw <- data.frame(pathway = c("a", "b"), overlap_size = c(2, 1),
                   pathway_size = c(5, 4), `p-value` = c(0.02, 0.005),
                   check.names = FALSE)
  tbl <- build_mummichog_pathway_table(pw)
  tbl[["adj.p Holm FET"]] <- c(0.04, 0.01)
  ov <- build_mummichog_overview(tbl, p_cutoff = 0.01, engine = "metaboanalyst")
  a <- ov_row(ov, "a")      # FET 0.02, Holm 0.04: not nominal at 0.01
  expect_identical(a[["ORA status"]], "Not nominal")
  expect_identical(a$Interpretation, "No nominal Mummichog signal")
  expect_true(is.na(a[["ORA overlap ECs (Sig.EC)"]]))
  b <- ov_row(ov, "b")      # FET 0.005, Holm 0.01: nominal and adjusted
  expect_identical(b[["ORA status"]], "Holm-adjusted p < 0.05")
  expect_false(identical(b$Interpretation, "No nominal Mummichog signal"))
})

test_that("p-value and count columns are named for the report's table helper", {
  tc <- ma_two_contrast_sections()
  s  <- tc$sections[["HL_vs_LL"]]
  f  <- mmc_dt_column_formats(s$ora_detail)
  expect_true(all(c("ORA raw p (FET)", "ORA raw p (EASE)", "ORA p (Gamma-adjusted FET)",
                    "ORA p (empirical)", "Holm-adjusted FET", "Holm-adjusted EASE",
                    "Holm-adjusted Gamma") %in% f$pval))
  expect_true(all(c("ORA overlap ECs (Sig.EC)", "ECs in pathway (Total.EC)",
                    "Putative compound matches, overlap ECs (Hits.sig)") %in% f$int))
  expect_false("Expected overlap ECs" %in% c(f$pval, f$int))   # not a whole number
  fo <- mmc_dt_column_formats(s$overview)
  expect_setequal(fo$pval, c("ORA raw p", "ORA Holm-adjusted p", "GSEA raw p", "GSEA BH FDR"))
  expect_true(all(c("ORA overlap ECs (Sig.EC)", "Features in overlap ECs") %in% fo$int))
  fa <- mmc_dt_column_formats(build_mummichog_across_contrasts(tc$sections)$wide)
  expect_true(all(c("Contrasts with adjusted evidence", "HL_vs_LL: GSEA BH FDR") %in%
                    c(fa$int, fa$pval)))
  expect_true("HL_vs_LL: GSEA BH FDR" %in% fa$pval)
})

# ---- detailed ORA table ----------------------------------------------------------

test_that("the detailed ORA table copies every engine column, Hits.* named as putative", {
  fx <- ma_evidence_fixture()
  s  <- suppressMessages(build_mummichog_report_sections(
    list(HL_vs_LL = fx$pathways), fx$config, files = fx$files,
    row_data = fx$row_data, de_res = fx$de_res, model = fx$model))[["HL_vs_LL"]]
  od <- s$ora_detail
  pw <- fx$pathways[order(fx$pathways$ma_rank), ]
  expect_identical(od$Pathway, pw$pathway)
  expect_equal(od[["Putative compound matches, overlap ECs (Hits.sig)"]], pw$Hits.sig)
  expect_equal(od[["Putative compound matches, all ECs (Hits.total)"]], pw$Hits.total)
  expect_equal(od[["ORA overlap ECs (Sig.EC)"]], pw$Sig.EC)
  expect_equal(od[["ORA raw p (FET)"]], pw$FET)
  expect_equal(od[["Holm-adjusted FET"]], pw$Holm_adj_FET)
  expect_identical(od[["Overlap EC IDs"]], pw$ec_hits)
})

# ---- input counts ----------------------------------------------------------------

test_that("input counts are the run record's, and pass_any_contrast is only reported", {
  fx <- ma_evidence_fixture()
  info <- attr(fx$pathways, "ma_info")
  info$counts <- list(features_given = 5, features_kept = 4,
                      removed_out_of_range = 0, removed_duplicate_mz = 1,
                      significant_features = 3, empirical_compounds = 3,
                      significant_ecs = 3, pathways_tested = 3, pathways_reported = 2)
  attr(fx$pathways, "ma_info") <- info
  de_res <- list(de_tables = list(HL_vs_LL = fx$de),
                 summary_df = data.frame(feature_id = fx$de$feature_id,
                                         pass_any_contrast = c(1, 0, 0, 1, 1, 0)))
  ic <- mmc_input_counts(fx$pathways, fx$de, de_res)
  expect_equal(ic[["Features at DE (this contrast)"]], 6)
  expect_equal(ic[["Passed to the engine"]], 5)
  expect_equal(ic[["Not passed (missing m/z, RT, p or statistic)"]], 1)
  expect_equal(ic[["Removed: duplicate m/z"]], 1)
  expect_equal(ic[["Features retained for Mummichog analysis"]], 4)
  expect_equal(ic[["ORA significant features (raw DE p ≤ cutoff)"]], 3)
  expect_identical(ic[["Significance rule"]], "raw DE p ≤ 0.05")
  expect_false(any(c("Eligible input features", "Significant input features") %in% names(ic)))
  # The rule is the cutoff the run recorded, whatever the config says now.
  info$params$p_cutoff <- 0.01
  attr(fx$pathways, "ma_info") <- info
  expect_identical(mmc_input_counts(fx$pathways, fx$de, de_res)[["Significance rule"]],
                   "raw DE p ≤ 0.01")
  expect_equal(ic[["Features in DE results (all contrasts)"]], 6)
  expect_equal(ic[["Passing overall DE selection (pass_any_contrast; not the Mummichog input)"]], 3)
  # No run record (python engine): engine counts are NA, nothing is invented.
  py <- mmc_input_counts(data.frame(pathway = "a"), fx$de, NULL)
  expect_true(is.na(py[["Features retained for Mummichog analysis"]]))
  expect_true(is.na(py[["Significance rule"]]))
  expect_equal(py[["Features at DE (this contrast)"]], 6)
})

# ---- the Shiny contract is untouched -----------------------------------------------

test_that("the Shiny entry, Summary and ORA table are unchanged by the presentation layer", {
  fx <- ma_evidence_fixture()
  s <- suppressMessages(build_mummichog_report_sections(
    list(HL_vs_LL = fx$pathways), fx$config, files = fx$files,
    row_data = fx$row_data, de_res = fx$de_res, model = fx$model))[["HL_vs_LL"]]
  # The Summary is still build_mummichog_summary() with the stage's arguments.
  base <- build_mummichog_summary(
    s$table, gsea_table = s$gsea$table, evidence = s$evidence,
    pathway_ids = fx$model$pathway_ids,
    identified = mmc_identified_support(NULL, "HL_vs_LL"),
    library_ids = character(0), p_cutoff = 0.05)
  expect_identical(s$summary, base)
  expect_identical(s$table, build_mummichog_pathway_table(fx$pathways))
  entry <- metab_shiny_mummichog_entry(s)
  expect_named(entry, c(.METAB_SHINY_MUMMICHOG_KEYS, .METAB_SHINY_MUMMICHOG_ADDED_KEYS))
  expect_identical(entry$summary_table, flatten_list_columns(base))
  expect_identical(entry$ora_table, flatten_list_columns(s$table))
  expect_false(any(c("overview", "ora_detail", "input_counts") %in% names(entry)))
})
