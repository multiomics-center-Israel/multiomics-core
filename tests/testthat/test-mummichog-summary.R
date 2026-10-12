# tests/testthat/test-mummichog-summary.R
#
# Unit tests for the mummichog Summary layer (06h): pathway selection, the
# exact-id join to identified-metabolite enrichment, the joinable / tested /
# supported labels, the contribution category and the compact cells.
#
# All fixtures are SYNTHETIC: invented pathway ids, names and p-values.

# Section-shaped inputs for one contrast.
#   P1  ORA nominal, GSEA nominal  -> evidence present
#   P2  ORA not nominal, GSEA nominal (GSEA-only), yet with an ORA evidence row
#   P3  ORA nominal only           -> evidence present
#   P4  nominal nowhere            -> ranked last, "No nominal Mummichog signal"
#   P5  GSEA nominal, no model id  -> not comparable
summary_fixture <- function() {
  ora <- data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    "Pathway" = c("Path one", "Path two", "Path three", "Path four"),
    "p.value" = c(0.01, 0.40, 0.03, 0.60))
  gsea <- data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    "Pathway" = c("Path one", "Path two", "Path four", "Path five"),
    "P.Value" = c(0.02, 0.004, 0.70, 0.045),
    "padj"    = c(0.30, 0.03, 0.80, 0.30))
  evidence <- list(pathway_summary = data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    "Pathway"             = c("Path one", "Path three", "Path two"),
    "Supporting ECs"      = c(3, 1, 2),
    "Supporting features" = c(5, 1, 2),
    "ECs Match"           = c(1, 0, 2),
    "ECs Conflict"        = c(2, 0, 0),
    "ECs Mixed"           = c(0, 0, 0),
    "ECs Not assessed"    = c(0, 1, 0)))
  ids <- c("Path one" = "xx00010", "Path two" = "xx00020",
           "Path three" = "xx00030", "Path four" = "xx00040",
           "Path five" = NA_character_)
  list(ora = ora, gsea = gsea, evidence = evidence, ids = ids)
}

# Identified-enrichment results shaped like mod_metabolomics_enrichment().
identified_fixture <- function() {
  tbl <- function(pathway, fdr) data.frame(pathway = pathway, FDR = fdr,
                                           stringsAsFactors = FALSE)
  list(
    ora_by_contrast = list(A_vs_B = list(res = list(table = tbl(
      c("xx00010 - Path one", "xx00030 - Path three"), c(0.01, 0.50))))),
    qea_by_contrast = list(A_vs_B = list(res = list(table = tbl(
      c("xx00010 - Path one", "xx00030 - Path three"), c(0.04, 0.20))))),
    gsea_by_contrast = list(A_vs_B = list(res = list(table = tbl(
      "xx00030 - Path three", 0.90)))),
    # ssGSEA is deliberately not read, even when it would say "significant".
    ssgsea_by_contrast = list(A_vs_B = list(res = list(table = tbl(
      "xx00020 - Path two", 0.001))))
  )
}

build_fixture_summary <- function(library_ids = c("xx00010", "xx00020",
                                                  "xx00030", "xx00040")) {
  f <- summary_fixture()
  build_mummichog_summary(
    f$ora, gsea_table = f$gsea, evidence = f$evidence, pathway_ids = f$ids,
    identified  = mmc_identified_support(identified_fixture(), "A_vs_B"),
    library_ids = library_ids, p_cutoff = 0.05)
}

row_of <- function(s, pw) s[s$Pathway == pw, , drop = FALSE]


test_that("the pathway id is the label text before the first ' - '", {
  expect_identical(
    .mmc_identified_pathway_id(c("mng00010 - Glycolysis / Gluconeogenesis",
                                 "mng00020 - A - B", "mng00030")),
    c("mng00010", "mng00020", "mng00030"))
})

test_that("identified support reads ORA, QEA and GSEA only, at FDR < 0.05", {
  sup <- mmc_identified_support(identified_fixture(), "A_vs_B")
  expect_setequal(unique(sup$method), c("ORA", "QEA", "GSEA"))
  expect_false("xx00020" %in% sup$pathway_id)
  expect_identical(sup$supported[sup$method == "QEA" & sup$pathway_id == "xx00010"],
                   TRUE)
  expect_equal(nrow(mmc_identified_support(NULL, "A_vs_B")), 0)
  expect_equal(nrow(mmc_identified_support(identified_fixture(), "other")), 0)
})

test_that("every pathway is ranked by min(ORA p, GSEA p), missing as Inf", {
  s <- build_fixture_summary()
  expect_identical(s$Pathway,
                   c("Path two", "Path one", "Path three", "Path five", "Path four"))
  expect_identical(s$`Nominal Mummichog signal`, c(TRUE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(names(s)[seq_along(.MMC_SUMMARY_DISPLAY_COLS)],
                   .MMC_SUMMARY_DISPLAY_COLS)
})

test_that("only the top_n ranked pathways are kept", {
  f <- summary_fixture()
  s <- build_mummichog_summary(f$ora, f$gsea, f$evidence, f$ids,
                               library_ids = "xx00010", top_n = 2)
  expect_identical(s$Pathway, c("Path two", "Path one"))
  # Pathways missing from both tables' p-values still rank, last.
  f$ora$p.value[f$ora$Pathway == "Path four"] <- NA
  f$gsea$P.Value[f$gsea$Pathway == "Path four"] <- NA
  s <- build_mummichog_summary(f$ora, f$gsea, f$evidence, f$ids,
                               library_ids = "xx00010")
  expect_identical(tail(s$Pathway, 1), "Path four")
})

test_that("joinable, tested and supported are kept apart", {
  s <- build_fixture_summary()

  # supported by identified ORA and QEA (fixed display order)
  expect_identical(row_of(s, "Path one")$`Identified-metabolite support`, "ORA + QEA")
  expect_identical(row_of(s, "Path one")$`Mummichog contribution`,
                   "Confirms identified signal")

  # tested by all three methods, significant in none
  expect_identical(row_of(s, "Path three")$`Identified-metabolite support`, "None")
  expect_identical(row_of(s, "Path three")$`Mummichog contribution`,
                   "Mummichog-only candidate")

  # in the shared library but absent from every identified table (ssGSEA ignored)
  expect_identical(row_of(s, "Path two")$`Identified-metabolite support`, "Not tested")
  expect_identical(row_of(s, "Path two")$`Mummichog contribution`,
                   "Mummichog-only candidate")

  # no model id at all
  expect_identical(row_of(s, "Path five")$`Identified-metabolite support`,
                   "Not comparable")
  expect_identical(row_of(s, "Path five")$`Mummichog contribution`, "Not comparable")

  # ranked but not nominal: never a candidate; identified support is unchanged
  four <- row_of(s, "Path four")
  expect_identical(four$`Mummichog contribution`, "No nominal Mummichog signal")
  expect_identical(four$`Identified-metabolite support`, "Not tested")
  expect_identical(four$`Supporting evidence`, "—")
})

test_that("an id outside the identified pathway library is not comparable", {
  s <- build_fixture_summary(library_ids = c("xx00020", "xx00030"))
  expect_identical(row_of(s, "Path one")$`Identified-metabolite support`,
                   "Not comparable")
  expect_identical(row_of(s, "Path one")$`Mummichog contribution`, "Not comparable")
})

test_that("the compact p-value, evidence and agreement cells", {
  s <- build_fixture_summary()
  one <- row_of(s, "Path one")
  expect_identical(one$`Mummichog ORA`, "p = 0.01")
  expect_identical(one$`Mummichog GSEA`, "p = 0.02; FDR = 0.3")
  expect_identical(one$`Supporting evidence`, "3 ECs / 5 features")
  expect_identical(one$`Annotation agreement`, "1 Match · 2 Conflict")
  expect_identical(row_of(s, "Path three")$`Supporting evidence`, "1 EC / 1 feature")
  expect_identical(row_of(s, "Path three")$`Annotation agreement`, "1 Not assessed")

  # GSEA-only: its ORA evidence row is not shown as support for the signal.
  two <- row_of(s, "Path two")
  expect_identical(two$`Supporting evidence`, "Not available — GSEA-only signal")
  expect_identical(two$`Annotation agreement`, "—")
  expect_true(is.na(two$`Supporting ECs`))
  expect_true(is.na(two$`Supporting features`))
  expect_identical(row_of(s, "Path five")$`Mummichog ORA`, "Not tested")
  expect_identical(row_of(s, "Path three")$`Mummichog GSEA`, "Not tested")

  # the values behind the labels travel with the export
  expect_equal(one$`GSEA padj`, 0.30)
  expect_identical(one$`Pathway ID`, "xx00010")
  expect_identical(one$`Identified methods tested`, "ORA + QEA")
})

test_that("missing GSEA reads 'Not run' and the ORA side still works", {
  f <- summary_fixture()
  s <- build_mummichog_summary(f$ora, gsea_table = NULL, evidence = f$evidence,
                               pathway_ids = f$ids, library_ids = "xx00010")
  expect_identical(s$Pathway, c("Path one", "Path three", "Path two", "Path four"))
  expect_true(all(s$`Mummichog GSEA` == "Not run"))
  # no identified results supplied: comparable pathways are untested
  expect_identical(s$`Identified-metabolite support`,
                   c("Not tested", "Not comparable", "Not comparable", "Not comparable"))
  expect_identical(s$`Mummichog contribution`,
                   c("Mummichog-only candidate", "Not comparable",
                     "No nominal Mummichog signal", "No nominal Mummichog signal"))
})

test_that("with no nominal pathway every row says so and none is a candidate", {
  f <- summary_fixture()
  f$ora$p.value <- 0.5
  f$gsea$P.Value <- 0.5
  s <- build_mummichog_summary(f$ora, f$gsea, f$evidence, f$ids,
                               library_ids = "xx00010")
  expect_equal(nrow(s), 5)
  expect_false(any(s$`Nominal Mummichog signal`))
  expect_true(all(s$`Mummichog contribution` == "No nominal Mummichog signal"))
  expect_true(all(s$`Supporting evidence` == "—"))
})

test_that("no pathway results at all gives a zero-row summary with the full schema", {
  s <- build_mummichog_summary(NULL, NULL)
  expect_equal(nrow(s), 0)
  expect_true(all(.MMC_SUMMARY_DISPLAY_COLS %in% names(s)))
})

test_that("agreement labels omit zero categories", {
  ps <- data.frame(check.names = FALSE,
                   "ECs Match" = c(0, 0), "ECs Conflict" = c(2, 0),
                   "ECs Mixed" = c(1, 0), "ECs Not assessed" = c(0, 0))
  expect_identical(mmc_format_agreement(ps), c("2 Conflict · 1 Mixed", NA))
})

test_that("plot-tab sentences separate nominal from BH-adjusted significance", {
  expect_identical(mmc_nominal_sentence(c(0.2, 0.3), 0.05, "Mummichog ORA"),
                   "No pathway reached nominal p < 0.05 in Mummichog ORA.")
  expect_identical(mmc_nominal_sentence(c(0.01, 0.3), 0.05, "Mummichog ORA"),
                   "One pathway reached nominal p < 0.05 in Mummichog ORA.")
  expect_identical(
    mmc_nominal_sentence(c(0.01, 0.02, 0.03, 0.04, 0.5), 0.05, "Mummichog GSEA",
                         padj = c(0.6, 0.6, 0.6, 0.6, 0.9)),
    "4 pathways reached nominal p < 0.05 in Mummichog GSEA; none remained significant after BH correction (FDR < 0.05).")
  expect_identical(
    mmc_nominal_sentence(c(0.001, 0.02), 0.05, "Mummichog GSEA",
                         padj = c(0.002, 0.2)),
    "2 pathways reached nominal p < 0.05 in Mummichog GSEA; 1 remained significant after BH correction (FDR < 0.05).")
})

test_that("the identified pathway library is read from the configured GMT", {
  dir <- withr::local_tempdir()
  gmt <- file.path(dir, "lib.gmt")
  writeLines(c("xx00010\tPath one\tC00001\tC00002",
               "xx00020\tPath two\tC00003\tC00004"), gmt)
  cfg <- list(modes = list(metabolomics = list(enrichment = list(gmt_file = gmt))))
  expect_setequal(mmc_identified_pathway_library(cfg), c("xx00010", "xx00020"))
  expect_length(mmc_identified_pathway_library(list()), 0)
})

test_that("save_mummichog_exports writes the summary and clears a stale one", {
  out <- withr::local_tempdir()
  s <- build_fixture_summary()
  written <- save_mummichog_exports(NULL, NULL, out, contrast_label = "A_vs_B",
                                    summary = s)
  expect_setequal(basename(written),
                  c("mummichog_summary_A_vs_B.tsv", "mummichog_summary_A_vs_B.csv"))
  back <- readr::read_tsv(written[grepl("tsv$", written)], show_col_types = FALSE)
  expect_equal(nrow(back), nrow(s))

  # A rerun without a summary must not leave the old one reading as current.
  save_mummichog_exports(NULL, NULL, out, contrast_label = "A_vs_B")
  expect_false(any(file.exists(written)))
})
