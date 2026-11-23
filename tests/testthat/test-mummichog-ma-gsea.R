# tests/testthat/test-mummichog-ma-gsea.R
#
# MetaboAnalyst-compatible GSEA on the engine's own joint-polarity ECs (06o).
# Unit tests on a small synthetic EC set, plus parity against the pinned
# MetaboAnalystR peaks GSEA run end to end
# (fixtures/mummichog_ma_parity/generate_gsea_reference.R). All data synthetic.

skip_if_not_installed("fgsea")
skip_if_not_installed("fastmatch")
skip_if_not_installed("data.table")
skip_if_not_installed("BiocParallel")

GSEA_TOL <- 1e-8

gs_ecs <- function() {
  ec <- sprintf("E%02d", 1:12)
  list(ec_exp = stats::setNames(c(5, 4.2, 3.1, 2.5, 1.8, 1.1, 0.4, -0.3, -1.2,
                                  -2.2, -3.3, -4.6), ec),
       cpd_exp_dict = stats::setNames(as.list(seq_along(ec)), sprintf("C%02d", 1:12)))
}
gs_pathways <- list(
  cpds = list(sprintf("C%02d", 1:5), sprintf("C%02d", 8:12),
              sprintf("C%02d", c(3, 6, 9, 12)), c("C01", "C99")),
  name = c("top", "bottom", "spread", "small"),
  emp_cpds = list(sprintf("E%02d", 1:5), sprintf("E%02d", 8:12),
                  sprintf("E%02d", c(3, 6, 9, 12)), "E01"),
  id = c("p1", "p2", "p3", "p4"))

test_that("GSEA runs on the ECs' scores and the ORA's pathway EC sets", {
  g <- suppressWarnings(mmc_ma_gsea(gs_ecs(), gs_pathways, n_perm = 100))
  expect_setequal(g$table$Pathway, gs_pathways$name)
  expect_identical(g$table$Pathway, rownames(g$resmat))
  expect_identical(colnames(g$resmat), c("Pathway_Total", "Hits", "P_val", "P_adj", "NES"))
  expect_false(is.unsorted(g$resmat[, "P_val"]))
  # upstream's sizes: model compounds, and those among the matched compounds
  expect_identical(unname(g$resmat["small", "Pathway_Total"]), 2)
  expect_identical(unname(g$resmat["small", "Hits"]), 1)
  # BH over the tested pathways, then rounded as upstream rounds
  f <- g$fgsea
  expect_equal(f$padj, stats::p.adjust(f$pval, "fdr"))
  expect_equal(unname(g$resmat[f$pathway, "P_adj"]), signif(f$padj, 4))
  # each row carries its own pathway's ECs
  expect_identical(g$table$EmpiricalCompounds[g$table$Pathway == "spread"],
                   "E03;E06;E09;E12")
  expect_length(g$deviations, 0L)
})

test_that("non-default GSEA settings are reported as deviations", {
  g <- suppressWarnings(mmc_ma_gsea(gs_ecs(), gs_pathways, n_perm = 50, seed = 7))
  expect_true(any(grepl("^n_perm = 50", g$deviations)))
  expect_true(any(grepl("^seed = 7", g$deviations)))
})

test_that("GSEA does not change the ORA, and travels through the result files", {
  ref_file <- file.path(ma_fixture_dir, "stats_reference.rds")
  skip_if_not(file.exists(ref_file), "stats reference not generated yet")
  skip_if_not_installed("fitdistrplus")
  ref <- readRDS(ref_file)
  inp <- ref$inputs
  lib <- mmc_ma_library_from_model(ma_write_model(inp$compounds, inp$pathways,
                                                  pathway_names = inp$pathway_names))
  args <- list(inp$peaks, lib, "mixed", ppm = inp$ppm, cutoff = inp$cutoff,
               n_perm = 20)
  without <- suppressWarnings(do.call(mmc_ma_run_contrast, args))
  with <- suppressWarnings(do.call(mmc_ma_run_contrast,
                                   c(args, list(gsea = list(metric_label = "x")))))
  expect_null(without$gsea)
  expect_false(is.null(with$gsea))
  expect_identical(with$stats, without$stats)
  expect_identical(with$perms, without$perms)

  dir <- tempfile("ma_gsea_io_")
  files <- write_mummichog_ma_results(with, inp$peaks, dir, info = list(),
                                      pathway_ids = NULL)
  expect_true("ma_gsea_results.tsv" %in% basename(files))
  tbl <- read_mummichog_ma_pathways(files)
  g <- attr(tbl, "ma_gsea")
  expect_identical(g$table$Pathway, with$gsea$table$Pathway)
  expect_equal(g$table$P.Value, with$gsea$table$P.Value)
  expect_equal(g$table$padj, with$gsea$table$padj)
  expect_identical(g$engine, "metaboanalyst")
})

# ---- parity with the pinned MetaboAnalystR peaks GSEA --------------------------

ma_gsea_ref <- function() {
  ref_file <- file.path(ma_fixture_dir, "gsea_reference.rds")
  skip_if_not(file.exists(ref_file), "GSEA parity reference not generated yet")
  ref <- readRDS(ref_file)
  expect_identical(ref$metaboanalystr_commit,
                   "398476ae2a0c996e390925fa62adc46b5ecde334")
  skip_if_not(identical(Sys.getlocale("LC_COLLATE"), ref$locale_collate),
              "GSEA reference generated under another LC_COLLATE")
  skip_if_not(identical(as.character(utils::packageVersion("fgsea")), ref$fgsea_version),
              paste("GSEA reference generated under fgsea", ref$fgsea_version))
  ref
}

for (case in c("mixed", "positive", "negative", "mixed_no_primary")) {
  test_that(sprintf("GSEA matches MetaboAnalystR end to end (%s)", case), {
    ref <- ma_gsea_ref()
    up <- ref$cases[[case]]
    inp <- ref$inputs
    lib <- mmc_ma_library_from_model(ma_write_model(
      inp$compounds, inp$pathways, pathway_names = inp$pathway_names))
    prep <- mmc_ma_prepare_peaks(inp$peaks, up$mode, rt_frac = inp$rt_frac)
    sig <- mmc_ma_input_mzlist(prep, inp$cutoff)
    ecs <- mmc_ma_build_ecs(mmc_ma_match_peaks(prep, lib, ppm = inp$ppm), prep, sig,
                            force_primary_ion = identical(up$primary_ion, "yes"))
    pw <- mmc_ma_pathway_ecs(lib, ecs, min_lib = ref$min_lib)

    # the inputs: EC scores and pathway EC sets
    expect_identical(ecs$ec_exp, up$ec_exp)
    expect_identical(pw$emp_cpds, up$emp_cpds)
    expect_identical(pw$name, up$pathway_names)

    g <- suppressWarnings(mmc_ma_gsea(ecs, pw, n_perm = ref$perm_num))
    f <- g$fgsea
    e <- up$fgsea
    # tested universe, order and sizes: exact
    expect_identical(f$pathway, e$pathway)
    expect_identical(as.numeric(f$size), e$size)
    lead <- ifelse(is.na(f$leadingEdge), "", f$leadingEdge)
    expect_identical(lead, e$leadingEdgeMatched)
    # statistics: tight tolerance (same primitives, same seeds)
    for (col in c("ES", "NES", "pval", "padj", "nMoreExtreme")) {
      expect_equal(as.numeric(f[[col]]), e[[col]], tolerance = GSEA_TOL, label = col)
    }
    # MetaboAnalyst's result table: values, rounding and row order
    expect_identical(rownames(g$resmat), rownames(up$resmat))
    expect_identical(colnames(g$resmat), colnames(up$resmat))
    expect_equal(unname(g$resmat), unname(up$resmat), tolerance = GSEA_TOL)
    # each row's ECs are that pathway's own (upstream's pathwaysFiltered)
    own <- vapply(rownames(g$resmat), function(p)
      paste(names(up$pathways_filtered[[p]]), collapse = ";"), character(1),
      USE.NAMES = FALSE)
    expect_identical(g$table$EmpiricalCompounds, own)
  })
}

test_that("the GSEA fixture is informative", {
  ref <- ma_gsea_ref()
  m <- ref$cases$mixed$fgsea
  expect_gte(nrow(m), 3L)
  expect_true(any(m$pval < 1))
  expect_false(isTRUE(all.equal(stats::p.adjust(m$pval, "fdr"),
                                stats::p.adjust(m$pval, "holm"))))
})
