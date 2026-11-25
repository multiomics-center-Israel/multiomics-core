# tests/testthat/test-mummichog-ma-gsea.R
#
# GSEA on the MetaboAnalyst-compatible engine's own joint-polarity ECs (06o).
# The default method is standard fgsea (mmc_gsea_fgsea(), shared with the
# python engine). The LEGACY method (gsea_method = "metaboanalyst_legacy") is
# tested for historical compatibility: unit tests on a small synthetic EC set,
# plus parity against the pinned MetaboAnalystR peaks GSEA run end to end
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

test_that("the default GSEA is standard fgsea on the ECs' scores and the ORA's pathway EC sets", {
  ecs <- gs_ecs()
  g <- suppressWarnings(mmc_ma_gsea(ecs, gs_pathways))
  expect_identical(g$method, "fgsea")
  expect_match(g$method_label, "fgsea::fgsea, adaptive multilevel", fixed = TRUE)
  expect_false(grepl("fgseaSimple", g$method_label, fixed = TRUE))
  expect_identical(g$reference$method, "fgsea::fgsea (fgseaMultilevel)")
  expect_identical(g$params$nperm_simple, 10000L)
  expect_null(g$params[["n_perm"]])           # no fixed permutation count
  expect_null(g$reference$commit)
  expect_null(g$resmat)                       # MetaboAnalyst's matrix is legacy-only
  expect_length(g$deviations, 0L)
  expect_identical(g$nes_note, mmc_gsea_nes_note("fgsea"))
  expect_identical(g$n_positions, 12L)        # one ranked position per EC
  expect_identical(g$ties$n_tied_ec, 0L)
  expect_identical(g$n_p_unavailable, 0L)

  # The same function the python adapter calls, on this engine's own inputs.
  direct <- suppressWarnings(mmc_gsea_fgsea(
    stats::setNames(gs_pathways$emp_cpds, gs_pathways$name), ecs$ec_exp))
  expect_identical(g$log2err_summary, mmc_gsea_log2err_summary(direct))
  # Per pathway too, in the python adapter's shape: named by the table's
  # Pathway. Beside the table, not in it.
  expect_setequal(names(g$log2err), g$table$Pathway)
  expect_equal(g$log2err[direct$pathway], stats::setNames(direct$log2err,
                                                          direct$pathway))
  expect_false("log2err" %in% names(g$table))
  k <- match(g$table$Pathway, direct$pathway)
  expect_setequal(g$table$Pathway, gs_pathways$name)
  expect_equal(g$table[["P.Value"]], direct$pval[k])
  expect_equal(g$table$padj, direct$padj[k])
  expect_equal(g$table$ES, round(direct$ES[k], 4))
  expect_false(is.unsorted(g$table[["P.Value"]], na.rm = TRUE))

  # signed ranking: the five highest-scoring ECs give +1, the five lowest -1
  es <- stats::setNames(g$table$ES, g$table$Pathway)
  expect_equal(es[["top"]], 1)
  expect_equal(es[["bottom"]], -1)
  row <- function(p) g$table[g$table$Pathway == p, , drop = FALSE]
  expect_identical(row("top")[["Leading-edge EmpiricalCompounds"]],
                   "E01; E02; E03; E04; E05")
  # the size columns and each row's own ECs, as the legacy table has them
  expect_identical(row("small")[["Pathway size (model compounds)"]], 2)
  expect_identical(row("small")[["Matched compounds"]], 1)
  expect_identical(row("spread")$EmpiricalCompounds, "E03;E06;E09;E12")
  # one schema for both methods
  legacy <- suppressWarnings(mmc_ma_gsea(ecs, gs_pathways, n_perm = 100,
                                         method = "metaboanalyst_legacy"))
  expect_identical(names(g$table), names(legacy$table))
  expect_false(any(g$table[["ES defaulted by reference implementation"]]))
})

test_that("legacy: the method warns, is recorded, and is never selected by a partial name", {
  expect_warning(g <- mmc_ma_gsea(gs_ecs(), gs_pathways, n_perm = 100,
                                  method = "metaboanalyst_legacy"),
                 "known defects")
  expect_identical(g$method, "metaboanalyst_legacy")
  expect_identical(g$reference$commit, .MMC_GSEA_REF_COMMIT)
  expect_identical(g$nes_note, mmc_gsea_nes_note("metaboanalyst_legacy"))
  expect_null(g$ties)
  expect_error(mmc_ma_gsea(gs_ecs(), gs_pathways, method = "metaboanalyst"),
               "gsea_method")
})

test_that("legacy: GSEA runs on the ECs' scores and the ORA's pathway EC sets", {
  g <- suppressWarnings(mmc_ma_gsea(gs_ecs(), gs_pathways, n_perm = 100,
                                    method = "metaboanalyst_legacy"))
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

test_that("legacy: non-default GSEA settings are reported as deviations", {
  g <- suppressWarnings(mmc_ma_gsea(gs_ecs(), gs_pathways, n_perm = 50, seed = 7,
                                    method = "metaboanalyst_legacy"))
  expect_true(any(grepl("^n_perm = 50", g$deviations)))
  expect_true(any(grepl("^seed = 7", g$deviations)))

  # Under the default method a seed is a setting, not a deviation, and there
  # is no permutation count: n_perm is refused rather than given a meaning.
  std <- suppressWarnings(mmc_ma_gsea(gs_ecs(), gs_pathways, seed = 7))
  expect_length(std$deviations, 0L)
  expect_identical(std$params$seed, 7)
  expect_error(mmc_ma_gsea(gs_ecs(), gs_pathways, n_perm = 50),
               "n_perm is the legacy method")
  expect_null(g$log2err_summary)              # legacy: no fgsea error estimate
  expect_null(g$log2err)
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

  # The method that produced the table travels with it: the default here.
  expect_identical(with$gsea$method, "fgsea")
  expect_identical(g$method, "fgsea")
  expect_identical(g$method_label, with$gsea$method_label)
  expect_identical(g$nes_note, mmc_gsea_nes_note("fgsea"))
  expect_identical(attr(tbl, "ma_info")$gsea$method, "fgsea")
  # ... with its tie summary and its count of tested pathways without a p.
  expect_equal(g$ties$n_ec, with$gsea$ties$n_ec)
  expect_equal(g$ties$n_tied_ec, with$gsea$ties$n_tied_ec)
  expect_identical(g$ties$tie_sensitive_pathways,
                   with$gsea$ties$tie_sensitive_pathways)
  expect_equal(g$n_p_unavailable, with$gsea$n_p_unavailable)
  # ... and how its p-values were estimated, with fgsea's error summary.
  expect_identical(g$reference$method, "fgsea::fgsea (fgseaMultilevel)")
  expect_identical(mmc_gsea_permutation_label(g$params),
                   "adaptive multilevel p-values, nPermSimple = 10000")
  expect_equal(g$log2err_summary, with$gsea$log2err_summary)
  # ... and fgsea's error estimate for every tested pathway, each with its
  # pathway name, restored in the shape the engine returned it.
  expect_identical(names(g$log2err), names(with$gsea$log2err))
  expect_equal(as.numeric(g$log2err), as.numeric(with$gsea$log2err))
  expect_setequal(names(g$log2err), g$table$Pathway)
  expect_length(attr(g$log2err, "pathway_id"), length(g$log2err))
  # The table that the exports and the Shiny payload carry is unchanged.
  expect_false("log2err" %in% names(g$table))
  expect_identical(names(g$table), names(with$gsea$table))

  # A record written while the default method ran fgsea::fgseaSimple() names
  # the same method, "fgsea". It is described by what it stores -- its label,
  # its reference and its permutation count -- and never as multilevel.
  new_file <- files[basename(files) == "ma_run_info.json"]
  rec <- jsonlite::fromJSON(new_file, simplifyVector = FALSE)
  simple_label <- paste0("standard preranked GSEA (fgsea::fgseaSimple, fgsea ",
                         "0.0.0) on the signed EmpiricalCompound scores, one ",
                         "ranking position per EC")
  rec$gsea$method_label <- simple_label
  rec$gsea$reference$method <- "fgsea::fgseaSimple"
  rec$gsea$params$nperm_simple <- NULL
  rec$gsea$params$n_perm <- 100L
  rec$gsea$log2err <- NULL
  jsonlite::write_json(rec, new_file, auto_unbox = TRUE, pretty = TRUE,
                       digits = NA, null = "null")
  was <- attr(read_mummichog_ma_pathways(files), "ma_gsea")
  expect_identical(was$method, "fgsea")
  expect_identical(was$method_label, simple_label)
  expect_identical(was$reference$method, "fgsea::fgseaSimple")
  expect_identical(mmc_gsea_permutation_label(was$params), "100 permutations")
  expect_null(was$log2err_summary)
  # It stored no per-pathway estimate: absent, not "present with gaps".
  rec$gsea$log2err_by_pathway <- NULL
  jsonlite::write_json(rec, new_file, auto_unbox = TRUE, pretty = TRUE,
                       digits = NA, null = "null")
  was <- attr(read_mummichog_ma_pathways(files), "ma_gsea")
  expect_null(was$log2err)
  expect_null(was$log2err_summary)
  # Even without its stored label, such a record is not given today's.
  rec$gsea$method_label <- NULL
  jsonlite::write_json(rec, new_file, auto_unbox = TRUE, pretty = TRUE,
                       digits = NA, null = "null")
  bare <- attr(read_mummichog_ma_pathways(files), "ma_gsea")
  expect_false(grepl("multilevel", bare$method_label, fixed = TRUE))
  expect_match(bare$method_label, "fgsea::fgseaSimple", fixed = TRUE)

  # Legacy, when asked for by name, is recorded and read back as legacy.
  old <- suppressWarnings(do.call(mmc_ma_run_contrast, c(args, list(gsea = list(
    metric_label = "x", method = "metaboanalyst_legacy")))))
  expect_identical(old$stats, without$stats)          # the ORA is still untouched
  old_dir <- tempfile("ma_gsea_legacy_")
  old_files <- write_mummichog_ma_results(old, inp$peaks, old_dir, info = list(),
                                          pathway_ids = NULL)
  og <- attr(read_mummichog_ma_pathways(old_files), "ma_gsea")
  expect_identical(og$method, "metaboanalyst_legacy")
  expect_identical(og$nes_note, mmc_gsea_nes_note("metaboanalyst_legacy"))
  expect_null(og$ties)
  expect_identical(mmc_gsea_permutation_label(og$params), "100 permutations")
  expect_null(og$log2err_summary)
  expect_null(og$log2err)

  # A run record written before the method was recorded names no method. It
  # can only be the legacy method, and must be described as such, whatever the
  # current default is.
  info_file <- old_files[basename(old_files) == "ma_run_info.json"]
  info <- jsonlite::fromJSON(info_file, simplifyVector = FALSE)
  info$gsea$method <- NULL
  info$gsea$method_label <- NULL
  jsonlite::write_json(info, info_file, auto_unbox = TRUE, pretty = TRUE,
                       digits = NA, null = "null")
  pre <- attr(read_mummichog_ma_pathways(old_files), "ma_gsea")
  expect_identical(pre$method, "metaboanalyst_legacy")
  expect_identical(pre$nes_note, mmc_gsea_nes_note("metaboanalyst_legacy"))
  expect_match(pre$method_label, "legacy MetaboAnalyst-compatible", fixed = TRUE)

  # An unknown recorded method is an error, not a silent default.
  info$gsea$method <- "something_else"
  jsonlite::write_json(info, info_file, auto_unbox = TRUE, pretty = TRUE,
                       digits = NA, null = "null")
  expect_error(read_mummichog_ma_pathways(old_files), "unknown GSEA method")
})

test_that("per-pathway log2err survives the run record, missing values and ids included", {
  # Three tested pathways: an estimate, a p-value without one, and no p-value.
  # A pathway name that needs JSON escaping, and one the model has no id for.
  err <- c("Alpha pathway" = 0.125, "Beta \"quoted\", pathway" = NA_real_,
           "Gamma pathway" = NA_real_)
  ids <- c("Alpha pathway" = "map00010", "Beta \"quoted\", pathway" = "map00020")

  entries <- .mmc_ma_log2err_entries(err, ids)
  expect_length(entries, 3L)
  expect_identical(entries[[1]],
                   list(pathway = "Alpha pathway", pathway_id = "map00010",
                        log2err = 0.125))
  expect_null(entries[[2]]$log2err)           # NA is written as JSON null
  expect_null(entries[[3]]$pathway_id)        # no id for this pathway

  # Through the JSON the writer produces and the reader parses.
  path <- withr::local_tempfile(fileext = ".json")
  jsonlite::write_json(list(gsea = list(log2err_by_pathway = entries)), path,
                       auto_unbox = TRUE, pretty = TRUE, digits = NA,
                       null = "null")
  expect_false(any(grepl("\"NA\"", readLines(path), fixed = TRUE)))
  back <- .mmc_ma_log2err_restore(
    jsonlite::fromJSON(path, simplifyVector = TRUE)$gsea[["log2err_by_pathway"]])
  expect_identical(names(back), names(err))
  expect_equal(as.numeric(back), as.numeric(err))
  expect_identical(attr(back, "pathway_id"), c("map00010", "map00020", NA))

  # The same entries unsimplified, as a list of records.
  as_list <- .mmc_ma_log2err_restore(
    jsonlite::fromJSON(path, simplifyVector = FALSE)$gsea[["log2err_by_pathway"]])
  expect_identical(names(as_list), names(err))
  expect_equal(as.numeric(as_list), as.numeric(err))
  expect_identical(attr(as_list, "pathway_id"), c("map00010", "map00020", NA))

  # All values missing, a single pathway, and no ids at all.
  one <- .mmc_ma_log2err_entries(c(Only = NA_real_))
  jsonlite::write_json(list(x = one), path, auto_unbox = TRUE, null = "null")
  back1 <- .mmc_ma_log2err_restore(jsonlite::fromJSON(path)$x)
  expect_identical(names(back1), "Only")
  expect_true(is.numeric(back1) && is.na(back1[[1]]))
  expect_identical(attr(back1, "pathway_id"), NA_character_)

  # Absent is absent: the legacy method, and every record written before the
  # values were stored.
  expect_null(.mmc_ma_log2err_entries(NULL))
  expect_null(.mmc_ma_log2err_restore(NULL))
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
  test_that(sprintf("legacy: GSEA matches MetaboAnalystR end to end (%s)", case), {
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

    # Historical compatibility: the reference is MetaboAnalystR's own output,
    # so it is compared with the legacy method, named explicitly.
    g <- suppressWarnings(mmc_ma_gsea(ecs, pw, n_perm = ref$perm_num,
                                      method = "metaboanalyst_legacy"))
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
