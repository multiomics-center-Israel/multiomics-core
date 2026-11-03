# tests/testthat/test-mummichog-ma-integration.R
#
# MetaboAnalyst-compatible mummichog, component 5: the engine switch, the stage
# module (05d), the result files (06n), the report sections and the DAG
# wiring. The end-to-end tests run the stage on the shared synthetic fixture
# and compare what it WRITES with the pinned MetaboAnalystR's own output
# (fixtures/mummichog_ma_parity/stats_reference.rds). All data synthetic.

ma_cfg <- function(...) {
  list(modes = list(metabolomics = list(organism = "Synthetic organism",
    enrichment = list(mummichog = utils::modifyList(list(
      enabled = TRUE, engine = "metaboanalyst", ionization_mode = "mixed",
      polarity_column = "Reference Ion", p_cutoff = 0.05, tolerance_ppm = 10,
      n_permutations = 100), list(...))))),
    params = list(seed = 1))
}

# pre / de_res / model built from the fixture inputs, in fixture feature order.
ma_stage_inputs <- function(inp) {
  pk <- inp$peaks
  row_data <- data.frame(
    feature_id = pk$feature_id, `m/z` = pk$mz, `RT [min]` = pk$rt,
    `Reference Ion` = ifelse(pk$polarity == "positive", "[M+H]+1", "[M-H]-1"),
    check.names = FALSE, stringsAsFactors = FALSE)
  de <- data.frame(feature_id = pk$feature_id, P.Value = pk$p_value,
                   logFC = pk$statistic, stringsAsFactors = FALSE)
  # A shuffled DE table: the stage must restore feature-table order itself.
  de <- de[rev(seq_len(nrow(de))), ]
  model <- ma_write_model(inp$compounds, inp$pathways,
                          pathway_names = inp$pathway_names)
  list(pre = list(row_data = row_data),
       de_res = list(de_tables = list(HL_vs_LL = de)), model = model)
}

ma_int_ref <- function() {
  ref_file <- file.path(ma_fixture_dir, "stats_reference.rds")
  skip_if_not(file.exists(ref_file), "stats parity reference not generated yet")
  ref <- readRDS(ref_file)
  skip_if_not(identical(Sys.getlocale("LC_COLLATE"), ref$locale_collate),
              "stats reference generated under another LC_COLLATE")
  skip_if_not_installed("fitdistrplus")
  skip_if_not(identical(as.character(utils::packageVersion("fitdistrplus")),
                        ref$fitdistrplus_version),
              "stats reference generated under another fitdistrplus")
  ref
}

ma_run_stage <- function(ref, mode) {
  si <- ma_stage_inputs(ref$inputs)
  cfg <- ma_cfg(ionization_mode = mode, model_json = si$model)
  out <- tempfile("ma_stage_")
  files <- suppressWarnings(
    mod_mummichog_metaboanalyst(si$pre, si$de_res, cfg, out))
  list(files = files, cfg = cfg, si = si, out = out)
}

# ---- engine switch -------------------------------------------------------------

test_that("the engine defaults to metaboanalyst and accepts only the two engines", {
  expect_identical(mmc_engine(list()), "metaboanalyst")
  expect_identical(mmc_engine(ma_cfg(engine = NULL)), "metaboanalyst")
  expect_identical(mmc_engine(ma_cfg()), "metaboanalyst")
  expect_identical(mmc_engine(ma_cfg(engine = "Python")), "python")
  expect_error(mmc_engine(ma_cfg(engine = "mummichog3")), "must be")
})

test_that("an engine-less config written for the python engine stops, not switches", {
  expect_error(mmc_engine(ma_cfg(engine = NULL, ionization_mode = "pos_default")),
               "engine: python")
  expect_error(mmc_engine(ma_cfg(engine = NULL, positive_mode = "positive")),
               "positive_mode")
  # named explicitly, the same settings keep running on the python engine
  expect_identical(mmc_engine(ma_cfg(engine = "python",
                                     ionization_mode = "pos_default")), "python")
})

test_that("a joint mixed run is labelled by its contrast, never split", {
  lab <- "A_vs_B, positive ions"
  expect_identical(mmc_run_label_parts(lab, ma_cfg())$contrast, lab)
  expect_true(is.na(mmc_run_label_parts(lab, ma_cfg())$polarity))
  py <- ma_cfg(engine = "python")
  expect_identical(mmc_run_label_parts(lab, py)$polarity, "positive")
})

test_that("the stage refuses settings with no MetaboAnalyst meaning", {
  si <- list(pre = list(row_data = data.frame(feature_id = "F1")),
             de_res = list(de_tables = list(a = data.frame())))
  expect_error(mod_mummichog_metaboanalyst(
    si$pre, si$de_res, ma_cfg(ionization_mode = "pos_default"), tempdir()),
    "no\\s+MetaboAnalyst equivalent|pos_default")
  human <- ma_cfg()
  human$modes$metabolomics$organism <- "Homo sapiens"   # built-in human_mfn
  expect_error(mod_mummichog_metaboanalyst(si$pre, si$de_res, human, tempdir()),
               "readable metabolic model")
  expect_null(suppressMessages(mod_mummichog_metaboanalyst(
    si$pre, si$de_res, ma_cfg(enabled = FALSE), tempdir())))
})

# ---- DAG wiring ------------------------------------------------------------------

test_that("the DAG runs the MA stage by default and the python stage on request", {
  skip_if_not_installed("targets")
  # pipe_metabolomics() calls tar_target() unqualified, as _targets.R does
  # after attaching targets; attach it for this test only, so the file runs
  # on its own.
  withr::local_package("targets")
  cmd <- function(tl) {
    nm <- vapply(tl, function(t) t$settings$name, character(1))
    paste(deparse(tl[[which(nm == "metab_mummichog_pinned_files")]]$command$expr),
          collapse = " ")
  }
  py <- pipe_metabolomics("pqn", mummichog_enabled = TRUE,
                          mummichog_engine = "python")
  ma <- pipe_metabolomics("pqn", mummichog_enabled = TRUE,
                          mummichog_engine = "metaboanalyst")
  # the default plan is the MetaboAnalyst-compatible stage
  expect_identical(cmd(pipe_metabolomics("pqn", mummichog_enabled = TRUE)), cmd(ma))
  expect_match(cmd(py), "mod_mummichog_pinned\\(")
  expect_false(grepl("metaboanalyst", cmd(py)))
  expect_match(cmd(ma), "mod_mummichog_metaboanalyst\\(")
  names_of <- function(tl) sort(vapply(tl, function(t) t$settings$name, character(1)))
  expect_identical(names_of(ma), names_of(py))
})

# ---- the stage end to end, against MetaboAnalystR's own output -------------------

for (mode in c("mixed", "positive", "negative")) {
  test_that(sprintf("the stage reproduces MetaboAnalystR end to end (%s)", mode), {
    ref <- ma_int_ref()
    up <- ref$cases[[mode]]
    st <- ma_run_stage(ref, mode)
    f <- st$files
    expect_true(all(file.exists(f)))
    expect_true(any(grepl("mummichog_pinned/contrasts\\.tsv$", f)))
    ma_dir <- file.path(st$out, "mummichog_pinned", "HL_vs_LL", "metaboanalyst")
    expect_setequal(basename(f[dirname(f) == ma_dir]),
                    c("ma_pathway_results.tsv",
                      "mummichog_pathway_enrichment_mummichog.csv",
                      "ma_matched_features.tsv", "ma_input.tsv",
                      "ma_run_info.json", "ma_gsea_results.tsv"))

    # MetaboAnalystR's own result CSV, cell for cell.
    csv <- utils::read.csv(file.path(ma_dir, "mummichog_pathway_enrichment_mummichog.csv"),
                           check.names = FALSE, row.names = 1,
                           colClasses = "character")
    expect_identical(rownames(csv), rownames(up$restable))
    expect_identical(unname(as.matrix(csv)), unname(up$restable))

    # The report-facing table: same pathways, order and values.
    tbl <- read_mummichog_pathways_by_contrast(f)[["HL_vs_LL"]]
    res <- up$mummi.resmat
    expect_identical(tbl$pathway, rownames(res))
    for (col in colnames(res)) {
      expect_equal(tbl[[col]], unname(res[, col]), label = col)
    }
    expect_equal(tbl$Holm_adj_FET, suppressWarnings(as.numeric(up$restable[, "AdjP.Fisher"])))
    expect_equal(tbl$Holm_adj_Gamma, suppressWarnings(as.numeric(up$restable[, "AdjP.Gamma"])))
    expect_identical(tbl[["p-value"]], tbl$FET)
    expect_identical(tbl$overlap_size, tbl$Sig.EC)
    expect_identical(tbl$pathway_size, tbl$Total.EC)
    expect_identical(attr(tbl, "ma_info")$params$ionization_mode, mode)
  })
}

test_that("the matched-feature table carries feature ids and polarity", {
  ref <- ma_int_ref()
  st <- ma_run_stage(ref, "mixed")
  feats <- readr::read_tsv(st$files[basename(st$files) == "ma_matched_features.tsv"],
                           show_col_types = FALSE)
  expect_true(all(c("feature_id", "polarity", "Query.Mass", "Matched.Compound",
                    "Matched.Form", "Empirical.Compound") %in% names(feats)))
  expect_true(all(feats$feature_id %in% ref$inputs$peaks$feature_id))
  # one EC with members of both polarities: the joint grouping is visible here
  both <- tapply(feats$polarity, feats$Empirical.Compound,
                 function(p) length(unique(p)) == 2)
  expect_true(any(both))
})

# ---- report sections and the Shiny contract --------------------------------------

test_that("MA sections show FET with Holm columns, the MA GSEA, and no evidence", {
  ref <- ma_int_ref()
  st <- ma_run_stage(ref, "mixed")
  pw <- read_mummichog_pathways_by_contrast(st$files)
  secs <- suppressMessages(build_mummichog_report_sections(
    pw, st$cfg, files = st$files, row_data = st$si$pre$row_data,
    de_res = st$si$de_res))
  s <- secs[["HL_vs_LL"]]
  expect_identical(s$engine, "metaboanalyst")
  # GSEA: the engine's own, on its joint-polarity ECs; never python's
  expect_false(is.null(s$gsea))
  expect_identical(s$gsea$engine, "metaboanalyst")
  expect_true(all(c("Pathway", "NES", "P.Value", "padj", "EmpiricalCompounds") %in%
                    names(s$gsea$table)))
  expect_false(is.null(s$gsea_plot))
  expect_true(all(s$summary$Pathway %in% c(s$table$Pathway, s$gsea$table$Pathway)))
  expect_null(s$evidence)
  expect_true(is.na(s$polarity))
  expect_false(is.null(s$ma_info))
  expect_true(all(c("p.value", "p.value EASE", "p.value Gamma",
                    "p.value empirical", "adj.p Holm FET", "adj.p Holm EASE",
                    "adj.p Holm Gamma") %in% names(s$table)))
  expect_false(any(grepl("BH|FDR|fdr", names(s$table))))
  expect_identical(s$table$Pathway, pw[["HL_vs_LL"]]$pathway)   # MA order kept
  expect_identical(s$table$p.value, pw[["HL_vs_LL"]]$FET)
  # the keys the Shiny payload keeps are all present
  expect_true(all(.METAB_SHINY_MUMMICHOG_KEYS %in% names(s)))
})

test_that("python pathway tables keep their table shape and order", {
  pw <- data.frame(pathway = c("b", "a"), overlap_size = c(1, 2),
                   pathway_size = c(4, 5), `p-value` = c(0.2, 0.01),
                   check.names = FALSE)
  tbl <- build_mummichog_pathway_table(pw)
  expect_named(tbl, c("Pathway", "Overlap", "Pathway size", "Enrichment ratio",
                      "p.value"))
  expect_identical(tbl$Pathway, c("a", "b"))
})
