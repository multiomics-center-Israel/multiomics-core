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

# ---- engine selection and migration -------------------------------------------
# Configs written before the `engine` key ran the python engine, so an enabled
# config must name its engine; it is never guessed or defaulted.

test_that("an explicit engine is honoured whatever the other settings", {
  for (m in c("pos_default", "positive", "negative", "mixed")) {
    expect_identical(mmc_engine(ma_cfg(engine = "python", ionization_mode = m)),
                     "python", label = m)
  }
  expect_identical(mmc_engine(ma_cfg(engine = "python", positive_mode = "positive")),
                   "python")
  for (m in c("positive", "negative", "mixed")) {
    expect_identical(mmc_engine(ma_cfg(engine = "metaboanalyst", ionization_mode = m)),
                     "metaboanalyst", label = m)
  }
  expect_identical(mmc_engine(ma_cfg(engine = " Python ")), "python")
  expect_identical(mmc_engine(ma_cfg(engine = "MetaboAnalyst")), "metaboanalyst")
  expect_error(mmc_engine(ma_cfg(engine = "mummichog3")), "must be")
})

test_that("an enabled config without an engine stops, whatever its polarity", {
  for (m in list("positive", "negative", "mixed", "pos_default", NULL)) {
    cfg <- ma_cfg(engine = NULL)
    cfg$modes$metabolomics$enrichment$mummichog$ionization_mode <- m
    err <- tryCatch(mmc_engine(cfg), error = conditionMessage)
    lab <- if (is.null(m)) "no ionization_mode" else m
    expect_match(err, "no engine is chosen", label = lab)
    expect_match(err, "engine: python", fixed = TRUE, label = lab)
    expect_match(err, "engine: metaboanalyst", fixed = TRUE, label = lab)
    expect_match(err, "reproduces", label = lab)
  }
})

test_that("a disabled config needs no engine", {
  expect_identical(mmc_engine(list()), "metaboanalyst")
  expect_identical(mmc_engine(ma_cfg(engine = NULL, enabled = FALSE)), "metaboanalyst")
  expect_identical(mmc_engine(ma_cfg(engine = NULL, enabled = FALSE,
                                     ionization_mode = "pos_default")), "metaboanalyst")
})

test_that("the template selects the MetaboAnalyst-compatible engine explicitly", {
  tpl <- file.path(normalizePath(if (dir.exists("R")) "." else "../..", mustWork = FALSE),
                   "config", "templates", "metabolomics_template.yaml")
  skip_if_not(file.exists(tpl), "template not found")
  mm <- yaml::read_yaml(tpl)$modes$metabolomics$enrichment$mummichog
  expect_identical(mm$engine, "metaboanalyst")
  mm$enabled <- TRUE
  cfg <- list(modes = list(metabolomics = list(enrichment = list(mummichog = mm))))
  expect_identical(mmc_engine(cfg), "metaboanalyst")
})

test_that("the wizard names the python engine for its built-in human model", {
  run_r <- file.path(normalizePath(if (dir.exists("R")) "." else "../..", mustWork = FALSE),
                     "run.R")
  skip_if_not(file.exists(run_r), "run.R not found")
  src <- paste(readLines(run_r, warn = FALSE), collapse = "\n")
  block <- regmatches(src, regexpr("mummichog_yaml <- '[^']*'", src))
  expect_length(block, 1L)
  expect_match(block, "engine: python", fixed = TRUE)
  expect_match(block, "enabled: true", fixed = TRUE)
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

# ---- failures: GSEA is optional, a failed contrast leaves nothing stale ---------

test_that("a GSEA failure keeps the contrast's ORA and records why", {
  ref <- ma_int_ref()
  skip_if_not_installed("fgsea")
  # Sourced functions, not a package: stub by assignment in their environment.
  target <- environment(mmc_ma_run_contrast)
  orig <- get("mmc_ma_gsea", envir = target)
  assign("mmc_ma_gsea", function(...) stop("synthetic GSEA failure"),
         envir = target)
  withr::defer(assign("mmc_ma_gsea", orig, envir = target))

  si <- ma_stage_inputs(ref$inputs)
  cfg <- ma_cfg(model_json = si$model)
  out <- withr::local_tempdir()
  expect_warning(
    files <- mod_mummichog_metaboanalyst(si$pre, si$de_res, cfg, out),
    "GSEA failed, ORA kept: synthetic GSEA failure")

  expect_true("ma_pathway_results.tsv" %in% basename(files))
  expect_false(any(basename(files) %in% c("ma_gsea_results.tsv", "ma_failed.json")))
  tbl <- read_mummichog_pathways_by_contrast(files)[["HL_vs_LL"]]
  expect_identical(tbl$pathway, rownames(ref$cases$mixed$mummi.resmat))
  expect_null(attr(tbl, "ma_gsea"))
  gi <- attr(tbl, "ma_info")$gsea
  expect_false(gi$ran)
  expect_identical(gi$error, "synthetic GSEA failure")
})

test_that("a failed contrast's earlier results and exports do not survive a rerun", {
  ref <- ma_int_ref()
  si <- ma_stage_inputs(ref$inputs)
  cfg <- ma_cfg(model_json = si$model)
  out <- withr::local_tempdir()
  good <- si$de_res$de_tables$HL_vs_LL
  run <- function(second) {
    de_res <- list(de_tables = list(HL_vs_LL = good, LL_vs_DD = second))
    suppressWarnings(mod_mummichog_metaboanalyst(si$pre, de_res, cfg, out))
  }
  dd_dir <- file.path(out, "mummichog_pinned", "LL_vs_DD", "metaboanalyst")

  run(good)
  expect_true(file.exists(file.path(dd_dir, "ma_pathway_results.tsv")))

  # The rerun breaks LL_vs_DD only; the stage completes on HL_vs_LL.
  broken <- good[, setdiff(names(good), "logFC")]
  expect_warning(
    files <- mod_mummichog_metaboanalyst(
      si$pre, list(de_tables = list(HL_vs_LL = good, LL_vs_DD = broken)), cfg, out),
    "contrast 'LL_vs_DD' failed")
  expect_identical(list.files(dd_dir), "ma_failed.json")
  expect_true(file.path(dd_dir, "ma_failed.json") %in% files)
  rec <- jsonlite::fromJSON(file.path(dd_dir, "ma_failed.json"))
  expect_identical(rec$contrast, "LL_vs_DD")
  expect_match(rec$error, "logFC")
  # It gets no report section; the healthy contrast is untouched.
  expect_identical(names(read_mummichog_pathways_by_contrast(files)), "HL_vs_LL")

  # The previous run's presentation exports of the failed contrast are all
  # cleared -- ORA and evidence too -- and the healthy contrast's are kept.
  dest <- file.path(out, "mummichog_pinned")
  kinds <- c("mummichog_pathway_bubble_%s.png", "mummichog_pathway_table_%s.tsv",
             "mummichog_evidence_pathways_%s.csv", "mummichog_gsea_table_%s.tsv",
             "mummichog_summary_%s.tsv")
  stale <- file.path(dest, sprintf(kinds, "LL_vs_DD"))
  keep  <- file.path(dest, sprintf(kinds, "HL_vs_LL"))
  file.create(c(stale, keep))
  gone <- clear_orphan_mummichog_gsea_exports(
    files, list(HL_vs_LL = list(slug = "HL_vs_LL")), out)
  expect_setequal(gone, stale)
  expect_true(all(file.exists(keep)))

  # Once the contrast succeeds again, its failure record goes.
  run(good)
  expect_false(file.exists(file.path(dd_dir, "ma_failed.json")))
  expect_true(file.exists(file.path(dd_dir, "ma_pathway_results.tsv")))
})

test_that("when every contrast fails, earlier exports are cleared before the stop", {
  ref <- ma_int_ref()
  si <- ma_stage_inputs(ref$inputs)
  cfg <- ma_cfg(model_json = si$model)
  out <- withr::local_tempdir()
  suppressWarnings(mod_mummichog_metaboanalyst(si$pre, si$de_res, cfg, out))

  # The stop below keeps the downstream exports target from running, so the
  # stage itself must clear the previous run's exports.
  dest <- file.path(out, "mummichog_pinned")
  stale <- file.path(dest, c("mummichog_pathway_bubble_HL_vs_LL.pdf",
                             "mummichog_pathway_table_HL_vs_LL.csv",
                             "mummichog_evidence_features_HL_vs_LL.tsv",
                             "mummichog_gsea_scatter_HL_vs_LL.png",
                             "mummichog_summary_HL_vs_LL.tsv"))
  keep <- file.path(dest, "mummichog_pathway_table_OTHER.tsv")
  file.create(c(stale, keep))

  good <- si$de_res$de_tables$HL_vs_LL
  broken <- list(de_tables = list(HL_vs_LL = good[, setdiff(names(good), "logFC")]))
  expect_error(suppressWarnings(
    mod_mummichog_metaboanalyst(si$pre, broken, cfg, out)),
    "every contrast failed")
  expect_false(any(file.exists(stale)))
  expect_true(file.exists(keep))
  ma_dir <- file.path(dest, "HL_vs_LL", "metaboanalyst")
  expect_identical(list.files(ma_dir), "ma_failed.json")
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
