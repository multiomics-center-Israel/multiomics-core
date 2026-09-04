# tests/testthat/test-mummichog-gsea.R
#
# Unit tests for the MetaboAnalyst-style MS peaks-to-pathways GSEA layer (06g):
# the ranking-metric choice, the EmpiricalCompound score construction (signed
# max over member features, with same-m/z features merged first), the
# pathway -> detected-EC gene sets, the fgsea run, the NES direction wording,
# and the summary scatter's encodings.
#
# The EC-score semantics asserted here are MetaboAnalystR's, verified in
# peaks_to_function.R (`ec.exp.vec <- unlist(lapply(ec_exp_dict, max))` in the
# mumRT/v2 branch) — deliberately NOT mean or max(abs()).
#
# All fixtures are synthetic.

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Per-(EC, feature) frame in the shape read_mummichog_ec_features() returns.
ec_features_fixture <- function() {
  data.frame(
    stringsAsFactors = FALSE,
    EID        = c("E1", "E1", "E2", "E2", "E3"),
    input_row  = paste0("row", 1:5),
    feature_id = paste0("f", 1:5),
    mz         = c(100, 200, 300, 300, 400),
    retention_time = c(1, 1, 2, 2, 3),
    p_value    = c(0.01, 0.02, 0.03, 0.04, 0.20),
    statistic  = c(3, -5, 1, 3, -2),
    adduct     = rep("M+H[1+]", 5)
  )
}

feature_scores_fixture <- function() {
  stats::setNames(c(3, -5, 1, 3, -2), paste0("f", 1:5))
}

# A larger synthetic mummichog tree + model, big enough for a real fgsea run:
# 24 EmpiricalCompounds, each with one feature and one candidate compound, and
# three pathways of 8 compounds each. Pathway "PW_pos" holds the 8 highest
# statistics, "PW_neg" the 8 lowest, "PW_mid" the middle 8.
build_gsea_fixture <- function(root, contrast_dir = "HL_vs_LL") {
  n     <- 24L
  cpds  <- sprintf("C%05d", seq_len(n))
  eids  <- sprintf("E%03d", seq_len(n))
  feats <- sprintf("feat_%02d", seq_len(n))
  # statistics from +12 down to -12 (never 0), so the ranking is unambiguous
  stat  <- c(seq(12, 1), seq(-1, -12))
  mzs   <- 100 + seq_len(n)

  tables <- file.path(root, "mummichog_pinned", contrast_dir, "v2",
                      "1700000000.1.run", "tables")
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)

  writeLines(c(
    "EID\tmassfeature_rows\tstr_row_ion\tcompounds\tcompound_names",
    sprintf("%s\trow%d\trow%d_M+H[1+]\t%s\tcpd %s", eids, seq_len(n),
            seq_len(n), cpds, cpds)
  ), file.path(tables, "ListOfEmpiricalCompounds.tsv"))

  writeLines(c(
    "input_row\tEID\tstr_row_ion\tcompounds\tcompound_names\tinput_row\tm/z\tretention_time\tp_value\tstatistic\tCompoundID_from_user",
    sprintf("row%d\t%s\trow%d_M+H[1+]\t%s\tcpd %s\trow%d\t%s\t1.0\t0.05\t%s\t%s",
            seq_len(n), eids, seq_len(n), cpds, cpds, seq_len(n),
            format(mzs), format(stat), feats)
  ), file.path(tables, "userInput_to_EmpiricalCompounds.tsv"))

  writeLines(c(
    "pathway\toverlap_size\tpathway_size\tp-value\toverlap_EmpiricalCompounds (id)\toverlap_features (id)\toverlap_features (name)",
    sprintf("PW_pos\t3\t8\t0.01\t%s\t\t", paste(eids[1:3], collapse = ",")),
    sprintf("PW_neg\t3\t8\t0.02\t%s\t\t", paste(eids[22:24], collapse = ","))
  ), file.path(tables, "mcg_pathwayanalysis_HL_vs_LL.tsv"))

  readr::write_tsv(
    data.frame(dir = contrast_dir, contrast = "HL_vs_LL",
               stringsAsFactors = FALSE),
    file.path(root, "mummichog_pinned", "contrasts.tsv"))

  # Model: three disjoint pathways over the same 24 compounds, in native
  # mummichog-2 shape (the compact one).
  model_path <- file.path(root, "model.json")
  jsonlite::write_json(list(
    metabolic_pathways = list(
      list(id = "p1", name = "PW_pos", cpds = as.list(cpds[1:8])),
      list(id = "p2", name = "PW_mid", cpds = as.list(cpds[9:16])),
      list(id = "p3", name = "PW_neg", cpds = as.list(cpds[17:24]))
    ),
    dict_cpds_def = stats::setNames(as.list(paste("cpd", cpds)), cpds)
  ), model_path, auto_unbox = TRUE)

  list(
    files    = list_mummichog_files(root),
    model    = read_mummichog_model_pathways(model_path),
    de_table = data.frame(feature_id = feats, logFC = stat / 4,
                          statistic = stat, P.Value = rep(0.05, n),
                          stringsAsFactors = FALSE),
    eids     = eids,
    stat     = stat
  )
}

# A GSEA result table in the shape run_mummichog_gsea()$table returns.
gsea_table_fixture <- function() {
  data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    "Pathway"                         = c("PW_pos", "PW_neg", "PW_mid"),
    "Pathway size (model compounds)"  = c(8, 8, 8),
    "Detected ECs"                    = c(8, 8, 8),
    "Hits used in GSEA"               = c(8, 8, 8),
    "ES"                              = c(0.80, -0.75, 0.10),
    "NES"                             = c(2.10, -1.90, 0.00),
    "P.Value"                         = c(0.001, 0.004, 0.700),
    "padj"                            = c(0.003, 0.006, 0.700),
    "Direction"                       = c("toward HL", "toward LL", "toward HL"),
    "Leading-edge EmpiricalCompounds" = c("E001;E002", "E024;E023", NA)
  )
}


# ---------------------------------------------------------------------------
# ranking statistic
# ---------------------------------------------------------------------------

test_that("the GSEA ranking metric prefers the moderated t statistic", {
  de <- data.frame(feature_id = c("f1", "f2"), logFC = c(1, -1),
                   statistic = c(4, -6), stringsAsFactors = FALSE)
  r <- mmc_gsea_rank_metric(de)

  expect_identical(r$metric, "moderated_t")
  expect_equal(unname(r$values[c("f1", "f2")]), c(4, -6))
  expect_match(r$label, "moderated t")
})

test_that("the GSEA ranking metric falls back to logFC without a statistic", {
  de <- data.frame(feature_id = c("f1", "f2"), logFC = c(1, -1),
                   stringsAsFactors = FALSE)
  r <- mmc_gsea_rank_metric(de)

  expect_identical(r$metric, "logFC")
  expect_equal(unname(r$values[c("f1", "f2")]), c(1, -1))
  expect_match(r$label, "no t statistic")

  # an all-NA statistic column is not usable either
  de2 <- transform(de, statistic = c(NA_real_, NA_real_))
  expect_identical(mmc_gsea_rank_metric(de2)$metric, "logFC")
})

test_that("the GSEA ranking metric is NULL when nothing is usable", {
  expect_null(mmc_gsea_rank_metric(NULL))
  expect_null(mmc_gsea_rank_metric(data.frame()))
  expect_null(mmc_gsea_rank_metric(data.frame(feature_id = "f1",
                                              AveExpr = 1)))
})

test_that("a literal 't' column is accepted as the moderated statistic", {
  de <- data.frame(feature_id = c("f1", "f2"), logFC = c(1, -1),
                   t = c(2, -3), stringsAsFactors = FALSE)
  r <- mmc_gsea_rank_metric(de)
  expect_identical(r$metric, "moderated_t")
  expect_equal(unname(r$values[["f2"]]), -3)
})


# ---------------------------------------------------------------------------
# EmpiricalCompound scores (MetaboAnalyst v2 semantics)
# ---------------------------------------------------------------------------

test_that("an EC score is the SIGNED maximum of its member features", {
  ec <- mmc_ec_scores(ec_features_fixture(), feature_scores_fixture())

  # E1 = max(3, -5) = 3 — not the mean (-1) and not max(abs()) (5 / -5)
  expect_equal(unname(ec[["E1"]]), 3)
  expect_false(isTRUE(all.equal(unname(ec[["E1"]]), -1)))
  expect_false(isTRUE(all.equal(unname(ec[["E1"]]), 5)))
  # E3 keeps its negative score (no absolute value anywhere)
  expect_equal(unname(ec[["E3"]]), -2)
})

test_that("features sharing one m/z are merged before the EC max", {
  # f3 and f4 both sit at m/z 300 with scores 1 and 3 -> merged to mean 2,
  # so E2 scores 2 rather than 3.
  ec <- mmc_ec_scores(ec_features_fixture(), feature_scores_fixture())
  expect_equal(unname(ec[["E2"]]), 2)

  # the merge rule is selectable, matching MetaboAnalyst's rank.metric
  expect_equal(unname(mmc_ec_scores(ec_features_fixture(),
                                    feature_scores_fixture(),
                                    rank_metric = "max")[["E2"]]), 3)
  expect_equal(unname(mmc_ec_scores(ec_features_fixture(),
                                    feature_scores_fixture(),
                                    rank_metric = "min")[["E2"]]), 1)
})

test_that("the EC score vector is sorted decreasing and skips unscored features", {
  ec <- mmc_ec_scores(ec_features_fixture(), feature_scores_fixture())
  expect_identical(names(ec), c("E1", "E2", "E3"))
  expect_true(all(diff(unname(ec)) <= 0))

  # a feature with no DE score simply does not contribute
  partial <- feature_scores_fixture()[c("f1", "f5")]
  ec2 <- mmc_ec_scores(ec_features_fixture(), partial)
  expect_setequal(names(ec2), c("E1", "E3"))
  expect_equal(unname(ec2[["E1"]]), 3)
})

test_that("EC scores are empty when there is nothing to score", {
  expect_length(mmc_ec_scores(NULL, feature_scores_fixture()), 0)
  expect_length(mmc_ec_scores(ec_features_fixture(),
                              stats::setNames(numeric(0), character(0))), 0)
  expect_length(mmc_ec_scores(ec_features_fixture()[0, ],
                              feature_scores_fixture()), 0)
})


# ---------------------------------------------------------------------------
# pathway -> detected EC sets
# ---------------------------------------------------------------------------

test_that("pathway EC sets collect every EC with a candidate in the pathway", {
  cand <- data.frame(
    stringsAsFactors = FALSE,
    EID         = c("E1", "E1", "E2", "E3"),
    compound_id = c("C1", "C2", "C2", "C9"),
    compound_name = c("a", "b", "b", "i")
  )
  sets <- mmc_pathway_ec_sets(list(P1 = c("C1", "C3"), P2 = c("C2"),
                                   P3 = c("C7")),
                              cand, ec_ids = c("E1", "E2", "E3"))

  expect_setequal(names(sets), c("P1", "P2"))       # P3 has no detected EC
  expect_setequal(sets$P1, "E1")
  expect_setequal(sets$P2, c("E1", "E2"))           # both ECs carry C2
})

test_that("pathway EC sets are restricted to the ranked universe", {
  cand <- data.frame(EID = c("E1", "E2"), compound_id = c("C1", "C1"),
                     compound_name = c("a", "a"), stringsAsFactors = FALSE)
  sets <- mmc_pathway_ec_sets(list(P1 = "C1"), cand, ec_ids = "E1")
  expect_setequal(sets$P1, "E1")
  expect_length(mmc_pathway_ec_sets(list(P1 = "C1"), cand, character(0)), 0)
  expect_length(mmc_pathway_ec_sets(list(), cand, "E1"), 0)
})


# ---------------------------------------------------------------------------
# direction wording
# ---------------------------------------------------------------------------

test_that("the NES direction note carries the contrast orientation", {
  note <- mmc_gsea_direction_note("HL_vs_LL", "moderated t statistic")
  expect_match(note, "positive NES")
  expect_match(note, "higher moderated t statistic in HL relative to LL")
  # never a bare "up"/"down"
  expect_false(grepl("\\b(up|down)regulated\\b", note))
})

test_that("an unparseable contrast gets neutral statistic-based wording", {
  note <- mmc_gsea_direction_note(NULL, "log2 fold change")
  expect_match(note, "positive end of the ranked log2 fold change")
  expect_match(mmc_gsea_direction_note("weird-label", "logFC"),
               "positive end of the ranked")
})


# ---------------------------------------------------------------------------
# end-to-end GSEA run
# ---------------------------------------------------------------------------

test_that("GSEA runs over the full ranked EC list and preserves NES sign", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  res <- run_mummichog_gsea(f$files, f$de_table, f$model,
                            contrast = "HL_vs_LL", n_perm = 500, seed = 42)
  expect_false(is.null(res))

  # the FULL EC universe is ranked, not only the ORA-significant overlap
  expect_equal(res$n_ec, 24L)
  expect_identical(res$metric, "moderated_t")
  expect_setequal(res$table$Pathway, c("PW_pos", "PW_mid", "PW_neg"))

  expect_true(all(c("Pathway", "Pathway size (model compounds)", "Detected ECs",
                    "Hits used in GSEA", "ES", "NES", "P.Value", "padj",
                    "Direction", "Leading-edge EmpiricalCompounds") %in%
                    names(res$table)))

  # PW_pos holds the 8 highest statistics -> positive NES; PW_neg the lowest.
  nes <- stats::setNames(res$table$NES, res$table$Pathway)
  expect_gt(nes[["PW_pos"]], 0)
  expect_lt(nes[["PW_neg"]], 0)
  # direction is expressed in the contrast's own terms
  dir <- stats::setNames(res$table$Direction, res$table$Pathway)
  expect_identical(dir[["PW_pos"]], "toward HL")
  expect_identical(dir[["PW_neg"]], "toward LL")

  # leading edge is reported in EC space
  le <- res$table[["Leading-edge EmpiricalCompounds"]][
    res$table$Pathway == "PW_pos"]
  expect_true(grepl("^E0", le))

  # sorted by ascending p, and the size columns are distinct concepts
  expect_equal(res$table[["P.Value"]], sort(res$table[["P.Value"]]))
  expect_true(all(res$table[["Pathway size (model compounds)"]] == 8))
})

test_that("GSEA is reproducible for a fixed seed", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  a <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "HL_vs_LL",
                          n_perm = 500, seed = 7)
  b <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "HL_vs_LL",
                          n_perm = 500, seed = 7)
  expect_equal(a$table[["P.Value"]], b$table[["P.Value"]])
  expect_equal(a$table$NES, b$table$NES)
})

test_that("GSEA skips gracefully when a prerequisite is missing", {
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  # no model
  expect_message(expect_null(
    run_mummichog_gsea(f$files, f$de_table, NULL, contrast = "c")),
    "no metabolic model")
  # no usable ranking statistic
  expect_message(expect_null(
    run_mummichog_gsea(f$files, data.frame(feature_id = "x", AveExpr = 1),
                       f$model, contrast = "c")),
    "no usable ranking statistic")
  # no EC tables among the files
  empty <- withr::local_tempdir()
  expect_message(expect_null(
    run_mummichog_gsea(list_mummichog_files(empty), f$de_table, f$model,
                       contrast = "c")),
    "no EmpiricalCompound tables")
})


# ---------------------------------------------------------------------------
# the summary scatter
# ---------------------------------------------------------------------------

test_that("the GSEA scatter maps NES on x, -log10(p) on y and NES to colour", {
  g <- gsea_table_fixture()
  p <- plot_mummichog_gsea_scatter(g, title = "T", subtitle = "S",
                                   p_cutoff = 0.05)

  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "T")
  expect_identical(p$labels$subtitle, "S")
  expect_match(p$labels$x, "NES")
  expect_match(p$labels$y, "-log10\\(GSEA p-value\\)")
  expect_identical(p$labels$fill, "NES")

  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomPoint" %in% geoms)
  expect_true("GeomHline" %in% geoms)   # the p = 0.05 reference line

  b  <- ggplot2::ggplot_build(p)
  pd <- b$data[[which(geoms == "GeomPoint")[1]]]
  expect_equal(sort(pd$x), sort(g$NES))
  expect_equal(sort(pd$y), sort(-log10(g[["P.Value"]])), tolerance = 1e-6)
})

test_that("the GSEA scatter draws the p = 0.05 reference line", {
  p <- plot_mummichog_gsea_scatter(gsea_table_fixture(), title = "T",
                                   p_cutoff = 0.05)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  b  <- ggplot2::ggplot_build(p)
  hl <- b$data[[which(geoms == "GeomHline")[1]]]

  expect_equal(unique(hl$yintercept), -log10(0.05))
  expect_identical(unique(hl$linetype), "dashed")

  # a different cutoff moves the line
  p2 <- plot_mummichog_gsea_scatter(gsea_table_fixture(), title = "T",
                                    p_cutoff = 0.01)
  geoms2 <- vapply(p2$layers, function(l) class(l$geom)[1], character(1))
  hl2 <- ggplot2::ggplot_build(p2)$data[[which(geoms2 == "GeomHline")[1]]]
  expect_equal(unique(hl2$yintercept), -log10(0.01))
})

test_that("the GSEA scatter's colour scale is diverging and centred at NES = 0", {
  g <- gsea_table_fixture()
  p <- plot_mummichog_gsea_scatter(g, title = "T")
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  pd <- ggplot2::ggplot_build(p)$data[[which(geoms == "GeomPoint")[1]]]

  # NES == 0 lands exactly on the scale's midpoint colour
  mid <- pd$fill[which(abs(pd$x) < 1e-12)]
  expect_equal(toupper(substr(mid, 1, 7)), "#FFFEE0")
  # positive and negative NES sit on opposite arms of the scale
  expect_false(identical(pd$fill[which.max(pd$x)], pd$fill[which.min(pd$x)]))
})

test_that("the GSEA scatter's point size is MetaboAnalyst's sqrt(-log10 p)", {
  g <- gsea_table_fixture()
  p <- plot_mummichog_gsea_scatter(g, title = "T")

  expect_equal(p$data$significance_size, sqrt(-log10(g[["P.Value"]])),
               tolerance = 1e-9)
  # the legend names the quantity honestly and never leaks the source's
  # internal variable name
  size_name <- p$scales$get_scales("size")$name
  expect_match(size_name, "Significance")
  expect_match(size_name, "log10 p")
  expect_false(grepl("radi", size_name, ignore.case = TRUE))
})

test_that("the GSEA scatter returns NULL on empty / unusable input", {
  expect_null(plot_mummichog_gsea_scatter(NULL, title = "T"))
  expect_null(plot_mummichog_gsea_scatter(data.frame(), title = "T"))
  expect_null(plot_mummichog_gsea_scatter(
    data.frame(Pathway = "x", NES = 1, check.names = FALSE), title = "T"))
})


# ---------------------------------------------------------------------------
# report-section integration
# ---------------------------------------------------------------------------

test_that("report sections carry the GSEA scatter as the pathway plot", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)
  cfg  <- list(modes = list(metabolomics = list(
    enrichment = list(mummichog = list(p_cutoff = 0.05,
                                       gsea_permutations = 500)))))
  pw   <- read_mummichog_pathways_by_contrast(f$files)

  secs <- build_mummichog_report_sections(
    pw, cfg,
    files    = f$files,
    de_res   = list(de_tables = list("HL_vs_LL" = f$de_table)),
    row_data = data.frame(feature_id = f$de_table$feature_id,
                          Name = paste("cpd", f$de_table$feature_id),
                          stringsAsFactors = FALSE),
    model    = f$model
  )

  expect_named(secs, "HL_vs_LL")
  s <- secs[["HL_vs_LL"]]
  expect_identical(s$plot_kind, "gsea_scatter")
  expect_s3_class(s$plot, "ggplot")
  expect_match(s$plot$labels$x, "NES")
  # the ORA table and plot are kept alongside, never replaced
  expect_s3_class(s$table, "data.frame")
  expect_true("p.value" %in% names(s$table))
  expect_s3_class(s$ora_plot, "ggplot")
  expect_false(is.null(s$gsea))
  expect_false(is.null(s$evidence))
})

test_that("report sections fall back to the ORA plot when GSEA cannot run", {
  # No mummichog files at all: the evidence and GSEA layers have nothing to read,
  # so the section must still render the ORA plot + table (this is also the
  # "evidence files missing" case the report has to survive).
  cfg <- list(modes = list(metabolomics = list(
    enrichment = list(mummichog = list(p_cutoff = 0.05)))))
  pw  <- data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    pathway      = c("A", "B"),
    overlap_size = c(2, 1),
    pathway_size = c(4, 5),
    "p-value"    = c(0.01, 0.30),
    "overlap_EmpiricalCompounds (id)" = c("E1,E2", "E3")
  )
  secs <- build_mummichog_report_sections(list("HL_vs_LL" = pw), cfg)

  s <- secs[["HL_vs_LL"]]
  expect_identical(s$plot_kind, "ora_bubble")
  expect_s3_class(s$plot, "ggplot")
  expect_match(s$plot$labels$x, "Enrichment ratio")
  expect_null(s$gsea)
  expect_null(s$gsea_plot)
  expect_null(s$evidence)
  expect_s3_class(s$table, "data.frame")
})

test_that("exports write the GSEA and evidence artefacts alongside the ORA ones", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)
  out  <- withr::local_tempdir()

  gsea <- run_mummichog_gsea(f$files, f$de_table, f$model,
                             contrast = "HL_vs_LL", n_perm = 200, seed = 42)
  pw   <- read_mummichog_pathways(f$files)
  annot <- normalize_metab_annotation(
    data.frame(feature_id = f$de_table$feature_id,
               Name = paste("cpd", sprintf("C%05d", seq_len(24))),
               stringsAsFactors = FALSE))
  ev <- build_mummichog_pathway_evidence(pw, f$files, f$model, annot)

  paths <- save_mummichog_exports(
    plot = NULL, table = build_mummichog_pathway_table(pw), out_dir = out,
    contrast_label = "HL_vs_LL", gsea = gsea, evidence = ev)

  base <- basename(paths)
  expect_true("mummichog_pathway_table_HL_vs_LL.tsv" %in% base)
  expect_true("mummichog_gsea_table_HL_vs_LL.tsv" %in% base)
  expect_true("mummichog_evidence_pathways_HL_vs_LL.tsv" %in% base)
  expect_true("mummichog_evidence_empirical_compounds_HL_vs_LL.tsv" %in% base)
  expect_true("mummichog_evidence_features_HL_vs_LL.tsv" %in% base)
  expect_true(all(file.exists(paths)))
})
