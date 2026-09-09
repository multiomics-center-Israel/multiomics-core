# tests/testthat/test-mummichog-gsea.R
#
# Unit tests for the MetaboAnalyst-style MS peaks-to-pathways GSEA layer (06g):
# the ranking-metric choice, the EmpiricalCompound score construction (signed
# max over member features, with same-m/z features merged first), the
# pathway -> detected-EC gene sets, the MetaboAnalyst-equivalent engine, the
# NES interpretation, and the summary scatter's encodings.
#
# Parity against the pinned upstream implementation lives in
# test-mummichog-gsea-parity.R; this file covers the surrounding contract.
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

# A GSEA result table in the shape run_mummichog_gsea()$table returns.
gsea_table_fixture <- function() {
  data.frame(
    check.names = FALSE, stringsAsFactors = FALSE,
    "Pathway"                         = c("PW_pos", "PW_neg", "PW_mid"),
    "Pathway size (model compounds)"  = c(8, 8, 8),
    "Detected ECs"                    = c(8, 8, 8),
    "Tested size"                     = c(8, 8, 8),
    "ES"                              = c(0.80, -0.75, 0.10),
    "NES"                             = c(2.10, -1.90, 0.00),
    "P.Value"                         = c(0.001, 0.004, 0.700),
    "padj"                            = c(0.003, 0.006, 0.700),
    "Leading-edge EmpiricalCompounds" = c("E001; E002", "E024; E023", NA)
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
# NES interpretation (must NOT be a biological direction)
# ---------------------------------------------------------------------------

test_that("the NES note is conservative: parity only, no direction claim", {
  note <- mmc_gsea_nes_note()

  # what it MUST say
  expect_match(note, "retained to reproduce the pinned MetaboAnalystR result",
               fixed = TRUE)
  expect_match(note, "transforms and reorders the score vector after pathway positions are constructed",
               fixed = TRUE)
  expect_match(note, "should not be interpreted as biological up/down direction",
               fixed = TRUE)
  expect_match(note, "magnitude of the pathway members' original scores",
               fixed = TRUE)

  # what it must NOT say: pathway indices come from the SIGNED ordering and are
  # then scored against a transformed, re-sorted vector, so even a claim about
  # which |score| END a pathway sits at is unsupported.
  expect_false(grepl("high-|score| end", note, fixed = TRUE))
  expect_false(grepl("low-|score| end", note, fixed = TRUE))
  expect_false(grepl("enrichment toward", note, fixed = TRUE))
})

test_that("no biological-direction wording survives anywhere in the GSEA layer", {
  # Regression guard: MetaboAnalyst ranks on |score|, so NES sign must never be
  # presented as treatment direction. This locks out the wording that used to be
  # here (a `Direction` column and a numerator/denominator note) and anything
  # equivalent creeping back in.
  expect_false(exists("mmc_gsea_direction_note"))

  root <- normalizePath(if (dir.exists("R")) "." else "../..", mustWork = FALSE)
  targets <- c(
    file.path(root, "R", "domain", "metabolomics", "06g_mummichog_gsea.R"),
    file.path(root, "R", "domain", "metabolomics", "06e_mummichog_plots.R"),
    file.path(root, "R", "pipeline", "metabolomics", "templates",
              "report_metabolomics.Rmd")
  )
  banned <- c("toward numerator", "toward denominator",
              "up-regulated pathway", "down-regulated pathway",
              "upregulated", "downregulated",
              "direction_note", "mmc_gsea_direction_note",
              # over-strong magnitude claims: pathway positions come from the
              # signed ordering, so even an |score|-end reading is unsupported
              "high-|score| end", "low-|score| end",
              "enrichment toward the high", "enrichment toward the low")
  for (f in targets) {
    skip_if_not(file.exists(f), paste("missing", basename(f)))
    txt <- tolower(paste(readLines(f, warn = FALSE), collapse = "\n"))
    for (b in banned) {
      expect_false(grepl(tolower(b), txt, fixed = TRUE),
                   info = paste0("banned direction wording '", b,
                                 "' found in ", basename(f)))
    }
    # the only sanctioned NES sentence, or none at all
    expect_false(grepl("toward hl|toward ll", txt, fixed = FALSE),
                 info = paste("group-direction wording in", basename(f)))
  }
})

test_that("the GSEA result table carries no Direction column", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)
  res  <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "HL_vs_LL")

  expect_false("Direction" %in% names(res$table))
  expect_false(any(grepl("direction", tolower(names(res$table)))))
  expect_null(res$direction_note)
  expect_identical(res$nes_note, mmc_gsea_nes_note())
  # the contrast is carried as provenance only
  expect_identical(res$contrast, "HL_vs_LL")
})


# ---------------------------------------------------------------------------
# end-to-end GSEA run
# ---------------------------------------------------------------------------

test_that("GSEA runs over the full ranked EC list and preserves NES sign", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  res <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "HL_vs_LL")
  expect_false(is.null(res))

  # the FULL EC universe is ranked, not only the ORA-significant overlap
  expect_equal(res$n_ec, 24L)
  expect_identical(res$metric, "moderated_t")
  expect_setequal(res$table$Pathway, c("PW_pos", "PW_mid", "PW_neg"))

  expect_true(all(c("Pathway", "Pathway size (model compounds)", "Detected ECs",
                    "Tested size", "ES", "NES", "P.Value", "padj",
                    "Leading-edge EmpiricalCompounds") %in% names(res$table)))

  # signed ES/NES are preserved (not flattened to absolute values)
  expect_true(any(res$table$NES < 0) || any(res$table$ES < 0))
  expect_true(all(is.finite(res$table$NES)))

  # engine settings are the MetaboAnalyst-equivalent defaults
  expect_identical(res$params$n_perm, 100L)
  expect_equal(res$params$gsea_param, 1)
  expect_equal(res$params$min_size, 1)
  expect_equal(res$params$max_size, Inf)
  expect_identical(res$params$seed, 123L)

  # sorted by ascending raw p, and the size columns are distinct concepts
  expect_equal(res$table[["P.Value"]], sort(res$table[["P.Value"]]))
  expect_true(all(res$table[["Pathway size (model compounds)"]] == 8))
  expect_true(all(res$table[["Detected ECs"]] <= 8))

  # BH is taken over the tested pathways only
  expect_equal(res$table$padj,
               stats::p.adjust(res$table[["P.Value"]], method = "fdr"),
               tolerance = 1e-12)
})

test_that("GSEA is reproducible for a fixed seed", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  a <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "HL_vs_LL")
  b <- run_mummichog_gsea(f$files, f$de_table, f$model, contrast = "HL_vs_LL")
  expect_equal(a$table[["P.Value"]], b$table[["P.Value"]])
  expect_equal(a$table$NES, b$table$NES)
  expect_equal(a$table$ES, b$table$ES)
  # and the default seed is upstream's, not an invented one
  expect_identical(a$params$seed, 123L)
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
  expect_identical(
    p$scales$get_scales("fill")$name,
    "NES"
  )

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
    enrichment = list(mummichog = list(p_cutoff = 0.05)))))
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
  # the ORA TABLE and its evidence are kept alongside, never replaced
  expect_s3_class(s$table, "data.frame")
  expect_true("p.value" %in% names(s$table))
  expect_false(is.null(s$gsea))
  expect_false(is.null(s$evidence))
  # ...but the ORA bubble does NOT compete as a second primary plot
  expect_null(s$ora_plot)
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
  # in the fallback case the ORA bubble IS the primary plot
  expect_s3_class(s$ora_plot, "ggplot")
})

test_that("exports write the GSEA and evidence artefacts alongside the ORA ones", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)
  out  <- withr::local_tempdir()

  gsea <- run_mummichog_gsea(f$files, f$de_table, f$model,
                             contrast = "HL_vs_LL")
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
