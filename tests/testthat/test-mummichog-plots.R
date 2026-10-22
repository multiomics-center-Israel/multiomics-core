# tests/testthat/test-mummichog-plots.R
#
# Unit tests for the mummichog report presentation builders (06e):
# plot_mummichog_bubble(), build_mummichog_pathway_table(), and
# mummichog_report_titles(). Object-introspection style (no vdiffr / ggsave),
# matching test-plot-volcano.R; fixture carries the literal hyphenated `p-value`
# column that read_mummichog_pathways() produces.

# A small pathway table shaped like read_mummichog_pathways() output.
make_pw <- function() {
  data.frame(
    check.names      = FALSE,
    stringsAsFactors = FALSE,
    pathway        = c("Tryptophan metabolism", "Ether lipid metabolism",
                       "Aminoacyl-tRNA biosynthesis"),
    overlap_size   = c(5, 3, 8),
    pathway_size   = c(10, 6, 12),
    "p-value"      = c(0.02, 0.20, 0.60),
    "overlap_EmpiricalCompounds (id)" = c("E1", "E2", "E3")
  )
}

test_that("plot_mummichog_bubble builds a ggplot with the mockup mappings", {
  pw <- make_pw()
  p  <- plot_mummichog_bubble(pw, title = "T", subtitle = "S", p_cutoff = 0.05)

  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "T")
  expect_identical(p$labels$subtitle, "S")
  expect_match(p$labels$x, "Enrichment ratio")
  expect_match(p$labels$y, "-log10")

  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomHline" %in% geoms)   # the p = 0.05 reference line
  expect_true("GeomPoint" %in% geoms)   # the bubbles

  b  <- ggplot2::ggplot_build(p)
  pd <- b$data[[which(geoms == "GeomPoint")[1]]]
  expect_equal(sort(pd$x), sort(pw$overlap_size / pw$pathway_size))
  expect_equal(sort(pd$y), sort(-log10(pw[["p-value"]])), tolerance = 1e-6)
})

test_that("plot_mummichog_bubble returns NULL on empty / unusable input", {
  expect_null(plot_mummichog_bubble(NULL, title = "T"))
  expect_null(plot_mummichog_bubble(data.frame(), title = "T"))
  # missing the required overlap_size / pathway_size columns
  expect_null(plot_mummichog_bubble(data.frame(pathway = "x", `p-value` = 0.1,
                                               check.names = FALSE), title = "T"))
})

test_that("build_mummichog_pathway_table sorts by p and adds enrichment ratio", {
  pw <- make_pw()
  tb <- build_mummichog_pathway_table(pw)

  expect_s3_class(tb, "data.frame")
  expect_identical(colnames(tb),
                   c("Pathway", "Overlap", "Pathway size",
                     "Enrichment ratio", "p.value"))
  expect_equal(tb[["p.value"]], sort(pw[["p-value"]]))         # ascending by p
  expect_equal(tb[["Enrichment ratio"]][1], round(5 / 10, 3)) # smallest-p row first
  expect_equal(nrow(tb), nrow(pw))
})

test_that("build_mummichog_pathway_table returns NULL on empty input", {
  expect_null(build_mummichog_pathway_table(NULL))
  expect_null(build_mummichog_pathway_table(data.frame()))
})

test_that("mummichog_report_titles composes title/subtitle by model precedence", {
  de <- list(de_tables = list("LL_vs_HL" = data.frame(x = 1)))

  # model_ref URL -> leading token of the basename ("cre")
  cfg_ref <- list(modes = list(metabolomics = list(
    organism   = "Coelastrella",
    enrichment = list(mummichog = list(model_ref = list(
      url    = "https://x/releases/cre_kegg_20260711.json",
      sha256 = paste(rep("a", 64), collapse = "")))))))
  t1 <- mummichog_report_titles(cfg_ref, de)
  expect_identical(t1$title, "Mummichog pathway analysis — Coelastrella (cre model)")
  expect_identical(t1$subtitle, "LL vs HL, all features")

  # model_json -> basename without extension
  cfg_json <- list(modes = list(metabolomics = list(
    organism   = "Mus musculus",
    enrichment = list(mummichog = list(model_json = "/data/mmu_model.json")))))
  t2 <- mummichog_report_titles(cfg_json, de)
  expect_match(t2$title, "Mus musculus")
  expect_match(t2$title, "\\(mmu_model model\\)")

  # no organism, no custom model -> built-in human_mfn, no "— organism"
  cfg_builtin <- list(modes = list(metabolomics = list(
    enrichment = list(mummichog = list(enabled = TRUE)))))
  t3 <- mummichog_report_titles(cfg_builtin, de)
  expect_identical(t3$title, "Mummichog pathway analysis (human_mfn model)")

  # no contrast -> NULL subtitle
  expect_null(mummichog_report_titles(cfg_builtin, list())$subtitle)
})

test_that("save_mummichog_exports writes the table as TSV + CSV", {
  tmp <- withr::local_tempdir()
  tb  <- build_mummichog_pathway_table(make_pw())

  paths <- save_mummichog_exports(NULL, tb, tmp, contrast_label = "LL vs HL")

  tsv <- file.path(tmp, "mummichog_pinned", "mummichog_pathway_table_LL_vs_HL.tsv")
  csv <- file.path(tmp, "mummichog_pinned", "mummichog_pathway_table_LL_vs_HL.csv")
  expect_setequal(paths, c(tsv, csv))
  expect_true(file.exists(tsv))
  expect_true(file.exists(csv))
  # round-trips with the same rows
  back <- readr::read_tsv(tsv, show_col_types = FALSE)
  expect_equal(nrow(back), nrow(tb))
})

test_that("save_mummichog_exports writes the plot as PNG + PDF (when a device exists)", {
  skip_if_not(isTRUE(capabilities("png")), "no PNG graphics device")
  tmp <- withr::local_tempdir()
  p   <- plot_mummichog_bubble(make_pw(), title = "T")

  paths <- save_mummichog_exports(p, NULL, tmp, contrast_label = "LL_vs_HL")

  png <- file.path(tmp, "mummichog_pinned", "mummichog_pathway_bubble_LL_vs_HL.png")
  pdf <- file.path(tmp, "mummichog_pinned", "mummichog_pathway_bubble_LL_vs_HL.pdf")
  expect_true(all(c(png, pdf) %in% paths))
  expect_true(file.exists(png))
  expect_true(file.exists(pdf))
})

test_that("save_mummichog_exports is a no-op when there is nothing to save", {
  tmp <- withr::local_tempdir()
  expect_length(save_mummichog_exports(NULL, NULL, tmp), 0)
})

test_that("mummichog_report_titles uses an explicit contrast for the subtitle", {
  cfg <- list(modes = list(metabolomics = list(
    organism   = "Coelastrella",
    enrichment = list(mummichog = list(model_json = "/m/cre.json")))))
  t <- mummichog_report_titles(cfg, contrast = "HL_48h_vs_LL_48h")
  expect_identical(t$subtitle, "HL_48h vs LL_48h, all features")
  expect_match(t$title, "Coelastrella")
})

test_that("build_mummichog_report_sections builds one section per usable contrast", {
  cfg <- list(modes = list(metabolomics = list(
    enrichment = list(mummichog = list(p_cutoff = 0.05)))))
  by_contrast <- list(
    "HL_24h_vs_LL_24h" = make_pw(),
    "HL_48h_vs_LL_48h" = make_pw(),
    "empty_contrast"   = NULL             # no result -> dropped
  )
  secs <- build_mummichog_report_sections(by_contrast, cfg)

  expect_named(secs, c("HL_24h_vs_LL_24h", "HL_48h_vs_LL_48h"))
  for (nm in names(secs)) {
    expect_s3_class(secs[[nm]]$plot, "ggplot")
    expect_s3_class(secs[[nm]]$table, "data.frame")
    expect_true(nzchar(secs[[nm]]$title))
    expect_true(nzchar(secs[[nm]]$slug))
  }
  expect_match(secs[["HL_24h_vs_LL_24h"]]$subtitle, "HL_24h vs LL_24h")
})

test_that("build_mummichog_report_sections gives colliding labels distinct slugs", {
  cfg <- list(modes = list(metabolomics = list(
    enrichment = list(mummichog = list(p_cutoff = 0.05)))))
  # two distinct labels that sanitise to the same token must not share a slug,
  # else their standalone exports would overwrite each other
  secs  <- build_mummichog_report_sections(
    list("A-B" = make_pw(), "A_B" = make_pw()), cfg)
  slugs <- vapply(secs, `[[`, character(1), "slug")

  expect_named(secs, c("A-B", "A_B"))                    # display keys keep originals
  expect_length(unique(slugs), 2L)                       # slugs are unique
  expect_true(all(grepl("^[A-Za-z0-9_]+$", slugs)))      # filesystem-safe
})

test_that("build_mummichog_report_sections returns empty list on NULL/empty input", {
  cfg <- list(modes = list(metabolomics = list(
    enrichment = list(mummichog = list()))))
  expect_length(build_mummichog_report_sections(NULL, cfg), 0)
  expect_length(build_mummichog_report_sections(list(), cfg), 0)
  expect_length(build_mummichog_report_sections(list(a = NULL), cfg), 0)
})


# ---- GSEA integration (06g wired into the 06e section/export machinery) ----
#
# These exercise the wiring, not the statistics: that the engine is called with
# the right arguments, that its results travel on the section without changing
# what was already there, and that each way of not running it is distinguishable
# from the others. The engine's own behaviour is covered by test-mummichog-gsea.R.

# The shared fixture builds a faithful v2 result tree plus a matching model and
# DE table; read the pathway tables back the way the pipeline does.
gsea_wiring_fixture <- function(root, method = "limma") {
  f <- build_gsea_fixture(root)
  list(
    pathways = read_mummichog_pathways_by_contrast(f$files),
    files    = f$files,
    model    = f$model,
    de_res   = list(de_tables = list(HL_vs_LL = f$de_table), method = method),
    config   = list(modes = list(metabolomics = list(
      enrichment = list(mummichog = list(enabled = TRUE, p_cutoff = 0.05)))))
  )
}

skip_unless_gsea_deps <- function() {
  for (pkg in c("fgsea", "fastmatch", "data.table", "BiocParallel")) {
    skip_if_not_installed(pkg)
  }
}

test_that("GSEA is attached per contrast without disturbing the ORA fields", {
  skip_unless_gsea_deps()
  fx <- gsea_wiring_fixture(withr::local_tempdir())

  secs <- build_mummichog_report_sections(
    fx$pathways, fx$config, files = fx$files, de_res = fx$de_res,
    model = fx$model)

  s <- secs[["HL_vs_LL"]]
  expect_false(is.null(s$gsea))
  expect_false(is.null(s$gsea_plot))
  expect_s3_class(s$gsea_plot, "ggplot")

  # `plot` is still the ORA bubble, not the GSEA scatter: the Shiny payload
  # carries this field and its meaning must not change.
  expect_s3_class(s$plot, "ggplot")
  bubble <- plot_mummichog_bubble(fx$pathways[["HL_vs_LL"]], title = s$title,
                                  subtitle = s$subtitle, p_cutoff = 0.05)
  expect_identical(s$plot$labels, bubble$labels)
  expect_false(identical(s$plot$labels, s$gsea_plot$labels))

  # And the ORA table is untouched by GSEA running.
  expect_identical(s$table, build_mummichog_pathway_table(fx$pathways[["HL_vs_LL"]]))
})

test_that("the DE method decides the ranking statistic, not the column name", {
  skip_unless_gsea_deps()
  root <- withr::local_tempdir()

  # limma puts a SIGNED moderated t in `statistic`, so GSEA ranks on it.
  lim <- gsea_wiring_fixture(file.path(root, "limma"), method = "limma")
  s_lim <- build_mummichog_report_sections(
    lim$pathways, lim$config, files = lim$files, de_res = lim$de_res,
    model = lim$model)[["HL_vs_LL"]]
  expect_identical(s_lim$gsea$metric, "t_statistic")
  expect_identical(s_lim$gsea$de_method, "limma")

  # wilcoxon puts an UNSIGNED rank-sum W in the SAME column, so ranking on it
  # would order by magnitude with no direction while claiming a t statistic.
  wil <- gsea_wiring_fixture(file.path(root, "wilcox"), method = "wilcoxon")
  s_wil <- build_mummichog_report_sections(
    wil$pathways, wil$config, files = wil$files, de_res = wil$de_res,
    model = wil$model)[["HL_vs_LL"]]
  expect_identical(s_wil$gsea$metric, "logFC")
  expect_match(s_wil$gsea$metric_label, "unsigned rank-sum")

  # The fixtures differ only in de_res$method, so this is the method deciding.
  expect_false(identical(s_lim$gsea$metric, s_wil$gsea$metric))
})

test_that("omitting de_res skips GSEA but keeps ORA and evidence", {
  skip_unless_gsea_deps()
  fx <- gsea_wiring_fixture(withr::local_tempdir())

  # The unwired default: no DE table at all means no trustworthy ranking.
  s <- build_mummichog_report_sections(
    fx$pathways, fx$config, files = fx$files, model = fx$model)[["HL_vs_LL"]]

  expect_null(s$gsea)
  expect_null(s$gsea_plot)
  expect_s3_class(s$plot, "ggplot")
  expect_false(is.null(s$table))
})

test_that("gsea_enabled = FALSE runs the ORA alone", {
  skip_unless_gsea_deps()
  fx <- gsea_wiring_fixture(withr::local_tempdir())
  fx$config$modes$metabolomics$enrichment$mummichog$gsea_enabled <- FALSE

  s <- build_mummichog_report_sections(
    fx$pathways, fx$config, files = fx$files, de_res = fx$de_res,
    model = fx$model)[["HL_vs_LL"]]

  expect_null(s$gsea)
  expect_s3_class(s$plot, "ggplot")

  # Absent key means ON, so configs written before GSEA existed keep working.
  fx$config$modes$metabolomics$enrichment$mummichog$gsea_enabled <- NULL
  on <- build_mummichog_report_sections(
    fx$pathways, fx$config, files = fx$files, de_res = fx$de_res,
    model = fx$model)[["HL_vs_LL"]]
  expect_false(is.null(on$gsea))
})

test_that("a contrast with no matching DE table is named, not silently skipped", {
  skip_unless_gsea_deps()
  fx <- gsea_wiring_fixture(withr::local_tempdir())
  # Two tables, neither named for this contrast: picking one would rank the
  # contrast on another contrast's statistics, so the engine is not called.
  fx$de_res$de_tables <- list(OTHER_A = fx$de_res$de_tables[[1]],
                              OTHER_B = fx$de_res$de_tables[[1]])

  expect_message(
    secs <- build_mummichog_report_sections(
      fx$pathways, fx$config, files = fx$files, de_res = fx$de_res,
      model = fx$model),
    "HL_vs_LL")

  s <- secs[["HL_vs_LL"]]
  expect_null(s$gsea)
  expect_s3_class(s$plot, "ggplot")   # ORA unaffected
  expect_false(is.null(s$table))
})

test_that("save_mummichog_exports writes and then clears the GSEA files", {
  skip_unless_gsea_deps()
  root <- withr::local_tempdir()
  fx   <- gsea_wiring_fixture(root)
  out  <- withr::local_tempdir()

  s <- build_mummichog_report_sections(
    fx$pathways, fx$config, files = fx$files, de_res = fx$de_res,
    model = fx$model)[["HL_vs_LL"]]

  written <- save_mummichog_exports(s$plot, s$table, out, contrast_label = s$slug,
                                    gsea = s$gsea, gsea_plot = s$gsea_plot)
  gsea_tsv <- file.path(out, "mummichog_pinned",
                        "mummichog_gsea_table_HL_vs_LL.tsv")
  expect_true(file.exists(gsea_tsv))
  expect_true(gsea_tsv %in% written)

  # A rerun that comes back without GSEA must not leave the old table behind
  # reading as current -- and must not touch the pathway exports beside it.
  pathway_tsv <- file.path(out, "mummichog_pinned",
                           "mummichog_pathway_table_HL_vs_LL.tsv")
  expect_true(file.exists(pathway_tsv))
  save_mummichog_exports(s$plot, s$table, out, contrast_label = s$slug)
  expect_false(file.exists(gsea_tsv))
  expect_true(file.exists(pathway_tsv))
})

test_that("the GSEA scatter is drawn at the configured p cutoff", {
  skip_unless_gsea_deps()
  fx <- gsea_wiring_fixture(withr::local_tempdir())
  fx$config$modes$metabolomics$enrichment$mummichog$p_cutoff <- 0.2

  s <- build_mummichog_report_sections(
    fx$pathways, fx$config, files = fx$files, de_res = fx$de_res,
    model = fx$model)[["HL_vs_LL"]]

  # The report's figure legend quotes the cutoff from config, so a scatter
  # drawn at the function default would caption a line it did not draw.
  ref <- plot_mummichog_gsea_scatter(s$gsea$table, title = s$title,
                                     subtitle = s$subtitle, p_cutoff = 0.2)
  other <- plot_mummichog_gsea_scatter(s$gsea$table, title = s$title,
                                       subtitle = s$subtitle, p_cutoff = 0.05)
  yint <- function(p) {
    vapply(p$layers, function(l) {
      v <- l$data$yintercept
      if (is.null(v)) NA_real_ else as.numeric(v)[1]
    }, numeric(1))
  }
  expect_equal(yint(s$gsea_plot), yint(ref))
  expect_false(isTRUE(all.equal(yint(s$gsea_plot), yint(other))))
})

test_that("the evidence and GSEA cleanups coexist in one call", {
  out <- withr::local_tempdir()
  dest <- file.path(out, "mummichog_pinned")

  tbl  <- data.frame(Pathway = "P", stringsAsFactors = FALSE)
  evid <- list(pathway_summary = tbl, ec_table = tbl, feature_table = tbl)
  gsea <- list(table = tbl)

  first <- save_mummichog_exports(NULL, tbl, out, contrast_label = "C",
                                  gsea = gsea, evidence = evid)
  ev_tsv   <- file.path(dest, "mummichog_evidence_pathways_C.tsv")
  gs_tsv   <- file.path(dest, "mummichog_gsea_table_C.tsv")
  path_tsv <- file.path(dest, "mummichog_pathway_table_C.tsv")
  expect_true(all(file.exists(ev_tsv, gs_tsv, path_tsv)))
  expect_true(all(c(ev_tsv, gs_tsv) %in% first))

  # A rerun with neither must clear BOTH sets. Adding the GSEA cleanup once
  # replaced the evidence one instead of sitting beside it, which left stale
  # evidence tables on disk reading as current; this pins that they coexist.
  save_mummichog_exports(NULL, tbl, out, contrast_label = "C")
  expect_false(file.exists(ev_tsv))
  expect_false(file.exists(gs_tsv))
  expect_true(file.exists(path_tsv))
})

test_that("a contrast with no section still loses its stale GSEA exports", {
  root <- withr::local_tempdir()
  out  <- withr::local_tempdir()
  dest <- file.path(out, "mummichog_pinned")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  f <- build_gsea_fixture(root)          # one contrast: HL_vs_LL

  # A previous run's GSEA exports for that contrast, left behind because this
  # run produced no section for it (no plottable ORA table, so the section
  # loop never visits it and never clears them by name).
  stale <- file.path(dest, c("mummichog_gsea_table_HL_vs_LL.tsv",
                             "mummichog_gsea_table_HL_vs_LL.csv",
                             "mummichog_gsea_scatter_HL_vs_LL.png",
                             "mummichog_gsea_scatter_HL_vs_LL.pdf"))
  file.create(stale)
  # Exports of other kinds, and of other contrasts, must survive.
  keep <- file.path(dest, c("mummichog_pathway_table_HL_vs_LL.tsv",
                            "mummichog_evidence_pathways_HL_vs_LL.tsv",
                            "mummichog_gsea_table_OTHER.tsv"))
  file.create(keep)

  gone <- clear_orphan_mummichog_gsea_exports(f$files, list(), out)

  expect_setequal(gone, stale)
  expect_false(any(file.exists(stale)))
  expect_true(all(file.exists(keep)))
})

test_that("a contrast with no section also loses its stale Summary exports", {
  root <- withr::local_tempdir()
  out  <- withr::local_tempdir()
  dest <- file.path(out, "mummichog_pinned")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  f <- build_gsea_fixture(root)          # one contrast: HL_vs_LL

  stale <- file.path(dest, c("mummichog_summary_HL_vs_LL.tsv",
                             "mummichog_summary_HL_vs_LL.csv"))
  file.create(stale)
  keep <- file.path(dest, "mummichog_summary_OTHER.tsv")
  file.create(keep)

  gone <- clear_orphan_mummichog_gsea_exports(f$files, list(), out)
  expect_setequal(gone, stale)
  expect_false(any(file.exists(stale)))
  expect_true(file.exists(keep))

  # A live section's Summary is the section's to manage, never an orphan.
  file.create(stale)
  gone <- clear_orphan_mummichog_gsea_exports(
    f$files, list(HL_vs_LL = list(slug = "HL_vs_LL")), out)
  expect_length(gone, 0L)
  expect_true(all(file.exists(stale)))
})

test_that("a live section's GSEA exports are never cleared as orphans", {
  root <- withr::local_tempdir()
  out  <- withr::local_tempdir()
  dest <- file.path(out, "mummichog_pinned")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  f <- build_gsea_fixture(root)
  live <- file.path(dest, "mummichog_gsea_table_HL_vs_LL.tsv")
  file.create(live)

  # The same contrast DOES have a section this run, so its exports are the
  # section's to manage -- deleting them here would remove a current result.
  gone <- clear_orphan_mummichog_gsea_exports(
    f$files, list(HL_vs_LL = list(slug = "HL_vs_LL")), out)

  expect_length(gone, 0L)
  expect_true(file.exists(live))
})

test_that("orphan cleanup follows the de-duplicated slug, not the bare sanitise", {
  out  <- withr::local_tempdir()
  root <- withr::local_tempdir()
  dest <- file.path(out, "mummichog_pinned")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  # Two contrast labels that collapse to ONE token. The stage therefore names
  # its directories A_B and A_B_1, and a run in which both produced sections
  # wrote exports under both slugs.
  mp <- file.path(root, "mummichog_pinned")
  for (d in c("A_B", "A_B_1")) {
    dir.create(file.path(mp, d, "v2", "1.run", "tables"), recursive = TRUE,
               showWarnings = FALSE)
    file.create(file.path(mp, d, "v2", "1.run", "tables",
                          sprintf("mcg_pathwayanalysis_%s.tsv", d)))
  }
  readr::write_tsv(
    data.frame(dir = c("A_B", "A_B_1"), contrast = c("A-B", "A_B"),
               stringsAsFactors = FALSE),
    file.path(mp, "contrasts.tsv"))
  files <- list.files(mp, recursive = TRUE, full.names = TRUE)

  live_f  <- file.path(dest, "mummichog_gsea_table_A_B.tsv")
  stale_f <- file.path(dest, "mummichog_gsea_table_A_B_1.tsv")
  file.create(live_f, stale_f)

  # Only one contrast has a section this run, so its slug is the bare token.
  # Sanitising each label independently would yield c("A_B", "A_B") and the
  # "_1" file would never be named -- the case this function exists to catch.
  gone <- clear_orphan_mummichog_gsea_exports(
    files, list(A_B = list(slug = "A_B")), out)

  expect_true(stale_f %in% gone)
  expect_false(file.exists(stale_f))
  expect_true(file.exists(live_f))       # never the live section's own file
})
