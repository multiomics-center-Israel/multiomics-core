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
