# tests/testthat/test-mummichog-across.R
#
# The across-contrast view (06u): a lookup of each contrast's existing results,
# aligned on the model's pathway id. Expected values are worked out by hand
# from ma_two_contrast_sections() (helper-mummichog-ma.R). All data SYNTHETIC.

across_row <- function(wide, pathway) {
  r <- wide[wide$Pathway == pathway, , drop = FALSE]
  expect_identical(nrow(r), 1L, label = pathway)
  r
}

test_that("each contrast's evidence is read from its own overview", {
  ac <- build_mummichog_across_contrasts(ma_two_contrast_sections()$sections)
  expect_identical(ac$contrasts, c("HL_vs_LL", "LL_vs_DD"))
  w <- ac$wide
  g <- across_row(w, "Glycolysis")
  expect_identical(g[["HL_vs_LL"]], "ORA Holm · GSEA nominal")
  expect_identical(g[["LL_vs_DD"]], "GSEA FDR")
  expect_identical(g[["Contrasts with adjusted evidence"]], 2L)
  expect_identical(g[["Contrasts with nominal evidence"]], 2L)
  p <- across_row(w, "Pyruvate metabolism")
  expect_identical(p[["HL_vs_LL"]], "—")
  expect_identical(p[["LL_vs_DD"]], "ORA Holm")
  expect_identical(p[["Contrasts with adjusted evidence"]], 1L)
  tca <- across_row(w, "TCA cycle")
  expect_identical(tca[["HL_vs_LL"]], "GSEA FDR")
  expect_identical(tca[["LL_vs_DD"]], "No result")
  expect_identical(tca[["Contrasts with a result"]], 1L)
})

test_that("ORA Holm and GSEA BH FDR are carried separately, per contrast", {
  tc <- ma_two_contrast_sections()
  w <- build_mummichog_across_contrasts(tc$sections)$wide
  g <- across_row(w, "Glycolysis")
  expect_equal(g[["HL_vs_LL: ORA raw p"]], 0.01)
  expect_equal(g[["HL_vs_LL: ORA Holm-adjusted p"]], 0.02)
  expect_equal(g[["HL_vs_LL: GSEA raw p"]], 0.03)
  expect_equal(g[["HL_vs_LL: GSEA BH FDR"]], 0.06)
  expect_equal(g[["LL_vs_DD: ORA raw p"]], 0.3)
  expect_equal(g[["LL_vs_DD: GSEA BH FDR"]], 0.04)
  tca <- across_row(w, "TCA cycle")
  expect_true(is.na(tca[["LL_vs_DD: GSEA raw p"]]))
  # Every value is the overview's own.
  for (ct in c("HL_vs_LL", "LL_vs_DD")) {
    ov <- tc$sections[[ct]]$overview
    at <- match(ov$Pathway, w$Pathway)
    expect_equal(w[[paste0(ct, ": ORA raw p")]][at], ov[["ORA raw p"]])
    expect_equal(w[[paste0(ct, ": GSEA BH FDR")]][at], ov[["GSEA BH FDR"]])
  }
})

test_that("rows are ordered by adjusted, then nominal recurrence, then smallest p", {
  w <- build_mummichog_across_contrasts(ma_two_contrast_sections()$sections)$wide
  # Glycolysis (2 adjusted), then Pyruvate (best p 0.004) before TCA (0.005).
  expect_identical(w$Pathway, c("Glycolysis", "Pyruvate metabolism", "TCA cycle"))
})

test_that("pathways align on the model id, and on the name only without one", {
  mk <- function(id, p) {
    list(p_cutoff = 0.05, overview = data.frame(
      check.names = FALSE, stringsAsFactors = FALSE,
      Pathway = "Same pathway", "Pathway ID" = id,
      "ORA raw p" = p, "ORA Holm-adjusted p" = NA_real_,
      "GSEA raw p" = NA_real_, "GSEA BH FDR" = NA_real_))
  }
  # Same name, different ids: two pathways, never merged by name.
  w <- build_mummichog_across_contrasts(list(A = mk("x00010", 0.01),
                                             B = mk("x00020", 0.01)))$wide
  expect_identical(nrow(w), 2L)
  # No id anywhere: aligned by name.
  w2 <- build_mummichog_across_contrasts(list(A = mk(NA_character_, 0.01),
                                              B = mk(NA_character_, 0.02)))$wide
  expect_identical(nrow(w2), 1L)
  expect_identical(w2[["Contrasts with nominal evidence"]], 2L)
})

test_that("per-polarity runs of one contrast count as one contrast", {
  run <- function(de_contrast, p) {
    list(p_cutoff = 0.05, de_contrast = de_contrast, overview = data.frame(
      check.names = FALSE, stringsAsFactors = FALSE,
      Pathway = "P", "Pathway ID" = "x00010",
      "ORA raw p" = p, "ORA Holm-adjusted p" = NA_real_,
      "GSEA raw p" = NA_real_, "GSEA BH FDR" = NA_real_))
  }
  # python mixed mode: one section per polarity, both of contrast A_vs_B.
  secs <- list("A_vs_B, positive ions" = run("A_vs_B", 0.01),
               "A_vs_B, negative ions" = run("A_vs_B", 0.02))
  ac <- build_mummichog_across_contrasts(secs)
  expect_identical(ac$de_contrasts, "A_vs_B")
  expect_identical(ac$wide[["Contrasts with nominal evidence"]], 1L)
  # Each polarity keeps its own evidence cell.
  expect_identical(ac$wide[["A_vs_B, positive ions"]], "ORA nominal")
  expect_identical(ac$wide[["A_vs_B, negative ions"]], "ORA nominal")
  # The real sections carry their DE contrast.
  tc <- ma_two_contrast_sections()
  expect_identical(tc$sections[["HL_vs_LL"]]$de_contrast, "HL_vs_LL")
})

test_that("nothing directional or combined is claimed", {
  ac <- build_mummichog_across_contrasts(ma_two_contrast_sections()$sections)
  all_names <- c(names(ac$wide), names(ac$long))
  expect_false(any(grepl("NES|direction|up|down|combined|pooled", all_names,
                         ignore.case = TRUE)))
})

test_that("the report lists nominal-in-two OR adjusted-in-one; the workbook lists all", {
  w <- data.frame(
    check.names = FALSE,
    Pathway = c("nominal in 2", "adjusted in 1", "nominal in 1 only", "nothing", "NA counts"),
    "Contrasts with nominal evidence"  = c(2L, 1L, 1L, 0L, NA),
    "Contrasts with adjusted evidence" = c(0L, 1L, 0L, 0L, NA))
  expect_identical(mmc_across_report_rows(w), c(TRUE, TRUE, FALSE, FALSE, FALSE))

  # On the two-contrast fixture every pathway qualifies: Glycolysis (nominal in
  # both), Pyruvate (ORA Holm in LL_vs_DD), TCA (GSEA FDR in HL_vs_LL).
  ac <- build_mummichog_across_contrasts(ma_two_contrast_sections()$sections)
  expect_identical(ac$wide$Pathway[mmc_across_report_rows(ac$wide)],
                   c("Glycolysis", "Pyruvate metabolism", "TCA cycle"))
})

test_that("sections without an overview give no across-contrast view", {
  expect_null(build_mummichog_across_contrasts(list()))
  expect_null(build_mummichog_across_contrasts(list(A = list(table = data.frame()))))
})
