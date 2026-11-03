# tests/testthat/test-mummichog-mixed-polarity.R
#
# ionization_mode: mixed -- one mummichog run per contrast x ion polarity
# (05b + 06i), the report wiring that keys a polarity run back to its contrast
# (06e), and the single-polarity behaviour that must not change. The engine is
# stubbed: these tests check what the stage SUBMITS, not mummichog itself.
# All data are synthetic.


# ---- fixtures ----------------------------------------------------------------

.mixed_pre <- function(ions = c("[M+H]+1", "[M+Na]+1", "[M+2H]+2",
                                "[M-H]-1", "[M-2H]-2")) {
  list(row_data = data.frame(
    check.names     = FALSE,
    feature_id      = c("P1", "P2", "P3", "N1", "N2")[seq_along(ions)],
    "m/z"           = c(181.07, 203.05, 342.12, 179.06, 133.01)[seq_along(ions)],
    RT              = c(1, 2, 3, 4, 5)[seq_along(ions)],
    "Reference Ion" = ions,
    stringsAsFactors = FALSE))
}

.mixed_de <- function(ids = c("P1", "P2", "P3", "N1", "N2")) {
  tab <- function(p) data.frame(feature_id = ids, P.Value = p[seq_along(ids)],
                                logFC = c(1, -1, 0.5, 2, -0.3)[seq_along(ids)],
                                stringsAsFactors = FALSE)
  list(de_tables = list(A_vs_B = tab(c(0.01, 0.2, 0.5, 0.03, 0.6)),
                        C_vs_D = tab(c(0.4, 0.02, 0.04, 0.7, 0.01))))
}

# Runs the real stage with the interpreter check and the engine call stubbed.
# Each stubbed engine call records the ids it was handed and the -m mode.
.mixed_run_stage <- function(pre, de_res, mummi_cfg, calls = new.env()) {
  calls$runs <- list()
  env  <- environment(mod_mummichog_pinned)
  orig <- list(pf  = get(".mmc_preflight_runner", envir = env),
               run = get("run_mummichog_v2", envir = env))
  withr::defer({
    assign(".mmc_preflight_runner", orig$pf, envir = env)
    assign("run_mummichog_v2", orig$run, envir = env)
  })
  assign(".mmc_preflight_runner", function(python) python, envir = env)
  assign("run_mummichog_v2", function(infile, out_dir, project = "mummichog_run", ...) {
    args <- list(...)
    sent <- utils::read.delim(infile, colClasses = "character", check.names = FALSE)
    calls$runs[[project]] <- list(ids = sent$CompoundID, mode = args$mode)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    stub <- file.path(out_dir, "stub.txt")
    writeLines("stub", stub)
    stub
  }, envir = env)

  out <- withr::local_tempdir()
  cfg <- list(modes = list(metabolomics = list(enrichment = list(
    mummichog = utils::modifyList(list(enabled = TRUE), mummi_cfg)))))
  mod_mummichog_pinned(pre, de_res, cfg, out_dir = out, python = "python")
  map <- utils::read.delim(file.path(out, "mummichog_pinned", "contrasts.tsv"),
                           colClasses = "character")
  list(runs = calls$runs, map = map)
}

# The python engine's mixed mode; named explicitly since metaboanalyst is the
# default engine.
.mixed_cfg <- list(engine = "python", ionization_mode = "mixed",
                   polarity_column = "Reference Ion", positive_mode = "positive")


# ---- polarity parsing and labels ----------------------------------------------

test_that("mmc_feature_polarity reads ion labels and plain words", {
  rd <- data.frame(ion = c("[M+H]+1", "[M-H]-1", "[M+2H]2+", "[M-2H]-2",
                           "[M]+", " [M+Cl]- ", "Positive", "neg"),
                   stringsAsFactors = FALSE)
  expect_identical(mmc_feature_polarity(rd, "ion"),
                   c("positive", "negative", "positive", "negative",
                     "positive", "negative", "positive", "negative"))
})

test_that("mmc_feature_polarity stops instead of assuming a polarity", {
  ok <- data.frame(ion = c("[M+H]+1", "[M-H]-1"), stringsAsFactors = FALSE)
  expect_error(mmc_feature_polarity(ok, NULL), "needs polarity_column")
  expect_error(mmc_feature_polarity(ok, "Polarity"), "not a column")
  for (bad in list(NA_character_, "", "[M+H]", "M+H", "unknown")) {
    rd <- data.frame(ion = c("[M+H]+1", bad), stringsAsFactors = FALSE)
    expect_error(mmc_feature_polarity(rd, "ion"), "no readable polarity")
  }
})

test_that("polarity labels round-trip; a plain contrast keeps polarity NA", {
  lab <- mmc_polarity_label("A_vs_B", c("positive", "negative"))
  expect_identical(lab, c("A_vs_B, positive ions", "A_vs_B, negative ions"))
  back <- mmc_split_polarity_label(c(lab, "A_vs_B"))
  expect_identical(back$contrast, c("A_vs_B", "A_vs_B", "A_vs_B"))
  expect_identical(back$polarity, c("positive", "negative", NA))
})


# ---- the stage: mixed mode ----------------------------------------------------

test_that("mixed mode runs each contrast once per polarity, each on its own features", {
  res <- .mixed_run_stage(.mixed_pre(), .mixed_de(), .mixed_cfg)

  expect_named(res$runs, c("A_vs_B_positive_ions", "A_vs_B_negative_ions",
                           "C_vs_D_positive_ions", "C_vs_D_negative_ions"))
  for (ct in c("A_vs_B", "C_vs_D")) {
    pos <- res$runs[[paste0(ct, "_positive_ions")]]
    neg <- res$runs[[paste0(ct, "_negative_ions")]]
    # every eligible feature of the polarity is its background, and none other
    expect_setequal(pos$ids, c("P1", "P2", "P3"))
    expect_setequal(neg$ids, c("N1", "N2"))
    expect_identical(pos$mode, "positive")
    expect_identical(neg$mode, "negative")
  }
  expect_identical(res$map$contrast,
                   c("A_vs_B, positive ions", "A_vs_B, negative ions",
                     "C_vs_D, positive ions", "C_vs_D, negative ions"))
  expect_identical(res$map$dir, names(res$runs))
})

test_that("positive_mode selects the mode of the positive-ion runs", {
  res <- .mixed_run_stage(.mixed_pre(), .mixed_de(),
                          utils::modifyList(.mixed_cfg, list(positive_mode = "pos_default")))
  expect_identical(res$runs[["A_vs_B_positive_ions"]]$mode, "pos_default")
  expect_identical(res$runs[["A_vs_B_negative_ions"]]$mode, "negative")
})

test_that("mixed mode stops before any engine call when polarity is not reliable", {
  calls <- new.env()
  no_col <- utils::modifyList(.mixed_cfg, list(polarity_column = NULL))
  expect_error(.mixed_run_stage(.mixed_pre(), .mixed_de(), no_col, calls),
               "needs polarity_column")
  expect_length(calls$runs, 0)

  no_mode <- utils::modifyList(.mixed_cfg, list(positive_mode = NULL))
  expect_error(.mixed_run_stage(.mixed_pre(), .mixed_de(), no_mode, calls),
               "needs positive_mode")
  expect_length(calls$runs, 0)

  unreadable <- .mixed_pre(c("[M+H]+1", "[M+Na]+1", "[M+2H]", "[M-H]-1", NA))
  expect_error(.mixed_run_stage(unreadable, .mixed_de(), .mixed_cfg, calls),
               "2 of 5 feature\\(s\\) have no readable polarity")
  expect_length(calls$runs, 0)
})

test_that("a polarity with no features is skipped with a message; the other runs", {
  neg_only <- .mixed_pre(c("[M-H]-1", "[M-2H]-2"))
  expect_message(
    res <- .mixed_run_stage(neg_only, .mixed_de(c("P1", "P2")), .mixed_cfg),
    "no positive-ion features")
  expect_named(res$runs, c("A_vs_B_negative_ions", "C_vs_D_negative_ions"))
})


# ---- the stage: single polarity is unchanged -----------------------------------

test_that("single-polarity modes keep one run per contrast over all features", {
  for (mode in c("negative", "positive", "pos_default")) {
    res <- .mixed_run_stage(.mixed_pre(), .mixed_de(), list(ionization_mode = mode))
    expect_named(res$runs, c("A_vs_B", "C_vs_D"))
    expect_setequal(res$runs[["A_vs_B"]]$ids, c("P1", "P2", "P3", "N1", "N2"))
    expect_identical(res$runs[["A_vs_B"]]$mode, mode)
    expect_identical(res$map$contrast, c("A_vs_B", "C_vs_D"))
    expect_identical(res$map$dir, c("A_vs_B", "C_vs_D"))
  }
  # unset -> pos_default, and no polarity column is needed or read
  pre <- .mixed_pre()
  pre$row_data[["Reference Ion"]] <- NULL
  res <- .mixed_run_stage(pre, .mixed_de(), list())
  expect_identical(res$runs[["A_vs_B"]]$mode, "pos_default")
})


# ---- report wiring and GSEA isolation -------------------------------------------

test_that("a polarity section keys DE and identified support back to its contrast", {
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root, contrast_dir = "HL_vs_LL_positive_ions")
  label <- mmc_polarity_label("HL_vs_LL", "positive")
  readr::write_tsv(data.frame(dir = "HL_vs_LL_positive_ions", contrast = label),
                   file.path(root, "mummichog_pinned", "contrasts.tsv"))

  gmt <- file.path(root, "library.gmt")
  writeLines(c("p1\tPW_pos\tC00001", "p2\tPW_mid\tC00009", "p3\tPW_neg\tC00017"), gmt)
  cfg <- list(modes = list(metabolomics = list(enrichment = list(
    gmt_file  = gmt,
    mummichog = list(engine = "python", p_cutoff = 0.05, ionization_mode = "mixed",
                     polarity_column = "Reference Ion", positive_mode = "positive")))))
  # Two DE tables, so a lookup by the run label (not the contrast) cannot fall
  # back to "the only table".
  de_res <- list(method = "limma",
                 de_tables = list(HL_vs_LL = f$de_table, Other = f$de_table[0, ]))
  enr <- list(qea_by_contrast = list(HL_vs_LL = list(res = list(
    table = data.frame(pathway = "p1 - PW_pos", FDR = 0.01)))))

  secs <- build_mummichog_report_sections(
    stats::setNames(list(read_mummichog_pathways(f$files)), label), cfg,
    files = f$files, de_res = de_res, model = f$model, enrichment_res = enr)

  s <- secs[[label]]
  expect_identical(s$polarity, "positive")
  expect_match(s$subtitle, "HL vs LL, positive-ion features only", fixed = TRUE)
  row <- s$summary[s$summary$Pathway == "PW_pos", , drop = FALSE]
  expect_identical(row[["Identified-metabolite support"]], "QEA")
  if (requireNamespace("fgsea", quietly = TRUE)) expect_false(is.null(s$gsea))

  # single-polarity labels are untouched
  expect_match(mummichog_report_titles(cfg, contrast = "HL_vs_LL")$subtitle,
               "HL vs LL, all features", fixed = TRUE)
})

test_that("a run label is split into contrast and polarity only in mixed mode", {
  mode_cfg <- function(m) list(modes = list(metabolomics = list(enrichment = list(
    mummichog = list(engine = "python", ionization_mode = m)))))
  lab <- "X, positive ions"

  mixed <- mmc_run_label_parts(lab, mode_cfg("mixed"))
  expect_identical(mixed$contrast, "X")
  expect_identical(mixed$polarity, "positive")

  for (m in list("negative", "positive", "pos_default", NULL)) {
    single <- mmc_run_label_parts(lab, mode_cfg(m))
    expect_identical(single$contrast, lab)
    expect_identical(single$polarity, NA_character_)
  }
})

test_that("a single-polarity contrast named like a polarity label stays one contrast", {
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root, contrast_dir = "HL_vs_LL_positive_ions")
  label <- "HL_vs_LL, positive ions"          # a user's own contrast name
  readr::write_tsv(data.frame(dir = "HL_vs_LL_positive_ions", contrast = label),
                   file.path(root, "mummichog_pinned", "contrasts.tsv"))

  gmt <- file.path(root, "library.gmt")
  writeLines(c("p1\tPW_pos\tC00001", "p2\tPW_mid\tC00009", "p3\tPW_neg\tC00017"), gmt)
  cfg <- list(modes = list(metabolomics = list(enrichment = list(
    gmt_file = gmt, mummichog = list(p_cutoff = 0.05, ionization_mode = "negative")))))
  # Two DE tables: one under the full name, a decoy under the stripped name.
  de_res <- list(method = "limma",
                 de_tables = stats::setNames(list(f$de_table, f$de_table[0, ]),
                                             c(label, "HL_vs_LL")))
  enr <- list(qea_by_contrast = stats::setNames(list(list(res = list(
    table = data.frame(pathway = "p1 - PW_pos", FDR = 0.01)))), label))

  secs <- build_mummichog_report_sections(
    stats::setNames(list(read_mummichog_pathways(f$files)), label), cfg,
    files = f$files, de_res = de_res, model = f$model, enrichment_res = enr)

  s <- secs[[label]]
  expect_identical(s$polarity, NA_character_)
  expect_match(s$subtitle, "HL vs LL, positive ions, all features", fixed = TRUE)
  row <- s$summary[s$summary$Pathway == "PW_pos", , drop = FALSE]
  expect_identical(row[["Identified-metabolite support"]], "QEA")
  if (requireNamespace("fgsea", quietly = TRUE)) expect_false(is.null(s$gsea))
})

test_that("opposite-polarity features in the DE table cannot change mummichog GSEA", {
  skip_if_not_installed("fgsea")
  root <- withr::local_tempdir()
  f    <- build_gsea_fixture(root)

  # The report passes the contrast's full DE table to each polarity's GSEA.
  # Features of the other polarity are absent from this run's EmpiricalCompounds;
  # give them extreme statistics so any leak would move the result.
  n <- nrow(f$de_table)
  other <- data.frame(feature_id = sprintf("neg_%02d", seq_len(n)),
                      logFC = rep(c(25, -25), length.out = n),
                      statistic = rep(c(50, -50), length.out = n),
                      P.Value = 1e-8, stringsAsFactors = FALSE)

  own  <- run_mummichog_gsea(f$files, f$de_table, f$model, de_method = "limma")
  full <- run_mummichog_gsea(f$files, rbind(f$de_table, other), f$model,
                             de_method = "limma")
  expect_false(is.null(own))
  expect_identical(full$table, own$table)
})
