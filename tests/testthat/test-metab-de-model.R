# tests/testthat/test-metab-de-model.R
#
# The limma model metabolomics DE hands downstream (de_res$de_model, exported
# unchanged as the Shiny payload's de_model):
#   - one fit covers every contrast, with coefficients named by Contrast_name;
#   - the per-contrast tables are read from that same fit, and agree with
#     fitting each contrast on its own (the previous implementation);
#   - the Shiny app's own topTable() call works for every contrast and returns
#     the statistics the tables hold;
#   - QC/blank samples are not part of the fit.
# All data are synthetic and seeded.

STAT_COLS <- c("logFC", "AveExpr", "P.Value", "adj.P.Val")

.model_contrasts <- data.frame(
  Factor        = "grp",
  Numerator     = c("B", "C", "C"),
  Denominator   = c("A", "A", "B"),
  Contrast_name = c("B_vs_A", "C_over_A", "late.vs.mid"),
  stringsAsFactors = FALSE
)

# 12 biological samples (A, B, C x 4) and, optionally, 3 pooled QC samples.
# Features m1-m5 differ between groups; a few values are missing.
.model_fixture <- function(with_qc = TRUE, qc_shift = 0) {
  groups <- c(rep(c("A", "B", "C"), each = 4), if (with_qc) rep("QC", 3))
  ids    <- c(sprintf("S%02d", 1:12), if (with_qc) sprintf("QC%d", 1:3))
  n_feat <- 30
  mat <- withr::with_seed(433, matrix(stats::rnorm(n_feat * length(ids), 12, 0.4),
                                      nrow = n_feat))
  dimnames(mat) <- list(sprintf("m%d", seq_len(n_feat)), ids)
  mat[1:5, groups == "B"] <- mat[1:5, groups == "B"] + 1.5
  mat[1:5, groups == "C"] <- mat[1:5, groups == "C"] + 3
  mat[cbind(c(7, 12, 20), c(2, 6, 11))] <- NA
  if (with_qc) mat[, groups == "QC"] <- mat[, groups == "QC"] + qc_shift

  meta <- data.frame(sample_id = ids, grp = groups, stringsAsFactors = FALSE)
  # Metadata deliberately in a different order from the matrix columns.
  meta <- meta[rev(seq_len(nrow(meta))), , drop = FALSE]
  list(expr_unscaled = mat, expr_work = mat, expr_filt = 2^mat, expr_raw = 2^mat,
       meta = meta)
}

.model_config <- function(method = "limma") {
  list(modes = list(metabolomics = list(
    preprocessing = list(chosen_norm = "pqn", transform = "log2", scaling = "none"),
    effects = list(samples = "sample_id", color = "grp"),
    de      = list(method = method, condition_column = "grp",
                   p_cutoff = 0.05, linear_fc_cutoff = 1.5),
    excel   = list(group_cv = FALSE)
  )))
}

.model_de <- function(pre = .model_fixture(), method = "limma",
                      contrasts = .model_contrasts) {
  suppressMessages(run_metabolomics_de(pre, .model_config(method), contrasts))
}

# The call the Shiny app makes for its volcano and MA plots.
.shiny_top_table <- function(fit, contrast) {
  limma::topTable(fit, coef = contrast, number = Inf, sort.by = "none")
}

.bio_matrix <- function(pre) {
  pre$expr_unscaled[, !grepl("^QC", colnames(pre$expr_unscaled)), drop = FALSE]
}


# --- one model, every contrast -------------------------------------------------

test_that("de_model has one coefficient per contrast, named by Contrast_name", {
  res <- .model_de()

  # The class test the Shiny app applies before calling topTable().
  expect_true(inherits(res$de_model, "MArrayLM"))
  expect_identical(colnames(res$de_model$coefficients),
                   .model_contrasts$Contrast_name)
  expect_identical(names(res$de_tables), .model_contrasts$Contrast_name)
  expect_identical(rownames(res$de_model$coefficients),
                   rownames(.model_fixture()$expr_unscaled))
})

test_that("a single contrast is still a named coefficient", {
  res <- .model_de(contrasts = .model_contrasts[2, , drop = FALSE])

  expect_identical(colnames(res$de_model$coefficients), "C_over_A")
  expect_identical(names(res$de_tables), "C_over_A")
})

test_that("contrasts without a Contrast_name are named <numerator>_vs_<denominator>", {
  res <- .model_de(contrasts = .model_contrasts[, c("Numerator", "Denominator")])

  expected <- c("B_vs_A", "C_vs_A", "C_vs_B")
  expect_identical(colnames(res$de_model$coefficients), expected)
  expect_identical(names(res$de_tables), expected)
})


# --- the Shiny app's topTable() call -------------------------------------------

test_that("the Shiny topTable() call works for every contrast and returns the tables' statistics", {
  res <- .model_de()

  for (ctr in .model_contrasts$Contrast_name) {
    tt <- .shiny_top_table(res$de_model, ctr)
    expect_true(all(STAT_COLS %in% colnames(tt)), info = ctr)

    tab <- res$de_tables[[ctr]]
    expect_false(anyDuplicated(tab$feature_id) > 0, info = ctr)
    expect_setequal(rownames(tt), tab$feature_id)

    # Compare by feature id, in an order that is neither table's own.
    ids <- rev(tab$feature_id)
    for (col in STAT_COLS) {
      expect_identical(tt[ids, col], tab[[col]][match(ids, tab$feature_id)],
                       info = paste(ctr, col))
    }
  }
})

test_that("adjusted p-values are Benjamini-Hochberg within each contrast, over every tested feature", {
  res <- .model_de()

  for (ctr in names(res$de_tables)) {
    tab <- res$de_tables[[ctr]]
    expect_equal(tab$adj.P.Val, stats::p.adjust(tab$P.Value, method = "BH"),
                 tolerance = 1e-12, info = ctr)
  }
})


# --- statistical equivalence with one fit per contrast --------------------------

test_that("the joint fit reproduces fitting each contrast on its own", {
  pre <- .model_fixture()
  mat <- .bio_matrix(pre)
  # A feature never observed in one group: its contrasts against that group
  # cannot be estimated, the remaining one can.
  mat["m9", 1:4] <- NA
  condition <- factor(rep(c("A", "B", "C"), each = 4))
  formulas  <- paste(.model_contrasts$Numerator, .model_contrasts$Denominator,
                     sep = " - ")

  joint <- fit_metab_limma_contrasts(mat, condition, formulas,
                                     .model_contrasts$Contrast_name)

  for (i in seq_along(formulas)) {
    ctr <- .model_contrasts$Contrast_name[i]
    alone <- de_limma(mat, condition, formulas[i])
    tab   <- joint$de_tables[[ctr]]

    expect_setequal(tab$feature_id, alone$feature_id)
    idx <- match(alone$feature_id, tab$feature_id)
    for (col in c(STAT_COLS, "statistic")) {
      expect_identical(is.na(tab[[col]][idx]), is.na(alone[[col]]),
                       info = paste(ctr, col, "NA pattern"))
      expect_equal(tab[[col]][idx], alone[[col]], tolerance = 1e-10,
                   info = paste(ctr, col))
    }
  }
})

test_that("run_metabolomics_de() tables match one fit per contrast on the biological samples", {
  pre <- .model_fixture()
  res <- .model_de(pre)
  mat <- .bio_matrix(pre)
  condition <- factor(rep(c("A", "B", "C"), each = 4))

  for (i in seq_len(nrow(.model_contrasts))) {
    ctr   <- .model_contrasts$Contrast_name[i]
    alone <- de_limma(mat, condition,
                      paste(.model_contrasts$Numerator[i],
                            .model_contrasts$Denominator[i], sep = " - "))
    tab <- res$de_tables[[ctr]]
    idx <- match(alone$feature_id, tab$feature_id)
    for (col in STAT_COLS) {
      expect_equal(tab[[col]][idx], alone[[col]], tolerance = 1e-10,
                   info = paste(ctr, col))
    }
  }
})


# --- output structure ------------------------------------------------------------

test_that("the DE result keeps its structure", {
  res <- .model_de()

  expect_identical(res$method, "limma")
  for (tab in res$de_tables) {
    expect_identical(colnames(tab), c("feature_id", "logFC", "AveExpr",
                                      "statistic", "P.Value", "adj.P.Val"))
    expect_identical(tab$feature_id, rownames(.model_fixture()$expr_unscaled))
  }
  per_contrast <- as.vector(outer(c("linearFC.", "AveExpr.", "pvalue.", "padj.", "pass."),
                                  .model_contrasts$Contrast_name, paste0))
  expect_true(all(c("feature_id", per_contrast, "pass_any_contrast") %in%
                    colnames(res$summary_df)))
  expect_identical(res$padj_provenance$contrast, .model_contrasts$Contrast_name)
  expect_true(all(res$padj_provenance$padj_source == "computed_bh"))
})

test_that("the t-test and Wilcoxon methods carry no limma model", {
  for (method in c("t_test", "t_test_equal", "wilcoxon")) {
    res <- suppressWarnings(.model_de(method = method))
    expect_null(res$de_model, info = method)
    expect_identical(names(res$de_tables), .model_contrasts$Contrast_name,
                     info = method)
  }
})


# --- sample filtering --------------------------------------------------------------

test_that("QC samples are not part of the fit", {
  with_qc    <- .model_de(.model_fixture(with_qc = TRUE))
  shifted_qc <- .model_de(.model_fixture(with_qc = TRUE, qc_shift = 5))
  without_qc <- .model_de(.model_fixture(with_qc = FALSE))

  for (ctr in .model_contrasts$Contrast_name) {
    expect_identical(shifted_qc$de_tables[[ctr]], with_qc$de_tables[[ctr]], info = ctr)
    expect_identical(without_qc$de_tables[[ctr]], with_qc$de_tables[[ctr]], info = ctr)
    expect_identical(.shiny_top_table(without_qc$de_model, ctr),
                     .shiny_top_table(with_qc$de_model, ctr), info = ctr)
  }

  # limma's AveExpr is the row mean of the matrix it fitted.
  pre <- .model_fixture(with_qc = TRUE)
  tab <- with_qc$de_tables[[1]]
  expect_equal(tab$AveExpr,
               unname(rowMeans(.bio_matrix(pre)[tab$feature_id, ], na.rm = TRUE)))
})


# --- guards ------------------------------------------------------------------------

test_that("a repeated Contrast_name is refused", {
  dup <- .model_contrasts
  dup$Contrast_name[3] <- dup$Contrast_name[1]

  expect_error(.model_de(contrasts = dup), "used more than once")
})

test_that("metadata that does not cover every matrix column is refused", {
  pre <- .model_fixture(with_qc = FALSE)
  pre$meta <- pre$meta[pre$meta$sample_id != "S03", , drop = FALSE]

  expect_error(.model_de(pre), "does not line up")
})

test_that("fit_metab_limma_contrasts() refuses names it cannot use as coefficients", {
  mat <- .bio_matrix(.model_fixture())
  condition <- factor(rep(c("A", "B", "C"), each = 4))

  expect_error(fit_metab_limma_contrasts(mat, condition, c("B - A", "C - A"), "only_one"),
               "one name per contrast")
  expect_error(fit_metab_limma_contrasts(mat, condition, c("B - A", "C - A"), c("x", "x")),
               "unique")
  expect_error(fit_metab_limma_contrasts(mat, condition[-1], "B - A", "x"),
               "non-missing")
})


# --- the Shiny payload ---------------------------------------------------------------

test_that("the payload's de_model answers topTable() for every payload contrast", {
  pre <- .model_fixture()
  pre$row_data <- data.frame(feature_id = rownames(pre$expr_unscaled),
                             row.names = rownames(pre$expr_unscaled),
                             stringsAsFactors = FALSE)
  res <- .model_de(pre)

  payload <- suppressMessages(suppressWarnings(
    build_shiny_payload_metabolomics(
      pre = pre, de_res = res, inputs = list(contrasts = .model_contrasts),
      config = .model_config(), include_legacy = FALSE
    )
  ))

  expect_identical(colnames(payload$de_model$coefficients),
                   as.character(payload$contrasts$Contrast_name))
  for (ctr in payload$contrasts$Contrast_name) {
    tt <- .shiny_top_table(payload$de_model, ctr)
    expect_identical(tt, .shiny_top_table(res$de_model, ctr), info = ctr)
    expect_setequal(rownames(tt), payload$de_stats$feature_id)
  }
})
