# tests/testthat/test-metab-de-percontrast.R
#
# de.method = "limma_percontrast" for metabolomics:
#   - each contrast is normalized and fitted on the samples of its two groups
#     only, and equals doing that by hand on the subset;
#   - samples of other groups and QC pools change nothing in a contrast;
#   - the per-contrast matrices are returned for export;
#   - unsupported normalizations, drift correction and under-sized groups stop
#     with a message instead of running something else.
# All data are synthetic and seeded.

STAT_COLS_PC <- c("feature_id", "logFC", "AveExpr", "statistic", "P.Value", "adj.P.Val")

.pc_contrasts <- data.frame(
  Factor        = "grp",
  Numerator     = c("B", "C"),
  Denominator   = c("A", "B"),
  Contrast_name = c("B_vs_A", "C_vs_B"),
  stringsAsFactors = FALSE
)

# 12 biological samples (A, B, C x 4) plus 3 pooled QCs, linear scale, with
# per-sample loading differences so normalization matters.
.pc_fixture <- function(c_shift = 0, qc_shift = 0) {
  groups <- c(rep(c("A", "B", "C"), each = 4), rep("QC", 3))
  ids    <- c(sprintf("S%02d", 1:12), sprintf("QC%d", 1:3))
  n_feat <- 25
  log_mat <- withr::with_seed(501, matrix(stats::rnorm(n_feat * length(ids), 10, 0.5),
                                          nrow = n_feat))
  dimnames(log_mat) <- list(sprintf("m%d", seq_len(n_feat)), ids)
  log_mat[1:4, groups == "B"] <- log_mat[1:4, groups == "B"] + 1.5
  log_mat[1:4, groups == "C"] <- log_mat[1:4, groups == "C"] + 2.5 + c_shift
  log_mat[, groups == "QC"]   <- log_mat[, groups == "QC"] + qc_shift
  loading <- withr::with_seed(502, stats::runif(length(ids), -0.6, 0.6))
  lin <- 2^sweep(log_mat, 2, loading, FUN = "+")
  lin[cbind(c(6, 15), c(3, 9))] <- NA

  meta <- data.frame(sample_id = ids, grp = groups, stringsAsFactors = FALSE)
  # A whole-run log matrix stands in for expr_unscaled; per-contrast DE must
  # not read it.
  whole <- log2(lin + 1)
  list(expr_unscaled = whole, expr_work = whole, expr_filt = lin, expr_raw = lin,
       meta = meta)
}

.pc_config <- function(chosen_norm = "pqn") {
  list(modes = list(metabolomics = list(
    preprocessing = list(chosen_norm = chosen_norm, transform = "log2",
                         scaling = "none", pseudocount = 1),
    effects = list(samples = "sample_id", color = "grp"),
    de      = list(method = "limma_percontrast", condition_column = "grp",
                   p_cutoff = 0.05, linear_fc_cutoff = 1.5)
  )))
}

.pc_de <- function(pre = .pc_fixture(), chosen_norm = "pqn",
                   contrasts = .pc_contrasts) {
  suppressMessages(run_metabolomics_de(pre, .pc_config(chosen_norm), contrasts))
}

# The same contrast done by hand: subset, normalize, log2(x + 1), limma.
.pc_by_hand <- function(lin, groups, num, den, normalize) {
  keep <- groups %in% c(num, den)
  sub  <- transform_metab(normalize(lin[, keep, drop = FALSE]), "log2", 1)
  cond <- factor(groups[keep], levels = c(den, num))
  fit_metab_limma_contrasts(sub, cond, paste(num, "-", den),
                            paste0(num, "_vs_", den))$de_tables[[1]]
}


test_that("each contrast equals normalizing and fitting its two groups by hand", {
  pre <- .pc_fixture()
  res <- .pc_de(pre)
  groups <- pre$meta$grp

  expect_identical(res$method, "limma_percontrast")
  expect_identical(names(res$de_tables), .pc_contrasts$Contrast_name)
  expect_null(res$de_model)

  for (i in seq_len(nrow(.pc_contrasts))) {
    num <- .pc_contrasts$Numerator[i]
    den <- .pc_contrasts$Denominator[i]
    expected <- .pc_by_hand(pre$expr_filt, groups, num, den,
                            function(m) norm_pqn(m))
    got <- res$de_tables[[.pc_contrasts$Contrast_name[i]]]
    expect_equal(got[, STAT_COLS_PC], expected[, STAT_COLS_PC])
  }
})

test_that("tss and none normalize the subset the same way", {
  pre <- .pc_fixture()
  groups <- pre$meta$grp
  for (nm in c("tss", "none")) {
    res <- .pc_de(pre, chosen_norm = nm)
    normalize <- if (nm == "tss") norm_total_sum else function(m) m
    expected <- .pc_by_hand(pre$expr_filt, groups, "B", "A", normalize)
    expect_equal(res$de_tables$B_vs_A[, STAT_COLS_PC], expected[, STAT_COLS_PC],
                 info = nm)
    expect_identical(res$percontrast_norm$norm_method, nm)
  }
})

test_that("other groups and QC pools do not change a contrast", {
  base    <- .pc_de(.pc_fixture())
  shifted <- .pc_de(.pc_fixture(c_shift = 3, qc_shift = 4))
  expect_equal(shifted$de_tables$B_vs_A, base$de_tables$B_vs_A)
  # The C contrast does move, so the shift was real.
  expect_false(isTRUE(all.equal(shifted$de_tables$C_vs_B$logFC,
                                base$de_tables$C_vs_B$logFC)))
})

test_that("the matrix each contrast was tested on is returned, two groups only", {
  pre <- .pc_fixture()
  res <- .pc_de(pre)
  pc  <- res$percontrast_norm
  expect_identical(names(pc$matrices), .pc_contrasts$Contrast_name)
  expect_identical(colnames(pc$matrices$B_vs_A), sprintf("S%02d", 1:8))
  expect_identical(colnames(pc$matrices$C_vs_B), sprintf("S%02d", 5:12))
  expect_identical(rownames(pc$matrices$B_vs_A), rownames(pre$expr_filt))
  expect_identical(pc$transform, "log2")
  expect_equal(pc$pseudocount, 1)
})

test_that("normalizations fitted another way are refused, not approximated", {
  for (nm in c("median", "eigenms", "bio_factor")) {
    expect_error(.pc_de(chosen_norm = nm), "supports only", info = nm)
  }
})

test_that("drift correction is refused", {
  pre <- .pc_fixture()
  pre$info <- list(drift = list(drift_applied = TRUE))
  expect_error(.pc_de(pre), "drift correction")
})

test_that("a group with fewer than 2 samples stops and names the contrast", {
  ctr <- data.frame(Factor = "grp", Numerator = "B", Denominator = "Z",
                    Contrast_name = "B_vs_Z", stringsAsFactors = FALSE)
  expect_error(.pc_de(contrasts = ctr), "B_vs_Z")
})
