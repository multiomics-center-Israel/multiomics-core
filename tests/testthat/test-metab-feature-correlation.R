# tests/testthat/test-metab-feature-correlation.R
#
# correlate_feature_vs_all() answers "which metabolites behave like this one?".
# The invariants worth pinning are mostly about honesty under bad input: a
# correlation resting on too few shared samples, or on a feature that happens to
# be constant across just those samples, must surface as NA rather than as a
# confident 1.0 at the top of the table.
#
# All synthetic data is seeded so a red test is a real failure, never the RNG.

# ---- helpers ----------------------------------------------------------------

make_mat <- function(rows, samples = NULL) {
  m <- do.call(rbind, rows)
  rownames(m) <- names(rows)
  colnames(m) <- samples %||% paste0("S", seq_len(ncol(m)))
  m
}


# ---- core behaviour ----------------------------------------------------------

test_that("a perfect copy ranks first and a mirrored feature reads as negative", {
  m <- make_mat(list(
    F1 = c(1, 2, 3, 4, 5, 6),
    F2 = c(2, 4, 6, 8, 10, 12),   # perfectly proportional
    F3 = c(6, 5, 4, 3, 2, 1)      # perfectly mirrored
  ))

  res <- correlate_feature_vs_all(m, "F1", min_n = 3)

  expect_equal(res$feature_id[1], "F2")
  expect_equal(res$pearson_r[res$feature_id == "F2"], 1)
  expect_equal(res$pearson_r[res$feature_id == "F3"], -1)
  expect_equal(res$direction[res$feature_id == "F2"], "positive")
  expect_equal(res$direction[res$feature_id == "F3"], "negative")
})

test_that("the query feature never appears among its own results", {
  m <- withr::with_seed(42, make_mat(list(
    F1 = rnorm(10), F2 = rnorm(10), F3 = rnorm(10)
  )))

  res <- correlate_feature_vs_all(m, "F1", min_n = 3)

  expect_false("F1" %in% res$feature_id)
  expect_equal(nrow(res), 2L)
})

test_that("uncorrelated noise does not reach significance", {
  m <- withr::with_seed(7, make_mat(list(
    F1 = rnorm(30), F2 = rnorm(30)
  )))

  res <- correlate_feature_vs_all(m, "F1", min_n = 5)

  expect_gt(res$pearson_pvalue[res$feature_id == "F2"], 0.05)
})

test_that("row order is deterministic and tie-broken by feature_id", {
  # F2 and F3 are both perfectly correlated with F1, so |r| cannot separate them.
  m <- make_mat(list(
    F1 = c(1, 2, 3, 4, 5, 6),
    Fb = c(2, 4, 6, 8, 10, 12),
    Fa = c(3, 6, 9, 12, 15, 18)
  ))

  res <- correlate_feature_vs_all(m, "F1", min_n = 3)
  expect_equal(res$feature_id[1:2], c("Fa", "Fb"))

  # Permuting the input rows must not change the output order.
  res_perm <- correlate_feature_vs_all(m[c(3, 1, 2), , drop = FALSE], "F1", min_n = 3)
  expect_equal(res_perm$feature_id, res$feature_id)
})

test_that("untestable pairs sort last rather than being dropped", {
  m <- make_mat(list(
    F1   = c(1, 2, 3, 4, 5, 6),
    Good = c(2, 4, 6, 8, 10, 12),
    Flat = rep(3, 6)
  ))

  res <- correlate_feature_vs_all(m, "F1", min_n = 3)

  expect_equal(nrow(res), 2L)
  expect_equal(res$feature_id[2], "Flat")
  expect_true(is.na(res$pearson_r[res$feature_id == "Flat"]))
})


# ---- agreement with stats::cor.test ------------------------------------------

test_that("pearson_pvalue matches stats::cor.test(method = 'pearson')", {
  m <- withr::with_seed(11, make_mat(list(
    F1 = rnorm(15),
    F2 = rnorm(15)
  )))

  res <- correlate_feature_vs_all(m, "F1", min_n = 5)
  expected <- stats::cor.test(m["F1", ], m["F2", ], method = "pearson")

  expect_equal(res$pearson_r[1], unname(expected$estimate))
  expect_equal(res$pearson_pvalue[1], expected$p.value)
})

test_that("spearman_pvalue matches cor.test(method='spearman', exact=FALSE), ties included", {
  # F2 is tie-free; F3 is deliberately full of ties, which is where the
  # asymptotic approximation is most likely to drift from the rank statistic.
  m <- withr::with_seed(3, make_mat(list(
    F1 = rnorm(14),
    F2 = rnorm(14),
    F3 = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7)
  )))

  res <- correlate_feature_vs_all(m, "F1", min_n = 5)

  for (fid in c("F2", "F3")) {
    expected <- suppressWarnings(
      stats::cor.test(m["F1", ], m[fid, ], method = "spearman", exact = FALSE)
    )
    row <- res[res$feature_id == fid, ]
    expect_equal(row$spearman_rho, unname(expected$estimate),
                 info = paste("rho for", fid))
    expect_equal(row$spearman_pvalue, expected$p.value,
                 info = paste("p-value for", fid))
  }
})

test_that("a monotonic but non-linear partner separates Spearman from Pearson", {
  x <- 1:8
  m <- make_mat(list(F1 = as.numeric(x), F2 = exp(x)))

  res <- correlate_feature_vs_all(m, "F1", min_n = 5)

  expect_equal(res$spearman_rho[1], 1)   # perfectly monotonic
  expect_lt(res$pearson_r[1], 0.9)       # but clearly not linear
})


# ---- missing data ------------------------------------------------------------

test_that("n_used counts each pair's own shared samples", {
  m <- make_mat(list(
    F1 = c(1, 2, 3, 4, 5, 6, NA, NA),
    F2 = c(1, 2, 3, 4, 5, 6, 7, 8),     # shares 6
    F3 = c(NA, NA, 3, 4, 5, 6, 7, 8)    # shares 4
  ))

  res <- correlate_feature_vs_all(m, "F1", min_n = 3)

  expect_equal(res$n_used[res$feature_id == "F2"], 6L)
  expect_equal(res$n_used[res$feature_id == "F3"], 4L)
})

test_that("min_n blanks an under-powered pair instead of ranking it", {
  m <- make_mat(list(
    F1   = c(1, 2, 3, NA, NA, NA, NA, NA),
    Thin = c(1, 2, 3, 9, 9, 9, 9, 9),   # shares only 3 samples, correlates 1.0
    Wide = c(1, 2, 3, 4, 5, 6, 7, 8)
  ))

  res <- correlate_feature_vs_all(m, "F1", min_n = 5)
  thin <- res[res$feature_id == "Thin", ]

  expect_equal(thin$n_used, 3L)
  expect_true(is.na(thin$pearson_r))
  expect_true(is.na(thin$pearson_padj))
  expect_true(is.na(thin$spearman_rho))
  expect_true(is.na(thin$direction))
})


# ---- the zero-variance guard, pairwise-complete aware ------------------------

test_that("a partner constant only on its shared samples is NA for both coefficients", {
  # FlatShared varies across the full matrix but is constant on S1-S5, which is
  # exactly the window it shares with F1. A global SD check would miss this.
  m <- make_mat(list(
    F1         = c(1, 2, 3, 4, 5, NA, NA, NA),
    FlatShared = c(7, 7, 7, 7, 7, 1, 2, 3),
    Control    = c(1, 2, 3, 4, 5, 9, 9, 9)
  ))

  expect_gt(stats::sd(m["FlatShared", ]), 0)  # not constant overall

  res <- expect_silent(correlate_feature_vs_all(m, "F1", min_n = 5))
  flat <- res[res$feature_id == "FlatShared", ]
  ctrl <- res[res$feature_id == "Control", ]

  expect_equal(flat$n_used, 5L)   # enough samples, so min_n is not what blanked it
  expect_true(is.na(flat$pearson_r))
  expect_true(is.na(flat$pearson_pvalue))
  expect_true(is.na(flat$pearson_padj))
  expect_true(is.na(flat$spearman_rho))
  expect_true(is.na(flat$spearman_pvalue))
  expect_true(is.na(flat$spearman_padj))

  # The guard is per-pair: the control in the same call is unaffected.
  expect_equal(ctrl$pearson_r, 1)
  expect_true(is.finite(ctrl$spearman_rho))
})

test_that("a query constant only on the shared samples is NA for both coefficients", {
  # Mirror of the case above: this time it is the *query* that flatlines on the
  # window it shares with the partner.
  m <- make_mat(list(
    Query   = c(5, 5, 5, 5, 5, 1, 2, 3),
    Partner = c(1, 2, 3, 4, 5, NA, NA, NA),
    Control = c(1, 2, 3, 4, 5, 6, 7, 8)
  ))

  expect_gt(stats::sd(m["Query", ]), 0)

  res <- expect_silent(correlate_feature_vs_all(m, "Query", min_n = 5))
  partner <- res[res$feature_id == "Partner", ]
  ctrl <- res[res$feature_id == "Control", ]

  expect_equal(partner$n_used, 5L)
  expect_true(is.na(partner$pearson_r))
  expect_true(is.na(partner$spearman_rho))

  # Control shares all 8 samples, where the query does vary.
  expect_equal(ctrl$n_used, 8L)
  expect_true(is.finite(ctrl$pearson_r))
})

test_that("a planted -Inf is excluded from Pearson, Spearman and n_used alike", {
  # transform_metab() can emit -Inf by log-transforming a non-positive value.
  # is.na() does not catch it and the payload fill does not replace it, so the
  # engine must treat non-finite as unobserved in BOTH coefficient paths --
  # asserting Pearson alone would let an Inf slip into the Spearman rank fast
  # path, where rank() would score it as a legitimate extreme value.
  m <- make_mat(list(
    Q    = c(1, 2, 3, 4, 5, 6, 7, 8),
    Bad  = c(2, 4, -Inf, 8, 10, 12, 14, 16),
    Good = c(2, 4, 6, 8, 10, 12, 14, 16)
  ))

  res <- expect_silent(correlate_feature_vs_all(m, "Q", min_n = 3))
  bad  <- res[res$feature_id == "Bad", ]
  good <- res[res$feature_id == "Good", ]

  # The -Inf sample is dropped from the pair, not counted and not ranked.
  expect_equal(bad$n_used, 7L)
  expect_equal(good$n_used, 8L)

  # Both coefficients are computed on the 7 usable samples, and both are finite.
  expect_true(is.finite(bad$pearson_r))
  expect_true(is.finite(bad$spearman_rho))
  expect_equal(bad$pearson_r, 1)
  expect_equal(bad$spearman_rho, 1)

  keep <- is.finite(m["Bad", ])
  expect_equal(bad$pearson_r,
               unname(stats::cor(m["Q", keep], m["Bad", keep])))
})

test_that("a -Inf in the query itself is excluded from every partner", {
  m <- make_mat(list(
    Q  = c(1, 2, -Inf, 4, 5, 6, 7, 8),
    P1 = c(2, 4, 6, 8, 10, 12, 14, 16),
    P2 = c(8, 7, 6, 5, 4, 3, 2, 1)
  ))

  res <- expect_silent(correlate_feature_vs_all(m, "Q", min_n = 3))

  expect_true(all(res$n_used == 7L))
  expect_true(all(is.finite(res$pearson_r)))
  expect_true(all(is.finite(res$spearman_rho)))
})

test_that("+Inf routes Spearman to the loop path rather than being ranked", {
  # If the fast path were gated on !anyNA() instead of all(is.finite()), the Inf
  # would reach rank() and be scored as the largest observation, silently
  # inflating rho. Compare against the correlation on the usable samples only.
  m <- make_mat(list(
    Q   = c(1, 2, 3, 4, 5, 6, 7, 8),
    Big = c(1, 2, 3, Inf, 5, 6, 7, 8)
  ))

  res <- expect_silent(correlate_feature_vs_all(m, "Q", min_n = 3))
  keep <- is.finite(m["Big", ])
  expected <- suppressWarnings(
    stats::cor(m["Q", keep], m["Big", keep], method = "spearman")
  )

  expect_equal(res$n_used, 7L)
  expect_equal(res$spearman_rho, unname(expected))
})

test_that("a globally constant feature is NA and raises no warning", {
  m <- make_mat(list(
    F1   = c(1, 2, 3, 4, 5, 6),
    Flat = rep(2, 6)
  ))

  res <- expect_silent(correlate_feature_vs_all(m, "F1", min_n = 3))

  expect_true(is.na(res$pearson_r))
  expect_true(is.na(res$spearman_rho))
  expect_true(is.na(res$direction))
})


# ---- Spearman code paths -----------------------------------------------------

test_that("the complete-matrix fast path and the NA loop path agree", {
  complete <- withr::with_seed(21, make_mat(list(
    F1 = rnorm(12), F2 = rnorm(12), F3 = rnorm(12)
  )))

  # Adding one feature with a missing value forces every pair through the loop,
  # but F2/F3 still share all 12 samples with F1, so their numbers must not move.
  with_na <- rbind(complete, FX = c(NA, withr::with_seed(5, rnorm(11))))

  fast <- correlate_feature_vs_all(complete, "F1", min_n = 5)
  loop <- correlate_feature_vs_all(with_na, "F1", min_n = 5)

  for (fid in c("F2", "F3")) {
    a <- fast[fast$feature_id == fid, ]
    b <- loop[loop$feature_id == fid, ]
    expect_equal(a$pearson_r, b$pearson_r, info = fid)
    expect_equal(a$spearman_rho, b$spearman_rho, info = fid)
    expect_equal(a$spearman_pvalue, b$spearman_pvalue, info = fid)
    expect_equal(a$n_used, b$n_used, info = fid)
  }
})


# ---- multiple-testing correction ---------------------------------------------

test_that("BH is computed over testable pairs only, with the query excluded", {
  m <- withr::with_seed(99, make_mat(list(
    F1   = c(rnorm(9), NA, NA, NA),
    A    = rnorm(12),
    B    = rnorm(12),
    C    = rnorm(12),
    Flat = rep(4, 12),                                  # zero variance
    Thin = c(rnorm(2), rep(NA, 7), rnorm(3))            # too few shared samples
  )))

  res <- correlate_feature_vs_all(m, "F1", min_n = 5)

  testable <- res[!is.na(res$pearson_padj), ]
  expect_setequal(testable$feature_id, c("A", "B", "C"))

  # Exactly p.adjust() over those three raw p-values -- no self, no untested rows.
  expect_equal(
    testable$pearson_padj,
    stats::p.adjust(testable$pearson_pvalue, method = "BH")
  )

  expect_true(is.na(res$pearson_padj[res$feature_id == "Flat"]))
  expect_true(is.na(res$pearson_padj[res$feature_id == "Thin"]))
  expect_equal(nrow(res), nrow(m) - 1L)
})

test_that("a perfect correlation yields a finite p-value, not NaN or Inf", {
  m <- make_mat(list(F1 = c(1, 2, 3, 4, 5, 6), F2 = c(2, 4, 6, 8, 10, 12)))

  res <- correlate_feature_vs_all(m, "F1", min_n = 3)

  expect_equal(abs(res$pearson_r[1]), 1)
  expect_true(is.finite(res$pearson_pvalue[1]))
  expect_equal(res$pearson_pvalue[1], 0)
  expect_true(is.finite(res$pearson_padj[1]))
})


# ---- direction and top_n -----------------------------------------------------

test_that("direction distinguishes zero from positive and from untestable", {
  # F1 and F2 are orthogonal by construction: both are mean-zero and their
  # cross-products cancel exactly, so r is 0 rather than merely small.
  m <- make_mat(list(
    F1   = c(-1, 1, -1, 1),
    F2   = c(1, 1, -1, -1),
    Flat = rep(0, 4)
  ))

  res <- correlate_feature_vs_all(m, "F1", min_n = 3)

  expect_equal(res$pearson_r[res$feature_id == "F2"], 0)
  expect_equal(res$direction[res$feature_id == "F2"], "none")
  expect_true(is.na(res$direction[res$feature_id == "Flat"]))
})

test_that("top_n is applied after self-removal and sorting", {
  m <- make_mat(list(
    F1 = c(1, 2, 3, 4, 5, 6),
    A  = c(2, 4, 6, 8, 10, 12),      # r = 1
    B  = c(6, 5, 4, 3, 2, 1),        # r = -1
    C  = c(1, 1, 2, 2, 9, 1)         # weak
  ))

  res <- correlate_feature_vs_all(m, "F1", min_n = 3, top_n = 2)

  expect_equal(nrow(res), 2L)
  expect_false("F1" %in% res$feature_id)
  expect_setequal(res$feature_id, c("A", "B"))
})


# ---- input validation --------------------------------------------------------

test_that("an unknown feature errors informatively and suggests near matches", {
  m <- make_mat(list(Glucose = 1:6, Fructose = 6:1, Lactate = c(1, 3, 2, 4, 3, 5)))

  expect_error(correlate_feature_vs_all(m, "Glucoze", min_n = 3), "Glucoze")
  expect_error(correlate_feature_vs_all(m, "Glucoze", min_n = 3), "Glucose")
  expect_error(correlate_feature_vs_all(m, "totally_absent", min_n = 3),
               "not in the matrix")
})

test_that("degenerate inputs are rejected with actionable messages", {
  m <- make_mat(list(F1 = 1:6, F2 = 6:1))

  expect_error(correlate_feature_vs_all(m, c("F1", "F2")), "single")
  expect_error(correlate_feature_vs_all(m, "F1", min_n = 2), "min_n")
  expect_error(correlate_feature_vs_all(m, "F1", top_n = 0), "top_n")
  expect_error(correlate_feature_vs_all(m["F1", , drop = FALSE], "F1"),
               "at least 2 features")

  unnamed <- matrix(1:12, nrow = 2)
  expect_error(correlate_feature_vs_all(unnamed, "F1"), "rownames")
})


# ---- restore_missing_values --------------------------------------------------

test_that("a mask reinstates NA exactly where TRUE", {
  m <- make_mat(list(F1 = c(1, 2, 3, 4), F2 = c(5, 6, 7, 8)))
  mask <- matrix(c(FALSE, TRUE, FALSE, FALSE,
                   FALSE, FALSE, TRUE, FALSE),
                 nrow = 2, byrow = TRUE, dimnames = dimnames(m))

  out <- restore_missing_values(m, mask)

  expect_true(is.na(out["F1", "S2"]))
  expect_true(is.na(out["F2", "S3"]))
  expect_equal(sum(is.na(out)), 2L)
  expect_equal(out["F1", "S1"], 1)   # untouched cells unchanged
})

test_that("a NULL mask is a no-op, so the GUI can call it unconditionally", {
  m <- make_mat(list(F1 = c(1, 2, 3, 4), F2 = c(5, 6, 7, 8)))
  expect_identical(restore_missing_values(m, NULL), m)
})

test_that("an all-FALSE mask restores nothing (verified-complete payload)", {
  m <- make_mat(list(F1 = c(1, 2, 3, 4), F2 = c(5, 6, 7, 8)))
  mask <- matrix(FALSE, nrow = 2, ncol = 4, dimnames = dimnames(m))
  expect_identical(restore_missing_values(m, mask), m)
})

test_that("a malformed mask is rejected rather than silently mis-indexing", {
  m <- make_mat(list(F1 = c(1, 2, 3, 4), F2 = c(5, 6, 7, 8)))
  good <- matrix(FALSE, nrow = 2, ncol = 4, dimnames = dimnames(m))

  numeric_mask <- matrix(0, nrow = 2, ncol = 4, dimnames = dimnames(m))
  expect_error(restore_missing_values(m, numeric_mask), "logical")

  na_mask <- good; na_mask[1, 1] <- NA
  expect_error(restore_missing_values(m, na_mask), "contains NA")

  expect_error(restore_missing_values(m, t(good)), "mask is")

  renamed <- good; colnames(renamed) <- paste0("X", 1:4)
  expect_error(restore_missing_values(m, renamed), "dimnames")
})

test_that("mask -> fill -> restore round-trips to the un-filled statistics", {
  # The test that actually pins the fix: it fails if the mask is captured after
  # the payload builder's row-median fill rather than before it.
  original <- withr::with_seed(4, make_mat(list(
    Q  = c(rnorm(6), NA, NA),
    P1 = rnorm(8),
    P2 = c(rnorm(5), NA, rnorm(2))
  )))

  mask <- !is.finite(original)

  filled <- original
  for (i in seq_len(nrow(filled))) {
    row_i <- filled[i, ]; na_i <- is.na(row_i)
    if (any(na_i)) {
      md <- stats::median(row_i, na.rm = TRUE)
      filled[i, na_i] <- if (is.na(md)) 0 else md
    }
  }
  expect_false(anyNA(filled))                       # as the payload ships it

  restored <- restore_missing_values(filled, mask)
  expect_equal(restored, original)

  from_restored <- correlate_feature_vs_all(restored, "Q", min_n = 3)
  from_original <- correlate_feature_vs_all(original, "Q", min_n = 3)
  expect_equal(from_restored, from_original)

  # And the filled matrix on its own gives different, wrong-n statistics.
  from_filled <- correlate_feature_vs_all(filled, "Q", min_n = 3)
  expect_true(all(from_filled$n_used == ncol(original)))
  expect_false(isTRUE(all.equal(from_filled$n_used, from_original$n_used)))
})


# ---- prepare_correlation_matrix ----------------------------------------------

test_that("QC, blank and pool samples are removed before correlating", {
  meta <- data.frame(
    sample_id = c("S1", "S2", "S3", "S4", "QC1", "Blank1"),
    condition = c("ctrl", "ctrl", "trt", "trt", "QC", "blank"),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample_id
  m <- make_mat(list(F1 = 1:6, F2 = 6:1), samples = meta$sample_id)

  out <- prepare_correlation_matrix(m, meta, condition_col = "condition")

  expect_equal(colnames(out), c("S1", "S2", "S3", "S4"))
  expect_equal(nrow(out), 2L)
})

test_that("an all-biological matrix passes through untouched", {
  meta <- data.frame(
    sample_id = paste0("S", 1:5),
    condition = c("a", "a", "b", "b", "b"),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample_id
  m <- make_mat(list(F1 = 1:5, F2 = 5:1), samples = meta$sample_id)

  out <- prepare_correlation_matrix(m, meta, condition_col = "condition")

  expect_equal(dim(out), dim(m))
  expect_equal(colnames(out), colnames(m))
})

test_that("sample IDs are taken from rownames when the column is absent", {
  meta <- data.frame(condition = c("a", "a", "b", "b"), stringsAsFactors = FALSE)
  rownames(meta) <- paste0("S", 1:4)
  m <- make_mat(list(F1 = 1:4, F2 = 4:1), samples = rownames(meta))

  out <- prepare_correlation_matrix(m, meta, condition_col = "condition")
  expect_equal(colnames(out), paste0("S", 1:4))
})

test_that("a configured qc_flag_column is honoured when supplied", {
  meta <- data.frame(
    sample_id = paste0("S", 1:5),
    condition = c("a", "a", "b", "b", "b"),
    run_type  = c("sample", "sample", "sample", "sample", "pool"),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample_id
  m <- make_mat(list(F1 = 1:5, F2 = 5:1), samples = meta$sample_id)

  kept <- prepare_correlation_matrix(m, meta, condition_col = "condition")
  expect_equal(ncol(kept), 5L)  # not detectable without the flag column

  filtered <- prepare_correlation_matrix(m, meta, condition_col = "condition",
                                         qc_flag_column = "run_type")
  expect_equal(colnames(filtered), paste0("S", 1:4))
})

test_that("misaligned or over-filtered inputs fail loudly", {
  meta <- data.frame(
    sample_id = c("S1", "S2", "QC1", "QC2"),
    condition = c("a", "b", "QC", "QC"),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample_id
  m <- make_mat(list(F1 = 1:4, F2 = 4:1), samples = meta$sample_id)

  # Only 2 biological samples survive -- too few to correlate.
  expect_error(prepare_correlation_matrix(m, meta, condition_col = "condition"),
               "not enough to")

  expect_error(
    prepare_correlation_matrix(m, meta, condition_col = "nope"),
    "not found"
  )

  orphan <- m
  colnames(orphan) <- c("S1", "S2", "S9", "S10")
  expect_error(prepare_correlation_matrix(orphan, meta, condition_col = "condition"),
               "no row in")
})


# ---- sample-ID resolution (both public functions) ----------------------------

# The shipped metabolomics and proteomics templates set effects$samples to
# "SampleID", and the payload does not expose that name -- it exposes rownames.
# A "sample_id" default therefore silently mismatched the standard template
# (and hard-errored inside the plot). These pin the NULL-default behaviour.

sampleid_fixture <- function() {
  meta <- data.frame(
    SampleID  = paste0("S", 1:6),
    treatment = c("ctrl", "ctrl", "mid", "mid", "high", "high"),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$SampleID           # as the payload builder sets it
  m <- make_mat(list(
    F1 = c(1, 2, 3, 4, 5, 6),
    F2 = c(2, 3, 4, 5, 6, 7)
  ), samples = meta$SampleID)
  list(m = m, meta = meta)
}

test_that("a SampleID-style payload works with sample_col omitted", {
  fx <- sampleid_fixture()

  out <- prepare_correlation_matrix(fx$m, fx$meta, condition_col = "treatment")
  expect_equal(colnames(out), paste0("S", 1:6))

  p <- plot_feature_correlation_profiles(
    fx$m, fx$meta, feature_id = "F1", partner_ids = "F2",
    group_col = "treatment"
  )
  expect_s3_class(p, "ggplot")
})

test_that("an explicit sample_col is honoured by both functions", {
  fx <- sampleid_fixture()

  out <- prepare_correlation_matrix(fx$m, fx$meta, condition_col = "treatment",
                                    sample_col = "SampleID")
  expect_equal(colnames(out), paste0("S", 1:6))

  p <- plot_feature_correlation_profiles(
    fx$m, fx$meta, feature_id = "F1", partner_ids = "F2",
    group_col = "treatment", sample_col = "SampleID"
  )
  expect_s3_class(p, "ggplot")
})

test_that("a bogus sample_col errors naming the available columns", {
  fx <- sampleid_fixture()
  expect_error(
    prepare_correlation_matrix(fx$m, fx$meta, condition_col = "treatment",
                               sample_col = "nope"),
    "SampleID"
  )
})

test_that("integer-defaulted rownames error telling the caller to pass sample_col", {
  meta <- data.frame(SampleID = paste0("S", 1:4),
                     treatment = c("a", "a", "b", "b"),
                     stringsAsFactors = FALSE)   # rownames are "1".."4"
  m <- make_mat(list(F1 = 1:4, F2 = 4:1), samples = meta$SampleID)

  expect_error(
    prepare_correlation_matrix(m, meta, condition_col = "treatment"),
    "Pass `sample_col`"
  )
})

test_that("duplicate and NA IDs are rejected via BOTH routes", {
  # The explicit path must carry the same guarantees as the rownames fallback:
  # a duplicated ID mis-aligns samples just as silently either way.
  m <- make_mat(list(F1 = 1:4, F2 = 4:1), samples = paste0("S", 1:4))

  dup <- data.frame(SampleID = c("S1", "S2", "S3", "S3"),
                    treatment = c("a", "a", "b", "b"),
                    stringsAsFactors = FALSE)
  expect_error(
    prepare_correlation_matrix(m, dup, condition_col = "treatment",
                               sample_col = "SampleID"),
    "not unique"
  )

  na_ids <- data.frame(SampleID = c("S1", "S2", NA, "S4"),
                       treatment = c("a", "a", "b", "b"),
                       stringsAsFactors = FALSE)
  expect_error(
    prepare_correlation_matrix(m, na_ids, condition_col = "treatment",
                               sample_col = "SampleID"),
    "missing or"
  )

  # And via rownames: make.unique-free duplicates are impossible in rownames, so
  # the rownames route is exercised by the coverage check instead.
  mismatched <- data.frame(treatment = c("a", "a", "b", "b"),
                           stringsAsFactors = FALSE)
  rownames(mismatched) <- c("S1", "S2", "S9", "S10")
  expect_error(
    prepare_correlation_matrix(m, mismatched, condition_col = "treatment"),
    "no row in"
  )
})


# ---- plot_feature_correlation_profiles ---------------------------------------

plot_fixture <- function() {
  meta <- data.frame(
    sample_id = paste0("S", 1:6),
    condition = c("ctrl", "ctrl", "mid", "mid", "high", "high"),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample_id
  m <- make_mat(list(
    F1   = c(1, 2, 3, 4, 5, 6),
    F2   = c(2, 3, 4, 5, 6, 7),
    F3   = c(6, 5, 4, 3, 2, 1),
    Flat = rep(3, 6)
  ), samples = meta$sample_id)
  list(m = m, meta = meta)
}

test_that("the profile plot returns a ggplot object", {
  fx <- plot_fixture()

  p <- plot_feature_correlation_profiles(
    fx$m, fx$meta, feature_id = "F1", partner_ids = c("F2", "F3"),
    group_col = "condition"
  )

  expect_s3_class(p, "ggplot")
})

test_that("unknown partner IDs are dropped with a message, not an error", {
  fx <- plot_fixture()

  expect_message(
    p <- plot_feature_correlation_profiles(
      fx$m, fx$meta, feature_id = "F1",
      partner_ids = c("F2", "does_not_exist"),
      group_col = "condition"
    ),
    "dropping 1 partner"
  )
  expect_s3_class(p, "ggplot")
})

test_that("an unknown query feature errors in the plot too", {
  fx <- plot_fixture()

  expect_error(
    plot_feature_correlation_profiles(
      fx$m, fx$meta, feature_id = "nope", partner_ids = "F2",
      group_col = "condition"
    ),
    "not in the matrix"
  )
})

test_that("label_map degrades safely when NULL, partial or unnamed", {
  fx <- plot_fixture()

  args <- list(fx$m, fx$meta, feature_id = "F1", partner_ids = c("F2", "F3"),
               group_col = "condition")

  expect_s3_class(do.call(plot_feature_correlation_profiles, args), "ggplot")

  partial <- do.call(plot_feature_correlation_profiles,
                     c(args, list(label_map = c(F1 = "Glucose"))))
  expect_s3_class(partial, "ggplot")

  # An unnamed map would index by position and mislabel everything; it must be
  # ignored rather than trusted.
  unnamed <- do.call(plot_feature_correlation_profiles,
                     c(args, list(label_map = c("Glucose", "Lactate"))))
  expect_s3_class(unnamed, "ggplot")
})

test_that("a group with no observed values becomes a gap, never a fabricated point", {
  # Reachable only once missingness is restored: every sample in one group is
  # unobserved for this feature, so its group mean is undefined. It must be
  # drawn as a break in the line -- interpolating would show the reader a
  # measurement nobody took.
  fx <- plot_fixture()
  m <- fx$m
  m["F2", c("S1", "S2")] <- NA          # the whole "ctrl" group

  expect_message(
    p <- plot_feature_correlation_profiles(
      m, fx$meta, feature_id = "F1", partner_ids = "F2",
      group_col = "condition"
    ),
    "no observed values"
  )
  expect_s3_class(p, "ggplot")

  built <- ggplot2::ggplot_build(p)
  ys <- unlist(lapply(built$data, function(d) d$y))
  expect_false(any(is.nan(ys)))
  expect_false(any(is.infinite(ys)))

  # F2 still contributes its other two groups -- the gap is one point, not the
  # whole series.
  f2_layer <- built$data[[1]]
  expect_true(sum(!is.na(f2_layer$y)) >= 2)
})

test_that("a -Inf in the matrix does not reach the plot as a coordinate", {
  fx <- plot_fixture()
  m <- fx$m
  m["F2", "S3"] <- -Inf

  p <- plot_feature_correlation_profiles(
    m, fx$meta, feature_id = "F1", partner_ids = "F2",
    group_col = "condition"
  )
  built <- ggplot2::ggplot_build(p)
  ys <- unlist(lapply(built$data, function(d) d$y))
  expect_false(any(is.nan(ys)))
  expect_false(any(is.infinite(ys)))
})

test_that("a flat profile is drawn rather than crashing or producing NaN", {
  fx <- plot_fixture()

  p <- plot_feature_correlation_profiles(
    fx$m, fx$meta, feature_id = "F1", partner_ids = c("F2", "Flat"),
    group_col = "condition"
  )

  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  ys <- unlist(lapply(built$data, function(d) d$y))
  expect_false(any(is.nan(ys)))
})
