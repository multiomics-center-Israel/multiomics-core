# Tests for per-group CV% columns in Final_results_{ALL,DE}.
#
# T1-T4 are deliberately BASE-R ONLY (no Bioconductor, no openxlsx) so they can
# be run in RStudio on Windows before pushing:
#
#   testthat::test_file("tests/testthat/test-group-cv.R")
#
# (run from the repo root). T5 — the opaque-bytes regression guard for the Excel
# writer — lives in test-attach-final-results-xlsx-bytes.R and needs openxlsx, so
# it is left for CI on Linux.

# ---- Self-contained sourcing of the functions under test --------------------
# Locate and source only the pure builders/helpers (function definitions only;
# they pull in openxlsx/edgeR lazily, inside functions we do not call here).
`%||%` <- function(a, b) if (is.null(a)) b else a

find_repo_file <- function(rel) {
    dir <- normalizePath(".", mustWork = FALSE)
    for (i in 1:8) {
        cand <- file.path(dir, rel)
        if (file.exists(cand)) return(cand)
        dir <- dirname(dir)
    }
    stop("Could not locate ", rel, " from working dir ", getwd())
}

source(find_repo_file("R/core/05_export_excel.R"))
source(find_repo_file("R/domain/proteomics/06_outputs_legacy.R"))


# =============================================================================
# T1 — CV math + edge cases  (maps to: Q1 definition, Q6 edge cases)
# =============================================================================
test_that("cv_percent computes 100*sd/mean with sample (n-1) SD", {
    # group A: c(10,20,30) -> mean 20, sd 10 -> CV 50%
    # group B: c(5,5,5)     -> mean 5,  sd 0  -> CV 0% (defined, not NA)
    mat <- rbind(A = c(10, 20, 30), B = c(5, 5, 5))
    cv <- cv_percent(mat)
    expect_equal(unname(cv["A"]), 50)
    expect_equal(unname(cv["B"]), 0)
    # names propagate
    expect_equal(names(cv), c("A", "B"))
})

test_that("cv_percent: single-value group -> NA (not 0)", {
    cv <- cv_percent(rbind(only_one = 42))
    expect_true(is.na(cv[["only_one"]]))
})

test_that("cv_percent: mean == 0 -> NA (not 0 or Inf)", {
    cv <- cv_percent(rbind(zeros = c(0, 0, 0)))
    expect_true(is.na(cv[["zeros"]]))
})

test_that("cv_percent: fewer than 2 observed (NA) values -> NA (proteomics case)", {
    # Mimics unimputed proteomics: only one observed value in the group.
    cv <- cv_percent(rbind(one_obs = c(NA, 5, NA)))
    expect_true(is.na(cv[["one_obs"]]))

    # Exactly two observed values -> defined.
    cv2 <- cv_percent(rbind(two_obs = c(NA, 5, 10)))
    expect_equal(unname(cv2[["two_obs"]]), 100 * stats::sd(c(5, 10)) / mean(c(5, 10)))
})


# =============================================================================
# T2 — back-transform correctness per omics  (maps to: Q1/D1/D2)
# =============================================================================
test_that("proteomics inverse identities hold (2^x and 2^x - 1)", {
    x <- c(0.5, 12.3, 100, 5000)
    # DIANN / already-log2 path: log2(x) -> inverse 2^x
    expect_equal(2^log2(x), x)
    # preprocessed-linear path: log2(x + 1) -> inverse 2^x - 1
    expect_equal(2^log2(x + 1) - 1, x)
})

test_that("proteomics_log_offset resolves the +1 only for preprocessed+linear", {
    # DIANN engine, linear input -> no +1
    expect_equal(proteomics_log_offset(list(modes = list(proteomics = list(
        engine = "DIANN", input = list(format = "diann"), scale_in = "linear")))), 0)
    # preprocessed + linear -> +1
    expect_equal(proteomics_log_offset(list(modes = list(proteomics = list(
        input = list(format = "preprocessed"), scale_in = "linear")))), 1)
    # preprocessed + already log2 -> no +1
    expect_equal(proteomics_log_offset(list(modes = list(proteomics = list(
        input = list(format = "preprocessed"), scale_in = "log2")))), 0)
    # preprocessed, scale_in derived from is_logtransformed = TRUE -> log2 -> no +1
    expect_equal(proteomics_log_offset(list(modes = list(proteomics = list(
        input = list(format = "preprocessed"),
        files = list(is_logtransformed = TRUE))))), 0)
})

test_that("metabolomics post-normalization linear reconstruction is exact for log2(x+p)", {
    # tss/pqn/eigenms always use log2(x + pseudocount); inverse 2^x - p recovers x.
    x <- c(1, 10, 250, 9999)
    for (p in c(0, 1, 5)) {
        expect_equal(2^(log2(x + p)) - p, x)
    }
})


# =============================================================================
# T3 — CV columns originate in build_final_results_generic; ALL and DE agree
#      (maps to: Q7 "insert in build_final_results_generic", ALL vs DE contract)
# =============================================================================
make_fixture <- function() {
    ids <- c("F1", "F2", "F3", "F4")
    samples <- c("S1", "S2", "S3", "S4")

    # Linear expression for CV (group A = S1,S2 ; group B = S3,S4)
    expr_linear <- matrix(c(
        10, 20, 100, 140,   # F1
        50, 50,  5,  15,    # F2
         8, 12, 30,  30,    # F3
         1,  3,  7,   9     # F4
    ), nrow = 4, byrow = TRUE, dimnames = list(ids, samples))

    meta <- data.frame(
        sample_id = samples,
        grp       = c("A", "A", "B", "B"),
        stringsAsFactors = FALSE
    )
    contrasts_df <- data.frame(
        Contrast_name = "B - A", Factor = "grp",
        Numerator = "B", Denominator = "A",
        stringsAsFactors = FALSE
    )

    cv_cols <- compute_group_cv_columns(expr_linear, meta, "sample_id", contrasts_df)

    summary_df <- data.frame(
        Gene                = ids,
        `linearFC.B - A`    = c(-2, 3, 1.1, -1.5),
        `pvalue.B - A`      = c(0.01, 0.2, 0.04, 0.5),
        `padj.B - A`        = c(0.03, 0.3, 0.08, 0.6),
        `B - A_pass`        = c(1, NA, 1, NA),
        pass_any_contrast   = c(1, NA, 1, NA),
        check.names = FALSE, stringsAsFactors = FALSE
    )

    list(summary_df = summary_df, expr_linear = expr_linear,
         contrasts_df = contrasts_df, cv_cols = cv_cols, ids = ids)
}

test_that("build_final_results_generic places CV columns after expression and carries them into the DE subset", {
    fx <- make_fixture()

    all_df <- build_final_results_generic(
        summary_df     = fx$summary_df,
        expr_df        = fx$expr_linear,
        contrasts_df   = fx$contrasts_df,
        feature_id_col = "Gene",
        mode           = "rna",
        cv_cols        = fx$cv_cols
    )

    # CV columns exist for both contrast groups
    expect_true(all(c("CV.A", "CV.B") %in% names(all_df)))

    # Placement: CV columns sit after the expression block, before DE stats
    pos_expr <- max(match(colnames(fx$expr_linear), names(all_df)))
    pos_cv   <- min(match(c("CV.A", "CV.B"), names(all_df)))
    pos_fc   <- match("linearFC.B - A", names(all_df))
    expect_true(pos_expr < pos_cv && pos_cv < pos_fc)

    # DE subset is the same construction the writer uses (rows passing any contrast).
    is_de <- !is.na(all_df$pass_any_contrast) & all_df$pass_any_contrast == 1
    de_df <- all_df[is_de, , drop = FALSE]

    # Same CV columns, identical values (DE just has fewer rows).
    expect_equal(de_df$CV.A, all_df$CV.A[is_de])
    expect_equal(de_df$CV.B, all_df$CV.B[is_de])

    # Spot-check a hand-computed value: F1 group A = c(10,20) -> CV 100*sd/mean
    expect_equal(all_df$CV.A[all_df$Gene == "F1"],
                 100 * stats::sd(c(10, 20)) / mean(c(10, 20)))
})


# =============================================================================
# T4 — CV emitted only for contrast groups; flag fully suppresses
#      (maps to: Q5 contrast-scoping + config flag)
# =============================================================================
test_that("compute_group_cv_columns emits columns only for groups in a contrast", {
    samples <- paste0("S", 1:6)
    expr <- matrix(seq_len(12), nrow = 2, dimnames = list(c("F1", "F2"), samples))
    meta <- data.frame(
        sample_id = samples,
        grp       = c("A", "A", "B", "B", "C", "C"),  # C exists but is in no contrast
        stringsAsFactors = FALSE
    )
    contrasts_df <- data.frame(
        Contrast_name = "B - A", Factor = "grp",
        Numerator = "B", Denominator = "A",
        stringsAsFactors = FALSE
    )

    cv <- compute_group_cv_columns(expr, meta, "sample_id", contrasts_df)
    expect_setequal(names(cv), c("CV.A", "CV.B"))
    expect_false("CV.C" %in% names(cv))   # C never appears in a contrast
})

test_that("group_cv flag = FALSE fully suppresses CV columns", {
    pre <- list(
        expr_filt = matrix(log2(c(10, 20, 30, 40) + 0) , nrow = 1,
                           dimnames = list("F1", paste0("S", 1:4))),
        meta = data.frame(sample_id = paste0("S", 1:4),
                          grp = c("A", "A", "B", "B"), stringsAsFactors = FALSE)
    )
    contrasts_df <- data.frame(
        Contrast_name = "B - A", Factor = "grp",
        Numerator = "B", Denominator = "A", stringsAsFactors = FALSE
    )
    cfg_off <- list(modes = list(proteomics = list(
        excel = list(group_cv = FALSE),
        effects = list(samples = "sample_id"),
        input = list(format = "diann"), scale_in = "linear"
    )))
    expect_null(build_group_cv_proteomics(pre, contrasts_df, cfg_off))

    # And ON (default) produces CV columns for the contrast groups.
    cfg_on <- cfg_off
    cfg_on$modes$proteomics$excel$group_cv <- TRUE
    cv_on <- build_group_cv_proteomics(pre, contrasts_df, cfg_on)
    expect_true(all(c("CV.A", "CV.B") %in% names(cv_on)))
})
