# Tests for per-group CV% columns in Final_results_{ALL,DE}.
#
# T1-T4 + the wiring tests (W1-W3) are deliberately BASE-R ONLY (no
# Bioconductor, no openxlsx) so they can be run in RStudio on Windows before
# pushing:
#
#   testthat::test_file("tests/testthat/test-group-cv.R")
#
# (run from the repo root). T5 — the opaque-bytes regression guard for the Excel
# writer — lives in test-attach-final-results-xlsx-bytes.R and needs openxlsx, so
# it is left for CI on Linux.
#
# Test map:
#   T1 cv_percent math + edge cases (single-sample/mean-0/<2-observed -> NA)
#   T2 back-transform identities + proteomics_log_offset resolution
#   T3 CV originates in build_final_results_generic; ALL and DE agree
#   T4 CV emitted only for contrast groups; group_cv flag fully suppresses
#   W1 happy-path real-config wiring per omics (sample_id_col + Factor actually
#      resolve against the matrix/meta columns -> finite CV, not all-NA)
#   W2 sample-ID mismatch -> graceful all-NA (the symptom of a misconfigured ID
#      column), never an error
#   W3 multi-Factor contrasts -> warning + NULL (documented limitation)

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
source(find_repo_file("R/domain/rnaseq/01_expression.R"))       # compute_cpm (base R)
source(find_repo_file("R/domain/rnaseq/05_outputs_legacy.R"))   # build_group_cv_rnaseq
source(find_repo_file("R/domain/metabolomics/05_outputs_legacy.R")) # build_group_cv_metabolomics


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


# =============================================================================
# W1 — Happy-path real-config wiring per omics
#      (maps to: the silent failure mode where a misresolved sample_id_col /
#       Factor makes col_grp all-NA and every CV column comes out all-NA)
#
# Each omics resolves sample_id_col from a DIFFERENT config path
# (modes.<omics>.effects$samples) and the group from contrasts_df$Factor. These
# tests build a fixture that mimics a real config — a realistic sample-ID column
# name, a separate condition/group column, an expression matrix whose colnames
# ARE those sample IDs, and a structured contrasts_df referencing the group
# column — then run the actual per-omics builder and assert the CV columns carry
# real finite numbers (not all-NA). A sample-ID/Factor mismatch would surface
# here as all-NA, so these assertions are the wiring guard.
# =============================================================================

# Shared wiring fixture: 6 samples, two 3-replicate groups (ctrl/trt).
wiring_meta <- function(sample_id_col = "SampleID", group_col = "condition") {
    m <- data.frame(
        x_id  = paste0("sample_", 1:6),
        x_grp = c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"),
        stringsAsFactors = FALSE
    )
    names(m) <- c(sample_id_col, group_col)
    m
}
wiring_contrasts <- function(group_col = "condition") {
    data.frame(
        Contrast_name = "trt - ctrl", Factor = group_col,
        Numerator = "trt", Denominator = "ctrl",
        stringsAsFactors = FALSE
    )
}
# A linear feature x sample matrix with within-group variation (so CV is finite
# and > 0 for the multi-replicate groups).
wiring_linear_matrix <- function(sample_ids) {
    m <- matrix(c(
        100, 120, 110, 300, 330, 360,   # F1
         50,  55,  45, 200, 180, 220,   # F2
         10,  11,   9,  40,  42,  38    # F3
    ), nrow = 3, byrow = TRUE,
    dimnames = list(c("F1", "F2", "F3"), sample_ids))
    m
}

test_that("W1 rnaseq: real-config wiring yields finite CV (linear CPM)", {
    meta <- wiring_meta()                          # SampleID / condition
    counts <- wiring_linear_matrix(meta$SampleID)  # raw counts (compute_cpm input)
    # source_type "matrix" => genuine raw counts, on the CPM whitelist.
    pre <- list(expr_filt = counts, meta = meta, info = list(source_type = "matrix"))
    config <- list(modes = list(rna = list(
        effects = list(samples = "SampleID")        # resolved sample_id_col
        # excel$group_cv defaults to TRUE
    )))

    cv <- build_group_cv_rnaseq(pre, wiring_contrasts(), config)
    expect_false(is.null(cv))
    expect_setequal(names(cv), c("CV.trt", "CV.ctrl"))
    # The wiring is correct only if real numbers come out (NOT all-NA).
    expect_true(all(is.finite(cv$CV.ctrl)))
    expect_true(all(is.finite(cv$CV.trt)))
    expect_true(any(cv$CV.ctrl > 0))
})

test_that("W1 proteomics: real-config wiring yields finite CV (2^expr_filt)", {
    meta <- wiring_meta()
    lin  <- wiring_linear_matrix(meta$SampleID)
    pre  <- list(expr_filt = log2(lin), meta = meta)  # log2 intensities (DIANN path)
    config <- list(modes = list(proteomics = list(
        effects = list(samples = "SampleID"),
        input = list(format = "diann"), scale_in = "linear"
    )))

    cv <- build_group_cv_proteomics(pre, wiring_contrasts(), config)
    expect_false(is.null(cv))
    expect_setequal(names(cv), c("CV.trt", "CV.ctrl"))
    expect_true(all(is.finite(cv$CV.ctrl)))
    expect_true(all(is.finite(cv$CV.trt)))
    # Back-transform round-trips to the linear matrix, so CV matches linear CV.
    expect_equal(cv$CV.ctrl[1], 100 * stats::sd(lin[1, 1:3]) / mean(lin[1, 1:3]))
})

test_that("W1 metabolomics: real-config wiring yields finite CV (2^expr_work - p)", {
    meta <- wiring_meta()
    lin  <- wiring_linear_matrix(meta$SampleID)
    pseudo <- 1
    pre  <- list(expr_work = log2(lin + pseudo), meta = meta)  # post-norm log2(x+p)
    config <- list(modes = list(metabolomics = list(
        effects = list(samples = "SampleID"),
        preprocessing = list(chosen_norm = "pqn"),   # always-log2 -> provably invertible
        normalization = list(scaling = "none", pseudocount = pseudo)
    )))

    cv <- build_group_cv_metabolomics(pre, wiring_contrasts(), config)
    expect_false(is.null(cv))
    expect_setequal(names(cv), c("CV.trt", "CV.ctrl"))
    expect_true(all(is.finite(cv$CV.ctrl)))
    expect_true(all(is.finite(cv$CV.trt)))
    # Reconstruction is exact here -> CV equals linear CV.
    expect_equal(cv$CV.trt[1], 100 * stats::sd(lin[1, 4:6]) / mean(lin[1, 4:6]))
})


# =============================================================================
# W2 — Sample-ID mismatch is graceful (all-NA), never an error
#      (maps to: documenting that an all-NA CV column == a misconfigured ID col)
# =============================================================================
test_that("W2 sample-ID/matrix mismatch -> all-NA CV columns, no error", {
    meta <- wiring_meta()                       # IDs: sample_1..sample_6
    lin  <- wiring_linear_matrix(meta$SampleID)
    # Break the wiring: matrix colnames no longer match meta[[sample_id_col]].
    colnames(lin) <- paste0("WRONG_", seq_len(ncol(lin)))

    # IMPORTANT: when the resolved sample_id_col does not match the matrix
    # column names, every sample maps to group NA, so each group's CV is all-NA.
    # This is the documented symptom of a misconfigured ID column — readers who
    # see an entirely-NA CV column should suspect a sample_id_col mismatch.
    cv <- NULL
    expect_error(
        cv <- compute_group_cv_columns(lin, meta, "SampleID", wiring_contrasts()),
        NA  # NA => assert that NO error is thrown
    )
    expect_setequal(names(cv), c("CV.trt", "CV.ctrl"))
    expect_true(all(is.na(cv$CV.ctrl)))
    expect_true(all(is.na(cv$CV.trt)))
})


# =============================================================================
# W3 — Multiple distinct Factors in the contrasts table is unsupported
#      (maps to: documented limitation — one grouping Factor per run)
# =============================================================================
test_that("W3 multi-Factor contrasts -> warning + NULL (no CV columns)", {
    meta <- data.frame(
        SampleID  = paste0("sample_", 1:4),
        condition = c("ctrl", "ctrl", "trt", "trt"),
        genotype  = c("wt", "ko", "wt", "ko"),
        stringsAsFactors = FALSE
    )
    lin <- matrix(c(
        100, 120, 300, 330,   # F1
         50,  55, 200, 180    # F2
    ), nrow = 2, byrow = TRUE,
    dimnames = list(c("F1", "F2"), meta$SampleID))

    # Two contrasts using two DIFFERENT Factor columns (both present in meta).
    contrasts_df <- data.frame(
        Contrast_name = c("trt - ctrl", "ko - wt"),
        Factor        = c("condition", "genotype"),
        Numerator     = c("trt", "ko"),
        Denominator   = c("ctrl", "wt"),
        stringsAsFactors = FALSE
    )

    expect_warning(
        res <- compute_group_cv_columns(lin, meta, "SampleID", contrasts_df),
        "single grouping column"
    )
    expect_null(res)
})


# =============================================================================
# G1 — Metabolomics log2-workspace guard (maps to: review finding #1)
#      The 2^x reconstruction is only valid on a provably-log2 workspace. The
#      only reachable non-log2 expr_work is chosen_norm="median" with a non-log2
#      transform; that (and any unknown chosen_norm) must SKIP + warn, never
#      produce wrong CV by inverting a non-log2 matrix.
# =============================================================================
test_that("G1 metab: median + non-log2 transform -> NULL + warning (no wrong CV)", {
    meta <- wiring_meta()
    lin  <- wiring_linear_matrix(meta$SampleID)
    # expr_work happens to be on a log10 base here; the guard must refuse to 2^ it.
    pre  <- list(expr_work = log10(lin + 1), meta = meta)
    config <- list(modes = list(metabolomics = list(
        effects = list(samples = "SampleID"),
        preprocessing = list(chosen_norm = "median"),
        normalization = list(scaling = "none", transform = "log10", pseudocount = 1)
    )))

    expect_warning(
        cv <- build_group_cv_metabolomics(pre, wiring_contrasts(), config),
        "provably-log2"
    )
    expect_null(cv)
})

test_that("G1 metab: median + log2 transform IS invertible -> finite CV", {
    meta <- wiring_meta()
    lin  <- wiring_linear_matrix(meta$SampleID)
    pseudo <- 1
    pre  <- list(expr_work = log2(lin + pseudo), meta = meta)
    config <- list(modes = list(metabolomics = list(
        effects = list(samples = "SampleID"),
        preprocessing = list(chosen_norm = "median"),
        normalization = list(scaling = "none", transform = "log2", pseudocount = pseudo)
    )))

    cv <- build_group_cv_metabolomics(pre, wiring_contrasts(), config)
    expect_false(is.null(cv))
    expect_setequal(names(cv), c("CV.trt", "CV.ctrl"))
    expect_true(all(is.finite(cv$CV.ctrl)))
    expect_equal(cv$CV.trt[1], 100 * stats::sd(lin[1, 4:6]) / mean(lin[1, 4:6]))
})

test_that("G1 metab: unknown chosen_norm -> NULL + warning (fail-safe)", {
    meta <- wiring_meta()
    lin  <- wiring_linear_matrix(meta$SampleID)
    pre  <- list(expr_work = log2(lin + 1), meta = meta)
    config <- list(modes = list(metabolomics = list(
        effects = list(samples = "SampleID"),
        preprocessing = list(chosen_norm = "some_future_method"),
        normalization = list(scaling = "none", transform = "log2", pseudocount = 1)
    )))

    expect_warning(
        cv <- build_group_cv_metabolomics(pre, wiring_contrasts(), config),
        "provably-log2"
    )
    expect_null(cv)
})


# =============================================================================
# G2 — RNA-seq raw-counts source-type whitelist (maps to: review finding #3)
#      compute_cpm() (counts / colSums) is only valid on genuine raw counts.
#      Only "matrix"/"tximport" are whitelisted; everything else (incl.
#      "preprocessed" and any unknown/missing source_type) SKIPS + warns.
# =============================================================================
test_that("G2 rnaseq: preprocessed source_type -> NULL + warning (CPM meaningless)", {
    meta <- wiring_meta()
    vals <- wiring_linear_matrix(meta$SampleID)   # already-processed values, not counts
    pre <- list(expr_filt = vals, meta = meta,
                info = list(source_type = "preprocessed"))
    config <- list(modes = list(rna = list(effects = list(samples = "SampleID"))))

    expect_warning(
        cv <- build_group_cv_rnaseq(pre, wiring_contrasts(), config),
        "not a known raw-count type"
    )
    expect_null(cv)
})

test_that("G2 rnaseq: missing source_type -> NULL + warning (fail-safe whitelist)", {
    meta <- wiring_meta()
    counts <- wiring_linear_matrix(meta$SampleID)
    pre <- list(expr_filt = counts, meta = meta)   # no info$source_type at all
    config <- list(modes = list(rna = list(effects = list(samples = "SampleID"))))

    expect_warning(
        cv <- build_group_cv_rnaseq(pre, wiring_contrasts(), config),
        "not a known raw-count type"
    )
    expect_null(cv)
})

test_that("G2 rnaseq: tximport source_type IS whitelisted -> finite CV", {
    meta <- wiring_meta()
    counts <- wiring_linear_matrix(meta$SampleID)
    pre <- list(expr_filt = counts, meta = meta,
                info = list(source_type = "tximport"))
    config <- list(modes = list(rna = list(effects = list(samples = "SampleID"))))

    cv <- build_group_cv_rnaseq(pre, wiring_contrasts(), config)
    expect_false(is.null(cv))
    expect_setequal(names(cv), c("CV.trt", "CV.ctrl"))
    expect_true(all(is.finite(cv$CV.ctrl)))
})
