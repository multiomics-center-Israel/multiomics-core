# tests/testthat/test-de-reconciliation.R
#
# The DE_reconciliation sheet: build_de_reconciliation_proteomics() in
# R/domain/proteomics/06_outputs_legacy.R, and the writer that puts it in the
# workbook (R/core/05_export_excel.R).
#
# Multi-imputation reports log2FC.imputs = log2( mean( 2^logFC ) ) across
# independently imputed fits. That rule is invisible in the results table, so
# this sheet shows it: the per-run coefficients, their ratios, the mean of
# those ratios, and the reported value beside it.
#
# The expected values here are NOT a second implementation of the pooling
# rule. Every assertion compares the reconciliation against what
# summarize_limma_mult_imputation() actually produced from the same runs, so
# if the pooling rule ever changes these tests follow it instead of pinning a
# copy that would silently drift.
#
# All fixtures are synthetic (p1..p4, S_1/S_2/NS_1/NS_2).

recon_meta <- function() {
    data.frame(
        SampleID  = c("S_1", "S_2", "NS_1", "NS_2"),
        condition = c("S", "S", "NS", "NS"),
        stringsAsFactors = FALSE
    )
}

recon_contrasts <- function(names = "S_vs_NS") {
    data.frame(
        Contrast_name = names,
        Factor        = "condition",
        Numerator     = "S",
        Denominator   = "NS",
        stringsAsFactors = FALSE
    )
}

recon_config <- function(n_reps = 3, multi = TRUE) {
    list(modes = list(proteomics = list(
        de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5, use_adj_for_pass1 = TRUE),
        imputation = list(multi_imputation = multi, no_repetitions = n_reps,
                          min_no_passed = 1),
        de_table = list(id_col = "FeatureID"),
        effects = list(samples = "SampleID")
    )))
}

# One DE table per run per contrast. `spread` controls how far the runs
# disagree and is recycled over features, so a vector gives a mixed fixture:
# 0 for a feature every run agreed on (what a fully measured feature looks
# like, since there was nothing to impute) and non-zero for one they did not.
recon_runs <- function(n_runs = 3, contrasts = "S_vs_NS", spread = 0.05,
                       features = c("p1", "p2", "p3", "p4")) {
    base_lfc <- c(-0.7, 2, 0.1, -1.4)[seq_along(features)]
    lapply(seq_len(n_runs), function(i) {
        per_contrast <- lapply(seq_along(contrasts), function(ci) {
            data.frame(
                FeatureID = features,
                # Offset per run and per contrast, so no two columns of the
                # reconciliation can coincide by accident.
                logFC     = base_lfc + spread * (i - 1) + 0.3 * (ci - 1),
                P.Value   = rep(1e-4, length(features)),
                adj.P.Val = rep(1e-3, length(features)),
                stringsAsFactors = FALSE
            )
        })
        names(per_contrast) <- contrasts
        per_contrast
    })
}

# A de_res as the pipeline builds it: the runs, and the summary the production
# summariser derives from exactly those runs.
recon_de_res <- function(n_runs = 3, contrasts = "S_vs_NS", spread = 0.05,
                         config = NULL, features = c("p1", "p2", "p3", "p4")) {
    cfg <- config %||% recon_config(n_reps = n_runs)
    runs <- recon_runs(n_runs, contrasts, spread, features)
    list(
        runs_de_tables = runs,
        summary_df     = summarize_limma_mult_imputation(runs, cfg)
    )
}


# =============================================================================
# The reconciliation reproduces what the summariser reported
# =============================================================================

test_that("mean of the per-run ratios reproduces the reported linearRatio.imputs", {
    de_res <- recon_de_res()
    rec <- build_de_reconciliation_proteomics(de_res, recon_config())

    ratio_cols <- grep("^run[0-9]+\\.ratio$", names(rec), value = TRUE)
    expect_length(ratio_cols, 3L)

    by_hand <- rowMeans(as.matrix(rec[, ratio_cols, drop = FALSE]))
    expect_equal(rec$mean.ratio, by_hand)
    # ...and that mean is what the summariser reported, not merely what this
    # sheet recomputed.
    expect_equal(rec$mean.ratio, rec$linearRatio.imputs)
})

test_that("log2 of that mean reproduces the reported log2FC.imputs", {
    rec <- build_de_reconciliation_proteomics(recon_de_res(), recon_config())

    expect_equal(rec$log2FC.from_mean_ratio, log2(rec$mean.ratio))
    expect_equal(rec$log2FC.from_mean_ratio, rec$log2FC.imputs)
})

test_that("both deltas are zero to floating-point tolerance", {
    rec <- build_de_reconciliation_proteomics(recon_de_res(), recon_config())

    # These are reconciliation checks against the production summariser. A
    # non-zero value means the sheet and the results table disagree.
    expect_equal(rec$delta.linearRatio, rep(0, nrow(rec)), tolerance = 1e-12)
    expect_equal(rec$delta.log2FC, rep(0, nrow(rec)), tolerance = 1e-12)
})

test_that("run columns carry the per-run coefficients and their ratios", {
    de_res <- recon_de_res()
    rec <- build_de_reconciliation_proteomics(de_res, recon_config())
    p1 <- rec[rec$FeatureID == "p1", ]

    for (i in 1:3) {
        expected <- de_res$runs_de_tables[[i]][["S_vs_NS"]]
        expected <- expected$logFC[expected$FeatureID == "p1"]
        expect_equal(p1[[paste0("run", i, ".log2FC")]], expected)
        expect_equal(p1[[paste0("run", i, ".ratio")]], 2^expected)
    }
})


# =============================================================================
# jensen_gap: why the mean of the logs is not the reported value
# =============================================================================

test_that("the pooled value is at or above the mean of the per-run log2FCs", {
    rec <- build_de_reconciliation_proteomics(recon_de_res(spread = 0.4), recon_config())

    expect_equal(rec$jensen_gap, rec$log2FC.imputs - rec$mean.log2FC.runs)
    # Jensen: log2(mean(2^x)) >= mean(x). Never the other way round.
    expect_true(all(rec$jensen_gap >= -1e-12))
    # With the runs genuinely disagreeing the gap is real, not a rounding
    # artefact -- which is the whole reason the column exists.
    expect_true(all(rec$jensen_gap > 1e-6))
})

test_that("the gap is exactly zero for a feature whose runs agree", {
    # p1's runs are identical, as they are for a feature measured in every
    # sample: nothing was imputed, so every draw saw the same matrix entries.
    # The other three disagree.
    cfg <- recon_config()
    rec <- build_de_reconciliation_proteomics(
        recon_de_res(spread = c(0, 0.4, 0.4, 0.4), config = cfg), cfg)

    agreed <- rec[rec$FeatureID == "p1", ]
    expect_equal(agreed$jensen_gap, 0, tolerance = 1e-12)
    expect_equal(agreed$mean.log2FC.runs, agreed$log2FC.imputs)

    # ...and the feature is not simply being skipped: the others still carry a
    # real gap in the same sheet.
    expect_true(all(rec$jensen_gap[rec$FeatureID != "p1"] > 1e-6))
})


# =============================================================================
# Shape: contrasts add rows, not run-column blocks
# =============================================================================

test_that("a second contrast adds rows and leaves the run columns alone", {
    cfg <- recon_config()
    one <- build_de_reconciliation_proteomics(
        recon_de_res(contrasts = "S_vs_NS", config = cfg), cfg)
    two <- build_de_reconciliation_proteomics(
        recon_de_res(contrasts = c("S_vs_NS", "S_vs_REF"), config = cfg), cfg)

    expect_equal(nrow(two), 2L * nrow(one))
    # The column set is identical: the run block is per row, not per contrast.
    expect_identical(names(two), names(one))
    expect_setequal(unique(two$Contrast), c("S_vs_NS", "S_vs_REF"))
    # Every feature appears once per contrast.
    expect_true(all(table(two$FeatureID) == 2L))
})

test_that("rows are feature-major, so one feature's contrasts sit together", {
    cfg <- recon_config()
    rec <- build_de_reconciliation_proteomics(
        recon_de_res(contrasts = c("S_vs_NS", "S_vs_REF"), config = cfg), cfg)

    # Not merely present in some order: adjacent. Someone checking one protein
    # by hand should not have to scroll past every other feature.
    first_seen <- match(unique(rec$FeatureID), rec$FeatureID)
    expect_equal(rec$FeatureID, rep(unique(rec$FeatureID), each = 2L))
    expect_equal(first_seen, seq(1L, nrow(rec), by = 2L))
})

test_that("the column order is the documented one", {
    rec <- build_de_reconciliation_proteomics(recon_de_res(n_runs = 2), recon_config(2))

    expect_identical(
        names(rec),
        c("FeatureID", "Contrast",
          "run1.log2FC", "run1.ratio", "run2.log2FC", "run2.ratio",
          "mean.ratio", "linearRatio.imputs", "delta.linearRatio",
          "log2FC.from_mean_ratio", "log2FC.imputs", "delta.log2FC",
          "mean.log2FC.runs", "jensen_gap")
    )
})


# =============================================================================
# Availability: only when two or more runs were really pooled
# =============================================================================

test_that("identical runs produce no sheet, however many of them there are", {
    # Several supported imputation methods return the same matrix for every
    # repetition: none, minval and DEP2 MinDet are deterministic, and both
    # impute_proteomics_qrilc() and impute_proteomics_dep2() call set.seed()
    # with a fixed configured seed inside each call, overwriting the per-run
    # seed. The runs are then not independent draws, and a sheet of identical
    # columns would read as agreement between draws that never differed.
    cfg <- recon_config()
    de_res <- recon_de_res(spread = 0, config = cfg)

    expect_message(
        res <- build_de_reconciliation_proteomics(de_res, cfg),
        "no pooling to reconcile"
    )
    expect_null(res)
})

test_that("a single run produces no sheet rather than a zero-delta one", {
    cfg <- recon_config(n_reps = 1, multi = FALSE)
    de_res <- recon_de_res(n_runs = 1, config = cfg)

    expect_null(build_de_reconciliation_proteomics(de_res, cfg))
})

test_that("precomputed DE produces no sheet", {
    # load_precomputed_proteomics_de() wraps the loaded tables as one pseudo-run
    # and carries no imputations at all, so there is no pooling to reconcile.
    cfg <- recon_config(n_reps = 1, multi = FALSE)
    de_res <- recon_de_res(n_runs = 1, config = cfg)
    de_res$imputations <- NULL

    expect_null(build_de_reconciliation_proteomics(de_res, cfg))
})

test_that("missing inputs return NULL rather than a malformed table", {
    cfg <- recon_config()
    expect_null(build_de_reconciliation_proteomics(list(), cfg))
    expect_null(build_de_reconciliation_proteomics(
        list(runs_de_tables = recon_runs(), summary_df = NULL), cfg))

    # Runs present, but the summary has no reported columns to check against.
    de_res <- recon_de_res()
    de_res$summary_df <- de_res$summary_df[, "FeatureID", drop = FALSE]
    expect_message(
        res <- build_de_reconciliation_proteomics(de_res, cfg),
        "no reported columns"
    )
    expect_null(res)
})


# =============================================================================
# The sheet in the workbook
# =============================================================================

test_that("the sheet note says this is reconciliation, not a fitted matrix", {
    skip_if_not_installed("openxlsx")
    rec <- build_de_reconciliation_proteomics(recon_de_res(), recon_config())

    wb <- openxlsx::createWorkbook()
    expect_true(add_de_reconciliation_sheet(wb, rec))
    expect_true("DE_reconciliation" %in% names(wb))

    note <- openxlsx::readWorkbook(wb, sheet = "DE_reconciliation",
                                   colNames = FALSE, rows = 1:8,
                                   skipEmptyRows = FALSE)
    note <- paste(unlist(note, use.names = FALSE), collapse = " ")

    expect_match(note, "independently imputed DE fit", fixed = TRUE)
    # Method-neutral: limma is not the only DE method mod_proteomics_de()
    # supports, so the note must not call every per-run estimate a limma
    # coefficient.
    expect_false(grepl("limma coefficient", note, fixed = TRUE))
    expect_match(note, "modes.proteomics.de.method", fixed = TRUE)
    expect_match(note, "written only when the runs genuinely differ", fixed = TRUE)
    expect_match(note, "2^run<N>.log2FC", fixed = TRUE)
    expect_match(note, "arithmetic mean", fixed = TRUE)
    expect_match(note, "log2(mean.ratio)", fixed = TRUE)
    expect_match(note, "same pooling rule", fixed = TRUE)
    # The two claims that matter most: nothing here was fitted, and the
    # existing column semantics are untouched.
    expect_match(note, "none of them was used to fit anything", fixed = TRUE)
    expect_match(note, "reconciliation checks, not results", fixed = TRUE)
    expect_match(note, "keep the matrix semantics they have always had", fixed = TRUE)
})

test_that("nothing is written when there is no reconciliation to show", {
    skip_if_not_installed("openxlsx")
    wb <- openxlsx::createWorkbook()

    expect_false(add_de_reconciliation_sheet(wb, NULL))
    expect_false(add_de_reconciliation_sheet(wb, data.frame()))
    expect_false("DE_reconciliation" %in% names(wb))
})

test_that("the How to read sheet points at the reconciliation only when pooled", {
    multi <- build_provenance_notes("proteomics", multi_imputation = TRUE)
    single <- build_provenance_notes("proteomics", multi_imputation = FALSE)

    expect_true(any(grepl("DE_reconciliation", multi$notes, fixed = TRUE)))
    expect_false(any(grepl("DE_reconciliation", single$notes, fixed = TRUE)))
    # Worded for the case the notes cannot see: multi_imputation: true with
    # no_repetitions: 1 pools nothing, so the sheet is absent although the flag
    # is on. An unconditional sentence would point at a missing sheet.
    expect_true(any(grepl("two or more runs were pooled and their estimates actually differ",
                          multi$notes, fixed = TRUE)))
})
