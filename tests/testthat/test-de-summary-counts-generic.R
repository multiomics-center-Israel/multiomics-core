# tests/testthat/test-de-summary-counts-generic.R
#
# Unit tests for build_de_summary_counts_generic() and the three domain
# wrappers.  All tests use the explicit-contrasts signature — no regex over
# names(de_stats) — so user-supplied annotation columns merged from row_data
# cannot corrupt the counts.

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

make_de_stats_metab <- function(include_any = TRUE) {
    df <- data.frame(
        feature_id       = paste0("F", 1:10),
        pass.A_vs_B      = c(1, 1, 1, NA, 0, 1, NA, 0, 1, 1),
        pass.C_vs_D      = c(1, NA, 0, 1, 1, 0, NA, 1, 0, 0),
        linearFC.A_vs_B  = c(0.5, -0.3, 1.2, 0.1, -0.8, 0.4, 0.2, -1.0, -0.6, 0.9),
        linearFC.C_vs_D  = c(-0.2, 0.7, 0.3, -0.5, 0.1, -0.4, 0.6, 0.8, -0.3, 0.2),
        stringsAsFactors = FALSE
    )
    if (include_any) {
        df$pass_any_contrast <- ifelse(
            (!is.na(df$pass.A_vs_B) & df$pass.A_vs_B == 1) |
            (!is.na(df$pass.C_vs_D) & df$pass.C_vs_D == 1),
            1, NA
        )
    }
    df
}

make_de_stats_prot <- function() {
    data.frame(
        feature_id             = paste0("P", 1:6),
        pass.imputs.X_vs_Y     = c(1, 1, 0, 1, NA, 1),
        linearFC.imputs.X_vs_Y = c(0.5, -0.3, 0.1, -0.8, 0.4, 0.2),
        pass_any_contrast      = c(1, 1, NA, 1, NA, 1),
        stringsAsFactors       = FALSE
    )
}

make_de_stats_prot_logfc <- function() {
    data.frame(
        feature_id             = paste0("P", 1:4),
        pass.imputs.A_vs_B     = c(1, 1, 0, 1),
        logFC_A_vs_B           = c(0.5, -0.3, 0.1, -0.8),
        linearFC.imputs.A_vs_B = c(0.9, -0.9, 0.1, -0.1),
        pass_any_contrast      = c(1, 1, NA, 1),
        stringsAsFactors       = FALSE
    )
}

make_de_stats_rna <- function() {
    data.frame(
        feature_id        = paste0("G", 1:8),
        A_vs_B_pass       = c(1, 1, 1, NA, 0, 1, 0, 1),
        linearFC.A_vs_B   = c(0.5, -0.3, 1.2, 0.1, -0.8, 0.4, -0.6, 0.9),
        pass_any_contrast = c(1, 1, 1, NA, NA, 1, NA, 1),
        stringsAsFactors  = FALSE
    )
}

# Callback factories so each test reads cleanly.
metab_pass <- function(cn) paste0("pass.", cn)
metab_fc   <- function(cn, cols) {
    fc <- paste0("linearFC.", cn)
    if (fc %in% cols) fc else NULL
}
rna_pass <- function(cn) paste0(cn, "_pass")
rna_fc   <- metab_fc

# ===========================================================================
# Tests for build_de_summary_counts_generic()
# ===========================================================================

test_that("generic: NULL input returns NULL", {
    result <- build_de_summary_counts_generic(
        NULL, contrasts = "A_vs_B",
        pass_col_for = metab_pass, fc_col_for = metab_fc
    )
    expect_null(result)
})

test_that("generic: empty contrasts returns NULL", {
    df <- make_de_stats_metab(include_any = FALSE)
    expect_null(build_de_summary_counts_generic(
        df, contrasts = character(0),
        pass_col_for = metab_pass, fc_col_for = metab_fc
    ))
    expect_null(build_de_summary_counts_generic(
        df, contrasts = NULL,
        pass_col_for = metab_pass, fc_col_for = metab_fc
    ))
})

test_that("generic: no matching pass columns returns NULL", {
    df <- data.frame(feature_id = 1:3, foo = 1:3)
    result <- build_de_summary_counts_generic(
        df, contrasts = "A_vs_B",
        pass_col_for = metab_pass, fc_col_for = metab_fc
    )
    expect_null(result)
})

test_that("generic: correct per-contrast counts with FC column", {
    df <- make_de_stats_metab(include_any = FALSE)
    result <- build_de_summary_counts_generic(
        df, contrasts = c("A_vs_B", "C_vs_D"),
        pass_col_for = metab_pass, fc_col_for = metab_fc
    )

    expect_equal(names(result), c("contrast", "up", "down", "total"))
    expect_true("A_vs_B" %in% result$contrast)
    expect_true("C_vs_D" %in% result$contrast)

    row_a <- result[result$contrast == "A_vs_B", ]
    # pass.A_vs_B == 1 at rows 1,2,3,6,9,10 => total = 6
    expect_equal(row_a$total, 6)
    # FC > 0 among sig: rows 1(0.5), 3(1.2), 6(0.4), 10(0.9) => up = 4
    expect_equal(row_a$up, 4)
    # FC < 0 among sig: rows 2(-0.3), 9(-0.6) => down = 2
    expect_equal(row_a$down, 2)
})

test_that("generic: up/down are NA when FC column is missing", {
    df <- data.frame(
        feature_id = 1:3,
        pass.X     = c(1, 1, 0),
        stringsAsFactors = FALSE
    )
    result <- build_de_summary_counts_generic(
        df, contrasts = "X",
        pass_col_for = metab_pass,
        fc_col_for   = function(cn, cols) NULL
    )
    expect_equal(result$total, 2)
    expect_true(is.na(result$up))
    expect_true(is.na(result$down))
})

test_that("generic: 'any' row appended when pass_any_contrast exists", {
    df <- make_de_stats_metab(include_any = TRUE)
    result <- build_de_summary_counts_generic(
        df, contrasts = c("A_vs_B", "C_vs_D"),
        pass_col_for = metab_pass, fc_col_for = metab_fc
    )
    any_row <- result[result$contrast == "any", ]
    expect_equal(nrow(any_row), 1)
    expect_equal(any_row$up, 0)
    expect_equal(any_row$down, 0)
    expected_any <- sum(!is.na(df$pass_any_contrast) & df$pass_any_contrast == 1)
    expect_equal(any_row$total, expected_any)
    expect_equal(result$contrast[nrow(result)], "any")
})

test_that("generic: no 'any' row when pass_any_contrast is absent", {
    df <- make_de_stats_metab(include_any = FALSE)
    result <- build_de_summary_counts_generic(
        df, contrasts = c("A_vs_B", "C_vs_D"),
        pass_col_for = metab_pass, fc_col_for = metab_fc
    )
    expect_false("any" %in% result$contrast)
})

test_that("generic: contrast with missing pass column is silently skipped", {
    df <- data.frame(
        feature_id      = 1:3,
        pass.A          = c(1, 1, 0),
        linearFC.A      = c(0.5, -0.3, 0.1),
        stringsAsFactors = FALSE
    )
    # "B" has no pass column — should be skipped, not error.
    result <- build_de_summary_counts_generic(
        df, contrasts = c("A", "B"),
        pass_col_for = metab_pass, fc_col_for = metab_fc
    )
    expect_equal(result$contrast, "A")
})

test_that("generic: custom is_significant callback is honoured", {
    df <- data.frame(
        feature_id = 1:5,
        pass.X     = c(TRUE, TRUE, FALSE, NA, TRUE),
        linearFC.X = c(0.5, -0.3, 0.1, 0.2, -0.8),
        stringsAsFactors = FALSE
    )
    result <- build_de_summary_counts_generic(
        df, contrasts = "X",
        pass_col_for  = metab_pass,
        fc_col_for    = metab_fc,
        is_significant = function(x) !is.na(x) & x == TRUE
    )
    expect_equal(result$total[result$contrast == "X"], 3)
    expect_equal(result$up[result$contrast == "X"], 1)
    expect_equal(result$down[result$contrast == "X"], 2)
})

test_that("generic: is_significant default matches == 1 semantics", {
    df <- data.frame(
        feature_id = 1:4,
        pass.X     = c(1, 0, NA, 1),
        linearFC.X = c(0.5, -0.3, 0.1, -0.8),
        stringsAsFactors = FALSE
    )
    result <- build_de_summary_counts_generic(
        df, contrasts = "X",
        pass_col_for = metab_pass, fc_col_for = metab_fc
    )
    expect_equal(result$total[result$contrast == "X"], 2)
})


# ===========================================================================
# Tests for domain wrappers
# ===========================================================================

test_that("metabolomics wrapper: correct output with pass_any_contrast", {
    df <- make_de_stats_metab(include_any = TRUE)
    result <- build_de_summary_counts_metabolomics(
        df, contrasts = c("A_vs_B", "C_vs_D")
    )

    expect_equal(names(result), c("contrast", "up", "down", "total"))
    expect_true("A_vs_B" %in% result$contrast)
    expect_true("C_vs_D" %in% result$contrast)
    expect_true("any" %in% result$contrast)

    row_a <- result[result$contrast == "A_vs_B", ]
    expect_equal(row_a$total, 6)
    expect_equal(row_a$up, 4)
    expect_equal(row_a$down, 2)
    expect_equal(result$contrast[nrow(result)], "any")
})

test_that("metabolomics wrapper: NULL input returns NULL", {
    expect_null(build_de_summary_counts_metabolomics(NULL, contrasts = "A_vs_B"))
})

test_that("proteomics wrapper: pass.imputs. schema with correct contrast names", {
    df <- make_de_stats_prot()
    result <- build_de_summary_counts_proteomics(df, contrasts = "X_vs_Y")

    expect_equal(names(result), c("contrast", "up", "down", "total"))
    row_x <- result[result$contrast == "X_vs_Y", ]
    expect_equal(nrow(row_x), 1)
    # pass.imputs.X_vs_Y == 1 at rows 1,2,4,6 => total = 4
    expect_equal(row_x$total, 4)
    # linearFC.imputs.X_vs_Y > 0 among sig: rows 1(0.5), 6(0.2) => up = 2
    expect_equal(row_x$up, 2)
    # linearFC.imputs.X_vs_Y < 0 among sig: rows 2(-0.3), 4(-0.8) => down = 2
    expect_equal(row_x$down, 2)
    expect_equal(result$contrast[nrow(result)], "any")
    expect_equal(result[result$contrast == "any", "total"], 4)
})

test_that("proteomics wrapper: FC candidate priority (logFC_ before linearFC.imputs.)", {
    df <- make_de_stats_prot_logfc()
    result <- build_de_summary_counts_proteomics(df, contrasts = "A_vs_B")

    row <- result[result$contrast == "A_vs_B", ]
    expect_equal(row$total, 3)
    # Should use logFC_A_vs_B (higher priority), not linearFC.imputs.A_vs_B
    expect_equal(row$up, 1)
    expect_equal(row$down, 2)
})

test_that("proteomics wrapper: NULL input returns NULL", {
    expect_null(build_de_summary_counts_proteomics(NULL, contrasts = "X_vs_Y"))
})

test_that("rnaseq wrapper: <contrast>_pass suffix handling", {
    df <- make_de_stats_rna()
    result <- build_de_summary_counts_rnaseq(df, contrasts = "A_vs_B")

    expect_equal(names(result), c("contrast", "up", "down", "total"))
    row_a <- result[result$contrast == "A_vs_B", ]
    expect_equal(row_a$total, 5)
    expect_equal(row_a$up, 4)
    expect_equal(row_a$down, 1)
    expect_equal(result$contrast[nrow(result)], "any")
    expect_equal(result[result$contrast == "any", "total"], 5)
})

test_that("rnaseq wrapper: mid-string '_pass_' columns are ignored", {
    # A row_data-style 'foo_pass_bar' column must never be counted, even if
    # the regex-era code would have matched _pass$. With explicit contrasts,
    # only paste0(cn, '_pass') is consulted.
    df <- data.frame(
        feature_id       = 1:3,
        A_vs_B_pass      = c(1, 0, 1),
        foo_pass_bar     = c(1, 1, 1),
        linearFC.A_vs_B  = c(0.5, -0.3, 0.8),
        stringsAsFactors = FALSE
    )
    result <- build_de_summary_counts_rnaseq(df, contrasts = "A_vs_B")
    per_contrast <- result[result$contrast != "any", ]
    expect_equal(nrow(per_contrast), 1)
    expect_equal(per_contrast$contrast, "A_vs_B")
})

test_that("rnaseq wrapper: non-exact FC column name yields NA up/down (no grep fallback)", {
    # The old wrapper grepped `linearFC.*<cn>` and would pick up
    # `linearFC.alt.contrast1`. The new wrapper requires exact match
    # `linearFC.<cn>` and reports NA when it's missing.
    df <- data.frame(
        feature_id             = paste0("G", 1:4),
        contrast1_pass         = c(1, 1, 0, 1),
        linearFC.alt.contrast1 = c(0.5, -0.3, 0.1, 0.8),
        pass_any_contrast      = c(1, 1, NA, 1),
        stringsAsFactors       = FALSE
    )
    result <- build_de_summary_counts_rnaseq(df, contrasts = "contrast1")
    row <- result[result$contrast == "contrast1", ]
    expect_equal(row$total, 3)
    expect_true(is.na(row$up))
    expect_true(is.na(row$down))
})

test_that("rnaseq wrapper: NULL input returns NULL", {
    expect_null(build_de_summary_counts_rnaseq(NULL, contrasts = "A_vs_B"))
})
