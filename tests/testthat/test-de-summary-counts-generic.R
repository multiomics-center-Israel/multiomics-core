# tests/testthat/test-de-summary-counts-generic.R
#
# Unit tests for build_de_summary_counts_generic() and the three domain wrappers.

# ---------------------------------------------------------------------------
# Helper: build a minimal de_stats data.frame
# ---------------------------------------------------------------------------
# Generic test helper — uses "pass." and "linearFC." dot-prefix naming
# for testing build_de_summary_counts_generic() with custom patterns.
make_de_stats_generic <- function(include_any = TRUE) {
    df <- data.frame(
        feature_id          = paste0("F", 1:10),
        pass.A_vs_B         = c(1, 1, 1, NA, 0, 1, NA, 0, 1, 1),
        pass.C_vs_D         = c(1, NA, 0, 1, 1, 0, NA, 1, 0, 0),
        linearFC.A_vs_B     = c(0.5, -0.3, 1.2, 0.1, -0.8, 0.4, 0.2, -1.0, -0.6, 0.9),
        linearFC.C_vs_D     = c(-0.2, 0.7, 0.3, -0.5, 0.1, -0.4, 0.6, 0.8, -0.3, 0.2),
        stringsAsFactors    = FALSE
    )
    if (include_any) {
        # features significant in at least one contrast
        df$pass_any_contrast <- ifelse(
            (!is.na(df$pass.A_vs_B) & df$pass.A_vs_B == 1) |
            (!is.na(df$pass.C_vs_D) & df$pass.C_vs_D == 1),
            1, NA
        )
    }
    df
}

# Metabolomics wrapper helper — uses "pass_" and "logFC_" underscore naming
# matching the current build_de_summary_counts_metabolomics() expectations.
make_de_stats_metab <- function(include_any = TRUE) {
    df <- data.frame(
        feature_id          = paste0("F", 1:10),
        pass_A_vs_B         = c(1, 1, 1, NA, 0, 1, NA, 0, 1, 1),
        pass_C_vs_D         = c(1, NA, 0, 1, 1, 0, NA, 1, 0, 0),
        logFC_A_vs_B        = c(0.5, -0.3, 1.2, 0.1, -0.8, 0.4, 0.2, -1.0, -0.6, 0.9),
        logFC_C_vs_D        = c(-0.2, 0.7, 0.3, -0.5, 0.1, -0.4, 0.6, 0.8, -0.3, 0.2),
        stringsAsFactors    = FALSE
    )
    if (include_any) {
        # features significant in at least one contrast
        df$pass_any_contrast <- ifelse(
            (!is.na(df$pass_A_vs_B) & df$pass_A_vs_B == 1) |
            (!is.na(df$pass_C_vs_D) & df$pass_C_vs_D == 1),
            1, NA
        )
    }
    df
}

make_de_stats_prot <- function() {
    # Proteomics schema: pass.imputs.<contrast>, linearFC.imputs.<contrast>
    data.frame(
        feature_id                = paste0("P", 1:6),
        pass.imputs.X_vs_Y       = c(1, 1, 0, 1, NA, 1),
        linearFC.imputs.X_vs_Y   = c(0.5, -0.3, 0.1, -0.8, 0.4, 0.2),
        pass_any_contrast         = c(1, 1, NA, 1, NA, 1),
        stringsAsFactors          = FALSE
    )
}

make_de_stats_prot_logfc <- function() {
    # Proteomics schema with logFC_ candidate (higher priority than linearFC.imputs.)
    data.frame(
        feature_id                = paste0("P", 1:4),
        pass.imputs.A_vs_B       = c(1, 1, 0, 1),
        logFC_A_vs_B             = c(0.5, -0.3, 0.1, -0.8),
        linearFC.imputs.A_vs_B   = c(0.9, -0.9, 0.1, -0.1),
        pass_any_contrast         = c(1, 1, NA, 1),
        stringsAsFactors          = FALSE
    )
}

make_de_stats_rna <- function() {
    data.frame(
        feature_id          = paste0("G", 1:8),
        A_vs_B_pass         = c(1, 1, 1, NA, 0, 1, 0, 1),
        linearFC.A_vs_B     = c(0.5, -0.3, 1.2, 0.1, -0.8, 0.4, -0.6, 0.9),
        pass_any_contrast   = c(1, 1, 1, NA, NA, 1, NA, 1),
        stringsAsFactors    = FALSE
    )
}


# ===========================================================================
# Tests for build_de_summary_counts_generic()
# ===========================================================================

test_that("generic: NULL input returns NULL", {
    result <- build_de_summary_counts_generic(
        NULL, "^pass\\.", function(x) x, function(cn, cols) NULL
    )
    expect_null(result)
})

test_that("generic: no matching pass columns returns NULL", {
    df <- data.frame(feature_id = 1:3, foo = 1:3)
    result <- build_de_summary_counts_generic(
        df, "^pass\\.", function(x) x, function(cn, cols) NULL
    )
    expect_null(result)
})

test_that("generic: correct per-contrast counts with FC column", {
    df <- make_de_stats_generic(include_any = FALSE)
    result <- build_de_summary_counts_generic(
        de_stats         = df,
        pass_pattern     = "^pass\\.",
        extract_contrast = function(col) sub("^pass\\.", "", col),
        find_fc_col      = function(cn, cols) {
            fc <- paste0("linearFC.", cn)
            if (fc %in% cols) fc else NULL
        }
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
        df, "^pass\\.",
        function(col) sub("^pass\\.", "", col),
        function(cn, cols) NULL
    )
    expect_equal(result$total, 2)
    expect_true(is.na(result$up))
    expect_true(is.na(result$down))
})

test_that("generic: pass_any_contrast excluded from per-contrast rows", {
    df <- data.frame(
        feature_id        = 1:3,
        pass.A            = c(1, 0, 1),
        pass_any_contrast = c(1, 0, 1),
        linearFC.A        = c(0.5, -0.3, -0.1),
        stringsAsFactors  = FALSE
    )
    result <- build_de_summary_counts_generic(
        df, "^pass\\.",
        function(col) sub("^pass\\.", "", col),
        function(cn, cols) {
            fc <- paste0("linearFC.", cn)
            if (fc %in% cols) fc else NULL
        }
    )
    # Should have per-contrast row for "A" + "any" row, nothing for "any_contrast"
    expect_equal(result$contrast, c("A", "any"))
})

test_that("generic: 'any' row appended when pass_any_contrast exists", {
    df <- make_de_stats_generic(include_any = TRUE)
    result <- build_de_summary_counts_generic(
        df, "^pass\\.",
        function(col) sub("^pass\\.", "", col),
        function(cn, cols) {
            fc <- paste0("linearFC.", cn)
            if (fc %in% cols) fc else NULL
        }
    )
    any_row <- result[result$contrast == "any", ]
    expect_equal(nrow(any_row), 1)
    expect_equal(any_row$up, 0)
    expect_equal(any_row$down, 0)

    # Verify total = number of features with pass_any_contrast == 1
    expected_any <- sum(!is.na(df$pass_any_contrast) & df$pass_any_contrast == 1)
    expect_equal(any_row$total, expected_any)

    # "any" is the last row
    expect_equal(result$contrast[nrow(result)], "any")
})

test_that("generic: no 'any' row when pass_any_contrast is absent", {
    df <- make_de_stats_generic(include_any = FALSE)
    result <- build_de_summary_counts_generic(
        df, "^pass\\.",
        function(col) sub("^pass\\.", "", col),
        function(cn, cols) {
            fc <- paste0("linearFC.", cn)
            if (fc %in% cols) fc else NULL
        }
    )
    expect_false("any" %in% result$contrast)
})

test_that("generic: custom is_significant callback is honoured", {
    df <- data.frame(
        feature_id = 1:5,
        pass.X     = c(TRUE, TRUE, FALSE, NA, TRUE),
        linearFC.X = c(0.5, -0.3, 0.1, 0.2, -0.8),
        stringsAsFactors = FALSE
    )
    # Use a boolean-based significance test instead of == 1
    result <- build_de_summary_counts_generic(
        de_stats         = df,
        pass_pattern     = "^pass\\.",
        extract_contrast = function(col) sub("^pass\\.", "", col),
        find_fc_col      = function(cn, cols) {
            fc <- paste0("linearFC.", cn)
            if (fc %in% cols) fc else NULL
        },
        is_significant   = function(x) !is.na(x) & x == TRUE
    )
    # TRUE at rows 1, 2, 5 => total = 3
    expect_equal(result$total[result$contrast == "X"], 3)
    # FC > 0: row 1(0.5) => up = 1; FC < 0: rows 2(-0.3), 5(-0.8) => down = 2
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
        df, "^pass\\.",
        function(col) sub("^pass\\.", "", col),
        function(cn, cols) {
            fc <- paste0("linearFC.", cn)
            if (fc %in% cols) fc else NULL
        }
        # is_significant omitted — uses default
    )
    # == 1 at rows 1, 4 => total = 2
    expect_equal(result$total[result$contrast == "X"], 2)
})


# ===========================================================================
# Tests for domain wrappers
# ===========================================================================

test_that("metabolomics wrapper: correct output with pass_any_contrast excluded", {
    df <- make_de_stats_metab(include_any = TRUE)
    result <- build_de_summary_counts_metabolomics(df)

    expect_equal(names(result), c("contrast", "up", "down", "total"))
    expect_true("A_vs_B" %in% result$contrast)
    expect_true("C_vs_D" %in% result$contrast)
    # build_de_summary_counts_metabolomics excludes pass_any_contrast
    # and does not append an "any" row
    expect_false("any" %in% result$contrast)
    expect_false("any_contrast" %in% result$contrast)
})

test_that("metabolomics wrapper: NULL input returns NULL", {
    expect_null(build_de_summary_counts_metabolomics(NULL))
})

test_that("proteomics wrapper: pass.imputs. schema with correct contrast names", {
    df <- make_de_stats_prot()
    result <- build_de_summary_counts_proteomics(df)

    expect_equal(names(result), c("contrast", "up", "down", "total"))
    # Contrast name should be "X_vs_Y", NOT "imputs.X_vs_Y"
    row_x <- result[result$contrast == "X_vs_Y", ]
    expect_equal(nrow(row_x), 1)
    # pass.imputs.X_vs_Y == 1 at rows 1,2,4,6 => total = 4
    expect_equal(row_x$total, 4)
    # linearFC.imputs.X_vs_Y > 0 among sig: rows 1(0.5), 6(0.2) => up = 2
    expect_equal(row_x$up, 2)
    # linearFC.imputs.X_vs_Y < 0 among sig: rows 2(-0.3), 4(-0.8) => down = 2
    expect_equal(row_x$down, 2)
    # "any" row
    expect_equal(result$contrast[nrow(result)], "any")
    expect_equal(result[result$contrast == "any", "total"], 4)
})

test_that("proteomics wrapper: FC candidate priority (logFC_ before linearFC.imputs.)", {
    df <- make_de_stats_prot_logfc()
    result <- build_de_summary_counts_proteomics(df)

    row <- result[result$contrast == "A_vs_B", ]
    # pass.imputs.A_vs_B == 1 at rows 1,2,4 => total = 3
    expect_equal(row$total, 3)
    # Should use logFC_A_vs_B (higher priority), not linearFC.imputs.A_vs_B
    # logFC_A_vs_B > 0 among sig: row 1(0.5) => up = 1
    expect_equal(row$up, 1)
    # logFC_A_vs_B < 0 among sig: rows 2(-0.3), 4(-0.8) => down = 2
    expect_equal(row$down, 2)
})

test_that("proteomics wrapper: NULL input returns NULL", {
    expect_null(build_de_summary_counts_proteomics(NULL))
})

test_that("rnaseq wrapper: _pass$ suffix handling (anchored)", {
    df <- make_de_stats_rna()
    result <- build_de_summary_counts_rnaseq(df)

    expect_equal(names(result), c("contrast", "up", "down", "total"))
    row_a <- result[result$contrast == "A_vs_B", ]
    # A_vs_B_pass == 1 at rows 1,2,3,6,8 => total = 5
    expect_equal(row_a$total, 5)
    # FC > 0 among sig: rows 1(0.5), 3(1.2), 6(0.4), 8(0.9) => up = 4
    expect_equal(row_a$up, 4)
    # FC < 0 among sig: row 2(-0.3) => down = 1
    expect_equal(row_a$down, 1)
    # "any" row
    expect_equal(result$contrast[nrow(result)], "any")
    expect_equal(result[result$contrast == "any", "total"], 5)
})

test_that("rnaseq wrapper: _pass$ does not match mid-string _pass_", {
    # Ensure a column like "foo_pass_bar" is NOT picked up
    df <- data.frame(
        feature_id          = 1:3,
        A_vs_B_pass         = c(1, 0, 1),
        foo_pass_bar        = c(1, 1, 1),
        linearFC.A_vs_B     = c(0.5, -0.3, 0.8),
        stringsAsFactors    = FALSE
    )
    result <- build_de_summary_counts_rnaseq(df)
    # Only A_vs_B should appear, not "foo" or "foo_pass_bar"
    per_contrast <- result[result$contrast != "any", ]
    expect_equal(nrow(per_contrast), 1)
    expect_equal(per_contrast$contrast, "A_vs_B")
})

test_that("rnaseq wrapper: grep fallback for FC column", {
    df <- data.frame(
        feature_id             = paste0("G", 1:4),
        contrast1_pass         = c(1, 1, 0, 1),
        linearFC.alt.contrast1 = c(0.5, -0.3, 0.1, 0.8),
        pass_any_contrast      = c(1, 1, NA, 1),
        stringsAsFactors       = FALSE
    )
    result <- build_de_summary_counts_rnaseq(df)
    row <- result[result$contrast == "contrast1", ]
    expect_equal(row$total, 3)
    # grep fallback finds linearFC.alt.contrast1 => up/down should work
    expect_equal(row$up, 2)
    expect_equal(row$down, 1)
})

test_that("rnaseq wrapper: NULL input returns NULL", {
    expect_null(build_de_summary_counts_rnaseq(NULL))
})
