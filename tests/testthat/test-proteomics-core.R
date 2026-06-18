# Tests for proteomics core DE summary functions
# test-proteomics-core.R

# --- Mock data helpers ---

create_mock_de_table <- function(n = 20) {
    data.frame(
        FeatureID = paste0("PROT", seq_len(n)),
        logFC = rnorm(n, 0, 2),
        adj.P.Val = runif(n, 0, 0.1),
        P.Value = runif(n, 0, 0.05),
        t = rnorm(n, 0, 3),
        stringsAsFactors = FALSE
    )
}

create_mock_summary <- function(n = 50) {
    df <- data.frame(
        FeatureID = paste0("PROT", seq_len(n)),
        pass.imputs.A_vs_B = sample(c(1, NA), n, replace = TRUE),
        pass.imputs.C_vs_D = sample(c(1, NA), n, replace = TRUE),
        linearFC.imputs.A_vs_B = runif(n, 0.5, 3),
        linearFC.imputs.C_vs_D = runif(n, 0.5, 3),
        padj.imputs.A_vs_B = runif(n, 0, 0.2),
        padj.imputs.C_vs_D = runif(n, 0, 0.2),
        stringsAsFactors = FALSE
    )
    df
}

# --- Tests for mark_pass1 ---

test_that("mark_pass1 correctly marks passing features", {
    de_table <- data.frame(
        adj.P.Val = c(0.01, 0.1, 0.04, NA, 0.001),
        logFC = c(2, 0.1, -1.5, 1, 3),
        P.Value = c(0.001, 0.05, 0.02, NA, 0.0001),
        stringsAsFactors = FALSE
    )
    result <- mark_pass1(de_table, p_cutoff = 0.05, lfc_cutoff = log2(1.5), use_adj = TRUE)
    # Feature 1: padj=0.01<0.05, |logFC|=2>0.585 -> pass (1)
    expect_equal(result[1], 1)
    # Feature 2: padj=0.1>0.05 -> fail (NA)
    expect_true(is.na(result[2]))
    # Feature 3: padj=0.04<0.05, |logFC|=1.5>0.585 -> pass (1)
    expect_equal(result[3], 1)
    # Feature 4: padj=NA -> fail (NA)
    expect_true(is.na(result[4]))
    # Feature 5: padj=0.001<0.05, |logFC|=3>0.585 -> pass (1)
    expect_equal(result[5], 1)
})

test_that("mark_pass1 uses P.Value when use_adj is FALSE", {
    de_table <- data.frame(
        adj.P.Val = c(0.1, 0.2),
        logFC = c(2, 2),
        P.Value = c(0.01, 0.01),
        stringsAsFactors = FALSE
    )
    result_adj <- mark_pass1(de_table, p_cutoff = 0.05, lfc_cutoff = 0, use_adj = TRUE)
    result_raw <- mark_pass1(de_table, p_cutoff = 0.05, lfc_cutoff = 0, use_adj = FALSE)
    # With adj: both fail (adj.P.Val > 0.05)
    expect_true(all(is.na(result_adj)))
    # Without adj: both pass (P.Value < 0.05)
    expect_true(all(result_raw == 1))
})

test_that("mark_pass1 returns all NA when no features pass", {
    de_table <- data.frame(
        adj.P.Val = c(0.5, 0.8, 0.9),
        logFC = c(0.1, 0.05, 0.01),
        P.Value = c(0.3, 0.5, 0.7),
        stringsAsFactors = FALSE
    )
    result <- mark_pass1(de_table, p_cutoff = 0.05, lfc_cutoff = log2(1.5), use_adj = TRUE)
    expect_true(all(is.na(result)))
})

test_that("mark_pass1 requires logFC and p-value columns", {
    de_table <- data.frame(
        some_col = c(1, 2, 3),
        stringsAsFactors = FALSE
    )
    expect_error(mark_pass1(de_table, p_cutoff = 0.05, lfc_cutoff = 0, use_adj = TRUE))
})

# --- Tests for add_pass_any_contrast ---

test_that("add_pass_any_contrast adds correct columns", {
    df <- data.frame(
        FeatureID = c("A", "B", "C"),
        pass.imputs.X = c(1, NA, 1),
        pass.imputs.Y = c(NA, NA, 1),
        stringsAsFactors = FALSE
    )
    result <- add_pass_any_contrast(df)
    expect_true("pass_any_contrast" %in% names(result))
    expect_true("n_pass_contrasts" %in% names(result))
    # A passes in 1 contrast
    expect_equal(result$pass_any_contrast[1], 1)
    expect_equal(result$n_pass_contrasts[1], 1L)
    # B passes in 0 contrasts
    expect_true(is.na(result$pass_any_contrast[2]))
    expect_equal(result$n_pass_contrasts[2], 0L)
    # C passes in 2 contrasts
    expect_equal(result$pass_any_contrast[3], 1)
    expect_equal(result$n_pass_contrasts[3], 2L)
})

test_that("add_pass_any_contrast handles all-NA edge case", {
    df <- data.frame(
        FeatureID = c("A", "B"),
        pass.imputs.X = c(NA, NA),
        stringsAsFactors = FALSE
    )
    result <- add_pass_any_contrast(df)
    expect_true(all(is.na(result$pass_any_contrast)))
    expect_true(all(result$n_pass_contrasts == 0L))
})

test_that("add_pass_any_contrast handles no pass columns gracefully", {
    df <- data.frame(
        FeatureID = c("A", "B"),
        some_other_col = c(1, 2),
        stringsAsFactors = FALSE
    )
    result <- add_pass_any_contrast(df)
    expect_true("pass_any_contrast" %in% names(result))
    expect_true("n_pass_contrasts" %in% names(result))
    expect_true(all(is.na(result$pass_any_contrast)))
    expect_true(all(result$n_pass_contrasts == 0L))
})

# --- Tests for build_de_contrast_summary ---
# Note: 07_shiny_export.R also defines build_de_contrast_summary with different
# signature. Re-source 05_de_summary.R to test the correct version.
local({
    root <- if (dir.exists("R")) "." else "../.."
    source(file.path(root, "R", "domain", "proteomics", "05_de_summary.R"), local = TRUE)
    build_de_contrast_summary_05 <<- build_de_contrast_summary
})

test_that("build_de_contrast_summary computes correct counts", {
    df <- data.frame(
        FeatureID = paste0("P", 1:5),
        pass.imputs.A_vs_B = c(1, 1, 1, NA, NA),
        linearFC.imputs.A_vs_B = c(2, -0.4, 3, 1, 1),
        padj.imputs.A_vs_B = c(0.01, 0.02, 0.03, 0.5, 0.8),
        stringsAsFactors = FALSE
    )
    result <- build_de_contrast_summary_05(df)
    expect_true(is.data.frame(result))
    expect_true("Contrast" %in% names(result))
    expect_equal(nrow(result), 1)
    expect_equal(result$Contrast[1], "A_vs_B")
    expect_equal(result$total_genes[1], 5)
    expect_equal(result$n_significant[1], 3)
    # P1: pass, FC=2>0 -> up; P2: pass, FC=-0.4<0 -> down;
    # P3: pass, FC=3>0 -> up
    expect_equal(result$n_up[1], 2)
    expect_equal(result$n_down[1], 1)
})

test_that("build_de_contrast_summary returns NULL for empty input", {
    result <- build_de_contrast_summary_05(data.frame())
    expect_null(result)
})

test_that("build_de_contrast_summary returns NULL for NULL input", {
    result <- build_de_contrast_summary_05(NULL)
    expect_null(result)
})

test_that("build_de_contrast_summary handles multiple contrasts", {
    df <- data.frame(
        FeatureID = paste0("P", 1:3),
        pass.imputs.X_vs_Y = c(1, NA, 1),
        pass.imputs.A_vs_B = c(NA, 1, NA),
        linearFC.imputs.X_vs_Y = c(2, 1, -3),
        linearFC.imputs.A_vs_B = c(1, 2, 1),
        stringsAsFactors = FALSE
    )
    result <- build_de_contrast_summary_05(df)
    expect_equal(nrow(result), 2)
    expect_true(all(c("X_vs_Y", "A_vs_B") %in% result$Contrast))
})

# --- Tests for get_de_features ---

test_that("get_de_features returns correct IDs", {
    de_res <- list(
        summary_df = data.frame(
            FeatureID = c("P1", "P2", "P3", "P4"),
            pass_any_contrast = c(1, NA, 1, NA),
            stringsAsFactors = FALSE
        )
    )
    cfg <- list(
        de_table = list(
            id_col = "FeatureID",
            pass_any_col = "pass_any_contrast"
        )
    )
    result <- get_de_features(de_res, cfg)
    expect_equal(sort(result), c("P1", "P3"))
})

test_that("get_de_features returns empty for no DE features", {
    de_res <- list(
        summary_df = data.frame(
            FeatureID = c("P1", "P2"),
            pass_any_contrast = c(NA, NA),
            stringsAsFactors = FALSE
        )
    )
    cfg <- list(
        de_table = list(
            id_col = "FeatureID",
            pass_any_col = "pass_any_contrast"
        )
    )
    result <- get_de_features(de_res, cfg)
    expect_equal(length(result), 0)
})

test_that("get_de_features errors on missing summary_df", {
    de_res <- list(other_field = "value")
    cfg <- list(de_table = list(id_col = "FeatureID", pass_any_col = "pass_any_contrast"))
    expect_error(get_de_features(de_res, cfg))
})
