# tests/testthat/test-report-log2fc-rounding.R
#
# The proteomics report reconstructed log2 fold changes from the linearFC
# column, which is written with signif(x, 3). Anything in [1.4950, 1.5049] is
# stored as 1.50, and log2(1.50) = 0.5849625 is EXACTLY the log2(1.5) threshold,
# so features genuinely below 1.5-fold were counted as passing.
#
# On a real run this inflated the interactive volcano from 222 to 227:
# PLS3, SERINC3, BLTP3A, CKAP2 and ISYNA1, whose true ratios were 1.4957 to
# 1.4996. The tables, which use log2FC directly, reported 222 throughout.

signed_fc_to_log2 <- function(fc) {
    fc <- as.numeric(fc)
    ifelse(is.na(fc) | fc == 0, NA_real_, log2(abs(fc)) * sign(fc))
}

resolve_log2fc <- function(df, contrast) {
    for (nm in c(paste0("log2FC.imputs.", contrast), paste0("log2FC.", contrast))) {
        if (nm %in% names(df)) return(as.numeric(df[[nm]]))
    }
    for (nm in c(paste0("linearFC.imputs.", contrast), paste0("linearFC.", contrast))) {
        if (nm %in% names(df)) return(signed_fc_to_log2(as.numeric(df[[nm]])))
    }
    rep(NA_real_, nrow(df))
}

# The five real features, with their true log2FC and their stored linearFC.
borderline <- data.frame(
    Gene                      = c("PLS3", "SERINC3", "BLTP3A", "CKAP2", "ISYNA1"),
    `log2FC.imputs.C_vs_V`    = c(-0.58329, -0.58461, 0.58085, 0.58300, 0.58202),
    `linearFC.imputs.C_vs_V`  = c(-1.5, -1.5, 1.5, 1.5, 1.5),
    check.names = FALSE, stringsAsFactors = FALSE
)

test_that("the stored log2FC is preferred over the rounded linearFC", {
    got <- resolve_log2fc(borderline, "C_vs_V")
    expect_equal(got, borderline$`log2FC.imputs.C_vs_V`)
})

test_that("all five borderline features fall BELOW the 1.5-fold threshold", {
    cut <- log2(1.5)
    expect_true(all(abs(resolve_log2fc(borderline, "C_vs_V")) < cut))
    # and their true ratios really are under 1.5
    expect_true(all(2^abs(borderline$`log2FC.imputs.C_vs_V`) < 1.5))
})

test_that("reconstructing from linearFC would wrongly pass all five", {
    # This is the bug, pinned so it cannot come back.
    cut <- log2(1.5)
    from_linear <- signed_fc_to_log2(borderline$`linearFC.imputs.C_vs_V`)
    expect_true(all(abs(from_linear) >= cut))
    expect_equal(sum(abs(from_linear) >= cut), 5L)
})

test_that("the fix changes the count by exactly the borderline features", {
    cut <- log2(1.5)
    correct <- sum(abs(resolve_log2fc(borderline, "C_vs_V")) >= cut)
    buggy <- sum(abs(signed_fc_to_log2(borderline$`linearFC.imputs.C_vs_V`)) >= cut)
    expect_equal(correct, 0L)
    expect_equal(buggy, 5L)
})

test_that("the linearFC fallback still works when log2FC is absent", {
    df <- borderline
    df$`log2FC.imputs.C_vs_V` <- NULL
    got <- resolve_log2fc(df, "C_vs_V")
    expect_equal(got, signed_fc_to_log2(df$`linearFC.imputs.C_vs_V`))
})

test_that("un-suffixed column names are also resolved", {
    df <- data.frame(`log2FC.C_vs_V` = c(1, -1), check.names = FALSE)
    expect_equal(resolve_log2fc(df, "C_vs_V"), c(1, -1))
})

test_that("a contrast with no fold-change column yields NA, not an error", {
    df <- data.frame(Gene = c("A", "B"), stringsAsFactors = FALSE)
    got <- resolve_log2fc(df, "C_vs_V")
    expect_length(got, 2)
    expect_true(all(is.na(got)))
})

test_that("features comfortably past the threshold are unaffected", {
    df <- data.frame(
        `log2FC.imputs.C_vs_V`   = c(2.0, -2.0),
        `linearFC.imputs.C_vs_V` = c(4.0, -4.0),
        check.names = FALSE
    )
    expect_true(all(abs(resolve_log2fc(df, "C_vs_V")) >= log2(1.5)))
})

test_that("the template defines the helper and no longer converts linearFC inline", {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "proteomics",
                            "report_template_proteomics.Rmd"),
        "R/domain/proteomics/report_template_proteomics.Rmd"
    )
    f <- candidates[file.exists(candidates)][1]
    skip_if(is.na(f), "report template not found from the test working directory")

    src <- readLines(f, warn = FALSE)
    expect_gte(length(grep("resolve_log2fc <- function", src, fixed = TRUE)), 1)
    expect_length(grep("signed_fc_to_log2(lfc_vals)", src, fixed = TRUE), 0)
    expect_length(grep("signed_fc_to_log2(as.numeric(df[[lfc_col]]))", src, fixed = TRUE), 0)
})
