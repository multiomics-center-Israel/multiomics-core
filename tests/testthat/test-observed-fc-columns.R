# Tests for compute_observed_fc_columns(): the measured-only fold change that
# sits beside the pipeline's imputation-based linearFC in the final table.

make_meta <- function() {
    data.frame(
        SampleID = paste0("S", 1:6),
        Group    = rep(c("A", "B"), each = 3),
        stringsAsFactors = FALSE
    )
}

make_contrasts <- function(name = "A_vs_B") {
    data.frame(
        Contrast_name = name,
        Factor        = "Group",
        Numerator     = "A",
        Denominator   = "B",
        stringsAsFactors = FALSE
    )
}

make_expr <- function(values, features = "F1") {
    matrix(values,
           nrow = length(features), byrow = TRUE,
           dimnames = list(features, paste0("S", 1:6)))
}

test_that("observed log2FC is the difference of complete group means", {
    expr <- make_expr(c(10, 12, 14, 4, 6, 8))  # A mean 12, B mean 6
    res <- compute_observed_fc_columns(expr, make_meta(), "SampleID", make_contrasts())

    expect_equal(res[["obs.log2FC.A_vs_B"]], 6)
    expect_equal(res[["obs.linearFC.A_vs_B"]], 64)
    expect_equal(res[["n_obs.A"]], 3L)
    expect_equal(res[["n_obs.B"]], 3L)
})

test_that("missing values are dropped from the group mean, not treated as zero", {
    # A has one missing value; the mean must use the two observed values only.
    expr <- make_expr(c(10, 12, NA, 4, 6, 8))
    res <- compute_observed_fc_columns(expr, make_meta(), "SampleID", make_contrasts())

    expect_equal(res[["obs.log2FC.A_vs_B"]], 11 - 6)
    expect_equal(res[["n_obs.A"]], 2L)
    expect_equal(res[["n_obs.B"]], 3L)
})

test_that("a group with nothing observed yields NA, not a one-armed difference", {
    expr <- make_expr(c(NA, NA, NA, 4, 6, 8))
    res <- compute_observed_fc_columns(expr, make_meta(), "SampleID", make_contrasts())

    expect_true(is.na(res[["obs.log2FC.A_vs_B"]]))
    expect_true(is.na(res[["obs.linearFC.A_vs_B"]]))
    expect_equal(res[["n_obs.A"]], 0L)
})

test_that("linearFC is signed: ratios below 1 are written as -1/ratio", {
    expr <- make_expr(c(4, 4, 4, 6, 6, 6))  # log2FC = -2, ratio = 0.25
    res <- compute_observed_fc_columns(expr, make_meta(), "SampleID", make_contrasts())

    expect_equal(res[["obs.log2FC.A_vs_B"]], -2)
    expect_equal(res[["obs.linearFC.A_vs_B"]], -4)
})

test_that("contrast names are normalized for proteomics but left alone elsewhere", {
    expr <- make_expr(c(10, 12, 14, 4, 6, 8))
    ctr  <- make_contrasts("A vs B")

    prot <- compute_observed_fc_columns(expr, make_meta(), "SampleID", ctr,
                                        mode = "proteomics")
    metab <- compute_observed_fc_columns(expr, make_meta(), "SampleID", ctr,
                                         mode = "metabolomics")

    expect_true(paste0("obs.log2FC.", normalize_contrast_name("A vs B")) %in% names(prot))
    expect_true("obs.log2FC.A vs B" %in% names(metab))
})

test_that("multiple features keep their identity and row order", {
    expr <- matrix(
        c(10, 12, 14, 4, 6, 8,
          1, 1, 1, 1, 1, 1),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("F1", "F2"), paste0("S", 1:6))
    )
    res <- compute_observed_fc_columns(expr, make_meta(), "SampleID", make_contrasts())

    expect_equal(rownames(res), c("F1", "F2"))
    expect_equal(res[["obs.log2FC.A_vs_B"]], c(6, 0))
})

test_that("insufficient inputs return NULL rather than a malformed table", {
    expr <- make_expr(c(10, 12, 14, 4, 6, 8))

    expect_null(compute_observed_fc_columns(NULL, make_meta(), "SampleID", make_contrasts()))
    expect_null(compute_observed_fc_columns(expr, make_meta(), "NotAColumn", make_contrasts()))
    expect_null(compute_observed_fc_columns(expr, make_meta(), "SampleID",
                                            data.frame(x = 1)))
})

test_that("an unresolvable grouping column warns and skips instead of guessing", {
    expr <- make_expr(c(10, 12, 14, 4, 6, 8))
    ctr  <- make_contrasts()
    ctr$Factor <- "NoSuchColumn"

    expect_warning(
        res <- compute_observed_fc_columns(expr, make_meta(), "SampleID", ctr),
        "could not resolve a single grouping column"
    )
    expect_null(res)
})

# --- Column_guide note must describe the workbook it is written into ---------
# Regression tests for review feedback on PR #154: the note claimed blanks were
# filled in an imp.<sample> block even when that block was not exported.

make_guide_df <- function(with_imp) {
    df <- data.frame(
        FeatureID = c("F1", "F2"),
        S1 = c(10, 11), S2 = c(12, 13),
        `CV.A` = c(3, 4), `n_obs.A` = c(2L, 2L),
        `obs.log2FC.A_vs_B` = c(0.5, -0.5),
        `linearFC.A_vs_B` = c(1.4, -1.4),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    if (with_imp) {
        df[["imp.S1"]] <- c(10, 11)
        df[["imp.S2"]] <- c(12, 13)
    }
    df
}

make_cfg <- function(imp_method, reps = 1) {
    list(modes = list(proteomics = list(
        imputation = list(method = imp_method, no_repetitions = reps)
    )))
}

guide_note <- function(df, cfg) {
    wb <- openxlsx::createWorkbook()
    add_column_guide_sheet(wb, df, cfg, mode = "proteomics")
    openxlsx::readWorkbook(wb, sheet = "Column_guide", colNames = FALSE,
                           rows = 1, cols = 1)[1, 1]
}

test_that("the guide note mentions the imp block only when that block is exported", {
    note <- guide_note(make_guide_df(TRUE), make_cfg("perseus"))
    expect_true(grepl("imp.<sample> block", note, fixed = TRUE))
})

test_that("with imputation on but the block withheld, the note says so instead", {
    note <- guide_note(make_guide_df(FALSE), make_cfg("perseus"))
    expect_false(grepl("imp.<sample> block", note, fixed = TRUE))
    expect_true(grepl("not exported in this workbook", note, fixed = TRUE))
    expect_true(grepl("perseus", note, fixed = TRUE))
})

test_that("with no imputation the note does not claim anything was filled in", {
    note <- guide_note(make_guide_df(FALSE), make_cfg("none"))
    expect_true(grepl("No imputation was performed", note, fixed = TRUE))
    expect_true(grepl("nothing was filled in", note, fixed = TRUE))
    expect_false(grepl("imp.<sample> block", note, fixed = TRUE))
    expect_false(grepl("were filled", note, fixed = TRUE))
})

test_that("build_final_results_proteomics tolerates a NULL config", {
    # config has always been optional. The imputed-block flag lookup must not
    # depend on NULL$a$b semantics; without a config the block is exported.
    expr <- matrix(c(10, 12, 14, 4, 6, 8), nrow = 1,
                   dimnames = list("F1", paste0("S", 1:6)))
    pre <- list(
        expr_filt = expr, expr_imp_single = expr,
        meta = data.frame(SampleID = paste0("S", 1:6),
                          Group = rep(c("A", "B"), each = 3),
                          stringsAsFactors = FALSE),
        row_data = NULL
    )
    summary_df <- data.frame(
        FeatureID = "F1",
        `linearFC.imputs.A_vs_B` = 1.5,
        `pvalue.imputs.A_vs_B` = 0.01,
        `padj.imputs.A_vs_B` = 0.02,
        `pass.imputs.A_vs_B` = 1,
        pass_any_contrast = 1,
        check.names = FALSE, stringsAsFactors = FALSE
    )
    ctr <- data.frame(Contrast_name = "A_vs_B", Factor = "Group",
                      Numerator = "A", Denominator = "B",
                      stringsAsFactors = FALSE)

    expect_no_error(
        res <- build_final_results_proteomics(pre, summary_df, ctr, config = NULL)
    )
    expect_true(any(grepl("^imp\\.", names(res))))
})
