# tests/testthat/test-proteomics-qc-boxplot-title.R
#
# The post-imputation QC boxplot used the generic "Normalized expression
# boxplots" title, while the report tab called it "after imputation". The title
# now names each step the plotted matrix went through, and only those that ran.
#
# proteomics_qc_boxplot_title() lives in R/modules/proteomics/01_mod_qc_pre.R
# and is loaded by helper.R along with the rest of R/.

cfg_with <- function(norm = "none", imp = "none", ...) {
    c(list(normalization = list(method = norm),
           imputation    = list(method = imp)), list(...))
}

test_that("pre-normalized input is not called normalized", {
    expect_identical(
        proteomics_qc_boxplot_title(list(), cfg_with(imp = "perseus_like")),
        "Expression boxplots (after perseus_like imputation)"
    )
})

test_that("imputation method none is not called imputed", {
    # The boxplot is still drawn under method none (the imputer returns
    # missingness flags), so its existence says nothing about imputation.
    expect_identical(
        proteomics_qc_boxplot_title(list(), cfg_with(norm = "median", imp = "none")),
        "Expression boxplots (after median normalization)"
    )
    expect_identical(
        proteomics_qc_boxplot_title(list(), cfg_with()),
        "Expression boxplots (no normalization, imputation or batch correction)"
    )
})

test_that("imputation is resolved the way impute_proteomics() resolves it", {
    # Missing keys: normalization defaults to none, imputation to perseus_like.
    expect_identical(
        proteomics_qc_boxplot_title(list(), list()),
        "Expression boxplots (after perseus_like imputation)"
    )
    # "perseus" is an alias of perseus_like.
    expect_identical(
        proteomics_qc_boxplot_title(list(), cfg_with(imp = "perseus")),
        "Expression boxplots (after perseus_like imputation)"
    )
})

test_that("batch correction is named only when it ran", {
    cfg <- cfg_with(norm = "median", imp = "minval",
                    batch_correction = list(method = "combat"))

    # Configured but not applied: the runtime flag decides, not the config.
    expect_identical(
        proteomics_qc_boxplot_title(list(batch_corrected = NULL), cfg),
        "Expression boxplots (after median normalization + minval imputation)"
    )
    expect_identical(
        proteomics_qc_boxplot_title(list(batch_corrected = TRUE, batch_method = "combat"), cfg),
        "Expression boxplots (after median normalization + minval imputation + combat batch correction)"
    )
    expect_identical(
        proteomics_qc_boxplot_title(list(batch_corrected = TRUE), cfg),
        "Expression boxplots (after median normalization + minval imputation + batch correction)"
    )
})
