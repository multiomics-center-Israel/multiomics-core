# tests/testthat/test-proteomics-qc-boxplot-title.R
#
# The post-imputation QC boxplot used the generic "Normalized expression
# boxplots" title, while the report tab called it "after imputation". The title
# now names each step the plotted matrix went through, and only those that ran.
#
# proteomics_qc_boxplot_title() lives in R/modules/proteomics/01_mod_qc_pre.R
# and is loaded by helper.R along with the rest of R/.

test_that("pre-normalized input is not called normalized", {
    expect_identical(
        proteomics_qc_boxplot_title(list(), list(normalization = list(method = "none"))),
        "Expression boxplots (after imputation)"
    )
    # A missing normalization block resolves to "none", as preprocess_proteomics() does.
    expect_identical(
        proteomics_qc_boxplot_title(list(), list()),
        "Expression boxplots (after imputation)"
    )
})

test_that("configured normalization is named", {
    expect_identical(
        proteomics_qc_boxplot_title(list(), list(normalization = list(method = "median"))),
        "Expression boxplots (after median normalization + imputation)"
    )
})

test_that("batch correction is named only when it ran", {
    cfg <- list(normalization = list(method = "median"),
                batch_correction = list(method = "combat"))

    # Configured but not applied: the runtime flag decides, not the config.
    expect_identical(
        proteomics_qc_boxplot_title(list(batch_corrected = NULL), cfg),
        "Expression boxplots (after median normalization + imputation)"
    )
    expect_identical(
        proteomics_qc_boxplot_title(list(batch_corrected = TRUE, batch_method = "combat"), cfg),
        "Expression boxplots (after median normalization + imputation + combat batch correction)"
    )
    expect_identical(
        proteomics_qc_boxplot_title(list(batch_corrected = TRUE), cfg),
        "Expression boxplots (after median normalization + imputation + batch correction)"
    )
})
