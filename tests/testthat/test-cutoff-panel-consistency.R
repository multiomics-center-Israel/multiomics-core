# The interactive cutoff panel overwrites the up/down/total counts printed above
# the volcano tabs, so its client-side rule has to agree with the pipeline's own
# significance call. Two defects made it disagree by one feature:
#
#   1. p-values were serialised at jsonlite's default 4 decimals, so an adjusted
#      p of 0.05004524 reached the browser as 0.0500 and satisfied "<= 0.05".
#   2. classify() skipped its pass-flag gate when the flag was null. A null
#      element is a real NA from the pipeline and means "did not pass", not
#      "no flag supplied".
#
# Reproduced on a real run where the static table said 222 and the panel said 223.

test_that("p-values are serialised at full precision, not rounded to 4 decimals", {
    df <- data.frame(
        name     = c("Lrp1", "Ipo7", "Pfkm"),
        logFC    = c(-1.5705, -1.5410, 0.9900),
        pval     = c(0.0027485658957283657, 0.0104, 0.0011),
        adjPval  = c(0.05004524, 0.05029277, 0.04957949),
        passFlag = c(NA, NA, 1),
        stringsAsFactors = FALSE
    )
    js <- cutoff_register_plot("volcano_test", df, plot_type = "volcano",
                               entity_label = "Protein", contrast = "A_vs_B")

    # The value that caused the off-by-one must survive verbatim.
    expect_true(grepl("0.05004524", js, fixed = TRUE))
    expect_true(grepl("0.05029277", js, fixed = TRUE))

    # And must not appear pre-rounded onto the cutoff.
    expect_false(grepl("[0.05,", js, fixed = TRUE))
    expect_false(grepl(",0.05,", js, fixed = TRUE))
})

test_that("logFC is serialised at full precision too", {
    df <- data.frame(
        name    = "F1",
        logFC   = 0.5849625007211562,  # exactly log2(1.5), the default cutoff
        pval    = 0.01,
        adjPval = 0.02,
        stringsAsFactors = FALSE
    )
    js <- cutoff_register_plot("volcano_lfc", df)
    expect_true(grepl("0.58496250072", js, fixed = TRUE))
})

test_that("the emitted classify() treats a null pass flag as not significant", {
    panel <- cutoff_panel_html(default_lfc = log2(1.5), default_p = 0.05)

    # The old guard let null through; the fixed guard must not mention null.
    expect_false(grepl("passFlag !== null", panel, fixed = TRUE))
    expect_true(grepl("passFlag !== undefined && passFlag !== 1", panel, fixed = TRUE))
})

test_that("a missing pass-flag column still leaves the gate open", {
    # Pipelines without an imputation-consensus vote supply no flag at all.
    # Those features must stay classifiable on p-value and fold change alone.
    df <- data.frame(
        name    = "F1",
        logFC   = 2,
        pval    = 0.001,
        adjPval = 0.002,
        stringsAsFactors = FALSE
    )
    js <- cutoff_register_plot("volcano_noflag", df)
    expect_true(grepl("passFlag: null", js, fixed = TRUE))

    panel <- cutoff_panel_html()
    # The call site must pass undefined only when the whole array is absent.
    expect_true(grepl("data.passFlag ? data.passFlag[i] : undefined", panel, fixed = TRUE))
})

test_that("NA pass flags are serialised as null, not dropped or coerced to 0", {
    df <- data.frame(
        name     = c("F1", "F2", "F3"),
        logFC    = c(1, -1, 2),
        pval     = c(0.01, 0.02, 0.03),
        adjPval  = c(0.01, 0.02, 0.03),
        passFlag = c(1, NA, 0),
        stringsAsFactors = FALSE
    )
    js <- cutoff_register_plot("volcano_flags", df)
    expect_true(grepl("passFlag: [1,null,0]", js, fixed = TRUE))
})
