# tests/testthat/test-metab-de-summary-pval-gate.R
#
# build_de_summary() gates significance on either the BH-adjusted p-value or the
# raw p-value, controlled by `use_adjusted_pval` (default TRUE = adjusted, the
# pipeline default). The metabolomics pipeline sets it from de$use_adjusted_pval.
# These tests pin both branches using a feature (F2) that is significant under
# the raw p-value but not under the adjusted one, with |log2FC| large enough that
# the fold-change gate never binds.

make_de_tables <- function() {
    tbl <- data.frame(
        feature_id = c("F1", "F2", "F3"),
        logFC      = c(1, 1, 1),              # |log2FC| >= log2(1.5); FC gate never binds
        AveExpr    = c(10, 10, 10),
        P.Value    = c(0.001, 0.010, 0.500),
        adj.P.Val  = c(0.001, 0.200, 0.600),  # F2: raw-significant, adjusted not
        stringsAsFactors = FALSE
    )
    list(A_vs_B = tbl)
}

test_that("use_adjusted_pval = TRUE gates on the BH-adjusted p-value", {
    out  <- build_de_summary(make_de_tables(), padj_cutoff = 0.05,
                             log2fc_cut = log2(1.5), use_adjusted_pval = TRUE)
    pass <- setNames(out$pass.A_vs_B, out$feature_id)
    expect_equal(pass[["F1"]], 1L)
    expect_equal(pass[["F2"]], 0L)   # adj p = 0.200 >= 0.05 -> dropped
    expect_equal(pass[["F3"]], 0L)
})

test_that("use_adjusted_pval = FALSE gates on the raw p-value", {
    out  <- build_de_summary(make_de_tables(), padj_cutoff = 0.05,
                             log2fc_cut = log2(1.5), use_adjusted_pval = FALSE)
    pass <- setNames(out$pass.A_vs_B, out$feature_id)
    expect_equal(pass[["F1"]], 1L)
    expect_equal(pass[["F2"]], 1L)   # raw p = 0.010 < 0.05
    expect_equal(pass[["F3"]], 0L)
})

test_that("the two gates disagree only on the raw-only feature, and the default is adjusted", {
    adj <- build_de_summary(make_de_tables(), 0.05, log2(1.5), use_adjusted_pval = TRUE)
    raw <- build_de_summary(make_de_tables(), 0.05, log2(1.5), use_adjusted_pval = FALSE)

    expect_equal(sum(adj$pass_any_contrast), 1L)  # F1 only
    expect_equal(sum(raw$pass_any_contrast), 2L)  # F1 + F2

    # The function default (no arg) must match the adjusted branch, so callers
    # that omit the argument (e.g. lipidomics) keep main's adjusted-p behavior.
    dflt <- build_de_summary(make_de_tables(), 0.05, log2(1.5))
    expect_equal(dflt$pass.A_vs_B, adj$pass.A_vs_B)
})
