# tests/testthat/test-proteomics-clustering-de-cutoffs.R
#
# Rows of the proteomics clustering heatmap:
#   - get_clustering_de_cutoffs(), apply_clustering_de_cutoffs() (05_de_summary.R)
#
# Synthetic data only.

test_that("clustering cut-offs are off by default and fall back to the de block", {
    expect_null(get_clustering_de_cutoffs(list(clustering = list(enabled = TRUE))))
    co <- get_clustering_de_cutoffs(list(
        clustering = list(de_log2fc_cutoff = 0.6),
        de = list(p_cutoff = 0.01, linear_fc_cutoff = 1.5)))
    expect_equal(co$padj_cutoff, 0.01)
    expect_equal(co$log2fc_cutoff, 0.6)
    co2 <- get_clustering_de_cutoffs(list(clustering = list(de_padj_cutoff = 0.05),
                                          de = list(linear_fc_cutoff = 1.5)))
    expect_equal(co2$log2fc_cutoff, log2(1.5))
})

test_that("clustering cut-offs re-mark the pass columns with padj < and |log2FC| >=", {
    df <- data.frame(
        FeatureID = paste0("P", 1:5),
        padj.imputs.a = c(0.01, 0.01, 0.05, 0.01, NA),
        log2FC.imputs.a = c(0.59, -0.6, 1.0, 0.2, 2),
        pass.imputs.a = c(1, 1, 1, NA, NA),
        padj.imputs.b = c(0.2, 0.2, 0.2, 0.001, 0.2),
        log2FC.imputs.b = c(0, 0, 0, -0.7, 0),
        pass.imputs.b = c(NA, NA, NA, 1, NA),
        stringsAsFactors = FALSE)
    expect_identical(apply_clustering_de_cutoffs(df, NULL), df)

    out <- apply_clustering_de_cutoffs(df, list(padj_cutoff = 0.05, log2fc_cutoff = 0.6))
    # P1: |0.59| < 0.6; P2: exactly 0.6 passes; P3: padj 0.05 is not < 0.05.
    expect_equal(out$pass.imputs.a, c(NA, 1, NA, NA, NA))
    expect_equal(out$pass.imputs.b, c(NA, NA, NA, 1, NA))
    expect_equal(out$pass_any_contrast, c(NA, 1, NA, 1, NA))
    expect_equal(get_de_features(list(summary_df = out), list()), c("P2", "P4"))
})
