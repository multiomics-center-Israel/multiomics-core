# tests/testthat/test-proteomics-gsea-rank-direction.R
#
# The default GSEA ranking in extract_de_table_for_pathway() is
# sign(log2FC) * -log10(pvalue). Two things have to hold at once:
#
#   - the sign must survive linearFC's signif(x, 3) rounding, which collapses
#     any ratio in [0.995, 1.005) to exactly 1 and would otherwise zero the rank
#     whatever the p-value;
#   - the sign must come from a SIGNED column. linearRatio.imputs is not one:
#     the precomputed-input path writes it as 2^abs(logFC)
#     (R/domain/proteomics/05_de_summary.R), so every feature looks upregulated.

# Mirrors the precomputed-input shape: an unsigned linearRatio beside a signed
# log2FC, which is what makes the column choice observable.
rank_summary_df <- function(lfc, pvals, contrast = "B_vs_A", with_log2fc = TRUE) {
    cn <- normalize_contrast_name(contrast)
    linear_ratio <- 2^abs(lfc)                                   # always >= 1
    linear_fc <- ifelse(lfc >= 0, linear_ratio, -linear_ratio)

    df <- data.frame(
        FeatureID = paste0("PROT", seq_along(lfc)),
        stringsAsFactors = FALSE
    )
    df[[paste0("linearFC.imputs.", cn)]]    <- signif(linear_fc, 3)
    df[[paste0("linearRatio.imputs.", cn)]] <- linear_ratio
    df[[paste0("pvalue.imputs.", cn)]]      <- pvals
    df[[paste0("padj.imputs.", cn)]]        <- pvals
    if (with_log2fc) df[[paste0("log2FC.imputs.", cn)]] <- lfc
    df
}

rank_config <- function() {
    list(modes = list(proteomics = list(de_table = list(id_col = "FeatureID"))))
}

test_that("downregulated proteins keep a negative GSEA rank", {
    # PROT2 is down two-fold. linearRatio.imputs is 4 for it, exactly as for the
    # up-regulated PROT1, so taking direction from that column would rank both
    # as up.
    lfc   <- c(2, -2)
    pvals <- c(1e-4, 1e-6)

    res <- extract_de_table_for_pathway(
        rank_summary_df(lfc, pvals), "B_vs_A", rank_config()
    )

    expect_gt(res$stat[1], 0)
    expect_lt(res$stat[2], 0)
    expect_equal(res$stat[2], -(-log10(pvals[2] + 1e-300)), tolerance = 1e-8)
})

test_that("a fold change that rounds to 1.00 still ranks by its p-value", {
    # 2^0.001 = 1.00069, which signif(x, 3) stores as 1.00, so sign(linearFC)
    # is 0 here. The rank must come from the p-value, not collapse to zero.
    res <- extract_de_table_for_pathway(
        rank_summary_df(c(0.001), c(1e-5)), "B_vs_A", rank_config()
    )

    expect_gt(abs(res$stat[1]), 0)
    expect_equal(res$stat[1], -log10(1e-5 + 1e-300), tolerance = 1e-8)
})

test_that("a protein with no fold-change estimate gets no rank", {
    # A precomputed DE row can carry a p-value and no fold-change estimate.
    # Inventing a direction for it would seat it in the upregulated tail on the
    # strength of that p-value; run_pathway_analysis() drops an NA rank instead.
    res <- extract_de_table_for_pathway(
        rank_summary_df(c(2, NA), c(1e-4, 1e-6)), "B_vs_A", rank_config()
    )

    expect_gt(res$stat[1], 0)
    expect_true(is.na(res$stat[2]))
})

test_that("a genuine zero fold change ranks neutrally", {
    # Not the [0.995, 1.005) rounding case: log2FC.imputs is exactly zero here,
    # so there is a direction to read and it is neither up nor down.
    res <- extract_de_table_for_pathway(
        rank_summary_df(c(0), c(1e-8)), "B_vs_A", rank_config()
    )

    expect_equal(res$stat[1], 0)
})

test_that("direction falls back to the rounded linearFC when log2FC is absent", {
    lfc   <- c(2, -2)
    pvals <- c(1e-4, 1e-6)

    res <- extract_de_table_for_pathway(
        rank_summary_df(lfc, pvals, with_log2fc = FALSE), "B_vs_A", rank_config()
    )

    # Rounding does not reach these two, so the fallback still gets the signs
    # right; only the [0.995, 1.005) window is lost without log2FC.imputs.
    expect_gt(res$stat[1], 0)
    expect_lt(res$stat[2], 0)
})
