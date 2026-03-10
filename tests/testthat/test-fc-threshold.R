# Tests for FC threshold counting logic in pipeline summary
# Regression tests for the log2FC threshold bug fixed in PR #48

test_that("log2FC up/down counts use actual threshold, not zero", {
    # Simulate a DE table with log2FoldChange values around the threshold
    fc_cutoff <- 1.5
    log2_fc <- log2(fc_cutoff)  # ~0.585

    tbl <- data.frame(
        padj = c(0.01, 0.01, 0.01, 0.01, 0.5),
        log2FoldChange = c(0.6, -0.7, 0.3, -0.2, 1.0)
    )

    p_cut <- 0.05
    sig <- !is.na(tbl$padj) & tbl$padj <= p_cut
    if (log2_fc > 0) sig <- sig & (abs(tbl$log2FoldChange) >= log2_fc)

    up <- sum(sig & tbl$log2FoldChange >= log2_fc, na.rm = TRUE)
    dn <- sum(sig & tbl$log2FoldChange <= -log2_fc, na.rm = TRUE)

    # 0.6 >= 0.585 -> up; -0.7 <= -0.585 -> down
    # 0.3 < 0.585 -> filtered out; -0.2 > -0.585 -> filtered out
    expect_equal(up, 1)
    expect_equal(dn, 1)
})

test_that("log2FC = 0.3 is NOT counted when threshold is log2(1.5)", {
    tbl <- data.frame(
        padj = c(0.01),
        log2FoldChange = c(0.3)
    )

    fc_cutoff <- 1.5
    log2_fc <- log2(fc_cutoff)

    sig <- !is.na(tbl$padj) & tbl$padj <= 0.05
    if (log2_fc > 0) sig <- sig & (abs(tbl$log2FoldChange) >= log2_fc)

    up <- sum(sig & tbl$log2FoldChange >= log2_fc, na.rm = TRUE)
    expect_equal(up, 0)
})

test_that("linearFC threshold filters correctly", {
    # linearFC: >1 means up, <1 means down.
    # Current filter: abs(fc) >= 1.5 OR 1/abs(fc) <= 1/1.5
    # FC=2.0 passes (2.0 >= 1.5), counted as up
    # FC=1.1 fails (1.1 < 1.5 and 1/1.1=0.91 > 0.667)
    tbl <- data.frame(
        padj = c(0.01, 0.01),
        linearFC = c(2.0, 1.1)
    )

    fc_lin <- 1.5
    sig <- !is.na(tbl$padj) & tbl$padj <= 0.05
    sig <- sig & (abs(tbl$linearFC) >= fc_lin | (1 / abs(tbl$linearFC)) <= (1 / fc_lin))

    up <- sum(sig & tbl$linearFC > 1, na.rm = TRUE)
    dn <- sum(sig & tbl$linearFC < 1, na.rm = TRUE)

    expect_equal(up, 1)
    expect_equal(dn, 0)  # 1.1 filtered out
})
