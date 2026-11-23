# tests/testthat/test-metab-qc-subsets.R
#
# build_qc_subsets() adds a "_noQC" subset when QC samples exist. Metadata often
# marks pools only through an is_QC flag or a "Pool" group, with no treatment
# column. Given a condition column, those samples must still be split off, by
# the same filter_to_biological() rule that DE and clustering use. Without one,
# the older treatment-only rule applies unchanged.

qc_fixture <- function(meta) {
    mat <- matrix(seq_len(3 * nrow(meta)), nrow = 3,
                  dimnames = list(paste0("F", 1:3), meta$sample_id))
    list(mat = mat, meta = meta)
}

test_that("an is_QC flag alone gives a noQC subset when condition_col is set", {
    fx <- qc_fixture(data.frame(
        sample_id = c("S1", "S2", "S3", "S4", "S5"),
        condition = c("A", "A", "B", "B", "Mix"),
        is_QC     = c(FALSE, FALSE, FALSE, FALSE, TRUE),
        stringsAsFactors = FALSE
    ))
    subsets <- suppressMessages(build_qc_subsets(
        fx$mat, fx$mat, fx$meta, "sample_id", expr_log = fx$mat,
        condition_col = "condition"
    ))

    expect_length(subsets, 2)
    expect_identical(subsets[[2]]$tag, "_noQC")
    expect_identical(colnames(subsets[[2]]$expr_work), c("S1", "S2", "S3", "S4"))
    expect_identical(colnames(subsets[[2]]$expr_filt), c("S1", "S2", "S3", "S4"))
    expect_identical(colnames(subsets[[2]]$expr_log),  c("S1", "S2", "S3", "S4"))
    expect_identical(subsets[[2]]$meta$sample_id,      c("S1", "S2", "S3", "S4"))
    # The first subset always keeps every sample
    expect_identical(colnames(subsets[[1]]$expr_work), fx$meta$sample_id)
})

test_that("a Pool group is split off without is_QC or treatment columns", {
    fx <- qc_fixture(data.frame(
        sample_id = c("S1", "S2", "S3", "S4", "S5"),
        condition = c("A", "A", "B", "B", "Pool"),
        stringsAsFactors = FALSE
    ))
    subsets <- suppressMessages(build_qc_subsets(
        fx$mat, fx$mat, fx$meta, "sample_id", condition_col = "condition"
    ))

    expect_length(subsets, 2)
    expect_identical(subsets[[2]]$meta$sample_id, c("S1", "S2", "S3", "S4"))
    expect_null(subsets[[2]]$expr_log)
})

test_that("without condition_col only the treatment column is used", {
    fx <- qc_fixture(data.frame(
        sample_id = c("S1", "S2", "S3"),
        condition = c("A", "B", "Pool"),
        is_QC     = c(FALSE, FALSE, TRUE),
        stringsAsFactors = FALSE
    ))
    expect_length(build_qc_subsets(fx$mat, fx$mat, fx$meta, "sample_id"), 1)

    fx$meta$treatment <- c("ctrl", "ctrl", "QC")
    subsets <- build_qc_subsets(fx$mat, fx$mat, fx$meta, "sample_id")
    expect_length(subsets, 2)
    expect_identical(subsets[[2]]$meta$sample_id, c("S1", "S2"))
})

test_that("no subset is added when every sample is biological", {
    fx <- qc_fixture(data.frame(
        sample_id = c("S1", "S2", "S3", "S4"),
        condition = c("A", "A", "B", "B"),
        is_QC     = FALSE,
        stringsAsFactors = FALSE
    ))
    subsets <- build_qc_subsets(fx$mat, fx$mat, fx$meta, "sample_id",
                                condition_col = "condition")
    expect_length(subsets, 1)
})
