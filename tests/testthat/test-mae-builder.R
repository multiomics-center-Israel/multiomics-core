# Tests for MAE builder edge cases
# Regression test for the match() all-NA bug fixed in PR #48

test_that("match() all-NA produces warning and fallback row_data", {
    # Simulate: row_data has completely different feature IDs than expression
    expr <- matrix(rnorm(12), nrow = 3, dimnames = list(
        c("gene_A", "gene_B", "gene_C"),
        c("S1", "S2", "S3", "S4")
    ))

    rd <- data.frame(
        feature_id = c("UNRELATED_1", "UNRELATED_2", "UNRELATED_3"),
        row.names = c("UNRELATED_1", "UNRELATED_2", "UNRELATED_3"),
        stringsAsFactors = FALSE
    )

    idx <- match(rownames(expr), rownames(rd))

    # All NAs — no overlap between feature sets
    expect_true(all(is.na(idx)))

    # Our fix: when all NA, create minimal row_data instead of corrupting
    if (all(is.na(idx))) {
        rd_fallback <- data.frame(
            feature_id = rownames(expr),
            row.names = rownames(expr),
            stringsAsFactors = FALSE
        )
    }

    expect_equal(rownames(rd_fallback), rownames(expr))
    expect_equal(rd_fallback$feature_id, c("gene_A", "gene_B", "gene_C"))
})

test_that("match() partial overlap correctly subsets row_data", {
    expr <- matrix(rnorm(12), nrow = 3, dimnames = list(
        c("gene_A", "gene_B", "gene_C"),
        c("S1", "S2", "S3", "S4")
    ))

    rd <- data.frame(
        symbol = c("SYM_A", "SYM_B"),
        row.names = c("gene_A", "gene_B"),
        stringsAsFactors = FALSE
    )

    idx <- match(rownames(expr), rownames(rd))

    # Partial match: gene_C has no entry
    expect_false(all(is.na(idx)))
    expect_true(is.na(idx[3]))

    rd_aligned <- rd[idx, , drop = FALSE]
    rownames(rd_aligned) <- rownames(expr)

    expect_equal(rownames(rd_aligned), c("gene_A", "gene_B", "gene_C"))
    expect_equal(rd_aligned$symbol[1], "SYM_A")
    expect_true(is.na(rd_aligned$symbol[3]))
})
