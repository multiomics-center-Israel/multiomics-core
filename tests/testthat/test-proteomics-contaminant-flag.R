# tests/testthat/test-proteomics-contaminant-flag.R
#
# Tests for filter_contaminants() honouring filtering$remove_contaminants
# (R/domain/proteomics/02_filtering.R). The flag was validated but never read,
# so contaminants were removed unconditionally. These tests pin the three cases
# that matter: flag absent (remove, the historical default), flag TRUE (remove),
# flag FALSE (keep).

make_expr <- function() {
    ids <- c("P1", "cRAP-KRT1", "P2", "cRAP-TRYP", "CXCL8;cRAP-IL8")
    m <- matrix(seq_len(length(ids) * 4), nrow = length(ids),
                dimnames = list(ids, paste0("S", 1:4)))
    storage.mode(m) <- "double"
    m
}

make_row_data <- function(expr) {
    data.frame(FeatureID = rownames(expr), stringsAsFactors = FALSE)
}

cfg_with <- function(...) list(filtering = list(...))

test_that("contaminants are removed when the flag is absent", {
    e <- make_expr()
    res <- filter_contaminants(e, make_row_data(e), cfg_with(contaminant_prefix = "cRAP-"))

    expect_equal(res$n_removed, 2)
    expect_equal(nrow(res$expr_mat), 3)
    expect_false(any(startsWith(rownames(res$expr_mat), "cRAP-")))
})

test_that("contaminants are removed when the flag is TRUE", {
    e <- make_expr()
    res <- filter_contaminants(
        e, make_row_data(e),
        cfg_with(remove_contaminants = TRUE, contaminant_prefix = "cRAP-"))

    expect_equal(res$n_removed, 2)
    expect_equal(nrow(res$expr_mat), 3)
})

test_that("contaminants are kept when the flag is FALSE", {
    e <- make_expr()
    res <- filter_contaminants(
        e, make_row_data(e),
        cfg_with(remove_contaminants = FALSE, contaminant_prefix = "cRAP-"))

    expect_equal(res$n_removed, 0L)
    expect_equal(nrow(res$expr_mat), nrow(e))
    expect_identical(rownames(res$expr_mat), rownames(e))
})

test_that("disabling reports how many features were kept", {
    e <- make_expr()
    expect_message(
        filter_contaminants(e, make_row_data(e),
                            cfg_with(remove_contaminants = FALSE,
                                     contaminant_prefix = "cRAP-")),
        "DISABLED.*keeping 2 features"
    )
})

test_that("row_data stays aligned to expr_mat in both modes", {
    e <- make_expr()
    rd <- make_row_data(e)

    kept <- filter_contaminants(e, rd, cfg_with(remove_contaminants = FALSE))
    expect_identical(kept$row_data$FeatureID, rownames(kept$expr_mat))

    dropped <- filter_contaminants(e, rd, cfg_with(remove_contaminants = TRUE))
    expect_identical(dropped$row_data$FeatureID, rownames(dropped$expr_mat))
})

test_that("the prefix match stays anchored in both modes", {
    e <- make_expr()
    # "CXCL8;cRAP-IL8" contains the prefix but does not start with it.
    dropped <- filter_contaminants(e, make_row_data(e),
                                   cfg_with(remove_contaminants = TRUE))
    expect_true("CXCL8;cRAP-IL8" %in% rownames(dropped$expr_mat))
})

test_that("a custom prefix is honoured", {
    e <- make_expr()
    res <- filter_contaminants(e, make_row_data(e),
                               cfg_with(contaminant_prefix = "P"))
    expect_equal(res$n_removed, 2)
    expect_false(any(startsWith(rownames(res$expr_mat), "P")))
})

# The prefix is config-supplied and validated only as a non-empty string, so it
# must be matched literally rather than compiled as a regular expression.
make_pipe_expr <- function() {
    ids <- c("sp|P12345", "sp|Q67890", "TRYP_PIG", "KRT1_HUMAN")
    m <- matrix(seq_len(length(ids) * 4), nrow = length(ids),
                dimnames = list(ids, paste0("S", 1:4)))
    storage.mode(m) <- "double"
    m
}

test_that("a prefix containing '|' removes only the features that start with it", {
    e <- make_pipe_expr()
    res <- suppressMessages(
        filter_contaminants(e, make_row_data(e),
                            cfg_with(contaminant_prefix = "sp|"))
    )

    # As a pattern, "^sp|" is an alternation with the empty string and matches
    # every id, silently emptying the matrix.
    expect_equal(res$n_removed, 2)
    expect_identical(rownames(res$expr_mat), c("TRYP_PIG", "KRT1_HUMAN"))
})

test_that("a prefix that is not a valid regex does not abort either mode", {
    e <- make_pipe_expr()

    # A bare "[" fails to compile as a pattern, so reaching grepl() at all would
    # abort preprocessing — including on the disabled path, which promises to
    # skip filtering entirely.
    kept <- suppressMessages(
        filter_contaminants(e, make_row_data(e),
                            cfg_with(remove_contaminants = FALSE,
                                     contaminant_prefix = "["))
    )
    expect_equal(nrow(kept$expr_mat), nrow(e))
    expect_equal(kept$n_removed, 0L)

    res <- suppressMessages(
        filter_contaminants(e, make_row_data(e), cfg_with(contaminant_prefix = "["))
    )
    expect_equal(res$n_removed, 0)
    expect_equal(nrow(res$expr_mat), nrow(e))
})

test_that("a NULL row_data is tolerated", {
    e <- make_expr()
    res <- filter_contaminants(e, NULL, cfg_with(remove_contaminants = FALSE))
    expect_null(res$row_data)
    expect_equal(res$n_removed, 0L)
})
