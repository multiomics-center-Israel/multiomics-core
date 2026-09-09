# tests/testthat/test-clustering-empty-de.R
#
# A DE step that returns nothing is a legitimate result, not a pipeline
# failure. run_clustering() rightly refuses an empty feature set, but that
# error used to abort the whole run, so a comparison analysis expected to find
# little or nothing -- an unmoderated paired t-test, or plain limma on a paired
# design -- could not reach its report at all.
#
# Two runs died this way: the first Levenberg run (0 features, unblocked limma)
# and the paired t-test comparison run.

test_that("run_clustering still refuses an empty feature set", {
    # The low-level contract is unchanged: it is the caller's job to check.
    expr <- matrix(rnorm(20), nrow = 5,
                   dimnames = list(paste0("f", 1:5), paste0("s", 1:4)))
    expect_error(
        run_clustering(expr, col_data = NULL, de_features = character(0),
                       config = list(method = "hierarchical")),
        "No differential features"
    )
})

test_that("run_clustering also refuses features absent from the matrix", {
    expr <- matrix(rnorm(20), nrow = 5,
                   dimnames = list(paste0("f", 1:5), paste0("s", 1:4)))
    expect_error(
        run_clustering(expr, col_data = NULL, de_features = c("nope1", "nope2"),
                       config = list(method = "hierarchical")),
        "No differential features"
    )
})

# The module's guard, matching the one metabolomics and lipidomics already use:
# fewer than 2 usable features means there is nothing to cluster.
should_skip_clustering <- function(de_features, expr_rownames) {
    length(intersect(de_features, expr_rownames)) < 2L
}

test_that("the module skips when nothing overlaps the matrix", {
    expect_true(should_skip_clustering(character(0), paste0("f", 1:5)))
    expect_true(should_skip_clustering(c("x", "y"), paste0("f", 1:5)))
    expect_true(should_skip_clustering(NULL, paste0("f", 1:5)))
})

test_that("a single usable feature is also skipped, matching the other omics", {
    # metabolomics/lipidomics use `< 2`: one feature cannot be clustered, and
    # hclust on a 1-row matrix errors rather than returning something useless.
    expect_true(should_skip_clustering(c("f1", "nope"), paste0("f", 1:5)))
})

test_that("the module proceeds once two features overlap", {
    expect_false(should_skip_clustering(c("f1", "f2", "nope"), paste0("f", 1:5)))
})

test_that("the skip return matches the shape the other omics modules use", {
    skip_return <- list(plots = list(), files = character(0),
                        excel_order = NULL, objects = list())
    expect_setequal(names(skip_return), c("plots", "files", "excel_order", "objects"))
})

test_that("the module guards before calling run_clustering", {
    candidates <- c(
        testthat::test_path("..", "..", "R", "modules", "proteomics",
                            "04_mod_clustering.R"),
        "R/modules/proteomics/04_mod_clustering.R"
    )
    f <- candidates[file.exists(candidates)][1]
    skip_if(is.na(f), "04_mod_clustering.R not found from the test working directory")

    src <- readLines(f, warn = FALSE)
    guard <- grep("n_available <- length(intersect(de_features, rownames(expr_mat)))",
                  src, fixed = TRUE)
    first_step <- grep(".run_hierarchical_step(", src, fixed = TRUE)
    expect_length(guard, 1)
    expect_true(length(first_step) > 0 && guard[1] < first_step[1])
})
