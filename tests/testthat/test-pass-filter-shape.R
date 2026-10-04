# tests/testthat/test-pass-filter-shape.R
#
# Shape contract of pass_filter().
#
# passes_per_group was built with sapply(), where each group contributes one
# logical per feature. With exactly one feature every group returns a single
# value, sapply() simplifies to a plain vector, and the rowSums() that counts
# groups per feature aborted with "'x' must be an array of at least two
# dimensions". One surviving feature after contaminant removal reaches it.
#
# pass_filter() is defined twice while #142 is open -- in
# R/domain/proteomics/02_filtering.R and R/domain/rnaseq/02_filtering.R -- and
# the domain files are sourced into one environment in sorted order, so the
# rnaseq definition shadows the proteomics one. Testing whatever is in scope
# would therefore exercise only one of the two copies. Every behavioural test
# below runs against both definitions, read from their own files, and the last
# test pins that the two remain the same function.
#
# All fixtures are synthetic (f1..f4, S1..S6, groups A/B/C).

# Read one top-level `<name> <- function(...)` definition out of a source file.
# keep.source = FALSE so deparse() below compares code rather than srcrefs,
# whatever the session's keep.source option happens to be.
read_definition <- function(rel_parts, name) {
    f <- do.call(testthat::test_path, c(list("..", ".."), as.list(rel_parts)))
    if (!file.exists(f)) {
        testthat::skip(paste(basename(f), "not found from the test working directory"))
    }
    exprs <- as.list(parse(f, keep.source = FALSE))
    hits <- Filter(function(e) {
        is.call(e) && length(e) == 3L &&
            as.character(e[[1]]) %in% c("<-", "=") &&
            identical(as.character(e[[2]]), name)
    }, exprs)
    # A second definition in the same file would make "the two copies" ambiguous.
    expect_length(hits, 1L)
    # globalenv() as parent: the function looks its helpers up exactly as it
    # does when the pipeline sources it.
    eval(hits[[1]][[3]], envir = new.env(parent = globalenv()))
}

pass_filter_defs <- function() {
    list(
        proteomics = read_definition(c("R", "domain", "proteomics", "02_filtering.R"), "pass_filter"),
        rnaseq     = read_definition(c("R", "domain", "rnaseq", "02_filtering.R"), "pass_filter")
    )
}

# Run one expectation against both copies, naming the copy when it fails.
for_each_def <- function(f) {
    defs <- pass_filter_defs()
    for (nm in names(defs)) f(defs[[nm]], nm)
}

pf_matrix <- function(values, features, samples) {
    matrix(values, nrow = length(features), byrow = TRUE,
           dimnames = list(features, samples))
}

samples6 <- c("S1", "S2", "S3", "S4", "S5", "S6")
groups6  <- c("A", "A", "B", "B", "C", "C")


# =============================================================================
# One feature
# =============================================================================

test_that("one feature and several groups returns one verdict", {
    # The regression case: under sapply() this aborted in rowSums().
    expr <- pf_matrix(c(1, 2, 3, NA, NA, NA), "f1", samples6)

    for_each_def(function(pass_filter, nm) {
        expect_no_error(keep <- pass_filter(expr, groups6, min_per_group = 2, min_groups = 1))
        expect_length(keep, 1L)
        expect_type(keep, "logical")
        expect_named(keep, "f1")
        expect_true(keep[["f1"]], info = nm)   # A has 2 observations, B has 1, C none
    })
})

test_that("one feature and one group returns one verdict", {
    expr <- pf_matrix(c(1, 2, NA), "f1", c("S1", "S2", "S3"))
    grp <- c("A", "A", "A")

    for_each_def(function(pass_filter, nm) {
        expect_true(pass_filter(expr, grp, min_per_group = 2)[["f1"]], info = nm)
        expect_false(pass_filter(expr, grp, min_per_group = 3)[["f1"]], info = nm)
    })
})

test_that("min_groups counts groups, not observations, for one feature", {
    # f1 is observed twice in A, once in B and not at all in C. At a threshold
    # of one observation it clears exactly two of the three groups, so the
    # min_groups cutoff falls between 2 and 3.
    expr <- pf_matrix(c(1, 2, 3, NA, NA, NA), "f1", samples6)

    for_each_def(function(pass_filter, nm) {
        expect_true(pass_filter(expr, groups6, min_per_group = 1, min_groups = 2)[["f1"]], info = nm)
        expect_false(pass_filter(expr, groups6, min_per_group = 1, min_groups = 3)[["f1"]], info = nm)
        # At a threshold of two observations only A qualifies, so asking for two
        # groups now fails on the same matrix.
        expect_true(pass_filter(expr, groups6, min_per_group = 2, min_groups = 1)[["f1"]], info = nm)
        expect_false(pass_filter(expr, groups6, min_per_group = 2, min_groups = 2)[["f1"]], info = nm)
    })
})

test_that("a per-group threshold applies to the right group for one feature", {
    # Named thresholds must not be reassigned by position when there is a
    # single feature to simplify over.
    expr <- pf_matrix(c(1, 2, 3, NA, NA, NA), "f1", samples6)
    thr <- c(A = 3, B = 1, C = 1)   # A now fails, B passes, C fails

    for_each_def(function(pass_filter, nm) {
        expect_true(pass_filter(expr, groups6, min_per_group = thr, min_groups = 1)[["f1"]], info = nm)
        expect_false(pass_filter(expr, groups6, min_per_group = thr, min_groups = 2)[["f1"]], info = nm)
    })
})


# =============================================================================
# The shapes that already worked
# =============================================================================

test_that("several features and one group are unchanged", {
    expr <- pf_matrix(c(1, 2, 3,
                        1, NA, NA,
                        NA, NA, NA),
                      c("f1", "f2", "f3"), c("S1", "S2", "S3"))
    grp <- c("A", "A", "A")

    for_each_def(function(pass_filter, nm) {
        expect_equal(unname(pass_filter(expr, grp, min_per_group = 2)),
                     c(TRUE, FALSE, FALSE), info = nm)
    })
})

test_that("several features and several groups are unchanged", {
    expr <- pf_matrix(c(1, 2, 3, 4, 5, 6,        # f1: observed everywhere
                        1, 2, 3, NA, NA, NA,     # f2: A fully, B once, C never
                        1, NA, NA, NA, NA, NA,   # f3: A once only
                        NA, NA, NA, NA, NA, NA), # f4: never observed
                      c("f1", "f2", "f3", "f4"), samples6)

    for_each_def(function(pass_filter, nm) {
        keep <- pass_filter(expr, groups6, min_per_group = 2, min_groups = 1)
        expect_equal(unname(keep), c(TRUE, TRUE, FALSE, FALSE), info = nm)
        expect_named(keep, c("f1", "f2", "f3", "f4"))

        # Asking for two groups at that threshold keeps only the complete one.
        expect_equal(unname(pass_filter(expr, groups6, min_per_group = 2, min_groups = 2)),
                     c(TRUE, FALSE, FALSE, FALSE), info = nm)
    })
})

test_that("group order follows first appearance, not sorting", {
    # min_groups only counts, so the verdict cannot depend on group order --
    # this pins that the explicit shape did not reorder the columns either.
    expr <- pf_matrix(c(1, 2, 3, 4,
                        NA, NA, 3, 4),
                      c("f1", "f2"), c("S1", "S2", "S3", "S4"))
    grp <- c("B", "B", "A", "A")

    for_each_def(function(pass_filter, nm) {
        expect_equal(unname(pass_filter(expr, grp, min_per_group = 2, min_groups = 2)),
                     c(TRUE, FALSE), info = nm)
        expect_equal(unname(pass_filter(expr, grp, min_per_group = 2, min_groups = 1)),
                     c(TRUE, TRUE), info = nm)
    })
})


# =============================================================================
# Alone and embedded must agree
# =============================================================================

test_that("a feature filtered alone gets the verdict it gets in a larger matrix", {
    # The strongest statement of the contract, and the one only the sapply()
    # simplification could have broken: a feature's verdict depends on that
    # feature, never on how many others were filtered alongside it.
    expr <- pf_matrix(c(1, 2, 3, 4, 5, 6,
                        1, 2, 3, NA, NA, NA,
                        1, NA, NA, NA, NA, NA,
                        NA, NA, NA, NA, NA, NA),
                      c("f1", "f2", "f3", "f4"), samples6)

    for_each_def(function(pass_filter, nm) {
        for (min_groups in 1:2) {
            many <- pass_filter(expr, groups6, min_per_group = 2, min_groups = min_groups)
            alone <- vapply(rownames(expr), function(id) {
                pass_filter(expr[id, , drop = FALSE], groups6,
                            min_per_group = 2, min_groups = min_groups)[[id]]
            }, logical(1))
            expect_equal(alone, many, info = sprintf("%s, min_groups %d", nm, min_groups))
        }
    })
})


# =============================================================================
# The two copies stay in step
# =============================================================================

test_that("the proteomics and rnaseq definitions are the same function", {
    # While #142 is open pass_filter() lives in two files and the rnaseq copy
    # wins by sort order, so a fix applied to one alone would either do nothing
    # or split the two modes' filtering. Comparing the deparsed body and formals
    # compares code: comments are already gone at parse time, so the two copies
    # may still be documented differently, and neither closure's environment
    # enters the comparison.
    defs <- pass_filter_defs()
    expect_identical(deparse(body(defs$proteomics)), deparse(body(defs$rnaseq)))
    expect_identical(deparse(formals(defs$proteomics)), deparse(formals(defs$rnaseq)))
})
