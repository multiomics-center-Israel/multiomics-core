# tests/testthat/test-report-explorer-missing.R
#
# Three defects in the proteomics report, all in code that only exists inside
# report_template_proteomics.Rmd:
#
# 1. The protein explorer averaged missing values as zeros. The explorer reads
#    the measured `<sample>` block of final_results.tsv (not the imputed
#    `<sample>.norm` block), so missing intensities are NA in R and null in the
#    JSON the template emits. The info line filtered them with `!isNaN(v)`, and
#    isNaN(null) is false because Number(null) is 0 -- so every missing value
#    survived the filter, was added into the sum as a real zero, and was counted
#    toward the reported sample count.
#
# 2. A group with no measured values disappeared. Its trace held only nulls,
#    which Plotly drops, and an empty y-only box trace contributes no category,
#    so the group vanished from the axis instead of being shown as unmeasured.
#
# 3. ORA section headings stripped `_up`/`_down` along with the file suffix, so
#    a database with both files produced two sections headed identically.
#
# These are assertions against the template source, not behavioural tests. The
# first two live in browser JavaScript that this suite cannot execute, and the
# third runs only inside a knitr chunk during a render. Source-level pinning is
# the same idiom used by test-imputation-seeding.R for the preprocessing seed.
# The render checklist in the PR is what actually verifies the output.

# Locate a repo file from either the repo root or tests/testthat.
repo_file <- function(...) {
    candidates <- c(testthat::test_path("..", "..", ...), file.path(...))
    candidates[file.exists(candidates)][1]
}

report_src <- function() {
    f <- repo_file("R", "domain", "proteomics", "report_template_proteomics.Rmd")
    skip_if(is.na(f) || !file.exists(f), "proteomics report template not found")
    paste(readLines(f, warn = FALSE), collapse = "\n")
}


# =============================================================================
# 1. Missing values are not measured values
# =============================================================================

test_that("the protein explorer has one guard for what counts as measured", {
    src <- report_src()

    expect_match(src, "function bpIsMeasured(v)", fixed = TRUE)
    # The null test is the whole point: !isNaN() alone is what let nulls in.
    expect_match(src, "v !== null", fixed = TRUE)
    expect_match(src, "v !== undefined", fixed = TRUE)
})

test_that("the mean and the sample count exclude missing values", {
    src <- report_src()

    # The exact filter that produced the bug, gone.
    expect_false(grepl("filter(function(v) { return !isNaN(v); })", src, fixed = TRUE))
    expect_match(src, "values.filter(bpIsMeasured)", fixed = TRUE)
    # The count is reported as measured samples, so a reader can tell that a
    # protein with gaps is summarised over fewer samples than the run has.
    expect_match(src, "measured samples", fixed = TRUE)
    # A protein with nothing measured must not print a mean of NaN.
    expect_match(src, "no measured values", fixed = TRUE)
})

test_that("both explorer traces drop missing values with that same guard", {
    src <- report_src()

    # Single-protein branch: value and sample name stay in lockstep, because the
    # hover text is positional and would otherwise mislabel points.
    expect_match(src, "if (groups[i] === g && bpIsMeasured(values[i])) {", fixed = TRUE)
    # Every use: the definition, the two trace loops, the empty-group scan, and
    # the info-line filter -- which names the function without calling it, so
    # this counts the bare name. Catches a guard being dropped again without
    # blocking a new one being added.
    n_uses <- lengths(regmatches(src, gregexpr("bpIsMeasured", src, fixed = TRUE)))
    expect_gte(n_uses, 5L)
})

test_that("the peptide explorer keeps the guard this fix was modelled on", {
    # The protein explorer now matches what the peptide explorer already did.
    # If that reference is ever weakened, this fix loses its precedent.
    expect_match(report_src(), "pep.expr[i] !== null && !isNaN(pep.expr[i])", fixed = TRUE)
})


# =============================================================================
# 2. A group with nothing measured stays visible
# =============================================================================

test_that("every group keeps its place on the axis", {
    src <- report_src()

    # An empty trace contributes no category, so the axis is seeded from the
    # group list itself -- in first-appearance order, which is the order
    # uniqueGroups is built in and the order the plot has always used.
    expect_match(src, "categoryorder: \"array\"", fixed = TRUE)
    expect_match(src, "categoryarray: uniqueGroups", fixed = TRUE)
})

test_that("a group with nothing measured is named rather than left blank", {
    src <- report_src()

    expect_match(src, "no measured values in: ", fixed = TRUE)
    expect_match(src, "emptyGroups", fixed = TRUE)
})


# =============================================================================
# 3. ORA headings keep their direction
# =============================================================================

test_that("ORA sections name their direction", {
    src <- report_src()

    expect_match(src, "_ora_up\\\\.csv$", fixed = TRUE)
    expect_match(src, "_ora_down\\\\.csv$", fixed = TRUE)
    expect_match(src, "(up-regulated)", fixed = TRUE)
    expect_match(src, "(down-regulated)", fixed = TRUE)
})

test_that("ORA file discovery and statistics are untouched", {
    src <- report_src()

    # The fix renames sections only: same files read, same significance rule.
    expect_match(src, "pattern = \"pathway_.*_ora.*\\\\.csv$\"", fixed = TRUE)
    expect_match(src, "ora_df[[padj_col]] < 0.05", fixed = TRUE)
})
