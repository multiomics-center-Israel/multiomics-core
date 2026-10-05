# tests/testthat/test-de-subset-contrasts.R
#
# Unit tests for the combined group x subset-factor design helpers used by
# run_limma_proteomics() when de$subset_factor is set:
#   - build_group_subset_factor()
#   - build_subset_contrast_formulas()
#
# Synthetic data mirroring the Noam_Ziv layout: Group {perfused, non_perfused}
# crossed with Month {Feb2026 (2/2), July2025 (1/1), May2026 (1/1)}.

make_meta <- function() {
    data.frame(
        SampleName = paste0("S", 1:8),
        Group = c("perfused", "non_perfused", "perfused", "perfused",
                  "non_perfused", "non_perfused", "non_perfused", "perfused"),
        Month = c("July2025", "July2025", "Feb2026", "Feb2026",
                  "Feb2026", "Feb2026", "May2026", "May2026"),
        stringsAsFactors = FALSE
    )
}

test_that("build_group_subset_factor produces one make.names-safe level per observed cell", {
    meta <- make_meta()
    combo <- build_group_subset_factor(meta, "Group", "Month")

    expect_s3_class(combo, "factor")
    expect_length(combo, nrow(meta))
    # 2 groups x 3 months, all combinations observed
    expect_setequal(
        levels(combo),
        make.names(c("perfused__Feb2026", "perfused__July2025", "perfused__May2026",
                     "non_perfused__Feb2026", "non_perfused__July2025", "non_perfused__May2026"))
    )
    # Levels are valid R names (usable as design column names)
    expect_identical(levels(combo), make.names(levels(combo)))
})

test_that("blank Months averages over all subset levels", {
    meta <- make_meta()
    levs <- levels(build_group_subset_factor(meta, "Group", "Month"))
    contrasts_df <- data.frame(
        Contrast_name = "All", Factor = "Group",
        Numerator = "perfused", Denominator = "non_perfused",
        Months = "", stringsAsFactors = FALSE
    )

    f <- build_subset_contrast_formulas(contrasts_df, meta, "Group", "Month", levs)
    expect_named(f, "All")
    # Each side averages 3 cells
    expect_match(f[["All"]], "/3 - .*/3")
    expect_match(f[["All"]], "perfused__Feb2026")
    expect_false(grepl("non_perfused__Feb2026.*-.*non_perfused", f[["All"]]) &&
                     !grepl("non_perfused__Feb2026", f[["All"]]))
})

test_that("specific Months restricts and averages only those cells", {
    meta <- make_meta()
    levs <- levels(build_group_subset_factor(meta, "Group", "Month"))
    contrasts_df <- data.frame(
        Contrast_name = c("Feb", "FebJuly"),
        Factor = "Group",
        Numerator = "perfused", Denominator = "non_perfused",
        Months = c("Feb2026", "Feb2026;July2025"),
        stringsAsFactors = FALSE
    )

    f <- build_subset_contrast_formulas(contrasts_df, meta, "Group", "Month", levs)
    # Feb-only: single cell per side
    expect_equal(unname(f[["Feb"]]),
                 "(perfused__Feb2026)/1 - (non_perfused__Feb2026)/1")
    # Feb+July: two cells per side, averaged
    expect_match(f[["FebJuly"]], "/2 - .*/2")
    expect_match(f[["FebJuly"]], "perfused__July2025")
    expect_false(grepl("May2026", f[["FebJuly"]]))
})

test_that("the generated expressions parse via limma::makeContrasts", {
    meta <- make_meta()
    combo <- build_group_subset_factor(meta, "Group", "Month")
    design <- stats::model.matrix(~ 0 + combo)
    colnames(design) <- levels(combo)

    contrasts_df <- data.frame(
        Contrast_name = c("All", "Feb", "MayJuly"),
        Factor = "Group",
        Numerator = "perfused", Denominator = "non_perfused",
        Months = c("", "Feb2026", "May2026;July2025"),
        stringsAsFactors = FALSE
    )
    f <- build_subset_contrast_formulas(contrasts_df, meta, "Group", "Month",
                                        colnames(design))

    cm <- limma::makeContrasts(contrasts = f, levels = design)
    expect_equal(ncol(cm), 3L)
    # Each contrast sums to zero (balanced numerator vs denominator weights)
    expect_true(all(abs(colSums(cm)) < 1e-8))
})

test_that("unknown subset level and empty cell raise clear errors", {
    meta <- make_meta()
    levs <- levels(build_group_subset_factor(meta, "Group", "Month"))

    bad_month <- data.frame(
        Contrast_name = "X", Factor = "Group",
        Numerator = "perfused", Denominator = "non_perfused",
        Months = "Dec2099", stringsAsFactors = FALSE
    )
    expect_error(
        build_subset_contrast_formulas(bad_month, meta, "Group", "Month", levs),
        "not present in data"
    )

    bad_group <- data.frame(
        Contrast_name = "Y", Factor = "Group",
        Numerator = "ghost", Denominator = "non_perfused",
        Months = "Feb2026", stringsAsFactors = FALSE
    )
    expect_error(
        build_subset_contrast_formulas(bad_group, meta, "Group", "Month", levs),
        "no samples for numerator"
    )
})
