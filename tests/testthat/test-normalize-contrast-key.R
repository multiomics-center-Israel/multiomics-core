# tests/testthat/test-normalize-contrast-key.R
#
# normalize_contrast_key() decides which contrast labels name the same
# comparison, across modes: the multiomics per-contrast joins and the
# DE-integration pairing both rely on it. Synthetic labels only.

test_that("the documented spellings of one comparison share a key", {
    expect_identical(normalize_contrast_key(c("A vs. B", "A_vs_B", "A - B", "a - b",
                                              "A vs B", "A-B")),
                     rep("avsb", 6L))
    # The space-stripped proteomics/ORA spelling agrees with the spaced one.
    expect_identical(normalize_contrast_key("1.56ppmvs.0ppm"),
                     normalize_contrast_key("1.56ppm vs. 0ppm"))
    expect_identical(normalize_contrast_key("1.56ppm - 0ppm"),
                     normalize_contrast_key("1.56ppm_vs_0ppm"))
})

test_that("a decimal point inside a number is content", {
    expect_false(identical(normalize_contrast_key("1.56ppm vs 0ppm"),
                           normalize_contrast_key("15.6ppm vs 0ppm")))
    # make.names padding and "vs." stay punctuation.
    expect_identical(normalize_contrast_key("Day.1_vs_Day.0"),
                     normalize_contrast_key("Day 1 vs Day 0"))
})

test_that("a hyphen inside a group name is not read as the separator", {
    k <- normalize_contrast_key(c("A-B_vs_C", "A_vs_B-C"))
    expect_false(identical(k[1], k[2]))
    # A hyphenated group written with or without the spaces around "vs" agrees.
    expect_identical(normalize_contrast_key("WT-1_vs_KO"), normalize_contrast_key("WT-1vsKO"))
    # A bare hyphen is still a comparison, not a group name with a space.
    expect_false(identical(normalize_contrast_key("A-B"), normalize_contrast_key("A B")))
})

test_that("a leading decimal point keeps its number", {
    expect_false(identical(normalize_contrast_key(".5ppm_vs_0ppm"),
                           normalize_contrast_key("5ppm_vs_0ppm")))
    expect_identical(normalize_contrast_key(".5ppm_vs_0ppm"),
                     normalize_contrast_key("0.5ppm_vs_0ppm"))
    expect_identical(normalize_contrast_key("A vs .5ppm"),
                     normalize_contrast_key("A vs 0.5ppm"))
})

test_that("missing and empty input keep their shape", {
    expect_identical(normalize_contrast_key(character(0)), character(0))
    expect_identical(normalize_contrast_key(c("A vs B", NA)), c("avsb", NA))
})
