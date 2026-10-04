# tests/testthat/test-load-gene-sets-null-gmt.R
#
# Regression tests for load_gene_sets() with no custom GMT (R/core/09_enrichment.R).
# `gmt_file: null` is the normal config for a model organism using built-in
# GO/KEGG. unlist(NULL) stays NULL and file.exists(NULL) raises
# "invalid 'file' argument", so the GMT-path branch used to abort before the
# built-in databases were ever reached. The caller wraps load_gene_sets() in
# tryCatch, so this surfaced as a silent "No gene sets loaded" rather than a
# failure — hence a test on the path-normalisation itself.

normalize_gmt_paths <- function(gmt_file) {
    # Mirrors the first two lines of the GMT branch in load_gene_sets().
    p <- as.character(unlist(gmt_file, use.names = FALSE))
    p[nzchar(p)]
}

test_that("NULL gmt_file normalises to an empty character vector", {
    expect_identical(normalize_gmt_paths(NULL), character(0))
})

test_that("file.exists tolerates the normalised empty result", {
    # The bug: file.exists(NULL) errors instead of returning logical(0).
    expect_error(file.exists(NULL), "invalid 'file' argument")
    expect_identical(file.exists(normalize_gmt_paths(NULL)), logical(0))
})

test_that("empty and blank gmt_file entries are dropped, not passed through", {
    expect_identical(normalize_gmt_paths(list()), character(0))
    expect_identical(normalize_gmt_paths(""), character(0))
    expect_identical(normalize_gmt_paths(c("", "a.gmt")), "a.gmt")
})

test_that("real gmt paths survive normalisation, scalar or list", {
    expect_identical(normalize_gmt_paths("go.gmt"), "go.gmt")
    expect_identical(normalize_gmt_paths(list("go.gmt", "kegg.gmt")),
                     c("go.gmt", "kegg.gmt"))
})

test_that("load_gene_sets returns built-in collections when gmt_file is NULL", {
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("AnnotationDbi")
    skip_on_cran()

    gs <- load_gene_sets(organism = "Homo sapiens",
                         pathway_database = "GO_BP",
                         gmt_file = NULL,
                         annotation = NULL,
                         target_id_type = "symbol")

    expect_gt(length(gs), 0)
    expect_true(any(grepl("^GO", names(gs))))
})
