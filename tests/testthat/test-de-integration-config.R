# tests/testthat/test-de-integration-config.R
#
# The DE-integration config contract: every mistake is reported at config time,
# naming the key to fix, and the defaults the readers rely on are filled in.

two_layers <- function() {
    list(
        list(name = "cells", omics_type = "proteomics",
             format = "proteomics_summary", path = "/tmp/cells.tsv"),
        list(name = "media", omics_type = "proteomics",
             format = "proteomics_final", path = "/tmp/media.tsv"))
}

test_that("defaults are filled and labels default to the layer name", {
    cfg <- validate_de_integration_config(list(layers = two_layers()))
    expect_true(cfg$hits$use_table_flag)
    expect_equal(cfg$hits$p_cutoff, 0.05)
    expect_true(cfg$hits$use_adjusted)
    expect_equal(cfg$hits$linear_fc_cutoff, 1.5)
    expect_equal(cfg$concordance$well_observed_min, 2)
    expect_identical(cfg$layers[[1]]$label, "cells")
    expect_length(cfg$comparisons, 0)
})

test_that("one layer is not enough to compare", {
    expect_error(validate_de_integration_config(list(layers = two_layers()[1])),
                 "at least two layers")
})

test_that("layer names must be one lower-case word, and unique", {
    ly <- two_layers(); ly[[1]]$name <- "Cells media"
    expect_error(validate_de_integration_config(list(layers = ly)), "one lower-case word")
    ly <- two_layers(); ly[[2]]$name <- "cells"
    expect_error(validate_de_integration_config(list(layers = ly)), "unique")
})

test_that("omics type and format are checked against each other", {
    ly <- two_layers(); ly[[1]]$omics_type <- "metabolomics"
    expect_error(validate_de_integration_config(list(layers = ly)), "omics_type")
    ly <- two_layers(); ly[[1]]$format <- "rnaseq_summary"
    expect_error(validate_de_integration_config(list(layers = ly)), "not a proteomics export")
    ly <- two_layers(); ly[[1]]$format <- "excel"
    expect_error(validate_de_integration_config(list(layers = ly)), "format must be one of")
})

test_that("a generic layer must map its id, p-value and fold change", {
    ly <- two_layers()
    ly[[1]]$format <- "generic"
    ly[[1]]$columns <- list(id = "protein")
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "pvalue, log2fc \\(or linear_fc\\)")
    ly[[1]]$columns <- list(id = "protein", pvalue = "P", linear_fc = "FC")
    expect_silent(validate_de_integration_config(list(layers = ly)))
})

test_that("an observed block must be complete", {
    ly <- two_layers()
    ly[[1]]$observed <- list(matrix = "m.tsv", samplesheet = "s.csv")
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "sample_col, contrasts_file")
})

test_that("comparisons must name known layers, and flip only their members", {
    cfg <- list(layers = two_layers(), comparisons = list(
        list(name = "c1", members = list(cells = "A_vs_B", tissue = "A_vs_B"))))
    expect_error(validate_de_integration_config(cfg), "unknown layer\\(s\\): tissue")

    cfg$comparisons[[1]]$members <- list(cells = "A_vs_B", media = "A_vs_B")
    cfg$comparisons[[1]]$flip <- list("tissue")
    expect_error(validate_de_integration_config(cfg), "flip names layer\\(s\\) not in")

    cfg$comparisons[[1]]$flip <- list("media")
    expect_identical(validate_de_integration_config(cfg)$comparisons[[1]]$flip, "media")
})

test_that("hit cutoffs are range-checked", {
    expect_error(validate_de_integration_config(
        list(layers = two_layers(), hits = list(p_cutoff = 0))), "p_cutoff")
    expect_error(validate_de_integration_config(
        list(layers = two_layers(), hits = list(linear_fc_cutoff = 0.5))), ">= 1")
})

test_that("validate_config() hands the section to the validator and keeps its defaults", {
    cfg <- validate_config(list(paths = list(raw = "data"), params = list(seed = 1),
                                modes = list(de_integration = list(layers = two_layers()))))
    expect_equal(cfg$modes$de_integration$hits$p_cutoff, 0.05)
})

test_that("the shipped template parses and validates", {
    f <- file.path(root_dir, "config", "templates", "de_integration_config.yaml")
    skip_if_not(file.exists(f))
    cfg <- yaml::read_yaml(f)
    expect_true(all(c("project", "paths", "params", "modes") %in% names(cfg)))
    out <- validate_de_integration_config(cfg$modes$de_integration)
    expect_length(out$layers, 2)
})
