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

test_that("a comparison name that is only spaces is rejected", {
    cfg <- list(layers = two_layers(), comparisons = list(
        list(name = "   ", members = list(cells = "A_vs_B", media = "A_vs_B"))))
    expect_error(validate_de_integration_config(cfg), "comparisons\\[\\[1\\]\\].name is required")
})

test_that("comparison members must map each layer to one contrast label", {
    cfg <- list(layers = two_layers(), comparisons = list(
        list(name = "c1", members = list("A_vs_B", "A_vs_B"))))
    expect_error(validate_de_integration_config(cfg), "entries without a layer name")

    cfg$comparisons[[1]]$members <- list(cells = "A_vs_B", "A_vs_B")
    expect_error(validate_de_integration_config(cfg), "entries without a layer name")

    cfg$comparisons[[1]]$members <- list(cells = "A_vs_B", cells = "C_vs_D")
    expect_error(validate_de_integration_config(cfg), "more than once: cells")

    cfg$comparisons[[1]]$members <- list(cells = c("A_vs_B", "C_vs_D"), media = "A_vs_B")
    expect_error(validate_de_integration_config(cfg), "single non-empty label for: cells")

    cfg$comparisons[[1]]$members <- list(cells = list("A_vs_B"), media = "")
    expect_error(validate_de_integration_config(cfg), "label for: cells, media")

    cfg$comparisons[[1]]$members <- list(cells = "A_vs_B", media = NA_character_)
    expect_error(validate_de_integration_config(cfg), "label for: media")

    cfg$comparisons[[1]]$members <- c(cells = "A_vs_B", media = "A_vs_B")
    expect_error(validate_de_integration_config(cfg), "must name a contrast")

    cfg$comparisons[[1]]$members <- list(cells = "A_vs_B", media = "A vs B")
    expect_silent(validate_de_integration_config(cfg))
})

test_that("id_col is one column name, for native formats only", {
    ly <- two_layers(); ly[[1]]$id_col <- "Protein.Group"
    expect_identical(validate_de_integration_config(list(layers = ly))$layers[[1]]$id_col,
                     "Protein.Group")
    ly[[1]]$id_col <- c("a", "b")
    expect_error(validate_de_integration_config(list(layers = ly)), "id_col must be one")
    ly[[1]]$id_col <- ""
    expect_error(validate_de_integration_config(list(layers = ly)), "id_col must be one")

    ly <- two_layers()
    ly[[1]]$format <- "generic"
    ly[[1]]$columns <- list(id = "protein", pvalue = "P", log2fc = "logFC")
    ly[[1]]$id_col <- "protein"
    expect_error(validate_de_integration_config(list(layers = ly)), "columns.id")
})

test_that("a generic layer maps both observed-count columns, or neither", {
    ly <- two_layers()
    ly[[1]]$format <- "generic"
    ly[[1]]$columns <- list(id = "protein", pvalue = "P", log2fc = "logFC", n_obs_num = "nA")
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "only one of n_obs_num and n_obs_den")
    ly[[1]]$columns$n_obs_den <- "nB"
    expect_silent(validate_de_integration_config(list(layers = ly)))
})

test_that("hit cutoffs are range-checked", {
    expect_error(validate_de_integration_config(
        list(layers = two_layers(), hits = list(p_cutoff = 0))), "p_cutoff")
    expect_error(validate_de_integration_config(
        list(layers = two_layers(), hits = list(linear_fc_cutoff = 0.5))), ">= 1")
})

test_that("a misspelt key anywhere in the section is an error, not a silent default", {
    base <- list(layers = two_layers())
    expect_error(validate_de_integration_config(c(base, list(hits = list(p_cutof = 0.01)))),
                 "de_integration.hits has unknown key\\(s\\): p_cutof")
    expect_error(validate_de_integration_config(c(base, list(comparison = list()))),
                 "de_integration has unknown key\\(s\\): comparison")
    expect_error(validate_de_integration_config(c(base, list(
        concordance = list(well_observed = 3)))),
        "concordance has unknown key\\(s\\): well_observed")

    ly <- two_layers(); ly[[2]]$annotaton_file <- "ann.csv"
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "layers\\[\\[2\\]\\] has unknown key\\(s\\): annotaton_file")

    ly <- two_layers()
    ly[[1]]$observed <- list(matrix = "m.tsv", samplesheet = "s.csv", sample_col = "S",
                             contrasts_file = "c.csv", idcol = "FeatureID")
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "observed has unknown key\\(s\\): idcol")
    ly[[1]]$observed <- list(matrix = "m.tsv", samplesheet = "s.csv", sample_col = "",
                             contrasts_file = "c.csv")
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "observed.sample_col must each be one non-empty string")

    cfg <- list(layers = two_layers(), comparisons = list(list(
        name = "c1", members = list(cells = "A_vs_B", media = "A_vs_B"), flipp = "media")))
    expect_error(validate_de_integration_config(cfg), "unknown key\\(s\\): flipp")
    cfg$comparisons[[1]]$flipp <- NULL
    cfg$comparisons[[1]]$flip <- list(1)
    expect_error(validate_de_integration_config(cfg), "flip must list layer names")
})

test_that("contrast and columns belong to generic tables, and contrast is one label", {
    ly <- two_layers(); ly[[1]]$contrast <- "A_vs_B"
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "contrast is only read for a \"generic\" table")
    ly <- two_layers(); ly[[1]]$columns <- list(id = "x")
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "columns is only read for a \"generic\" table")

    ly <- two_layers()
    ly[[1]]$format <- "generic"
    ly[[1]]$columns <- list(id = "protein", pvalue = "P", log2fc = "logFC")
    ly[[1]]$contrast <- c("A_vs_B", "C_vs_D")
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "contrast must be one non-empty label")
    ly[[1]]$contrast <- ""
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "contrast must be one non-empty label")
    ly[[1]]$contrast <- "   "
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "contrast must be one non-empty label")
    ly[[1]]$contrast <- "A_vs_B"
    expect_silent(validate_de_integration_config(list(layers = ly)))
})

test_that("well_observed_min is a number", {
    expect_error(validate_de_integration_config(list(layers = two_layers(),
        concordance = list(well_observed_min = "2"))), "well_observed_min must be a finite number")
    expect_error(validate_de_integration_config(list(layers = two_layers(),
        concordance = list(well_observed_min = -1))), "well_observed_min")
})

test_that("cutoffs and thresholds must be finite", {
    base <- list(layers = two_layers())
    expect_error(validate_de_integration_config(c(base, list(hits = list(linear_fc_cutoff = Inf)))),
                 "linear_fc_cutoff .*finite number >= 1")
    expect_error(validate_de_integration_config(c(base, list(
        concordance = list(well_observed_min = Inf)))), "well_observed_min must be a finite")
    y <- yaml::yaml.load("hits: {linear_fc_cutoff: .inf}")
    expect_error(validate_de_integration_config(c(base, y)), "linear_fc_cutoff")
})

test_that("a layer's hit override is checked as the settings it will use", {
    ly <- two_layers(); ly[[1]]$hits <- list(linear_fc_cutoff = 0.5)
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "layers\\[\\[1\\]\\] \\('cells'\\).hits.linear_fc_cutoff .*>= 1")
    ly[[1]]$hits <- list(p_cutoff = 2)
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "\\('cells'\\).hits.p_cutoff must be a number in \\(0, 1\\]")
    ly[[1]]$hits <- list(use_adjusted = "yes")
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "\\('cells'\\).hits.use_adjusted must be true or false")
    ly[[1]]$hits <- list(p_cutof = 0.01)
    expect_error(validate_de_integration_config(list(layers = ly)), "unknown key\\(s\\): p_cutof")
    ly[[1]]$hits <- list(p_cutoff = 0.01, use_table_flag = FALSE)
    expect_silent(validate_de_integration_config(list(layers = ly)))
})

test_that("every column-map entry of a generic layer is one column name", {
    ly <- two_layers()
    ly[[1]]$format <- "generic"
    ly[[1]]$columns <- list(id = "protein", pvalue = NULL, log2fc = "logFC")
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "columns.pvalue must each be one column name")
    ly[[1]]$columns <- list(id = "protein", pvalue = "P", log2fc = "logFC", padj = "")
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "columns.padj must each be one column name")
    ly[[1]]$columns <- list(id = "protein", pvalue = "P", log2fc = c("a", "b"))
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "columns.log2fc must each be one column name")
    ly[[1]]$columns <- list(id = "protein", pvalue = "P", log2fc = "logFC", padjj = "q")
    expect_error(validate_de_integration_config(list(layers = ly)),
                 "unknown key\\(s\\): padjj")
})

test_that("a blank YAML mapping is caught, not dropped", {
    y <- yaml::yaml.load("columns: {id: protein, pvalue: P, log2fc: logFC, padj: }")
    ly <- two_layers()
    ly[[1]]$format <- "generic"
    ly[[1]]$columns <- y$columns
    expect_error(validate_de_integration_config(list(layers = ly)), "columns.padj")
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
