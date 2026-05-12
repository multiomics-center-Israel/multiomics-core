# Tests for config template validity
# Ensures YAML templates parse without errors and contain expected sections

cfg_dir <- file.path(root_dir, "config", "templates")

test_that("proteomics config template is valid YAML with required keys", {
    f <- file.path(cfg_dir, "proteins_config.yaml")
    skip_if_not(file.exists(f))

    cfg <- yaml::read_yaml(f)
    expect_true(is.list(cfg))
    expect_true("project" %in% names(cfg) || "modes" %in% names(cfg))
})

test_that("rna config template is valid YAML with required keys", {
    f <- file.path(cfg_dir, "rna_config.yaml")
    skip_if_not(file.exists(f))

    cfg <- yaml::read_yaml(f)
    expect_true(is.list(cfg))
    expect_true("project" %in% names(cfg) || "modes" %in% names(cfg))
})

test_that("metabolomics config template is valid YAML", {
    f <- file.path(cfg_dir, "metabolomics_template.yaml")
    skip_if_not(file.exists(f))

    cfg <- yaml::read_yaml(f)
    expect_true(is.list(cfg))
})

test_that("multiomics config template is valid YAML with integration section", {
    f <- file.path(cfg_dir, "multiomics_config.yaml")
    skip_if_not(file.exists(f))

    cfg <- yaml::read_yaml(f)
    expect_true(is.list(cfg))

    multi <- cfg$modes$multiomics %||% cfg$multiomics
    if (!is.null(multi)) {
        expect_true("integration" %in% names(multi) || "feature_selection" %in% names(multi))
    }
})
