# tests/testthat/test-pipeline-orchestration.R
# Tests for pipeline orchestration functions (pipe_rnaseq, pipe_proteomics, etc.)
# Validates that each pipe function returns correctly structured lists of tar_target objects
# without evaluating the underlying data (testing DAG creation only).

skip_if_not_installed("targets")
library(targets)

# =============================================================================
# Helper: Extract target names from a list of tar_target objects
# =============================================================================
get_target_names <- function(target_list) {
    vapply(target_list, function(t) t$settings$name, character(1))
}

# =============================================================================
# Helper: Validate that every element in a target list is a proper target object
# =============================================================================
assert_valid_targets <- function(target_list) {
    expect_type(target_list, "list")
    expect_true(length(target_list) > 0, info = "Target list should not be empty")

    for (i in seq_along(target_list)) {
        tgt <- target_list[[i]]
        expect_true(
            !is.null(tgt$settings$name),
            info = sprintf("Element %d should have $settings$name", i)
        )
    }
}


# =============================================================================
# Tests: pipe_rnaseq()
# =============================================================================

test_that("pipe_rnaseq returns a valid target list", {
    result <- pipe_rnaseq()
    assert_valid_targets(result)
})

test_that("pipe_rnaseq contains expected target names", {
    result <- pipe_rnaseq()
    names <- get_target_names(result)
    expect_true("rna_inputs" %in% names)
    expect_true("rna_pre" %in% names)
    expect_true("rna_de_res" %in% names)
    expect_true("rna_out_dir" %in% names)
})

test_that("pipe_rnaseq returns expected number of targets", {
    result <- pipe_rnaseq()
    expect_gte(length(result), 15)
})


# =============================================================================
# Tests: pipe_proteomics()
# =============================================================================

test_that("pipe_proteomics returns a valid target list", {
    result <- pipe_proteomics()
    assert_valid_targets(result)
})

test_that("pipe_proteomics contains expected target names", {
    result <- pipe_proteomics()
    names <- get_target_names(result)
    expect_true("prot_out_dir" %in% names)
    expect_true("prot_pre" %in% names)
    expect_true("prot_de_res" %in% names)
})

test_that("pipe_proteomics returns expected number of targets", {
    result <- pipe_proteomics()
    expect_gte(length(result), 10)
})


# =============================================================================
# Tests: pipe_metabolomics()
# =============================================================================

test_that("pipe_metabolomics returns a valid target list", {
    result <- pipe_metabolomics(chosen_norm = "tss")
    assert_valid_targets(result)
})

test_that("pipe_metabolomics contains expected target names", {
    result <- pipe_metabolomics(chosen_norm = "tss")
    names <- get_target_names(result)
    expect_true("metab_out_dir" %in% names)
    expect_true("metab_pre" %in% names)
    expect_true("metab_de_res" %in% names)
})

test_that("pipe_metabolomics returns expected number of targets", {
    result <- pipe_metabolomics(chosen_norm = "tss")
    expect_gte(length(result), 8)
})


# =============================================================================
# Tests: pipe_multiomics()
# =============================================================================

# Minimal cfg_raw with all three omics modes active, used by pipe_multiomics()
# to determine gateway target wiring at plan-definition time.
mock_cfg_raw <- list(
    modes = list(
        rna          = list(enabled = TRUE),
        proteomics   = list(enabled = TRUE),
        metabolomics = list(enabled = TRUE)
    )
)

test_that("pipe_multiomics returns a valid target list", {
    result <- pipe_multiomics(cfg_raw = mock_cfg_raw)
    assert_valid_targets(result)
})

test_that("pipe_multiomics contains expected target names", {
    result <- pipe_multiomics(cfg_raw = mock_cfg_raw)
    names <- get_target_names(result)
    expect_true("multiomics_out_dir" %in% names)
    expect_true("multiomics_harmonization" %in% names)
    expect_true("multiomics_integration" %in% names)
    expect_true("multiomics_concordance" %in% names)
    expect_true("multiomics_cross_enrichment" %in% names)
    expect_true("multiomics_report" %in% names)
})

test_that("pipe_multiomics returns expected number of targets", {
    result <- pipe_multiomics(cfg_raw = mock_cfg_raw)
    expect_gte(length(result), 10)
})

test_that("multiomics_report target has file format", {
    result <- pipe_multiomics(cfg_raw = mock_cfg_raw)
    names <- get_target_names(result)
    report_idx <- which(names == "multiomics_report")
    expect_length(report_idx, 1)
    expect_equal(result[[report_idx]]$settings$format, "file")
})


# =============================================================================
# Tests: Target name uniqueness across pipes
# =============================================================================

test_that("target names are unique within each pipe", {
    pipes <- list(
        rna = pipe_rnaseq(),
        prot = pipe_proteomics(),
        metab = pipe_metabolomics(chosen_norm = "tss"),
        multi = pipe_multiomics(cfg_raw = mock_cfg_raw)
    )

    for (pipe_name in names(pipes)) {
        target_names <- get_target_names(pipes[[pipe_name]])
        expect_equal(
            length(target_names), length(unique(target_names)),
            info = sprintf("Duplicate target names in pipe_%s", pipe_name)
        )
    }
})
