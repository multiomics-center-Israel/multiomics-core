# tests/testthat/test-multiomics-harmonization.R
# Tests for mod_multiomics_harmonization() with mock data.
# Validates graceful handling of missing, partial, and correct data configurations.

skip_if_not_installed("MultiAssayExperiment")
skip_if_not_installed("SummarizedExperiment")

# =============================================================================
# Helper: Create mock preprocessed omics data (matches *_pre object structure)
# =============================================================================
create_mock_pre <- function(omics_type = "rna", n_features = 50, n_samples = 6,
                            sample_prefix = "Sample") {
    set.seed(42)

    sample_ids <- paste0(sample_prefix, seq_len(n_samples))
    feature_ids <- paste0(toupper(substr(omics_type, 1, 1)), "FEAT",
                          sprintf("%04d", seq_len(n_features)))

    # Expression matrix (features x samples)
    expr_work <- matrix(
        rnorm(n_features * n_samples, mean = 10, sd = 2),
        nrow = n_features, ncol = n_samples,
        dimnames = list(feature_ids, sample_ids)
    )

    expr_raw <- matrix(
        as.integer(rpois(n_features * n_samples, lambda = 100)),
        nrow = n_features, ncol = n_samples,
        dimnames = list(feature_ids, sample_ids)
    )

    # Metadata
    meta <- data.frame(
        SampleID = sample_ids,
        condition = factor(rep(c("Control", "Treatment"), length.out = n_samples)),
        row.names = sample_ids,
        stringsAsFactors = FALSE
    )

    # Row data (feature annotations)
    row_data <- data.frame(
        feature_id = feature_ids,
        gene_symbol = paste0("GENE", seq_len(n_features)),
        row.names = feature_ids,
        stringsAsFactors = FALSE
    )

    list(
        expr_raw  = expr_raw,
        expr_filt = expr_work,
        expr_work = expr_work,
        meta      = meta,
        row_data  = row_data
    )
}

# =============================================================================
# Helper: Create mock config object with multiomics structure
# =============================================================================
create_mock_config <- function() {
    list(
        project = list(
            name = "TestProject",
            analysis_round = "01",
            dir = tempdir(),
            analyst = "Test"
        ),
        paths = list(
            raw = "data",
            out = "outputs"
        ),
        params = list(
            seed = 42
        ),
        design = list(
            condition_column = "condition"
        ),
        modes = list(
            rna = list(
                effects = list(color = "condition", samples = "SampleID"),
                id_columns = list(gene_id_col = "gene_id", sample_col = "SampleID")
            ),
            proteomics = list(
                effects = list(color = "condition", samples = "SampleID"),
                id_columns = list(protein_id = "Protein.Group", sample_col = "SampleID")
            ),
            metabolomics = list(
                effects = list(color = "condition", samples = "SampleID"),
                id_columns = list(sample_col = "SampleID")
            ),
            multiomics = list(
                condition_column = "condition",
                require_one_to_one_mapping = FALSE,
                feature_selection = list(method = "variance", top_n = 500),
                integration = list(methods = c("DIABLO", "SNF"))
            )
        )
    )
}


# =============================================================================
# Tests: All NULL inputs should fail
# =============================================================================

test_that("harmonization fails when all inputs are NULL", {
    config <- create_mock_config()
    out_dir <- file.path(tempdir(), "test_harm_null")

    expect_error(
        mod_multiomics_harmonization(
            config = config,
            rna_pre = NULL,
            prot_pre = NULL,
            metab_pre = NULL,
            out_dir = out_dir
        ),
        "at least 2 omics"
    )
})


# =============================================================================
# Tests: Single omics should fail
# =============================================================================

test_that("harmonization fails with only 1 omics layer (rna only)", {
    rna <- create_mock_pre("rna")
    config <- create_mock_config()
    out_dir <- file.path(tempdir(), "test_harm_1omics_rna")

    expect_error(
        mod_multiomics_harmonization(
            config = config,
            rna_pre = rna,
            prot_pre = NULL,
            metab_pre = NULL,
            out_dir = out_dir
        ),
        "at least 2 omics"
    )
})

test_that("harmonization fails with only 1 omics layer (proteomics only)", {
    prot <- create_mock_pre("prot")
    config <- create_mock_config()
    out_dir <- file.path(tempdir(), "test_harm_1omics_prot")

    expect_error(
        mod_multiomics_harmonization(
            config = config,
            rna_pre = NULL,
            prot_pre = prot,
            metab_pre = NULL,
            out_dir = out_dir
        ),
        "at least 2 omics"
    )
})


# =============================================================================
# Tests: Non-overlapping samples should fail
# =============================================================================

test_that("harmonization fails with no sample overlap", {
    rna  <- create_mock_pre("rna", n_features = 50, n_samples = 6,
                            sample_prefix = "RNA_S")
    prot <- create_mock_pre("prot", n_features = 30, n_samples = 6,
                            sample_prefix = "PROT_S")
    config <- create_mock_config()
    out_dir <- file.path(tempdir(), "test_harm_no_overlap")

    expect_error(
        mod_multiomics_harmonization(
            config = config,
            rna_pre = rna,
            prot_pre = prot,
            metab_pre = NULL,
            out_dir = out_dir
        ),
        "Insufficient sample overlap|Common samples: 0"
    )
})


# =============================================================================
# Tests: Two omics with matching samples should succeed
# =============================================================================

test_that("harmonization succeeds with 2 omics and overlapping samples", {
    rna  <- create_mock_pre("rna", n_features = 50, n_samples = 6)
    prot <- create_mock_pre("prot", n_features = 30, n_samples = 6)
    config <- create_mock_config()
    out_dir <- file.path(tempdir(), paste0("test_harm_2omics_", Sys.getpid()))

    result <- mod_multiomics_harmonization(
        config = config,
        rna_pre = rna,
        prot_pre = prot,
        metab_pre = NULL,
        out_dir = out_dir
    )

    # Validate return structure
    expect_type(result, "list")
    expect_true("mae" %in% names(result))
    expect_true("gene_protein_mapping" %in% names(result))
    expect_true("inputs" %in% names(result))
    expect_true("omics_available" %in% names(result))

    # MAE should be a MultiAssayExperiment
    expect_s4_class(result$mae, "MultiAssayExperiment")

    # Should have 2 omics
    expect_equal(length(result$omics_available), 2)
    expect_true("transcriptomics" %in% result$omics_available)
    expect_true("proteomics" %in% result$omics_available)
})


# =============================================================================
# Tests: Three omics should succeed
# =============================================================================

test_that("harmonization succeeds with 3 omics layers", {
    rna   <- create_mock_pre("rna", n_features = 50, n_samples = 6)
    prot  <- create_mock_pre("prot", n_features = 30, n_samples = 6)
    metab <- create_mock_pre("metab", n_features = 20, n_samples = 6)
    config <- create_mock_config()
    out_dir <- file.path(tempdir(), paste0("test_harm_3omics_", Sys.getpid()))

    result <- mod_multiomics_harmonization(
        config = config,
        rna_pre = rna,
        prot_pre = prot,
        metab_pre = metab,
        out_dir = out_dir
    )

    expect_equal(length(result$omics_available), 3)
    expect_s4_class(result$mae, "MultiAssayExperiment")
    expect_true("metabolomics" %in% result$omics_available)
})


# =============================================================================
# Tests: Output directory is created
# =============================================================================

test_that("harmonization creates output directory", {
    rna  <- create_mock_pre("rna", n_features = 50, n_samples = 6)
    prot <- create_mock_pre("prot", n_features = 30, n_samples = 6)
    config <- create_mock_config()
    out_dir <- file.path(tempdir(), paste0("test_harm_dir_", Sys.getpid(), "_", sample(1e6, 1)))

    result <- mod_multiomics_harmonization(
        config = config,
        rna_pre = rna,
        prot_pre = prot,
        metab_pre = NULL,
        out_dir = out_dir
    )

    expect_true(dir.exists(out_dir))
})


# =============================================================================
# Tests: load_multiomics_inputs validation
# =============================================================================

test_that("load_multiomics_inputs fails with 0 omics", {
    config <- create_mock_config()
    expect_error(
        load_multiomics_inputs(config, rna_pre = NULL, prot_pre = NULL, metab_pre = NULL),
        "at least 2 omics"
    )
})

test_that("load_multiomics_inputs returns correct structure with 2 omics", {
    rna  <- create_mock_pre("rna")
    prot <- create_mock_pre("prot")
    config <- create_mock_config()

    result <- load_multiomics_inputs(config, rna_pre = rna, prot_pre = prot, metab_pre = NULL)

    expect_type(result, "list")
    expect_true("omics_available" %in% names(result))
    expect_true("n_omics" %in% names(result))
    expect_equal(result$n_omics, 2)
    expect_true("transcriptomics" %in% result$omics_available)
    expect_true("proteomics" %in% result$omics_available)
})


# =============================================================================
# Tests: validate_multiomics_inputs
# =============================================================================

test_that("validate_multiomics_inputs passes with sufficient overlap", {
    rna  <- create_mock_pre("rna", n_samples = 6)
    prot <- create_mock_pre("prot", n_samples = 6)
    config <- create_mock_config()

    inputs <- load_multiomics_inputs(config, rna_pre = rna, prot_pre = prot)

    # validate_multiomics_inputs prints informational messages on success
    expect_no_error(
        validate_multiomics_inputs(inputs, config)
    )
})

test_that("validate_multiomics_inputs fails when sample column missing from metadata", {
    rna  <- create_mock_pre("rna", n_samples = 6)
    prot <- create_mock_pre("prot", n_samples = 6)

    # Remove the SampleID column from RNA metadata
    rna$meta$SampleID <- NULL

    config <- create_mock_config()
    inputs <- load_multiomics_inputs(config, rna_pre = rna, prot_pre = prot)

    expect_error(
        validate_multiomics_inputs(inputs, config),
        "Sample column.*not found"
    )
})
