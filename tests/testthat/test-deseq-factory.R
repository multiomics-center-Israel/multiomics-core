# tests/testthat/test-deseq-factory.R
# Integration tests for DESeq2 factory with dual input support (matrix and tximport)

# =============================================================================
# Helper: Create mock tximport object
# =============================================================================
create_mock_tximport <- function(n_genes = 100, n_samples = 6,
                                  counts_from_abundance = "no") {
    set.seed(42)

    gene_ids <- paste0("ENSG", sprintf("%05d", 1:n_genes))
    sample_ids <- paste0("Sample", 1:n_samples)

    # Create realistic-ish counts
    counts <- matrix(
        rpois(n_genes * n_samples, lambda = 100),
        nrow = n_genes, ncol = n_samples,
        dimnames = list(gene_ids, sample_ids)
    )

    # Create abundance (TPM-like values)
    abundance <- matrix(
        runif(n_genes * n_samples, min = 0.1, max = 1000),
        nrow = n_genes, ncol = n_samples,
        dimnames = list(gene_ids, sample_ids)
    )
    # Normalize columns to sum to 1e6 (TPM)
    abundance <- sweep(abundance, 2, colSums(abundance), "/") * 1e6

    # Create length matrix (effective lengths)
    gene_lengths <- runif(n_genes, min = 500, max = 5000)
    length_mat <- matrix(
        rep(gene_lengths, n_samples),
        nrow = n_genes, ncol = n_samples,
        dimnames = list(gene_ids, sample_ids)
    )

    list(
        counts = counts,
        abundance = abundance,
        length = length_mat,
        countsFromAbundance = counts_from_abundance
    )
}

# =============================================================================
# Helper: Create mock count matrix
# =============================================================================
create_mock_counts <- function(n_genes = 100, n_samples = 6) {
    set.seed(42)

    gene_ids <- paste0("ENSG", sprintf("%05d", 1:n_genes))
    sample_ids <- paste0("Sample", 1:n_samples)

    counts <- matrix(
        rpois(n_genes * n_samples, lambda = 100),
        nrow = n_genes, ncol = n_samples,
        dimnames = list(gene_ids, sample_ids)
    )

    storage.mode(counts) <- "integer"
    counts
}

# =============================================================================
# Helper: Create mock metadata
# =============================================================================
create_mock_meta <- function(sample_ids, groups = NULL) {
    n_samples <- length(sample_ids)
    if (is.null(groups)) {
        groups <- rep(c("Control", "Treatment"), length.out = n_samples)
    }

    data.frame(
        SampleID = sample_ids,
        Condition = factor(groups),
        row.names = sample_ids,
        stringsAsFactors = FALSE
    )
}

# =============================================================================
# Tests: Source Type Detection
# =============================================================================
test_that("detect_source_type correctly identifies tximport objects", {
    txi <- create_mock_tximport()
    expect_equal(detect_source_type(txi), "tximport")
})

test_that("detect_source_type correctly identifies matrix input", {
    counts <- create_mock_counts()
    expect_equal(detect_source_type(counts), "matrix")
})

test_that("detect_source_type respects and validates annotation", {
    txi <- create_mock_tximport()
    attr(txi, "source_type") <- "tximport"

    # Should work when annotation matches structure
    expect_equal(detect_source_type(txi), "tximport")

    # Should error when annotation doesn't match structure
    counts <- create_mock_counts()
    attr(counts, "source_type") <- "tximport"  # Wrong annotation
    expect_error(detect_source_type(counts), "does not match object structure")
})

# =============================================================================
# Tests: tximport Validation
# =============================================================================
test_that("is_valid_tximport_structure validates correctly", {
    txi <- create_mock_tximport()

    # Valid tximport

    expect_true(is_valid_tximport_structure(txi, validate_only = TRUE))

    # Missing field
    txi_bad <- txi
    txi_bad$length <- NULL
    expect_false(is_valid_tximport_structure(txi_bad, validate_only = TRUE))

    # Dimension mismatch
    txi_bad <- txi
    txi_bad$length <- txi_bad$length[1:50, , drop = FALSE]
    expect_false(is_valid_tximport_structure(txi_bad, validate_only = TRUE))
})

test_that("validate_tximport stops on invalid structure", {
    txi_bad <- list(counts = matrix(1:10, 5, 2))  # Missing abundance and length
    expect_error(validate_tximport(txi_bad), "Missing required fields")
})

# =============================================================================
# Tests: Count Matrix Validation
# =============================================================================
test_that("validate_count_matrix accepts valid integer counts", {
    counts <- create_mock_counts()
    expect_silent(validate_count_matrix(counts))
})
test_that("validate_count_matrix rejects non-integer values", {
    counts <- create_mock_counts()
    counts[1, 1] <- 10.5  # Non-integer value

    expect_error(
        validate_count_matrix(counts),
        "Non-integer counts detected"
    )
})

test_that("validate_count_matrix rejects negative values", {
    counts <- create_mock_counts()
    counts[1, 1] <- -5

    expect_error(validate_count_matrix(counts), "negative values")
})

# =============================================================================
# Tests: Sample Alignment
# =============================================================================
test_that("align_samples_strict fails on mismatch by default", {
    samples <- paste0("Sample", 1:6)
    meta <- create_mock_meta(paste0("Sample", 1:5))  # Missing Sample6

    expect_error(
        align_samples_strict(samples, meta, "SampleID", lenient = FALSE),
        "Sample mismatch detected"
    )
})

test_that("align_samples_strict allows intersection in lenient mode", {
    samples <- paste0("Sample", 1:6)
    meta <- create_mock_meta(paste0("Sample", 1:5))

    expect_warning(
        result <- align_samples_strict(samples, meta, "SampleID", lenient = TRUE),
        "Sample mismatch"
    )

    expect_equal(nrow(result), 5)
    expect_equal(rownames(result), paste0("Sample", 1:5))
})

# =============================================================================
# Tests: DESeqDataSet Factory
# =============================================================================
test_that("create_deseq_dataset works with raw count matrix", {
    skip_if_not_installed("DESeq2")

    counts <- create_mock_counts()
    meta <- create_mock_meta(colnames(counts))

    dds <- create_deseq_dataset(
        expr = counts,
        meta = meta,
        design = ~ Condition,
        sample_col = "SampleID"
    )

    expect_s4_class(dds, "DESeqDataSet")
    expect_equal(ncol(dds), ncol(counts))
    expect_equal(nrow(dds), nrow(counts))
    expect_equal(S4Vectors::metadata(dds)$source_type, "matrix")
})

test_that("create_deseq_dataset works with tximport object", {
    skip_if_not_installed("DESeq2")

    txi <- create_mock_tximport()
    meta <- create_mock_meta(colnames(txi$counts))

    dds <- create_deseq_dataset(
        expr = txi,
        meta = meta,
        design = ~ Condition,
        sample_col = "SampleID"
    )

    expect_s4_class(dds, "DESeqDataSet")
    expect_equal(ncol(dds), ncol(txi$counts))
    expect_equal(nrow(dds), nrow(txi$counts))
    expect_equal(S4Vectors::metadata(dds)$source_type, "tximport")
})

# =============================================================================
# Tests: tximport Subsetting
# =============================================================================
test_that("subset_tximport subsets all matrices together", {
    txi <- create_mock_tximport(n_samples = 6)
    samples_to_keep <- paste0("Sample", 1:3)

    txi_sub <- subset_tximport(txi, samples = samples_to_keep)

    expect_equal(ncol(txi_sub$counts), 3)
    expect_equal(ncol(txi_sub$abundance), 3)
    expect_equal(ncol(txi_sub$length), 3)
    expect_equal(colnames(txi_sub$counts), samples_to_keep)
})

test_that("subset_tximport_genes subsets all matrices together", {
    txi <- create_mock_tximport(n_genes = 100)

    # Using logical keep_vec
    keep_vec <- rep(c(TRUE, FALSE), 50)
    txi_sub <- subset_tximport_genes(txi, genes = keep_vec)

    expect_equal(nrow(txi_sub$counts), 50)
    expect_equal(nrow(txi_sub$abundance), 50)
    expect_equal(nrow(txi_sub$length), 50)
})

# =============================================================================
# Tests: Normalization with Dual Input
# =============================================================================
test_that("normalize_counts works with matrix input (TMMlogCPM)", {
    skip_if_not_installed("edgeR")

    counts <- create_mock_counts()
    meta <- create_mock_meta(colnames(counts))

    result <- normalize_counts(counts, meta, method = "TMMlogCPM")

    expect_true(is.matrix(result))
    expect_equal(dim(result), dim(counts))
    expect_equal(attr(result, "method"), "TMMlogCPM")
    expect_equal(attr(result, "source_type"), "matrix")
})

test_that("normalize_counts works with matrix input (VST)", {
    skip_if_not_installed("DESeq2")

    counts <- create_mock_counts()
    meta <- create_mock_meta(colnames(counts))

    result <- normalize_counts(counts, meta, method = "VST", sample_col = "SampleID")

    expect_true(is.matrix(result))
    expect_equal(attr(result, "method"), "VST")
    expect_equal(attr(result, "source_type"), "matrix")
})

test_that("normalize_counts works with tximport input (VST)", {
    skip_if_not_installed("DESeq2")

    txi <- create_mock_tximport()
    meta <- create_mock_meta(colnames(txi$counts))

    result <- normalize_counts(txi, meta, method = "VST", sample_col = "SampleID")

    expect_true(is.matrix(result))
    expect_equal(attr(result, "method"), "VST")
    expect_equal(attr(result, "source_type"), "tximport")
})

test_that("normalize_counts rejects tximport with TMMlogCPM", {
    txi <- create_mock_tximport()
    meta <- create_mock_meta(colnames(txi$counts))

    expect_error(
        normalize_counts(txi, meta, method = "TMMlogCPM"),
        "TMMlogCPM requires raw count matrix"
    )
})

# =============================================================================
# Tests: Integration - run_deseq2_de
# =============================================================================
test_that("run_deseq2_de works with raw count matrix", {
    skip_if_not_installed("DESeq2")

    counts <- create_mock_counts()
    meta <- create_mock_meta(colnames(counts))

    contrasts_df <- data.frame(
        Contrast_name = "Treatment_vs_Control",
        Factor = "Condition",
        Numerator = "Treatment",
        Denominator = "Control",
        stringsAsFactors = FALSE
    )

    de_cfg <- list(
        sample_col = "SampleID",
        p_cutoff = 0.05
    )

    result <- run_deseq2_de(counts, meta, contrasts_df, de_cfg)

    expect_type(result, "list")
    expect_s4_class(result$dds, "DESeqDataSet")
    expect_true("Treatment_vs_Control" %in% names(result$tables))
    expect_true("log2FoldChange" %in% colnames(result$tables$Treatment_vs_Control))
})

test_that("run_deseq2_de works with tximport object", {
    skip_if_not_installed("DESeq2")

    txi <- create_mock_tximport()
    meta <- create_mock_meta(colnames(txi$counts))

    contrasts_df <- data.frame(
        Contrast_name = "Treatment_vs_Control",
        Factor = "Condition",
        Numerator = "Treatment",
        Denominator = "Control",
        stringsAsFactors = FALSE
    )

    de_cfg <- list(
        sample_col = "SampleID",
        p_cutoff = 0.05
    )

    result <- run_deseq2_de(txi, meta, contrasts_df, de_cfg)

    expect_type(result, "list")
    expect_s4_class(result$dds, "DESeqDataSet")
    expect_equal(S4Vectors::metadata(result$dds)$source_type, "tximport")
    expect_true("Treatment_vs_Control" %in% names(result$tables))
})

# =============================================================================
# Tests: Anti-patterns (should fail)
# =============================================================================
test_that("non-integer counts in matrix path fails", {
    skip_if_not_installed("DESeq2")

    counts <- create_mock_counts()
    counts[1, 1] <- 10.5  # Make it non-integer
    meta <- create_mock_meta(colnames(counts))

    expect_error(
        create_deseq_dataset(counts, meta, ~ Condition, "SampleID"),
        "Non-integer counts"
    )
})

test_that("misannotated source_type fails validation", {
    txi <- create_mock_tximport()
    # Incorrectly annotate as matrix
    attr(txi, "source_type") <- "matrix"

    expect_error(
        detect_source_type(txi),
        "does not match object structure"
    )
})
