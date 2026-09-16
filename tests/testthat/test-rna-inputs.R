# tests/testthat/test-rna-inputs.R
# Tests for RNA input helpers (R/domain/rnaseq/00_inputs.R)

# Write a minimal RSEM gene-level results file. tximport's "rsem" importer
# reads the gene_id, effective_length, expected_count and TPM columns.
write_mock_rsem <- function(path, eff_length) {
    df <- data.frame(
        gene_id = c("ENSG00001", "ENSG00002", "ENSG00003"),
        `transcript_id(s)` = c("t1", "t2", "t3"),
        length = c(1000, 2000, 1500),
        effective_length = eff_length,
        expected_count = c(100, 0, 50),
        TPM = c(10.5, 0, 5.2),
        FPKM = c(8.1, 0, 4.0),
        check.names = FALSE
    )
    write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

test_that("load_rsem_as_tximport imports gene-level RSEM output", {
    skip_if_not_installed("tximport")

    dir <- withr::local_tempdir()
    write_mock_rsem(file.path(dir, "SampleA.genes.results"), c(800, 0, 1200))
    write_mock_rsem(file.path(dir, "SampleB.genes.results"), c(820, 0, 1180))

    txi <- load_rsem_as_tximport(dir)

    expect_true(all(c("counts", "abundance", "length") %in% names(txi)))
    # Sample IDs derived from filenames, in sorted order
    expect_equal(colnames(txi$counts), c("SampleA", "SampleB"))
    expect_equal(rownames(txi$counts), c("ENSG00001", "ENSG00002", "ENSG00003"))
    # Result must satisfy the structure the pipeline validates
    expect_true(is_valid_tximport_structure(txi, validate_only = TRUE))
})

test_that("load_rsem_as_tximport replaces zero effective lengths with 1", {
    skip_if_not_installed("tximport")

    dir <- withr::local_tempdir()
    # Gene 2 has effective_length 0 in both samples
    write_mock_rsem(file.path(dir, "SampleA.genes.results"), c(800, 0, 1200))
    write_mock_rsem(file.path(dir, "SampleB.genes.results"), c(820, 0, 1180))

    txi <- load_rsem_as_tximport(dir, fix_zero_length = TRUE)
    expect_true(all(txi$length > 0))

    txi_raw <- load_rsem_as_tximport(dir, fix_zero_length = FALSE)
    expect_true(any(txi_raw$length == 0))
})

test_that("load_rsem_as_tximport honours explicit sample_names", {
    skip_if_not_installed("tximport")

    dir <- withr::local_tempdir()
    write_mock_rsem(file.path(dir, "SampleA.genes.results"), c(800, 0, 1200))
    write_mock_rsem(file.path(dir, "SampleB.genes.results"), c(820, 0, 1180))

    txi <- load_rsem_as_tximport(dir, sample_names = c("ctrl_1", "trt_1"))
    expect_equal(colnames(txi$counts), c("ctrl_1", "trt_1"))
})

test_that("load_rsem_as_tximport errors on empty dir and bad inputs", {
    dir <- withr::local_tempdir()

    # No matching files
    expect_error(load_rsem_as_tximport(dir), "No files matching")

    # Nonexistent directory
    expect_error(
        load_rsem_as_tximport(file.path(dir, "does-not-exist")),
        "Directory does not exist"
    )

    skip_if_not_installed("tximport")
    write_mock_rsem(file.path(dir, "SampleA.genes.results"), c(800, 0, 1200))
    # sample_names length mismatch
    expect_error(
        load_rsem_as_tximport(dir, sample_names = c("a", "b")),
        "files were matched"
    )
})
