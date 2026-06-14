# tests/testthat/test-build-feature-annot.R
#
# Unit tests for build_feature_annot(): standardises a per-omics row_data
# feature-metadata table into the canonical feature_annot shape
# (feature IDs as rownames, ID column dropped, all other columns carried).

test_that("proteomics-style row_data keeps annotation columns, IDs become rownames", {
    rd <- data.frame(
        protein_id                = c("P1", "P2", "P3"),
        Protein.Names             = c("Alpha", "Beta", "Gamma"),
        Genes                     = c("AAA", "BBB", "CCC"),
        First.Protein.Description = c("d1", "d2", "d3"),
        stringsAsFactors = FALSE
    )

    out <- build_feature_annot(rd, "protein_id")

    expect_s3_class(out, "data.frame")
    expect_equal(rownames(out), c("P1", "P2", "P3"))
    expect_equal(colnames(out),
                 c("Protein.Names", "Genes", "First.Protein.Description"))
    expect_false("protein_id" %in% colnames(out))
})

test_that("metabolomics-style row_data drops feature_id and keeps all metadata", {
    rd <- data.frame(
        feature_id = c("m1", "m2"),
        name       = c("Glucose", "Lactate"),
        mz         = c(180.06, 90.03),
        RT         = c(1.2, 3.4),
        stringsAsFactors = FALSE
    )

    out <- build_feature_annot(rd, "feature_id")

    expect_equal(rownames(out), c("m1", "m2"))
    expect_equal(colnames(out), c("name", "mz", "RT"))
    expect_false("feature_id" %in% colnames(out))
})

test_that("rnaseq matrix-style row_data (ID only) yields a 0-column table with rownames", {
    rd <- data.frame(gene_id = c("g1", "g2", "g3"), stringsAsFactors = FALSE)

    out <- build_feature_annot(rd, "gene_id")

    expect_s3_class(out, "data.frame")
    expect_equal(ncol(out), 0L)
    expect_equal(nrow(out), 3L)
    expect_equal(rownames(out), c("g1", "g2", "g3"))
})

test_that("falls back to the first column when id_col is absent", {
    # Mirrors the rnaseq matrix/tximport branch, where the column is literally
    # "gene_id" even if config$id_columns$gene_id differs.
    rd <- data.frame(
        gene_id = c("g1", "g2"),
        symbol  = c("SYM1", "SYM2"),
        stringsAsFactors = FALSE
    )

    out <- build_feature_annot(rd, "ensembl_id")  # not present

    expect_equal(rownames(out), c("g1", "g2"))
    expect_equal(colnames(out), "symbol")
})

test_that("NULL row_data returns NULL", {
    expect_null(build_feature_annot(NULL, "feature_id"))
})

test_that("duplicate feature IDs are made unique for rownames", {
    rd <- data.frame(
        feature_id = c("m1", "m1"),
        name       = c("a", "b"),
        stringsAsFactors = FALSE
    )

    out <- build_feature_annot(rd, "feature_id")

    expect_equal(nrow(out), 2L)
    expect_equal(anyDuplicated(rownames(out)), 0L)
})
