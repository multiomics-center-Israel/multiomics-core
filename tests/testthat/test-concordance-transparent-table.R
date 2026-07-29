# tests/testthat/test-concordance-transparent-table.R
#
# Covers the transparent RNA-protein log2FC concordance contract:
#   - build_gene_protein_mapping_from_ids() rebuilds the mapping from scratch
#     (from a custom file + the original ID space), independent of any MAE.
#   - analyze_rna_protein_concordance() returns a checkable table carrying the
#     original gene_id / protein_id beside both log2FCs, and the reported
#     correlation is computed from that exact table.

make_mapping_file <- function(pairs) {
    f <- tempfile(fileext = ".csv")
    utils::write.csv(pairs, f, row.names = FALSE)
    f
}

test_that("build_gene_protein_mapping_from_ids reads the custom file and filters to present IDs", {
    pairs <- data.frame(
        gene_id = paste0("EHI_", sprintf("%03d", 1:5)),
        protein_id = paste0("XP_", 1:5),
        stringsAsFactors = FALSE
    )
    map_file <- make_mapping_file(pairs)
    on.exit(unlink(map_file), add = TRUE)

    config <- list(modes = list(multiomics = list(gene_protein_mapping_file = map_file)))

    res <- build_gene_protein_mapping_from_ids(
        gene_ids = c("EHI_001", "EHI_002", "EHI_999"),   # EHI_999 absent from file
        protein_ids = c("XP_1", "XP_2", "XP_999"),
        config = config
    )

    expect_true(all(c("gene_id", "protein_id", "mapping_source") %in% names(res)))
    expect_identical(unique(res$mapping_source), "custom_file")
    # Only the pairs whose gene AND protein are in the supplied ID space survive.
    expect_setequal(res$gene_id, c("EHI_001", "EHI_002"))
    expect_setequal(res$protein_id, c("XP_1", "XP_2"))
})

test_that("build_gene_protein_mapping_from_ids returns NULL when no mapping file is configured", {
    expect_null(build_gene_protein_mapping_from_ids("EHI_001", "XP_1", config = list()))
})

test_that("concordance table carries original IDs beside both log2FCs, and drives the correlation", {
    n <- 10
    genes <- paste0("EHI_", sprintf("%03d", seq_len(n)))
    prots <- paste0("XP_", seq_len(n))
    mapping <- data.frame(gene_id = genes, protein_id = prots,
                          mapping_source = "custom_file", stringsAsFactors = FALSE)

    # Deterministic log2FCs (no RNG): partially concordant in sign/magnitude.
    rna_lfc  <- c(-1.5, 1.2, -0.8, 0.3, 2.1, -2.0, 0.9, -1.1, 0.2, 1.7)
    prot_lfc <- c(-2.0, 0.7, -1.4, -0.1, 1.3, -0.5, 1.6, -1.9, 1.0, 0.4)

    rna_de <- data.frame(feature_id = genes, logFC = rna_lfc,
                         padj = rep(0.2, n), stringsAsFactors = FALSE)
    prot_de <- data.frame(feature_id = prots, logFC = prot_lfc,
                          padj = rep(0.2, n), stringsAsFactors = FALSE)

    res <- analyze_rna_protein_concordance(rna_de, prot_de, mapping,
                                           config = list(), out_dir = NULL)
    tbl <- res$concordance_table

    # Transparent: real IDs on both sides + both log2FCs.
    expect_true(all(c("gene_id", "protein_id", "logFC_rna", "logFC_prot") %in% names(tbl)))
    expect_equal(nrow(tbl), n)
    expect_setequal(tbl$gene_id, genes)
    expect_setequal(tbl$protein_id, prots)

    # The reported correlation must be the correlation OF THIS TABLE.
    expect_equal(
        unname(res$stats$correlation),
        unname(stats::cor(tbl$logFC_rna, tbl$logFC_prot)),
        tolerance = 1e-8
    )
})
