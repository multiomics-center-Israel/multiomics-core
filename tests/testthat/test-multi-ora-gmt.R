# tests/testthat/test-multi-ora-gmt.R
#
# Unit tests for the GMT-based multi-ORA fallback used when an organism has no
# KEGG/OrgDb (non-model organisms), in R/domain/multiomics/07b_multigsea_plots.R:
#   - gmt_to_term2gene()        reshapes a GMT into TERM2GENE / TERM2NAME
#   - run_multi_ora_enricher()  clusterProfiler::enricher wrapper, KEGG-parity shape

write_test_gmt <- function() {
    path <- tempfile(fileext = ".gmt")
    lines <- c(
        paste(c("T1", "Term One", paste0("g", 1:20)), collapse = "\t"),
        paste(c("T2", "Term Two", paste0("g", 21:40)), collapse = "\t"),
        # A too-short line (no genes) must be ignored.
        paste(c("T3", "Term Three"), collapse = "\t")
    )
    writeLines(lines, path)
    path
}

test_that("gmt_to_term2gene reshapes a GMT into term/gene and term/name frames", {
    gmt <- write_test_gmt()
    on.exit(unlink(gmt), add = TRUE)

    gs <- gmt_to_term2gene(gmt)

    expect_named(gs, c("t2g", "t2n"))
    expect_setequal(colnames(gs$t2g), c("term", "gene"))
    expect_setequal(colnames(gs$t2n), c("term", "name"))
    # T3 has no genes -> only T1, T2 survive in the membership table.
    expect_setequal(unique(gs$t2g$term), c("T1", "T2"))
    expect_equal(sum(gs$t2g$term == "T1"), 20L)
    expect_equal(gs$t2n$name[gs$t2n$term == "T1"], "Term One")
})

test_that("gmt_to_term2gene returns NULL for a GMT with no usable sets", {
    empty <- tempfile(fileext = ".gmt")
    on.exit(unlink(empty), add = TRUE)
    writeLines("Tonly\tno genes here", empty)  # < 3 fields -> dropped by read_gmt
    expect_null(gmt_to_term2gene(empty))
})

test_that("run_multi_ora_enricher returns the KEGG-parity shape and finds an enriched term", {
    skip_if_not_installed("clusterProfiler")

    term2gene <- rbind(
        data.frame(term = "T1", gene = paste0("g", 1:20), stringsAsFactors = FALSE),
        data.frame(term = "T2", gene = paste0("g", 21:40), stringsAsFactors = FALSE)
    )
    term2name <- data.frame(term = c("T1", "T2"),
                            name = c("Term One", "Term Two"),
                            stringsAsFactors = FALSE)
    universe  <- paste0("g", 1:100)
    sig_genes <- paste0("g", 1:15)   # heavily loaded onto T1

    res <- run_multi_ora_enricher(sig_genes, universe, term2gene, term2name, "pooled")

    expect_s3_class(res, "data.frame")
    expect_true(all(c("pathway", "ID", "pvalue", "padj",
                      "GeneRatio", "Count", "geneID") %in% colnames(res)))
    expect_true("T1" %in% res$ID)
    # Readable label came from TERM2NAME, not the raw term ID.
    expect_equal(res$pathway[res$ID == "T1"], "Term One")
})

test_that("run_multi_ora_enricher returns NULL when too few significant genes", {
    skip_if_not_installed("clusterProfiler")
    t2g <- data.frame(term = "T1", gene = paste0("g", 1:20), stringsAsFactors = FALSE)
    expect_null(run_multi_ora_enricher(c("g1", "g2"), paste0("g", 1:100), t2g))
})
