# Tests for the clustering of the batch-corrected RNA matrix
# (R/modules/rnaseq/04b_mod_clustering_batch.R)
# test-rna-batch-clustering.R
#
# All data here is synthetic; the fixtures are in helper-rna-batch.R.

bc_clustering_config <- function(enabled = TRUE) {
    cfg <- bc_config()
    cfg$modes$rna$clustering <- list(
        enabled   = enabled,
        group_col = "Condition",
        steps = list(
            hierarchical    = list(enabled = TRUE, k = 3),
            partition       = list(enabled = FALSE),
            binary_patterns = list(enabled = FALSE)
        )
    )
    cfg
}

test_that("compare_rna_cluster_assignments lines the two runs up by gene", {
    cmp <- compare_rna_cluster_assignments(c(g2 = 2, g1 = 1, g3 = 2), c(g1 = 1, g2 = 1, g4 = 3))
    expect_identical(cmp$gene, c("g1", "g2", "g3", "g4"))
    expect_identical(cmp$cluster_uncorrected, c(1, 2, 2, NA))
    expect_identical(cmp$cluster_batch_corrected, c(1, 1, NA, 3))

    expect_null(compare_rna_cluster_assignments(NULL, c(g1 = 1)))
    expect_null(compare_rna_cluster_assignments(c(g1 = 1), NULL))
    expect_null(compare_rna_cluster_assignments(c(1, 2), c(g1 = 1)))
})

test_that("no correction ran: nothing is clustered and earlier outputs are removed", {
    out_dir <- withr::local_tempdir()
    old_dir <- file.path(out_dir, "Batch_corrected_limma", "Clustering", "Hierarchical")
    dir.create(old_dir, recursive = TRUE)
    dir.create(file.path(out_dir, "Datasets"))
    stale <- c(file.path(old_dir, "Hierarchical_DE_heatmap.png"),
               file.path(out_dir, "Datasets", "clustering_batch_comparison.tsv"))
    file.create(stale)
    # The original clustering directory is not this module's to clean.
    dir.create(file.path(out_dir, "Clustering"))
    original <- file.path(out_dir, "Clustering", "Hierarchical_DE_heatmap.png")
    file.create(original)

    pre <- make_bc_pre()
    res <- suppressMessages(mod_rnaseq_clustering_batch_corr(
        pre, list(ran = FALSE), list(tables = list()), list(objects = list()),
        bc_clustering_config(), out_dir))
    expect_identical(res, list(ran = FALSE))
    expect_false(any(file.exists(stale)))
    expect_true(file.exists(original))
})

test_that("clustering disabled: the corrected run is skipped like the original", {
    out_dir <- withr::local_tempdir()
    pre <- make_bc_pre()
    cfg <- bc_clustering_config(enabled = FALSE)
    corrected <- compute_rna_batch_correction(pre, cfg)
    res <- suppressMessages(mod_rnaseq_clustering_batch_corr(
        pre, corrected, list(tables = list()), list(objects = list()), cfg, out_dir))
    expect_true(res$ran)
    expect_identical(res$method, "limma")
    expect_null(res$comparison)
    expect_length(res$files, 0)
})

test_that("the corrected matrix is what is clustered, next to an untouched original", {
    skip_if_not_installed("pheatmap")
    out_dir <- withr::local_tempdir()
    pre <- make_bc_pre()
    before <- pre
    cfg <- bc_clustering_config()
    de_res <- suppressMessages(mod_rnaseq_de(pre, list(contrasts = pre$contrasts), cfg))
    corrected <- compute_rna_batch_correction(pre, cfg)

    original <- suppressMessages(mod_rnaseq_clustering(pre, de_res, cfg, out_dir))
    original_files <- original$files[file.exists(original$files)]
    expect_gt(length(original_files), 0)
    md5_before <- tools::md5sum(original_files)

    res <- suppressMessages(mod_rnaseq_clustering_batch_corr(
        pre, corrected, de_res, original, cfg, out_dir))
    expect_identical(pre, before)
    expect_true(res$ran)

    # Its own directory, same file names; the original files are not rewritten.
    corrected_root <- file.path(out_dir, "Batch_corrected_limma", "Clustering")
    clustering_files <- setdiff(res$files,
                                file.path(out_dir, "Datasets", "clustering_batch_comparison.tsv"))
    expect_true(all(startsWith(clustering_files, corrected_root)))
    expect_setequal(basename(clustering_files), basename(original$files))
    expect_identical(tools::md5sum(original_files), md5_before)

    # The z-scores behind the corrected heatmap come from the corrected values.
    z_corrected <- res$clustering$excel_order$zscore_mat
    ids <- intersect(original$excel_order$ordered_ids, rownames(z_corrected))
    expect_gt(length(ids), 1)
    expect_equal(unname(z_corrected[ids, ]),
                 unname(zscore_rows(corrected$expr_work[ids, ])))
    expect_false(isTRUE(all.equal(unname(z_corrected[ids, ]),
                                  unname(original$excel_order$zscore_mat[ids, ]))))

    # Both assignments are kept, for the same genes.
    cmp <- res$comparison
    expect_named(cmp, c("gene", "cluster_uncorrected", "cluster_batch_corrected"))
    expect_setequal(cmp$gene, names(original$objects$clusters))
    expect_false(anyNA(cmp))
    written <- read.delim(file.path(out_dir, "Datasets", "clustering_batch_comparison.tsv"),
                          stringsAsFactors = FALSE)
    expect_identical(written$gene, cmp$gene)
})
