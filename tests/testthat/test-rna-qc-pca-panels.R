# The RNA-seq report's PCA section pairs the working matrix (log2 TMM CPM by
# default) with a blind DESeq2 VST for each gene set. These tests pin which
# files the dropdown pairs and in what order, and when the VST is computed.
# All fixtures are synthetic.

touch_pngs <- function(dir, names) {
    for (n in names) file.create(file.path(dir, n))
}

test_that("list_rna_pca_panels pairs the two transforms per gene set", {
    d <- withr::local_tempdir()
    touch_pngs(d, c("PCA_PC1.vs.PC2.png", "PCA_vst_all.png",
                    "PCA_top500.png", "PCA_vst_top500.png",
                    "PCA_top2000.png", "PCA_vst_top2000.png",
                    "PCA_PC1.vs.PC3.png"))

    panels <- list_rna_pca_panels(d)

    expect_identical(panels$key, c("all", "top2000", "top500"))
    expect_identical(panels$label, c("All genes", "Top 2,000 variable genes", "Top 500 variable genes"))
    expect_identical(basename(panels$cpm_path), c("PCA_PC1.vs.PC2.png", "PCA_top2000.png", "PCA_top500.png"))
    expect_identical(basename(panels$vst_path), c("PCA_vst_all.png", "PCA_vst_top2000.png", "PCA_vst_top500.png"))
})

test_that("list_rna_pca_panels keeps a gene set when only one transform exists", {
    d <- withr::local_tempdir()
    touch_pngs(d, c("PCA_PC1.vs.PC2.png", "PCA_top1000.png"))

    panels <- list_rna_pca_panels(d)

    expect_identical(panels$key, c("all", "top1000"))
    expect_true(all(is.na(panels$vst_path)))
})

test_that("list_rna_pca_panels returns zero rows when there are no PCA images", {
    expect_equal(nrow(list_rna_pca_panels(withr::local_tempdir())), 0)
})

synthetic_pre <- function(n_genes = 300, seed = 1) {
    withr::with_seed(seed, {
        counts <- matrix(rnbinom(n_genes * 6, mu = 200, size = 5), nrow = n_genes,
                         dimnames = list(paste0("g", seq_len(n_genes)), paste0("S", 1:6)))
    })
    meta <- data.frame(SampleID = paste0("S", 1:6), group = rep(c("A", "B"), each = 3))
    work <- log2(counts + 1)
    attr(work, "method") <- "TMMlogCPM"
    list(expr_filt = counts, expr_work = work,
         de_input = annotate_source_type(counts, "matrix"),
         meta = meta, info = list(source_type = "matrix"))
}

test_that("compute_rna_qc_vst returns a VST matrix on the filtered genes", {
    skip_if_not_installed("DESeq2")
    pre <- synthetic_pre()

    vst <- suppressMessages(compute_rna_qc_vst(pre, sample_col = "SampleID"))

    expect_true(is.matrix(vst))
    expect_identical(dim(vst), dim(pre$expr_filt))
    expect_identical(attr(vst, "method"), "VST")
})

test_that("compute_rna_qc_vst skips preprocessed input and an expr_work that is already VST", {
    pre <- synthetic_pre()

    pre_pp <- pre
    pre_pp$info$source_type <- "preprocessed"
    expect_null(compute_rna_qc_vst(pre_pp))

    pre_vst <- pre
    attr(pre_vst$expr_work, "method") <- "VST"
    expect_null(compute_rna_qc_vst(pre_vst))
})

test_that("write_rnaseq_datasets_legacy removes a normalized matrix left by another method", {
    d <- withr::local_tempdir()
    dirs <- create_legacy_output_dirs(d)
    file.create(file.path(dirs$datasets, "rna_norm_TMMlogCPM.tsv"))
    pre <- synthetic_pre(n_genes = 20)
    attr(pre$expr_work, "method") <- "VST"

    write_rnaseq_datasets_legacy(pre, config = list(), dirs = dirs)

    expect_identical(list.files(dirs$datasets, pattern = "^rna_norm_"), "rna_norm_VST.tsv")
    expect_true(file.exists(file.path(dirs$datasets, "rna_counts_filtered.tsv")))
})

test_that("write_rna_pca_transform_labels records the work and VST transforms", {
    d <- withr::local_tempdir()
    cfg <- list(normalization = list(method = "TMMlogCPM"))
    work <- matrix(0, 2, 2)
    attr(work, "method") <- "TMMlogCPM"

    write_rna_pca_transform_labels(d, rna_qc_transform_label(work, cfg),
                                   vst_label = "DESeq2 VST (blind)")

    labels <- read_rna_pca_transform_labels(d)
    expect_identical(labels[["work"]], "log2(TMM CPM + 1)")
    expect_identical(labels[["vst"]], "DESeq2 VST (blind)")
})

test_that("a rerun on VST rewrites the labels and drops the earlier VST row", {
    d <- withr::local_tempdir()
    write_rna_pca_transform_labels(d, "log2(TMM CPM + 1)", vst_label = "DESeq2 VST (blind)")
    work <- matrix(0, 2, 2)
    attr(work, "method") <- "VST"

    write_rna_pca_transform_labels(d, rna_qc_transform_label(work, list()))

    tbl <- utils::read.delim(file.path(d, "pca_transform_labels.tsv"), stringsAsFactors = FALSE)
    expect_identical(tbl$panel, "work")
    expect_identical(read_rna_pca_transform_labels(d)[["work"]], "DESeq2 VST (blind)")
})

test_that("read_rna_pca_transform_labels falls back when the file is missing", {
    labels <- read_rna_pca_transform_labels(withr::local_tempdir())
    expect_identical(labels[["work"]], "Normalized expression")
    expect_identical(labels[["vst"]], "DESeq2 VST (blind)")
})

test_that("heatmap sample IDs fall back when effects$samples is unset", {
    meta <- data.frame(SampleID = c("S1", "S2"), Lib = c("L1", "L2"),
                       Group = c("A", "B"), row.names = c("r1", "r2"))
    expect_identical(rna_heatmap_sample_ids(meta, list()), c("S1", "S2"))
    expect_identical(rna_heatmap_sample_ids(meta, list(effects = list(samples = "Lib"))),
                     c("L1", "L2"))
    expect_identical(rna_heatmap_sample_ids(meta, list(), de_sample_col = "Lib"), c("L1", "L2"))
    # No such column: rownames, as the DE runners use.
    expect_identical(rna_heatmap_sample_ids(meta[, c("Lib", "Group")], list()), c("r1", "r2"))
})
