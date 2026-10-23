# Tests for de$block_col in the RNA-seq DESeq2 path and describe_rna_de()
# test-rnaseq-de-block.R

# Six libraries, three blocks, both groups in every block. Synthetic only.
make_rna_block_meta <- function() {
    data.frame(
        SampleID = c("B1_ctl", "B1_trt", "B2_ctl", "B2_trt", "B3_ctl", "B3_trt"),
        Group    = rep(c("ctl", "trt"), times = 3),
        Block    = rep(c("B1", "B2", "B3"), each = 2),
        stringsAsFactors = FALSE
    )
}

# Large gene-specific block effects plus a modest, consistent treatment effect
# on the first 40 genes: the case pairing is meant for. The block effects vary
# by gene so that size factors cannot absorb them the way they absorb depth.
make_rna_block_counts <- function(meta, n_genes = 400, seed = 7) {
    withr::with_seed(seed, {
        base <- rgamma(n_genes, shape = 2, rate = 0.01)
        block_eff <- matrix(exp(rnorm(n_genes * 3, sd = 1)), nrow = n_genes,
                            dimnames = list(NULL, c("B1", "B2", "B3")))
        mu <- base * block_eff[, meta$Block]
        mu[1:40, meta$Group == "trt"] <- mu[1:40, meta$Group == "trt"] * 1.8
        counts <- matrix(rnbinom(length(mu), mu = mu, size = 20), nrow = n_genes)
    })
    rownames(counts) <- paste0("GENE", seq_len(n_genes))
    colnames(counts) <- meta$SampleID
    counts
}

rna_block_contrasts <- function() {
    data.frame(Contrast_name = "trt_vs_ctl", Factor = "Group",
               Numerator = "trt", Denominator = "ctl", stringsAsFactors = FALSE)
}

# ---- describe_rna_de ---------------------------------------------------------

test_that("describe_rna_de names the unpaired DESeq2 design that ran", {
    d <- describe_rna_de(list(method = "deseq2", deseq_mode = "default",
                              design = "~Group", block_col = NULL))
    expect_identical(d$label, "DESeq2")
    expect_match(d$short, "design ~Group", fixed = TRUE)
    expect_match(d$methods, "Wald test, default mode", fixed = TRUE)
    expect_false(grepl("fixed term", d$methods))
})

test_that("describe_rna_de states the block as a fixed term with its cost", {
    d <- describe_rna_de(list(method = "deseq2", deseq_mode = "default",
                              design = "~Block + Group", block_col = "Block",
                              n_blocks = 3L))
    expect_match(d$short, "~Block + Group", fixed = TRUE)
    expect_match(d$methods, "Block as a fixed term", fixed = TRUE)
    expect_match(d$methods, "each of its 3 levels", fixed = TRUE)
    expect_match(d$methods, "no random effect", fixed = TRUE)
})

test_that("describe_rna_de covers limma, pre-computed and missing info", {
    expect_match(describe_rna_de(list(method = "limma", design = "~0 + Group"))$short,
                 "limma", fixed = TRUE)
    expect_match(describe_rna_de(list(method = "precomputed"))$methods,
                 "did not fit a model", fixed = TRUE)
    none <- describe_rna_de(NULL)
    expect_identical(none$method, "unknown")
    expect_match(none$methods, "not recorded", fixed = TRUE)
})

# ---- run_deseq2_de with de$block_col ----------------------------------------

test_that("run_deseq2_de keeps ~ Group when no block_col is set", {
    skip_if_not_installed("DESeq2")
    meta <- make_rna_block_meta()
    res <- suppressMessages(run_deseq2_de(make_rna_block_counts(meta), meta,
                                          rna_block_contrasts(),
                                          list(sample_col = "SampleID")))
    expect_identical(res$info$method, "deseq2")
    expect_identical(res$info$design, "~Group")
    expect_null(res$info$block_col)
    expect_identical(as.character(DESeq2::design(res$dds))[2], "Group")
})

test_that("run_deseq2_de fits ~ block + Group and keeps the contrast table shape", {
    skip_if_not_installed("DESeq2")
    meta <- make_rna_block_meta()
    counts <- make_rna_block_counts(meta)
    cfg <- list(sample_col = "SampleID", block_col = "Block")

    paired <- suppressMessages(run_deseq2_de(counts, meta, rna_block_contrasts(), cfg))
    unpaired <- suppressMessages(run_deseq2_de(counts, meta, rna_block_contrasts(),
                                               list(sample_col = "SampleID")))

    expect_identical(paired$info$design, "~Block + Group")
    expect_identical(paired$info$block_col, "Block")
    expect_identical(paired$info$n_blocks, 3L)
    expect_true(any(grepl("^\\.block_", DESeq2::resultsNames(paired$dds))))
    expect_identical(paired$info$sample_col, "SampleID")

    tab <- paired$tables$trt_vs_ctl
    expect_true(all(c("log2FoldChange", "pvalue", "padj", "FeatureID") %in% colnames(tab)))
    expect_identical(rownames(tab), rownames(unpaired$tables$trt_vs_ctl))
    # The block offsets swamp the treatment effect unless the fit accounts for them.
    n_sig <- function(t) sum(t$padj[1:40] <= 0.05, na.rm = TRUE)
    expect_gt(n_sig(tab), n_sig(unpaired$tables$trt_vs_ctl))
})

test_that("run_deseq2_de refuses a block confounded with the group", {
    skip_if_not_installed("DESeq2")
    meta <- make_rna_block_meta()
    meta$Block <- ifelse(meta$Group == "ctl", "Bx", "By")
    counts <- make_rna_block_counts(make_rna_block_meta())
    expect_error(
        suppressMessages(run_deseq2_de(counts, meta, rna_block_contrasts(),
                                       list(sample_col = "SampleID", block_col = "Block"))),
        "cannot be fitted"
    )
})

test_that("run_deseq2_de refuses block_col equal to the contrast factor", {
    skip_if_not_installed("DESeq2")
    meta <- make_rna_block_meta()
    expect_error(
        suppressMessages(run_deseq2_de(make_rna_block_counts(meta), meta,
                                       rna_block_contrasts(),
                                       list(sample_col = "SampleID", block_col = "Group"))),
        "contrast factor itself"
    )
})

test_that("run_deseq2_de passes resolver errors through (missing column)", {
    skip_if_not_installed("DESeq2")
    meta <- make_rna_block_meta()
    expect_error(
        suppressMessages(run_deseq2_de(make_rna_block_counts(meta), meta,
                                       rna_block_contrasts(),
                                       list(sample_col = "SampleID", block_col = "Donor"))),
        "not in the sample metadata"
    )
})

test_that("run_limma_rna_de records its design in info", {
    skip_if_not_installed("limma")
    meta <- make_rna_block_meta()
    expr <- log2(make_rna_block_counts(meta) + 1)
    res <- suppressMessages(run_limma_rna_de(expr, meta, rna_block_contrasts(),
                                             list(sample_col = "SampleID")))
    expect_identical(res$info$method, "limma")
    expect_identical(res$info$design, "~0 + Group")
})

# ---- validate_rna_config ----------------------------------------------------

rna_block_valid_cfg <- function() {
    list(
        id_columns = list(gene_id = "gene_id"),
        files = list(counts = "counts.csv", metadata = "metadata.csv"),
        normalization = list(method = "TMMlogCPM"),
        filtering = list(group_col = "Group", threshold = 10),
        de = list(padj_cutoff = 0.05, linear_fc_cutoff = 1.5, block_col = "Block")
    )
}

test_that("validate_rna_config accepts block_col on count input", {
    expect_true(validate_rna_config(rna_block_valid_cfg()))
})

test_that("validate_rna_config checks de$block_pca is a flag", {
    cfg <- rna_block_valid_cfg()
    cfg$de$block_pca <- FALSE
    expect_true(validate_rna_config(cfg))
    cfg$de$block_pca <- "no"
    expect_error(validate_rna_config(cfg), "block_pca")
})

test_that("validate_rna_config rejects a non-string block_col", {
    cfg <- rna_block_valid_cfg()
    cfg$de$block_col <- 3
    expect_error(validate_rna_config(cfg), "de\\$block_col")
})

test_that("validate_rna_config rejects block_col where no DESeq2 fit runs", {
    cfg <- rna_block_valid_cfg()
    cfg$files$de_table <- list(trt_vs_ctl = "de.csv")
    expect_error(validate_rna_config(cfg), "pre-computed DE tables")
})

# ---- build_rna_block_pca / plot_rna_block_pca ----------------------------------

test_that("build_rna_block_pca centres every block after correction", {
    skip_if_not_installed("limma")
    meta <- make_rna_block_meta()
    expr <- log2(make_rna_block_counts(meta) + 1)
    bp <- build_rna_block_pca(expr, meta, "SampleID", "Block", "Group", n_top = 200)
    after <- bp$coords[bp$coords$panel == levels(bp$coords$panel)[2], ]
    # Block means are removed, so each block averages to the origin; with two
    # samples per block that makes every pair mirror-image.
    expect_equal(as.vector(tapply(after$PC1, after$block, mean)), rep(0, 3), tolerance = 1e-8)
    expect_equal(as.vector(tapply(after$PC2, after$block, mean)), rep(0, 3), tolerance = 1e-8)
    expect_equal(bp$n_top, 200)
    expect_equal(nrow(bp$coords), 12)
    expect_s3_class(plot_rna_block_pca(bp, "Block", "Group"), "ggplot")
})

test_that("build_rna_block_pca keeps the within-block group difference", {
    skip_if_not_installed("limma")
    meta <- make_rna_block_meta()
    expr <- log2(make_rna_block_counts(meta) + 1)
    corrected <- limma::removeBatchEffect(expr, batch = factor(meta$Block),
                                          design = stats::model.matrix(~ factor(meta$Group)))
    d <- function(m) m[, meta$Group == "trt"] - m[, meta$Group == "ctl"]
    expect_equal(d(corrected), d(expr), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("build_rna_block_pca refuses expression columns missing from metadata", {
    skip_if_not_installed("limma")
    meta <- make_rna_block_meta()
    expr <- log2(make_rna_block_counts(meta) + 1)
    colnames(expr)[1] <- "not_a_sample"
    expect_error(build_rna_block_pca(expr, meta, "SampleID", "Block", "Group"),
                 "no row in the metadata")
})

test_that("run_deseq2_de records the contrast factor for the block PCA", {
    skip_if_not_installed("DESeq2")
    meta <- make_rna_block_meta()
    res <- suppressMessages(run_deseq2_de(make_rna_block_counts(meta), meta,
                                          rna_block_contrasts(),
                                          list(sample_col = "SampleID", block_col = "Block")))
    expect_identical(res$info$factor_col, "Group")
})

test_that("run_deseq2_de accepts a block_col that is not a syntactic name", {
    skip_if_not_installed("DESeq2")
    meta <- make_rna_block_meta()
    names(meta)[names(meta) == "Block"] <- "Donor ID"
    res <- suppressMessages(run_deseq2_de(make_rna_block_counts(make_rna_block_meta()), meta,
                                          rna_block_contrasts(),
                                          list(sample_col = "SampleID", block_col = "Donor ID")))
    expect_identical(res$info$design, "~Donor ID + Group")
    expect_identical(res$info$block_col, "Donor ID")
})

test_that("qc_post removes a stale block PCA even when it returns early", {
    out <- withr::local_tempdir()
    dirs <- create_legacy_output_dirs(out, create = TRUE)
    stale <- file.path(dirs$diagnostic_plots, c("PCA_block_corrected.png", "PCA_block_corrected_coords.tsv"))
    file.create(stale)
    m <- matrix(1, 2, 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
    pre <- list(expr_raw = m, expr_filt = m, expr_work = m,
                meta = data.frame(SampleID = c("s1", "s2")), row_data = data.frame(id = c("g1", "g2")))
    # No tables, so the summary cannot be built and the module returns early.
    res <- suppressWarnings(mod_rnaseq_qc_post(pre, list(tables = NULL),
                                               list(modes = list(rna = list())), out))
    expect_length(res$files, 0)
    expect_false(any(file.exists(stale)))
})
