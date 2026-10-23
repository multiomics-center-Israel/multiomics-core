# Tests for de$method: limma_voom (run_voom_rna_de) and its description
# test-rnaseq-de-voom.R

# Six libraries, three blocks, both groups in every block. Synthetic only.
make_voom_meta <- function() {
    data.frame(
        SampleID = c("B1_ctl", "B1_trt", "B2_ctl", "B2_trt", "B3_ctl", "B3_trt"),
        Group    = rep(c("ctl", "trt"), times = 3),
        Block    = rep(c("B1", "B2", "B3"), each = 2),
        stringsAsFactors = FALSE
    )
}

# Gene-specific block effects (not absorbed by TMM) plus a consistent
# treatment effect on the first 40 genes.
make_voom_counts <- function(meta, n_genes = 400, seed = 11) {
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

voom_contrasts <- function() {
    data.frame(Contrast_name = "trt_vs_ctl", Factor = "Group",
               Numerator = "trt", Denominator = "ctl", stringsAsFactors = FALSE)
}

skip_unless_voom <- function() {
    skip_if_not_installed("limma")
    skip_if_not_installed("edgeR")
}

test_that("run_voom_rna_de without a block returns the RNA DE table contract", {
    skip_unless_voom()
    meta <- make_voom_meta()
    res <- suppressMessages(run_voom_rna_de(make_voom_counts(meta), meta, voom_contrasts(),
                                            list(sample_col = "SampleID")))
    expect_null(res$dds)
    expect_identical(res$info$method, "limma_voom")
    expect_identical(res$info$design, "~Group")
    expect_null(res$info$block_col)
    tab <- res$tables$trt_vs_ctl
    expect_true(all(c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                      "FeatureID") %in% colnames(tab)))
    expect_identical(tab$FeatureID, rownames(tab))
    expect_true(all(is.finite(tab$lfcSE)))
    # Treated up-regulated genes come out with a positive log2FC (trt over ctl).
    expect_gt(stats::median(tab$log2FoldChange[1:40]), 0)
})

test_that("run_voom_rna_de blocks via duplicateCorrelation and gains power", {
    skip_unless_voom()
    meta <- make_voom_meta()
    counts <- make_voom_counts(meta)
    blocked <- suppressMessages(run_voom_rna_de(counts, meta, voom_contrasts(),
                                                list(sample_col = "SampleID", block_col = "Block")))
    unblocked <- suppressMessages(run_voom_rna_de(counts, meta, voom_contrasts(),
                                                  list(sample_col = "SampleID")))
    expect_identical(blocked$info$block_col, "Block")
    expect_identical(blocked$info$n_blocks, 3L)
    expect_true(is.finite(blocked$info$block_correlation))
    expect_gt(blocked$info$block_correlation, 0)
    # The design has no Line term: the block is a random effect.
    expect_identical(blocked$info$design, "~Group")
    n_sig <- function(r) sum(r$tables$trt_vs_ctl$padj[1:40] <= 0.05, na.rm = TRUE)
    expect_gt(n_sig(blocked), n_sig(unblocked))
})

test_that("run_voom_rna_de equals a direct single-coefficient limma fit", {
    skip_unless_voom()
    meta <- make_voom_meta()
    counts <- make_voom_counts(meta)
    res <- suppressMessages(run_voom_rna_de(counts, meta, voom_contrasts(),
                                            list(sample_col = "SampleID", block_col = "Block")))
    # The limma user's guide recipe, written out by hand with trt as a coefficient.
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts), method = "TMM")
    des <- stats::model.matrix(~ factor(meta$Group, levels = c("ctl", "trt")))
    blk <- factor(meta$Block)
    v <- limma::voom(dge, des)
    cor <- limma::duplicateCorrelation(v, des, block = blk)$consensus
    v <- limma::voom(dge, des, block = blk, correlation = cor)
    cor <- limma::duplicateCorrelation(v, des, block = blk)$consensus
    fit <- limma::eBayes(limma::lmFit(v, des, block = blk, correlation = cor))
    tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
    expect_equal(res$info$block_correlation, cor, tolerance = 1e-8)
    expect_equal(res$tables$trt_vs_ctl$stat, tt$t, tolerance = 1e-8)
    expect_equal(res$tables$trt_vs_ctl$pvalue, tt$P.Value, tolerance = 1e-8)
})

test_that("run_voom_rna_de fits a contrast between two non-reference levels", {
    skip_unless_voom()
    meta <- rbind(make_voom_meta(), transform(make_voom_meta()[c(1, 3, 5), ],
                                               SampleID = paste0("B", 1:3, "_alt"), Group = "alt"))
    counts <- make_voom_counts(meta)
    cdf <- data.frame(Contrast_name = c("trt_vs_alt", "alt_vs_ctl"), Factor = "Group",
                      Numerator = c("trt", "alt"), Denominator = c("alt", "ctl"))
    res <- suppressMessages(run_voom_rna_de(counts, meta, cdf, list(sample_col = "SampleID")))
    expect_named(res$tables, c("trt_vs_alt", "alt_vs_ctl"))
    expect_true(all(is.finite(res$tables$trt_vs_alt$log2FoldChange)))
})

test_that("run_voom_rna_de aligns metadata to the count columns", {
    skip_unless_voom()
    meta <- make_voom_meta()
    counts <- make_voom_counts(meta)
    shuffled <- meta[c(4, 1, 6, 2, 5, 3), ]
    a <- suppressMessages(run_voom_rna_de(counts, meta, voom_contrasts(),
                                          list(sample_col = "SampleID")))
    b <- suppressMessages(run_voom_rna_de(counts, shuffled, voom_contrasts(),
                                          list(sample_col = "SampleID")))
    expect_equal(a$tables$trt_vs_ctl$log2FoldChange, b$tables$trt_vs_ctl$log2FoldChange)
})

test_that("run_voom_rna_de refuses tximport input", {
    skip_unless_voom()
    meta <- make_voom_meta()
    counts <- make_voom_counts(meta)
    txi <- list(counts = counts, abundance = counts, length = counts * 0 + 1000,
                countsFromAbundance = "no")
    expect_error(run_voom_rna_de(txi, meta, voom_contrasts(), list(sample_col = "SampleID")),
                 "tximport")
})

test_that("describe_rna_de describes limma-voom with and without a block", {
    plain <- describe_rna_de(list(method = "limma_voom", design = "~0 + Group"))
    expect_identical(plain$label, "limma-voom")
    expect_match(plain$methods, "voom precision weights", fixed = TRUE)
    expect_false(grepl("duplicateCorrelation", plain$methods))

    blk <- describe_rna_de(list(method = "limma_voom", design = "~0 + Group",
                                block_col = "Block", n_blocks = 3L,
                                block_correlation = 0.5, block_requested = "Block"))
    expect_match(blk$methods, "random effect", fixed = TRUE)
    expect_match(blk$methods, "consensus 0.500", fixed = TRUE)
    expect_match(blk$short, "block Block", fixed = TRUE)

    fell_back <- describe_rna_de(list(method = "limma_voom", design = "~0 + Group",
                                      block_col = NULL, block_requested = "Block"))
    expect_match(fell_back$methods, "fitted without blocking", fixed = TRUE)
})

voom_valid_cfg <- function() {
    list(
        id_columns = list(gene_id = "gene_id"),
        files = list(counts = "counts.csv", metadata = "metadata.csv"),
        normalization = list(method = "TMMlogCPM"),
        filtering = list(group_col = "Group", threshold = 10),
        de = list(method = "limma_voom", padj_cutoff = 0.05, linear_fc_cutoff = 1.5,
                  block_col = "Block")
    )
}

test_that("validate_rna_config accepts limma_voom with a block on counts", {
    expect_true(validate_rna_config(voom_valid_cfg()))
})

test_that("validate_rna_config rejects limma_voom on tximport input", {
    cfg <- voom_valid_cfg()
    cfg$files$txi <- "txi.rds"
    cfg$normalization$method <- "VST"   # tximport rejects TMMlogCPM earlier
    expect_error(validate_rna_config(cfg), "needs a raw count matrix")
})

test_that("validate_rna_config warns on a de$method that is not a count method", {
    cfg <- voom_valid_cfg()
    cfg$de$method <- "limma"
    expect_warning(validate_rna_config(cfg), "DESeq2 will be fitted")
})

# ---- de$block_model -----------------------------------------------------------

test_that("resolve_rna_block_model defaults by method and honours the key", {
    expect_null(resolve_rna_block_model(list(method = "limma_voom")))
    expect_identical(resolve_rna_block_model(list(block_col = "Block")), "fixed")
    expect_identical(resolve_rna_block_model(list(method = "limma_voom", block_col = "Block")), "random")
    expect_identical(resolve_rna_block_model(list(method = "limma_voom", block_col = "Block",
                                                  block_model = "Fixed")), "fixed")
})

test_that("run_voom_rna_de with block_model fixed equals limma ~ block + group by hand", {
    skip_unless_voom()
    meta <- make_voom_meta()
    counts <- make_voom_counts(meta)
    res <- suppressMessages(run_voom_rna_de(counts, meta, voom_contrasts(),
                                            list(sample_col = "SampleID", block_col = "Block",
                                                 block_model = "fixed")))
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts), method = "TMM")
    des <- stats::model.matrix(~ factor(meta$Block) + factor(meta$Group, levels = c("ctl", "trt")))
    fit <- limma::eBayes(limma::lmFit(limma::voom(dge, des), des))
    tt <- limma::topTable(fit, coef = ncol(des), number = Inf, sort.by = "none")
    expect_identical(res$info$block_model, "fixed")
    expect_identical(res$info$design, "~Block + Group")
    expect_true(is.na(res$info$block_correlation))
    expect_equal(res$tables$trt_vs_ctl$stat, tt$t, tolerance = 1e-8)
    expect_equal(res$tables$trt_vs_ctl$log2FoldChange, tt$logFC, tolerance = 1e-8)
})

test_that("run_voom_rna_de with block_model fixed refuses a confounded block", {
    skip_unless_voom()
    meta <- make_voom_meta()
    counts <- make_voom_counts(meta)
    meta$Block <- ifelse(meta$Group == "ctl", "Bx", "By")
    expect_error(run_voom_rna_de(counts, meta, voom_contrasts(),
                                 list(sample_col = "SampleID", block_col = "Block",
                                      block_model = "fixed")),
                 "relabelling")   # caught before the rank check, which would also stop it
})

test_that("describe_rna_de names a fixed block in limma-voom", {
    d <- describe_rna_de(list(method = "limma_voom", design = "~Block + Group",
                              block_col = "Block", block_model = "fixed", n_blocks = 3L))
    expect_match(d$methods, "Block as a fixed term", fixed = TRUE)
    expect_false(grepl("duplicateCorrelation", d$methods))
    expect_match(d$short, "fixed term", fixed = TRUE)
    expect_identical(d$block_model, "fixed")
})

test_that("validate_rna_config rejects a random block for DESeq2 and bad values", {
    cfg <- voom_valid_cfg()
    cfg$de$block_model <- "fixed"
    expect_true(validate_rna_config(cfg))
    cfg$de$method <- "deseq2"; cfg$de$block_model <- "random"
    expect_error(validate_rna_config(cfg), "no random effect")
    cfg$de$block_model <- "mixed"
    expect_error(validate_rna_config(cfg), "block_model")
})

# ---- review fixes ---------------------------------------------------------------

test_that("voom fixed tests the group, not a block column with a look-alike name", {
    skip_unless_voom()
    meta <- make_voom_meta()
    # Group level "block2" reads like the model.matrix name of block level "2".
    meta$Group <- ifelse(meta$Group == "trt", "block2", "ctl")
    meta$Block <- sub("^B", "", meta$Block)
    counts <- make_voom_counts(make_voom_meta())
    cdf <- data.frame(Contrast_name = "b2_vs_ctl", Factor = "Group",
                      Numerator = "block2", Denominator = "ctl")
    res <- suppressMessages(run_voom_rna_de(counts, meta, cdf,
                                            list(sample_col = "SampleID", block_col = "Block",
                                                 block_model = "fixed")))
    # The first 40 genes carry the treatment effect; the block coefficient would not.
    expect_gt(stats::median(res$tables$b2_vs_ctl$log2FoldChange[1:40]), 0.4)
})

test_that("voom random refuses a block that relabels the contrast factor", {
    skip_unless_voom()
    meta <- make_voom_meta(); counts <- make_voom_counts(meta)
    expect_error(run_voom_rna_de(counts, meta, voom_contrasts(),
                                 list(sample_col = "SampleID", block_col = "Group")),
                 "relabelling")
    meta$Pair <- ifelse(meta$Group == "ctl", "P1", "P2")
    expect_error(run_voom_rna_de(counts, meta, voom_contrasts(),
                                 list(sample_col = "SampleID", block_col = "Pair")),
                 "relabelling")
})

test_that("voom skips a self-contrast and refuses levels that collide after make.names", {
    skip_unless_voom()
    meta <- make_voom_meta(); counts <- make_voom_counts(meta)
    cdf <- rbind(voom_contrasts(),
                 data.frame(Contrast_name = "same", Factor = "Group", Numerator = "ctl", Denominator = "ctl"))
    res <- suppressMessages(expect_warning(
        run_voom_rna_de(counts, meta, cdf, list(sample_col = "SampleID")), "with itself"))
    expect_named(res$tables, "trt_vs_ctl")
    bad <- meta; bad$Group <- ifelse(bad$Group == "ctl", "A B", "A.B")
    expect_error(run_voom_rna_de(counts, bad, data.frame(Contrast_name = "x", Factor = "Group",
                                                         Numerator = "A B", Denominator = "A.B"),
                                 list(sample_col = "SampleID")),
                 "become identical")
})
