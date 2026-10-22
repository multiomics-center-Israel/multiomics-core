# Tests for de$method: edger_ql (run_edger_rna_de) and its description
# test-rnaseq-de-edger.R

make_edger_meta <- function() {
    data.frame(
        SampleID = c("B1_ctl", "B1_trt", "B2_ctl", "B2_trt", "B3_ctl", "B3_trt"),
        Group    = rep(c("ctl", "trt"), times = 3),
        Block    = rep(c("B1", "B2", "B3"), each = 2),
        stringsAsFactors = FALSE
    )
}

# Synthetic only: gene-specific block effects plus a treatment effect on 40 genes.
make_edger_counts <- function(meta, n_genes = 400, seed = 13) {
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

edger_contrasts <- function() {
    data.frame(Contrast_name = "trt_vs_ctl", Factor = "Group",
               Numerator = "trt", Denominator = "ctl", stringsAsFactors = FALSE)
}

# The edgeR user's guide QL recipe, written out by hand.
edger_by_hand <- function(counts, design) {
    dge <- edgeR::estimateDisp(edgeR::calcNormFactors(edgeR::DGEList(counts), method = "TMM"), design)
    qlf <- edgeR::glmQLFTest(edgeR::glmQLFit(dge, design), coef = ncol(design))
    edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
}

test_that("run_edger_rna_de unpaired equals the edgeR QL recipe and keeps the contract", {
    skip_if_not_installed("edgeR")
    meta <- make_edger_meta(); counts <- make_edger_counts(meta)
    res <- suppressMessages(run_edger_rna_de(counts, meta, edger_contrasts(),
                                             list(sample_col = "SampleID")))
    tt <- edger_by_hand(counts, stats::model.matrix(~ factor(meta$Group, levels = c("ctl", "trt"))))
    tab <- res$tables$trt_vs_ctl
    expect_true(all(c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                      "FeatureID") %in% colnames(tab)))
    expect_equal(tab$log2FoldChange, tt$logFC, tolerance = 1e-8)
    expect_equal(tab$pvalue, tt$PValue, tolerance = 1e-8)
    expect_equal(tab$stat, sign(tt$logFC) * sqrt(tt$F), tolerance = 1e-8)
    expect_identical(res$info$method, "edger_ql")
    expect_identical(res$info$design, "~Group")
    expect_null(res$info$block_col)
})

test_that("run_edger_rna_de with a block fits ~ block + group and gains power", {
    skip_if_not_installed("edgeR")
    meta <- make_edger_meta(); counts <- make_edger_counts(meta)
    paired <- suppressMessages(run_edger_rna_de(counts, meta, edger_contrasts(),
                                                list(sample_col = "SampleID", block_col = "Block")))
    unpaired <- suppressMessages(run_edger_rna_de(counts, meta, edger_contrasts(),
                                                  list(sample_col = "SampleID")))
    tt <- edger_by_hand(counts, stats::model.matrix(~ factor(meta$Block) +
                                                     factor(meta$Group, levels = c("ctl", "trt"))))
    expect_equal(paired$tables$trt_vs_ctl$pvalue, tt$PValue, tolerance = 1e-8)
    expect_identical(paired$info$design, "~Block + Group")
    expect_identical(paired$info$block_model, "fixed")
    n_sig <- function(r) sum(r$tables$trt_vs_ctl$padj[1:40] <= 0.05, na.rm = TRUE)
    expect_gt(n_sig(paired), n_sig(unpaired))
})

test_that("run_edger_rna_de refuses a confounded block and tximport input", {
    skip_if_not_installed("edgeR")
    meta <- make_edger_meta(); counts <- make_edger_counts(meta)
    bad <- meta; bad$Block <- ifelse(bad$Group == "ctl", "Bx", "By")
    expect_error(run_edger_rna_de(counts, bad, edger_contrasts(),
                                  list(sample_col = "SampleID", block_col = "Block")),
                 "cannot be fitted")
    txi <- list(counts = counts, abundance = counts, length = counts * 0 + 1000,
                countsFromAbundance = "no")
    expect_error(run_edger_rna_de(txi, meta, edger_contrasts(), list(sample_col = "SampleID")),
                 "tximport")
})

test_that("describe_rna_de describes edgeR QL with a fixed block", {
    d <- describe_rna_de(list(method = "edger_ql", design = "~Block + Group", block_col = "Block",
                              block_model = "fixed", n_blocks = 3L, edger_version = "9.9.9"))
    expect_identical(d$label, "edgeR")
    expect_match(d$methods, "quasi-likelihood F-test", fixed = TRUE)
    expect_match(d$methods, "version 9.9.9", fixed = TRUE)
    expect_match(d$methods, "Block as a fixed term", fixed = TRUE)
})

test_that("validate_rna_config handles edger_ql", {
    cfg <- list(id_columns = list(gene_id = "gene_id"),
                files = list(counts = "counts.csv", metadata = "metadata.csv"),
                normalization = list(method = "TMMlogCPM"),
                filtering = list(group_col = "Group", threshold = 10),
                de = list(method = "edger_ql", block_col = "Block"))
    expect_true(validate_rna_config(cfg))
    cfg$de$block_model <- "random"
    expect_error(validate_rna_config(cfg), "no random effect")
    cfg$de$block_model <- NULL
    cfg$input <- list(format = "preprocessed")
    expect_error(validate_rna_config(cfg), "needs a raw count matrix")
})

test_that("edgeR tests the group, not a block column with a look-alike name", {
    skip_if_not_installed("edgeR")
    meta <- make_edger_meta()
    meta$Group <- ifelse(meta$Group == "trt", "block2", "ctl")
    meta$Block <- sub("^B", "", meta$Block)
    counts <- make_edger_counts(make_edger_meta())
    cdf <- data.frame(Contrast_name = "b2_vs_ctl", Factor = "Group",
                      Numerator = "block2", Denominator = "ctl")
    res <- suppressMessages(run_edger_rna_de(counts, meta, cdf,
                                             list(sample_col = "SampleID", block_col = "Block")))
    expect_gt(stats::median(res$tables$b2_vs_ctl$log2FoldChange[1:40]), 0.4)
    expect_identical(res$info$sample_col, "SampleID")
})

test_that("edgeR skips a contrast that compares a level with itself", {
    skip_if_not_installed("edgeR")
    meta <- make_edger_meta(); counts <- make_edger_counts(meta)
    cdf <- rbind(edger_contrasts(),
                 data.frame(Contrast_name = "same", Factor = "Group", Numerator = "trt", Denominator = "trt"))
    res <- suppressMessages(expect_warning(
        run_edger_rna_de(counts, meta, cdf, list(sample_col = "SampleID")), "with itself"))
    expect_named(res$tables, "trt_vs_ctl")
})
