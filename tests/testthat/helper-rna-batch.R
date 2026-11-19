# Synthetic fixtures for the RNA batch-correction tests
# (test-rna-batch-correction.R, test-rna-batch-qc.R).
#
# All data here is synthetic.

# Twelve samples, six per condition, in two batches that hold the conditions
# unevenly (B1: 4 Ctrl + 2 Treat, B2: 2 Ctrl + 4 Treat). `Run` is a second
# batch factor crossed with the first; `Pair` is a numeric label.
bc_meta <- function() {
    meta <- data.frame(
        SampleID  = sprintf("S%02d", 1:12),
        Condition = rep(c("Ctrl", "Treat"), each = 6),
        Batch     = c("B1", "B1", "B1", "B1", "B2", "B2",
                      "B1", "B1", "B2", "B2", "B2", "B2"),
        Run       = rep(c("R1", "R2"), times = 6),
        Pair      = rep(c(1, 2, 3), times = 4),
        stringsAsFactors = FALSE
    )
    rownames(meta) <- meta$SampleID
    meta
}

# The simulated truth on the log2 scale, 200 genes:
#   * group: genes 1-20 are +2 in Treat;
#   * Batch: genes 31-90 move 1.5 in B2, alternating up and down;
#   * Run:   genes 101-140 move 1.5 in R2, alternating up and down.
# A minority of the genes carries each effect, as in a real experiment.
bc_ids <- function() sprintf("G%03d", 1:200)
bc_group_genes <- function() bc_ids()[1:20]
bc_sign <- function(idx) {
    s <- rep(0, 200)
    s[idx] <- rep(c(1, -1), length.out = length(idx))
    s
}
bc_batch_sign <- function() bc_sign(31:90)
bc_run_sign <- function() bc_sign(101:140)

# A `pre` object as preprocess_rna() returns it, with a log-scale analysis
# matrix. `method` is what is recorded as having produced expr_work.
make_bc_pre <- function(method = "VST", seed = 21) {
    meta <- bc_meta()
    expr <- withr::with_seed(seed, {
        base <- stats::runif(200, min = 6, max = 12)
        base +
            outer(c(rep(2, 20), rep(0, 180)), as.numeric(meta$Condition == "Treat")) +
            outer(1.5 * bc_batch_sign(), as.numeric(meta$Batch == "B2")) +
            outer(1.5 * bc_run_sign(), as.numeric(meta$Run == "R2")) +
            matrix(stats::rnorm(200 * 12, sd = 0.2), nrow = 200)
    })
    dimnames(expr) <- list(bc_ids(), meta$SampleID)
    counts <- round(2^expr)
    storage.mode(counts) <- "integer"
    attr(expr, "method") <- method
    pre <- list(
        expr_raw  = counts,
        expr_filt = counts,
        expr_work = expr,
        de_input  = counts,
        row_data  = data.frame(gene = bc_ids(), stringsAsFactors = FALSE),
        meta      = meta,
        info      = list(expr_work_method = method),
        contrasts = data.frame(Contrast_name = "Treat_vs_Ctrl", Factor = "Condition",
                               Numerator = "Treat", Denominator = "Ctrl",
                               stringsAsFactors = FALSE)
    )
    attr(pre, "source_type") <- "preprocessed"
    pre
}

bc_config <- function(enabled = TRUE, method = "limma", cols = "Batch", ...) {
    list(
        params = list(seed = 1),
        modes = list(rna = list(
            id_columns = list(sample_col = "SampleID"),
            effects = list(samples = "SampleID", color = "Condition"),
            # The 3D PCA needs pandoc to save; keep it out of the unit tests.
            qc = list(thresholds = list(min_samples_for_pca3d = 100)),
            de = list(p_cutoff = 0.05),
            batch_correction = list(enabled = enabled, method = method,
                                    batch_cols = cols, ...)
        ))
    )
}

# Mean difference between the two levels of a metadata column, over `genes`,
# each gene taken along `sign` (its simulated direction).
bc_shift <- function(expr, col, genes, sign = rep(1, length(genes))) {
    lv <- sort(unique(bc_meta()[[col]]))
    hi <- bc_meta()[[col]] == lv[2]
    mean(sign * (rowMeans(expr[genes, hi]) - rowMeans(expr[genes, !hi])))
}
# Genes that carry a batch effect and no group effect.
bc_batch_genes <- function() bc_ids()[bc_batch_sign() != 0]
bc_run_genes <- function() bc_ids()[bc_run_sign() != 0]
bc_plain <- function(x) {
    attr(x, "method") <- NULL
    x
}
