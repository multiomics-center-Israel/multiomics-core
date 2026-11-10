# Tests for batch in the RNA-seq DE model (R/domain/rnaseq/04b_de_batch.R,
# run_deseq2_de() and run_limma_rna_de())
# test-rna-de-batch-model.R
#
# All data here is synthetic.

# Twelve samples, six per condition, in two batches that hold the conditions
# unevenly (B1: 4 Ctrl + 2 Treat, B2: 2 Ctrl + 4 Treat), so an unmodelled
# batch effect leaks into the condition contrast. `Run` is a second batch
# factor crossed with the first; `Pair` is a numeric label with three values.
make_batch_meta <- function() {
    data.frame(
        SampleID  = sprintf("S%02d", 1:12),
        Condition = rep(c("Ctrl", "Treat"), each = 6),
        Batch     = c("B1", "B1", "B1", "B1", "B2", "B2",
                      "B1", "B1", "B2", "B2", "B2", "B2"),
        Run       = rep(c("R1", "R2"), times = 6),
        Pair      = rep(c(1, 2, 3), times = 4),
        stringsAsFactors = FALSE
    )
}

# The simulated truth. DESeq2's median-of-ratios size factors assume that most
# genes do not differ between samples and that those which do are balanced; a
# fixture that breaks this biases every log2FC through the size factors, in
# models with and without batch alike. So:
#   * DE: genes 1-20 go up 4-fold in Treat, genes 21-40 down 4-fold.
#   * Batch: 90 of the 300 genes (31-120) move 2^1.5-fold in B2, alternating
#     up and down. Genes 31-40 are DE and batch-affected.
#   * The other 190 genes are unchanged and anchor the size factors.
gene_ids <- function(n_genes = 300) sprintf("G%03d", seq_len(n_genes))
truth_lfc <- function(n_genes = 300) c(rep(2, 20), rep(-2, 20), rep(0, n_genes - 40))
batch_sign <- function(n_genes = 300) {
    s <- rep(0, n_genes)
    s[31:120] <- rep(c(1, -1), length.out = 90)
    s
}

make_batch_counts <- function(meta, n_genes = 300, batch_lfc = 1.5, seed = 11) {
    withr::with_seed(seed, {
        base <- stats::rgamma(n_genes, shape = 4, rate = 0.004) + 200
        log2mu <- outer(log2(base), rep(0, nrow(meta)), "+") +
            outer(truth_lfc(n_genes), as.numeric(meta$Condition == "Treat")) +
            outer(batch_lfc * batch_sign(n_genes), as.numeric(meta$Batch == "B2"))
        counts <- matrix(stats::rnbinom(length(log2mu), mu = 2^log2mu, size = 100),
                         nrow = n_genes)
    })
    dimnames(counts) <- list(gene_ids(n_genes), meta$SampleID)
    counts
}

# Estimation error of the Treat vs Ctrl log2FC against the simulated truth.
lfc_error <- function(tab, genes = gene_ids()) {
    tab[genes, "log2FoldChange"] - truth_lfc()[match(genes, gene_ids())]
}
rmse <- function(e) sqrt(mean(e^2))
batch_genes <- function() gene_ids()[batch_sign() != 0]

# Diagnostic: the mean error of the batch-affected null genes along each gene's
# own batch direction. Unbalanced batches leak about 1.5 * (4/6 - 2/6) = 0.5
# when batch is not modelled. Blind to a shift common to all genes (it cancels
# between up and down batch genes), which is why the bias check below exists.
batch_leak <- function(tab) {
    genes <- gene_ids()[batch_sign() != 0 & truth_lfc() == 0]
    mean(lfc_error(tab, genes) * batch_sign()[match(genes, gene_ids())])
}

# The scientific claim, checked for one fit with and one without batch:
#   * on the batch-affected genes, modelling batch at least halves the RMSE
#     against the truth (expected: about 0.5 -> about 0.1, a 4-5x reduction);
#   * the batch-adjusted estimates have no global bias (expected about 0;
#     a size-factor or contrast error shows up here);
#   * every DE gene keeps the direction of its true effect (|truth| = 2,
#     about 20 standard errors);
#   * the leak diagnostic agrees.
expect_batch_adjustment_improves <- function(off_tab, on_tab, info) {
    expect_lt(rmse(lfc_error(on_tab, batch_genes())),
              0.5 * rmse(lfc_error(off_tab, batch_genes())), label = paste(info, "RMSE"))
    expect_lt(abs(mean(lfc_error(on_tab))), 0.1, label = paste(info, "global bias"))
    de <- gene_ids()[truth_lfc() != 0]
    expect_true(all(sign(on_tab[de, "log2FoldChange"]) == sign(truth_lfc()[truth_lfc() != 0])),
                label = paste(info, "DE directions"))
    expect_gt(batch_leak(off_tab), 0.3, label = paste(info, "leak without batch"))
    expect_lt(abs(batch_leak(on_tab)), 0.1, label = paste(info, "leak with batch"))
}

# Seeds the estimation checks are repeated over, so they are not tuned to one.
batch_seeds <- c(11, 23, 37, 101, 2024)

batch_contrasts <- function() {
    data.frame(Contrast_name = "Treat_vs_Ctrl", Factor = "Condition",
               Numerator = "Treat", Denominator = "Ctrl", stringsAsFactors = FALSE)
}

batch_on <- function(cols = "Batch") {
    get_rna_batch_model_config(list(batch_in_model = list(enabled = TRUE, batch_cols = cols)))
}

# --- Config resolution and validation ----------------------------------------

test_that("get_rna_batch_model_config is off unless enabled, and checks batch_cols", {
    off <- get_rna_batch_model_config(list())
    expect_false(off$enabled)
    expect_identical(off$batch_cols, character(0))
    # batch_cols are ignored, not validated, while disabled
    expect_false(get_rna_batch_model_config(list(batch_in_model = list(
        enabled = FALSE, batch_cols = c("a", "b", "c"))))$enabled)

    expect_identical(batch_on("Batch")$batch_cols, "Batch")
    expect_identical(batch_on(c("Batch", "Run"))$batch_cols, c("Batch", "Run"))
    for (bad in list(character(0), c("A", "B", "C"), c("Batch", "Batch"), "")) {
        expect_error(batch_on(bad), "batch_cols")
    }
})

test_that("validate_rna_config refuses batch in the model with precomputed DE", {
    cfg <- list(
        id_columns = list(gene_id = "gene"),
        files = list(counts = "counts.tsv", metadata = "meta.tsv", de_table = "de.csv"),
        batch_in_model = list(enabled = TRUE, batch_cols = "Batch")
    )
    expect_error(validate_rna_config(cfg), "precomputed")
    cfg$batch_in_model$enabled <- FALSE
    expect_true(validate_rna_config(cfg))
    cfg$files$de_table <- NULL
    cfg$batch_in_model$enabled <- TRUE
    expect_true(validate_rna_config(cfg))
})

# --- Design checks (no DE package needed) -------------------------------------

test_that("prepare_rna_batch_design makes numeric batch labels categorical", {
    meta <- make_batch_meta()
    out <- prepare_rna_batch_design(meta, meta$SampleID, "SampleID", "Condition", "Pair")
    expect_s3_class(out$Pair, "factor")
    expect_identical(levels(out$Pair), c("1", "2", "3"))
    # Three pair levels -> two coefficients, not one slope.
    x <- stats::model.matrix(rna_batch_design_formula("Condition", "Pair"), data = out)
    expect_identical(ncol(x), 4L)
})

test_that("prepare_rna_batch_design stops on unusable batch columns", {
    meta <- make_batch_meta()
    prep <- function(m, cols, ids = m$SampleID) {
        prepare_rna_batch_design(m, ids, "SampleID", "Condition", cols)
    }
    expect_error(prep(meta, "Plate"), "not found")
    expect_error(prep(meta, "Condition"), "group column")

    na_meta <- meta
    na_meta$Batch[3] <- NA
    expect_error(prep(na_meta, "Batch"), "no value")

    one <- meta
    one$Batch <- "B1"
    expect_error(prep(one, "Batch"), "single level")

    odd <- meta
    names(odd)[names(odd) == "Batch"] <- "my batch"
    expect_error(prep(odd, "my batch"), "formula")

    conf <- meta
    conf$Batch <- ifelse(conf$Condition == "Ctrl", "B1", "B2")
    expect_error(prep(conf, "Batch"), "confounded")

    # The second batch column is checked against the first as well.
    nested <- meta
    nested$Run <- nested$Batch
    expect_error(prep(nested, c("Batch", "Run")), "confounded")
})

test_that("prepare_rna_batch_design requires residual degrees of freedom", {
    # Full rank, but 4 coefficients for 4 samples.
    meta <- data.frame(SampleID = paste0("S", 1:4), Condition = c("A", "A", "B", "B"),
                       B1 = c("x", "y", "x", "y"), B2 = c("p", "q", "q", "p"))
    expect_error(prepare_rna_batch_design(meta, meta$SampleID, "SampleID", "Condition",
                                          c("B1", "B2")),
                 "residual degrees")
})

test_that("prepare_rna_batch_design judges only the samples that enter DE", {
    meta <- make_batch_meta()
    # An extra sample outside DE with no batch label must not stop the run.
    extra <- rbind(meta, data.frame(SampleID = "S99", Condition = "Ctrl", Batch = NA,
                                    Run = "R1", Pair = 1))
    expect_silent(prepare_rna_batch_design(extra, meta$SampleID, "SampleID",
                                           "Condition", "Batch"))
})

test_that("build_rna_limma_design keeps the group columns and adds batch columns", {
    meta <- prepare_rna_batch_design(make_batch_meta(), sprintf("S%02d", 1:12),
                                     "SampleID", "Condition", c("Batch", "Pair"))
    group <- factor(meta$Condition)
    d <- build_rna_limma_design(group, meta, c("Batch", "Pair"))
    expect_identical(colnames(d), c("Ctrl", "Treat", "Batch_B2", "Pair_2", "Pair_3"))
    expect_identical(nrow(d), 12L)
})

# --- limma (preprocessed) path -----------------------------------------------

limma_expr <- function(meta, seed = 11) {
    log2(make_batch_counts(meta, seed = seed) + 1)
}

test_that("run_limma_rna_de without batch matches the current limma fit", {
    meta <- make_batch_meta()
    expr <- limma_expr(meta)
    de_cfg <- list(sample_col = "SampleID")

    off   <- suppressMessages(run_limma_rna_de(expr, meta, batch_contrasts(), de_cfg))
    off_2 <- suppressMessages(run_limma_rna_de(expr, meta, batch_contrasts(), de_cfg,
                                               batch_cfg = get_rna_batch_model_config(list())))
    expect_identical(off$tables, off_2$tables)

    # The fit as written before batch_in_model existed.
    cond <- factor(meta$Condition)
    design <- stats::model.matrix(~ 0 + cond)
    colnames(design) <- levels(cond)
    fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr, design),
                         limma::makeContrasts(Treat_vs_Ctrl = Treat - Ctrl, levels = design)))
    ref <- limma::topTable(fit, coef = 1, sort.by = "none", number = Inf)
    tab <- off$tables$Treat_vs_Ctrl
    expect_equal(tab$log2FoldChange, ref$logFC)
    expect_equal(tab$pvalue, ref$P.Value)
    expect_equal(tab$padj, ref$adj.P.Val)

    expect_identical(off$design$formula, "~ 0 + Condition")
    expect_identical(off$design$batch_cols, "")
})

test_that("run_limma_rna_de with batch reduces the estimation error against the truth", {
    meta <- make_batch_meta()
    de_cfg <- list(sample_col = "SampleID")
    for (seed in batch_seeds) {
        expr <- limma_expr(meta, seed = seed)
        off <- suppressMessages(run_limma_rna_de(expr, meta, batch_contrasts(), de_cfg))
        on  <- suppressMessages(run_limma_rna_de(expr, meta, batch_contrasts(), de_cfg,
                                                 batch_cfg = batch_on("Batch")))
        expect_batch_adjustment_improves(off$tables$Treat_vs_Ctrl, on$tables$Treat_vs_Ctrl,
                                         info = paste("limma seed", seed))
    }

    expect_identical(on$design$formula, "~ 0 + Condition + Batch")
    expect_identical(on$design$batch_cols, "Batch")
    expect_identical(names(on$tables), names(off$tables))
    expect_identical(colnames(on$tables$Treat_vs_Ctrl), colnames(off$tables$Treat_vs_Ctrl))
})

# --- DESeq2 path ----------------------------------------------------------------

run_ds <- function(meta, counts, mode = "default", batch_cfg = NULL) {
    suppressWarnings(suppressMessages(run_deseq2_de(
        counts, meta, batch_contrasts(),
        de_cfg = list(sample_col = "SampleID", deseq_mode = mode),
        batch_cfg = batch_cfg)))
}

test_that("run_deseq2_de without batch matches the current DESeq2 fit", {
    skip_if_not_installed("DESeq2")
    meta <- make_batch_meta()
    counts <- make_batch_counts(meta)

    off   <- run_ds(meta, counts)
    off_2 <- run_ds(meta, counts, batch_cfg = get_rna_batch_model_config(list()))
    expect_equal(off$tables, off_2$tables)
    expect_identical(deparse(DESeq2::design(off$dds)), "~Condition")

    # The fit as written before batch_in_model existed.
    cd <- meta
    cd$Condition <- factor(cd$Condition)
    rownames(cd) <- cd$SampleID
    ref_dds <- suppressMessages(DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(
        counts, colData = cd, design = ~ Condition)))
    ref <- as.data.frame(DESeq2::results(ref_dds, contrast = c("Condition", "Treat", "Ctrl"),
                                         alpha = 0.05))
    tab <- off$tables$Treat_vs_Ctrl
    expect_equal(tab$log2FoldChange, ref$log2FoldChange)
    expect_equal(tab$pvalue, ref$pvalue)
    expect_equal(tab$padj, ref$padj)

    expect_identical(off$design$formula, "~Condition")
    expect_identical(off$design$batch_cols, "")
})

test_that("run_deseq2_de with batch fits ~ Batch + Condition and reduces the error", {
    skip_if_not_installed("DESeq2")
    meta <- make_batch_meta()
    for (seed in batch_seeds) {
        counts <- make_batch_counts(meta, seed = seed)
        off <- run_ds(meta, counts)
        on  <- run_ds(meta, counts, batch_cfg = batch_on("Batch"))
        expect_batch_adjustment_improves(off$tables$Treat_vs_Ctrl, on$tables$Treat_vs_Ctrl,
                                         info = paste("DESeq2 seed", seed))
    }

    expect_identical(deparse(DESeq2::design(on$dds)), "~Batch + Condition")
    expect_identical(on$design$formula, "~Batch + Condition")
    expect_identical(on$design$batch_cols, "Batch")
    expect_identical(on$design$residual_df, 9L)

    # Output contract unchanged: same contrasts, same columns.
    expect_identical(names(on$tables), names(off$tables))
    expect_identical(colnames(on$tables$Treat_vs_Ctrl), colnames(off$tables$Treat_vs_Ctrl))
})

test_that("results(dds) without a contrast reports the condition in both DESeq2 modes", {
    skip_if_not_installed("DESeq2")
    meta <- make_batch_meta()
    counts <- make_batch_counts(meta)
    for (mode in c("default", "legacy")) {
        on <- run_ds(meta, counts, mode = mode, batch_cfg = batch_on("Batch"))
        desc <- S4Vectors::mcols(suppressMessages(DESeq2::results(on$dds)))$description[2]
        expect_match(desc, "Condition Treat vs Ctrl", fixed = TRUE, info = mode)
        expect_identical(on$design$mode, mode)

        # The explicit contrast the pipeline uses is the same comparison.
        explicit <- suppressMessages(DESeq2::results(
            on$dds, contrast = c("Condition", "Treat", "Ctrl"), alpha = 0.05))
        expect_equal(on$tables$Treat_vs_Ctrl$log2FoldChange, explicit$log2FoldChange,
                     info = mode)
    }
})

test_that("run_deseq2_de accepts two batch columns and numeric batch labels", {
    skip_if_not_installed("DESeq2")
    meta <- make_batch_meta()
    counts <- make_batch_counts(meta)

    two <- run_ds(meta, counts, batch_cfg = batch_on(c("Batch", "Run")))
    expect_identical(deparse(DESeq2::design(two$dds)), "~Batch + Run + Condition")
    expect_identical(two$design$batch_cols, "Batch;Run")

    num <- run_ds(meta, counts, batch_cfg = batch_on("Pair"))
    expect_s3_class(SummarizedExperiment::colData(num$dds)$Pair, "factor")
    expect_identical(num$design$n_coefficients, 4L)
})

test_that("run_deseq2_de stops on a confounded batch before DESeq2 does", {
    skip_if_not_installed("DESeq2")
    meta <- make_batch_meta()
    counts <- make_batch_counts(meta)
    meta$Batch <- ifelse(meta$Condition == "Ctrl", "B1", "B2")
    expect_error(run_ds(meta, counts, batch_cfg = batch_on("Batch")), "confounded")
})

# --- Methods sentence from the recorded design ---------------------------------

test_that("rna_de_design_sentence describes only a batch that was fitted", {
    rec <- function(formula, cols, engine = "DESeq2", mode = "default") {
        data.frame(engine = engine, mode = mode, formula = formula, batch_cols = cols,
                   stringsAsFactors = FALSE)
    }
    expect_null(rna_de_design_sentence(NULL))
    expect_null(rna_de_design_sentence(rec("~Condition", "")))
    expect_null(rna_de_design_sentence(rec("~Condition", NA)))

    one <- rna_de_design_sentence(rec("~Batch + Condition", "Batch"))
    expect_match(one, "The DESeq2 model included batch as a covariate", fixed = TRUE)
    expect_match(one, "`~Batch + Condition`; batch column Batch)", fixed = TRUE)
    expect_false(grepl("legacy", one, fixed = TRUE))

    two <- rna_de_design_sentence(rec("~Batch + Run + Condition", "Batch;Run"))
    expect_match(two, "batch columns Batch and Run", fixed = TRUE)

    legacy <- rna_de_design_sentence(rec("~Batch + Condition", "Batch", mode = "legacy"))
    expect_match(legacy, "In legacy mode (betaPrior)", fixed = TRUE)

    lim <- rna_de_design_sentence(rec("~ 0 + Condition + Batch", "Batch",
                                      engine = "limma", mode = NA))
    expect_match(lim, "The limma model included batch", fixed = TRUE)
})

test_that("de_design.tsv always reflects the current run in a reused output directory", {
    dir <- withr::local_tempdir()
    path <- file.path(dir, "de_design.tsv")
    # What report_template.Rmd does: the sentence from the file if it exists.
    methods_text <- function() {
        if (file.exists(path)) rna_de_design_sentence(read.delim(path, stringsAsFactors = FALSE))
    }
    x <- stats::model.matrix(~ Batch + Condition, data = make_batch_meta())
    batch_run <- list(design = rna_de_design_record(
        "DESeq2", "default", "~Batch + Condition", "Batch", x))

    # Run 1 fits batch: the file is written and Methods describes it.
    expect_identical(write_rna_de_design(batch_run, dir), path)
    expect_match(methods_text(), "batch column Batch", fixed = TRUE)

    # Run 2, same directory, precomputed DE: no model, so the old file must go.
    expect_identical(write_rna_de_design(list(dds = NULL, tables = list()), dir), character(0))
    expect_false(file.exists(path))
    expect_null(methods_text())
    # Nothing to remove is not an error.
    expect_identical(write_rna_de_design(list(), dir), character(0))

    # Run 3 fits batch again, then run 4 fits a model without batch.
    write_rna_de_design(batch_run, dir)
    no_batch_run <- list(design = rna_de_design_record(
        "DESeq2", "default", "~Condition", character(0),
        stats::model.matrix(~ Condition, data = make_batch_meta())))
    write_rna_de_design(no_batch_run, dir)
    expect_true(file.exists(path))
    expect_identical(read.delim(path, stringsAsFactors = FALSE)$formula, "~Condition")
    expect_null(methods_text())
})

test_that("de_design.tsv round-trips into the Methods sentence", {
    dir <- withr::local_tempdir()
    meta <- make_batch_meta()
    x <- stats::model.matrix(~ Batch + Condition, data = meta)
    for (cols in list(character(0), "Batch")) {
        rec <- rna_de_design_record("DESeq2", "default", "~Batch + Condition", cols, x)
        path <- save_tsv(rec, dir, "de_design.tsv")
        back <- read.delim(path, stringsAsFactors = FALSE)
        if (length(cols) == 0) {
            expect_null(rna_de_design_sentence(back))
        } else {
            expect_match(rna_de_design_sentence(back), "batch column Batch", fixed = TRUE)
        }
    }
})
