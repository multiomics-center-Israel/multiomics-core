# tests/testthat/test-rna-top-variable-pca.R
#
# The top-variable PCA (PCA_top500/1000/2000/5000.png) ranks only genes whose
# mean TMM-normalized CPM reaches qc$top_variable_pca$min_mean_cpm (default 1).
# On log-CPM with prior.count = 1, genes with a handful of reads can reach a
# large variance from counting noise alone and crowd out expressed genes. The
# gate changes those plots only: expr_work, the main PCA and everything the
# module hands downstream must stay as they were, and min_mean_cpm = 0 must
# give back the previous selection exactly. All data here are synthetic.

# 700 expressed genes (half shifted between groups) and 500 genes with a few
# reads. Libraries are ~1M reads, so one read is ~1 CPM and the few-read genes
# sit well below a mean of 1 CPM while their log-CPM variance is among the
# highest. Sample S4 carries a block of very high counts so that TMM factors
# move away from 1.
make_top_var_counts <- function() {
    withr::with_seed(42, {
        hi <- matrix(rpois(700 * 6, lambda = rep(runif(700, 300, 3000), 6)), nrow = 700)
        hi[1:350, 4:6] <- round(hi[1:350, 4:6] * 1.8)
        lo <- matrix(rpois(500 * 6, lambda = 0.3), nrow = 500)
    })
    counts <- rbind(hi, lo)
    counts[1:40, 4] <- counts[1:40, 4] * 25
    counts <- counts[rowSums(counts) > 0, , drop = FALSE]
    dimnames(counts) <- list(paste0("G", seq_len(nrow(counts))), paste0("S", 1:6))
    counts
}

make_top_var_pre <- function(counts = make_top_var_counts(), method = "TMMlogCPM") {
    expr_work <- suppressMessages(normalize_counts(counts, method = "TMMlogCPM", prior.count = 1))
    list(
        expr_raw  = counts,
        expr_filt = counts,
        expr_work = expr_work,
        meta      = data.frame(SampleID = colnames(counts),
                               group = rep(c("ctrl", "trt"), each = 3),
                               stringsAsFactors = FALSE),
        row_data  = data.frame(gene_id = rownames(counts)),
        info      = list(expr_work_method = method)
    )
}

# Heatmaps and the 3D PCA are switched off: they are not under test and slow.
make_top_var_config <- function(min_mean_cpm = NULL) {
    qc <- list(thresholds = list(min_samples_for_pca3d = 100,
                                 max_samples_for_heatmaps = 1,
                                 max_samples_for_expr_heatmap = 1))
    if (!is.null(min_mean_cpm)) qc$top_variable_pca <- list(min_mean_cpm = min_mean_cpm)
    list(modes = list(rna = list(
        effects = list(samples = "SampleID", color = "group", shape = NULL),
        qc = qc
    )))
}

run_qc_pre <- function(pre, config, out_dir = tempfile("qc_pre_")) {
    res <- suppressMessages(mod_rnaseq_qc_pre(pre, config, out_dir))
    res$out_qc <- file.path(out_dir, "Diagnostic_plots")
    res
}

# The selection mod_rnaseq_qc_pre() made before the gate existed.
previous_top_idx <- function(mat, n_top) {
    gene_vars <- apply(mat, 1, var, na.rm = TRUE)
    order(gene_vars, decreasing = TRUE)[1:n_top]
}

top_genes_used <- function(p) rownames(attr(p, "pca_result")$rotation)

read_selection <- function(res) {
    utils::read.csv(file.path(res$out_qc, "PCA_top_variable_selection.csv"), stringsAsFactors = FALSE)
}

# ---- compute_tmm_cpm -----------------------------------------------------------

test_that("compute_tmm_cpm uses the effective library sizes of normalize_counts()", {
    counts <- make_top_var_counts()
    cpm <- compute_tmm_cpm(counts)
    # With a vanishing prior, 2^logCPM from normalize_counts() is linear CPM on
    # exactly the library sizes it used; compare where counts are positive.
    log_cpm <- suppressMessages(normalize_counts(counts, method = "TMMlogCPM", prior.count = 1e-9))
    pos <- counts > 0
    expect_equal(cpm[pos], 2^log_cpm[pos], tolerance = 1e-6)
    expect_identical(dimnames(cpm), dimnames(counts))

    # TMM scaling matters on this fixture, so raw-library CPM is a different answer.
    nf <- edgeR::calcNormFactors(edgeR::DGEList(counts), method = "TMM")$samples$norm.factors
    expect_gt(max(abs(log2(nf))), 0.05)
    expect_false(isTRUE(all.equal(cpm, compute_cpm(counts))))
})

# ---- rank_top_variable_genes ------------------------------------------------------

test_that("rank_top_variable_genes ranks only genes at or above the mean CPM", {
    mat <- matrix(c(0, 10,   # G1 var 50, low CPM
                    1, 7,    # G2 var 18
                    2, 3,    # G3 var 0.5
                    0, 4),   # G4 var 8
                  ncol = 2, byrow = TRUE,
                  dimnames = list(paste0("G", 1:4), c("S1", "S2")))
    mean_cpm <- c(0.4, 5, 1, 2)
    expect_identical(rank_top_variable_genes(mat, mean_cpm, 1), c(2L, 4L, 3L))
    expect_identical(rank_top_variable_genes(mat, mean_cpm, 3), 2L)
    expect_identical(rank_top_variable_genes(mat, c(NA, 5, 1, 2), 0), c(2L, 4L, 3L))
    expect_length(rank_top_variable_genes(mat, mean_cpm, 10), 0L)
})

test_that("rank_top_variable_genes with min_mean_cpm = 0 or no CPM is the previous ranking", {
    mat <- make_top_var_pre()$expr_work
    mean_cpm <- rowMeans(compute_tmm_cpm(make_top_var_counts()))
    old <- order(apply(mat, 1, var, na.rm = TRUE), decreasing = TRUE)
    expect_identical(rank_top_variable_genes(mat, mean_cpm, 0), old)
    expect_identical(rank_top_variable_genes(mat, NULL), old)
    expect_identical(rank_top_variable_genes(mat, NULL, 1), old)
})

test_that("rank_top_variable_genes keeps the previous order with NA, Inf and constant rows", {
    mat <- matrix(c(1, 5, 2,
                    NA, NA, NA,
                    3, 3, 3,
                    1, Inf, 2,
                    0, NA, 4,
                    2, 8, 1),
                  ncol = 3, byrow = TRUE,
                  dimnames = list(paste0("G", 1:6), paste0("S", 1:3)))
    old <- order(apply(mat, 1, var, na.rm = TRUE), decreasing = TRUE)
    expect_identical(rank_top_variable_genes(mat, rep(2, 6), 0), old)
    expect_identical(rank_top_variable_genes(mat, rep(2, 6), 1), old)
    expect_identical(rank_top_variable_genes(mat, c(2, 2, 2, 2, 0.5, 2), 1), setdiff(old, 5L))
})

test_that("rank_top_variable_genes rejects a mean_cpm that does not match the rows", {
    mat <- matrix(1:6, nrow = 3, dimnames = list(c("G1", "G2", "G3"), c("S1", "S2")))
    expect_error(rank_top_variable_genes(mat, c(1, 2), 1), "one per row")
})

# ---- mod_rnaseq_qc_pre --------------------------------------------------------------

test_that("the default gate ranks only genes with mean TMM CPM >= 1", {
    pre <- make_top_var_pre()
    mean_cpm <- rowMeans(compute_tmm_cpm(pre$expr_filt))
    n_pass <- sum(mean_cpm >= 1)
    # The fixture must put low-CPM genes among the most variable for this to test anything.
    expect_gt(sum(mean_cpm[previous_top_idx(pre$expr_work, 500)] < 1), 50)
    expect_gte(n_pass, 500)
    expect_lt(n_pass, 1000)

    res <- run_qc_pre(pre, make_top_var_config())
    used <- top_genes_used(res$plots$pca_top500)
    expected <- rownames(pre$expr_work)[rank_top_variable_genes(pre$expr_work, mean_cpm, 1)[1:500]]
    expect_setequal(used, expected)
    expect_true(all(mean_cpm[used] >= 1))

    sel <- read_selection(res)
    expect_identical(names(sel), c("min_mean_cpm", "n_candidates", "n_top", "n_genes_used", "status"))
    expect_equal(sel$n_top, c(500, 1000, 2000, 5000))
    expect_equal(sel$min_mean_cpm, rep(1, 4))
    expect_equal(sel$n_candidates, rep(n_pass, 4))
    expect_equal(sel$status, c("generated", rep("skipped_insufficient_genes", 3)))
    expect_equal(sel$n_genes_used, c(500L, NA, NA, NA))
})

test_that("a top-N plot is drawn only when N genes pass, and a stale one is removed", {
    pre <- make_top_var_pre()
    out_dir <- withr::local_tempdir()
    dir.create(file.path(out_dir, "Diagnostic_plots"), recursive = TRUE)
    stale <- file.path(out_dir, "Diagnostic_plots", "PCA_top1000.png")
    file.create(stale)

    res <- run_qc_pre(pre, make_top_var_config(), out_dir)
    expect_true(file.exists(file.path(res$out_qc, "PCA_top500.png")))
    expect_false(file.exists(stale))
    expect_named(res$plots[grep("^pca_top", names(res$plots))], "pca_top500")
    expect_false(any(grepl("PCA_top1000", res$files)))
})

test_that("min_mean_cpm = 0 reproduces the previous top-variable PCA", {
    pre <- make_top_var_pre()
    res <- run_qc_pre(pre, make_top_var_config(0))
    for (n_top in c(500, 1000)) {
        p <- res$plots[[sprintf("pca_top%d", n_top)]]
        old <- qc_pca_scatter(pre$expr_work[previous_top_idx(pre$expr_work, n_top), ], pre$meta,
                              make_top_var_config()$modes$rna, pcs = c(1, 2))
        expect_equal(attr(p, "scores"), attr(old, "scores"), info = n_top)
        expect_equal(attr(p, "var_expl"), attr(old, "var_expl"), info = n_top)
    }
    expect_true(file.exists(file.path(res$out_qc, "PCA_top1000.png")))
    expect_null(res$plots$pca_top2000)   # 2000 > number of genes, as before
    sel <- read_selection(res)
    expect_equal(sel$min_mean_cpm, rep(0, 4))
    expect_equal(sel$n_genes_used, c(500L, 1000L, NA, NA))
    expect_equal(sel$status, c("generated", "generated", "skipped_insufficient_genes",
                               "skipped_insufficient_genes"))
})

test_that("n_genes_used counts the genes the PCA kept, not the genes requested", {
    pre <- make_top_var_pre(method = "VST")
    # 450 varying genes and 600 constant ones: the top 500 and top 1000 must take
    # constant genes, which the PCA cannot use and drops.
    n_var <- 450
    pre$expr_work <- pre$expr_work[seq_len(n_var), ]
    pre$expr_work <- rbind(pre$expr_work,
                           matrix(1, 600, ncol(pre$expr_work),
                                  dimnames = list(paste0("C", 1:600), colnames(pre$expr_work))))
    res <- run_qc_pre(pre, make_top_var_config())
    sel <- read_selection(res)
    expect_equal(sel$status[1:2], c("generated", "generated"))
    expect_equal(sel$n_genes_used[1:2], c(n_var, n_var))
    expect_equal(nrow(attr(res$plots$pca_top500, "pca_result")$rotation), n_var)
})

test_that("a top-N PCA that errors is reported as failed with a warning", {
    pre <- make_top_var_pre()
    out_dir <- withr::local_tempdir()
    dir.create(file.path(out_dir, "Diagnostic_plots"), recursive = TRUE)
    stale <- file.path(out_dir, "Diagnostic_plots", "PCA_top1000.png")
    file.create(stale)

    # Fail only the top-1000 PCA; the stub lives where the module resolves it.
    orig <- get("qc_pca_scatter", envir = globalenv())
    assign("qc_pca_scatter", function(expr_mat, ...) {
        if (nrow(expr_mat) == 1000) stop("synthetic PCA failure")
        orig(expr_mat, ...)
    }, envir = globalenv())
    withr::defer(assign("qc_pca_scatter", orig, envir = globalenv()))

    expect_warning(res <- run_qc_pre(pre, make_top_var_config(0), out_dir),
                   "top 1000 variable genes was not generated: synthetic PCA failure")
    sel <- read_selection(res)
    expect_equal(sel$status[1:2], c("generated", "failed"))
    expect_equal(sel$n_genes_used[1:2], c(500L, NA))
    expect_null(res$plots$pca_top1000)
    expect_false(file.exists(stale))
})

test_that("VST or preprocessed expression is ranked without the gate", {
    for (method in c("VST", "preprocessed")) {
        pre <- make_top_var_pre(method = method)
        res <- run_qc_pre(pre, make_top_var_config())
        expect_setequal(top_genes_used(res$plots$pca_top500),
                        rownames(pre$expr_work)[previous_top_idx(pre$expr_work, 500)])
        sel <- read_selection(res)
        expect_true(all(is.na(sel$min_mean_cpm)), info = method)
        expect_equal(sel$n_candidates, rep(nrow(pre$expr_work), 4), info = method)
        expect_equal(sel$n_genes_used, c(500L, 1000L, NA, NA), info = method)
    }
})

test_that("the gate leaves the main PCA and the objects passed downstream unchanged", {
    pre <- make_top_var_pre()
    gated <- run_qc_pre(pre, make_top_var_config())
    ungated <- run_qc_pre(pre, make_top_var_config(0))
    for (key in c("objects", "norm_log_counts_pca", "pca_scores")) {
        expect_identical(gated[[key]], ungated[[key]], info = key)
    }
    main <- qc_pca_scatter(pre$expr_work, pre$meta, make_top_var_config()$modes$rna, pcs = c(1, 2))
    expect_equal(gated$pca_scores, attr(main, "scores"))
    expect_true(file.exists(file.path(gated$out_qc, "PCA_PC1.vs.PC2.png")))
})

test_that("expr_filt that does not line up with expr_work is an error, not a silent mismatch", {
    pre <- make_top_var_pre()
    pre$expr_filt <- pre$expr_filt[rev(rownames(pre$expr_filt)), ]
    expect_error(run_qc_pre(pre, make_top_var_config()), "different genes or samples")
})

# ---- report ---------------------------------------------------------------------------

rna_report_template <- function() {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "rnaseq", "report_template.Rmd"),
        "R/domain/rnaseq/report_template.Rmd"
    )
    candidates[file.exists(candidates)][1]
}

run_rule_chunk <- function(sidecar = NULL, top_pngs = "PCA_top500.png") {
    rmd <- rna_report_template()
    fence <- find_rmd_chunk_lines(rmd, "pca-var-genes-rule")
    code <- readLines(rmd, warn = FALSE)[(fence[1] + 1):(fence[2] - 1)]
    env <- new.env(parent = globalenv())
    env$diag_dir <- withr::local_tempdir(.local_envir = parent.frame())
    file.create(file.path(env$diag_dir, top_pngs))
    if (!is.null(sidecar)) {
        utils::write.csv(sidecar, file.path(env$diag_dir, "PCA_top_variable_selection.csv"), row.names = FALSE)
    }
    paste(utils::capture.output(eval(parse(text = code), envir = env)), collapse = "\n")
}

test_that("the report states the rule that ran, and nothing without the sidecar", {
    expect_identical(run_rule_chunk(), "")
    expect_match(run_rule_chunk(data.frame(min_mean_cpm = 1, n_candidates = 9812)),
                 "9,812 genes with a mean TMM-normalized CPM of at least 1", fixed = TRUE)
    expect_match(run_rule_chunk(data.frame(min_mean_cpm = NA, n_candidates = 300)),
                 "no mean-CPM restriction was applied", fixed = TRUE)
    expect_match(run_rule_chunk(data.frame(min_mean_cpm = 0, n_candidates = 300)),
                 "restriction disabled", fixed = TRUE)
    expect_false(grepl("Too few genes", run_rule_chunk(data.frame(min_mean_cpm = 1, n_candidates = 9812))))
})

test_that("the report says when no top-N plot qualified and the all-genes plot is shown", {
    out <- run_rule_chunk(data.frame(min_mean_cpm = 1, n_candidates = 320), top_pngs = character(0))
    expect_match(out, "320 genes with a mean TMM-normalized CPM", fixed = TRUE)
    expect_match(out, "the plot below uses all genes", fixed = TRUE)
})

test_that("the report reads the one-row-per-size sidecar once and names failed sizes", {
    sidecar <- data.frame(min_mean_cpm = 1, n_candidates = 9812,
                          n_top = c(500, 1000, 2000, 5000),
                          n_genes_used = c(500, NA, 2000, NA),
                          status = c("generated", "failed", "generated", "skipped_insufficient_genes"))
    out <- run_rule_chunk(sidecar)
    expect_length(gregexpr("Genes were ranked", out, fixed = TRUE)[[1]], 1L)
    expect_match(out, "9,812 genes with a mean TMM-normalized CPM of at least 1", fixed = TRUE)
    expect_match(out, "could not be generated (see the pipeline log): N = 1,000.", fixed = TRUE)

    sidecar$status <- c("generated", rep("skipped_insufficient_genes", 3))
    expect_false(grepl("could not be generated", run_rule_chunk(sidecar)))
})
