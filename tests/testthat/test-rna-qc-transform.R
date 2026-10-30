# tests/testthat/test-rna-qc-transform.R
#
# RNA-seq QC plots use their own matrix. With TMMlogCPM it is
# log2(TMM CPM + qc$cpm_pseudocount): edgeR's prior.count is in reads and, at
# typical depths, a far smaller offset than 1 CPM, so genes with a handful of
# reads dominate the sample-level structure of expr_work. The QC matrix feeds
# the PCAs (PC1 vs PC2, which Shiny receives, PC1 vs PC3, 3D, top-variable), the
# sample distance/correlation heatmaps and the pre-DE expression heatmaps.
# expr_work, DE, clustering and the exported tables are unchanged, and VST and
# preprocessed input are plotted as before. All data here are synthetic.

# 700 expressed genes (half shifted between groups) and 500 genes with a few
# reads. Libraries are ~1M reads, so one read is ~1 CPM. Sample S4 carries a
# block of very high counts so that TMM factors move away from 1.
make_qc_counts <- function() {
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

make_qc_pre <- function(counts = make_qc_counts(), method = "TMMlogCPM") {
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

# Heatmaps and the 3D PCA are off unless a test switches them on.
make_qc_config <- function(cpm_pseudocount = NULL, heatmaps = FALSE) {
    qc <- list(thresholds = list(min_samples_for_pca3d = if (heatmaps) 1 else 100,
                                 max_samples_for_heatmaps = if (heatmaps) 100 else 1,
                                 max_samples_for_expr_heatmap = if (heatmaps) 100 else 1))
    qc$cpm_pseudocount <- cpm_pseudocount
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

# The QC matrix written out by hand, from edgeR directly.
expected_qc_matrix <- function(counts, pseudocount = 1) {
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts), method = "TMM")
    log2(edgeR::cpm(dge, log = FALSE) + pseudocount)
}

scatter_on <- function(mat, pre, pcs = c(1, 2)) {
    qc_pca_scatter(mat, pre$meta, make_qc_config()$modes$rna, pcs = pcs)
}

top_genes_used <- function(p) rownames(attr(p, "pca_result")$rotation)

# Replace a function where the module resolves it (the sourced functions live
# in globalenv), restored when the calling test ends.
local_stub <- function(name, fun, env = parent.frame()) {
    orig <- get(name, envir = globalenv())
    assign(name, fun, envir = globalenv())
    withr::defer(assign(name, orig, envir = globalenv()), envir = env)
}

# ---- compute_tmm_cpm -----------------------------------------------------------

test_that("compute_tmm_cpm uses the effective library sizes of normalize_counts()", {
    counts <- make_qc_counts()
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

# ---- resolve_rna_qc_transform / rna_qc_expr -----------------------------------------

test_that("the QC transform applies to TMMlogCPM only, with a default pseudocount of 1", {
    tmm <- list(expr_work_method = "TMMlogCPM")
    expect_identical(resolve_rna_qc_transform(tmm)$label, "log2(TMM CPM + 1)")
    expect_identical(resolve_rna_qc_transform(tmm, list(cpm_pseudocount = 0.5))$pseudocount, 0.5)
    expect_identical(resolve_rna_qc_transform(tmm, list(cpm_pseudocount = 0.5))$label, "log2(TMM CPM + 0.5)")
    for (method in list("VST", "preprocessed", NULL)) {
        tr <- resolve_rna_qc_transform(list(expr_work_method = method), list(cpm_pseudocount = 1))
        expect_false(tr$applied)
        expect_true(is.na(tr$pseudocount))
    }
})

test_that("the QC matrix is log2(TMM CPM + pseudocount) and expr_work is untouched", {
    pre <- make_qc_pre()
    before <- pre
    qc <- rna_qc_expr(pre)
    expect_equal(qc, expected_qc_matrix(pre$expr_filt, 1))
    expect_identical(dimnames(qc), dimnames(pre$expr_work))
    expect_false(isTRUE(all.equal(qc, pre$expr_work, check.attributes = FALSE)))
    expect_equal(rna_qc_expr(pre, list(cpm_pseudocount = 0.5)), expected_qc_matrix(pre$expr_filt, 0.5))
    expect_identical(pre, before)
})

test_that("VST and preprocessed expression are used for QC as they are", {
    for (method in c("VST", "preprocessed")) {
        pre <- make_qc_pre(method = method)
        expect_identical(rna_qc_expr(pre, list(cpm_pseudocount = 1)), pre$expr_work, info = method)
    }
})

test_that("expr_filt that does not line up with expr_work is an error, not a silent mismatch", {
    pre <- make_qc_pre()
    pre$expr_filt <- pre$expr_filt[rev(rownames(pre$expr_filt)), ]
    expect_error(rna_qc_expr(pre), "different genes or samples")
})

# ---- mod_rnaseq_qc_pre ------------------------------------------------------------

test_that("the main and top-variable PCAs use the QC matrix, ranked by plain variance", {
    pre <- make_qc_pre()
    qc <- expected_qc_matrix(pre$expr_filt)
    res <- run_qc_pre(pre, make_qc_config())

    # PC1 vs PC2 is also what Shiny receives.
    main <- scatter_on(qc, pre)
    expect_equal(res$pca_scores, attr(main, "scores"))
    expect_equal(res$norm_log_counts_pca$x, attr(main, "pca_result")$x)
    expect_equal(res$objects$var_expl, attr(main, "var_expl"))
    expect_equal(attr(res$plots$pca_1_3, "scores"), attr(scatter_on(qc, pre, c(1, 3)), "scores"))

    # No CPM gate: the top 500 are the 500 most variable rows of the QC matrix.
    top500 <- rownames(qc)[order(apply(qc, 1, var), decreasing = TRUE)[1:500]]
    expect_setequal(top_genes_used(res$plots$pca_top500), top500)
    expect_equal(attr(res$plots$pca_top500, "scores"), attr(scatter_on(qc[top500, ], pre), "scores"))
    expect_false(file.exists(file.path(res$out_qc, "PCA_top_variable_selection.csv")))
})

test_that("the PCAs change with the transform, but the Shiny PCA objects keep their shape", {
    pre <- make_qc_pre()
    res <- run_qc_pre(pre, make_qc_config())
    old <- scatter_on(pre$expr_work, pre)
    expect_false(isTRUE(all.equal(res$pca_scores$PC1, attr(old, "scores")$PC1)))
    expect_identical(names(res$pca_scores), names(attr(old, "scores")))
    expect_identical(res$pca_scores$sample, attr(old, "scores")$sample)
    expect_s3_class(res$norm_log_counts_pca, "prcomp")
    expect_identical(res$objects$pca_scores, res$pca_scores)
    expect_identical(res$objects$norm_log_counts_pca, res$norm_log_counts_pca)
})

test_that("heatmaps and the 3D PCA use the QC matrix; distribution plots use expr_work", {
    pre <- make_qc_pre()
    qc <- expected_qc_matrix(pre$expr_filt)
    calls <- new.env()
    calls$log <- list()
    recorder <- function(name) {
        force(name)
        function(expr_mat, ...) {
            calls$log[[length(calls$log) + 1]] <- list(fn = name, expr = as.matrix(expr_mat),
                                                       title = list(...)$title)
            invisible(NULL)
        }
    }
    for (fn in c("qc_pca_3d", "qc_sample_distance_heatmap", "qc_sample_correlation_heatmap",
                 "wrap_qc_heatmap", "plot_heatmap_core", "qc_omic_density", "norm_boxplot")) {
        local_stub(fn, recorder(fn))
    }
    local_stub("save_heatmap_to_file", function(...) invisible(NULL))

    run_qc_pre(pre, make_qc_config(heatmaps = TRUE))
    fns <- vapply(calls$log, `[[`, "", "fn")
    same_values <- function(a, b) {
        identical(dimnames(a), dimnames(b)) && isTRUE(all.equal(c(a), c(b)))
    }

    for (fn in c("qc_pca_3d", "qc_sample_distance_heatmap", "qc_sample_correlation_heatmap",
                 "wrap_qc_heatmap", "plot_heatmap_core")) {
        hits <- calls$log[fns == fn]
        expect_gte(length(hits), 1L)
        for (h in hits) expect_true(same_values(h$expr, qc), info = fn)
    }
    expect_length(calls$log[fns == "plot_heatmap_core"], 2L)

    dens <- calls$log[fns == "qc_omic_density"][[1]]
    expect_true(same_values(dens$expr, pre$expr_work))
    expect_match(dens$title, "analysis matrix, TMMlogCPM", fixed = TRUE)
    boxes <- calls$log[fns == "norm_boxplot"]
    norm_box <- boxes[vapply(boxes, function(b) grepl("analysis matrix", b$title), logical(1))]
    expect_length(norm_box, 1L)
    expect_true(same_values(norm_box[[1]]$expr, pre$expr_work))
})

test_that("VST and preprocessed input keep the previous QC plots", {
    for (method in c("VST", "preprocessed")) {
        pre <- make_qc_pre(method = method)
        res <- run_qc_pre(pre, make_qc_config(cpm_pseudocount = 1))
        expect_equal(res$pca_scores, attr(scatter_on(pre$expr_work, pre), "scores"), info = method)
        old_top <- order(apply(pre$expr_work, 1, var, na.rm = TRUE), decreasing = TRUE)[1:500]
        expect_setequal(top_genes_used(res$plots$pca_top500), rownames(pre$expr_work)[old_top])
    }
})

test_that("PCA file names are unchanged and sizes beyond the gene count are skipped", {
    pre <- make_qc_pre()
    res <- run_qc_pre(pre, make_qc_config())
    for (f in c("PCA_PC1.vs.PC2.png", "PCA_PC1.vs.PC3.png", "PCA_top500.png", "PCA_top1000.png")) {
        expect_true(file.exists(file.path(res$out_qc, f)), info = f)
    }
    expect_named(res$plots[grep("^pca_top", names(res$plots))], c("pca_top500", "pca_top1000"))
    expect_false(file.exists(file.path(res$out_qc, "PCA_top2000.png")))
})

test_that("a stale top-N plot is removed when this run does not draw that size", {
    pre <- make_qc_pre()
    out_dir <- withr::local_tempdir()
    dir.create(file.path(out_dir, "Diagnostic_plots"), recursive = TRUE)
    stale <- file.path(out_dir, "Diagnostic_plots", "PCA_top2000.png")
    file.create(stale)
    run_qc_pre(pre, make_qc_config(), out_dir)
    expect_false(file.exists(stale))
})

test_that("a top-N PCA that errors gives a warning and leaves no plot behind", {
    pre <- make_qc_pre()
    out_dir <- withr::local_tempdir()
    dir.create(file.path(out_dir, "Diagnostic_plots"), recursive = TRUE)
    stale <- file.path(out_dir, "Diagnostic_plots", "PCA_top1000.png")
    file.create(stale)

    orig <- get("qc_pca_scatter", envir = globalenv())
    local_stub("qc_pca_scatter", function(expr_mat, ...) {
        if (nrow(expr_mat) == 1000) stop("synthetic PCA failure")
        orig(expr_mat, ...)
    })

    expect_warning(res <- run_qc_pre(pre, make_qc_config(), out_dir),
                   "top 1000 variable genes was not generated: synthetic PCA failure")
    expect_null(res$plots$pca_top1000)
    expect_false(is.null(res$plots$pca_top500))
    expect_false(file.exists(stale))
})

# ---- fact sheet -----------------------------------------------------------------------

test_that("the fact sheet quotes the QC PCA and QC-matrix correlations", {
    pre <- make_qc_pre()
    config <- make_qc_config()
    res <- run_qc_pre(pre, config)
    qc <- expected_qc_matrix(pre$expr_filt)
    claim <- function(rows, what) Filter(function(r) !is.null(r) && r$claim == what, rows)[[1]]

    rows <- .facts_rnaseq_structure(NULL, NA_character_, pre, config, res)
    pc <- claim(rows, "variance explained by PC1 and PC2")
    expect_identical(pc$value, sprintf("%.2f%% and %.2f%%", 100 * res$objects$var_expl[1],
                                       100 * res$objects$var_expl[2]))
    expect_match(pc$source_file, "PCA_PC1.vs.PC2.png (QC matrix, log2(TMM CPM + 1))", fixed = TRUE)

    # Without the QC result it recomputes with the QC PCA's genes and settings:
    # the same numbers, not just the same input matrix.
    rows_re <- .facts_rnaseq_structure(NULL, NA_character_, pre, config)
    expect_identical(claim(rows_re, "variance explained by PC1 and PC2")$value, pc$value)

    cm <- stats::cor(qc)
    expect_identical(claim(rows, "pearson correlation between libraries")$value,
                     fmt_range(cm[upper.tri(cm)], 3))
    expect_match(claim(rows, "pearson correlation between libraries")$source_file,
                 "log2(TMM CPM + 1)", fixed = TRUE)
})

# ---- report and methods text ------------------------------------------------------------

rna_report_template <- function() {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "rnaseq", "report_template.Rmd"),
        "R/domain/rnaseq/report_template.Rmd"
    )
    candidates[file.exists(candidates)][1]
}

# The lines of the report's setup chunk that define rna_qc_sentence.
qc_sentence_for <- function(rna_info, rna_cfg = list(normalization = list(prior.count = 1))) {
    src <- readLines(rna_report_template(), warn = FALSE)
    start <- grep("^rna_qc_tr <- ", src)
    stopifnot(length(start) == 1)
    sentence_start <- grep("^rna_qc_sentence <- ", src)
    end <- sentence_start + which(src[(sentence_start + 1):length(src)] == "}")[1]
    env <- new.env(parent = globalenv())
    env$params <- list(rna_info = rna_info)
    env$rna_cfg <- rna_cfg
    eval(parse(text = src[start:end]), envir = env)
    env$rna_qc_sentence
}

test_that("the report names the QC transform apart from the analysis matrix", {
    tmm <- qc_sentence_for(list(expr_work_method = "TMMlogCPM"))
    expect_match(tmm, "QC-only transformation, log2(TMM CPM + 1)", fixed = TRUE)
    expect_match(tmm, "edgeR log2 CPM with prior.count = 1", fixed = TRUE)
    custom <- qc_sentence_for(list(expr_work_method = "TMMlogCPM"),
                              list(normalization = list(prior.count = 2),
                                   qc = list(cpm_pseudocount = 0.5)))
    expect_match(custom, "log2(TMM CPM + 0.5)", fixed = TRUE)
    expect_match(custom, "prior.count = 2", fixed = TRUE)
    expect_match(qc_sentence_for(list(expr_work_method = "VST")),
                 "use the VST values of the analysis matrix, the same matrix", fixed = TRUE)
    expect_match(qc_sentence_for(NULL), "not recorded", fixed = TRUE)
})

test_that("the pipeline summary's normalization step names the QC transform", {
    config <- list(project = list(name = "P"),
                   modes = list(rna = list(normalization = list(method = "TMMlogCPM", prior.count = 1))))
    pre <- list(info = list(expr_work_method = "TMMlogCPM"),
                expr_raw = matrix(1, 5, 2), expr_filt = matrix(1, 4, 2),
                meta = data.frame(Condition = c("A", "B")))
    body <- generate_summary_body_r(collect_pipeline_stats(config, pre, de_res = list()))
    expect_match(body, "QC plots only (PCA, sample and expression heatmaps) use log2(TMM CPM + 1)", fixed = TRUE)
    expect_match(body, "prior count = 1", fixed = TRUE)

    # VST configured but fell back to TMMlogCPM: describe what ran.
    config$modes$rna$normalization$method <- "VST"
    body <- generate_summary_body_r(collect_pipeline_stats(config, pre, de_res = list()))
    expect_match(body, "log2 CPM (prior count = 1)", fixed = TRUE)
    expect_no_match(body, "Variance Stabilizing", fixed = TRUE)
    expect_match(body, "QC plots only", fixed = TRUE)

    pre$info$expr_work_method <- "VST"
    body <- generate_summary_body_r(collect_pipeline_stats(config, pre, de_res = list()))
    expect_no_match(body, "QC plots only", fixed = TRUE)
})

test_that("the report describes the normalization that ran, not a failed VST", {
    src <- readLines(rna_report_template(), warn = FALSE)
    lines <- src[grep("^rna_norm_ran <- ", src):grep("^rna_norm_method <- ", src)]
    method_for <- function(configured, ran) {
        env <- new.env(parent = globalenv())
        env$params <- list(rna_info = if (is.null(ran)) NULL else list(expr_work_method = ran))
        env$rna_cfg <- list(normalization = list(method = configured))
        eval(parse(text = lines), envir = env)
        env$rna_norm_method
    }
    expect_identical(method_for("VST", "TMMlogCPM"), "TMMlogCPM")
    expect_identical(method_for("VST", "VST"), "VST")
    expect_identical(method_for("TMMlogCPM", "preprocessed"), "TMMlogCPM")
    expect_identical(method_for("VST", NULL), "VST")
    expect_identical(method_for(NULL, NULL), "TMMlogCPM")
})
