# Regression tests for the early-return guards in pathway / qc_post modules.
# Before this fix, isFALSE(NULL) returned FALSE so the guard didn't fire when
# the relevant config section was entirely absent, and the module silently
# ran with default settings (network annotation, file I/O, etc.). These
# tests pin the skip-when-absent behaviour so it can't regress.

test_that("mod_rnaseq_pathway: skips cleanly when modes.rna.pathway is absent", {
    config  <- list(modes = list(rna = list()))
    out_dir <- tempfile("test-rna-pathway-absent-")
    dir.create(out_dir)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    res <- suppressMessages(
        mod_rnaseq_pathway(de_res = NULL, pre = NULL,
                           config = config, out_dir = out_dir)
    )

    expect_type(res, "list")
    expect_named(res, c("annotation", "pathway_results", "plot_files"),
                 ignore.order = TRUE)
    expect_null(res$annotation)
    expect_length(res$pathway_results, 0)
    expect_length(res$plot_files, 0)

    # No work done: the Enrichment output dir should not have been created.
    expect_false(dir.exists(file.path(out_dir, "Enrichment")))
})

test_that("mod_rnaseq_pathway: skips when both pathway and enrichment are absent, even with DE data", {
    # Regression (#128): the enable gate must fire on a legacy/minimal config that
    # omits BOTH sections BEFORE any annotation/online work — not only when de_res
    # is NULL. Previously the guard leaned on de_res == NULL to short-circuit.
    config  <- list(modes = list(rna = list()))
    out_dir <- tempfile("test-rna-pathway-both-absent-")
    dir.create(out_dir)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    de_res <- list(tables = list(
        cA = data.frame(FeatureID = c("g1", "g2"),
                        log2FoldChange = c(2, -2),
                        pvalue = c(1e-4, 1e-3),
                        padj = c(1e-3, 1e-2),
                        stringsAsFactors = FALSE)
    ))

    res <- suppressMessages(
        mod_rnaseq_pathway(de_res = de_res, pre = NULL,
                           config = config, out_dir = out_dir)
    )
    expect_length(res$pathway_results, 0)
    expect_false(dir.exists(file.path(out_dir, "Enrichment")))
})

test_that("mod_proteomics_qc_post: skips cleanly when modes.proteomics.qc_post is absent", {
    config  <- list(modes = list(proteomics = list()))
    out_dir <- tempfile("test-prot-qcpost-absent-")
    dir.create(out_dir)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    ids <- paste0("P", 1:10)
    set.seed(42)
    summary_df <- data.frame(
        FeatureID = ids,
        `padj.imputs.A_vs_B`    = runif(10, 0, 0.1),
        `pvalue.imputs.A_vs_B`  = runif(10, 0, 0.05),
        `linearFC.imputs.A_vs_B` = runif(10, -3, 3),
        check.names = FALSE, stringsAsFactors = FALSE
    )

    expr_mat <- matrix(rnorm(40), nrow = 10, ncol = 4,
                       dimnames = list(ids, paste0("S", 1:4)))
    meta <- data.frame(sample_id = paste0("S", 1:4),
                       condition = rep(c("A", "B"), each = 2))

    de_res <- list(summary_df = summary_df)
    config <- list(modes = list(proteomics = list(
        de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5),
        de_table = list(id_col = "FeatureID"),
        qc = list(top_de_heatmap_sizes = c())
    )))
    pre <- list(expr_imp_single = expr_mat, meta = meta)

    res <- suppressWarnings(suppressMessages(
        mod_proteomics_qc_post(pre = pre, de_res = de_res,
                               config = config, out_dir = out_dir)
    ))

    expect_type(res, "list")
    expect_named(res, c("plots", "files"), ignore.order = TRUE)
    expect_true(length(res$files) > 0)
})
