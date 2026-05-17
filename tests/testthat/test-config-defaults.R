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

test_that("mod_proteomics_qc_post: skips cleanly when modes.proteomics.qc_post is absent", {
    config  <- list(modes = list(proteomics = list()))
    out_dir <- tempfile("test-prot-qcpost-absent-")
    dir.create(out_dir)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    res <- suppressMessages(
        mod_proteomics_qc_post(pre = NULL, de_res = NULL,
                               config = config, out_dir = out_dir)
    )

    expect_type(res, "list")
    expect_named(res, c("plots", "files"), ignore.order = TRUE)
    expect_length(res$plots, 0)
    expect_length(res$files, 0)

    # No work done: create_legacy_output_dirs() should not have been called.
    expect_false(dir.exists(file.path(out_dir, "Diagnostic_plots")))
})
