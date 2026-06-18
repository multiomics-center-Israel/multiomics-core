# Regression tests for arg-drift between pathway modules and run_pathway_analysis().
# Each module runs end-to-end with a tiny synthetic DE table + an inline GMT,
# and asserts no error — catches "unused argument" and similar wiring breakage
# that CI's unit tests don't otherwise exercise (CI doesn't invoke tar_make()).

test_that("mod_rnaseq_pathway: end-to-end wiring with synthetic DE + tiny GMT", {
    skip_if_not_installed("fgsea")

    gmt <- tempfile(fileext = ".gmt")
    writeLines(c(
        "SET_A\tdesc A\tg1\tg2\tg3\tg4\tg5",
        "SET_B\tdesc B\tg6\tg7\tg8\tg9\tg10\tg11\tg12"
    ), gmt)
    on.exit(unlink(gmt), add = TRUE)

    de_res <- list(tables = list(
        contrast_A = data.frame(
            FeatureID      = paste0("g", 1:15),
            log2FoldChange = c(2.0, -1.5, 0.1, 1.2, -2.1, 0.5, 0.3, -0.7, 1.8, -1.1,
                               0.4, -0.2, 1.5, -1.8, 0.6),
            pvalue         = c(0.001, 0.01, 0.8, 0.04, 0.005, 0.2, 0.3, 0.15, 0.002, 0.03,
                               0.5, 0.7, 0.008, 0.001, 0.4),
            padj           = c(0.01, 0.05, 0.9, 0.1, 0.02, 0.3, 0.4, 0.25, 0.01, 0.08,
                               0.6, 0.8, 0.04, 0.01, 0.5),
            stat           = c(5, -3, 0.2, 2.5, -4, 1, 1.2, -1.5, 4.5, -2.8,
                               0.8, -0.3, 3.2, -4.5, 1.1),
            stringsAsFactors = FALSE
        )
    ))
    pre    <- list(expr_filt = NULL, meta = NULL)
    config <- list(modes = list(rna = list(
        annotation = list(skip_annotation = TRUE, organism = "Homo sapiens"),
        pathway    = list(
            enabled   = TRUE,
            databases = character(0),
            gmt_file  = gmt,
            min_size  = 3,
            max_size  = 50
        )
    )))
    out_dir <- tempfile("test-rna-pathway-")
    dir.create(out_dir)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    expect_no_error(suppressMessages(suppressWarnings(
        mod_rnaseq_pathway(de_res, pre, config, out_dir)
    )))
})

test_that("run_proteomics_pathway: end-to-end wiring with synthetic summary + tiny GMT", {
    skip_if_not_installed("fgsea")

    gmt <- tempfile(fileext = ".gmt")
    writeLines(c(
        "SET_A\tdesc A\tg1\tg2\tg3\tg4\tg5",
        "SET_B\tdesc B\tg6\tg7\tg8\tg9\tg10\tg11\tg12"
    ), gmt)
    on.exit(unlink(gmt), add = TRUE)

    summary_df <- data.frame(
        FeatureID                  = paste0("p", 1:15),
        Genes                      = paste0("g", 1:15),
        padj.imputs.contrast_A     = c(0.01, 0.05, 0.9, 0.1, 0.02, 0.3, 0.4, 0.25, 0.01, 0.08,
                                       0.6, 0.8, 0.04, 0.01, 0.5),
        pvalue.imputs.contrast_A   = c(0.001, 0.01, 0.8, 0.04, 0.005, 0.2, 0.3, 0.15, 0.002, 0.03,
                                       0.5, 0.7, 0.008, 0.001, 0.4),
        linearFC.imputs.contrast_A = c(4.0, -2.8, 1.1, 2.3, -4.2, 1.4, 1.2, -1.6, 3.5, -2.1,
                                       1.3, -1.1, 2.8, -3.5, 1.5),
        stringsAsFactors = FALSE,
        check.names      = FALSE
    )

    de_res <- list(summary_df = summary_df)
    pre    <- list()
    config <- list(modes = list(proteomics = list(
        de_table   = list(id_col = "FeatureID"),
        annotation = list(skip_annotation = TRUE, organism = "Homo sapiens"),
        pathway    = list(
            enabled   = TRUE,
            databases = character(0),
            gmt_file  = gmt,
            min_size  = 3,
            max_size  = 50
        )
    )))
    out_dir <- tempfile("test-prot-pathway-")
    dir.create(out_dir)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    expect_no_error(suppressMessages(suppressWarnings(
        run_proteomics_pathway(de_res, pre, config, out_dir)
    )))
})
