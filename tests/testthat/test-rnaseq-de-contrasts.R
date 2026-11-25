# mod_rnaseq_de() returns the contrast table it used (de_res$contrasts), so
# the omics bundle records Factor / Numerator / Denominator even when the
# contrasts were auto-generated. The DE results themselves must not change.
#
# Uses the limma route (make_bc_pre() is a "preprocessed" pre object) from
# helper-rna-batch.R. All data are synthetic.

contrast_de_config <- function() bc_config(enabled = FALSE)

# The DE result exactly as the route returned it before $contrasts existed.
direct_limma_de <- function(pre, contrasts_df, config) {
    suppressMessages(run_limma_rna_de(
        expr = pre$expr_work, meta = pre$meta, contrasts_df = contrasts_df,
        de_cfg = config$modes$rna$de,
        batch_cfg = get_rna_batch_model_config(config$modes$rna)))
}


test_that("supplied contrasts are returned unchanged and DE results are unchanged", {
    skip_if_not_installed("limma")
    pre <- make_bc_pre()
    cfg <- contrast_de_config()

    res <- suppressMessages(mod_rnaseq_de(pre, list(contrasts = pre$contrasts), cfg))

    expect_identical(res$contrasts, pre$contrasts)
    expect_named(res$tables, pre$contrasts$Contrast_name)
    # Everything else is what the DE route returned before this change.
    expect_identical(res[setdiff(names(res), "contrasts")],
                     direct_limma_de(pre, pre$contrasts, cfg))
})

test_that("auto-generated contrasts are returned as the DE step resolved them", {
    skip_if_not_installed("limma")
    pre <- make_bc_pre()
    pre["contrasts"] <- list(NULL)   # as preprocess_rna() leaves it with no contrasts file
    cfg <- contrast_de_config()

    res <- suppressMessages(mod_rnaseq_de(pre, list(contrasts = NULL), cfg))

    # The same resolver mod_rnaseq_de() uses, on the group column it reads.
    expected <- suppressMessages(auto_generate_contrasts(pre$meta, "Condition"))
    expect_identical(res$contrasts, expected)
    expect_true(all(c("Contrast_name", "Factor", "Numerator", "Denominator")
                    %in% names(res$contrasts)))
    expect_named(res$tables, expected$Contrast_name)
    expect_identical(res[setdiff(names(res), "contrasts")],
                     direct_limma_de(pre, expected, cfg))
})

# Pre-computed route: one small DE table per file, in a temp raw directory.
precomputed_de_config <- function(files) {
    dir <- withr::local_tempdir(.local_envir = parent.frame())
    for (f in files) {
        utils::write.table(
            data.frame(FeatureID = c("g1", "g2"), log2FoldChange = c(1, -1),
                       pvalue = c(0.01, 0.5), padj = c(0.02, 0.5)),
            file.path(dir, f), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    list(project = list(dir = dir), paths = list(raw = "."),
         modes = list(rna = list(files = list(de_table = as.list(files)))))
}

test_that("pre-computed DE keeps supplied contrasts that labelled its tables", {
    cfg <- precomputed_de_config(c("de_first.tsv", "de_second.tsv"))
    supplied <- data.frame(Contrast_name = c("B_vs_A", "C_vs_A"), Factor = "Group",
                           Numerator = c("B", "C"), Denominator = "A",
                           stringsAsFactors = FALSE)

    res <- suppressMessages(load_precomputed_rna_de(cfg, contrasts_df = supplied))

    expect_named(res$tables, supplied$Contrast_name)
    expect_identical(res$contrasts, supplied)
})

test_that("pre-computed DE records no contrasts when labels came from file names", {
    # Two files, one supplied contrast: the loader cannot pair them and names
    # the tables after the files, so the supplied table does not describe them.
    cfg <- precomputed_de_config(c("de_first.tsv", "de_second.tsv"))
    supplied <- data.frame(Contrast_name = "B_vs_A", Factor = "Group",
                           Numerator = "B", Denominator = "A",
                           stringsAsFactors = FALSE)

    res <- suppressMessages(load_precomputed_rna_de(cfg, contrasts_df = supplied))

    expect_named(res$tables, c("first", "second"))
    expect_null(res$contrasts)
})

test_that("the RNA bundle carries the resolved contrasts", {
    skip_if_not_installed("limma")
    pre <- make_bc_pre()
    pre["contrasts"] <- list(NULL)
    cfg <- contrast_de_config()
    res <- suppressMessages(mod_rnaseq_de(pre, list(contrasts = NULL), cfg))

    parts <- build_rnaseq_bundle_parts(pre = pre, de_res = res, config = cfg)
    expect_identical(parts$de$contrasts, res$contrasts)
    expect_identical(parts$de$contrasts$Numerator, res$contrasts$Numerator)
})
