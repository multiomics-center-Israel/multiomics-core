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
