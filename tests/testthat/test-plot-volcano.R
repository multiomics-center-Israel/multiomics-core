# Tests for plot_volcano p-value-type selection.
#
# Regression tests for the bug where `pvalue_type` was swallowed by `...`, so the
# "padj" and "pval" volcano variants were byte-identical (both drawn from padj)
# and only the title differed. plot_volcano must now drive the y-axis, dashed
# line, significance coloring, and subtitle from the *selected* p-value column.

# Minimal config: a 0.05 cutoff and a 1.5 linear-FC cutoff.
volcano_cfg <- list(de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5))

# A synthetic DE table where padj and raw P.Value diverge enough that the two
# variants must look different. Feature 3 is significant by raw p but not padj.
make_de_tbl <- function() {
    data.frame(
        FeatureID = paste0("F", 1:4),
        logFC     = c(2.0, -2.0, 1.0, 0.1),
        P.Value   = c(1e-6, 1e-5, 1e-3, 0.9),
        adj.P.Val = c(1e-4, 1e-3, 0.20, 0.95),
        stringsAsFactors = FALSE
    )
}

# Pull the y aesthetic (-log10 p) actually mapped into the built plot.
volcano_y <- function(p) {
    sort(ggplot2::ggplot_build(p)$data[[1]]$y)
}

test_that("padj and pval variants are not identical", {
    de <- make_de_tbl()
    p_padj <- plot_volcano(de, volcano_cfg, title = "x", pvalue_type = "padj")
    p_pval <- plot_volcano(de, volcano_cfg, title = "x", pvalue_type = "pval")

    # The y-axis data must differ, because they transform different columns.
    expect_false(isTRUE(all.equal(volcano_y(p_padj), volcano_y(p_pval))))

    # Subtitle and y-axis label must name the column actually used.
    expect_match(p_padj$labels$subtitle, "adj\\.P\\.Val")
    expect_match(p_padj$labels$y, "adj\\.P\\.Val")
    expect_match(p_pval$labels$subtitle, "P\\.Value")
    expect_match(p_pval$labels$y, "P\\.Value")
})

test_that("default pvalue_type is padj", {
    de <- make_de_tbl()
    expect_equal(
        volcano_y(plot_volcano(de, volcano_cfg)),
        volcano_y(plot_volcano(de, volcano_cfg, pvalue_type = "padj"))
    )
})

test_that("invalid pvalue_type is rejected by match.arg", {
    de <- make_de_tbl()
    expect_error(plot_volcano(de, volcano_cfg, pvalue_type = "raw"))
})

test_that("pvalue_type='pval' errors when no raw p-value column exists", {
    # padj-only table (e.g. a DESeq2-style result with no raw p column).
    de <- data.frame(
        FeatureID = paste0("F", 1:3),
        logFC     = c(1.0, -1.0, 0.2),
        padj      = c(0.01, 0.02, 0.9),
        stringsAsFactors = FALSE
    )
    # padj still works...
    expect_silent(plot_volcano(de, volcano_cfg, pvalue_type = "padj"))
    # ...but pval must fail loudly rather than fall back to padj.
    expect_error(
        plot_volcano(de, volcano_cfg, pvalue_type = "pval"),
        "raw p-value"
    )
})

test_that("qc_post caller pattern skips only the pval variant on all-NA raw p", {
    # The volcano emit loops (rnaseq/proteomics/metabolomics qc_post) wrap each
    # ptype in a tryCatch so a missing/all-NA raw p-value skips ONLY the pval
    # plot (warning, no error), leaving padj + downstream targets intact. This
    # replicates that wrapper to lock the contract.
    de <- data.frame(
        FeatureID = paste0("F", 1:3),
        logFC     = c(1.0, -1.0, 0.2),
        P.Value   = NA_real_,
        adj.P.Val = c(0.01, 0.02, 0.9),
        stringsAsFactors = FALSE
    )

    emit_one <- function(ptype) {
        tryCatch(
            plot_volcano(de, volcano_cfg, title = "x", pvalue_type = ptype),
            error = function(e) {
                if (ptype == "pval") {
                    warning("Skipping pval volcano: ", conditionMessage(e))
                    return(NULL)
                }
                stop(e)
            }
        )
    }

    # padj must succeed and yield a real plot (a thrown error fails the test).
    p_padj <- emit_one("padj")
    expect_s3_class(p_padj, "ggplot")

    # pval must warn-and-skip (NULL), not abort the loop.
    expect_warning(p_pval <- emit_one("pval"), "Skipping pval volcano")
    expect_null(p_pval)
})

test_that("pvalue_type='pval' errors on an all-NA P.Value column (RNA fallback)", {
    # Mirrors get_rna_de_tables_qc_post leaving P.Value as NA for contrasts that
    # have no genuine raw p-value column, instead of silently copying padj.
    de <- data.frame(
        FeatureID = paste0("F", 1:3),
        logFC     = c(1.0, -1.0, 0.2),
        P.Value   = NA_real_,
        adj.P.Val = c(0.01, 0.02, 0.9),
        stringsAsFactors = FALSE
    )
    # padj is intact, so the padj volcano is fine.
    expect_silent(plot_volcano(de, volcano_cfg, pvalue_type = "padj"))
    # The pval volcano must not silently become a copy of the padj plot.
    expect_error(
        plot_volcano(de, volcano_cfg, pvalue_type = "pval"),
        "entirely NA"
    )
})
