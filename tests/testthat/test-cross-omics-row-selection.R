# Which rows a cross-omics figure shows, and what a blank cell means.
#
# The figures exist to show where the layers agree. Ordering by combined p-value
# alone does not do that: the layer with the largest gene-set collection
# contributes many single-layer pathways with very small p-values and fills every
# slot, so the pathways several layers support never appear.
#
# Selection is for display only -- the meta-analysis table keeps its own
# ordering, and no p-value or adjustment is touched.
#
# All fixtures are synthetic.

meta_fixture <- function(norm_id, n_omics, combined_pval,
                         pval_a = NA_real_, pval_b = NA_real_) {
    data.frame(
        pathway = paste("Pathway", norm_id),
        norm_id = norm_id,
        pval_transcriptomics = pval_a,
        pval_proteomics = pval_b,
        combined_pval = combined_pval,
        n_omics = n_omics,
        stringsAsFactors = FALSE
    )
}


# ---- what leads the figure -------------------------------------------------

test_that("a two-layer pathway outranks a stronger single-layer one", {
    # The single-layer row has a combined p-value four orders of magnitude
    # smaller, which is exactly the situation that filled every slot.
    meta <- meta_fixture(
        norm_id      = c("00010", "00020"),
        n_omics      = c(1L, 2L),
        combined_pval = c(1e-9, 1e-5)
    )

    expect_identical(select_multi_omics_pathways(meta, top_n = 2)$norm_id,
                     c("00020", "00010"))
    # And with room for one row, the supported pathway is the one shown.
    expect_identical(select_multi_omics_pathways(meta, top_n = 1)$norm_id, "00020")
})

test_that("within the same support level, the combined p-value decides", {
    meta <- meta_fixture(
        norm_id       = c("00030", "00010", "00020"),
        n_omics       = c(2L, 2L, 2L),
        combined_pval = c(1e-3, 1e-7, 1e-5)
    )

    expect_identical(select_multi_omics_pathways(meta, top_n = 3)$norm_id,
                     c("00010", "00020", "00030"))
})

test_that("rows that tie on both keys keep the order they arrived in", {
    # The incoming table is already sorted by combined p-value, so falling back
    # to it makes the selection reproducible rather than dependent on sort
    # stability.
    meta <- meta_fixture(
        norm_id       = c("00010", "00020", "00030"),
        n_omics       = c(2L, 2L, 2L),
        combined_pval = c(1e-5, 1e-5, 1e-5)
    )

    expect_identical(select_multi_omics_pathways(meta, top_n = 3)$norm_id,
                     c("00010", "00020", "00030"))
})

test_that("selection copes with a table missing the columns it ranks on", {
    bare <- data.frame(norm_id = c("00010", "00020"), stringsAsFactors = FALSE)

    # Nothing to rank on is not a reason to fail or to invent a preference.
    expect_identical(select_multi_omics_pathways(bare, top_n = 2)$norm_id,
                     c("00010", "00020"))
    expect_equal(nrow(select_multi_omics_pathways(bare, top_n = 1)), 1L)
    expect_equal(nrow(select_multi_omics_pathways(
        bare[0, , drop = FALSE], top_n = 5)), 0L)
})

test_that("selection reads the meta table without reordering it", {
    meta <- meta_fixture(
        norm_id       = c("00010", "00020"),
        n_omics       = c(1L, 2L),
        combined_pval = c(1e-9, 1e-5)
    )
    before <- meta

    select_multi_omics_pathways(meta, top_n = 2)

    # The CSV the report links keeps its combined_pval ordering; only the figure
    # selects differently.
    expect_identical(meta, before)
})


# ---- what a blank cell means ------------------------------------------------

test_that("an untested layer stays NA through the heatmap preparation", {
    skip_if_not_installed("pheatmap")

    # 00020 was never tested by proteomics. Flattening that to 0 rendered it the
    # same white as a pathway the layer did test and found nothing in, and left
    # pheatmap's na_col as dead configuration.
    meta <- meta_fixture(
        norm_id       = c("00010", "00020"),
        n_omics       = c(2L, 1L),
        combined_pval = c(1e-5, 1e-3),
        pval_a        = c(1e-5, 1e-3),
        pval_b        = c(1e-4, NA_real_)
    )

    out <- withr::local_tempfile(fileext = ".png")
    grDevices::png(out)
    on.exit(grDevices::dev.off(), add = TRUE)

    # The cap step assigns through a logical subscript, and a logical subscript
    # carrying NA is an error in `[<-`. With the union candidate universe every
    # untested cell is NA, so this is the path a real run now takes.
    expect_no_error(plot_cross_omics_pathway_heatmap(meta, c("transcriptomics",
                                                             "proteomics")))
})

test_that("the heatmap no longer flattens missing values to zero", {
    body_src <- paste(deparse(body(plot_cross_omics_pathway_heatmap)),
                      collapse = " ")

    expect_false(grepl("is.na(log_pval_matrix)] <- 0", body_src, fixed = TRUE))
    # And the cap it does apply is guarded, which is what keeps NA survivable.
    expect_true(grepl("!is.na(log_pval_matrix)", body_src, fixed = TRUE))
})
