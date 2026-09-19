# Which rows the cross-omics meta figures lead with.
#
# Ordering by combined p-value alone lets the layer with the largest gene-set
# collection fill every slot with single-layer pathways. Counting the layers
# that merely *held* a p-value does not fix that -- it rewards coverage, so a
# pathway every layer measured and none found anything in outranks one
# genuinely significant in two.
#
# The meta heatmap and dot plot therefore rank on NOMINAL SUPPORT: the layers
# whose own raw p-value for the pathway is below 0.05. Ties go to the combined
# p-value, then to the pathway identity.
#
# Nominal support is a display-ranking definition and nothing more. It sets no
# significance threshold, joins no tested family, and touches no p-value,
# adjustment, membership or exported ordering. The rule is exactly p < 0.05 --
# a p-value sitting on the threshold is not support -- and a layer with no
# p-value contributes zero.
#
# The ORA FDR heatmap (#208) ranks on its own adjusted-p columns and is
# deliberately untouched; the last test here holds that boundary.
#
# All fixtures synthetic.

nominal_meta_fixture <- function(norm_id, combined_pval,
                         pval_transcriptomics = NA_real_,
                         pval_proteomics = NA_real_) {
    data.frame(
        norm_id = norm_id,
        pathway = paste("Pathway", norm_id),
        pval_transcriptomics = pval_transcriptomics,
        pval_proteomics = pval_proteomics,
        combined_pval = combined_pval,
        stringsAsFactors = FALSE
    )
}

# What the two display callers do, without drawing anything.
select_as_figure <- function(meta, top_n = 30) {
    meta$n_nominal_support <- .nominal_support(meta)
    select_multi_omics_pathways(meta, top_n,
                                count_col = "n_nominal_support",
                                id_col = "norm_id")
}


# ---- what nominal support counts -------------------------------------------

test_that("support counts the layers under the threshold, not the layers present", {
    # 00010 was measured by both layers and found by neither; 00020 is under
    # 0.05 in both. The old rule scored these equally.
    meta <- nominal_meta_fixture(
        norm_id = c("00010", "00020"),
        pval_transcriptomics = c(0.19, 0.01),
        pval_proteomics      = c(0.98, 0.02),
        combined_pval        = c(0.4, 0.001)
    )

    expect_equal(.nominal_support(meta), c(0L, 2L))
})

test_that("a missing p-value contributes no support", {
    meta <- nominal_meta_fixture(
        norm_id = c("00010", "00020"),
        pval_transcriptomics = c(0.01, NA_real_),
        pval_proteomics      = c(NA_real_, NA_real_),
        combined_pval        = c(0.01, 0.5)
    )

    expect_equal(.nominal_support(meta), c(1L, 0L))
})

test_that("a p-value exactly at the threshold is not support", {
    # The rule is p < 0.05, not <=. Pinned because the boundary is the part
    # most easily changed by accident.
    meta <- nominal_meta_fixture(
        norm_id = c("00010", "00020"),
        pval_transcriptomics = c(0.05, 0.049999),
        combined_pval        = c(0.05, 0.05)
    )

    expect_equal(.nominal_support(meta), c(0L, 1L))
})

test_that("a table with no per-layer p-value columns scores zero without erroring", {
    meta <- data.frame(norm_id = c("00010", "00020"),
                       combined_pval = c(0.01, 0.02),
                       stringsAsFactors = FALSE)

    expect_no_error(support <- .nominal_support(meta))
    expect_equal(support, c(0L, 0L))
})

test_that("the threshold is a local default, not a fixed constant", {
    meta <- nominal_meta_fixture(norm_id = "00010",
                         pval_transcriptomics = 0.08,
                         combined_pval = 0.08)

    expect_equal(.nominal_support(meta), 0L)
    expect_equal(.nominal_support(meta, alpha = 0.1), 1L)
})


# ---- what the figures lead with ---------------------------------------------

test_that("two nominally supportive layers outrank a stronger single-layer pathway", {
    # The regression. 00010's combined p is six orders of magnitude smaller,
    # and it is under 0.05 in one layer only.
    meta <- nominal_meta_fixture(
        norm_id = c("00010", "00020"),
        pval_transcriptomics = c(1e-9, 0.01),
        pval_proteomics      = c(0.9, 0.02),
        combined_pval        = c(1e-9, 1e-3)
    )

    expect_identical(select_as_figure(meta)$norm_id, c("00020", "00010"))
    # And with room for one row, it is the supported pathway that is shown.
    expect_identical(select_as_figure(meta, top_n = 1)$norm_id, "00020")
})

test_that("within equal support the combined p-value decides", {
    meta <- nominal_meta_fixture(
        norm_id = c("00010", "00020", "00030"),
        pval_transcriptomics = c(0.01, 0.01, 0.01),
        pval_proteomics      = c(0.01, 0.01, 0.01),
        combined_pval        = c(0.03, 0.001, 0.02)
    )

    expect_identical(select_as_figure(meta)$norm_id,
                     c("00020", "00030", "00010"))
})

test_that("a complete tie resolves on identity and survives a reversed table", {
    # Same support, same combined p. With arrival order as the last key this
    # flips when the rows are bound the other way round.
    meta <- nominal_meta_fixture(
        norm_id = c("00020", "00010"),
        pval_transcriptomics = c(0.01, 0.01),
        pval_proteomics      = c(0.02, 0.02),
        combined_pval        = c(1e-4, 1e-4)
    )
    reversed <- meta[rev(seq_len(nrow(meta))), , drop = FALSE]

    expect_identical(select_as_figure(meta)$norm_id, c("00010", "00020"))
    expect_identical(select_as_figure(meta)$norm_id,
                     select_as_figure(reversed)$norm_id)
})


# ---- the boundaries this must not cross -------------------------------------

test_that("id_col = NULL leaves the selector ordering exactly as it was", {
    # Every caller predating this argument relies on the incoming order as the
    # final key. Two rows tied on both ranking keys must come back in the order
    # they arrived, not sorted by anything else.
    ranking <- data.frame(row = 1:2,
                          n_ora_layers = c(2L, 2L),
                          best_ora_padj = c(0.01, 0.01),
                          stringsAsFactors = FALSE)

    got <- select_multi_omics_pathways(ranking, top_n = 2,
                                       count_col = "n_ora_layers",
                                       score_col = "best_ora_padj")

    expect_identical(got$row, 1:2)
    expect_identical(
        select_multi_omics_pathways(ranking[2:1, , drop = FALSE], top_n = 2,
                                    count_col = "n_ora_layers",
                                    score_col = "best_ora_padj")$row,
        2:1)
})

test_that("the ORA figure's own ranking is untouched by the new argument", {
    # #208 ranks on the adjusted-p matrix, which nominal support must not
    # reach. Its frame carries no pval_* columns and no identity, so it keeps
    # both its columns and its incoming-order tie-break.
    ranking <- data.frame(row = 1:3,
                          n_ora_layers = c(1L, 2L, 2L),
                          best_ora_padj = c(1e-9, 0.02, 0.01),
                          stringsAsFactors = FALSE)

    got <- select_multi_omics_pathways(ranking, top_n = 2,
                                       count_col = "n_ora_layers",
                                       score_col = "best_ora_padj")

    expect_identical(got$row, c(3L, 2L))
    expect_equal(.nominal_support(ranking), c(0L, 0L, 0L))
})

test_that("an unknown id_col falls back rather than erroring", {
    meta <- nominal_meta_fixture(
        norm_id = c("00010", "00020"),
        pval_transcriptomics = c(0.01, 0.01),
        combined_pval        = c(1e-4, 1e-4)
    )
    meta$n_nominal_support <- .nominal_support(meta)

    expect_no_error(
        got <- select_multi_omics_pathways(meta, top_n = 2,
                                           count_col = "n_nominal_support",
                                           id_col = "no_such_column"))
    expect_identical(got$norm_id, c("00010", "00020"))
})

test_that("display selection leaves the caller's table alone", {
    # The figures attach the support count to their own copy. The exported
    # meta-analysis table keeps its columns, its rows and its own ordering.
    meta <- nominal_meta_fixture(
        norm_id = c("00010", "00020", "00030"),
        pval_transcriptomics = c(1e-9, 0.01, 0.4),
        pval_proteomics      = c(0.9, 0.02, 0.3),
        combined_pval        = c(1e-9, 1e-3, 0.5)
    )
    before <- meta

    selected <- select_as_figure(meta, top_n = 2)

    expect_identical(meta, before)
    expect_false("n_nominal_support" %in% names(meta))
    # The selection reordered its own result, not the source.
    expect_identical(meta$norm_id, c("00010", "00020", "00030"))
    expect_identical(selected$norm_id, c("00020", "00010"))
})

test_that("an empty or single-row table is handled without special-casing", {
    empty <- nominal_meta_fixture(character(0), numeric(0),
                          numeric(0), numeric(0))
    expect_equal(nrow(select_as_figure(empty)), 0L)

    one <- nominal_meta_fixture("00010", 0.01, 0.01, 0.02)
    expect_identical(select_as_figure(one)$norm_id, "00010")
})
