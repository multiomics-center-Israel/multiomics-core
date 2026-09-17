# The cross-omics candidate universe is the union of what the layers enriched.
#
# use_pathways is not a display choice: it is what merge_pathway_pvalues()
# assembles and stouffer_combined_pvalues() scores, so it decides which pathways
# reach the meta-analysis. It used to collapse to the intersection whenever five
# or more pathways were shared, which let the narrowest layer decide eligibility
# for every other.
#
# All fixtures are synthetic.

# Two layers sharing five pathways -- enough to have triggered the old collapse
# -- and each carrying one the other does not.
shared_ids <- c("00010", "00020", "00030", "00040", "00050")

omics_frame <- function(ids, pvalues) {
    data.frame(ID = ids, pathway = paste("Pathway", ids),
               pvalue = pvalues, padj = pvalues,
               stringsAsFactors = FALSE)
}

two_layers <- function() {
    list(
        transcriptomics = omics_frame(c(shared_ids, "00099"),
                                      c(0.001, 0.002, 0.003, 0.004, 0.005, 0.006)),
        proteomics      = omics_frame(c(shared_ids, "00088"),
                                      c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06))
    )
}

minimal_config <- function() {
    list(global = list(organism = "Unlisted nonmodel species"))
}


test_that("a pathway only one layer enriched still reaches the meta-analysis", {
    res <- suppressMessages(
        analyze_cross_omics_enrichment(two_layers(), minimal_config(), out_dir = NULL))

    skip_if(is.null(res) || is.null(res$meta_analysis), "no meta-analysis produced")
    ids <- res$meta_analysis$norm_id

    # The five shared pathways would have been the whole candidate set under the
    # old rule, because five is exactly where it collapsed.
    expect_true(all(shared_ids %in% ids))
    # And these are the ones the collapse discarded: present in one layer only.
    expect_true("00099" %in% ids)
    expect_true("00088" %in% ids)
})

test_that("the meta-analysis, not the candidate set, decides multi-omics evidence", {
    res <- suppressMessages(
        analyze_cross_omics_enrichment(two_layers(), minimal_config(), out_dir = NULL))
    skip_if(is.null(res) || is.null(res$meta_analysis), "no meta-analysis produced")
    meta <- res$meta_analysis

    # Union-only rows arrive with one layer's p-value and are counted as such,
    # which is what stouffer_combined_pvalues() exists to do. They are no longer
    # removed before it can count them.
    expect_equal(meta$n_omics[meta$norm_id == "00099"], 1L)
    expect_equal(meta$n_omics[meta$norm_id == "00088"], 1L)
    expect_true(all(meta$n_omics[meta$norm_id %in% shared_ids] == 2L))
})

test_that("the shared pathways behave exactly as they did before", {
    layers <- two_layers()
    res <- suppressMessages(
        analyze_cross_omics_enrichment(layers, minimal_config(), out_dir = NULL))
    skip_if(is.null(res) || is.null(res$meta_analysis), "no meta-analysis produced")
    meta <- res$meta_analysis

    # Widening the candidate set must not disturb the rows that were already in
    # it: each shared pathway still carries the per-layer p-values it came with.
    for (id in shared_ids) {
        row <- meta[meta$norm_id == id, ]
        expect_equal(row[["pval_transcriptomics"]],
                     layers$transcriptomics$pvalue[layers$transcriptomics$ID == id])
        expect_equal(row[["pval_proteomics"]],
                     layers$proteomics$pvalue[layers$proteomics$ID == id])
    }
})

test_that("the overlap is still reported, as a description rather than a gate", {
    res <- suppressMessages(
        analyze_cross_omics_enrichment(two_layers(), minimal_config(), out_dir = NULL))
    skip_if(is.null(res), "no result produced")

    expect_setequal(res$common_pathways, shared_ids)
    expect_setequal(res$union_pathways, c(shared_ids, "00099", "00088"))
})

test_that("fewer than five shared pathways behaves the same way it always did", {
    # The old rule already used the union here, so this case is a control: it
    # must be unchanged, and it pins that the union is now unconditional rather
    # than the other side of a threshold.
    layers <- list(
        transcriptomics = omics_frame(c("00010", "00099"), c(0.001, 0.002)),
        proteomics      = omics_frame(c("00010", "00088"), c(0.01, 0.02))
    )
    res <- suppressMessages(
        analyze_cross_omics_enrichment(layers, minimal_config(), out_dir = NULL))
    skip_if(is.null(res) || is.null(res$meta_analysis), "no meta-analysis produced")

    expect_setequal(res$meta_analysis$norm_id, c("00010", "00099", "00088"))
})
