# tests/testthat/test-pathway-volcano-duplicate-ids.R
#
# build_pathway_volcano_data() (R/domain/rnaseq/07_pathway.R) keys each
# significant gene set by its pathway label. Gene-set collections are now scored
# one per GMT file, so two of them can carry the same pathway id and label with
# different memberships; the later one used to replace the earlier one under
# that single key and its leading-edge genes disappeared from the export.
#
# All data are synthetic.

volcano_de <- function(ids) {
    data.frame(FeatureID = ids, log2FoldChange = seq_along(ids) / 10,
               pvalue = 0.001, padj = 0.01, stringsAsFactors = FALSE)
}

gsea_result <- function(pathway, label, genes) {
    data.frame(pathway = pathway, pathway_name = label, padj = 0.001,
               leadingEdge = paste(genes, collapse = ","),
               stringsAsFactors = FALSE)
}

test_that("the same pathway in two collections keeps both memberships", {
    de <- volcano_de(paste0("g", 1:6))

    # One id, one label, two collections, different leading edges.
    pathway_results <- list(contrast_A = list(
        KEGG_source_one_fgsea = gsea_result("map00010", "Glycolysis", c("g1", "g2")),
        KEGG_source_two_fgsea = gsea_result("map00010", "Glycolysis", c("g5", "g6"))
    ))

    out <- build_pathway_volcano_data(de, pathway_results)

    # g1 belongs only to the first collection's set, g5 only to the second's.
    # If one had replaced the other, one of these would be empty.
    expect_true(nzchar(out$pathways[out$FeatureID == "g1"]))
    expect_true(nzchar(out$pathways[out$FeatureID == "g5"]))

    # The second is qualified by the result it came from, so the two are
    # distinguishable in the export rather than merged.
    expect_match(out$pathways[out$FeatureID == "g5"], "KEGG_source_two_fgsea",
                 fixed = TRUE)
})

test_that("a label that does not collide is not qualified", {
    de <- volcano_de(paste0("g", 1:4))

    pathway_results <- list(contrast_A = list(
        KEGG_source_one_fgsea = gsea_result("map00010", "Glycolysis", c("g1", "g2")),
        PFAM_source_fgsea     = gsea_result("PF00089", "Trypsin", c("g3", "g4"))
    ))

    out <- build_pathway_volcano_data(de, pathway_results)

    expect_equal(out$pathways[out$FeatureID == "g1"], "Glycolysis")
    expect_equal(out$pathways[out$FeatureID == "g3"], "Trypsin")
})

test_that("the same pathway with identical membership is not duplicated", {
    de <- volcano_de(paste0("g", 1:4))

    # Two collections agreeing on the set: there is nothing to lose, so the
    # label stays as it was rather than gaining a source suffix.
    pathway_results <- list(contrast_A = list(
        source_one_fgsea = gsea_result("map00010", "Glycolysis", c("g1", "g2")),
        source_two_fgsea = gsea_result("map00010", "Glycolysis", c("g1", "g2"))
    ))

    out <- build_pathway_volcano_data(de, pathway_results)

    expect_equal(out$pathways[out$FeatureID == "g1"], "Glycolysis")
})
