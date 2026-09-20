# One cross-omics heatmap per gene-set collection, not one combined figure.
#
# The collections cover neither the same omics layers nor the same identifier
# space -- GO is annotated for the gene layers only, while the KEGG map space is
# the one every layer can share -- and they differ by an order of magnitude in
# size. A single figure with a flat top-N therefore gave every row to whichever
# collection was largest, and read as an analysis of that one alone.
#
# Two properties this holds to:
#
#   A figure shows only the layers that scored something in its collection. An
#   empty column would invite "tested here and found nothing", which is exactly
#   the inference these figures must not support.
#
#   The output contract changed from two fixed filenames to one pair per
#   collection, so a previous run's files -- the old combined ones, or a
#   collection this run does not produce -- would still be found by the report's
#   glob and shown as current. They are cleared, not merely overwritten.
#
# All fixtures synthetic.

meta_slice <- function(norm_id, ...) {
    layers <- list(...)
    n <- length(norm_id)
    out <- data.frame(norm_id = norm_id,
                      pathway = paste("Pathway", norm_id),
                      stringsAsFactors = FALSE)
    for (nm in names(layers)) {
        out[[paste0("pval_", nm)]] <- rep(layers[[nm]], length.out = n)
    }
    out
}


# ---- splitting the table by collection --------------------------------------

test_that("the meta table is classified on its accession, not on its label", {
    # By the time these figures are drawn, `pathway` holds a readable name that
    # attach_pathway_display_names() put there. Classifying on it would put
    # every row in "Other" and produce one figure named after nothing, so the
    # identity column is passed explicitly.
    meta <- meta_slice(c("00010", "GO:0006915", "PF00069"),
                       transcriptomics = 0.01)
    meta$pathway <- c("Glycolysis / Gluconeogenesis", "Apoptosis",
                      "Protein kinase domain")

    on_labels <- classify_pathway_collection(meta, "hsa")
    on_ids <- classify_pathway_collection(meta, "hsa",
                                          keys = as.character(meta$norm_id))

    # The label-driven default is exactly the trap: everything is "Other".
    expect_true(all(on_labels == "Other"))
    expect_identical(on_ids, c("KEGG", "GO", "Pfam"))
})

test_that("a collection nothing belongs to simply does not appear", {
    meta <- meta_slice(c("00010", "00020"), transcriptomics = 0.01)

    expect_setequal(
        unique(classify_pathway_collection(meta, "hsa",
                                           keys = as.character(meta$norm_id))),
        "KEGG")
})


# ---- which layers a per-collection figure may show --------------------------

test_that("a layer with no p-value in this collection is left out of its figure", {
    # GO is annotated for the gene layers only. Metabolomics has a column, and
    # nothing in it.
    sub <- meta_slice(c("GO:0006915", "GO:0006412"),
                      transcriptomics = 0.01, proteomics = 0.04,
                      metabolomics = NA_real_)

    expect_identical(
        .layers_with_values(sub, c("transcriptomics", "proteomics", "metabolomics")),
        c("transcriptomics", "proteomics"))
})

test_that("a layer with even one value is kept, and layer order is preserved", {
    sub <- meta_slice(c("00010", "00020"),
                      transcriptomics = c(NA_real_, 0.02),
                      metabolomics = 0.01)

    expect_identical(
        .layers_with_values(sub, c("metabolomics", "transcriptomics")),
        c("metabolomics", "transcriptomics"))
})

test_that("a layer with no column at all is not an error", {
    sub <- meta_slice("00010", transcriptomics = 0.01)

    expect_identical(.layers_with_values(sub, c("transcriptomics", "proteomics")),
                     "transcriptomics")
    expect_length(.layers_with_values(sub, "proteomics"), 0L)
})


# ---- the filename a collection maps to --------------------------------------

test_that("distinct collections cannot collapse onto one filename", {
    # Two collections sharing a slug would write one file, and the second would
    # be read as the first. The classifier's vocabulary is fixed, so this pins
    # injectivity over exactly that vocabulary.
    known <- c("KEGG", "GO", "Pfam", "InterPro", "Other")
    slugs <- vapply(known, .collection_slug, character(1))

    expect_equal(length(unique(slugs)), length(known))
    expect_false(any(grepl("[^A-Za-z0-9_]", slugs)))
})

test_that("a slug survives the round trip the report uses to label a figure", {
    # The report recovers the collection name from the filename by turning
    # underscores back into spaces, so a slug must not begin or end with one.
    expect_identical(.collection_slug("InterPro"), "InterPro")
    expect_identical(.collection_slug("Gene Ontology"), "Gene_Ontology")
    expect_false(grepl("^_|_$", .collection_slug(" KEGG ")))
})


# ---- what a previous run leaves behind ---------------------------------------

test_that("the old combined figures are removed when the contract changes", {
    dir <- withr::local_tempdir()
    old <- file.path(dir, c("cross_omics_pathway_heatmap.png",
                            "cross_omics_ora_heatmap.png"))
    file.create(old)

    .clear_collection_heatmaps(dir)

    expect_false(any(file.exists(old)))
})

test_that("a collection this run no longer produces does not survive it", {
    dir <- withr::local_tempdir()
    stale <- file.path(dir, c("cross_omics_pathway_heatmap_Pfam.png",
                              "cross_omics_ora_heatmap_Pfam.png"))
    file.create(stale)

    .clear_collection_heatmaps(dir)

    expect_false(any(file.exists(stale)))
})

test_that("nothing else in the directory is touched", {
    # Deliberately narrow: this clears two figure families in one directory and
    # is not a general output-lifecycle sweep.
    dir <- withr::local_tempdir()
    keep <- file.path(dir, c("cross_omics_enrichment_dotplot.png",
                             "cross_omics_pathways_meta_analysis.csv",
                             "transcriptomics_top_pathways.png",
                             "cross_omics_pathway_heatmap_KEGG.csv"))
    file.create(keep)
    file.create(file.path(dir, "cross_omics_pathway_heatmap_KEGG.png"))

    .clear_collection_heatmaps(dir)

    expect_true(all(file.exists(keep)))
})

test_that("clearing an empty or absent directory is silent and harmless", {
    dir <- withr::local_tempdir()

    expect_length(.clear_collection_heatmaps(dir), 0L)
    expect_length(.clear_collection_heatmaps(file.path(dir, "no_such_dir")), 0L)
})
