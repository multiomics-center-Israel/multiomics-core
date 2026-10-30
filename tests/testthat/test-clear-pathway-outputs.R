# tests/testthat/test-clear-pathway-outputs.R
#
# save_pathway_results() and generate_pathway_plots() never delete, so a
# collection that failed in a later run (e.g. a KEGG download error) kept the
# earlier run's CSV and dotplot, and the RNA report showed them as current.
# mod_rnaseq_pathway() now calls clear_pathway_outputs() first.

make_enrich_dir <- function() {
    d <- file.path(tempfile("enrich_"), "Enrichment")
    dir.create(file.path(d, "plots"), recursive = TRUE)
    dir.create(file.path(d, "GSEA"))
    file.create(file.path(d, c("pathway_A_vs_B_KEGG_fgsea.csv",
                               "pathway_volcano_data_A_vs_B.csv",
                               "gene_annotation.csv",
                               "notes.txt")))
    file.create(file.path(d, "plots", c("pathway_A_vs_B_KEGG_fgsea.png", "other.png")))
    file.create(file.path(d, "GSEA", "pathway_kept.csv"))
    d
}

test_that("clear_pathway_outputs removes only the pathway writers' files", {
    d <- make_enrich_dir()
    removed <- clear_pathway_outputs(d)

    expect_setequal(basename(removed), c("pathway_A_vs_B_KEGG_fgsea.csv",
                                         "pathway_volcano_data_A_vs_B.csv",
                                         "pathway_A_vs_B_KEGG_fgsea.png"))
    expect_false(any(file.exists(removed)))
    expect_true(file.exists(file.path(d, "notes.txt")))
    # The report's gene-symbol source stays until the offline route writes its
    # own (#338); removing it now would leave offline reports with raw IDs.
    expect_true(file.exists(file.path(d, "gene_annotation.csv")))
    expect_true(file.exists(file.path(d, "plots", "other.png")))
    # The offline route's subtrees are not touched
    expect_true(file.exists(file.path(d, "GSEA", "pathway_kept.csv")))
})

test_that("clear_pathway_outputs is a no-op on a missing or clean directory", {
    expect_length(clear_pathway_outputs(file.path(tempdir(), "no_such_enrichment_dir")), 0L)
    d <- file.path(tempfile("enrich_"), "Enrichment")
    dir.create(d, recursive = TRUE)
    expect_length(clear_pathway_outputs(d), 0L)
})

test_that("mod_rnaseq_pathway clears old pathway files even when enrichment is off", {
    out_dir <- dirname(make_enrich_dir())
    cfg <- list(modes = list(rna = list(pathway = list(enabled = FALSE))))
    res <- suppressMessages(mod_rnaseq_pathway(de_res = list(tables = list()), pre = NULL,
                                               config = cfg, out_dir = out_dir))
    expect_length(res$pathway_results, 0L)
    expect_false(file.exists(file.path(out_dir, "Enrichment", "pathway_A_vs_B_KEGG_fgsea.csv")))
})
