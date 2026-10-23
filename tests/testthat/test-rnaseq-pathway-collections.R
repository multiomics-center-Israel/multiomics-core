# tests/testthat/test-rnaseq-pathway-collections.R
#
# A KEGG download error only warned inside load_gene_sets(), so the RNA report
# omitted every KEGG section while its methods text still named KEGG (taken
# from the config). mod_rnaseq_pathway() now records what loaded in
# Enrichment/pathway_collections.tsv and the report reads it through
# describe_rna_pathway_collections().

test_that("summarise_pathway_collections marks a requested collection that did not load", {
    tab <- summarise_pathway_collections(c("GO", "KEGG"),
                                         list(GO_BP = list(), GO_CC = list(), GO_MF = list()))
    expect_equal(tab$collection, c("GO", "KEGG"))
    expect_equal(tab$status, c("loaded", "not loaded"))
})

test_that("summarise_pathway_collections expands 'all' and keeps custom GMTs", {
    tab <- summarise_pathway_collections("all", list(GO_BP = list(), KEGG = list(),
                                                     custom = list()))
    expect_equal(tab$collection, c("GO", "KEGG", "Reactome", "custom"))
    expect_equal(tab$status, c("loaded", "loaded", "not loaded", "loaded"))
})

test_that("summarise_pathway_collections accepts a single GO ontology request", {
    tab <- summarise_pathway_collections("GO_BP", list(GO_BP = list()))
    expect_equal(tab$collection, "GO_BP")
    expect_equal(tab$status, "loaded")
})

test_that("summarise_pathway_collections handles no gene sets at all", {
    tab <- summarise_pathway_collections(c("GO", "KEGG"), list())
    expect_equal(tab$status, c("not loaded", "not loaded"))
})

test_that("describe_rna_pathway_collections reads the recorded outcome", {
    d <- tempfile("enrich_")
    dir.create(d)
    utils::write.table(data.frame(collection = c("GO", "KEGG"),
                                  status = c("loaded", "not loaded")),
                       file.path(d, "pathway_collections.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    desc <- describe_rna_pathway_collections(d, list(method = "both",
                                                     databases = c("GO", "KEGG")))
    expect_true(desc$recorded)
    expect_equal(desc$method, "fGSEA and ORA")
    expect_equal(desc$collections, "GO")
    expect_equal(desc$missing, "KEGG")
    expect_match(desc$missing_note, "KEGG was requested but could not be loaded")
})

test_that("describe_rna_pathway_collections falls back to the config without a record", {
    desc <- describe_rna_pathway_collections(tempfile("no_enrich_"), list(method = "fgsea"))
    expect_false(desc$recorded)
    expect_equal(desc$method, "fGSEA")
    # Same default as mod_rnaseq_pathway()
    expect_equal(desc$collections, "GO and KEGG")
    expect_equal(desc$missing_note, "")
})

test_that("the RNA report reads pathway collections through one resolver", {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "rnaseq", "report_template.Rmd"),
        "R/domain/rnaseq/report_template.Rmd"
    )
    f <- candidates[file.exists(candidates)][1]
    if (is.na(f)) stop("R/domain/rnaseq/report_template.Rmd not found from the test directory")
    src <- readLines(f, warn = FALSE)
    expect_false(any(grepl("rna_cfg$pathway$databases", src, fixed = TRUE)))
    expect_gte(sum(grepl("describe_rna_pathway_collections(enrich_dir, rna_cfg$pathway)",
                         src, fixed = TRUE)), 3L)
})
