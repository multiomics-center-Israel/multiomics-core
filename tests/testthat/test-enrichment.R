# Tests for enrichment and pathway functions
# test-enrichment.R

# --- Tests for read_gmt ---

test_that("read_gmt reads GMT file correctly", {
    # Create temp GMT file
    tmp <- tempfile(fileext = ".gmt")
    writeLines(c(
        "PATHWAY_A\tDescription A\tGENE1\tGENE2\tGENE3",
        "PATHWAY_B\tDescription B\tGENE2\tGENE4",
        "PATHWAY_C\tNA\tGENE5\tGENE6\tGENE7\tGENE8"
    ), tmp)

    result <- read_gmt(tmp)
    unlink(tmp)

    expect_true(is.list(result))
    expect_equal(length(result), 3)
    expect_equal(names(result), c("PATHWAY_A", "PATHWAY_B", "PATHWAY_C"))
    expect_equal(result$PATHWAY_A, c("GENE1", "GENE2", "GENE3"))
    expect_equal(result$PATHWAY_B, c("GENE2", "GENE4"))
    expect_equal(length(result$PATHWAY_C), 4)
})

test_that("read_gmt stores descriptions as attribute", {
    tmp <- tempfile(fileext = ".gmt")
    writeLines(c(
        "PATHWAY_A\tFirst pathway\tGENE1\tGENE2",
        "PATHWAY_B\tSecond pathway\tGENE3"
    ), tmp)

    result <- read_gmt(tmp)
    unlink(tmp)

    descs <- attr(result, "descriptions")
    expect_false(is.null(descs))
    expect_equal(descs["PATHWAY_A"], c(PATHWAY_A = "First pathway"))
    expect_equal(descs["PATHWAY_B"], c(PATHWAY_B = "Second pathway"))
})

test_that("read_gmt handles empty genes gracefully", {
    tmp <- tempfile(fileext = ".gmt")
    writeLines(c(
        "PATHWAY_A\tDesc\tGENE1",
        "PATHWAY_EMPTY\tDesc\t"
    ), tmp)

    result <- read_gmt(tmp)
    unlink(tmp)

    expect_true(is.list(result))
    # At minimum PATHWAY_A should be there
    expect_true("PATHWAY_A" %in% names(result))
})

test_that("read_gmt errors on missing file", {
    expect_error(read_gmt("/nonexistent/path/fake.gmt"), "GMT file not found")
})

# --- Tests for load_gene_sets ---

test_that("load_gene_sets loads custom GMT file", {
    tmp <- tempfile(fileext = ".gmt")
    writeLines(c(
        "SET1\tDescription 1\tGENE1\tGENE2\tGENE3",
        "SET2\tDescription 2\tGENE4\tGENE5"
    ), tmp)

    result <- load_gene_sets(organism = "Homo sapiens", gmt_file = tmp)
    unlink(tmp)

    expect_true(is.list(result))
    # Custom GMT should be loaded under the "custom" key
    expect_true("custom" %in% names(result))
    expect_true(length(result$custom) > 0)
})

# --- Tests for add_pathway_names (RNA-seq) ---

test_that("add_pathway_names adds KEGG pathway names", {
    df <- data.frame(
        pathway = c("hsa04110", "hsa04151", "hsa00010"),
        padj = c(0.01, 0.02, 0.03),
        stringsAsFactors = FALSE
    )
    result <- add_pathway_names(df, database = "KEGG")
    expect_true("pathway_name" %in% names(result))
    # KEGG pathway_name is set to the pathway ID itself (no lookup)
    expect_equal(result$pathway_name, result$pathway)
})

test_that("add_pathway_names uses GMT descriptions for custom database", {
    df <- data.frame(
        pathway = c("SET1", "SET2"),
        padj = c(0.01, 0.02),
        stringsAsFactors = FALSE
    )
    gene_sets <- list(SET1 = c("G1", "G2"), SET2 = c("G3"))
    attr(gene_sets, "descriptions") <- c(SET1 = "First set", SET2 = "Second set")

    result <- add_pathway_names(df, database = "custom", gene_sets = gene_sets)
    expect_true("pathway_name" %in% names(result))
    expect_equal(result$pathway_name[1], "First set")
    expect_equal(result$pathway_name[2], "Second set")
})

test_that("add_pathway_names returns input unchanged when no pathway column", {
    df <- data.frame(
        gene_set = c("SET1", "SET2"),
        padj = c(0.01, 0.02),
        stringsAsFactors = FALSE
    )
    result <- add_pathway_names(df, database = "custom")
    expect_false("pathway_name" %in% names(result))
    expect_equal(ncol(result), ncol(df))
})

test_that("add_pathway_names returns NULL/empty input unchanged", {
    expect_null(add_pathway_names(NULL, database = "GO"))
    empty_df <- data.frame(pathway = character(0), padj = numeric(0), stringsAsFactors = FALSE)
    result <- add_pathway_names(empty_df, database = "GO")
    expect_equal(nrow(result), 0)
})

test_that("lookup_go_term_names returns names for valid GO IDs", {
    skip_if_not_installed("GO.db")
    result <- lookup_go_term_names(c("GO:0008150", "GO:0003674"))
    expect_true(is.character(result))
    expect_equal(length(result), 2)
    # GO:0008150 is "biological_process"
    expect_true(nzchar(result[1]))
})

test_that("lookup_go_term_names handles invalid GO IDs gracefully", {
    skip_if_not_installed("GO.db")
    result <- lookup_go_term_names(c("GO:9999999", "NOT_A_GO_ID"))
    expect_true(is.character(result))
    expect_equal(length(result), 2)
})

test_that("lookup_go_term_names returns empty for empty input", {
    result <- lookup_go_term_names(character(0))
    expect_true(is.character(result))
    expect_equal(length(result), 0)
})

test_that("add_pathway_names handles GO terms", {
    skip_if_not_installed("GO.db")
    df <- data.frame(
        pathway = c("GO:0008150", "GO:0003674"),
        padj = c(0.01, 0.02),
        stringsAsFactors = FALSE
    )
    result <- add_pathway_names(df, database = "GO")
    expect_true("pathway_name" %in% names(result))
    expect_true(nzchar(result$pathway_name[1]))
})

# --- Tests for map_compounds_for_enrichment (metabolomics ID mapping) ---

test_that("map_compounds_for_enrichment prefers KEGG IDs and aligns feature_map", {
    expr <- matrix(seq_len(12), nrow = 4,
                   dimnames = list(c("f1", "f2", "f3", "f4"), paste0("S", 1:3)))
    row_data <- data.frame(
        feature_id = c("f1", "f2", "f3", "f4"),
        Name       = c("Glucose", "", "Citrate", "Pyruvate"),
        KEGG       = c("C00031", "", "", "cpd:C00022"),
        stringsAsFactors = FALSE
    )

    res <- map_compounds_for_enrichment(row_data, expr)

    # KEGG IDs win where present (any prefix stripped); Name is the fallback.
    expect_equal(unname(res$feature_map["f1"]), "C00031")
    expect_equal(unname(res$feature_map["f4"]), "C00022")
    expect_equal(unname(res$feature_map["f3"]), "Citrate")

    # f2 has neither a name nor a KEGG id, so it drops out of the background —
    # but the map must still resolve the LATER features correctly. This is the
    # regression guard for the previous positional-misalignment bug, where
    # dropping f2 shifted every subsequent feature's compound by one.
    expect_false("" %in% res$compound_names)
    expect_true(all(c("C00031", "Citrate", "C00022") %in% res$compound_names))
})
