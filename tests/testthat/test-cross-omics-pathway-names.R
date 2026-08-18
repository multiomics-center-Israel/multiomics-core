# tests/testthat/test-cross-omics-pathway-names.R
#
# Unit tests for the cross-omics enrichment labelling and row-selection helpers
# in R/domain/multiomics/07_enrichment.R:
#   - attach_pathway_names_from_tables()
#   - format_pathway_labels()
#   - pathway_name_column()
#   - select_multi_omics_pathways()
#   - get_kegg_organism()  (defensive input handling)
#
# All data here is synthetic: made-up accessions and made-up term names.

make_pathway_tables <- function() {
    list(
        transcriptomics = data.frame(
            pathway      = c("map00001", "map00002", "GO:0000001"),
            pathway_name = c("Alpha metabolism", "Beta metabolism", "gamma process"),
            pvalue       = c(0.01, 0.20, 0.30),
            stringsAsFactors = FALSE
        ),
        proteomics = data.frame(
            pathway      = c("map00002", "GO:0000002"),
            pathway_name = c("Beta metabolism", "delta process"),
            pvalue       = c(0.02, 0.04),
            stringsAsFactors = FALSE
        )
    )
}

make_meta <- function() {
    data.frame(
        pathway          = c("map00001", "map00002", "GO:0000001", "GO:0000002"),
        pval_transcriptomics = c(0.01, 0.20, 0.30, NA),
        pval_proteomics      = c(NA,   0.02, NA,   0.04),
        combined_pval    = c(0.01, 0.05, 0.30, 0.04),
        n_omics          = c(1L, 2L, 1L, 1L),
        stringsAsFactors = FALSE
    )
}


# --- attach_pathway_names_from_tables() --------------------------------------

test_that("pathway_name propagates from the per-omics tables into the meta table", {
    meta <- attach_pathway_names_from_tables(make_meta(), make_pathway_tables())

    expect_true("pathway_name" %in% names(meta))
    expect_equal(
        meta$pathway_name,
        c("Alpha metabolism", "Beta metabolism", "gamma process", "delta process")
    )
    # The join must not reorder or drop rows.
    expect_equal(meta$pathway, make_meta()$pathway)
    expect_equal(nrow(meta), nrow(make_meta()))
})

test_that("accessions named by no source table get NA, not a wrong name", {
    meta <- make_meta()
    meta <- rbind(meta, data.frame(
        pathway = "map09999",
        pval_transcriptomics = 0.5, pval_proteomics = NA,
        combined_pval = 0.5, n_omics = 1L,
        stringsAsFactors = FALSE
    ))

    out <- attach_pathway_names_from_tables(meta, make_pathway_tables())
    expect_true(is.na(out$pathway_name[out$pathway == "map09999"]))
})

test_that("attach_pathway_names_from_tables is a no-op when no table carries names", {
    tables <- lapply(make_pathway_tables(), function(df) df[, c("pathway", "pvalue")])
    meta <- make_meta()

    expect_identical(attach_pathway_names_from_tables(meta, tables), meta)
})

test_that("attach_pathway_names_from_tables tolerates an empty meta table", {
    empty <- make_meta()[0, , drop = FALSE]
    expect_equal(nrow(attach_pathway_names_from_tables(empty, make_pathway_tables())), 0)
})


# --- format_pathway_labels() -------------------------------------------------

test_that("labels fall back to the bare accession when the name is NA or empty", {
    ids   <- c("GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004")
    names <- c("gamma process", NA, "", "   ")

    labels <- format_pathway_labels(ids, names)

    expect_equal(labels[1], "gamma process")
    expect_equal(labels[2], "GO:0000002")
    expect_equal(labels[3], "GO:0000003")
    expect_equal(labels[4], "GO:0000004")
})

test_that("labels fall back to accessions when there is no name column at all", {
    ids <- c("map00001", "map00002")
    expect_equal(format_pathway_labels(ids, NULL), ids)
    expect_null(pathway_name_column(data.frame(pathway = ids, stringsAsFactors = FALSE)))
})

test_that("append_id keeps the accession visible and makes labels unique", {
    ids   <- c("GO:0000001", "GO:0000002")
    names <- c("shared name", "shared name")

    labels <- format_pathway_labels(ids, names, append_id = TRUE)

    expect_equal(labels, c("shared name (GO:0000001)", "shared name (GO:0000002)"))
    expect_equal(anyDuplicated(labels), 0L)
})

test_that("long names are truncated to max_chars, accession appended after", {
    long <- strrep("x", 80)
    label <- format_pathway_labels("GO:0000001", long, max_chars = 20)

    expect_equal(nchar(label), 20)
    expect_true(endsWith(label, "..."))

    with_id <- format_pathway_labels("GO:0000001", long, max_chars = 20, append_id = TRUE)
    expect_true(endsWith(with_id, " (GO:0000001)"))
})


# --- select_multi_omics_pathways() -------------------------------------------

test_that("row selection puts higher n_omics first, then lower combined_pval", {
    meta <- data.frame(
        pathway       = c("P1", "P2", "P3", "P4", "P5"),
        combined_pval = c(1e-12, 1e-11, 0.30, 0.01, 0.20),
        n_omics       = c(1L, 1L, 3L, 2L, 3L),
        stringsAsFactors = FALSE
    )

    out <- select_multi_omics_pathways(meta, top_n = 5)

    # 3-omics rows first (P5 before P3 on p-value), then the 2-omics row,
    # then the single-omics rows — despite those having the smallest p-values.
    expect_equal(out$pathway, c("P5", "P3", "P4", "P1", "P2"))
    expect_true(!is.unsorted(rev(out$n_omics)))
})

test_that("row selection honours top_n and does not mutate the input order", {
    meta <- data.frame(
        pathway       = c("P1", "P2", "P3"),
        combined_pval = c(1e-9, 0.5, 0.4),
        n_omics       = c(1L, 2L, 2L),
        stringsAsFactors = FALSE
    )

    out <- select_multi_omics_pathways(meta, top_n = 2)

    expect_equal(nrow(out), 2)
    expect_equal(out$pathway, c("P3", "P2"))
    # The meta table itself keeps its combined_pval ordering — the plot reorders,
    # the CSV contract does not.
    expect_equal(meta$pathway, c("P1", "P2", "P3"))
})

test_that("row selection handles empty input and a missing n_omics column", {
    meta <- data.frame(
        pathway       = c("P1", "P2"),
        combined_pval = c(0.5, 0.1),
        stringsAsFactors = FALSE
    )

    # Without n_omics every row ties, so ordering falls through to combined_pval.
    expect_equal(select_multi_omics_pathways(meta, 2)$pathway, c("P2", "P1"))
    expect_equal(nrow(select_multi_omics_pathways(meta[0, ], 5)), 0)
})

test_that("NA combined p-values sort last rather than dropping rows", {
    meta <- data.frame(
        pathway       = c("P1", "P2", "P3"),
        combined_pval = c(NA_real_, 0.2, 0.1),
        n_omics       = c(2L, 2L, 2L),
        stringsAsFactors = FALSE
    )

    out <- select_multi_omics_pathways(meta, 3)
    expect_equal(out$pathway, c("P3", "P2", "P1"))
})


# --- get_kegg_organism() -----------------------------------------------------

test_that("get_kegg_organism returns NULL for absent or unusable organisms", {
    expect_null(get_kegg_organism(NULL))
    expect_null(get_kegg_organism(NA))
    expect_null(get_kegg_organism(NA_character_))
    expect_null(get_kegg_organism(""))
    expect_null(get_kegg_organism("   "))
    expect_null(get_kegg_organism(character(0)))
    expect_null(get_kegg_organism(c("human", "mouse")))
    expect_null(get_kegg_organism("Spalangia endius"))
})

test_that("get_kegg_organism still resolves known organisms", {
    expect_equal(get_kegg_organism("human"), "hsa")
    expect_equal(get_kegg_organism("Homo sapiens"), "hsa")
    expect_equal(get_kegg_organism("c_elegans"), "cel")
})
