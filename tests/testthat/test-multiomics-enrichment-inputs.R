# tests/testthat/test-multiomics-enrichment-inputs.R
#
# How the cross-omics enrichment step reads each layer's results:
#   - extract_enrichment_df() accepts the RNA ($pathway_results) and proteomics
#     (bare contrast-keyed list) shapes, but leaves the metabolomics wrapper to
#     the KEGG-from-DE fallback.
#   - MultiGSEA aligns omics on the `pathway` ID and labels with `pathway_name`,
#     with no row-name fallback.
#   - get_kegg_organism() returns NULL for a missing organism and errors on
#     several.
# All data are synthetic.

pw_table <- function(ids, padj, names = NULL) {
    df <- data.frame(pathway = ids, padj = padj, stringsAsFactors = FALSE)
    if (!is.null(names)) df$pathway_name <- names
    df
}

# --- extract_enrichment_df ----------------------------------------------------

test_that("extract_enrichment_df binds a bare contrast-keyed proteomics result", {
    res <- list(
        A_vs_B = list(GO_fgsea  = pw_table(c("GO:1", "GO:2"), c(0.01, 0.2)),
                      GO_ora_up = pw_table("GO:3", 0.03)),
        C_vs_B = list(GO_fgsea  = pw_table("GO:4", 0.04))
    )
    df <- extract_enrichment_df(res)
    expect_s3_class(df, "data.frame")
    expect_setequal(df$pathway, c("GO:1", "GO:2", "GO:3", "GO:4"))
})

test_that("extract_enrichment_df still reads the RNA $pathway_results slot", {
    res <- list(
        annotation      = NULL,
        pathway_results = list(A_vs_B = list(
            KEGG_fgsea = pw_table(c("hsa1", "hsa2"), c(0.01, 0.2))
        )),
        plot_files      = list()
    )
    expect_setequal(extract_enrichment_df(res)$pathway, c("hsa1", "hsa2"))
})

test_that("extract_enrichment_df leaves the metabolomics wrapper to the KEGG fallback", {
    # Shape of mod_metabolomics_enrichment(): method tables carry raw_p / p_value
    # and FDR, which merge_pathway_pvalues() cannot read.
    metab <- list(
        qea                = list(table = data.frame(ID = "map00010", raw_p = 0.01, FDR = 0.05)),
        ssgsea             = list(table = data.frame(ID = "map00020", p_value = 0.02, FDR = 0.1)),
        ora                = data.frame(ID = "map00030", p_value = 0.03, FDR = 0.2),
        gsea               = data.frame(ID = "map00040", p_value = 0.04, FDR = 0.3),
        qea_by_contrast    = list(),
        ssgsea_by_contrast = list(),
        ora_by_contrast    = list(),
        gsea_by_contrast   = list(),
        plots              = list(),
        files              = "enrichment_qea_results.tsv"
    )
    expect_null(extract_enrichment_df(metab))
})

# --- MultiGSEA term helpers ---------------------------------------------------

test_that(".multigsea_term_ids keys on pathway when there is no term or ID column", {
    df <- pw_table(c("GO:1", "GO:2"), c(0.01, 0.2), names = c("name 1", "name 2"))
    expect_identical(.multigsea_term_ids(df), c("GO:1", "GO:2"))
})

test_that(".multigsea_term_ids returns NULL instead of falling back to row names", {
    expect_null(.multigsea_term_ids(data.frame(padj = c(0.1, 0.2))))
})

test_that(".multigsea_term_names reads pathway_name and the older pathway/ID shape", {
    current <- pw_table(c("GO:1", "GO:2"), c(0.01, 0.2), names = c("name 1", "name 2"))
    older   <- data.frame(ID = "hsa00010", pathway = "Glycolysis", stringsAsFactors = FALSE)
    nm <- .multigsea_term_names(list(current, older))
    expect_identical(nm[["GO:1"]], "name 1")
    expect_identical(nm[["hsa00010"]], "Glycolysis")
})

test_that(".save_multigsea_pair_plot aligns two omics on the pathway ID, not row order", {
    skip_if_not_installed("ggplot2")

    ids  <- c("GO:1", "GO:2", "GO:3", "GO:4")
    padj <- c(0.001, 0.01, 0.2, 0.9)
    rna  <- pw_table(ids, padj, names = paste("name", 1:4))
    # Same pathways and values, stored in a different row order.
    ord  <- c(3, 1, 4, 2)
    prot <- pw_table(ids[ord], padj[ord], names = paste("name", 1:4)[ord])

    out_dir <- tempfile("test-multigsea-pair-")
    dir.create(out_dir)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

    suppressMessages(suppressWarnings(
        .save_multigsea_pair_plot(rna, prot, "transcriptomics", "proteomics",
                                  out_dir = out_dir)
    ))
    written <- utils::read.csv(
        file.path(out_dir, "multigsea_transcriptomics_vs_proteomics.csv"),
        stringsAsFactors = FALSE
    )
    expect_setequal(written$term, ids)
    expect_equal(written$x, written$y)
    expect_setequal(written$label, paste("name", 1:4))
})

# --- get_kegg_organism --------------------------------------------------------

test_that("get_kegg_organism returns NULL for a missing or blank organism", {
    expect_null(get_kegg_organism(NULL))
    expect_null(get_kegg_organism(character(0)))
    expect_null(get_kegg_organism(NA_character_))
    expect_null(get_kegg_organism(""))
})

test_that("get_kegg_organism maps a known organism name", {
    expect_identical(get_kegg_organism("human"), "hsa")
    expect_identical(get_kegg_organism("Homo sapiens"), "hsa")
})

test_that("get_kegg_organism errors on several organisms rather than skipping KEGG", {
    expect_error(get_kegg_organism(c("human", "mouse")), "single organism")
})
