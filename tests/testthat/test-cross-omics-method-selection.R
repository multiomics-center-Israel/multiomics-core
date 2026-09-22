# Each layer feeds the cross-omics meta-analysis with one kind of test.
#
# merge_pathway_pvalues() used to pick a p-value column for the whole table,
# `pvalue` before `pval`. A layer holding ORA rows in `pvalue` beside fgsea rows
# in `pval` therefore contributed ORA only, while a layer scored by GSEA alone
# contributed GSEA, and Stouffer combined unlike tests. The layer now
# contributes its rank-based rows where it has any, otherwise its ORA rows, and
# the table records which.
#
# All fixtures are synthetic.

gene_layer <- function() {
    # ORA rows carry `pvalue`, fgsea rows `pval`, as run_pathway_analysis()
    # writes them; bind_rows() NA-fills the other column.
    dplyr::bind_rows(
        data.frame(ID = c("map00010", "map00020"), pvalue = c(1e-6, 0.3),
                   padj = c(1e-5, 0.5), method = "ora",
                   stringsAsFactors = FALSE),
        data.frame(ID = c("map00010", "map00020"), pval = c(0.04, 0.2),
                   padj = c(0.1, 0.4), NES = c(1.8, -0.9), method = "fgsea",
                   stringsAsFactors = FALSE)
    )
}

ora_only_layer <- function() {
    data.frame(ID = c("map00010", "map00030"), pvalue = c(0.01, 0.02),
               padj = c(0.05, 0.05), method = "ora", stringsAsFactors = FALSE)
}

gsea_table <- function() {
    data.frame(pathway = c("Glycolysis", "Citrate cycle"),
               ID = c("map00010", "map00020"), pvalue = c(0.03, 0.6),
               padj = c(0.2, 0.8), NES = c(-1.5, 0.4), method = "fgsea",
               contrast = "A_vs_B", stringsAsFactors = FALSE)
}

nonmodel_config <- function() {
    list(global = list(organism = "Unlisted nonmodel species"))
}


# ---- the per-row p-value --------------------------------------------------

test_that("raw p-values are read per row across pvalue and pval", {
    df <- data.frame(pvalue = c(0.01, NA), pval = c(NA, 0.02))
    expect_equal(.raw_p_values(df), c(0.01, 0.02))
})

test_that("an adjusted column is used only when no raw column exists", {
    expect_equal(.raw_p_values(data.frame(padj = c(0.1, 0.2))), c(0.1, 0.2))
    expect_equal(.raw_p_values(data.frame(pvalue = c(0.01, NA), padj = c(0.5, 0.6))),
                 c(0.01, NA))
})

test_that("a table with no p-value column at all gives NULL", {
    expect_null(.raw_p_values(data.frame(ID = "map00010")))
})


# ---- choosing the method --------------------------------------------------

test_that("rank-based rows win over ORA rows in the same layer", {
    df <- gene_layer()
    sel <- select_layer_method_rows(df, .raw_p_values(df))
    expect_identical(sel$method, "fgsea")
    expect_identical(df$method[sel$keep], c("fgsea", "fgsea"))
})

test_that("a layer with only ORA rows contributes ORA", {
    df <- ora_only_layer()
    expect_identical(select_layer_method_rows(df, .raw_p_values(df))$method, "ora")
})

test_that("rank-based rows with no p-value do not win", {
    df <- gene_layer()
    df$pval <- NA_real_
    sel <- select_layer_method_rows(df, .raw_p_values(df))
    expect_identical(sel$method, "ora")
})

test_that("rows of unknown method are left out once a method is chosen", {
    df <- rbind(ora_only_layer(),
                data.frame(ID = "map00040", pvalue = 1e-9, padj = 1e-8,
                           method = NA, stringsAsFactors = FALSE))
    sel <- select_layer_method_rows(df, .raw_p_values(df))
    expect_identical(sel$method, "ora")
    expect_false(sel$keep[3])
})

test_that("a table with no method column is used whole and labelled unspecified", {
    df <- data.frame(ID = c("map00010", "map00020"), pvalue = c(0.01, 0.02))
    sel <- select_layer_method_rows(df, .raw_p_values(df))
    expect_identical(sel$method, "unspecified")
    expect_true(all(sel$keep))
})


# ---- the merge ------------------------------------------------------------

test_that("the merge takes the minimum within one method, never across two", {
    tables <- list(proteomics = gene_layer(), metabolomics = ora_only_layer())

    merged <- merge_pathway_pvalues(tables, c("00010", "00020", "00030"),
                                    names(tables), kegg_org = NULL)
    row <- merged[merged$norm_id == "00010", ]

    # The ORA row for 00010 is far smaller; it must not be what the gene layer
    # contributes once that layer has fgsea rows.
    expect_equal(row$pval_proteomics, 0.04)
    expect_identical(row$method_proteomics, "fgsea")
    expect_identical(row$method_metabolomics, "ora")
})

test_that("the method column is NA where the layer has no p-value", {
    tables <- list(proteomics = gene_layer(), metabolomics = ora_only_layer())

    merged <- merge_pathway_pvalues(tables, c("00020", "00030"),
                                    names(tables), kegg_org = NULL)

    expect_true(is.na(merged$method_metabolomics[merged$norm_id == "00020"]))
    expect_true(is.na(merged$method_proteomics[merged$norm_id == "00030"]))
})

test_that("the method columns do not enter the Stouffer combination", {
    tables <- list(proteomics = gene_layer(), metabolomics = ora_only_layer())
    merged <- merge_pathway_pvalues(tables, c("00010", "00020"),
                                    names(tables), kegg_org = NULL)

    meta <- stouffer_combined_pvalues(merged)

    expect_equal(meta$n_omics[meta$norm_id == "00010"], 2)
    expect_false(anyNA(meta$combined_pval))
})


# ---- supplementary rank-based tables --------------------------------------

test_that("a rank table without a method column is refused with a warning", {
    tab <- gsea_table()
    tab$method <- NULL
    expect_warning(out <- .usable_rank_tables(list(metabolomics = tab)),
                   "metabolomics")
    expect_length(out, 0)
})

test_that("only rank-based rows of a supplementary table are kept", {
    tab <- rbind(gsea_table(),
                 transform(gsea_table()[1, ], method = "ora"))
    out <- .usable_rank_tables(list(metabolomics = tab))
    expect_identical(unique(out$metabolomics$method), "fgsea")
    expect_equal(nrow(out$metabolomics), 2)
})

test_that("a layer supplied only as rank-based rows counts as a layer", {
    res <- suppressMessages(analyze_cross_omics_enrichment(
        list(proteomics = gene_layer()), nonmodel_config(), out_dir = NULL,
        rank_tables = list(metabolomics = gsea_table())))

    expect_false(is.null(res))
    meta <- res$meta_analysis
    expect_true(all(c("pval_proteomics", "pval_metabolomics") %in% names(meta)))
    methods <- meta$method_metabolomics
    expect_identical(unique(methods[!is.na(methods)]), "fgsea")
    expect_identical(unname(res$layer_methods[c("proteomics", "metabolomics")]),
                     c("fgsea", "fgsea"))
})

test_that("rank rows reach the meta-analysis but not the per-layer tables", {
    res <- suppressMessages(analyze_cross_omics_enrichment(
        list(proteomics = gene_layer(), metabolomics = ora_only_layer()),
        nonmodel_config(), out_dir = NULL,
        rank_tables = list(metabolomics = gsea_table())))

    # pathway_tables drives the per-layer figures and CSVs: still ORA only.
    expect_identical(unique(res$pathway_tables$metabolomics$method), "ora")

    # The meta-analysis reads the layer's GSEA rows, not its smaller ORA p.
    row <- res$meta_analysis[res$meta_analysis$norm_id == "00010", ]
    expect_equal(row$pval_metabolomics, 0.03)
    expect_identical(row$method_metabolomics, "fgsea")
})

test_that("one enriched layer and no rank table is still too few", {
    expect_null(suppressMessages(analyze_cross_omics_enrichment(
        list(proteomics = gene_layer()), nonmodel_config(), out_dir = NULL)))
})
