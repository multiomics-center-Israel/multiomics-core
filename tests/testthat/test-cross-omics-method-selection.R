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

# fgsea rows for 00010 and 00020, plus an ORA-only pathway 00099 that the
# fgsea rows never scored.
gene_layer_with_ora_only_pathway <- function() {
    dplyr::bind_rows(
        gene_layer(),
        data.frame(ID = "map00099", pvalue = 1e-4, padj = 1e-3, method = "ora",
                   stringsAsFactors = FALSE))
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


# ---- candidates come from contributing rows only --------------------------

test_that("a row contributes only when it belongs to the chosen method", {
    cb <- .layer_contribution(gene_layer_with_ora_only_pathway())
    expect_identical(cb$method, "fgsea")
    df <- gene_layer_with_ora_only_pathway()
    expect_identical(cb$keep, df$method == "fgsea")
})

test_that("a pathway only an unselected method supports is not a candidate", {
    res <- suppressMessages(analyze_cross_omics_enrichment(
        list(proteomics = gene_layer_with_ora_only_pathway(),
             metabolomics = ora_only_layer()),
        nonmodel_config(), out_dir = NULL))

    # 00099 exists only in the gene layer's ORA rows, and that layer's selected
    # method is fgsea: no layer can contribute it.
    expect_false("00099" %in% res$union_pathways)
    expect_false("00099" %in% res$meta_analysis$norm_id)
    expect_true(all(res$meta_analysis$n_omics >= 1))
})

test_that("that pathway is a candidate once another layer contributes it", {
    metab <- rbind(ora_only_layer(),
                   data.frame(ID = "map00099", pvalue = 0.03, padj = 0.1,
                              method = "ora", stringsAsFactors = FALSE))
    res <- suppressMessages(analyze_cross_omics_enrichment(
        list(proteomics = gene_layer_with_ora_only_pathway(), metabolomics = metab),
        nonmodel_config(), out_dir = NULL))

    row <- res$meta_analysis[res$meta_analysis$norm_id == "00099", ]
    expect_equal(nrow(row), 1)
    # Only the metabolomics layer is behind it; the gene layer's ORA row is not.
    expect_true(is.na(row$pval_proteomics))
    expect_equal(row$pval_metabolomics, 0.03)
    expect_equal(row$n_omics, 1)
})


# ---- disclosure: what choosing one method leaves out -----------------------

# Shaped like run_kegg_enrichment_for_omics() output: it picks GSEA or ORA per
# contrast on the mapped-feature count, so one layer can hold GSEA for one
# contrast and ORA for another, bound with .rbind_fill().
kegg_from_de_layer <- function() {
    .rbind_fill(list(
        data.frame(pathway = c("Glycolysis", "Citrate cycle"),
                   ID = c("map00010", "map00020"), pvalue = c(0.02, 0.3),
                   padj = c(0.1, 0.5), NES = c(1.7, -0.6), setSize = c(40L, 30L),
                   method = "gsea", contrast = "A_vs_B", omics = "proteomics",
                   stringsAsFactors = FALSE),
        data.frame(pathway = c("Glycolysis", "Pentose phosphate"),
                   ID = c("map00010", "map00030"), pvalue = c(1e-5, 0.001),
                   padj = c(1e-4, 0.01), Fold_enrichment = c(3.1, 2.4),
                   Count = c(8L, 5L), method = "ora", contrast = "C_vs_D",
                   omics = "proteomics", stringsAsFactors = FALSE)
    ))
}

test_that("a layer scored by GSEA for one contrast and ORA for another contributes GSEA", {
    tables <- list(proteomics = kegg_from_de_layer(), metabolomics = ora_only_layer())
    merged <- merge_pathway_pvalues(tables, c("00010", "00020", "00030"),
                                    names(tables), kegg_org = NULL)

    # The ORA contrast's far smaller p for 00010 is not what the layer gives.
    expect_equal(merged$pval_proteomics[merged$norm_id == "00010"], 0.02)
    expect_true(all(merged$method_proteomics[!is.na(merged$pval_proteomics)] == "gsea"))
    # 00030 only the ORA contrast scored: not this layer's evidence at run level.
    expect_true(is.na(merged$pval_proteomics[merged$norm_id == "00030"]))

    lost <- .method_selection_losses(kegg_from_de_layer(),
                                     .layer_contribution(kegg_from_de_layer()))
    expect_equal(lost$n_rows, 2)
    expect_identical(lost$contrasts, "C_vs_D")
})

test_that("the run log names the contrast the method choice left out", {
    expect_message(
        analyze_cross_omics_enrichment(
            list(proteomics = kegg_from_de_layer(), metabolomics = ora_only_layer()),
            nonmodel_config(), out_dir = NULL),
        "proteomics: the meta-analysis uses its gsea rows only; 2 row\\(s\\) of other methods are left out, including every row for contrast\\(s\\) C_vs_D")
})

test_that("a collection scored only by the other method is named too", {
    # As run_pathway_analysis() writes with method "both" when fgsea skipped a
    # collection but ORA did not.
    df <- dplyr::bind_rows(
        data.frame(pathway = "map00010", pval = 0.04, padj = 0.1, NES = 1.2,
                   method = "fgsea", database = "KEGG", contrast = "A_vs_B",
                   stringsAsFactors = FALSE),
        data.frame(ID = "GO:0006096", pvalue = 0.001, padj = 0.01,
                   method = "ora", database = "GO_BP", contrast = "A_vs_B",
                   stringsAsFactors = FALSE))
    lost <- .method_selection_losses(df, .layer_contribution(df))
    expect_identical(lost$collections, "GO_BP")
    # Its contrast still has fgsea rows, so no contrast is reported lost.
    expect_identical(lost$contrasts, character(0))

    expect_message(
        analyze_cross_omics_enrichment(
            list(proteomics = df, metabolomics = ora_only_layer()),
            nonmodel_config(), out_dir = NULL),
        "including every row for collection\\(s\\) GO_BP")
})

test_that("nothing is reported when no rows of another method are left out", {
    lost <- .method_selection_losses(ora_only_layer(),
                                     .layer_contribution(ora_only_layer()))
    expect_equal(lost$n_rows, 0)
})


# ---- per-contrast discovery reads the rank-based tables too ----------------

test_that("a contrast only a rank-based source scored is discovered", {
    po <- list(proteomics = data.frame(ID = "map00010", pvalue = 0.01,
                                       method = "ora", contrast = "A_vs_B",
                                       stringsAsFactors = FALSE))
    rk <- list(metabolomics = data.frame(ID = c("map00010", "map00010"),
                                         pval = c(0.02, 0.03), method = "fgsea",
                                         contrast = c("A vs B", "C_vs_D"),
                                         stringsAsFactors = FALSE))
    out <- canonicalize_enrichment_contrasts(po, rk)

    # One canonical name per contrast key, per_omics spelling first.
    expect_identical(out$contrast_names, c("A_vs_B", "C_vs_D"))
    expect_identical(out$per_omics$proteomics$contrast, "A_vs_B")
})

test_that("an unusable rank table does not add contrasts", {
    po <- list(proteomics = data.frame(ID = "map00010", pvalue = 0.01,
                                       method = "ora", contrast = "A_vs_B",
                                       stringsAsFactors = FALSE))
    rk <- list(metabolomics = data.frame(ID = "map00010", pval = 0.02,
                                         method = "ora", contrast = "C_vs_D",
                                         stringsAsFactors = FALSE))
    expect_identical(canonicalize_enrichment_contrasts(po, rk)$contrast_names,
                     "A_vs_B")
})

test_that("the per-contrast orchestration discovers contrasts with the rank tables", {
    src <- paste(deparse(body(mod_multiomics_enrichment)), collapse = " ")
    expect_true(grepl("canonicalize_enrichment_contrasts(per_omics, rank_tables)",
                      src, fixed = TRUE))
})
