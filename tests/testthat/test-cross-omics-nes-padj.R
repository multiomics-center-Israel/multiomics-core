# tests/testthat/test-cross-omics-nes-padj.R
#
# Unit tests for the significance-statistic selection and the signed-NES heatmap
# path in R/domain/multiomics/07_enrichment.R:
#   - resolve_pvalue_type()
#   - pick_pvalue_column()
#   - describe_pvalue_axis()
#   - merge_pathway_pvalues()      (padj vs raw p, and the per-layer fallback)
#   - build_layer_nes_matrix()
#   - plot_cross_omics_pathway_heatmap()  (NES mode vs -log10 fallback)
#   - enrich_feature_list_gmt()    (non-model GMT fallback for loadings)
#
# All data here is synthetic: made-up accessions, made-up term names and
# made-up feature ids. Nothing is taken from a real experiment.

# --- fixtures ----------------------------------------------------------------

# Two layers that both carry raw and adjusted p-values, and NES on one of them.
make_tables_with_nes <- function() {
    list(
        transcriptomics = data.frame(
            pathway      = c("GO:0000001", "GO:0000002", "GO:0000003"),
            pathway_name = c("alpha process", "beta process", "gamma process"),
            pvalue       = c(0.001, 0.010, 0.400),
            padj         = c(0.020, 0.300, 0.900),
            NES          = c(2.10, -1.60, 0.20),
            stringsAsFactors = FALSE
        ),
        proteomics = data.frame(
            pathway      = c("GO:0000001", "GO:0000002"),
            pathway_name = c("alpha process", "beta process"),
            pvalue       = c(0.005, 0.030),
            padj         = c(0.040, 0.100),
            NES          = c(1.30, -2.40),
            stringsAsFactors = FALSE
        )
    )
}

# Same layers, but the transcriptomics layer is over-representation output:
# it has no NES column at all and no usable padj.
make_tables_mixed <- function() {
    tabs <- make_tables_with_nes()
    tabs$transcriptomics$NES <- NULL
    tabs$transcriptomics$padj <- NA_real_
    tabs
}

make_meta_two_layers <- function() {
    data.frame(
        pathway              = c("GO:0000001", "GO:0000002", "GO:0000003"),
        pathway_name         = c("alpha process", "beta process", "gamma process"),
        pval_transcriptomics = c(0.020, 0.300, 0.900),
        pval_proteomics      = c(0.040, 0.100, NA),
        combined_pval        = c(0.001, 0.150, 0.900),
        n_omics              = c(2L, 2L, 1L),
        stringsAsFactors = FALSE
    )
}


# --- resolve_pvalue_type() ---------------------------------------------------

test_that("pvalue_type defaults to padj and honours an explicit choice", {
    expect_equal(resolve_pvalue_type(list()), "padj")
    expect_equal(
        resolve_pvalue_type(list(modes = list(multiomics = list(
            enrichment = list(pvalue_type = "pvalue"))))),
        "pvalue"
    )
    expect_equal(
        resolve_pvalue_type(list(modes = list(multiomics = list(
            enrichment = list(pvalue_type = " PADJ "))))),
        "padj"
    )
})

test_that("an unrecognised pvalue_type warns and falls back to padj", {
    cfg <- list(modes = list(multiomics = list(
        enrichment = list(pvalue_type = "bonferroni"))))
    expect_warning(res <- resolve_pvalue_type(cfg), "pvalue_type")
    expect_equal(res, "padj")
})


# --- pick_pvalue_column() ----------------------------------------------------

test_that("pick_pvalue_column prefers the requested statistic", {
    df <- data.frame(pvalue = c(0.01, 0.02), padj = c(0.10, 0.20))

    expect_equal(pick_pvalue_column(df, "padj"),
                 list(column = "padj", type = "padj"))
    expect_equal(pick_pvalue_column(df, "pvalue"),
                 list(column = "pvalue", type = "pvalue"))
})

test_that("pick_pvalue_column falls back and reports the statistic it used", {
    only_raw <- data.frame(pvalue = c(0.01, 0.02))
    expect_equal(pick_pvalue_column(only_raw, "padj"),
                 list(column = "pvalue", type = "pvalue"))

    # A column that exists but is entirely NA is not usable.
    all_na <- data.frame(pvalue = c(0.01, 0.02), padj = c(NA_real_, NA_real_))
    expect_equal(pick_pvalue_column(all_na, "padj"),
                 list(column = "pvalue", type = "pvalue"))

    expect_null(pick_pvalue_column(data.frame(pathway = "GO:0000001"), "padj"))
    expect_null(pick_pvalue_column(NULL, "padj"))
})

test_that("alternative adjusted-p spellings are recognised", {
    expect_equal(pick_pvalue_column(data.frame(p.adjust = 0.01), "padj")$column,
                 "p.adjust")
    expect_equal(pick_pvalue_column(data.frame(pval = 0.01), "pvalue")$column,
                 "pval")
})


# --- describe_pvalue_axis() --------------------------------------------------

test_that("the axis label names the statistic and any layer that fell back", {
    expect_equal(describe_pvalue_axis("padj"), "-log10(adjusted p)")
    expect_equal(describe_pvalue_axis("pvalue"), "-log10(p-value)")

    label <- describe_pvalue_axis(
        "padj", c(transcriptomics = "pvalue", proteomics = "padj"))
    expect_match(label, "^-log10\\(adjusted p\\)")
    expect_match(label, "raw p for: transcriptomics")
    expect_false(grepl("proteomics", label))
})


# --- merge_pathway_pvalues() -------------------------------------------------

test_that("merge_pathway_pvalues takes padj by default, raw p on request", {
    tabs <- make_tables_with_nes()
    targets <- c("GO:0000001", "GO:0000002", "GO:0000003")

    adj <- merge_pathway_pvalues(tabs, targets, names(tabs), pvalue_type = "padj")
    adj <- adj[match(targets, adj$pathway), ]
    expect_equal(adj$pval_transcriptomics, c(0.020, 0.300, 0.900))
    expect_equal(adj$pval_proteomics, c(0.040, 0.100, NA))
    expect_equal(attr(
        merge_pathway_pvalues(tabs, targets, names(tabs), pvalue_type = "padj"),
        "pvalue_type_used"),
        c(transcriptomics = "padj", proteomics = "padj"))

    raw <- merge_pathway_pvalues(tabs, targets, names(tabs), pvalue_type = "pvalue")
    raw <- raw[match(targets, raw$pathway), ]
    expect_equal(raw$pval_transcriptomics, c(0.001, 0.010, 0.400))
})

test_that("a layer without padj falls back to raw p and is flagged", {
    tabs <- make_tables_mixed()
    targets <- c("GO:0000001", "GO:0000002", "GO:0000003")

    merged <- merge_pathway_pvalues(tabs, targets, names(tabs), pvalue_type = "padj")
    used <- attr(merged, "pvalue_type_used")

    expect_equal(used[["transcriptomics"]], "pvalue")
    expect_equal(used[["proteomics"]], "padj")

    merged <- merged[match(targets, merged$pathway), ]
    expect_equal(merged$pval_transcriptomics, c(0.001, 0.010, 0.400))
    expect_equal(merged$pval_proteomics, c(0.040, 0.100, NA))
})

test_that("pathways absent from a layer stay NA rather than becoming zero", {
    tabs <- make_tables_with_nes()
    merged <- merge_pathway_pvalues(tabs, c("GO:0000003"), names(tabs))
    expect_true(is.na(merged$pval_proteomics))
})


# --- build_layer_nes_matrix() ------------------------------------------------

test_that("build_layer_nes_matrix aligns NES to the requested rows and columns", {
    tabs <- make_tables_with_nes()
    pathways <- c("GO:0000001", "GO:0000002", "GO:0000003")

    nes <- build_layer_nes_matrix(tabs, pathways, names(tabs))

    expect_equal(dim(nes), c(3L, 2L))
    expect_equal(rownames(nes), pathways)
    expect_equal(colnames(nes), c("transcriptomics", "proteomics"))
    expect_equal(nes[, "transcriptomics"],
                 c(`GO:0000001` = 2.10, `GO:0000002` = -1.60, `GO:0000003` = 0.20))
    # Not scored in proteomics -> NA, never 0.
    expect_true(is.na(nes["GO:0000003", "proteomics"]))
})

test_that("duplicate rows resolve to the most significant one", {
    tabs <- list(proteomics = data.frame(
        pathway = c("GO:0000001", "GO:0000001"),
        pvalue  = c(0.400, 0.001),
        padj    = c(0.500, 0.010),
        NES     = c(0.30, -2.80),
        stringsAsFactors = FALSE
    ))
    nes <- build_layer_nes_matrix(tabs, "GO:0000001", "proteomics")
    expect_equal(unname(nes[1, 1]), -2.80)
})

test_that("build_layer_nes_matrix returns NULL when no layer reports NES", {
    tabs <- make_tables_with_nes()
    tabs$transcriptomics$NES <- NULL
    tabs$proteomics$NES <- NULL
    expect_null(build_layer_nes_matrix(tabs, "GO:0000001", names(tabs)))
    expect_null(build_layer_nes_matrix(NULL, "GO:0000001", "proteomics"))
})


# --- plot_cross_omics_pathway_heatmap() --------------------------------------

# The plotters draw on the active device, so these tests only assert that the
# right branch is taken and that nothing errors; the colours themselves are
# checked by eye on the rendered PNG.

draw_heatmap_to_temp <- function(...) {
    f <- tempfile(fileext = ".png")
    png(f, width = 900, height = 700, res = 110)
    on.exit(dev.off(), add = TRUE)
    plot_cross_omics_pathway_heatmap(...)
    f
}

test_that("the heatmap uses signed NES when every layer reports one", {
    expect_silent(f <- draw_heatmap_to_temp(
        make_meta_two_layers(), c("transcriptomics", "proteomics"),
        pathway_tables = make_tables_with_nes(),
        value_label = "-log10(adjusted p)"
    ))
    expect_true(file.exists(f))
    expect_gt(file.size(f), 0)
})

test_that("a layer without NES keeps its column but is marked, not silently grey", {
    meta <- make_meta_two_layers()
    nes <- build_layer_nes_matrix(make_tables_mixed(), meta$pathway,
                                  c("transcriptomics", "proteomics"))

    # NES still available from the other layer, so the figure stays on the NES
    # scale; the ORA layer is all-NA and gets a "(no NES)" column header.
    expect_false(is.null(nes))
    expect_true(all(is.na(nes[, "transcriptomics"])))
    expect_false(all(is.na(nes[, "proteomics"])))

    expect_silent(f <- draw_heatmap_to_temp(
        meta, c("transcriptomics", "proteomics"),
        pathway_tables = make_tables_mixed()
    ))
    expect_true(file.exists(f))
})

test_that("the heatmap falls back to -log10 when no layer reports NES", {
    tabs <- make_tables_with_nes()
    tabs$transcriptomics$NES <- NULL
    tabs$proteomics$NES <- NULL

    expect_null(build_layer_nes_matrix(tabs, make_meta_two_layers()$pathway,
                                       names(tabs)))
    expect_silent(f <- draw_heatmap_to_temp(
        make_meta_two_layers(), c("transcriptomics", "proteomics"),
        pathway_tables = tabs
    ))
    expect_true(file.exists(f))
})

test_that("row selection and labels survive the NES change", {
    meta <- make_meta_two_layers()
    sel <- select_multi_omics_pathways(meta, top_n = 3)

    # Two layers below 0.05 outrank one, whatever the combined p says.
    expect_equal(sel$pathway[1], "GO:0000001")
    expect_equal(
        format_pathway_labels(sel$pathway, sel$pathway_name, append_id = TRUE)[1],
        "alpha process (GO:0000001)"
    )
})


# --- enrich_feature_list_gmt() -----------------------------------------------

write_synthetic_gmt <- function(path) {
    lines <- c(
        paste(c("GO:0000001", "alpha process",
                paste0("gene", 1:12)), collapse = "\t"),
        paste(c("GO:0000002", "beta process",
                paste0("gene", 13:30)), collapse = "\t")
    )
    writeLines(lines, path)
    path
}

make_gmt_config <- function(gmt_path) {
    list(
        project = list(dir = dirname(gmt_path)),
        paths = list(raw = "."),
        modes = list(rna = list(pathway = list(gmt_file = gmt_path)))
    )
}

test_that("the GMT fallback runs ORA on a loadings feature list", {
    skip_if_not_installed("clusterProfiler")

    gmt <- write_synthetic_gmt(tempfile(fileext = ".gmt"))
    harm <- list(inputs = list(transcriptomics = list(
        expr_work = matrix(0, nrow = 30, ncol = 2,
                           dimnames = list(paste0("gene", 1:30), c("S1", "S2")))
    )))

    res <- suppressMessages(enrich_feature_list_gmt(
        resolved_ids = paste0("gene", 1:11),   # 11 of the 12 alpha members
        omics_type = "transcriptomics",
        harmonization_res = harm,
        config = make_gmt_config(gmt)
    ))

    expect_s3_class(res, "data.frame")
    expect_true(all(c("pathway", "ID", "pvalue", "padj", "GeneRatio", "setSize")
                    %in% names(res)))
    expect_true("GO:0000001" %in% res$ID)
})

test_that("the GMT fallback declines cleanly when it has nothing to work with", {
    harm <- list(inputs = list(transcriptomics = list(expr_work = NULL)))

    expect_null(enrich_feature_list_gmt("gene1", "transcriptomics", harm, NULL))
    expect_null(enrich_feature_list_gmt(
        "gene1", "metabolomics", harm,
        make_gmt_config(tempfile(fileext = ".gmt"))))
    expect_null(suppressMessages(enrich_feature_list_gmt(
        "gene1", "transcriptomics", harm,
        make_gmt_config(file.path(tempdir(), "does-not-exist.gmt")))))
})


# --- resolve_gene_n_ids() ----------------------------------------------------

test_that("MOFA's view-suffixed GENE_N ids resolve like the bare ones", {
    harm <- list(gene_protein_mapping = data.frame(
        gene_id    = c("geneA", "geneB", "geneC"),
        protein_id = c("protA", "protB", "protC"),
        stringsAsFactors = FALSE
    ))

    expect_equal(
        suppressMessages(resolve_gene_n_ids(
            c("GENE_1", "GENE_2_transcriptomics", "GENE_3"),
            harm, "transcriptomics")),
        c("geneA", "geneB", "geneC")
    )
    expect_equal(
        suppressMessages(resolve_gene_n_ids("GENE_2_proteomics", harm, "proteomics")),
        "protB"
    )
    # Native ids pass straight through.
    expect_equal(resolve_gene_n_ids(c("geneA", "geneB"), harm, "transcriptomics"),
                 c("geneA", "geneB"))
})
