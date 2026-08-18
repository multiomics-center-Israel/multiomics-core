# tests/testthat/test-cross-omics-ora-heatmap.R
#
# Unit tests for the over-representation companion to the signed-NES heatmap in
# R/domain/multiomics/07_enrichment.R:
#   - pathway_id_column()
#   - filter_ora_rows()
#   - build_layer_ora_matrix()
#   - plot_cross_omics_ora_heatmap()
#
# All data here is synthetic: made-up accessions, made-up term names and
# made-up significance values. Nothing is taken from a real experiment.

# --- fixtures ----------------------------------------------------------------

# A layer that ran both fgsea and ORA, with the two ORA directions of the same
# pathway on separate rows — the shape every gene-based layer arrives in.
make_mixed_method_table <- function() {
    data.frame(
        pathway      = c("GO:0000001", "GO:0000002",
                         "GO:0000001", "GO:0000001", "GO:0000002"),
        pathway_name = c("alpha process", "beta process",
                         "alpha process", "alpha process", "beta process"),
        method       = c("fgsea", "fgsea", "ora", "ora", "ora"),
        direction    = c(NA, NA, "up", "down", "up"),
        NES          = c(2.10, -1.60, NA, NA, NA),
        pvalue       = c(0.001, 0.010, 0.004, 0.300, 0.200),
        padj         = c(0.020, 0.100, 0.010, 0.600, 0.400),
        stringsAsFactors = FALSE
    )
}

# A compound-ORA layer: no `method` column, no NES, just accessions and p-values.
make_compound_ora_table <- function() {
    data.frame(
        pathway      = c("map00001", "GO:0000001"),
        pathway_name = c("synthetic map", "alpha process"),
        pvalue       = c(0.030, 0.070),
        padj         = c(0.250, 0.400),
        stringsAsFactors = FALSE
    )
}

make_ora_tables <- function() {
    list(
        transcriptomics = filter_ora_rows(make_mixed_method_table()),
        proteomics      = filter_ora_rows(make_compound_ora_table())
    )
}


# --- pathway_id_column() -----------------------------------------------------

test_that("pathway_id_column prefers pathway, then ID, then Description", {
    expect_equal(pathway_id_column(data.frame(pathway = "GO:0000001",
                                              ID = "GO:0000001")), "pathway")
    expect_equal(pathway_id_column(data.frame(ID = "GO:0000001",
                                              Description = "alpha")), "ID")
    expect_equal(pathway_id_column(data.frame(Description = "alpha")), "Description")
    expect_null(pathway_id_column(data.frame(padj = 0.01)))
    expect_null(pathway_id_column(NULL))
})


# --- filter_ora_rows() -------------------------------------------------------

test_that("filter_ora_rows keeps the ORA rows of a mixed-method table", {
    ora <- filter_ora_rows(make_mixed_method_table())

    expect_equal(nrow(ora), 3L)
    expect_true(all(ora$method == "ora"))
    expect_true(all(is.na(ora$NES)))
})

test_that("the method label is matched case- and whitespace-insensitively", {
    df <- make_mixed_method_table()
    df$method <- c("fgsea", "FGSEA", " ORA ", "Ora", "ora")
    expect_equal(nrow(filter_ora_rows(df)), 3L)
})

test_that("a table with neither method nor NES is taken as over-representation", {
    # Compound ORA writes no method label; dropping it on that would hide the
    # only layer this figure exists for.
    expect_equal(nrow(filter_ora_rows(make_compound_ora_table())), 2L)
})

test_that("a GSEA-only table without a method label contributes nothing", {
    gsea_only <- data.frame(pathway = "GO:0000001", NES = 1.4, padj = 0.01)
    expect_equal(nrow(filter_ora_rows(gsea_only)), 0L)
})

test_that("filter_ora_rows passes NULL and empty tables straight through", {
    expect_null(filter_ora_rows(NULL))
    expect_equal(nrow(filter_ora_rows(data.frame(pathway = character(0)))), 0L)
})


# --- build_layer_ora_matrix() ------------------------------------------------

test_that("build_layer_ora_matrix aligns FDRs to the requested rows and columns", {
    tabs <- make_ora_tables()
    pathways <- c("GO:0000001", "GO:0000002", "map00001")

    pmat <- build_layer_ora_matrix(tabs, pathways, names(tabs))

    expect_equal(dim(pmat), c(3L, 2L))
    expect_equal(rownames(pmat), pathways)
    expect_equal(colnames(pmat), c("transcriptomics", "proteomics"))
    # Two ORA directions for GO:0000001; the more significant one wins.
    expect_equal(unname(pmat["GO:0000001", "transcriptomics"]), 0.010)
    expect_equal(unname(pmat["GO:0000002", "transcriptomics"]), 0.400)
    expect_equal(unname(pmat["map00001", "proteomics"]), 0.250)
})

test_that("a pathway a layer never tested stays NA rather than becoming zero", {
    tabs <- make_ora_tables()
    pmat <- build_layer_ora_matrix(tabs, c("GO:0000002", "map00001"), names(tabs))

    expect_true(is.na(pmat["GO:0000002", "proteomics"]))
    expect_true(is.na(pmat["map00001", "transcriptomics"]))
    expect_false(any(pmat == 0, na.rm = TRUE))
})

test_that("a layer that ran no ORA still gets a column, all NA", {
    tabs <- make_ora_tables()
    tabs$metabolomics <- data.frame(pathway = character(0), padj = numeric(0))

    pmat <- build_layer_ora_matrix(tabs, "GO:0000001", names(tabs))

    expect_true("metabolomics" %in% colnames(pmat))
    expect_true(all(is.na(pmat[, "metabolomics"])))
    expect_true(is.na(attr(pmat, "pvalue_type_used")[["metabolomics"]]))
})

test_that("the matrix reports which statistic each layer supplied", {
    tabs <- make_ora_tables()
    tabs$transcriptomics$padj <- NA_real_

    pmat <- build_layer_ora_matrix(tabs, "GO:0000001", names(tabs))
    used <- attr(pmat, "pvalue_type_used")

    expect_equal(used[["transcriptomics"]], "pvalue")
    expect_equal(used[["proteomics"]], "padj")
    expect_equal(unname(pmat["GO:0000001", "transcriptomics"]), 0.004)
})

test_that("raw p can be requested instead of the FDR", {
    tabs <- make_ora_tables()
    pmat <- build_layer_ora_matrix(tabs, "GO:0000001", names(tabs),
                                   pvalue_type = "pvalue")
    expect_equal(unname(pmat["GO:0000001", "transcriptomics"]), 0.004)
})


# --- plot_cross_omics_ora_heatmap() ------------------------------------------

# The plotters draw on the active device, so these tests assert that the right
# branch is taken and that nothing errors; the colours themselves are checked by
# eye on the rendered PNG.

draw_ora_heatmap_to_temp <- function(...) {
    f <- tempfile(fileext = ".png")
    png(f, width = 900, height = 700, res = 110)
    on.exit(dev.off(), add = TRUE)
    suppressMessages(plot_cross_omics_ora_heatmap(...))
    f
}

make_meta_two_layers <- function() {
    data.frame(
        pathway              = c("GO:0000001", "GO:0000002", "map00001"),
        pathway_name         = c("alpha process", "beta process", "synthetic map"),
        pval_transcriptomics = c(0.010, 0.400, NA),
        pval_proteomics      = c(0.400, NA, 0.250),
        combined_pval        = c(0.020, 0.400, 0.250),
        n_omics              = c(2L, 1L, 1L),
        stringsAsFactors = FALSE
    )
}

test_that("the ORA heatmap draws from the over-representation rows", {
    tabs <- list(transcriptomics = make_mixed_method_table(),
                 proteomics      = make_compound_ora_table())

    f <- draw_ora_heatmap_to_temp(make_meta_two_layers(), names(tabs),
                                  pathway_tables = tabs)
    expect_true(file.exists(f))
    expect_gt(file.size(f), 0)
})

test_that("a layer that ran no ORA keeps its column instead of vanishing", {
    tabs <- list(transcriptomics = make_mixed_method_table(),
                 proteomics      = make_compound_ora_table(),
                 metabolomics    = data.frame(pathway = "GO:0000009",
                                              NES = 1.1, padj = 0.02))

    pmat <- build_layer_ora_matrix(lapply(tabs, filter_ora_rows),
                                   c("GO:0000001", "map00001"), names(tabs))
    expect_true(all(is.na(pmat[, "metabolomics"])))

    f <- draw_ora_heatmap_to_temp(make_meta_two_layers(), names(tabs),
                                  pathway_tables = tabs)
    expect_true(file.exists(f))
})

test_that("the heatmap declines cleanly when no layer ran ORA", {
    gsea_only <- list(
        transcriptomics = data.frame(pathway = "GO:0000001", NES = 1.4, padj = 0.01),
        proteomics      = data.frame(pathway = "GO:0000002", NES = -1.2, padj = 0.03)
    )

    expect_silent(f <- draw_ora_heatmap_to_temp(make_meta_two_layers(),
                                                names(gsea_only),
                                                pathway_tables = gsea_only))
    expect_true(file.exists(f))

    expect_silent(f2 <- draw_ora_heatmap_to_temp(make_meta_two_layers(),
                                                 c("transcriptomics", "proteomics"),
                                                 pathway_tables = NULL))
    expect_true(file.exists(f2))
})

test_that("rows come from the ORA evidence, not from GSEA-only pathways", {
    # GO:0000002 is scored by fgsea in the meta table but its ORA rows exist too,
    # while GO:0000099 is fgsea-only and must not reach the figure.
    tabs <- list(transcriptomics = rbind(
        make_mixed_method_table(),
        data.frame(pathway = "GO:0000099", pathway_name = "omega process",
                   method = "fgsea", direction = NA, NES = 3.0,
                   pvalue = 1e-06, padj = 1e-04, stringsAsFactors = FALSE)
    ))

    ora <- lapply(tabs, filter_ora_rows)
    expect_false("GO:0000099" %in% ora$transcriptomics$pathway)

    pmat <- build_layer_ora_matrix(ora, c("GO:0000001", "GO:0000099"),
                                   names(tabs))
    expect_true(is.na(pmat["GO:0000099", "transcriptomics"]))
})
