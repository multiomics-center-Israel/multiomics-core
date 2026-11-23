# tests/testthat/test-metab-pathway-volcano.R
#
# The report's pathway overlay: a volcano whose dropdown colors the measured
# metabolites of one significant GSEA or QEA pathway. Checks the pathway ->
# feature mapping, the selector order (GSEA by |NES|, then QEA by FDR), and
# that an empty selection gives a dropdown with nothing to choose.
# All data here are synthetic.

synthetic_pre <- function() {
    ids <- paste0("F", 1:5)
    list(
        row_data = data.frame(
            feature_id = ids,
            Name       = c("cpd one", "cpd two", "cpd three", "cpd four", "cpd two isomer"),
            KEGG       = c("C00031", "C00022", "C00158", "", "C00022"),
            stringsAsFactors = FALSE
        ),
        expr_raw = matrix(1:20, nrow = 5, dimnames = list(ids, paste0("S", 1:4)))
    )
}

synthetic_gmt_config <- function() {
    gmt <- tempfile(fileext = ".gmt")
    writeLines(c(
        "map00010\tSynthetic set A\tC00031\tC00022",
        "map00020\tSynthetic set B\tC00022\tC00158\tC99999",
        "map00030\tSynthetic set C\tC88888\tC77777"
    ), gmt)
    list(
        project = list(dir = tempdir()),
        modes = list(metabolomics = list(enrichment = list(gmt_file = gmt)))
    )
}

test_that("pathway members map compounds to every feature that carries them", {
    res <- suppressMessages(
        build_pathway_feature_members(synthetic_pre(), synthetic_gmt_config())
    )

    expect_setequal(names(res$members),
                    c("map00010 - Synthetic set A", "map00020 - Synthetic set B"))
    expect_setequal(res$members[["map00010 - Synthetic set A"]], c("F1", "F2", "F5"))
    expect_setequal(res$members[["map00020 - Synthetic set B"]], c("F2", "F3", "F5"))
    expect_identical(unname(res$feature_compound["F1"]), "C00031")
    expect_false("F4" %in% names(res$feature_compound))
})

test_that("pathway members are empty without a GMT file", {
    cfg <- list(project = list(dir = tempdir()),
                modes = list(metabolomics = list(enrichment = list())))
    res <- build_pathway_feature_members(synthetic_pre(), cfg)
    expect_length(res$members, 0)
})

test_that("selector lists GSEA by |NES| then QEA by FDR, significant only", {
    members <- list(ga = "F1", gb = "F2", gc = "F3", qa = "F1", qb = "F2")
    gsea <- data.frame(
        pathway = c("ga", "gb", "gc", "gd"),
        FDR     = c(0.01, 0.04, 0.20, 0.03),
        NES     = c(1.5, -2.8, 2.0, 3.0),
        leading_edge = c("C00031", "C00022", "", ""),
        stringsAsFactors = FALSE
    )
    qea <- data.frame(
        pathway = c("qa", "qb", "qc"),
        FDR     = c(0.03, 0.001, 0.5),
        stringsAsFactors = FALSE
    )
    sel <- select_overlay_pathways(gsea, qea, members, fdr_cutoff = 0.05)

    # gc fails the FDR cutoff; gd and qc have no measured features
    expect_identical(sel$pathway, c("gb", "ga", "qb", "qa"))
    expect_identical(sel$method,  c("GSEA", "GSEA", "QEA", "QEA"))
    expect_true(all(grepl("^(GSEA|QEA) \\| ", sel$label)))
})

test_that("selector is empty when no pathway is significant", {
    gsea <- data.frame(pathway = "ga", FDR = 0.2, NES = 2.5, stringsAsFactors = FALSE)
    sel <- select_overlay_pathways(gsea, NULL, list(ga = "F1"))
    expect_equal(nrow(sel), 0)
    expect_true(all(c("method", "pathway", "FDR", "NES", "label") %in% colnames(sel)))
})

synthetic_de <- function() {
    data.frame(
        feature_id = paste0("F", 1:5),
        Name       = c("cpd one", "cpd two", "cpd three", "cpd four", "cpd two isomer"),
        logFC      = c(1.2, -0.8, 0.3, 2.1, -1.5),
        P.Value    = c(0.001, 0.02, 0.4, 0.0005, 0.01),
        adj.P.Val  = c(0.01, 0.05, 0.5, 0.008, 0.04),
        stringsAsFactors = FALSE
    )
}

test_that("overlay adds one hidden trace per pathway and a restyle button each", {
    skip_if_not_installed("plotly")
    members <- list(ga = c("F1", "F2"), qa = "F3")
    sel <- data.frame(
        method = c("GSEA", "QEA"), pathway = c("ga", "qa"), FDR = c(0.01, 0.02),
        NES = c(2.2, NA), leading_edge = c("C00031", NA), n_measured = c(2L, 1L),
        label = c("GSEA | ga", "QEA | qa"), stringsAsFactors = FALSE
    )
    fig <- plot_volcano_pathway_overlay(
        synthetic_de(), sel, members,
        feature_compound = c(F1 = "C00031", F2 = "C00022", F3 = "C00158")
    )
    built <- plotly::plotly_build(fig)$x

    expect_length(built$data, 3)
    buttons <- built$layout$updatemenus[[1]]$buttons
    expect_length(buttons, 2)
    expect_identical(vapply(buttons, function(b) b$method, character(1)),
                     c("restyle", "restyle"))
    expect_identical(unlist(buttons[[2]]$args[[2]]), c(TRUE, FALSE, TRUE))
    expect_false(isTRUE(built$data[[3]]$visible))
})

test_that("overlay with no significant pathway offers nothing to select", {
    skip_if_not_installed("plotly")
    fig <- plot_volcano_pathway_overlay(
        synthetic_de(), select_overlay_pathways(members = list()), list()
    )
    built <- plotly::plotly_build(fig)$x

    expect_length(built$data, 1)
    buttons <- built$layout$updatemenus[[1]]$buttons
    expect_length(buttons, 1)
    expect_identical(buttons[[1]]$method, "skip")
})
