# What a coloured node on a KEGG map is allowed to mean.
#
# Colouring every measured feature makes a map show coverage rather than
# signal: a saturated picture that says only "these genes were assayed". The
# rule is therefore AND, not OR -- a feature must have moved by more than the
# fold-change threshold AND have a p-value supporting it before it may
# contribute to a node.
#
# The three thresholds live in one list, `.PATHVIEW_THRESHOLDS`, because three
# things have to agree about them: which pathways are selected
# (.kegg_hits_by_contrast), which features may colour a node
# (filter_changed_features), and what the caption under the figure claims
# (pathview_significance_caption). Kept apart, a caption could describe a rule
# the filter was not applying, and nothing would ever catch it.
#
# All fixtures synthetic; no pathview call, no network.

changed_fixture <- function() {
    data.frame(
        feature_id = c("big_sig", "down_big_sig", "big_ns", "small_sig",
                       "small_ns", "borderline_sig"),
        log2fc = c(1.5, -1.5, 1.5, 0.1, 0.1, log2(1.5)),
        pvalue = c(0.01, 0.01, 0.90, 0.01, 0.90, 0.01),
        stringsAsFactors = FALSE
    )
}


# ---- the node-evidence rule -------------------------------------------------

test_that("a node needs both a large change and a supported one", {
    kept <- suppressMessages(filter_changed_features(changed_fixture()))

    # A big but unsupported change is usually a noisy low-abundance feature; a
    # confident but tiny one is not what a pathway map exists to highlight.
    expect_setequal(kept$feature_id, c("big_sig", "down_big_sig"))
})

test_that("the fold-change rule is strict, so sitting on the threshold is not enough", {
    # borderline_sig is exactly log2(1.5) and is significant. The comparison is
    # ">", so significance does not rescue it. Pinned because a boundary is the
    # part most easily changed by accident.
    kept <- suppressMessages(filter_changed_features(changed_fixture()))

    expect_false("borderline_sig" %in% kept$feature_id)
})

test_that("a feature with no usable fold change stays off the map", {
    # NaN is what log2() of a signed linear fold change leaves behind for every
    # down-regulated feature. Such a node must stay uncoloured rather than be
    # coloured on its p-value alone.
    de <- data.frame(
        feature_id = c("nan_fc", "na_fc", "inf_fc", "clean"),
        log2fc = c(NaN, NA_real_, Inf, 2),
        pvalue = rep(0.001, 4),
        stringsAsFactors = FALSE
    )

    expect_equal(suppressMessages(filter_changed_features(de))$feature_id, "clean")
})

test_that("an adjusted p column is used only where there is no raw one, and it says so", {
    # Fold changes clear the bar on their own, so this isolates which p column
    # the filter picked up.
    de <- data.frame(feature_id = c("a", "b"), log2fc = c(2, 2),
                     adj.P.Val = c(0.01, 0.90), stringsAsFactors = FALSE)

    expect_message(kept <- filter_changed_features(de), "adj\\.P\\.Val")
    expect_equal(kept$feature_id, "a")

    # With both present the raw p decides, and nothing is announced.
    both <- data.frame(feature_id = c("a", "b"), log2fc = c(2, 2),
                       pvalue = c(0.90, 0.01), adj.P.Val = c(0.01, 0.90),
                       stringsAsFactors = FALSE)
    expect_silent(kept_both <- filter_changed_features(both))
    expect_equal(kept_both$feature_id, "b")
})

test_that("a table the rule cannot be applied to colours nothing, rather than everything", {
    no_p <- data.frame(feature_id = "a", log2fc = 2, stringsAsFactors = FALSE)
    no_fc <- data.frame(feature_id = "a", pvalue = 0.001, stringsAsFactors = FALSE)

    expect_message(expect_equal(nrow(filter_changed_features(no_p)), 0L),
                   "no p-value column")
    expect_message(expect_equal(nrow(filter_changed_features(no_fc)), 0L),
                   "no fold-change column")
})

test_that("empty and absent tables pass through untouched", {
    empty <- data.frame(feature_id = character(0), log2fc = numeric(0),
                        pvalue = numeric(0), stringsAsFactors = FALSE)

    expect_equal(nrow(suppressMessages(filter_changed_features(empty))), 0L)
    expect_null(filter_changed_features(NULL))
})

test_that("the thresholds can be overridden as a whole", {
    # log2fc 0.3 is about a 1.23-fold change: under the default 1.5-fold bar,
    # over a 1.1-fold one.
    de <- data.frame(feature_id = c("a", "b"), log2fc = c(0.3, 0.3),
                     pvalue = c(0.02, 0.2), stringsAsFactors = FALSE)

    expect_equal(nrow(suppressMessages(filter_changed_features(de))), 0L)
    # Loosening the fold-change bar admits only the supported feature; the
    # p-value rule still applies.
    loose <- list(fdr_alpha = 0.05, node_fc = 1.1, node_p = 0.05)
    expect_equal(suppressMessages(filter_changed_features(de, loose))$feature_id,
                 "a")
})


# ---- the caption and the rules cannot drift ---------------------------------

test_that("the caption states the thresholds it was given", {
    cap <- pathview_significance_caption()

    expect_true(grepl("scored them below 0.05", cap, fixed = TRUE))
    expect_true(grepl("1.5-fold", cap, fixed = TRUE))
    expect_true(grepl("0.58", cap, fixed = TRUE))
    expect_true(grepl("raw p < 0.05", cap, fixed = TRUE))
})

test_that("the caption follows the thresholds rather than repeating constants", {
    cap <- pathview_significance_caption(
        list(fdr_alpha = 0.1, node_fc = 2, node_p = 0.01))

    expect_true(grepl("scored them below 0.1", cap, fixed = TRUE))
    expect_true(grepl("2-fold", cap, fixed = TRUE))
    expect_true(grepl("raw p < 0.01", cap, fixed = TRUE))
    expect_false(grepl("1.5-fold", cap, fixed = TRUE))
})

test_that("the caption does not promise a distinction the map cannot draw", {
    # pathview draws "measured but unchanged" and "never measured" identically.
    # A caption implying otherwise would invite the reader to read absence of
    # colour as evidence of absence.
    cap <- pathview_significance_caption()

    expect_true(grepl("no measured feature", cap, fixed = TRUE))
    expect_true(grepl("does not separate the two", cap, fixed = TRUE))
})

test_that("the pathway selector and the caption read the same cutoff", {
    # The whole reason the thresholds are one list. If a future change moves the
    # selector's default without the caption, this fails rather than shipping a
    # figure whose legend describes a rule it did not apply.
    expect_identical(deparse(formals(.kegg_hits_by_contrast)$alpha),
                     ".PATHVIEW_THRESHOLDS$fdr_alpha")
    expect_true(grepl(paste0("scored them below ",
                             .PATHVIEW_THRESHOLDS$fdr_alpha),
                      pathview_significance_caption(), fixed = TRUE))
})


# ---- both renderers, one caption --------------------------------------------

test_that("the caption does not claim an FDR where a selector fell back to raw p", {
    # Two different selectors put pathways on these maps and neither guarantees
    # an FDR: .kegg_hits_by_contrast() scores on padj where a layer's table
    # carries usable ones and on pvalue where it does not, while the supported
    # renderer's tiers drop to raw p when the adjusted ones select nothing --
    # even though usable adjusted values are present. A caption promising an
    # FDR throughout would overstate the evidence on exactly those runs.
    cap <- pathview_significance_caption()

    expect_true(grepl("adjusted p-values", cap, fixed = TRUE))
    expect_true(grepl("raw p-values", cap, fixed = TRUE))
    expect_true(grepl("floor rather than an FDR", cap, fixed = TRUE))
})

test_that("the caption does not name a selector that only one renderer uses", {
    # It previously attributed its fallback rule to .kegg_hits_by_contrast(),
    # which the supported-space renderer never calls -- so the caption
    # described a rule that did not produce the figure the report prefers.
    # The claim has to hold for whichever renderer wrote the PDF.
    expect_false(grepl(".kegg_hits_by_contrast",
                       pathview_significance_caption(), fixed = TRUE))
})

test_that("the supported renderer records whether compounds reached the map", {
    # Before the node rule, "metabolomics was in the run" was close enough to
    # "a compound node carries a value" for the report to assume it. It is not
    # any more: every metabolite can fail the rule while the gene layers still
    # produce a figure, and the report would then promise a colour the map does
    # not have. Pinned at the source -- reaching it needs a pathview call.
    body_src <- paste(deparse(body(generate_multi_ora_pathview)), collapse = " ")

    expect_true(grepl("any_compounds <- TRUE", body_src, fixed = TRUE))
    expect_true(grepl("multi_ora_pathview_supported.yaml", body_src,
                      fixed = TRUE))
})

test_that("both pathview renderers apply the node rule the caption states", {
    # The report attaches one caption to whichever PDF exists, and prefers the
    # supported-space one. Filtering only the KO-space renderer would leave that
    # figure colouring unchanged features under a caption saying otherwise.
    for (fn in list(generate_per_omic_union_pathview, generate_multi_ora_pathview)) {
        body_src <- paste(deparse(body(fn)), collapse = " ")
        expect_true(grepl("filter_changed_features", body_src, fixed = TRUE))
    }
})
