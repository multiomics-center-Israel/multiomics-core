# One layer's top pathways, looked up in another (07d_cross_lookup.R).
#
# All fixtures are synthetic: made-up accessions and p-values.

# A metabolite layer keyed "<accession> - <name>", as some compound GMTs name
# their sets, and scored by GSEA.
metab_gsea <- function() {
    data.frame(
        pathway = c("rno00010 - Pathway A", "rno00020 - Pathway B",
                    "rno00030 - Pathway C", "rno00040 - Pathway D"),
        pvalue  = c(0.01, 0.02, 0.04, 0.5),
        padj    = c(0.2, 0.2, 0.3, 0.6),
        NES     = c(-1.9, 1.7, -1.4, 0.6),
        setSize = c(6, 5, 4, 3),
        method  = "fgsea",
        stringsAsFactors = FALSE
    )
}

# A protein layer keyed on the bare accession, with fgsea and ORA rows side by
# side as run_pathway_analysis() writes them. It has nothing for rno00020.
prot_mixed <- function() {
    dplyr::bind_rows(
        data.frame(pathway = c("rno00010", "rno00030", "rno00050"),
                   pval = c(0.3, 0.001, 0.02), padj = c(0.6, 0.05, 0.2),
                   NES = c(1.2, 2.1, -1.5), size = c(40, 25, 30),
                   method = "fgsea", stringsAsFactors = FALSE),
        data.frame(pathway = "rno00010", pvalue = 0.04, padj = 0.3,
                   method = "ora", direction = "up", stringsAsFactors = FALSE)
    )
}

tables <- function() list(metabolomics = metab_gsea(), proteomics = prot_mixed())


# ---- joining and status ----------------------------------------------------

test_that("an accession carrying a name joins the bare accession", {
    lk <- build_cross_omics_lookup(tables(), "metabolomics", "proteomics",
                                   top_n = 10, kegg_org = "rno")
    row <- lk[lk$ranking == "all" & lk$norm_id == "00010", ]

    expect_identical(row$target_status, "tested")
    expect_equal(row$target_p, 0.3)
    expect_equal(row$target_NES, 1.2)
})

test_that("a pathway the target has no result for is reported, not dropped", {
    lk <- build_cross_omics_lookup(tables(), "metabolomics", "proteomics",
                                   top_n = 10, kegg_org = "rno")
    row <- lk[lk$ranking == "all" & lk$norm_id == "00020", ]

    expect_identical(row$target_status, "not in target results")
    expect_true(is.na(row$target_p))
})

test_that("the tested ranking holds only pathways the target has", {
    lk <- build_cross_omics_lookup(tables(), "metabolomics", "proteomics",
                                   top_n = 10, kegg_org = "rno")
    tested <- lk[lk$ranking == "tested_in_target", ]

    expect_setequal(tested$norm_id, c("00010", "00030"))
    expect_true(all(tested$target_status == "tested"))
    expect_identical(tested$rank, seq_len(nrow(tested)))
})

test_that("no adjustment is recomputed over the looked-up subset", {
    # A BH over pathways the source already selected reads as an ordinary FDR,
    # which it is not; only the target's own p and adjusted p are reported.
    lk <- build_cross_omics_lookup(tables(), "metabolomics", "proteomics",
                                   top_n = 10, kegg_org = "rno")
    expect_false("target_padj_within_lookup" %in% names(lk))
    row <- lk[lk$ranking == "all" & lk$norm_id == "00010", ]
    expect_equal(row$target_padj, 0.6)
})


# ---- which rows are read ----------------------------------------------------

test_that("the target is read from its rank-based rows, with ORA beside them", {
    lk <- build_cross_omics_lookup(tables(), "metabolomics", "proteomics",
                                   top_n = 10, kegg_org = "rno")
    row <- lk[lk$ranking == "all" & lk$norm_id == "00010", ]

    # The ORA p (0.04) is smaller than the fgsea p (0.3) but does not replace it.
    expect_identical(row$target_method, "fgsea")
    expect_equal(row$target_ora_p, 0.04)
    expect_equal(row$target_n_measured, 40)
})

test_that("a source layer with no rank-based rows is ranked on its ORA rows", {
    ora_only <- data.frame(pathway = c("rno00010", "rno00030"),
                           pvalue = c(0.02, 0.01), padj = c(0.1, 0.1),
                           method = "ora", stringsAsFactors = FALSE)
    lk <- build_cross_omics_lookup(list(proteomics = ora_only,
                                        metabolomics = metab_gsea()),
                                   "proteomics", "metabolomics",
                                   top_n = 10, kegg_org = "rno")
    blk <- lk[lk$ranking == "all", ]

    expect_identical(unique(blk$source_method), "ora")
    expect_true(all(is.na(blk$source_NES)))
    expect_true(all(is.na(blk$source_n_measured)))
    expect_identical(blk$norm_id, c("00030", "00010"))
})

test_that("a source layer with nothing usable gives NULL", {
    empty <- data.frame(pathway = character(0), pvalue = numeric(0))
    expect_null(build_cross_omics_lookup(list(proteomics = empty,
                                              metabolomics = metab_gsea()),
                                         "proteomics", "metabolomics"))
})


# ---- ranking ---------------------------------------------------------------

test_that("the source is ranked on its own raw p-value", {
    lk <- build_cross_omics_lookup(tables(), "metabolomics", "proteomics",
                                   top_n = 10, kegg_org = "rno")
    blk <- lk[lk$ranking == "all", ]
    expect_identical(blk$norm_id, c("00010", "00020", "00030", "00040"))
})

test_that("ties on p-value go to the larger |NES|, then to the pathway id", {
    # 00030 arrives before 00010 with the same p and |NES|: the pathway id, not
    # the incoming order, decides.
    tied <- data.frame(pathway = c("rno00030", "rno00020", "rno00010"),
                       pvalue = c(0.05, 0.05, 0.05), NES = c(1.1, -2.0, -1.1),
                       method = "fgsea", stringsAsFactors = FALSE)
    lk <- build_cross_omics_lookup(list(a = tied, b = metab_gsea()), "a", "b",
                                   top_n = 3, kegg_org = "rno")
    expect_identical(lk$norm_id[lk$ranking == "all"], c("00020", "00010", "00030"))

    reversed <- tied[rev(seq_len(nrow(tied))), ]
    lk2 <- build_cross_omics_lookup(list(a = reversed, b = metab_gsea()), "a", "b",
                                    top_n = 3, kegg_org = "rno")
    expect_identical(lk2$norm_id[lk2$ranking == "all"], c("00020", "00010", "00030"))
})

test_that("a top_n larger than the table returns every pathway", {
    lk <- build_cross_omics_lookup(tables(), "metabolomics", "proteomics",
                                   top_n = 500, kegg_org = "rno")
    expect_equal(sum(lk$ranking == "all"), 4)
})

test_that("an unusable top_n falls back to the default", {
    lk <- build_cross_omics_lookup(tables(), "metabolomics", "proteomics",
                                   top_n = NA, kegg_org = "rno")
    expect_equal(sum(lk$ranking == "all"), 4)
})

test_that("one pathway across contrasts keeps its best row and says which", {
    two <- rbind(transform(metab_gsea(), contrast = "A_vs_B"),
                 transform(metab_gsea(), contrast = "C_vs_D",
                           pvalue = c(0.001, 0.9, 0.9, 0.9)))
    lk <- build_cross_omics_lookup(list(metabolomics = two,
                                        proteomics = prot_mixed()),
                                   "metabolomics", "proteomics",
                                   top_n = 10, kegg_org = "rno")
    row <- lk[lk$ranking == "all" & lk$norm_id == "00010", ]

    expect_equal(row$source_p, 0.001)
    expect_identical(row$source_contrast, "C_vs_D")
    expect_equal(sum(lk$ranking == "all" & lk$norm_id == "00010"), 1)
})


# ---- config ----------------------------------------------------------------

test_that("the lookup is off by default, with fifteen pathways when on", {
    expect_identical(.cross_lookup_config(list()), list(enabled = FALSE, top_n = 15L))
})

test_that("the reader's defaults match the config validator's", {
    # Two paths to one answer: the validator fills the key, and the reader
    # defaults it for a config that bypassed the validator.
    validated <- suppressMessages(suppressWarnings(validate_multiomics_config(
        list(integration = list(methods = "SNF"))
    )))
    defaults <- .cross_lookup_config(list())
    expect_identical(isTRUE(validated$enrichment$cross_lookup$enabled),
                     defaults$enabled)
    expect_equal(validated$enrichment$cross_lookup$top_n, defaults$top_n)
})

test_that("the config can switch the lookup off and set top_n", {
    cfg <- list(modes = list(multiomics = list(enrichment = list(
        cross_lookup = list(enabled = FALSE, top_n = 5)))))
    expect_identical(.cross_lookup_config(cfg), list(enabled = FALSE, top_n = 5L))
})


# ---- files -----------------------------------------------------------------

test_that("both directions are written as a table and two figures each", {
    skip_if_not_installed("ggplot2")
    out_dir <- withr::local_tempdir()

    written <- suppressMessages(write_cross_omics_lookups(
        tables(), c("metabolomics", "proteomics"), out_dir,
        top_n = 10, kegg_org = "rno"))

    expect_setequal(names(written),
                    c("metabolomics_to_proteomics", "proteomics_to_metabolomics"))
    for (stem in c("metabolomics_to_proteomics", "proteomics_to_metabolomics")) {
        base <- file.path(out_dir, paste0("cross_lookup_", stem))
        expect_true(file.exists(paste0(base, ".tsv")))
        expect_true(file.exists(paste0(base, ".png")))
        expect_true(file.exists(paste0(base, "_tested.png")))
    }
})

test_that("the written table reads back with the columns the report shows", {
    skip_if_not_installed("ggplot2")
    out_dir <- withr::local_tempdir()
    suppressMessages(write_cross_omics_lookups(
        tables(), c("metabolomics", "proteomics"), out_dir,
        top_n = 10, kegg_org = "rno"))

    tbl <- read.delim(file.path(out_dir, "cross_lookup_metabolomics_to_proteomics.tsv"),
                      stringsAsFactors = FALSE)
    expect_true(all(c("ranking", "pathway", "source_NES", "target_status",
                      "source_contrast", "target_contrast") %in% names(tbl)))
    expect_false("target_padj_within_lookup" %in% names(tbl))
})

test_that("a previous run's lookup files are cleared", {
    out_dir <- withr::local_tempdir()
    stale <- file.path(out_dir, c("cross_lookup_a_to_b.tsv", "cross_lookup_a_to_b.png"))
    for (f in stale) writeLines("left over from an earlier run", f)
    bystander <- file.path(out_dir, "cross_omics_pathways_meta_analysis.csv")
    writeLines("kept", bystander)

    .clear_cross_lookup_outputs(out_dir)

    expect_false(any(file.exists(stale)))
    expect_true(file.exists(bystander))
})

test_that("switching the lookup off writes none and clears the old ones", {
    out_dir <- withr::local_tempdir()
    stale <- file.path(out_dir, "cross_lookup_metabolomics_to_proteomics.tsv")
    writeLines("left over from an earlier run", stale)
    cfg <- list(global = list(organism = "Unlisted nonmodel species"),
                modes = list(multiomics = list(enrichment = list(
                    cross_lookup = list(enabled = FALSE)))))

    res <- suppressWarnings(suppressMessages(analyze_cross_omics_enrichment(
        list(proteomics = prot_mixed(), metabolomics = metab_gsea()),
        cfg, out_dir = out_dir)))

    expect_length(list.files(out_dir, pattern = "^cross_lookup_"), 0)
    expect_length(res$lookup_files, 0)
})

test_that("the lookup is drawn as an NES heatmap, one column per layer", {
    skip_if_not_installed("ggplot2")
    lk <- build_cross_omics_lookup(tables(), "metabolomics", "proteomics", top_n = 5,
                                   kegg_org = "rno")
    blk <- lk[lk$ranking == "all", , drop = FALSE]
    blk <- blk[order(blk$rank), , drop = FALSE]
    p <- .lookup_heatmap(blk)

    expect_s3_class(p, "ggplot")
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true("GeomTile" %in% geoms)
    # A missing NES is left unfilled, so it cannot be read as a value near zero,
    # which is what the scale's own centre now looks like.
    fill_scale <- p$scales$scales[[which(vapply(p$scales$scales,
        function(sc) "fill" %in% sc$aesthetics, logical(1)))[1]]]
    expect_true(is.na(fill_scale$na.value))
    expect_true(grepl("NES", paste(deparse(p$mapping$fill), collapse = "")))
    expect_identical(levels(p$data$layer),
                     c("metabolomics\n(ranked)", "proteomics\n(looked up)"))
    # A pathway the target has no result for is a dash, not a number.
    untested <- blk$pathway[blk$target_status != "tested"]
    if (length(untested) > 0) {
        tgt <- p$data[p$data$layer == "proteomics\n(looked up)", ]
        expect_true(all(tgt$label[is.na(tgt$NES)] %in% c("–", "ORA")))
    }
})

test_that("the heatmap is written to the path it was given", {
    skip_if_not_installed("ggplot2")
    skip_if_not(isTRUE(capabilities("png")), "no png device")
    lk <- build_cross_omics_lookup(tables(), "metabolomics", "proteomics", top_n = 5)
    out <- file.path(withr::local_tempdir(), "lookup.png")
    expect_identical(plot_cross_omics_lookup(lk, "all", out), out)
    expect_true(file.exists(out))
})


# ---- the target is read in the source row's own contrast --------------------

# Per-contrast tables as the run-level merge sees them: each layer spells the
# comparison its own way.
metab_two_contrasts <- function() {
    rbind(transform(metab_gsea(), contrast = "A_vs_B"),
          transform(metab_gsea(), contrast = "C_vs_D",
                    pvalue = c(0.3, 0.9, 0.9, 0.9)))
}

prot_two_contrasts <- function() {
    data.frame(pathway = c("rno00010", "rno00010", "rno00030"),
               pval = c(0.4, 0.001, 0.02), padj = c(0.7, 0.01, 0.2),
               NES = c(1.2, -2.5, 1.6), size = c(40, 40, 25),
               method = "fgsea", contrast = c("A vs B", "C_vs_D", "C_vs_D"),
               stringsAsFactors = FALSE)
}

test_that("run-level source and target come from the same normalized contrast", {
    lk <- build_cross_omics_lookup(
        list(metabolomics = metab_two_contrasts(), proteomics = prot_two_contrasts()),
        "metabolomics", "proteomics", top_n = 10, kegg_org = "rno")
    row <- lk[lk$ranking == "all" & lk$norm_id == "00010", ]

    # The source's best 00010 row is A_vs_B. The target's C_vs_D row has the far
    # smaller p but is a different comparison; its A vs B row is the one read.
    expect_identical(row$source_contrast, "A_vs_B")
    expect_identical(row$target_contrast, "A vs B")
    expect_equal(row$target_p, 0.4)
    expect_equal(row$target_NES, 1.2)
    expect_identical(normalize_contrast_key(row$source_contrast),
                     normalize_contrast_key(row$target_contrast))
})

test_that("a target missing in the source's contrast is not filled from another", {
    lk <- build_cross_omics_lookup(
        list(metabolomics = metab_two_contrasts(), proteomics = prot_two_contrasts()),
        "metabolomics", "proteomics", top_n = 10, kegg_org = "rno")
    row <- lk[lk$ranking == "all" & lk$norm_id == "00030", ]

    # The target has 00030 only in C_vs_D; the source's 00030 row is A_vs_B.
    expect_identical(row$source_contrast, "A_vs_B")
    expect_identical(row$target_status, "not in target results for this contrast")
    expect_true(is.na(row$target_p))
    expect_true(is.na(row$target_contrast))
    expect_false("00030" %in% lk$norm_id[lk$ranking == "tested_in_target"])
})

test_that("a target that names no contrast is not matched to a source that does", {
    lk <- build_cross_omics_lookup(
        list(metabolomics = metab_two_contrasts(), proteomics = prot_mixed()),
        "metabolomics", "proteomics", top_n = 10, kegg_org = "rno")
    expect_true(all(lk$target_status[lk$ranking == "all"] ==
                        "not in target results for this contrast"))
})


# ---- which rows are read: .layer_contribution() is the path ----------------

# Sourced functions live in the global environment rather than a package
# namespace, so a stub is assigned there and restored on exit -- the pattern of
# the other enrichment tests.
local_stubs <- function(stubs, env = parent.frame()) {
    target <- environment(build_cross_omics_lookup)
    nms <- names(stubs)
    had <- vapply(nms, exists, logical(1), envir = target, inherits = FALSE)
    old <- lapply(nms[had], get, envir = target, inherits = FALSE)
    names(old) <- nms[had]
    withr::defer({
        for (nm in nms) {
            if (nm %in% names(old)) {
                assign(nm, old[[nm]], envir = target)
            } else if (exists(nm, envir = target, inherits = FALSE)) {
                rm(list = nm, envir = target)
            }
        }
    }, envir = env)
    for (nm in nms) assign(nm, stubs[[nm]], envir = target)
    invisible(NULL)
}

test_that("the rows read are the ones .layer_contribution() keeps", {
    real <- .layer_contribution
    # Keep only the source's rno00040 row: if the lookup chose its rows any
    # other way, the other three pathways would appear.
    local_stubs(list(.layer_contribution = function(df, kegg_org = NULL) {
        cb <- real(df, kegg_org)
        if (identical(df$pathway, metab_gsea()$pathway)) {
            cb$keep <- cb$keep & grepl("00040", df$pathway)
        }
        cb
    }))
    lk <- build_cross_omics_lookup(tables(), "metabolomics", "proteomics",
                                   top_n = 10, kegg_org = "rno")
    expect_identical(lk$norm_id[lk$ranking == "all"], "00040")
})

test_that("the lookup no longer selects rows on its own", {
    src <- paste(deparse(body(build_cross_omics_lookup)), collapse = " ")
    expect_true(grepl("\\.layer_contribution\\(", src, perl = TRUE))
    expect_false(grepl("select_layer_method_rows\\(", src, perl = TRUE))
})


# ---- run level: read from merge_tables after exclusions, rank tables in ----

test_that("the run-level lookup reads the merged tables after exclusion", {
    out_dir <- withr::local_tempdir()
    prot <- data.frame(ID = c("map00010", "map00020", "map00030"),
                       pval = c(0.01, 0.02, 0.03), padj = c(0.1, 0.1, 0.1),
                       NES = c(1.5, -1.4, 1.2), method = "fgsea",
                       contrast = "A_vs_B", stringsAsFactors = FALSE)
    compound_gsea <- data.frame(ID = c("map00010", "map00020"),
                                pvalue = c(0.04, 0.05), padj = c(0.2, 0.2),
                                NES = c(-1.3, 1.1), method = "fgsea",
                                contrast = "A vs B", stringsAsFactors = FALSE)
    # A class exclusion that drops 00020, without the KEGG classification.
    local_stubs(list(
        .excluded_pathway_classes = function(config) "Stubbed class",
        keep_kegg_pathways = function(ids, ...) !grepl("00020", ids)))
    cfg <- list(global = list(organism = "Unlisted nonmodel species"),
                modes = list(multiomics = list(enrichment = list(
                    cross_lookup = list(enabled = TRUE, top_n = 10)))))

    res <- suppressWarnings(suppressMessages(analyze_cross_omics_enrichment(
        list(proteomics = prot), cfg, out_dir = out_dir,
        rank_tables = list(metabolomics = compound_gsea))))

    tsv <- file.path(out_dir, "cross_lookup_proteomics_to_metabolomics.tsv")
    expect_true(file.exists(tsv))
    lk <- read.delim(tsv, stringsAsFactors = FALSE, colClasses = c(norm_id = "character"))
    # The excluded pathway is gone; the metabolomics layer is there only as the
    # rank table supplied beside the per-layer ones.
    expect_false("00020" %in% lk$norm_id)
    row <- lk[lk$ranking == "all" & lk$norm_id == "00010", ]
    expect_identical(row$target_method, "fgsea")
    expect_equal(row$target_p, 0.04)
})


# ---- the report shows where each value comes from ----------------------------

test_that("the report's lookup table carries both contrasts", {
    rmd <- paste(readLines(testthat::test_path(
        "..", "..", "R", "domain", "multiomics", "report_template_multiomics.Rmd")),
        collapse = "\n")
    cols_def <- regmatches(rmd, regexpr("lookup_cols <- c\\([^)]*\\)", rmd))
    expect_length(cols_def, 1)
    expect_true(grepl('"source_contrast"', cols_def, fixed = TRUE))
    expect_true(grepl('"target_contrast"', cols_def, fixed = TRUE))
    expect_false(grepl("target_padj_within_lookup", rmd, fixed = TRUE))
})
