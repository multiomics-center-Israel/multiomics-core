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
    expect_true(is.na(row$target_padj_within_lookup))
})

test_that("the tested ranking holds only pathways the target has", {
    lk <- build_cross_omics_lookup(tables(), "metabolomics", "proteomics",
                                   top_n = 10, kegg_org = "rno")
    tested <- lk[lk$ranking == "tested_in_target", ]

    expect_setequal(tested$norm_id, c("00010", "00030"))
    expect_true(all(tested$target_status == "tested"))
    expect_identical(tested$rank, seq_len(nrow(tested)))
})

test_that("the within-lookup adjustment covers the tested rows of one ranking", {
    lk <- build_cross_omics_lookup(tables(), "metabolomics", "proteomics",
                                   top_n = 10, kegg_org = "rno")
    blk <- lk[lk$ranking == "all", ]
    tested <- blk$target_status == "tested"

    expect_equal(blk$target_padj_within_lookup[tested],
                 stats::p.adjust(blk$target_p[tested], method = "BH"))
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

test_that("ties on p-value go to the larger |NES|, then to incoming order", {
    tied <- data.frame(pathway = c("rno00010", "rno00020", "rno00030"),
                       pvalue = c(0.05, 0.05, 0.05), NES = c(1.1, -2.0, 1.1),
                       method = "fgsea", stringsAsFactors = FALSE)
    lk <- build_cross_omics_lookup(list(a = tied, b = metab_gsea()), "a", "b",
                                   top_n = 3, kegg_org = "rno")
    expect_identical(lk$norm_id[lk$ranking == "all"], c("00020", "00010", "00030"))
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

test_that("the lookup is on by default with fifteen pathways", {
    expect_identical(.cross_lookup_config(list()), list(enabled = TRUE, top_n = 15L))
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
                      "target_padj_within_lookup") %in% names(tbl)))
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
