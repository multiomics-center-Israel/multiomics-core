# Excluding KEGG pathway classes a project does not report.
#
# KEGG's reference maps are pan-species, so a run with no KEGG organism code of
# its own is tested against maps of organs its organism does not have. The
# exclusion is a reporting layer: it drops finished rows, and never touches the
# tested universe, the p-values or the adjustment behind them.
#
# All fixtures are synthetic; nothing here reaches KEGG or the network.

# A miniature br08901, in the shape the real hierarchy has: A categories over B
# subcategories over C pathway lines, plus lines that are not pathways at all.
brite_fixture <- function() {
    c(
        "+D\tPathway",
        "#<h2>KEGG Pathway Maps</h2>",
        "!",
        # Headings carry their own BRITE hierarchy code before the readable
        # name, which the config never spells.
        "A<b>09100 Metabolism</b>",
        "B  09101 Carbohydrate metabolism",
        "C    00010  Glycolysis / Gluconeogenesis [PATH:map00010]",
        "C    00020  Citrate cycle (TCA cycle) [PATH:map00020]",
        "B  09102 Energy metabolism",
        "C    00190  Oxidative phosphorylation [PATH:map00190]",
        "A<b>09150 Organismal Systems</b>",
        "B  09153 Circulatory system",
        "C    04260  Cardiac muscle contraction [PATH:map04260]",
        "B  09151 Immune system",
        "C    04620  Toll-like receptor signaling pathway [PATH:map04620]",
        # ...and a heading with no code at all, since KEGG has shipped this
        # hierarchy in more than one shape.
        "A<b>Uncoded Category</b>",
        "B  Uncoded Subcategory",
        "C    09999  Placeholder pathway [PATH:map09999]",
        "!",
        "#Last updated"
    )
}


# ---- reading the hierarchy --------------------------------------------------

test_that("parse_kegg_brite_pathways reads A/B/C and ignores the rest", {
    cls <- parse_kegg_brite_pathways(brite_fixture())

    expect_named(cls, c("pathway_id", "category", "subcategory", "pathway_name"))
    expect_equal(nrow(cls), 6L)
    expect_identical(cls$pathway_name[cls$pathway_id == "00010"],
                     "Glycolysis / Gluconeogenesis [PATH:map00010]")
    # The "+D", "#" and "!" lines are not pathways.
    expect_false(any(grepl("^[+#!]", cls$pathway_id)))
})

test_that("headings expose the readable name, not its hierarchy code", {
    # The config says "Organismal Systems"; the file says
    # "A<b>09150 Organismal Systems</b>". If the code survived into the parsed
    # value, every exclusion a project writes would silently match nothing.
    cls <- parse_kegg_brite_pathways(brite_fixture())

    expect_identical(cls$category[cls$pathway_id == "00010"], "Metabolism")
    expect_identical(cls$subcategory[cls$pathway_id == "00010"],
                     "Carbohydrate metabolism")
    expect_identical(cls$category[cls$pathway_id == "04620"], "Organismal Systems")
    expect_identical(cls$subcategory[cls$pathway_id == "04620"], "Immune system")

    # Nothing parsed still carries a leading code.
    expect_false(any(grepl("^[0-9]{5}\\s", c(cls$category, cls$subcategory))))
})

test_that("a heading with no hierarchy code is kept as it is", {
    cls <- parse_kegg_brite_pathways(brite_fixture())

    expect_identical(cls$category[cls$pathway_id == "09999"], "Uncoded Category")
    expect_identical(cls$subcategory[cls$pathway_id == "09999"], "Uncoded Subcategory")
})

test_that("parse_kegg_brite_pathways returns NULL rather than an empty table", {
    # "Nothing is classified" and "exclude everything" must not look alike to a
    # caller, so an unusable hierarchy is absent, not empty.
    expect_null(parse_kegg_brite_pathways(NULL))
    expect_null(parse_kegg_brite_pathways(character(0)))
    expect_message(
        expect_null(parse_kegg_brite_pathways(c("A<b>Metabolism</b>", "B  Nothing here"))),
        "no pathway lines found"
    )
})

test_that("a cached object of the wrong shape is not trusted", {
    expect_false(.is_kegg_category_table(NULL))
    expect_false(.is_kegg_category_table(data.frame()))
    expect_false(.is_kegg_category_table(data.frame(pathway_id = "00010")))
    expect_true(.is_kegg_category_table(parse_kegg_brite_pathways(brite_fixture())))
})

test_that("kegg_pathway_categories reads a valid cache without fetching", {
    cache_dir <- withr::local_tempdir()
    saveRDS(parse_kegg_brite_pathways(brite_fixture()),
            file.path(cache_dir, "kegg_pathway_categories.rds"))

    cls <- kegg_pathway_categories(cache_dir = cache_dir)

    expect_equal(nrow(cls), 6L)
})


# ---- which pathways survive -------------------------------------------------

with_fixture_cache <- function() {
    cache_dir <- withr::local_tempdir(.local_envir = parent.frame())
    saveRDS(parse_kegg_brite_pathways(brite_fixture()),
            file.path(cache_dir, "kegg_pathway_categories.rds"))
    cache_dir
}

test_that("an excluded category drops its pathways and keeps the others", {
    cache <- with_fixture_cache()
    ids <- c("map00010", "map04260", "map04620")

    keep <- keep_kegg_pathways(ids, exclude = "Organismal Systems",
                               cache_dir = cache)

    expect_identical(keep, c(TRUE, FALSE, FALSE))
})

test_that("a subcategory can be excluded while its siblings stay", {
    cache <- with_fixture_cache()
    ids <- c("map04260", "map04620")

    # An insect has no circulatory system in the vertebrate sense, but it does
    # have innate immunity -- which is the whole reason subcategories matter.
    keep <- keep_kegg_pathways(ids, exclude = "Circulatory system",
                               cache_dir = cache)

    expect_identical(keep, c(FALSE, TRUE))
})

test_that("every KEGG spelling of one pathway classifies alike", {
    cache <- with_fixture_cache()
    ids <- c("04260", "map04260", "ko04260", "hsa04260",
             "hsa04260 Cardiac muscle contraction")

    keep <- keep_kegg_pathways(ids, exclude = "Organismal Systems",
                               kegg_org = "hsa", cache_dir = cache)

    expect_true(all(!keep))
})

test_that("identifiers that are not KEGG accessions are never excluded", {
    cache <- with_fixture_cache()
    # A GO term, a custom gene-set name, and an id shaped like an accession but
    # belonging to another organism than this run's.
    ids <- c("GO:0006915", "My custom set", "mmu04260")

    keep <- keep_kegg_pathways(ids, exclude = "Organismal Systems",
                               kegg_org = "hsa", cache_dir = cache)

    expect_true(all(keep))
})

test_that("an unclassified accession is kept", {
    cache <- with_fixture_cache()
    # Not in the hierarchy: a newer map, or one the fixture does not list.
    expect_true(keep_kegg_pathways("map08888", exclude = "Metabolism",
                                   cache_dir = cache))
})

test_that("no exclusion list means nothing is excluded", {
    cache <- with_fixture_cache()
    ids <- c("map00010", "map04260")

    expect_true(all(keep_kegg_pathways(ids, exclude = NULL, cache_dir = cache)))
    expect_true(all(keep_kegg_pathways(ids, exclude = character(0), cache_dir = cache)))
    expect_true(all(keep_kegg_pathways(ids, exclude = list(), cache_dir = cache)))
})

test_that("an unavailable classification keeps everything and says so", {
    # Fail open. The classification is passed in as absent rather than left to a
    # fetch, so this asserts the contract on every machine instead of only on
    # one with no network.
    ids <- c("map00010", "map04260")

    expect_message(
        keep <- keep_kegg_pathways(ids, exclude = "Metabolism",
                                   classification = NULL),
        "keeping all"
    )
    expect_true(all(keep))
})

test_that(".excluded_pathway_classes reads the key, and copes without it", {
    expect_identical(.excluded_pathway_classes(list()), character(0))

    cfg <- list(modes = list(multiomics = list(enrichment = list(
        exclude_pathway_classes = list("Human Diseases", "Circulatory system")
    ))))
    expect_identical(.excluded_pathway_classes(cfg),
                     c("Human Diseases", "Circulatory system"))

    cfg$modes$multiomics$enrichment$exclude_pathway_classes <- list()
    expect_identical(.excluded_pathway_classes(cfg), character(0))
})


# ---- the exclusion does not touch the statistics ----------------------------

test_that("excluding a class leaves the retained rows' p-values untouched", {
    cache <- with_fixture_cache()

    # Two metabolites per pathway, one of them significant, so the compound ORA
    # has something to test in each. Synthetic throughout.
    de <- data.frame(
        feature_id = paste0("m", 1:6),
        KEGG_ID = c("C00001", "C00002", "C00003", "C00004", "C00005", "C00006"),
        padj = c(0.001, 0.9, 0.001, 0.9, 0.001, 0.9),
        log2fc = c(2, 0.1, 2, 0.1, 2, 0.1),
        stringsAsFactors = FALSE
    )
    # The shape get_kegg_compound_pathways() caches: one row per
    # (pathway, compound), with a readable name alongside.
    sets <- data.frame(
        pathway  = rep(c("00010", "04260", "04620"), each = 2),
        compound = paste0("C0000", 1:6),
        name     = rep(c("Glycolysis", "Cardiac muscle contraction",
                         "Toll-like receptor signaling"), each = 2),
        stringsAsFactors = FALSE
    )
    saveRDS(sets, file.path(cache, "kegg_compound_pathways.rds"))

    unfiltered <- suppressMessages(
        run_compound_ora(de, cache, 1, 500, 1, universe = de$KEGG_ID))
    skip_if(is.null(unfiltered) || nrow(unfiltered) == 0,
            "compound ORA fixture produced no rows in this environment")

    filtered <- suppressMessages(
        run_compound_ora(de, cache, 1, 500, 1, universe = de$KEGG_ID,
                         exclude_classes = "Organismal Systems"))

    # The excluded class is gone...
    expect_false(any(filtered$ID %in% c("04260", "04620")))
    # ...and what remains carries exactly the values it had before, p-value and
    # adjustment alike: the tested family did not change, only the report did.
    kept <- unfiltered[unfiltered$ID %in% filtered$ID, ]
    expect_equal(filtered$pvalue, kept$pvalue)
    expect_equal(filtered$padj, kept$padj)
})
