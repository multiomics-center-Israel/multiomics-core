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
             "hsa04260 Cardiac muscle contraction",
             "hsa04260~Cardiac muscle contraction")

    keep <- keep_kegg_pathways(ids, exclude = "Organismal Systems",
                               kegg_org = "hsa", cache_dir = cache)

    expect_true(all(!keep))
})

test_that("identifiers that are not KEGG accessions are never excluded", {
    cache <- with_fixture_cache()
    # A GO term, a custom gene-set name, and an id shaped like an accession but
    # belonging to another organism than this run's -- bare and tilde-labelled.
    ids <- c("GO:0006915", "My custom set", "mmu04260",
             "GO:0006915~Apoptosis", "My~custom set", "mmu04260~Cardiac muscle")

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


# ---- the gene-based multi-ORA boundary --------------------------------------
#
# The gene-based path runs its own KEGG ORA and its tables feed the multi-ORA
# summary, the pooled/dot/support plots and the OrgDb pathview renderers. It
# used to bypass the exclusion entirely, so a class could vanish from the
# cross-omics tables and reappear here one section away. Filtered once, at the
# boundary, rather than in each renderer downstream.

gene_ora_table <- function() {
    data.frame(
        pathway = c("Glycolysis / Gluconeogenesis", "Cardiac muscle contraction"),
        ID      = c("00010", "04260"),
        pvalue  = c(0.001, 0.002),
        padj    = c(0.01, 0.02),
        Count   = c(5L, 4L),
        stringsAsFactors = FALSE
    )
}

test_that(".exclude_kegg_classes drops the class and changes nothing else", {
    cache <- with_fixture_cache()
    df <- gene_ora_table()

    out <- .exclude_kegg_classes(df, "Organismal Systems", kegg_org = "hsa",
                                 classification = kegg_pathway_categories(cache))

    expect_identical(out$ID, "00010")
    # The row that stays carries exactly the values it arrived with: the tested
    # universe and the adjustment behind them were never touched.
    expect_equal(out$pvalue, df$pvalue[df$ID == "00010"])
    expect_equal(out$padj, df$padj[df$ID == "00010"])
    # Column shape survives, since the summary and plots consume this frame.
    expect_named(out, names(df))
})

test_that(".exclude_kegg_classes is a no-op without an exclusion list", {
    df <- gene_ora_table()

    expect_identical(.exclude_kegg_classes(df, NULL), df)
    expect_identical(.exclude_kegg_classes(df, character(0)), df)
    expect_identical(.exclude_kegg_classes(df, list()), df)
    expect_null(.exclude_kegg_classes(NULL, "Metabolism"))
})

test_that(".exclude_kegg_classes reports nothing left as no pathways", {
    cache <- with_fixture_cache()

    # Emptying the table must look like every other "no pathways" path here,
    # or a zero-row frame reaches plotting code that expects rows or NULL.
    expect_null(.exclude_kegg_classes(
        gene_ora_table()[2, , drop = FALSE], "Organismal Systems",
        kegg_org = "hsa", classification = kegg_pathway_categories(cache)))
})

test_that("both gene-ORA wiring paths accept an exclusion list", {
    # The run-level and per-contrast callers reach run_multi_ora_kegg() through
    # different argument styles; neither may drop the exclusion on the floor.
    expect_true("exclude_classes" %in% names(formals(run_multi_ora_kegg)))

    run_src <- paste(deparse(body(run_multi_ora)), collapse = " ")
    grp_src <- paste(deparse(body(.run_multi_ora_contrast_group)), collapse = " ")

    # Run level: pooled and per-omics both read the config through one helper.
    expect_equal(
        lengths(regmatches(run_src, gregexpr("run_multi_ora_kegg", run_src)))[[1]], 2L)
    expect_true(grepl(".excluded_pathway_classes", run_src, fixed = TRUE))
    # Per contrast: both calls pass the list handed down to the helper.
    expect_equal(
        lengths(regmatches(grp_src, gregexpr("run_multi_ora_kegg", grp_src)))[[1]], 2L)
    expect_true(grepl("exclude_classes = exclude_classes", grp_src, fixed = TRUE))
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


# ---- run_multi_ora_kegg(): completion message and fallback (#222) -----------
#
# enrichKEGG is mocked in the clusterProfiler namespace, as test-loadings-gsea.R
# mocks gseKEGG; the Fisher fallback and the classification fetch are sourced
# globals, so they are stubbed by assignment into the function's environment.
# Nothing here reaches KEGG.

local_kegg_ora_stubs <- function(stubs, env = parent.frame()) {
    target <- environment(run_multi_ora_kegg)
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

# enrichKEGG's result, as as.data.frame() would give it.
fake_enrichkegg_table <- function(ids = c("hsa00010", "hsa04260"),
                                  pvalue = c(0.001, 0.002),
                                  p.adjust = c(0.01, 0.02)) {
    data.frame(ID = ids, Description = paste("Pathway", ids),
               pvalue = pvalue, p.adjust = p.adjust,
               GeneRatio = "3/10", Count = 3L, geneID = "hsa:1/hsa:2/hsa:3",
               stringsAsFactors = FALSE)
}

# Wires enrichKEGG to `enrich`, records whether the Fisher fallback ran, and
# serves the fixture classification. Returns the environment the Fisher stub
# writes to.
wire_multi_ora_kegg <- function(enrich, env = parent.frame()) {
    seen <- new.env()
    seen$fisher <- FALSE
    testthat::local_mocked_bindings(enrichKEGG = enrich,
                                    .package = "clusterProfiler", .env = env)
    local_kegg_ora_stubs(list(
        run_ora_kegg_fisher = function(...) {
            seen$fisher <- TRUE
            data.frame(pathway = "Fisher result", ID = "00190", pvalue = 0.001,
                       padj = 0.01, stringsAsFactors = FALSE)
        },
        kegg_pathway_categories = function(...) parse_kegg_brite_pathways(brite_fixture())
    ), env = env)
    seen
}

run_kegg_producer <- function(exclude_classes = NULL) {
    msgs <- character(0)
    res <- withCallingHandlers(
        run_multi_ora_kegg(c("hsa:1", "hsa:2", "hsa:3"), paste0("hsa:", 1:50),
                           kegg_org = "hsa", label = "pooled", pval_cutoff = 0.1,
                           exclude_classes = exclude_classes),
        message = function(m) {
            msgs <<- c(msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
        })
    list(res = res, msgs = msgs)
}

test_that("adjusted-p hits print the completion message and skip Fisher", {
    skip_if_not_installed("clusterProfiler")
    seen <- wire_multi_ora_kegg(function(...) fake_enrichkegg_table())

    out <- run_kegg_producer()

    expect_identical(out$res$ID, c("hsa00010", "hsa04260"))
    expect_equal(out$res$padj, c(0.01, 0.02))
    expect_true(any(grepl("pooled: 2 enriched pathways", out$msgs, fixed = TRUE)))
    expect_false(seen$fisher)
})

test_that("the raw-p fallback prints its own message and the completion message", {
    skip_if_not_installed("clusterProfiler")
    seen <- wire_multi_ora_kegg(function(...)
        fake_enrichkegg_table(p.adjust = c(0.5, 0.6)))

    out <- run_kegg_producer()

    expect_identical(out$res$ID, c("hsa00010", "hsa04260"))
    expect_true(any(grepl("padj too strict", out$msgs, fixed = TRUE)))
    expect_true(any(grepl("pooled: 2 enriched pathways", out$msgs, fixed = TRUE)))
    expect_false(seen$fisher)
})

test_that("hits that class exclusion removes entirely return NULL without Fisher", {
    skip_if_not_installed("clusterProfiler")
    # clusterProfiler found something; the project excluded all of it. That is
    # a finished answer, not a missing one, so the Fisher fallback must not run
    # and test the excluded pathways a second way.
    seen <- wire_multi_ora_kegg(function(...)
        fake_enrichkegg_table(ids = "hsa04260", pvalue = 0.001, p.adjust = 0.01))

    out <- run_kegg_producer(exclude_classes = "Organismal Systems")

    expect_null(out$res)
    expect_false(seen$fisher)
    expect_false(any(grepl("enriched pathways", out$msgs, fixed = TRUE)))
})

test_that("no clusterProfiler hits still fall back to Fisher", {
    skip_if_not_installed("clusterProfiler")
    seen <- wire_multi_ora_kegg(function(...)
        fake_enrichkegg_table(pvalue = c(0.5, 0.6), p.adjust = c(0.8, 0.9)))

    out <- run_kegg_producer()

    expect_true(seen$fisher)
    expect_identical(out$res$pathway, "Fisher result")
    expect_false(any(grepl("pooled: [0-9]+ enriched pathways", out$msgs)))
})

test_that("a clusterProfiler error still falls back to Fisher", {
    skip_if_not_installed("clusterProfiler")
    seen <- wire_multi_ora_kegg(function(...) stop("synthetic enrichKEGG failure"))

    out <- run_kegg_producer()

    expect_true(seen$fisher)
    expect_identical(out$res$pathway, "Fisher result")
    expect_true(any(grepl("clusterProfiler ORA failed: synthetic enrichKEGG failure",
                          out$msgs, fixed = TRUE)))
    expect_false(any(grepl("pooled: [0-9]+ enriched pathways", out$msgs)))
})
