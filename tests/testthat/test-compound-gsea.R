# Compound GSEA for metabolomics: what it ranks on, what it scores, and what it
# is deliberately kept away from.
#
# GSEA runs beside compound ORA, not instead of it, and its rows never enter
# `pathway_tables`, which drive the per-layer figures and the ORA figure. It
# reaches the cross-omics meta-analysis through its own `rank_tables` argument,
# where merge_pathway_pvalues() picks one method per layer instead of taking
# the minimum across methods (see test-cross-omics-method-selection.R).
#
# Nothing here touches the network. The compound-pathway table is injected, and
# the one function that reads it from disk reads a cache and never fetches.
#
# All fixtures are synthetic.

cpd_fixture <- function() {
    data.frame(
        pathway = c(rep("map00010", 6), rep("map00020", 5), rep("map00030", 4)),
        compound = c(sprintf("C%05d", 1:6),
                     sprintf("C%05d", 7:11),
                     sprintf("C%05d", 12:15)),
        name = c(rep("Glycolysis", 6), rep("Citrate cycle", 5),
                 rep("Pentose phosphate", 4)),
        stringsAsFactors = FALSE
    )
}

# Fifteen compounds, ranked so that map00010's members sit at the top: enough
# structure for fgsea to score, small enough to stay fast.
de_fixture <- function(n = 15, statistic = TRUE) {
    df <- data.frame(
        KEGG_ID = sprintf("C%05d", seq_len(n)),
        log2fc  = seq(2, -2, length.out = n),
        pvalue  = seq(0.001, 0.5, length.out = n),
        stringsAsFactors = FALSE
    )
    if (statistic) df$statistic <- seq(5, -5, length.out = n)
    df
}


# ---- the standardiser carries the moderated statistic ----------------------

metab_de_data <- function(method = "limma", statistic = c(4.2, -2.1, 0.4),
                          stat_col = "statistic") {
    tbl <- data.frame(
        feature_id = c("M1", "M2", "M3"),
        logFC      = c(1.5, -0.8, 0.2),
        P.Value    = c(0.001, 0.02, 0.5),
        adj.P.Val  = c(0.01, 0.08, 0.6),
        stringsAsFactors = FALSE
    )
    if (!is.null(statistic)) tbl[[stat_col]] <- statistic
    out <- list(de_tables = list(B_vs_A = tbl))
    if (!is.null(method)) out$method <- method
    out
}

test_that("a signed statistic survives extract_de_tables() unchanged", {
    std <- extract_de_tables(metab_de_data("limma"), "metabolomics",
                             harmonization_res = NULL)

    expect_true("statistic" %in% names(std$B_vs_A))
    expect_equal(std$B_vs_A$statistic, c(4.2, -2.1, 0.4))

    # Additive: what the standardiser produced before is untouched.
    expect_equal(std$B_vs_A$feature_id, c("M1", "M2", "M3"))
    expect_equal(std$B_vs_A$log2fc, c(1.5, -0.8, 0.2))
    expect_equal(std$B_vs_A$pvalue, c(0.001, 0.02, 0.5))
    expect_equal(std$B_vs_A$padj, c(0.01, 0.08, 0.6))
})

test_that("every t-based method carries its statistic through", {
    for (m in c("limma", "t_test", "t_test_equal")) {
        std <- extract_de_tables(metab_de_data(m), "metabolomics", NULL)
        expect_equal(std$B_vs_A$statistic, c(4.2, -2.1, 0.4),
                     info = paste("method:", m))
    }
})

test_that("a Wilcoxon result never exposes W as the statistic", {
    # wilcox.test() stores W in the same column the t-based methods use for a
    # signed t. W is non-negative and centred at n1*n2/2, so ranking on it would
    # give fgsea a list with no low end at all -- decreased metabolites weighted
    # slightly positive instead of strongly negative. It is not converted here,
    # it is withheld, and the ranker's signed fallback takes over.
    w_stats <- c(30, 2, 18)   # plausible rank sums: all non-negative
    std <- extract_de_tables(metab_de_data("wilcoxon", statistic = w_stats),
                             "metabolomics", NULL)

    expect_true("statistic" %in% names(std$B_vs_A))
    expect_true(all(is.na(std$B_vs_A$statistic)))
    # The rest of the row is untouched -- only the statistic is withheld.
    expect_equal(std$B_vs_A$log2fc, c(1.5, -0.8, 0.2))
    expect_equal(std$B_vs_A$pvalue, c(0.001, 0.02, 0.5))
})

test_that("a Wilcoxon DE result ranks on the signed fallback", {
    std <- extract_de_tables(
        metab_de_data("wilcoxon", statistic = c(30, 2, 18)), "metabolomics", NULL)
    de_mapped <- std$B_vs_A
    de_mapped$KEGG_ID <- c("C00001", "C00002", "C00003")

    ranks <- rank_compounds_for_gsea(de_mapped)

    # Signed by the fold change, not ordered by W: M2 fell, so its compound must
    # rank negative even though its W is the smallest of the three.
    expect_lt(ranks[["C00002"]], 0)
    expect_gt(ranks[["C00001"]], 0)
    expect_gt(ranks[["C00003"]], 0)
    # And the ordering is the fallback's, not W's ascending order.
    expect_equal(unname(ranks),
                 unname(sort(sign(c(1.5, -0.8, 0.2)) *
                             -log10(c(0.001, 0.02, 0.5) + 1e-300),
                             decreasing = TRUE)))
})

test_that("an unknown, missing or precomputed method withholds the statistic", {
    for (m in list("precomputed", "something_new", NULL)) {
        std <- extract_de_tables(metab_de_data(m), "metabolomics", NULL)
        expect_true(all(is.na(std$B_vs_A$statistic)),
                    info = paste("method:", m %||% "<missing>"))
    }
})

test_that("a table carrying t is accepted, and one carrying neither gives NA", {
    with_t <- metab_de_data("limma", statistic = c(3.3, -1.1, 0.2),
                            stat_col = "t")
    expect_equal(extract_de_tables(with_t, "metabolomics", NULL)$B_vs_A$statistic,
                 c(3.3, -1.1, 0.2))

    without <- metab_de_data("limma", statistic = NULL)
    got <- extract_de_tables(without, "metabolomics", NULL)$B_vs_A
    expect_true("statistic" %in% names(got))
    expect_true(all(is.na(got$statistic)))
})


# ---- what the ranking is built from ----------------------------------------

test_that("the moderated statistic is preferred when it carries usable values", {
    # Statistic and fallback disagree on sign for every row, so the ranking can
    # only match one of them.
    de <- data.frame(
        KEGG_ID   = c("C00001", "C00002"),
        statistic = c(3, -3),
        log2fc    = c(-1, 1),
        pvalue    = c(0.001, 0.001),
        stringsAsFactors = FALSE
    )

    ranks <- rank_compounds_for_gsea(de)

    expect_identical(names(ranks), c("C00001", "C00002"))
    expect_equal(unname(ranks), c(3, -3))
})

test_that("the fallback is used when the statistic column has no finite values", {
    # The column is present and full of NA, which is exactly what the
    # standardiser produces for a DE export that never had one. Testing the
    # column rather than its values would take this all-NA vector.
    de <- data.frame(
        KEGG_ID   = c("C00001", "C00002"),
        statistic = c(NA_real_, NA_real_),
        log2fc    = c(1, -1),
        pvalue    = c(0.01, 0.01),
        stringsAsFactors = FALSE
    )

    ranks <- rank_compounds_for_gsea(de)

    expect_equal(unname(ranks), c(2, -2))
})

test_that("the fallback is used when there is no statistic column at all", {
    de <- data.frame(KEGG_ID = c("C00001", "C00002"),
                     log2fc = c(1, -1), pvalue = c(0.01, 0.01),
                     stringsAsFactors = FALSE)

    expect_equal(unname(rank_compounds_for_gsea(de)), c(2, -2))
})

test_that("a p-value that underflowed to zero ranks high rather than infinite", {
    de <- data.frame(KEGG_ID = c("C00001", "C00002"),
                     log2fc = c(1, 1), pvalue = c(0, 0.01),
                     stringsAsFactors = FALSE)

    ranks <- rank_compounds_for_gsea(de)

    expect_true(all(is.finite(ranks)))
    expect_gt(ranks[["C00001"]], ranks[["C00002"]])
})

test_that("a zero fold change ranks neutrally instead of being dropped", {
    de <- data.frame(KEGG_ID = c("C00001", "C00002", "C00003"),
                     log2fc = c(1, 0, -1), pvalue = c(0.01, 0.01, 0.01),
                     stringsAsFactors = FALSE)

    ranks <- rank_compounds_for_gsea(de)

    expect_length(ranks, 3)
    expect_equal(unname(ranks[["C00002"]]), 0)
})

test_that("non-finite ranks and unusable ids are dropped", {
    de <- data.frame(
        KEGG_ID = c("C00001", NA, "", "C00004"),
        log2fc  = c(1, 1, 1, NA),
        pvalue  = c(0.01, 0.01, 0.01, 0.01),
        stringsAsFactors = FALSE
    )

    expect_identical(names(rank_compounds_for_gsea(de)), "C00001")
})

test_that("nothing rankable gives an empty vector rather than an error", {
    expect_length(rank_compounds_for_gsea(data.frame()), 0L)
    expect_length(rank_compounds_for_gsea(NULL), 0L)
    expect_length(rank_compounds_for_gsea(
        data.frame(log2fc = 1, pvalue = 0.01)), 0L)   # no KEGG_ID
})


# ---- duplicate compounds ---------------------------------------------------

test_that("one compound appears once, carrying its strongest rank", {
    # Three metabolites annotate to two compounds. fgsea would score C00001
    # twice without this.
    de <- data.frame(
        KEGG_ID   = c("C00001", "C00001", "C00002"),
        statistic = c(1.5, -4.0, 2.0),
        log2fc    = c(1, -1, 1),
        pvalue    = c(0.05, 0.001, 0.01),
        stringsAsFactors = FALSE
    )

    ranks <- rank_compounds_for_gsea(de)

    expect_length(ranks, 2L)
    # -4.0 beats 1.5 on absolute rank, and keeps its sign.
    expect_equal(unname(ranks[["C00001"]]), -4)
})

test_that("the collapse does not depend on the row order it arrives in", {
    de <- data.frame(
        KEGG_ID   = c("C00001", "C00001", "C00002"),
        statistic = c(1.5, -4.0, 2.0),
        log2fc    = c(1, -1, 1),
        pvalue    = c(0.05, 0.001, 0.01),
        stringsAsFactors = FALSE
    )

    forwards  <- rank_compounds_for_gsea(de)
    backwards <- rank_compounds_for_gsea(de[rev(seq_len(nrow(de))), , drop = FALSE])

    expect_identical(forwards, backwards)
})

test_that("two rows of ONE compound tying on magnitude resolve the same either way", {
    # The case the id cannot separate: same compound, same absolute rank,
    # opposite signs. Magnitude and name are both tied, so without a third key
    # order() falls back to the arrival index and the answer flips with the
    # input. Which sign wins is arbitrary; that it does not depend on row order
    # is the point.
    de <- data.frame(
        KEGG_ID   = c("C00001", "C00001"),
        statistic = c(2, -2),
        log2fc    = c(1, -1),
        pvalue    = c(0.01, 0.01),
        stringsAsFactors = FALSE
    )

    forwards  <- rank_compounds_for_gsea(de)
    backwards <- rank_compounds_for_gsea(de[c(2, 1), , drop = FALSE])

    expect_length(forwards, 1L)
    expect_identical(forwards, backwards)
})

test_that("compounds tying on absolute rank resolve deterministically", {
    de <- data.frame(
        KEGG_ID   = c("C00002", "C00001"),
        statistic = c(2, -2),
        log2fc    = c(1, -1),
        pvalue    = c(0.01, 0.01),
        stringsAsFactors = FALSE
    )

    expect_identical(rank_compounds_for_gsea(de),
                     rank_compounds_for_gsea(de[c(2, 1), , drop = FALSE]))
})


# ---- reading the cache, never fetching -------------------------------------

test_that("an absent or unusable cache yields NULL rather than a download", {
    empty <- withr::local_tempdir()
    expect_null(.cached_compound_pathways(empty))
    expect_null(.cached_compound_pathways(NULL))

    # Present but not the expected object.
    saveRDS(list(a = 1), file.path(empty, "kegg_compound_pathways.rds"))
    expect_null(.cached_compound_pathways(empty))
})

test_that("a populated cache is read back", {
    dir <- withr::local_tempdir()
    saveRDS(cpd_fixture(), file.path(dir, "kegg_compound_pathways.rds"))

    got <- .cached_compound_pathways(dir)

    expect_true(is.data.frame(got))
    expect_true(all(c("compound", "pathway") %in% names(got)))
})

test_that("compound GSEA declines to run without a compound-pathway table", {
    dir <- withr::local_tempdir()   # deliberately empty

    expect_message(
        res <- run_compound_gsea(de_fixture(), cache_dir = dir,
                                 min_gs = 2, max_gs = 500),
        "compound-pathway"
    )
    expect_null(res)
})

test_that("the orchestrator reaches for the cache and never the downloader", {
    # The guarantee is structural, not conditional: there is no branch in this
    # path that could call get_kegg_compound_pathways(), whose own cold-cache
    # behaviour is to fetch from KEGG.
    src <- paste(deparse(body(run_compound_gsea_for_contrasts)), collapse = " ")

    expect_true(grepl(".cached_compound_pathways(", src, fixed = TRUE))
    expect_false(grepl("get_kegg_compound_pathways(", src, fixed = TRUE))
})


# ---- scoring ---------------------------------------------------------------

test_that("a scored table carries the agreed schema and names its method", {
    skip_if_not_installed("fgsea")

    res <- suppressMessages(run_compound_gsea(
        de_fixture(), cache_dir = NULL, min_gs = 2, max_gs = 500,
        cpd_pathways = cpd_fixture()))

    expect_true(is.data.frame(res))
    expect_true(all(c("pathway", "ID", "pvalue", "padj", "NES", "ES",
                      "setSize", "leadingEdge", "database", "method")
                    %in% names(res)))
    expect_identical(unique(res$method), "fgsea")
    expect_identical(unique(res$database), "KEGG")
    # Accession in ID, readable name in pathway -- the way pathway_join_key()
    # and pathway_display_label() each read first.
    expect_true(all(grepl("^map[0-9]{5}$", res$ID)))
    expect_false(any(grepl("^map[0-9]{5}$", res$pathway)))
})

test_that("no significance cutoff is applied inside the producer", {
    skip_if_not_installed("fgsea")

    res <- suppressMessages(run_compound_gsea(
        de_fixture(), cache_dir = NULL, min_gs = 2, max_gs = 500,
        cpd_pathways = cpd_fixture()))

    # Every pathway large enough to score comes back, whatever its p-value --
    # which is the whole claim. Asserting that some row is non-significant would
    # be asserting a property of the fixture's numbers instead.
    expect_setequal(res$ID, c("map00010", "map00020", "map00030"))
})

test_that("the same seed gives the same scores", {
    skip_if_not_installed("fgsea")

    once  <- suppressMessages(run_compound_gsea(
        de_fixture(), cache_dir = NULL, min_gs = 2, max_gs = 500, seed = 42,
        cpd_pathways = cpd_fixture()))
    twice <- suppressMessages(run_compound_gsea(
        de_fixture(), cache_dir = NULL, min_gs = 2, max_gs = 500, seed = 42,
        cpd_pathways = cpd_fixture()))

    expect_equal(once$pvalue, twice$pvalue)
    expect_equal(once$NES, twice$NES)
    expect_identical(once$ID, twice$ID)
})

test_that("pathway membership comes from the supplied compound mapping", {
    skip_if_not_installed("fgsea")

    # One pathway removed from the mapping cannot be scored, and nothing else
    # changes: membership is not re-derived from the ranked list.
    partial <- cpd_fixture()
    partial <- partial[partial$pathway != "map00030", , drop = FALSE]

    res <- suppressMessages(run_compound_gsea(
        de_fixture(), cache_dir = NULL, min_gs = 2, max_gs = 500,
        cpd_pathways = partial))

    expect_false("map00030" %in% res$ID)
    expect_true(all(c("map00010", "map00020") %in% res$ID))
})

test_that("GSEA rows are refused by the ORA figure's membership rule", {
    skip_if_not_installed("fgsea")

    res <- suppressMessages(run_compound_gsea(
        de_fixture(), cache_dir = NULL, min_gs = 2, max_gs = 500,
        cpd_pathways = cpd_fixture()))

    m <- build_ora_adjusted_p_matrix(list(metabolomics = res), "00010",
                                     "metabolomics")

    # #208's figure shows ORA evidence; an fgsea row must not supply a cell.
    expect_true(is.na(m["00010", "metabolomics"]))
})


# ---- compound ORA is untouched ---------------------------------------------

test_that("scoring GSEA leaves the compound ORA result identical", {
    skip_if_not_installed("fgsea")

    cache <- withr::local_tempdir()
    saveRDS(cpd_fixture(), file.path(cache, "kegg_compound_pathways.rds"))

    universe <- sprintf("C%05d", 1:15)
    de_mapped <- data.frame(
        KEGG_ID = universe,
        pvalue  = c(rep(1e-4, 4), rep(0.5, 11)),
        padj    = c(rep(0.01, 4), rep(0.9, 11)),
        log2fc  = seq(2, -2, length.out = 15),
        statistic = seq(5, -5, length.out = 15),
        stringsAsFactors = FALSE
    )

    before <- suppressMessages(run_compound_ora(
        de_mapped, cache_dir = cache, min_gs = 2, max_gs = 500,
        pval_cutoff = 1, universe = universe))

    suppressMessages(run_compound_gsea(de_mapped, cache_dir = cache,
                                       min_gs = 2, max_gs = 500,
                                       cpd_pathways = cpd_fixture()))

    after <- suppressMessages(run_compound_ora(
        de_mapped, cache_dir = cache, min_gs = 2, max_gs = 500,
        pval_cutoff = 1, universe = universe))

    expect_identical(before, after)
})


# ---- degenerate inputs -----------------------------------------------------

test_that("too few rankable compounds gives NULL", {
    expect_message(
        res <- run_compound_gsea(de_fixture(n = 2), cache_dir = NULL,
                                 min_gs = 2, max_gs = 500,
                                 cpd_pathways = cpd_fixture()),
        "Too few"
    )
    expect_null(res)
})

test_that("no finite ranking value anywhere gives NULL", {
    de <- data.frame(KEGG_ID = sprintf("C%05d", 1:5),
                     statistic = NA_real_, log2fc = NA_real_,
                     pvalue = NA_real_, stringsAsFactors = FALSE)

    expect_message(
        res <- run_compound_gsea(de, cache_dir = NULL, min_gs = 2, max_gs = 500,
                                 cpd_pathways = cpd_fixture()),
        "Too few"
    )
    expect_null(res)
})

test_that("no pathway inside the size bounds gives NULL", {
    # Exercised through the sets rather than through min_gs: the floor is
    # clamped to at most 3 for compound pathways (see below), so a large
    # configured min_gs no longer excludes anything. One measured compound per
    # pathway is under the floor whatever the configuration says.
    thin <- data.frame(
        pathway  = c("map00010", "map00020"),
        compound = c("C00001", "C00002"),
        name     = c("Glycolysis", "Citrate cycle"),
        stringsAsFactors = FALSE
    )

    expect_message(
        res <- run_compound_gsea(de_fixture(), cache_dir = NULL,
                                 min_gs = 2, max_gs = 500,
                                 cpd_pathways = thin),
        "measured compounds"
    )
    expect_null(res)
})

test_that("a configured gene-set-scale floor does not disqualify compound pathways", {
    skip_if_not_installed("fgsea")

    # min_set_size ships at 10, which is modest for a gene set and far above
    # what a KEGG compound pathway carries in measured members. Clamping up to
    # it would leave compound GSEA testing nothing while reporting no error, so
    # the floor is clamped down to 3, exactly as run_compound_ora() does.
    three <- data.frame(
        pathway  = rep(c("map00010", "map00020"), each = 3),
        compound = c("C00001", "C00002", "C00003",
                     "C00004", "C00005", "C00006"),
        name     = rep(c("Glycolysis", "Citrate cycle"), each = 3),
        stringsAsFactors = FALSE
    )

    res <- suppressMessages(run_compound_gsea(
        de_fixture(), cache_dir = NULL, min_gs = 10, max_gs = 500,
        cpd_pathways = three))

    expect_true(is.data.frame(res))
    expect_setequal(res$ID, c("map00010", "map00020"))
})

test_that("the configured maximum is still honoured as given", {
    # Only the floor is compound-specific; max_gs is used as configured.
    expect_message(
        res <- run_compound_gsea(de_fixture(), cache_dir = NULL,
                                 min_gs = 2, max_gs = 3,
                                 cpd_pathways = cpd_fixture()),
        "measured compounds"
    )
    expect_null(res)
})

test_that("a single scorable pathway is not treated as a failure", {
    skip_if_not_installed("fgsea")

    only_one <- cpd_fixture()
    only_one <- only_one[only_one$pathway == "map00010", , drop = FALSE]

    res <- suppressMessages(run_compound_gsea(
        de_fixture(), cache_dir = NULL, min_gs = 2, max_gs = 500,
        cpd_pathways = only_one))

    expect_equal(nrow(res), 1L)
    expect_identical(res$ID, "map00010")
})

test_that("the orchestrator gives NULL when metabolomics has no DE results", {
    expect_null(run_compound_gsea_for_contrasts(
        de_results = list(transcriptomics = list()),
        harmonization_res = NULL,
        config = list(),
        out_dir = withr::local_tempdir()))
})


# ---- the export belongs to the run that produced it -------------------------

test_that("a run with GSEA results writes the export", {
    dir <- withr::local_tempdir()
    res <- data.frame(ID = "map00010", pathway = "Glycolysis", NES = 1.5,
                      pvalue = 0.01, padj = 0.03, method = "fgsea",
                      stringsAsFactors = FALSE)

    path <- write_compound_gsea_export(res, dir)

    expect_false(is.null(path))
    expect_true(file.exists(file.path(dir, "metabolomics_compound_gsea.csv")))
    expect_equal(nrow(read.csv(path, stringsAsFactors = FALSE)), 1L)
})

test_that("a run with no GSEA results removes the previous run's export", {
    # The report includes this file on file.exists() alone, so a leftover would
    # be read as this run's evidence. A run can legitimately score nothing --
    # no metabolomics DE, no mapping, a cold cache -- and that must leave no
    # file behind rather than the last run's.
    dir <- withr::local_tempdir()
    stale <- file.path(dir, "metabolomics_compound_gsea.csv")
    writeLines("left over from an earlier run", stale)

    expect_null(write_compound_gsea_export(NULL, dir))
    expect_false(file.exists(stale))
})

test_that("an empty result also clears rather than keeps", {
    dir <- withr::local_tempdir()
    stale <- file.path(dir, "metabolomics_compound_gsea.csv")
    writeLines("left over", stale)

    expect_null(write_compound_gsea_export(
        data.frame(ID = character(0), stringsAsFactors = FALSE), dir))
    expect_false(file.exists(stale))
})

test_that("writing replaces an earlier export rather than appending to it", {
    dir <- withr::local_tempdir()
    first  <- data.frame(ID = c("map00010", "map00020"), NES = c(1, -1),
                         stringsAsFactors = FALSE)
    second <- data.frame(ID = "map00030", NES = 2, stringsAsFactors = FALSE)

    write_compound_gsea_export(first, dir)
    path <- write_compound_gsea_export(second, dir)

    got <- read.csv(path, stringsAsFactors = FALSE)
    expect_identical(got$ID, "map00030")
})


# ---- class exclusion must not reach the network -----------------------------

brite_fixture <- function() {
    data.frame(
        pathway_id   = c("00010", "00020", "00030"),
        category     = c("Metabolism", "Metabolism", "Human Diseases"),
        subcategory  = c("Carbohydrate metabolism", "Carbohydrate metabolism",
                         "Cancer"),
        pathway_name = c("Glycolysis", "Citrate cycle", "Pentose phosphate"),
        stringsAsFactors = FALSE
    )
}

test_that("the BRITE classification is read from cache and never fetched", {
    dir <- withr::local_tempdir()
    expect_null(.cached_pathway_categories(dir))

    saveRDS(brite_fixture(), file.path(dir, "kegg_pathway_categories.rds"))
    got <- .cached_pathway_categories(dir)

    expect_true(is.data.frame(got))
    expect_true(all(c("pathway_id", "category", "subcategory") %in% names(got)))
})

test_that("an unusable BRITE cache is refused rather than half-read", {
    dir <- withr::local_tempdir()
    saveRDS(data.frame(nonsense = 1), file.path(dir, "kegg_pathway_categories.rds"))

    # Validated with the producer's own shape check, so a file that is not the
    # classification table reads as absent and the caller fails open.
    expect_null(.cached_pathway_categories(dir))
})

test_that("compound GSEA passes a cached classification rather than the fetching default", {
    # keep_kegg_pathways() defaults `classification` to kegg_pathway_categories(),
    # which downloads br08901 on a cold cache. Leaving that default in place
    # would make this path fetch indirectly whenever exclude_pathway_classes is
    # set -- the one hole in the cache-only guarantee.
    src <- paste(deparse(body(run_compound_gsea)), collapse = " ")

    expect_true(grepl("classification = .cached_pathway_categories(", src,
                      fixed = TRUE))
})

test_that("class exclusion works from cache alone, and fails open without one", {
    skip_if_not_installed("fgsea")

    cache <- withr::local_tempdir()
    saveRDS(brite_fixture(), file.path(cache, "kegg_pathway_categories.rds"))

    excluded <- suppressMessages(run_compound_gsea(
        de_fixture(), cache_dir = cache, min_gs = 2, max_gs = 500,
        exclude_classes = "Human Diseases", cpd_pathways = cpd_fixture()))

    # map00030 is the Human Diseases row of the fixture.
    expect_false("map00030" %in% excluded$ID)
    expect_true("map00010" %in% excluded$ID)

    # With no classification available, nothing is dropped: an exclusion that
    # cannot be resolved must not silently remove pathways.
    empty <- withr::local_tempdir()
    kept <- suppressMessages(run_compound_gsea(
        de_fixture(), cache_dir = empty, min_gs = 2, max_gs = 500,
        exclude_classes = "Human Diseases", cpd_pathways = cpd_fixture()))

    expect_true("map00030" %in% kept$ID)
})


# ---- GSEA survives a run where no ORA table was produced --------------------

test_that("compound GSEA still runs, returns and exports when per_omics is empty", {
    skip_if_not_installed("fgsea")

    tmp <- withr::local_tempdir()
    out_dir <- file.path(tmp, "cross_enrichment")
    dir.create(file.path(out_dir, "metabolomics"), recursive = TRUE)

    # The compound-pathway table the ORA step would have cached.
    saveRDS(cpd_fixture(),
            file.path(out_dir, "metabolomics", "kegg_compound_pathways.rds"))

    # HMDB -> KEGG mapping, seeded through the env var load_hmdb_to_kegg_map()
    # already honours, so nothing reaches for the repo's real data file.
    hmdb <- sprintf("HMDB%07d", 1:15)
    dir.create(file.path(tmp, "data"))
    writeLines(c("HMDB\tKEGG", paste0(hmdb, "\t", sprintf("C%05d", 1:15))),
               file.path(tmp, "data", "HMDB2kegg_cpd.Jan2026.v2.txt"))
    withr::local_envvar(c(PIPELINE_ROOT = tmp))

    harmonization_res <- list(inputs = list(metabolomics = list(
        row_data = data.frame(feature_id = hmdb, HMDB = hmdb,
                              stringsAsFactors = FALSE))))

    de_results <- list(metabolomics = list(de_tables = list(
        B_vs_A = data.frame(
            feature_id = hmdb,
            logFC      = seq(2, -2, length.out = 15),
            P.Value    = seq(0.001, 0.5, length.out = 15),
            adj.P.Val  = seq(0.01, 0.9, length.out = 15),
            statistic  = seq(5, -5, length.out = 15),
            stringsAsFactors = FALSE
        ))))

    # No omics layer is enriched: omics_present is empty, so
    # build_per_omics_enrichment() has nothing to loop over and per_omics comes
    # back empty. An organism with no KEGG code keeps the gene-side conversion
    # cache -- and its network calls -- out of the picture entirely.
    # Deliberately no min_set_size: this runs on the shipped default of 10, so
    # it also stands as the end-to-end check that a gene-set-scale floor does
    # not silence compound GSEA in a default configuration.
    config <- list(
        global = list(organism = "not a known organism",
                      omics_present = character(0)),
        modes = list(multiomics = list(enrichment = list(run_enrichment = TRUE)))
    )

    res <- suppressWarnings(suppressMessages(mod_multiomics_enrichment(
        enrichment_results = NULL,
        de_results = de_results,
        harmonization_res = harmonization_res,
        config = config,
        out_dir = out_dir)))

    # Returned rather than thrown away with the empty ORA result.
    expect_false(is.null(res))
    expect_true(is.data.frame(res$compound_gsea))
    expect_gt(nrow(res$compound_gsea), 0)
    expect_identical(unique(res$compound_gsea$method), "fgsea")

    # Reported as what it is: no ORA tables, no cross-omics analysis, no figures.
    expect_identical(res$per_omics, list())
    expect_null(res$cross_omics)
    expect_identical(res$plots, list())

    # And exported.
    expect_true(file.exists(file.path(out_dir, "metabolomics_compound_gsea.csv")))
})

test_that("a run with neither ORA nor GSEA still returns NULL", {
    tmp <- withr::local_tempdir()
    out_dir <- file.path(tmp, "cross_enrichment")
    dir.create(out_dir, recursive = TRUE)

    # No compound-pathway cache, so GSEA cannot score; no layers either.
    config <- list(
        global = list(organism = "not a known organism",
                      omics_present = character(0)),
        modes = list(multiomics = list(enrichment = list(run_enrichment = TRUE)))
    )

    res <- suppressWarnings(suppressMessages(mod_multiomics_enrichment(
        enrichment_results = NULL,
        de_results = list(),
        harmonization_res = NULL,
        config = config,
        out_dir = out_dir)))

    # Unchanged from before this PR: nothing to report is still NULL.
    expect_null(res)
})

test_that("the GSEA call precedes the empty-per_omics guard", {
    # Ordering is the whole fix: scored after the guard, it would be unreachable
    # in exactly the case the test above covers. Pinned because the two are far
    # apart in the function and easy to separate again.
    src <- paste(deparse(body(mod_multiomics_enrichment)), collapse = " ")
    gsea_at  <- regexpr("run_compound_gsea_for_contrasts(", src, fixed = TRUE)
    guard_at <- regexpr("No omics layers produced enrichment results", src,
                        fixed = TRUE)

    expect_gt(gsea_at, 0)
    expect_gt(guard_at, 0)
    expect_lt(gsea_at, guard_at)
})


test_that("disabling enrichment clears a previous run's GSEA export", {
    # The disabled guard returns before any of the enrichment work, so a rerun
    # with enrichment switched off would have left the earlier export for the
    # report to show as current. Cleared ahead of every return, not just the
    # ones that reach the scoring.
    tmp <- withr::local_tempdir()
    out_dir <- file.path(tmp, "cross_enrichment")
    dir.create(out_dir, recursive = TRUE)
    stale <- file.path(out_dir, "metabolomics_compound_gsea.csv")
    writeLines("left over from an earlier run", stale)

    config <- list(
        global = list(organism = "not a known organism",
                      omics_present = character(0)),
        modes = list(multiomics = list(enrichment = list(run_enrichment = FALSE)))
    )

    res <- suppressMessages(mod_multiomics_enrichment(
        enrichment_results = NULL, de_results = list(),
        harmonization_res = NULL, config = config, out_dir = out_dir))

    expect_null(res)
    expect_false(file.exists(stale))
})

test_that("the export cleanup precedes the disabled-enrichment guard", {
    # Ordering is the fix; the two sit close together and are easy to swap back.
    src <- paste(deparse(body(mod_multiomics_enrichment)), collapse = " ")
    clear_at <- regexpr("write_compound_gsea_export(NULL", src, fixed = TRUE)
    guard_at <- regexpr("Cross-omics enrichment disabled in config", src,
                        fixed = TRUE)

    expect_gt(clear_at, 0)
    expect_gt(guard_at, 0)
    expect_lt(clear_at, guard_at)
})
