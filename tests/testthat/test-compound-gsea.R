# tests/testthat/test-compound-gsea.R
#
# Unit tests for the compound-GSEA companion to the compound ORA in
# R/domain/multiomics/07_enrichment.R:
#   - rank_compounds_for_gsea()
#   - run_compound_gsea()
#
# All data is synthetic: made-up compound accessions, made-up map ids and
# made-up statistics. The two KEGG caches run_compound_gsea() consults are
# seeded on disk in a temp directory, so no test here touches the network.

# --- fixtures ----------------------------------------------------------------

# A DE table in the shape extract_de_tables() emits for metabolomics, already
# merged with KEGG compound ids: feature_id, log2fc, pvalue, padj, statistic,
# KEGG_ID.
make_compound_de <- function(n = 40, with_statistic = TRUE, seed = 7L) {
    withr::with_seed(seed, {
        stat <- stats::rnorm(n, sd = 2)
        data.frame(
            feature_id = paste0("feature_", seq_len(n)),
            log2fc     = stat / 4,
            pvalue     = 2 * stats::pnorm(-abs(stat)),
            padj       = pmin(1, 2 * stats::pnorm(-abs(stat)) * 2),
            statistic  = if (with_statistic) stat else NA_real_,
            KEGG_ID    = sprintf("C%05d", seq_len(n)),
            stringsAsFactors = FALSE
        )
    })
}

# Seed the two RDS caches run_compound_gsea() reads, so the KEGG REST calls
# never fire. Two synthetic maps overlap the measured compounds; one is filed
# under a class the exclusion test drops.
seed_kegg_caches <- function(dir, n_cpd = 40) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    cpds <- sprintf("C%05d", seq_len(n_cpd))

    links <- rbind(
        data.frame(pathway = "map90001", compound = cpds[1:12],
                   name = "synthetic core metabolism", stringsAsFactors = FALSE),
        data.frame(pathway = "map90002", compound = cpds[20:31],
                   name = "synthetic side metabolism", stringsAsFactors = FALSE),
        # Only two measured compounds: below any minimum this suite asks for.
        data.frame(pathway = "map90003", compound = cpds[1:2],
                   name = "synthetic tiny map", stringsAsFactors = FALSE)
    )
    saveRDS(links, file.path(dir, "kegg_compound_pathways.rds"))

    saveRDS(
        data.frame(
            pathway_id   = c("90001", "90002", "90003"),
            category     = c("Metabolism", "Synthetic Excluded Class", "Metabolism"),
            subcategory  = c("Carbohydrate metabolism", "Made-up subcategory",
                             "Carbohydrate metabolism"),
            pathway_name = c("synthetic core metabolism",
                             "synthetic side metabolism", "synthetic tiny map"),
            stringsAsFactors = FALSE
        ),
        file.path(dir, "kegg_pathway_categories.rds")
    )
    dir
}

# --- rank_compounds_for_gsea() -----------------------------------------------

test_that("ranking prefers the moderated t when it carries usable values", {
    de <- make_compound_de()
    ranks <- rank_compounds_for_gsea(de)

    expect_length(ranks, nrow(de))
    expect_equal(unname(ranks[de$KEGG_ID[1]]), de$statistic[1])
    expect_false(is.unsorted(rev(ranks)))
})

test_that("an all-NA statistic column falls back to sign(log2FC) * -log10(p)", {
    # The failure this guards: the standardised metabolomics DE table always
    # carries a `statistic` column, so a presence test would pick an all-NA
    # vector and hand fgsea nothing to score.
    de <- make_compound_de(with_statistic = FALSE)
    ranks <- rank_compounds_for_gsea(de)

    expect_length(ranks, nrow(de))
    expect_true(all(is.finite(ranks)))
    expected <- sign(de$log2fc) * -log10(de$pvalue + 1e-300)
    expect_equal(unname(ranks[de$KEGG_ID[1]]), expected[1])
})

test_that("a statistic column with a single usable value is still used", {
    de <- make_compound_de()
    de$statistic <- NA_real_
    de$statistic[3] <- 5.5
    ranks <- rank_compounds_for_gsea(de)

    # Only the one finite entry survives; the rest are dropped as non-finite
    # rather than silently re-derived from p, which would mix two scales.
    expect_length(ranks, 1L)
    expect_equal(unname(ranks), 5.5)
    expect_equal(names(ranks), de$KEGG_ID[3])
})

test_that("compounds annotated twice are collapsed to their strongest rank", {
    de <- make_compound_de(n = 6)
    de$KEGG_ID[2] <- de$KEGG_ID[1]
    de$statistic <- c(1, -4, 2, -2, 3, -3)

    ranks <- rank_compounds_for_gsea(de)

    expect_false(any(duplicated(names(ranks))))
    expect_equal(unname(ranks[de$KEGG_ID[1]]), -4)
})

test_that("collapsing does not depend on the DE table's row order", {
    de <- make_compound_de(n = 8)
    de$KEGG_ID[c(2, 5)] <- de$KEGG_ID[1]
    shuffled <- de[rev(seq_len(nrow(de))), , drop = FALSE]

    expect_equal(rank_compounds_for_gsea(de), rank_compounds_for_gsea(shuffled))
})

test_that("unrankable input yields an empty vector rather than an error", {
    expect_length(rank_compounds_for_gsea(NULL), 0L)
    expect_length(rank_compounds_for_gsea(data.frame()), 0L)

    de <- make_compound_de(n = 5)
    de$KEGG_ID <- NULL
    expect_length(rank_compounds_for_gsea(de), 0L)

    all_na <- make_compound_de(n = 5, with_statistic = FALSE)
    all_na$log2fc <- NA_real_
    all_na$pvalue <- NA_real_
    expect_length(rank_compounds_for_gsea(all_na), 0L)
})

# --- run_compound_gsea() -----------------------------------------------------

test_that("compound GSEA emits the fgsea column shape the gene layers emit", {
    skip_if_not_installed("fgsea")
    cache <- seed_kegg_caches(withr::local_tempdir())
    de <- make_compound_de()

    res <- suppressWarnings(
        run_compound_gsea(de, cache, min_gs = 5, max_gs = 500, seed = 1L,
                          out_dir = NULL)
    )

    expect_s3_class(res, "data.frame")
    expect_true(all(c("pathway", "pathway_name", "pval", "padj", "NES",
                      "size", "database", "method") %in% names(res)))
    # `pathway` must hold the map id: the cross-omics merge joins the layers on
    # this column and the gene layers put map ids there.
    expect_true(all(grepl("^map\\d+$", res$pathway)))
    expect_true(all(res$method == "fgsea"))
    expect_true(all(res$database == "KEGG"))
    expect_true(any(is.finite(res$NES)))
})

test_that("reported size counts measured compounds, not KEGG membership", {
    skip_if_not_installed("fgsea")
    # Only the first 20 compounds are measured, so map90001 (12 members) keeps
    # 12 and map90002 (members 20-31) keeps just one.
    cache <- seed_kegg_caches(withr::local_tempdir())
    de <- make_compound_de(n = 20)

    res <- suppressWarnings(
        run_compound_gsea(de, cache, min_gs = 5, max_gs = 500, seed = 1L,
                          out_dir = NULL)
    )

    expect_equal(res$pathway, "map90001")
    expect_equal(res$size, 12L)
})

test_that("the same seed reproduces the same scores", {
    skip_if_not_installed("fgsea")
    cache <- seed_kegg_caches(withr::local_tempdir())
    de <- make_compound_de()

    a <- suppressWarnings(run_compound_gsea(de, cache, 5, 500, seed = 42L, out_dir = NULL))
    b <- suppressWarnings(run_compound_gsea(de, cache, 5, 500, seed = 42L, out_dir = NULL))

    expect_equal(a, b)
})

test_that("excluded KEGG classes never reach the scored table", {
    skip_if_not_installed("fgsea")
    cache <- seed_kegg_caches(withr::local_tempdir())
    de <- make_compound_de()

    with_all <- suppressWarnings(
        run_compound_gsea(de, cache, 5, 500, seed = 1L, out_dir = NULL)
    )
    filtered <- suppressWarnings(
        run_compound_gsea(de, cache, 5, 500, seed = 1L,
                          exclude_classes = "Synthetic Excluded Class",
                          out_dir = NULL)
    )

    expect_true("map90002" %in% with_all$pathway)
    expect_false("map90002" %in% filtered$pathway)
    # Dropped before testing, so the surviving map is corrected against a
    # smaller family and its padj changes accordingly.
    expect_true("map90001" %in% filtered$pathway)
})

test_that("the min set size is honoured as given", {
    skip_if_not_installed("fgsea")
    cache <- seed_kegg_caches(withr::local_tempdir())
    de <- make_compound_de()

    lenient <- suppressWarnings(run_compound_gsea(de, cache, 2, 500, seed = 1L, out_dir = NULL))
    strict  <- suppressWarnings(run_compound_gsea(de, cache, 13, 500, seed = 1L, out_dir = NULL))

    # The 2-compound map is only testable under the lenient minimum.
    expect_true("map90003" %in% lenient$pathway)
    # Both real maps hold 12 measured compounds, so a minimum of 13 leaves none.
    expect_null(strict)
})

test_that("the scored table is written when an out_dir is given", {
    skip_if_not_installed("fgsea")
    cache <- seed_kegg_caches(withr::local_tempdir())
    out <- withr::local_tempdir()
    de <- make_compound_de()

    res <- suppressWarnings(
        run_compound_gsea(de, cache, 5, 500, seed = 1L, out_dir = out,
                          label = "cond_a_vs_cond_b")
    )

    written <- file.path(out, "cond_a_vs_cond_b_compound_gsea_all_tested.csv")
    expect_true(file.exists(written))
    expect_equal(nrow(utils::read.csv(written)), nrow(res))
})

test_that("too few rankable compounds returns NULL, not an error", {
    cache <- seed_kegg_caches(withr::local_tempdir())
    de <- make_compound_de(n = 2)

    expect_null(suppressWarnings(
        run_compound_gsea(de, cache, 5, 500, seed = 1L, out_dir = NULL)
    ))
})

test_that("missing compound-pathway associations return NULL", {
    empty_cache <- withr::local_tempdir()
    saveRDS(data.frame(pathway = character(0), compound = character(0),
                       name = character(0), stringsAsFactors = FALSE),
            file.path(empty_cache, "kegg_compound_pathways.rds"))

    expect_null(
        run_compound_gsea(make_compound_de(), empty_cache, 5, 500, seed = 1L,
                          out_dir = NULL)
    )
})
