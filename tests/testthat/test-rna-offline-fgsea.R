# tests/testthat/test-rna-offline-fgsea.R
#
# The offline RNA enrichment path (modes.rna.enrichment.annotation_dir) runs
# GSEA through the shared fgsea call by default, with the clusterProfiler
# engine kept as an option for reproducing earlier runs. These tests cover the
# TERM2GENE conversion, the settings resolver and validator, the engine
# routing, the stored and written table schema, and the readers of the new
# tables. Synthetic data only.

# Replace a function in globalenv (where helper.R sources R/) for one test.
local_global_stub <- function(name, value, env = parent.frame()) {
    orig <- get(name, envir = globalenv())
    assign(name, value, envir = globalenv())
    withr::defer(assign(name, orig, envir = globalenv()), envir = env)
}

# 60 ranked genes; "pTOP" holds the 12 highest-ranked genes, "pMIX" 12 genes
# spread evenly over the ranking (not enriched), "pTINY" only 3 genes.
make_fgsea_tables <- function() {
    genes <- paste0("g", 1:60)
    mix   <- genes[seq(3, 60, by = 5)]
    t2g <- data.frame(
        term = c(rep("pTOP", 12), rep("pMIX", 12), rep("pTINY", 3)),
        gene = c(genes[1:12], mix, genes[58:60]),
        stringsAsFactors = FALSE)
    t2n <- data.frame(term = c("pTOP", "pMIX", "pTINY"),
                      name = c("Top pathway", "Mixed pathway", "Tiny pathway"),
                      stringsAsFactors = FALSE)
    list(KEGG = list(TERM2GENE = t2g, TERM2NAME = t2n))
}

make_fgsea_ranked <- function() {
    genes <- paste0("g", 1:60)
    v <- setNames(seq(6, -6, length.out = 60), genes)
    list(pval_with_direction = list(c1 = v),
         fc                  = list(c1 = v))
}

# ---------------------------------------------------------------------------
# TERM2GENE -> gene-set list
# ---------------------------------------------------------------------------

test_that("term2gene_to_gene_sets builds one unique set per term, in file order", {
    t2g <- data.frame(term = c("b", "a", "b", "b", "a", NA, ""),
                      gene = c("g1", "g2", "g3", "g1", "g4", "g5", "g6"),
                      stringsAsFactors = FALSE)
    t2n <- data.frame(term = c("a", "b", "zzz"), name = c("Alpha", "Beta", "Unused"),
                      stringsAsFactors = FALSE)

    sets <- term2gene_to_gene_sets(t2g, t2n)

    expect_identical(names(sets), c("b", "a"))        # first appearance, not sorted
    expect_identical(sets$b, c("g1", "g3"))            # duplicates removed
    expect_identical(sets$a, c("g2", "g4"))
    expect_false(any(c("g5", "g6") %in% unlist(sets))) # NA / empty term dropped
    expect_identical(attr(sets, "descriptions"), c(b = "Beta", a = "Alpha"))
})

test_that("term2gene_to_gene_sets works without TERM2NAME", {
    t2g <- data.frame(term = c("t1", "t1"), gene = c("x", "y"), stringsAsFactors = FALSE)
    sets <- term2gene_to_gene_sets(t2g)
    expect_identical(sets, list(t1 = c("x", "y")))
})

# ---------------------------------------------------------------------------
# Settings resolver and validator
# ---------------------------------------------------------------------------

test_that("resolve_rna_gsea_settings defaults to fgsea, the earlier offline limits and params.seed", {
    s <- resolve_rna_gsea_settings(list(modes = list(rna = list())))
    expect_identical(s$engine, "fgsea")
    expect_equal(s$min_size, 4)
    expect_true("max_size" %in% names(s))
    expect_null(s$max_size)          # = all annotated genes, per database
    expect_equal(s$seed, 1L)

    s2 <- resolve_rna_gsea_settings(list(params = list(seed = 42L),
                                         modes = list(rna = list())))
    expect_equal(s2$seed, 42L)
})

test_that("resolve_rna_gsea_settings reads enrichment.gsea sizes and ignores the online pathway block", {
    cfg <- list(modes = list(rna = list(
        pathway    = list(min_size = 15, max_size = 300),
        enrichment = list(gsea_engine = "clusterProfiler"))))
    s <- resolve_rna_gsea_settings(cfg)
    expect_identical(s$engine, "clusterProfiler")
    expect_equal(s$min_size, 4)
    expect_null(s$max_size)

    cfg$modes$rna$enrichment$gsea <- list(min_size = 5, max_size = 200)
    s <- resolve_rna_gsea_settings(cfg)
    expect_equal(c(s$min_size, s$max_size), c(5, 200))
})

test_that("resolve_rna_gsea_settings resolves one GSEA significance cutoff", {
    cfg <- function(enr) list(modes = list(rna = list(enrichment = enr)))
    expect_equal(resolve_rna_gsea_settings(cfg(list()))$sig_cutoff, 0.05)
    expect_equal(resolve_rna_gsea_settings(cfg(list(pvalue_cutoff = 0.1)))$sig_cutoff, 0.1)
    expect_equal(resolve_rna_gsea_settings(
        cfg(list(pvalue_cutoff = 0.1, gsea_pvalue_cutoff = 0.2)))$sig_cutoff, 0.2)
    # The DE cutoff is not part of the chain (the module never used it).
    expect_equal(resolve_rna_gsea_settings(list(modes = list(rna = list(
        de = list(p_cutoff = 0.01)))))$sig_cutoff, 0.05)

    # ORA runs at enrichment.pvalue_cutoff only; gsea_pvalue_cutoff does not reach it.
    s <- resolve_rna_gsea_settings(cfg(list(pvalue_cutoff = 0.1, gsea_pvalue_cutoff = 0.2)))
    expect_equal(s$ora_cutoff, 0.1)
    expect_equal(resolve_rna_gsea_settings(cfg(list()))$ora_cutoff, 0.05)
})

test_that("the module's filter and the Shiny payload read the same GSEA cutoff, with <=", {
    cfg <- list(modes = list(rna = list(enrichment = list(pvalue_cutoff = 0.1))))
    expect_equal(.build_enrichment_config(list(), cfg)$gsea_pvalue_cutoff,
                 resolve_rna_gsea_settings(cfg)$sig_cutoff)

    # A set exactly at the cutoff counts as significant.
    tab <- data.frame(ID = c("a", "b"), padj = c(0.1, 0.10001))
    expect_identical(.gsea_significant_rows(tab, 0.1)$ID, "a")
    # The report's use of the same cutoff is tested by rendering its chunks
    # (see "the report filters GSEA tables at a non-default cutoff, inclusive").
})

test_that("the pipeline summary counts pathways with the module's cutoffs, not the DE cutoff", {
    gsea <- data.frame(pathway = paste0("p", 1:4), NES = c(2, -2, 1, 1),
                       padj = c(0.001, 0.03, 0.05, 0.2))
    ora  <- data.frame(pathway = paste0("o", 1:3), padj = c(0.001, 0.03, 0.05))
    pw   <- list(pathway_results = list(c1 = list(KEGG_gsea_fc = gsea),
                                        cluster_ora = list(KEGG_contrasts_c1_ora = ora)))
    stats <- function(rna) collect_pipeline_stats(
        list(modes = list(rna = rna)), pre = list(), de_res = list(), pathway_res = pw)

    # DE cutoff differs from the (unset) enrichment cutoffs: the module runs at
    # 0.05, so the summary must too. Before, both fell back to de.p_cutoff (0.01).
    s <- stats(list(de = list(p_cutoff = 0.01)))
    expect_equal(s$pathway$total, 3 + 2)   # GSEA padj <= 0.05: 3; ORA padj < 0.05: 2
    expect_equal(s$methods$pathway_p_cutoff, 0.05)

    # Separate GSEA and ORA cutoffs are each applied to their own tables.
    s <- stats(list(de = list(p_cutoff = 0.01),
                    enrichment = list(pvalue_cutoff = 0.04, gsea_pvalue_cutoff = 0.03)))
    expect_equal(s$pathway$total, 2 + 2)   # GSEA <= 0.03: 2; ORA < 0.04: 2
    expect_equal(c(s$pathway$up, s$pathway$down), c(1, 1))
    expect_equal(s$methods$pathway_p_cutoff, 0.03)
})

test_that("resolve_rna_gsea_settings accepts any case and rejects unknown engines", {
    cfg <- list(modes = list(rna = list(enrichment = list(gsea_engine = "FGSEA"))))
    expect_identical(resolve_rna_gsea_settings(cfg)$engine, "fgsea")

    cfg$modes$rna$enrichment$gsea_engine <- "gsea"
    expect_error(resolve_rna_gsea_settings(cfg), "gsea_engine")
})

test_that("validate_rna_config checks gsea_engine and the gsea size limits", {
    base <- list(id_columns = list(gene_id = "gene_id"),
                 files = list(counts = "counts.csv", metadata = "metadata.csv"))

    ok <- base
    ok$enrichment <- list(gsea_engine = "fgsea",
                          gsea = list(min_size = 10, max_size = 500))
    expect_true(validate_rna_config(ok))
    ok$enrichment$gsea_engine <- "clusterProfiler"
    expect_true(validate_rna_config(ok))

    bad <- base
    bad$enrichment <- list(gsea_engine = "camera")
    expect_error(validate_rna_config(bad), "gsea_engine")

    bad$enrichment <- list(gsea = list(min_size = 500, max_size = 10))
    expect_error(validate_rna_config(bad), "min_size")

    bad$enrichment <- list(gsea = list(min_size = 0))
    expect_error(validate_rna_config(bad), "min_size")
})

test_that("gsea_engine without a gsea block is not read as the gsea block", {
    # gsea_engine is the only gsea* key, so enrichment$gsea would partially
    # match it and hand back the string "fgsea".
    cfg <- list(id_columns = list(gene_id = "gene_id"),
                files = list(counts = "counts.csv", metadata = "metadata.csv"),
                enrichment = list(gsea_engine = "fgsea"))
    expect_true(validate_rna_config(cfg))

    s <- resolve_rna_gsea_settings(list(modes = list(rna = list(
        enrichment = list(gsea_engine = "fgsea")))))
    expect_identical(s$engine, "fgsea")
})

test_that("a non-BH gsea_padj_method warns under the fgsea engine only", {
    cfg <- list(id_columns = list(gene_id = "gene_id"),
                files = list(counts = "counts.csv", metadata = "metadata.csv"),
                enrichment = list(gsea_padj_method = "bonferroni"))
    expect_warning(validate_rna_config(cfg), "clusterProfiler engine only")
    cfg$enrichment$gsea_engine <- "fgsea"
    expect_warning(validate_rna_config(cfg), "clusterProfiler engine only")
    cfg$enrichment$gsea_engine <- "clusterProfiler"
    expect_no_warning(validate_rna_config(cfg))
    cfg$enrichment <- list(gsea_engine = "fgsea", gsea_padj_method = "BH")
    expect_no_warning(validate_rna_config(cfg))
})

test_that("run_fgsea_core tests an unsigned ranking with scoreType 'pos'", {
    skip_if_not_installed("fgsea")
    genes <- paste0("g", 1:200)
    pos <- setNames(seq(10, 0.05, length.out = 200), genes)   # all > 0, like -log10(p)
    sets <- list(top = genes[1:20], bottom = genes[181:200], mid = genes[seq(5, 200, by = 10)])

    expect_no_warning(res <- run_fgsea_core(pos, sets, min_size = 5, n_perm_simple = 1000))
    # Under "std" the bottom set is tested as a negative tail and gets pval NA.
    expect_false(anyNA(res$pval))
    expect_equal(nrow(res), 3L)
    expect_lt(res$pval[res$pathway == "top"], 0.05)
    expect_equal(res$pval[res$pathway == "bottom"], 1)

    # A signed ranking keeps "std" (a set at the bottom is negatively enriched).
    signed <- setNames(seq(5, -5, length.out = 200), genes)
    res_s <- suppressWarnings(run_fgsea_core(signed, sets, min_size = 5, n_perm_simple = 1000))
    expect_lt(res_s$NES[res_s$pathway == "bottom"], 0)
    expect_error(run_fgsea_core(pos, sets, score_type = "bad"))
})

test_that("the RNA config template declares the fgsea engine and leaves the size limits at their defaults", {
    tmpl <- yaml::read_yaml(file.path(root_dir, "config", "templates", "rna_config.yaml"))
    enr <- tmpl$modes$rna$enrichment
    expect_identical(enr$gsea_engine, "fgsea")
    # Set limits would override the earlier offline defaults (4 / all genes).
    expect_null(enr$gsea$min_size)
    expect_null(enr$gsea$max_size)
    s <- resolve_rna_gsea_settings(tmpl)
    expect_equal(s$min_size, 4)
    expect_null(s$max_size)
})

# ---------------------------------------------------------------------------
# Engine routing (no fgsea / clusterProfiler needed: the engines are stubbed)
# ---------------------------------------------------------------------------

test_that("run_gsea_all routes the fgsea engine through run_fgsea_core with gene-set lists", {
    skip_if_not_installed("BiocParallel")
    seen <- new.env()
    local_global_stub("run_fgsea_core", function(ranks, gene_sets, min_size, max_size,
                                                 seed, n_perm_simple = 10000L,
                                                 bpparam = NULL) {
        seen$gene_sets <- gene_sets
        seen$args <- list(min_size = min_size, max_size = max_size, seed = seed)
        data.frame()
    })

    out <- suppressMessages(run_gsea_all(
        make_fgsea_ranked(), make_fgsea_tables(), output_dir = NULL, workers = 1,
        seed = 7L, min_size = 3, max_size = 40))

    expect_true(is.list(seen$gene_sets))
    expect_setequal(names(seen$gene_sets), c("pTOP", "pMIX", "pTINY"))
    expect_identical(seen$args, list(min_size = 3, max_size = 40, seed = 7L))
    expect_true(all(out$manifest$status == "empty"))
})

test_that("run_gsea_all defaults to min 4 and max = unique genes of each TERM2GENE", {
    skip_if_not_installed("BiocParallel")
    seen <- new.env()
    seen$max <- list()
    local_global_stub("run_fgsea_core", function(ranks, gene_sets, min_size, max_size,
                                                 seed, n_perm_simple = 10000L,
                                                 bpparam = NULL) {
        seen$min <- min_size
        seen$max[[length(seen$max) + 1]] <- max_size
        data.frame()
    })
    tables <- make_fgsea_tables()
    # A second database with a different gene count, so the limit is per database.
    tables$GO_BP <- list(
        TERM2GENE = data.frame(term = c("go1", "go1", "go2"),
                               gene = c("g1", "g2", "g2"), stringsAsFactors = FALSE),
        TERM2NAME = data.frame(term = c("go1", "go2"), name = c("A", "B"),
                               stringsAsFactors = FALSE))
    ranked <- list(pval_with_direction = make_fgsea_ranked()$pval_with_direction)

    suppressMessages(run_gsea_all(ranked, tables, output_dir = NULL, workers = 1))

    expect_equal(seen$min, 4)
    n_kegg <- length(unique(tables$KEGG$TERM2GENE[, 2]))
    expect_setequal(unlist(seen$max), c(n_kegg, 2L))
    # An explicit max_size is passed through unchanged: see the routing test above.
})

test_that("the online path keeps pathway.min_size/max_size and ignores enrichment.gsea", {
    gmt <- withr::local_tempfile(fileext = ".gmt")
    writeLines("SET_A\tdesc A\tg1\tg2\tg3\tg4\tg5", gmt)
    seen <- new.env()
    local_global_stub("run_pathway_analysis", function(..., min_size, max_size) {
        seen$sizes <- c(min_size, max_size)
        list()
    })
    de_res <- list(tables = list(c1 = data.frame(
        FeatureID = paste0("g", 1:10), log2FoldChange = seq(2, -2, length.out = 10),
        pvalue = seq(0.001, 0.5, length.out = 10), padj = seq(0.01, 0.6, length.out = 10),
        stringsAsFactors = FALSE)))
    run <- function(rna) {
        cfg <- list(modes = list(rna = rna))
        out <- withr::local_tempdir()
        suppressMessages(suppressWarnings(
            mod_rnaseq_pathway(de_res, list(expr_filt = NULL, meta = NULL), cfg, out)))
        seen$sizes
    }
    base <- list(annotation = list(skip_annotation = TRUE, organism = "Homo sapiens"),
                 pathway = list(enabled = TRUE, databases = character(0), gmt_file = gmt))

    expect_equal(run(base), c(10, 500))

    base$pathway$min_size <- 3
    base$pathway$max_size <- 50
    base$enrichment <- list(gsea = list(min_size = 99, max_size = 999))
    expect_equal(run(base), c(3, 50))
})

test_that("engine = 'clusterProfiler' routes to clusterProfiler::GSEA on TERM2GENE", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("BiocParallel")
    seen <- new.env()
    testthat::local_mocked_bindings(
        GSEA = function(geneList, TERM2GENE, TERM2NAME, ...) {
            seen$TERM2GENE <- TERM2GENE
            seen$args <- list(...)
            data.frame()
        },
        .package = "clusterProfiler")
    local_global_stub("run_fgsea_core", function(...) stop("fgsea engine must not run"))

    out <- suppressMessages(run_gsea_all(
        make_fgsea_ranked(), make_fgsea_tables(), output_dir = NULL, workers = 1,
        engine = "clusterProfiler"))

    expect_true(is.data.frame(seen$TERM2GENE))
    expect_equal(seen$args$minGSSize, 4)
    expect_false(any(out$manifest$status == "failed"))
})

test_that("run_gsea_all rejects an unknown engine", {
    expect_error(run_gsea_all(make_fgsea_ranked(), make_fgsea_tables(), engine = "camera"))
})

# ---------------------------------------------------------------------------
# fgsea engine end to end (real fgsea)
# ---------------------------------------------------------------------------

test_that("fgsea engine keeps the clusterProfiler contract in memory and every tested set on disk", {
    skip_if_not_installed("fgsea")
    skip_if_not_installed("BiocParallel")
    out_dir <- withr::local_tempdir()

    out <- suppressWarnings(suppressMessages(run_gsea_all(
        make_fgsea_ranked(), make_fgsea_tables(), output_dir = out_dir, workers = 1,
        seed = 1L, min_size = 5, max_size = 500, dotplot = FALSE)))

    # Every ranking alike, under <db>_gsea_<ranking>; no primary ranking.
    c1 <- out$results$c1
    expect_true(all(grepl("^KEGG_gsea_", names(c1))))
    expect_false(any(names(c1) == "KEGG_fgsea"))
    expect_setequal(out$index$storage_key, paste0("KEGG_gsea_", out$index$group))

    for (tab in c1) {
        expect_false("method" %in% names(tab))            # no column main did not have
        expect_false("pathway_name" %in% names(tab))
        expect_true(all(tab$padj <= 0.05))                 # significant sets only
        expect_identical(tab$pathway, tab$Description)    # pathway = the name
        expect_identical(tab$p.adjust, tab$padj)
        expect_identical(tab$core_enrichment, gsub(",", "/", tab$leadingEdge, fixed = TRUE))
    }
    pwd <- c1$KEGG_gsea_pval_with_direction
    expect_identical(pwd$pathway[pwd$ID == "pTOP"], "Top pathway")

    # On disk, the full table: the non-enriched set too; the 3-gene set is
    # below min_size.
    disk <- read.csv(file.path(out_dir, "KEGG", "ranking_by_pval_with_direction", "c1",
                               "results.csv"), stringsAsFactors = FALSE)
    expect_setequal(disk$ID, c("pTOP", "pMIX"))
    expect_true(any(disk$padj > 0.05))
    expect_identical(disk$pathway[disk$ID == "pMIX"], "Mixed pathway")

    # The manifest counts the evaluated units whether or not they were stored.
    expect_true(all(out$manifest$evaluated))
    expect_equal(nrow(out$manifest), 2L)
})

test_that("the display gseaResult holds only significant rows", {
    skip_if_not_installed("DOSE")
    df <- data.frame(pathway = c("p1", "p2"), pval = c(1e-4, 0.5), padj = c(1e-3, 0.6),
                     log2err = NA_real_, ES = c(0.8, 0.1), NES = c(2.1, 0.3),
                     size = c(10L, 10L), leadingEdge = c("g1,g2", "g3"),
                     stringsAsFactors = FALSE)
    sets <- structure(list(p1 = paste0("g", 1:10), p2 = paste0("g", 11:20)),
                      descriptions = c(p1 = "One", p2 = "Two"))
    tab <- .format_offline_fgsea(df, sets, "KEGG", "c1", "pval_with_direction")
    sig <- .gsea_significant_rows(tab, 0.05)
    ranked <- setNames(as.numeric(20:1), paste0("g", 1:20))

    obj <- .fgsea_as_gsea_result(sig, ranked, sets, 0.05)
    skip_if(is.null(obj), "gseaResult class not available")
    expect_true(methods::is(obj, "gseaResult"))
    expect_identical(obj@result$ID, "p1")
    expect_identical(obj@result$core_enrichment, "g1/g2")
    expect_identical(obj@geneList, ranked)
})

test_that("an offline fgsea table has the clusterProfiler column contract", {
    # Expected values are written out here, not derived from the formatter, so
    # this is the check the formatter-built fixtures below rely on.
    raw <- data.frame(pathway = c("t2", "t1"), pval = c(0.3, 0.001), padj = c(0.4, 0.004),
                      log2err = NA_real_, ES = c(-0.2, 0.7), NES = c(-0.8, 2.2),
                      size = c(8L, 12L), leadingEdge = c("g9", "g1,g2"),
                      stringsAsFactors = FALSE)
    sets <- structure(list(t1 = paste0("g", 1:12), t2 = paste0("g", 13:20)),
                      descriptions = c(t1 = "Term one", t2 = ""))
    tab <- .format_offline_fgsea(raw, sets, "KEGG", "c1", "fc")

    expect_setequal(names(tab), c(
        "pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge",
        "contrast", "database", "ranking_method", "ID", "Description", "setSize",
        "enrichmentScore", "pvalue", "p.adjust", "core_enrichment"))
    expect_identical(tab$ID, c("t1", "t2"))                     # ordered by p-value
    expect_identical(tab$pathway, c("Term one", "t2"))          # name; the ID when unnamed
    expect_identical(tab$Description, tab$pathway)
    expect_identical(tab$core_enrichment, c("g1/g2", "g9"))
    expect_identical(tab$setSize, c(12L, 8L))
    expect_identical(tab$p.adjust, c(0.004, 0.4))
    expect_identical(unique(tab$ranking_method), "fc")
})

# ---------------------------------------------------------------------------
# .run_local_enrichment end to end
# ---------------------------------------------------------------------------

make_flat_fixture <- function(engine = NULL) {
    dir <- withr::local_tempdir(.local_envir = parent.frame())
    ann <- file.path(dir, "annotation"); dir.create(ann)
    lt <- make_fgsea_tables()$KEGG
    write.table(lt$TERM2GENE, file.path(ann, "KEGG_pathway2gene.tab"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(lt$TERM2NAME, file.path(ann, "KEGG_pathway2name.tab"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    genes <- paste0("g", 1:60)
    c1 <- data.frame(FeatureID = genes,
                     log2FoldChange = seq(3, -3, length.out = 60),
                     pvalue = c(10^-(seq(8, 1, length.out = 30)),
                                rev(10^-(seq(8, 1, length.out = 30)))),
                     stringsAsFactors = FALSE)
    c1$padj <- p.adjust(c1$pvalue, "BH")

    enr <- list(annotation_dir = ann, databases = "KEGG", workers = 1,
                pvalue_cutoff = 0.05,
                plots = list(dotplot = FALSE, ridgeplot = FALSE,
                             ridgeplot_all_genes = FALSE, shared_genes = FALSE,
                             pathway_heatmaps = FALSE),
                gsea = list(per_pathway_artifacts = FALSE, min_size = 5))
    if (!is.null(engine)) enr$gsea_engine <- engine
    config <- list(project = list(dir = dir), params = list(seed = 1L),
                   modes = list(rna = list(de = list(p_cutoff = 0.05),
                                           enrichment = enr)))
    list(de_tables = list(c1 = c1), genes = genes, config = config,
         out_dir = file.path(dir, "out"))
}

run_flat_fixture <- function(fx) {
    suppressWarnings(suppressMessages(.run_local_enrichment(
        de_tables = fx$de_tables, feature_ids = fx$genes,
        enr_cfg = fx$config$modes$rna$enrichment, config = fx$config,
        out_dir = fx$out_dir, clustering_res = NULL, pre = NULL)))
}

test_that("offline fgsea run writes no flat fgsea file and no primary-ranking key", {
    skip_if_not_installed("fgsea")
    skip_if_not_installed("BiocParallel")
    fx <- make_flat_fixture()
    res <- run_flat_fixture(fx)

    expect_length(list.files(file.path(fx$out_dir, "Enrichment"),
                             pattern = "^pathway_.*_fgsea\\.csv$"), 0L)
    keys <- unlist(lapply(res$pathway_results, names))
    expect_false(any(grepl("_fgsea$", keys)))
    # The full per-ranking tables are on disk.
    expect_true(file.exists(file.path(fx$out_dir, "Enrichment", "GSEA", "KEGG",
                                      "ranking_by_pval_with_direction", "c1",
                                      "results.csv")))
})

test_that(".run_local_enrichment runs with gsea_engine set and no gsea block", {
    skip_if_not_installed("BiocParallel")
    local_global_stub("run_fgsea_core", function(...) data.frame())
    fx <- make_flat_fixture(engine = "fgsea")
    # Leave gsea_engine as the only gsea* key: enr_cfg$gsea would partially
    # match it and fail on $per_pathway_artifacts.
    fx$config$modes$rna$enrichment$gsea <- NULL
    res <- run_flat_fixture(fx)
    gsea_rows <- res$enrichment_manifest[res$enrichment_manifest$analysis == "GSEA", ]
    expect_gt(nrow(gsea_rows), 0)
    expect_true(all(gsea_rows$evaluated))
})

test_that("offline enrichment records its manifest even when no unit has a table", {
    skip_if_not_installed("BiocParallel")
    # Every unit evaluates with no gene set in range, so no result table is written.
    local_global_stub("run_fgsea_core", function(...) data.frame())
    fx <- make_flat_fixture()
    run_flat_fixture(fx)

    f <- file.path(fx$out_dir, "Enrichment", "enrichment_manifest.csv")
    expect_true(file.exists(f))
    m <- read.csv(f, stringsAsFactors = FALSE)
    expect_gt(nrow(m), 0)
    expect_true(all(m$evaluated))
    expect_true(all(m$status == "empty"))
})

test_that("removing a stale manifest leaves every other output in place", {
    out <- withr::local_tempdir()
    enr <- file.path(out, "Enrichment")
    dir.create(file.path(enr, "GSEA", "KEGG"), recursive = TRUE)
    keep <- c(file.path(enr, "pathway_c1_KEGG_fgsea.csv"),
              file.path(enr, "gene_annotation.csv"),
              file.path(enr, "GSEA", "KEGG", "results.csv"),
              file.path(out, "enrichment_manifest.csv"))   # same name, outside Enrichment/
    for (f in keep) writeLines("x", f)
    writeLines("analysis,evaluated", file.path(enr, "enrichment_manifest.csv"))

    cfg <- list(modes = list(rna = list(enrichment = list(enabled = FALSE))))
    suppressMessages(mod_rnaseq_pathway(list(tables = list()), NULL, cfg, out))
    expect_false(file.exists(file.path(enr, "enrichment_manifest.csv")))
    expect_true(all(file.exists(keep)))
})

# ---------------------------------------------------------------------------
# RNA report: enrichment status, Methods and GSEA tables (rendered chunks)
# ---------------------------------------------------------------------------

# Body of one named chunk of the RNA report.
report_chunk <- function(label) {
    rmd <- readLines(file.path(root_dir, "R", "domain", "rnaseq", "report_template.Rmd"),
                     warn = FALSE)
    start <- grep(paste0("^```\\{r ", label, "[,}]"), rmd)
    end   <- start + which(grepl("^```\\s*$", rmd[(start + 1):length(rmd)]))[1]
    rmd[(start + 1):(end - 1)]
}

# The report's environment after its pathway-setup chunk, for an Enrichment
# directory and an RNA config.
report_env <- function(enrich_dir, rna_cfg = list()) {
    env <- new.env(parent = globalenv())
    env$enrich_dir      <- enrich_dir
    env$enrich_plot_dir <- file.path(enrich_dir, "plots")
    env$cfg             <- list(modes = list(rna = rna_cfg))
    env$rna_cfg         <- rna_cfg
    eval(parse(text = report_chunk("pathway-setup")), envir = env)
    env
}

# Printed output of report chunks, evaluated in order in `env`.
render_chunks <- function(env, labels) {
    paste(capture.output(for (l in labels) eval(parse(text = report_chunk(l)), envir = env)),
          collapse = "\n")
}

# One enrichment-manifest row, and the file of several.
manifest_row <- function(analysis, database, group, evaluated = TRUE, n_significant = 0L) {
    data.frame(analysis = analysis, database = database, group = group, item = "c1",
               evaluated = evaluated,
               status = if (!evaluated) "failed" else if (n_significant > 0) "significant" else "empty",
               n_significant = if (evaluated) n_significant else NA_integer_,
               stringsAsFactors = FALSE)
}
write_manifest <- function(d, rows) {
    write.csv(do.call(rbind, rows), file.path(d, "enrichment_manifest.csv"), row.names = FALSE)
}

# Enrichment status, as the report states it.
report_status <- function(enrich_dir) {
    env <- report_env(enrich_dir)
    chunk <- if (!env$has_enrichment) "pathway-skip" else if (env$offline_ran) "pathway-status"
    list(has_enrichment = env$has_enrichment, status = env$enr_status,
         text = if (is.null(chunk)) "" else render_chunks(env, chunk))
}

test_that("the RNA report tells not run, completed without significance, and failed apart", {
    d <- withr::local_tempdir()

    s <- report_status(d)                                           # nothing recorded
    expect_identical(s$status, "not_run")
    expect_false(s$has_enrichment)
    expect_match(s$text, "was not run", fixed = TRUE)

    write_manifest(d, list(manifest_row("GSEA", "KEGG", "fc"),
                           manifest_row("GSEA", "KEGG", "pval_with_direction")))
    s <- report_status(d)                                           # ran, nothing significant
    expect_identical(s$status, "completed")
    expect_true(s$has_enrichment)
    expect_match(s$text, "none reached significance", fixed = TRUE)

    write_manifest(d, list(manifest_row("GSEA", "KEGG", "fc", n_significant = 3L),
                           manifest_row("GSEA", "KEGG", "pval_with_direction",
                                        evaluated = FALSE)))
    s <- report_status(d)                                           # partly failed
    expect_identical(s$status, "completed")
    expect_match(s$text, "1 of 1 evaluated", fixed = TRUE)
    expect_match(s$text, "1 analysis unit(s) failed", fixed = TRUE)

    write_manifest(d, list(manifest_row("GSEA", "KEGG", "fc", evaluated = FALSE),
                           manifest_row("GSEA", "KEGG", "pval_with_direction",
                                        evaluated = FALSE)))
    s <- report_status(d)                                           # every unit failed
    expect_identical(s$status, "failed")
    expect_false(s$has_enrichment)
    expect_match(s$text, "no analysis unit completed", fixed = TRUE)
    expect_false(grepl("was not run", s$text, fixed = TRUE))
})

test_that("the RNA report counts the online flat files as a completed run", {
    d <- withr::local_tempdir()
    writeLines("pathway,padj", file.path(d, "pathway_c1_KEGG_fgsea.csv"))
    s <- report_status(d)
    expect_identical(s$status, "completed")
    expect_true(s$has_enrichment)
})

# The report's pathway Methods block (one if/else expression, grown until it parses).
report_methods_text <- function(enrich_dir, rna_cfg = list()) {
    env <- report_env(enrich_dir, rna_cfg)
    rmd <- readLines(file.path(root_dir, "R", "domain", "rnaseq", "report_template.Rmd"),
                     warn = FALSE)
    start <- grep("^pathway_methods_text <- ", rmd)
    for (end in start:length(rmd)) {
        expr <- tryCatch(parse(text = rmd[start:end]), error = function(e) NULL)
        if (!is.null(expr)) break
    }
    eval(expr, envir = env)
    env$pathway_methods_text
}

test_that("Methods state the GSEA and ORA cutoffs separately, for what ran", {
    d <- withr::local_tempdir()
    rna <- list(enrichment = list(pvalue_cutoff = 0.04, gsea_pvalue_cutoff = 0.03))

    write_manifest(d, list(manifest_row("GSEA", "KEGG", "fc"),
                           manifest_row("GSEA", "GO_BP", "pval_with_direction"),
                           manifest_row("ORA", "KEGG", "contrasts"),
                           manifest_row("GSEA", "GO_MF", "fc", evaluated = FALSE)))
    txt <- report_methods_text(d, rna)
    expect_match(txt, "fgsea on GO_BP, KEGG gene sets", fixed = TRUE)   # failed GO_MF left out
    expect_match(txt, "log2 fold change; signed -log10 p-value", fixed = TRUE)
    expect_match(txt, "not independent evidence", fixed = TRUE)
    expect_match(txt, "Benjamini-Hochberg adjusted p-value ≤ 0.03", fixed = TRUE)
    expect_match(txt, "enricher on KEGG gene sets", fixed = TRUE)
    expect_match(txt, "adjusted p-value cutoff of 0.04", fixed = TRUE)
    expect_false(grepl("0.05", txt, fixed = TRUE))

    # GSEA only: no ORA sentence.
    write_manifest(d, list(manifest_row("GSEA", "KEGG", "fc")))
    txt <- report_methods_text(d, rna)
    expect_false(grepl("Over-representation", txt, fixed = TRUE))

    # clusterProfiler engine: no Benjamini-Hochberg claim.
    rna$enrichment$gsea_engine <- "clusterProfiler"
    txt <- report_methods_text(d, rna)
    expect_match(txt, "clusterProfiler::GSEA", fixed = TRUE)
    expect_false(grepl("Benjamini", txt, fixed = TRUE))
})

test_that("Methods say nothing about pathways when enrichment did not run", {
    d <- withr::local_tempdir()
    expect_identical(report_methods_text(d), "")
    # Online flat fgsea files: the online sentence, with the GSEA cutoff.
    writeLines("pathway,padj", file.path(d, "pathway_c1_KEGG_fgsea.csv"))
    expect_match(report_methods_text(d), "adjusted p-value ≤ 0.05 are reported", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# Readers of the stored tables
# ---------------------------------------------------------------------------

# In-memory tables as run_gsea_all() stores them with the fgsea engine: the
# significant rows of .format_offline_fgsea(), one per ranking.
make_offline_like_results <- function() {
    sets <- structure(list(GO1 = paste0("g", 1:10), GO2 = paste0("g", 11:20)),
                      descriptions = c(GO1 = "Term one", GO2 = "Term two"))
    base <- data.frame(pathway = c("GO1", "GO2"), pval = c(1e-4, 0.4),
                       padj = c(1e-3, 0.5), log2err = NA_real_, ES = c(0.7, 0.2),
                       NES = c(2, 0.5), size = c(10L, 10L),
                       leadingEdge = c("g1,g2", "g11"), stringsAsFactors = FALSE)
    stored <- function(ranking) .gsea_significant_rows(
        .format_offline_fgsea(base, sets, "GO_BP", "c1", ranking), 0.05)
    list(pwd = stored("pval_with_direction"), fc = stored("fc"))
}

test_that("multiomics sees the main-branch shape: no method column, every ranking pooled", {
    # Pooling every ranking of the RNA layer (select_layer_method_rows() keeps
    # all rows without a method column) is kept for compatibility only, pending
    # the multiomics review; it is not an approved choice.
    tabs <- make_offline_like_results()
    rna_pathway_res <- list(pathway_results = list(
        c1 = list(GO_BP_gsea_pval_with_direction = tabs$pwd, GO_BP_gsea_fc = tabs$fc)))

    df <- extract_enrichment_df(rna_pathway_res)
    expect_equal(nrow(df), 2L)                      # significant rows only, both rankings
    expect_false("method" %in% names(df))
    sel <- select_layer_method_rows(df, .raw_p_values(df))
    expect_identical(sel$method, "unspecified")
    expect_true(all(sel$keep))
})

test_that("Shiny payload keeps pathway as the name and ID as the identifier", {
    tabs <- make_offline_like_results()
    pw <- list(
        pathway_results = list(c1 = list(GO_BP_gsea_pval_with_direction = tabs$pwd)),
        enrichment_index = data.frame(
            analysis = "GSEA", database = "GO_BP", group = "pval_with_direction",
            item = "c1", container = "c1", storage_key = "GO_BP_gsea_pval_with_direction",
            has_simplify = FALSE, simplify_key = NA_character_,
            stringsAsFactors = FALSE))
    gi <- paste0("g", 1:20)

    out <- .restructure_gsea(pw, gi, sig_cutoff = 0.05)
    tab <- out$gsea$GO_BP$pval_with_direction$c1
    expect_identical(tab$ID, "GO1")
    expect_identical(tab$pathway, "Term one")
    expect_false(any(c("core_enrichment", "leadingEdge", "method") %in% names(tab)))
    expect_identical(out$gsea_leading_edge$GO_BP$pval_with_direction$c1$GO1, 1:2)
})

test_that("pathway membership is keyed by term ID", {
    tabs <- make_offline_like_results()
    lt <- list(GO_BP = list(TERM2GENE = data.frame(
        term = c("GO1", "GO1", "GO2"), gene = c("g1", "g2", "g11"),
        stringsAsFactors = FALSE)))
    mem <- .build_gsea_pathway_membership(
        list(c1 = list(GO_BP_gsea_pval_with_direction = tabs$pwd)), lt, sig_cutoff = 0.05)
    expect_identical(names(mem$GO_BP), "GO1")
})

# ---------------------------------------------------------------------------
# Report: offline GSEA tables, one row per pathway and ranking
# ---------------------------------------------------------------------------

# One unit's results.csv, laid out as run_gsea_all() writes it.
write_gsea_unit <- function(enrich_dir, db, ranking, contrast, df) {
    d <- file.path(enrich_dir, "GSEA", db, paste0("ranking_by_", ranking), contrast)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    df$contrast <- contrast; df$database <- db; df$ranking_method <- ranking
    write.csv(df, file.path(d, "results.csv"), row.names = FALSE)
}

# An fgsea-engine unit for terms T1..Tn ("Term <id>"), built by the production
# formatter so the fixture cannot drift from it (its contract is pinned above).
write_fgsea_unit <- function(enrich_dir, db, ranking, contrast, ids, padj, nes) {
    raw <- data.frame(pathway = ids, pval = padj / 2, padj = padj, log2err = NA_real_,
                      ES = nes / 2, NES = nes, size = 10L, leadingEdge = "g1,g2,g3",
                      stringsAsFactors = FALSE)
    sets <- structure(setNames(rep(list(paste0("g", 1:10)), length(ids)), ids),
                      descriptions = setNames(paste("Term", ids), ids))
    write_gsea_unit(enrich_dir, db, ranking, contrast,
                    .format_offline_fgsea(raw, sets, db, contrast, ranking))
}

test_that("read_offline_gsea_tables keeps one row per pathway and ranking, significant only", {
    d <- withr::local_tempdir()
    write_fgsea_unit(d, "KEGG", "fc", "c1", c("T1", "T2"), c(0.01, 0.2), c(2, -1))
    write_fgsea_unit(d, "KEGG", "pval_with_direction", "c1", c("T1", "T3"),
                     c(0.02, 0.05), c(1.8, -1.5))
    write_fgsea_unit(d, "KEGG", "pval_wo_direction", "c1", "T1", 0.03, 1.2)

    tab <- read_offline_gsea_tables(file.path(d, "GSEA"), 0.05)
    # T1 in three rankings is three rows, each with its own p-values; T2 is
    # not significant; T3 sits exactly on the cutoff.
    expect_equal(nrow(tab), 4L)
    t1 <- tab[tab$ID == "T1", ]
    expect_equal(nrow(t1), 3L)
    expect_setequal(t1$padj, c(0.01, 0.02, 0.03))
    expect_identical(t1$Ranking, c("log2 fold change", "signed -log10 p-value",
                                   "unsigned -log10 p-value"))
    expect_false("T2" %in% tab$ID)
    expect_true("T3" %in% tab$ID)
    expect_identical(tab$Pathway[tab$ID == "T3"], "Term T3")
    # No direction for the unsigned ranking.
    expect_true(is.na(t1$Direction[t1$ranking_method == "pval_wo_direction"]))
    expect_identical(t1$Direction[t1$ranking_method == "fc"], "Up-regulated")
    expect_identical(tab$Direction[tab$ID == "T3"], "Down-regulated")
    expect_identical(tab$leading_edge[1], "g1,g2,g3")
})

test_that("read_offline_gsea_tables maps the clusterProfiler engine's columns", {
    d <- withr::local_tempdir()
    # Written out by hand: the clusterProfiler table as run_gsea_all() stores it.
    cp <- data.frame(ID = "map00010", Description = "Glycolysis", setSize = 12L,
                     enrichmentScore = -0.6, NES = -1.9, pvalue = 0.001, p.adjust = 0.01,
                     qvalue = 0.01, rank = 5L, leading_edge = "tags=40%",
                     core_enrichment = "g1/g2", padj = 0.01, pathway = "Glycolysis",
                     stringsAsFactors = FALSE)
    write_gsea_unit(d, "KEGG", "fc", "c1", cp)
    tab <- read_offline_gsea_tables(file.path(d, "GSEA"), 0.05)
    expect_equal(nrow(tab), 1L)
    expect_identical(tab$Pathway, "Glycolysis")
    expect_identical(tab$ID, "map00010")
    expect_equal(tab$pval, 0.001)
    expect_identical(tab$size, 12L)
    expect_identical(tab$Direction, "Down-regulated")
    expect_identical(tab$leading_edge, "g1,g2")
})

test_that("read_offline_gsea_tables returns the empty layout when nothing is there", {
    d <- withr::local_tempdir()
    expect_equal(nrow(read_offline_gsea_tables(file.path(d, "GSEA"), 0.05)), 0L)
    write_fgsea_unit(d, "KEGG", "fc", "c1", "T1", 0.5, 1)
    tab <- read_offline_gsea_tables(file.path(d, "GSEA"), 0.05)
    expect_equal(nrow(tab), 0L)
    expect_true(all(c("Pathway", "ID", "Ranking", "NES", "pval", "padj", "size") %in% names(tab)))
})

test_that("the report shows every ranking per contrast/database, and units with no result", {
    skip_if_not_installed("knitr")
    d <- withr::local_tempdir()
    write_fgsea_unit(d, "KEGG", "fc", "c1", c("T1", "T2"), c(0.01, 0.2), c(2, -1))
    write_fgsea_unit(d, "KEGG", "pval_with_direction", "c1", "T1", 0.02, 1.8)
    write_fgsea_unit(d, "GO_BP", "fc", "c1", "G1", 0.6, 0.5)    # nothing significant
    write_manifest(d, list(manifest_row("GSEA", "KEGG", "fc"),
                           manifest_row("GSEA", "KEGG", "pval_with_direction"),
                           manifest_row("GSEA", "GO_BP", "fc")))

    env <- report_env(d)
    tables <- render_chunks(env, c("fgsea-tables-prep", "fgsea-tables"))
    expect_match(tables, "not independent evidence", fixed = TRUE)
    expect_match(tables, "### c1 — KEGG", fixed = TRUE)
    expect_match(tables, "2 row(s): 1 distinct pathway(s) across 2 ranking(s)", fixed = TRUE)
    expect_match(tables, "log2 fold change", fixed = TRUE)
    expect_match(tables, "signed -log10 p-value", fixed = TRUE)
    expect_false(grepl("Term T2", tables, fixed = TRUE))       # not significant
    expect_match(tables, "### c1 — GO_BP", fixed = TRUE)
    expect_match(tables, "No significant pathways (padj ≤ 0.05) in any ranking", fixed = TRUE)
    # Leading edge: per ranking, the genes as reported.
    le <- render_chunks(env, "leading-edge")
    expect_match(le, "g1, g2, g3", fixed = TRUE)
    expect_match(le, "signed -log10 p-value", fixed = TRUE)
})

test_that("the report filters GSEA tables at a non-default cutoff, inclusive", {
    skip_if_not_installed("knitr")
    rna <- list(enrichment = list(gsea_pvalue_cutoff = 0.1))

    # Offline: per-ranking results.csv.
    d <- withr::local_tempdir()
    write_fgsea_unit(d, "KEGG", "fc", "c1", c("Tin", "Tat", "Tout"),
                     c(0.08, 0.1, 0.12), c(1, 1, 1))
    write_manifest(d, list(manifest_row("GSEA", "KEGG", "fc")))
    env <- report_env(d, rna)
    tables <- render_chunks(env, c("fgsea-tables-prep", "fgsea-tables"))
    le <- render_chunks(env, "leading-edge")
    for (out in list(tables, le)) {
        expect_match(out, "Term Tin", fixed = TRUE)
        expect_match(out, "Term Tat", fixed = TRUE)            # exactly on the cutoff
        expect_false(grepl("Term Tout", out, fixed = TRUE))
    }

    # Online: the flat pathway_<contrast>_<db>_fgsea.csv file.
    d2 <- withr::local_tempdir()
    write.csv(data.frame(pathway = c("Oin", "Oat", "Oout"), pathway_name = c("Oin", "Oat", "Oout"),
                         pval = c(0.04, 0.05, 0.06), padj = c(0.08, 0.1, 0.12), NES = 1,
                         size = 10L, leadingEdge = "g1,g2"),
              file.path(d2, "pathway_c1_KEGG_fgsea.csv"), row.names = FALSE)
    env2 <- report_env(d2, rna)
    tables2 <- render_chunks(env2, c("fgsea-tables-prep", "fgsea-tables"))
    le2 <- render_chunks(env2, "leading-edge")
    for (out in list(tables2, le2)) {
        expect_match(out, "Oin", fixed = TRUE)
        expect_match(out, "Oat", fixed = TRUE)
        expect_false(grepl("Oout", out, fixed = TRUE))
    }
})
