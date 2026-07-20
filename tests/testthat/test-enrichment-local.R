# tests/testthat/test-enrichment-local.R
#
# Phase 1 — local (offline) enrichment unit tests.
# Pure-R functions (loader, rankers, table processing, gene-list builder) are
# tested directly. clusterProfiler-dependent paths (GSEA/ORA) are skipped when
# the package is unavailable. Synthetic fixtures only — no real data.

fx_dir <- test_path("fixtures", "enrichment_local")

# ---------------------------------------------------------------------------
# load_local_pathway_tables: present DBs load; missing DBs skip + warn; overlap
# ---------------------------------------------------------------------------

test_that("load_local_pathway_tables loads present DBs and standardizes columns", {
    tabs <- load_local_pathway_tables(fx_dir, databases = c("KEGG", "GO_BP"))
    expect_true(all(c("KEGG", "GO_BP") %in% names(tabs)))
    expect_named(tabs$KEGG, c("TERM2GENE", "TERM2NAME"))
    expect_equal(colnames(tabs$KEGG$TERM2GENE), c("term", "gene"))
    expect_equal(colnames(tabs$KEGG$TERM2NAME), c("term", "name"))
    # header row was consumed, not treated as data
    expect_false(any(tabs$KEGG$TERM2GENE$gene == "gene"))
    expect_setequal(unique(tabs$KEGG$TERM2GENE$term), c("path1", "path2"))
})

test_that("missing database files are skipped with a warning, others still load", {
    expect_warning(
        tabs <- load_local_pathway_tables(fx_dir, databases = c("KEGG", "GO_BP", "GO_MF")),
        "GO_MF"
    )
    expect_true(all(c("KEGG", "GO_BP") %in% names(tabs)))
    expect_false("GO_MF" %in% names(tabs))
})

test_that("low feature overlap emits a warning (does not crash)", {
    expect_warning(
        load_local_pathway_tables(fx_dir, databases = "KEGG",
                                  feature_ids = paste0("zzz", 1:100)),
        "low overlap"
    )
})

test_that("non-existent annotation_dir errors clearly", {
    expect_error(load_local_pathway_tables(file.path(fx_dir, "nope")),
                 "Annotation directory not found")
})

# ---------------------------------------------------------------------------
# Ranked gene list builders
# ---------------------------------------------------------------------------

de_one <- data.frame(
    FeatureID      = c("g1", "g2", "g3"),
    log2FoldChange = c(1, -2, 0.5),
    pvalue         = c(1e-4, 1e-2, 1),
    stringsAsFactors = FALSE
)

test_that("rank_by_pval_wo_direction is unsigned -log10(p), sorted descending", {
    r <- rank_by_pval_wo_direction(de_one)
    expect_equal(names(r), c("g1", "g2", "g3"))
    expect_equal(unname(r["g1"]), 4)
    expect_equal(unname(r["g2"]), 2)
    expect_equal(unname(r["g3"]), 0)  # -log10(1) = 0
})

test_that("rank_by_pval_with_direction applies sign of log2FC", {
    r <- rank_by_pval_with_direction(de_one)
    expect_equal(unname(r["g1"]),  4)   # +sign, p=1e-4
    expect_equal(unname(r["g2"]), -2)   # -sign, p=1e-2
    expect_equal(unname(r["g3"]),  0)
    expect_true(which(names(r) == "g1") < which(names(r) == "g2"))
})

test_that("rank_by_fc recovers linear FC, applies legacy transform + signif(4)", {
    r <- rank_by_fc(de_one)
    # g1: 2^1=2 -> log2(2)=1 ; g2: -(2^2)=-4 -> -1/-4=0.25 -> log2=-2 ; g3: 2^0.5 -> log2=0.5
    expect_equal(unname(r["g1"]),  1)
    expect_equal(unname(r["g2"]), -2)
    expect_equal(unname(r["g3"]), 0.5)
})

test_that("rank_by_min_pval_any_contrast takes per-gene min p across contrasts", {
    de_tables <- list(
        c1 = data.frame(FeatureID = c("g1", "g2"), log2FoldChange = c(1, 1),
                        pvalue = c(1e-4, 1e-2), stringsAsFactors = FALSE),
        c2 = data.frame(FeatureID = c("g2", "g3"), log2FoldChange = c(1, 1),
                        pvalue = c(1e-6, 0.5), stringsAsFactors = FALSE)
    )
    r <- rank_by_min_pval_any_contrast(de_tables)
    expect_equal(unname(r["g2"]), 6)   # min(1e-2, 1e-6) -> -log10 = 6
    expect_equal(unname(r["g1"]), 4)
    expect_equal(names(r)[1], "g2")    # most significant first
})

test_that("build_ranked_gene_lists returns the four ranking variants + any_contrast", {
    de_tables <- list(c1 = de_one)
    rk <- build_ranked_gene_lists(de_tables)
    expect_setequal(names(rk), c("pval_wo_direction", "pval_with_direction", "fc"))
    expect_true("c1" %in% names(rk$fc))
    expect_true("any_contrast" %in% names(rk$pval_wo_direction))
})

# ---------------------------------------------------------------------------
# process_enrichment_table: GeneRatio/BgRatio expansion + fold enrichment
# ---------------------------------------------------------------------------

test_that("process_enrichment_table expands ratios and computes Fold_enrichment", {
    tbl <- data.frame(
        GeneRatio = "2/10",
        BgRatio   = "5/100",
        Count     = 2,
        geneID    = "g1/g2",
        stringsAsFactors = FALSE
    )
    out <- process_enrichment_table(tbl)
    expect_equal(out$in_cluster_in_term, 2)
    expect_equal(out$in_cluster, 10)
    expect_equal(out$in_term, 5)
    expect_equal(out$in_genome, 100)
    expect_equal(out$Fold_enrichment, signif((2 / 10) / (5 / 100), 2))  # = 4
})

# ---------------------------------------------------------------------------
# build_gene_lists: contrasts always; partition only when clustering provided
# ---------------------------------------------------------------------------

de_sig <- list(
    cA = data.frame(
        FeatureID      = c("g1", "g2", "g3"),
        log2FoldChange = c(2, -2, 0.1),       # g3 below lfc cutoff
        padj           = c(1e-3, 1e-3, 1e-3),
        stringsAsFactors = FALSE
    )
)

test_that("build_gene_lists builds up/down contrast lists and skips partition when NULL", {
    gl <- build_gene_lists(de_sig, clustering_res = NULL,
                           p_cutoff = 0.05, lfc_cutoff = log2(1.5))
    expect_true("contrasts" %in% names(gl))
    expect_equal(unname(gl$contrasts$cA["g1"]), "up")
    expect_equal(unname(gl$contrasts$cA["g2"]), "down")
    expect_false("g3" %in% names(gl$contrasts$cA))   # |lfc| below cutoff
    expect_false("partition" %in% names(gl))          # no clustering provided
})

test_that("build_gene_lists adds partition clusters when clustering_res is provided", {
    clustering_res <- list(objects = list(clusters = c(g1 = 1L, g2 = 2L, g3 = 1L)))
    gl <- build_gene_lists(de_sig, clustering_res = clustering_res)
    expect_true("partition" %in% names(gl))
    expect_equal(gl$partition$k, c(g1 = 1L, g2 = 2L, g3 = 1L))
})

test_that("build_gene_lists builds an all_DE collection (union of DE genes, single 'all' cluster)", {
    de_tables <- list(
        cA = data.frame(FeatureID = c("g1", "g2", "g3"), log2FoldChange = c(2, -2, 0.1),
                        padj = c(1e-3, 1e-3, 1e-3), stringsAsFactors = FALSE),
        cB = data.frame(FeatureID = c("g2", "g4"), log2FoldChange = c(2, 2),
                        padj = c(1e-3, 1e-3), stringsAsFactors = FALSE)
    )
    gl <- build_gene_lists(de_tables, clustering_res = NULL,
                           p_cutoff = 0.05, lfc_cutoff = log2(1.5))
    expect_true("all_DE" %in% names(gl))
    # union across contrasts; g3 excluded (|lfc| below cutoff)
    expect_setequal(names(gl$all_DE$any_contrast), c("g1", "g2", "g4"))
    expect_true(all(gl$all_DE$any_contrast == "all"))
})

# --- M2: partition sourced unambiguously (no mislabeling hierarchical) ---

test_that("partition prefers excel_order$partition_clusters over objects$clusters", {
    cr <- list(
        objects     = list(clusters = c(g1 = 9L, g2 = 9L)),          # would be wrong if used
        excel_order = list(partition_clusters = c(g1 = 1L, g2 = 2L)) # the real partition
    )
    gl <- build_gene_lists(list(), clustering_res = cr)
    expect_equal(gl$partition$k, c(g1 = 1L, g2 = 2L))
})

test_that("hierarchical-only run does NOT produce a partition collection", {
    # excel_order present (hierarchical ran) but partition_clusters NULL -> partition
    # did not run; objects$clusters then holds hierarchical cuts and must be ignored.
    cr <- list(
        objects     = list(clusters = c(g1 = 1L, g2 = 1L)),  # hierarchical cuts
        excel_order = list(ordered_ids = c("g1", "g2"), partition_clusters = NULL)
    )
    gl <- build_gene_lists(list(), clustering_res = cr)
    expect_false("partition" %in% names(gl))
})

test_that("partition-only run (no excel_order) uses objects$clusters", {
    cr <- list(objects = list(clusters = c(g1 = 1L, g2 = 2L)))  # excel_order absent
    gl <- build_gene_lists(list(), clustering_res = cr)
    expect_equal(gl$partition$k, c(g1 = 1L, g2 = 2L))
})

# --- M1: clustering IDs aligned to the DE / TERM2GENE space (Gene: prefix) ---

test_that("partition gene IDs are stripped of a leading 'Gene:' prefix", {
    cr <- list(objects = list(),
               excel_order = list(partition_clusters =
                                      setNames(c(1L, 2L), c("Gene:g1", "Gene:g2"))))
    gl <- build_gene_lists(list(), clustering_res = cr)
    expect_equal(names(gl$partition$k), c("g1", "g2"))
})

# --- Binary patterns: data.frame -> named vector, IDs stripped (M1 + M3) ---

test_that("binary_patterns read from the patterns data.frame, IDs stripped", {
    cr <- list(objects = list(patterns = data.frame(
        feature_id   = c("Gene:g1", "g2", "g3"),
        best_pattern = c("up", "down", NA),       # NA dropped
        best_corr    = c(0.9, 0.8, NA),
        stringsAsFactors = FALSE
    )))
    gl <- build_gene_lists(list(), clustering_res = cr)
    expect_true("binary_patterns" %in% names(gl))
    expect_equal(gl$binary_patterns$best, c(g1 = "up", g2 = "down"))
})

# ---------------------------------------------------------------------------
# clusterProfiler-dependent paths (skipped if package unavailable)
# ---------------------------------------------------------------------------

test_that("run_gsea_local returns NULL on empty input without error", {
    expect_null(run_gsea_local(numeric(0),
                               term2gene = data.frame(term = character(), gene = character()),
                               term2name = data.frame(term = character(), name = character())))
})

test_that("run_cluster_ora returns empty list for invalid clusters", {
    skip_if_not_installed("clusterProfiler")
    res <- run_cluster_ora(clusters = setNames(integer(0), character(0)),
                           TERM2GENE = data.frame(term = "p1", gene = "g1"),
                           TERM2NAME = data.frame(term = "p1", name = "P1"))
    expect_equal(length(res), 0)
})

# ---------------------------------------------------------------------------
# Phase 2 enrichment plots — fail-soft behavior (synthetic; no real data)
# ---------------------------------------------------------------------------

test_that("plot_gsea_ridgeplot is fail-soft on NULL/invalid input", {
    tmp <- withr::local_tempdir()
    out_dir <- file.path(tmp, "ridgeplot")
    # Second arg is now the output DIRECTORY (writes plot.png / data.csv inside).
    expect_warning(expect_null(plot_gsea_ridgeplot(NULL, out_dir)))
    expect_false(file.exists(file.path(out_dir, "plot.png")))
    expect_warning(expect_null(plot_gsea_ridgeplot(data.frame(a = 1), out_dir)))
    expect_false(file.exists(file.path(out_dir, "plot.png")))
})

test_that("plot_ora_shared_genes is fail-soft on empty / malformed input", {
    tmp <- withr::local_tempdir()
    # NULL and non-data.frame -> empty, no error, no files
    expect_equal(length(plot_ora_shared_genes(NULL, tmp)), 0)
    expect_equal(length(plot_ora_shared_genes(list(), tmp)), 0)
    # missing required columns -> warns + empty
    expect_warning(res <- plot_ora_shared_genes(data.frame(a = 1, b = 2), tmp))
    expect_equal(length(res), 0)
    # zero-row frame with the right columns -> empty, no files
    empty_df <- data.frame(Cluster = character(), ID = character(),
                           Description = character(), geneID = character())
    expect_equal(length(plot_ora_shared_genes(empty_df, tmp)), 0)
})

test_that("plot_ora_shared_genes writes genes_to_terms + terms_to_terms per cluster", {
    skip_if_not_installed("pheatmap")
    tmp <- withr::local_tempdir()
    # Two clusters, two terms each, overlapping genes -> non-degenerate matrices.
    df <- data.frame(
        Cluster     = c("1", "1", "2", "2"),
        ID          = c("GO:1", "GO:2", "GO:3", "GO:4"),
        Description = c("term one", "term two", "term three", "term four"),
        geneID      = c("g1/g2/g3", "g2/g3/g4", "g5/g6", "g6/g7/g8"),
        stringsAsFactors = FALSE
    )
    files <- plot_ora_shared_genes(df, tmp)
    # CSVs are always written (both views, both clusters) = 4 CSVs minimum.
    # Filenames carry the cluster label: cluster_<label>_{genes,terms}_to_terms.csv
    csvs <- list.files(tmp, pattern = "\\.csv$")
    expect_true(any(grepl("genes_to_terms", csvs)))
    expect_true(any(grepl("terms_to_terms", csvs)))
    expect_true(any(grepl("cluster_1_", csvs)))
    expect_true(all(file.info(file.path(tmp, csvs))$size > 0))
    # single-term cluster is skipped (legacy guard)
    df1 <- df[1, , drop = FALSE]
    expect_equal(length(plot_ora_shared_genes(df1, withr::local_tempdir())), 0)
})

test_that("plot_gsea_ridgeplot all-genes variant is fail-soft and uses distinct name", {
    tmp <- withr::local_tempdir()
    out_dir <- file.path(tmp, "ridgeplot")
    # core_enrichment = FALSE -> legacy "all genes" variant. Fail-soft on NULL,
    # and it must target plot_all_genes.png (not plot.png).
    expect_warning(expect_null(
        plot_gsea_ridgeplot(NULL, out_dir, core_enrichment = FALSE)))
    expect_false(file.exists(file.path(out_dir, "plot_all_genes.png")))
    expect_false(file.exists(file.path(out_dir, "plot.png")))
})

test_that("save_gsea_per_pathway_artifacts is fail-soft on NULL input", {
    tmp <- withr::local_tempdir()
    expect_null(save_gsea_per_pathway_artifacts(NULL, data.frame(), tmp))
    # No result, no output_dir -> no-op, no error
    expect_null(save_gsea_per_pathway_artifacts(NULL, data.frame(), NULL))
})

test_that(".go_semdata caches per ontology (no rebuild on hit)", {
    # Seed the cache with sentinels; a cache hit must return them WITHOUT calling
    # GOSemSim::godata() (keeps this fast and independent of GO.db). Seed "BP" too
    # so the NULL->"BP" mapping also hits the cache and never triggers a build.
    assign("CC", "SENTINEL_CC", envir = .enrich_semdata_cache)
    assign("BP", "SENTINEL_BP", envir = .enrich_semdata_cache)
    on.exit({
        rm("CC", envir = .enrich_semdata_cache)
        rm("BP", envir = .enrich_semdata_cache)
    }, add = TRUE)
    expect_identical(.go_semdata("CC"), "SENTINEL_CC")
    expect_identical(.go_semdata(NULL), "SENTINEL_BP")  # NULL -> "BP" key
})

test_that(".go_semdata builds GO Wang semData offline from GO.db (no OrgDb)", {
    skip_if_not_installed("GOSemSim")
    skip_if_not_installed("GO.db")
    # The crux of the legacy-parity fix: GO simplify's Wang semantic data is built
    # from the GO DAG in GO.db, with NO organism OrgDb and no network. Build the
    # smallest ontology (CC) to keep this affordable. (~tens of seconds.)
    if (exists("CC", envir = .enrich_semdata_cache)) rm("CC", envir = .enrich_semdata_cache)
    on.exit(if (exists("CC", envir = .enrich_semdata_cache))
                rm("CC", envir = .enrich_semdata_cache), add = TRUE)
    sem <- .go_semdata("CC")
    expect_false(is.null(sem))
    expect_s4_class(sem, "GOSemSimDATA")
    expect_identical(.go_semdata("CC"), sem)  # second call is a cache hit
})

test_that(".build_gene_context assembles DE stats + expression + z-scores", {
    skip_if_not(exists("build_rnaseq_summary_df"),
                "build_rnaseq_summary_df not available")
    genes <- paste0("g", 1:6)
    mk <- function(seed) {
        withr::with_seed(seed, data.frame(
            FeatureID      = genes,
            log2FoldChange = rnorm(6),
            pvalue         = runif(6),
            padj           = runif(6),
            stringsAsFactors = FALSE
        ))
    }
    de_tables <- list(A.vs.B = mk(1), A.vs.C = mk(2))
    expr <- matrix(withr::with_seed(3, rnorm(6 * 4)), nrow = 6,
                   dimnames = list(genes, paste0("S", 1:4)))
    pre <- list(expr_work = expr, meta = data.frame(SampleID = paste0("S", 1:4)))

    ctx <- .build_gene_context(de_tables, pre, list(p_cutoff = 0.05))
    expect_s3_class(ctx, "data.frame")
    expect_equal(names(ctx)[1], "FeatureID")           # keyed by feature id first
    expect_true(all(genes %in% ctx$FeatureID))
    # per-contrast DE columns present (from build_rnaseq_summary_df)
    expect_true(any(grepl("^linearFC\\.", names(ctx))))
    # per-sample expression + z-score columns appended
    expect_true(all(paste0("S", 1:4) %in% names(ctx)))
    expect_true(any(grepl("\\.zscore$", names(ctx))))
})

# ---------------------------------------------------------------------------
# Enrichment enable/disable gating (mod_rnaseq_pathway): the explicit
# enrichment.enabled switch + the existing pathway.enabled switch. Either being
# false must skip enrichment and return the empty-safe shape; an absent
# enrichment.enabled must keep running (backward compatibility).
# ---------------------------------------------------------------------------

# Minimal config builder for the gate. enrichment_enabled = NULL omits the key.
.mk_pw_cfg <- function(pathway_enabled = TRUE, enrichment_enabled = NULL) {
    enr <- list(annotation_dir = "")
    if (!is.null(enrichment_enabled)) enr$enabled <- enrichment_enabled
    list(modes = list(rna = list(
        pathway    = list(enabled = pathway_enabled),
        enrichment = enr
    )))
}
.empty_pw <- list(annotation = NULL, pathway_results = list(), plot_files = list())

test_that("enrichment.enabled: false skips enrichment (empty-safe return)", {
    res <- suppressMessages(mod_rnaseq_pathway(
        de_res = NULL, pre = NULL,
        config = .mk_pw_cfg(pathway_enabled = TRUE, enrichment_enabled = FALSE),
        out_dir = withr::local_tempdir()))
    expect_identical(res, .empty_pw)
})

test_that("pathway.enabled: false still skips enrichment (existing behavior)", {
    res <- suppressMessages(mod_rnaseq_pathway(
        de_res = NULL, pre = NULL,
        config = .mk_pw_cfg(pathway_enabled = FALSE, enrichment_enabled = NULL),
        out_dir = withr::local_tempdir()))
    expect_identical(res, .empty_pw)
})

test_that("absent enrichment.enabled still runs enrichment (backward compatible)", {
    # No `enabled` key + pathway on: the disable gate must NOT fire. With empty DE
    # tables the run proceeds past the gate to the downstream 'no DE tables' guard
    # (a different code path) — proving enrichment was NOT skipped by the switch.
    expect_warning(
        res <- suppressMessages(mod_rnaseq_pathway(
            de_res = list(tables = list()), pre = NULL,
            config = .mk_pw_cfg(pathway_enabled = TRUE, enrichment_enabled = NULL),
            out_dir = withr::local_tempdir())),
        "No DE tables")
    expect_identical(res, .empty_pw)
})

test_that("enrichment.enabled: true runs enrichment (reaches DE-table check)", {
    expect_warning(
        res <- suppressMessages(mod_rnaseq_pathway(
            de_res = list(tables = list()), pre = NULL,
            config = .mk_pw_cfg(pathway_enabled = TRUE, enrichment_enabled = TRUE),
            out_dir = withr::local_tempdir())),
        "No DE tables")
    expect_identical(res, .empty_pw)
})
