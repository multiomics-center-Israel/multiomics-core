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

# ===========================================================================
# Codex-review fixes (#1/#4/#5/#6/#7/#8/#9/#10)
# ===========================================================================

# ---- #1 + #4: annotation_dir resolution + missing-dir error ----------------

.mk_enr_cfg <- function(annotation_dir, project_dir = NULL) {
    list(
        project = list(dir = project_dir),
        modes = list(rna = list(
            pathway    = list(enabled = TRUE),
            de         = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5),
            enrichment = list(annotation_dir = annotation_dir)
        ))
    )
}
.call_local_enrich <- function(cfg) {
    suppressMessages(.run_local_enrichment(
        de_tables   = list(A.vs.B = data.frame(FeatureID = "g1")),
        feature_ids = "g1",
        enr_cfg     = cfg$modes$rna$enrichment,
        config      = cfg,
        out_dir     = withr::local_tempdir()))
}

test_that("#1 missing annotation dir raises a clear offline-enrichment error", {
    proj <- withr::local_tempdir()  # exists, but has no 'annot' subdir
    err <- tryCatch(.call_local_enrich(.mk_enr_cfg("annot", proj)),
                    error = function(e) conditionMessage(e))
    expect_match(err, "Local/offline enrichment was requested")
    expect_match(err, "does not exist")
    expect_match(err, 'annotation_dir: ""', fixed = TRUE)
    expect_match(err, "annot")  # the raw path is echoed
})

test_that("#4 relative annotation_dir resolves against project.dir, not CWD", {
    proj <- withr::local_tempdir()  # project dir WITHOUT the 'annot' subdir
    cwd  <- withr::local_tempdir()
    dir.create(file.path(cwd, "annot"))  # same-named dir under CWD must NOT win
    withr::local_dir(cwd)
    err <- tryCatch(.call_local_enrich(.mk_enr_cfg("annot", proj)),
                    error = function(e) conditionMessage(e))
    # Errors (proj/annot doesn't exist) instead of silently using cwd/annot,
    # and the reported resolved path is under the project dir.
    expect_match(err, "does not exist")
    expect_true(grepl(basename(proj), err, fixed = TRUE))
})

test_that("#4 absolute annotation_dir is used as-is (not prefixed by project.dir)", {
    proj <- withr::local_tempdir()
    abs_missing <- file.path(withr::local_tempdir(), "nope")  # absolute, missing
    err <- tryCatch(.call_local_enrich(.mk_enr_cfg(abs_missing, proj)),
                    error = function(e) conditionMessage(e))
    expect_match(err, "does not exist")
    expect_false(grepl(basename(proj), err, fixed = TRUE))  # project dir NOT prepended
})

# ---- #10: sequential fallback is seeded + reproducible ---------------------

test_that("#10 no-future fallback is seeded, reproducible, and non-leaking", {
    # Force the fallback by shadowing requireNamespace in the environment that
    # run_enrichment_jobs (sourced into globalenv) looks up from.
    had  <- exists("requireNamespace", envir = globalenv(), inherits = FALSE)
    orig <- if (had) get("requireNamespace", envir = globalenv()) else NULL
    assign("requireNamespace", function(...) FALSE, envir = globalenv())
    on.exit({
        if (had) assign("requireNamespace", orig, envir = globalenv())
        else rm("requireNamespace", envir = globalenv())
    }, add = TRUE)

    jobs <- as.list(1:6)
    fun  <- function(j) stats::runif(1)
    r1 <- run_enrichment_jobs(jobs, fun, workers = 1, seed = 123)
    r2 <- run_enrichment_jobs(jobs, fun, workers = 1, seed = 123)
    expect_identical(r1, r2)                       # reproducible across runs
    expect_false(identical(r1, run_enrichment_jobs(jobs, fun, 1, seed = 999)))

    # Does not leave the caller's global RNG state altered.
    set.seed(42); before <- .Random.seed
    invisible(run_enrichment_jobs(jobs, fun, workers = 1, seed = 7))
    expect_identical(.Random.seed, before)
})

# ---- #6: effective GSEA cutoff selection -----------------------------------

test_that("#6 .gsea_significant_rows uses the configured adjusted-p cutoff", {
    res_df <- data.frame(ID = c("a", "b", "c"), padj = c(0.03, 0.07, 0.20),
                         stringsAsFactors = FALSE)
    expect_equal(nrow(.gsea_significant_rows(res_df, 0.05)), 1L)   # only 0.03
    expect_setequal(.gsea_significant_rows(res_df, 0.10)$ID, c("a", "b"))  # +0.07
    expect_equal(nrow(.gsea_significant_rows(res_df, 0.01)), 0L)
    # NA padj rows never selected; empty/absent handled
    res_na <- data.frame(ID = "x", padj = NA_real_)
    expect_equal(nrow(.gsea_significant_rows(res_na, 0.10)), 0L)
    expect_equal(nrow(.gsea_significant_rows(data.frame(ID = "x"), 0.10)), 0L)
})

# ---- #8: Gene: prefix normalization for expression matching ----------------

test_that("#8 .build_gene_context matches DE ids to prefixed expr row names", {
    skip_if_not(exists("build_rnaseq_summary_df"))
    de <- list(A.vs.B = data.frame(
        FeatureID = c("ABC123", "DEF456"),
        log2FoldChange = c(1, -1), pvalue = c(0.01, 0.02), padj = c(0.02, 0.03),
        stringsAsFactors = FALSE))
    expr <- matrix(withr::with_seed(1, rnorm(2 * 4)), nrow = 2,
                   dimnames = list(c("Gene:ABC123", "Gene:DEF456"), paste0("S", 1:4)))
    ctx <- .build_gene_context(de, list(expr_work = expr), list(p_cutoff = 0.05))
    row <- ctx[ctx$FeatureID == "ABC123", , drop = FALSE]
    expect_false(any(is.na(row[, paste0("S", 1:4)])))   # expression was matched
    # No-prefix case still matches
    expr2 <- expr; rownames(expr2) <- c("ABC123", "DEF456")
    ctx2 <- .build_gene_context(de, list(expr_work = expr2), list(p_cutoff = 0.05))
    expect_false(any(is.na(ctx2[ctx2$FeatureID == "ABC123", paste0("S", 1:4)])))
})

test_that("#8 duplicate normalized expr keys are excluded (no silent wrong match)", {
    skip_if_not(exists("build_rnaseq_summary_df"))
    # "Gene:ABC123" and "ABC123" both normalize to ABC123 -> ambiguous. "DEF456" is
    # unique and must still match.
    de <- list(A.vs.B = data.frame(FeatureID = c("ABC123", "DEF456"),
                                   log2FoldChange = c(1, 1), pvalue = c(.01, .01),
                                   padj = c(.02, .02), stringsAsFactors = FALSE))
    expr <- matrix(withr::with_seed(1, rnorm(3 * 3)), nrow = 3,
                   dimnames = list(c("Gene:ABC123", "ABC123", "DEF456"),
                                   paste0("S", 1:3)))
    expect_warning(
        ctx <- .build_gene_context(de, list(expr_work = expr), list(p_cutoff = 0.05)),
        "ambiguous")
    # Ambiguous gene -> expression EXCLUDED (NA), not an arbitrary row.
    expect_true(all(is.na(ctx[ctx$FeatureID == "ABC123", paste0("S", 1:3)])))
    # Unique gene still matched.
    expect_false(any(is.na(ctx[ctx$FeatureID == "DEF456", paste0("S", 1:3)])))
})

# ---- #5 + #9: ORA dotplot toggle & targeted na-handling --------------------

# Small synthetic ORA result that is definitely significant (cluster1 == T1 genes).
.mk_ora_result <- function() {
    t2g <- data.frame(term = rep(c("T1", "T2", "T3"), each = 20),
                      gene = paste0("g", 1:60), stringsAsFactors = FALSE)
    t2n <- data.frame(term = c("T1", "T2", "T3"),
                      name = c("term one", "term two", "term three"),
                      stringsAsFactors = FALSE)
    clusters <- setNames(c(rep("1", 20), rep("2", 20)), paste0("g", 1:40))
    run_cluster_ora_compute(clusters, TERM2GENE = t2g, TERM2NAME = t2n,
                            type = "KEGG")   # KEGG: never simplified
}

test_that("#5 plots.dotplot=FALSE suppresses the ORA dotplot; TRUE writes it", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("enrichplot")
    res <- .mk_ora_result()
    skip_if(length(res) == 0, "fixture produced no ORA enrichment")

    d_off <- withr::local_tempdir()
    write_cluster_ora_outputs(res, unit_dir = d_off, type = "KEGG", dotplot = FALSE)
    expect_true(file.exists(file.path(d_off, "results.csv")))
    expect_false(file.exists(file.path(d_off, "dotplot.pdf")))   # suppressed

    d_on <- withr::local_tempdir()
    write_cluster_ora_outputs(res, unit_dir = d_on, type = "KEGG", dotplot = TRUE)
    expect_true(file.exists(file.path(d_on, "dotplot.pdf")))     # written
})

test_that("#9 shared-gene keeps enriched terms with NA in optional columns", {
    skip_if_not_installed("pheatmap")
    # 2 terms, overlapping genes; one term has qvalue = NA (an optional column).
    df <- data.frame(
        Cluster = c("1", "1"), ID = c("GO:1", "GO:2"),
        Description = c("term one", "term two"),
        geneID = c("g1/g2/g3", "g2/g3/g4"),
        qvalue = c(0.01, NA),   # unrelated NA must NOT drop the row
        stringsAsFactors = FALSE)
    tmp <- withr::local_tempdir()
    plot_ora_shared_genes(df, tmp)
    t2t <- list.files(tmp, pattern = "terms_to_terms.*\\.csv$", full.names = TRUE)
    expect_true(length(t2t) >= 1)
    m <- utils::read.csv(t2t[1], row.names = 1, check.names = FALSE)
    expect_equal(nrow(m), 2L)   # BOTH terms retained (old na.omit dropped NA row -> 1 -> skipped)
})

# ---- #7: per_pathway_artifacts gates core-gene CSV writes -------------------

.mk_gsea_result <- function() {
    genes <- paste0("g", 1:200)
    vals  <- c(seq(3, 2.1, length.out = 30), seq(2, -2, length.out = 140),
               seq(-2.1, -3, length.out = 30))
    stats <- setNames(vals, genes)
    t2g   <- data.frame(term = "T1", gene = paste0("g", 1:30), stringsAsFactors = FALSE)
    suppressWarnings(suppressMessages(clusterProfiler::GSEA(
        stats, TERM2GENE = t2g, minGSSize = 5, maxGSSize = 1000,
        pvalueCutoff = 1, verbose = FALSE, by = "fgsea")))
}

test_that("#7 core-gene CSVs written only when per_pathway_artifacts (plots) is on", {
    skip_if_not_installed("clusterProfiler")
    g <- tryCatch(.mk_gsea_result(), error = function(e) NULL)
    skip_if(is.null(g) || nrow(as.data.frame(g)) == 0, "no GSEA result produced")
    res_df <- as.data.frame(g); res_df$padj <- res_df$p.adjust
    expr <- matrix(withr::with_seed(1, rnorm(30 * 6)), nrow = 30,
                   dimnames = list(paste0("g", 1:30), paste0("S", 1:6)))
    core_csvs <- function(d) list.files(file.path(d, "core_genes"),
                                        pattern = "\\.csv$", full.names = TRUE)

    # artifacts FALSE, heatmaps FALSE -> nothing written
    d1 <- withr::local_tempdir()
    save_gsea_per_pathway_artifacts(g, res_df, d1, plots = FALSE, heatmaps = FALSE,
                                    sig_cutoff = 1)
    expect_equal(length(list.files(d1, recursive = TRUE)), 0L)

    # artifacts FALSE, heatmaps TRUE -> NO core-gene CSVs; heatmaps allowed
    skip_if_not_installed("pheatmap")
    d2 <- withr::local_tempdir()
    save_gsea_per_pathway_artifacts(g, res_df, d2, expr_mat = expr,
                                    plots = FALSE, heatmaps = TRUE, sig_cutoff = 1)
    expect_equal(length(core_csvs(d2)), 0L)                       # contract: no CSVs
    hm <- list.files(d2, recursive = TRUE)                        # relative paths
    expect_true(any(grepl("heatmaps_(all|core)_genes/.*\\.png$", hm)))  # heatmaps produced

    # artifacts TRUE, heatmaps FALSE -> core-gene CSVs written
    d3 <- withr::local_tempdir()
    save_gsea_per_pathway_artifacts(g, res_df, d3, plots = TRUE, heatmaps = FALSE,
                                    sig_cutoff = 1)
    expect_true(length(core_csvs(d3)) >= 1)

    # artifacts TRUE, heatmaps TRUE -> both
    d4 <- withr::local_tempdir()
    save_gsea_per_pathway_artifacts(g, res_df, d4, expr_mat = expr,
                                    plots = TRUE, heatmaps = TRUE, sig_cutoff = 1)
    expect_true(length(core_csvs(d4)) >= 1)
    expect_true(any(grepl("heatmaps_(all|core)_genes/.*\\.png$",
                          list.files(d4, recursive = TRUE))))
})

# ===========================================================================
# Mandatory GO simplification (config field removed; always applied for GO)
# ===========================================================================

test_that("GO simplify is mandatory for GO and never applied to KEGG (fail-soft)", {
    skip_if_not_installed("clusterProfiler")
    # Shadow .go_semdata so we can (a) observe whether simplify is attempted and
    # (b) force the fail-soft branch without a slow real godata() build.
    called <- new.env(); called$ont <- character(0)
    orig <- get(".go_semdata", envir = globalenv())
    assign(".go_semdata",
           function(ont) { called$ont <- c(called$ont, ont %||% "BP"); NULL },
           envir = globalenv())
    on.exit(assign(".go_semdata", orig, envir = globalenv()), add = TRUE)

    t2g <- data.frame(term = rep(c("T1", "T2", "T3"), each = 20),
                      gene = paste0("g", 1:60), stringsAsFactors = FALSE)
    t2n <- data.frame(term = c("T1", "T2", "T3"),
                      name = c("a", "b", "c"), stringsAsFactors = FALSE)
    clusters <- setNames(c(rep("1", 20), rep("2", 20)), paste0("g", 1:40))

    # GO: simplify is attempted unconditionally (no config flag exists). semData
    # NULL -> fail-soft: warn, keep unsimplified result (elements [[2]]/[[4]] NULL).
    called$ont <- character(0)
    expect_warning(
        res_go <- run_cluster_ora_compute(clusters, t2g, t2n, type = "GO", ont = "CC"),
        "GO simplify could not run")
    expect_true("CC" %in% called$ont)     # simplify attempted for GO (mandatory)
    expect_true(length(res_go) == 4)
    expect_false(is.null(res_go[[3]]))    # unsimplified result retained
    expect_null(res_go[[2]]); expect_null(res_go[[4]])   # no simplify output

    # KEGG: simplify is NEVER attempted.
    called$ont <- character(0)
    res_kegg <- run_cluster_ora_compute(clusters, t2g, t2n, type = "KEGG")
    expect_equal(length(called$ont), 0L)  # .go_semdata not called for KEGG
    expect_null(res_kegg[[2]]); expect_null(res_kegg[[4]])
})

test_that("run_cluster_ora_compute no longer accepts a go_simplify argument", {
    expect_false("go_simplify" %in% names(formals(run_cluster_ora_compute)))
    expect_false("go_simplify" %in% names(formals(run_cluster_ora)))
})

test_that("config template no longer exposes go_simplify", {
    tpl <- readLines(testthat::test_path("..", "..", "config", "templates",
                                         "rna_config.yaml"))
    # No ACTIVE (uncommented) go_simplify key remains.
    expect_false(any(grepl("^\\s*go_simplify\\s*:", tpl)))
})
