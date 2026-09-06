# tests/testthat/test-enrichment-globals-fix.R
#
# Runtime scalability fix (found during Stage-3E full Assaf validation):
# the ORA/GSEA workers used to capture the ENTIRE multi-database `local_tables`,
# and `workers=1` still routed through future_lapply (hitting the
# future.globals.maxSize guard). The fix: each job carries only its own
# database's TERM2GENE/TERM2NAME, the worker captures nothing large, and ORA at
# workers<=1 uses a true base::lapply path. GSEA stays on future_lapply for RNG
# reproducibility. These tests lock in the behaviour and prove results/RNG are
# unchanged.

smib <- function(x) length(serialize(x, NULL)) / 1024 / 1024

# --- deterministic synthetic annotation tables (3 DBs so "all" >> "one") ---
make_local_tables_fix <- function() {
    kegg_g <- data.frame(term = rep(c("p1","p2"), each = 10),
                         gene = paste0("g", 1:20), stringsAsFactors = FALSE)
    kegg_n <- data.frame(term = c("p1","p2"), name = c("P1","P2"), stringsAsFactors = FALSE)
    # two moderate "GO-like" DBs to make the multi-DB total clearly larger
    big <- function(pref) data.frame(
        term = rep(paste0(pref, 1:300), each = 60),
        gene = paste0("bg", rep(1:60, 300)), stringsAsFactors = FALSE)
    bign <- function(pref) data.frame(term = paste0(pref, 1:300),
                                      name = paste0(pref, " term ", 1:300), stringsAsFactors = FALSE)
    list(
        KEGG  = list(TERM2GENE = kegg_g,        TERM2NAME = kegg_n),
        GO_BP = list(TERM2GENE = big("GO:BP"),  TERM2NAME = bign("GO:BP")),
        GO_MF = list(TERM2GENE = big("GO:MF"),  TERM2NAME = bign("GO:MF"))
    )
}

# Build ORA jobs the way .run_local_enrichment() does (single-DB per job).
build_ora_jobs_fix <- function(gene_lists, lt) {
    jobs <- list()
    for (m in names(gene_lists)) for (r in names(gene_lists[[m]])) {
        clusters <- gene_lists[[m]][[r]]
        for (db in names(lt)) jobs[[length(jobs) + 1]] <- list(
            clust_method = m, clust_round = r, db_name = db,
            clusters = clusters, term2gene = lt[[db]]$TERM2GENE,
            term2name = lt[[db]]$TERM2NAME)
    }
    jobs
}

# ===========================================================================
# 1 / 2 — workers no longer capture the full multi-database local_tables
# ===========================================================================

test_that("ORA worker captures only scalar thresholds (no local_tables / gene_lists)", {
    w <- .make_ora_worker(pval_cutoff = 0.05, padj_method = "fdr", orgdb = NULL)
    expect_setequal(ls(environment(w)), c("pval_cutoff", "padj_method", "orgdb"))
    expect_false(any(c("local_tables", "gene_lists") %in% ls(environment(w))))
    expect_lt(smib(w), 0.1)   # worker closure is tiny
})

test_that("GSEA worker captures only scalar thresholds (no local_tables)", {
    w <- .make_gsea_worker(pvalueCutoff = 0.05, pAdjustMethod = "fdr")
    expect_setequal(ls(environment(w)), c("pvalueCutoff", "pAdjustMethod"))
    expect_false("local_tables" %in% ls(environment(w)))
    expect_lt(smib(w), 0.1)
})

test_that("per-job export is a single database, far smaller than all local_tables", {
    lt <- make_local_tables_fix()
    w  <- .make_ora_worker(0.05, "fdr", NULL)
    gl <- list(partition = list(k = setNames(rep(c("1","2"), each = 10), paste0("g",1:20))))
    jobs <- build_ora_jobs_fix(gl, lt)
    # a GO_BP job exports only GO_BP's tables, not all three DBs
    go_job <- Filter(function(j) j$db_name == "GO_BP", jobs)[[1]]
    expect_true(all(c("clusters","term2gene","term2name") %in% names(go_job)))
    expect_lt(smib(list(w, go_job)), smib(lt))               # one DB < all DBs
    expect_lt(smib(go_job$term2gene), smib(lt$GO_BP$TERM2GENE) + 0.01)
    # and the job carries ONLY its own DB (COW ref identical to source)
    expect_identical(go_job$term2gene, lt$GO_BP$TERM2GENE)
})

# ===========================================================================
# 3 / 6 — ORA results identical (worker == direct call; and across paths)
# ===========================================================================

test_that("ORA worker(job) equals a direct run_cluster_ora_compute() call (results unchanged)", {
    skip_if_not_installed("clusterProfiler")
    lt <- make_local_tables_fix()
    clusters <- setNames(rep(c("1","2"), each = 10), paste0("g", 1:20))  # c1->p1, c2->p2
    job <- list(clust_method = "partition", clust_round = "k", db_name = "KEGG",
                clusters = clusters, term2gene = lt$KEGG$TERM2GENE,
                term2name = lt$KEGG$TERM2NAME)
    w <- .make_ora_worker(pval_cutoff = 0.5, padj_method = "fdr", orgdb = NULL)

    res_worker <- suppressWarnings(suppressMessages(w(job)$result))
    res_direct <- suppressWarnings(suppressMessages(run_cluster_ora_compute(
        clusters = clusters, TERM2GENE = lt$KEGG$TERM2GENE, TERM2NAME = lt$KEGG$TERM2NAME,
        type = "KEGG", pvalueCutoff = 0.5, pAdjustMethod = "fdr", orgdb = NULL, ont = NULL)))
    expect_identical(res_worker, res_direct)
})

test_that("ORA results identical across lapply / future-sequential / future-multisession", {
    skip_if_not_installed("clusterProfiler")
    lt <- make_local_tables_fix()
    gl <- list(partition = list(k = setNames(rep(c("1","2"), each = 10), paste0("g",1:20))))
    jobs <- build_ora_jobs_fix(gl, lt)
    w <- .make_ora_worker(0.5, "fdr", NULL)

    r_lapply <- suppressWarnings(suppressMessages(
        run_enrichment_jobs(jobs, w, workers = 1L, seed = 1L, prefer_lapply_when_sequential = TRUE)))
    r_futseq <- suppressWarnings(suppressMessages(
        run_enrichment_jobs(jobs, w, workers = 1L, seed = 1L, prefer_lapply_when_sequential = FALSE)))
    r_multi  <- suppressWarnings(suppressMessages(
        run_enrichment_jobs(jobs, w, workers = 2L, seed = 1L, prefer_lapply_when_sequential = TRUE)))

    # extract just the results (order-preserving) for comparison
    res <- function(x) lapply(x, function(j) j$result)
    expect_identical(res(r_lapply), res(r_futseq))     # lapply == future sequential
    expect_identical(res(r_lapply), res(r_multi))      # == future multisession (workers=1 vs >1)
})

# ===========================================================================
# 5 / 7 — ORA workers<=1 is TRUE sequential (no future global discovery);
#         the future path (GSEA-style) DOES go through future
# ===========================================================================

test_that("prefer_lapply_when_sequential bypasses future global discovery at workers<=1", {
    withr::local_options(future.globals.maxSize = 1024^2)   # 1 MiB, test-local (NOT a fix)
    big <- runif(400000)                                    # ~3 MiB captured global
    fn  <- local({ b <- big; function(job) length(b) })

    # GSEA-style (future path) trips the (lowered) globals guard even at workers=1
    expect_error(
        run_enrichment_jobs(list(1, 2), fn, workers = 1L, seed = 1L,
                            prefer_lapply_when_sequential = FALSE),
        "globals|maximum allowed size|maxSize")
    # ORA-style (base::lapply path) runs with no global check
    res <- run_enrichment_jobs(list(1, 2), fn, workers = 1L, seed = 1L,
                               prefer_lapply_when_sequential = TRUE)
    expect_equal(length(res), 2L)
})

# ===========================================================================
# 4 / 8 — GSEA stays on future path; results reproducible & worker-count invariant
# ===========================================================================

make_gsea_fixture <- function() {
    # pathway genes p1 = g1..g8 sit at the TOP of the ranking -> strong ES
    ranked <- setNames(c(seq(10, 3, length.out = 8), seq(2, -8, length.out = 22)),
                       paste0("g", 1:30))
    lt <- list(KEGG = list(
        TERM2GENE = data.frame(term = "p1", gene = paste0("g", 1:8), stringsAsFactors = FALSE),
        TERM2NAME = data.frame(term = "p1", name = "P1", stringsAsFactors = FALSE)))
    list(ranked = list(fc = list(c1 = ranked)), lt = lt)
}

test_that("GSEA is reproducible and worker-count invariant under the seeded future path", {
    skip_if_not_installed("clusterProfiler")
    fx <- make_gsea_fixture()
    run1  <- function() suppressWarnings(suppressMessages(
        run_gsea_all(fx$ranked, fx$lt, output_dir = NULL, workers = 1L, seed = 111L, dotplot = FALSE)))
    o1  <- run1()
    o1b <- run1()
    o2  <- suppressWarnings(suppressMessages(
        run_gsea_all(fx$ranked, fx$lt, output_dir = NULL, workers = 2L, seed = 111L, dotplot = FALSE)))

    expect_identical(o1$results, o1b$results)   # same seed -> identical (reproducible)
    expect_identical(o1$results, o2$results)     # workers=1 vs workers=2 -> identical (RNG preserved)
})
