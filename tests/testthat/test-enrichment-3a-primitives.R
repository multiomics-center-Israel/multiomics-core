# tests/testthat/test-enrichment-3a-primitives.R
#
# Stage 3A — preserve enrichment metadata + biological primitives that were
# previously discarded (availability manifest, storage index, exact GSEA
# rankings, GSEA pathway membership). The novel logic lives in pure-R helpers
# tested directly here; a gated KEGG-only integration test exercises the full
# .run_local_enrichment() return contract on tiny synthetic data.

# ===========================================================================
# Empty-schema helpers
# ===========================================================================

test_that("empty manifest/index have the canonical column schema", {
    m <- .empty_enrichment_manifest()
    expect_s3_class(m, "data.frame")
    expect_equal(nrow(m), 0L)
    expect_identical(colnames(m), c("analysis", "database", "group", "item",
                                    "evaluated", "status", "n_significant",
                                    "has_simplify", "storage_key"))

    idx <- .empty_enrichment_index()
    expect_equal(nrow(idx), 0L)
    expect_identical(colnames(idx), c("analysis", "database", "group", "item",
                                      "container", "storage_key", "has_simplify",
                                      "simplify_key"))
})

# ===========================================================================
# .build_ora_manifest_rows — three-state granularity
# ===========================================================================

test_that("partition manifest has one row per evaluated cluster, incl. n_significant = 0", {
    # Evaluated clusters 1, 2, 3 (from the gene-list); the result table only has
    # significant terms for clusters 1 and 2 -> cluster 3 must appear with 0.
    clusters <- setNames(c(1, 1, 2, 2, 3, 3), paste0("g", 1:6))
    main_df  <- data.frame(
        Cluster     = c("1", "1", "2"),
        ID          = c("path1", "path2", "path1"),
        Description = c("a", "b", "a"),
        stringsAsFactors = FALSE
    )
    rows <- .build_ora_manifest_rows(
        collection = "partition", round = "k", db_name = "GO_BP",
        main_df = main_df, has_simplify = TRUE, clusters = clusters,
        storage_key = "GO_BP_partition_k_ora")

    expect_setequal(rows$item, c("1", "2", "3"))
    expect_true(all(rows$analysis == "ORA"))
    expect_true(all(rows$group == "partition"))
    expect_true(all(rows$evaluated))
    expect_true(all(rows$has_simplify))
    expect_equal(rows$n_significant[rows$item == "1"], 2L)
    expect_equal(rows$n_significant[rows$item == "2"], 1L)
    expect_equal(rows$n_significant[rows$item == "3"], 0L)  # evaluated-but-empty
    expect_true(all(rows$storage_key == "GO_BP_partition_k_ora"))
    # ORA is always a successful evaluation: "significant" / "empty", never "failed"
    expect_equal(rows$status[rows$item == "1"], "significant")
    expect_equal(rows$status[rows$item == "3"], "empty")
})

test_that("binary_patterns manifest has one row per evaluated pattern", {
    clusters <- setNames(c("01", "01", "10", "11"), paste0("g", 1:4))
    main_df  <- data.frame(Cluster = c("01"), ID = "path1", stringsAsFactors = FALSE)
    rows <- .build_ora_manifest_rows(
        collection = "binary_patterns", round = "best", db_name = "KEGG",
        main_df = main_df, has_simplify = FALSE, clusters = clusters,
        storage_key = "KEGG_binary_patterns_best_ora")

    expect_setequal(rows$item, c("01", "10", "11"))
    expect_equal(rows$n_significant[rows$item == "01"], 1L)
    expect_equal(rows$n_significant[rows$item == "10"], 0L)
    expect_equal(rows$n_significant[rows$item == "11"], 0L)
    expect_false(any(rows$has_simplify))
})

test_that("contrasts / all_DE manifest is a single row keyed by the round", {
    main_df <- data.frame(Cluster = c("up", "down", "up"),
                          ID = c("a", "b", "c"), stringsAsFactors = FALSE)
    r1 <- .build_ora_manifest_rows(
        collection = "contrasts", round = "M_HFD.vs.NC", db_name = "GO_BP",
        main_df = main_df, has_simplify = TRUE, clusters = NULL,
        storage_key = "GO_BP_contrasts_M_HFD.vs.NC_ora")
    expect_equal(nrow(r1), 1L)
    expect_equal(r1$item, "M_HFD.vs.NC")
    expect_equal(r1$n_significant, 3L)   # total rows across up/down

    r2 <- .build_ora_manifest_rows(
        collection = "all_DE", round = "any_contrast", db_name = "KEGG",
        main_df = NULL, has_simplify = FALSE, clusters = NULL,
        storage_key = "KEGG_all_DE_any_contrast_ora")
    expect_equal(nrow(r2), 1L)
    expect_equal(r2$item, "any_contrast")
    expect_equal(r2$n_significant, 0L)   # evaluated but empty
})

# ===========================================================================
# .build_gsea_pathway_membership — union per db, subset of TERM2GENE
# ===========================================================================

test_that("pathway membership is the significant-GSEA union per db, all in TERM2GENE", {
    gsea_results <- list(
        M_HFD.vs.NC = list(
            GO_BP_gsea_fc = data.frame(
                ID = c("GO:0001", "GO:0002"), database = "GO_BP",
                stringsAsFactors = FALSE),
            KEGG_gsea_fc = data.frame(
                ID = "path1", database = "KEGG", stringsAsFactors = FALSE)
        ),
        M_Rev.vs.NC = list(
            GO_BP_gsea_pval_wo_direction = data.frame(
                ID = c("GO:0002", "GO:0003"), database = "GO_BP",
                stringsAsFactors = FALSE)
        )
    )
    local_tables <- list(
        GO_BP = list(TERM2GENE = data.frame(
            term = c("GO:0001", "GO:0001", "GO:0002", "GO:0003", "GO:0099"),
            gene = c("g1", "g2", "g3", "g4", "g9"), stringsAsFactors = FALSE)),
        KEGG = list(TERM2GENE = data.frame(
            term = c("path1", "path1", "path2"),
            gene = c("g1", "g5", "g6"), stringsAsFactors = FALSE))
    )

    mem <- .build_gsea_pathway_membership(gsea_results, local_tables)

    # Union of significant GSEA pathways per db (GO:0099 / path2 are NOT significant)
    expect_setequal(names(mem$GO_BP), c("GO:0001", "GO:0002", "GO:0003"))
    expect_setequal(names(mem$KEGG), "path1")
    # Exact membership from TERM2GENE
    expect_setequal(mem$GO_BP[["GO:0001"]], c("g1", "g2"))
    expect_setequal(mem$KEGG[["path1"]], c("g1", "g5"))
    # Every membership pathway exists in the corresponding TERM2GENE source
    for (db in names(mem)) {
        expect_true(all(names(mem[[db]]) %in% local_tables[[db]]$TERM2GENE$term))
    }
})

test_that("membership is empty when no GSEA pathway is significant", {
    expect_length(.build_gsea_pathway_membership(list(), list()), 0L)
})

# ===========================================================================
# GSEA empty-vs-failure semantics (Part 1)
# ===========================================================================

# Temporarily replace the GSEA worker factory in the environment where
# run_gsea_all() resolves it (the sourced functions live in globalenv). Restored
# on test exit. Lets us drive the empty-vs-failure branches WITHOUT clusterProfiler.
local_stub_gsea_worker <- function(gsea_result_for, env = parent.frame()) {
    target <- globalenv()
    orig <- get(".make_gsea_worker", envir = target)
    stub <- function(pvalueCutoff, pAdjustMethod) {
        function(job) list(ranking_method = job$ranking_method,
                           contrast = job$contrast, db_name = job$db_name,
                           gsea_result = gsea_result_for)
    }
    assign(".make_gsea_worker", stub, envir = target)
    withr::defer(assign(".make_gsea_worker", orig, envir = target), envir = env)
}

test_that("run_gsea_all records a TECHNICAL FAILURE distinctly from a successful empty result", {
    # Every job returns a caught error (class "gsea_error") — a technical failure.
    local_stub_gsea_worker(structure(list(message = "stub failure"), class = "gsea_error"))

    ranked <- list(fc = list(c1 = setNames(c(2, 1), c("g1", "g2"))))
    lt <- list(KEGG = list(
        TERM2GENE = data.frame(term = "p1", gene = c("g1", "g2"), stringsAsFactors = FALSE),
        TERM2NAME = data.frame(term = "p1", name = "P1", stringsAsFactors = FALSE)))

    out <- suppressMessages(run_gsea_all(ranked, lt, output_dir = NULL, workers = 1))
    m <- out$manifest
    expect_equal(nrow(m), 1L)
    expect_equal(m$status, "failed")
    expect_false(m$evaluated)              # -> hidden, NOT greyed as "no significant results"
    expect_true(is.na(m$n_significant))    # unknown, never a fabricated 0
    expect_equal(nrow(out$index), 0L)      # nothing stored for a failed unit
    expect_length(out$results, 0L)
})

test_that("a successful empty GSEA unit is status='empty' with n_significant = 0 (not failed)", {
    # A NON-error, NON-null result whose as.data.frame() has 0 rows — the
    # successful-but-nothing-significant case (a 0-row data.frame stands in for an
    # empty gseaResult without needing the S4 class loaded).
    local_stub_gsea_worker(data.frame())

    ranked <- list(fc = list(c1 = setNames(c(2, 1), c("g1", "g2"))))
    lt <- list(KEGG = list(
        TERM2GENE = data.frame(term = "p1", gene = c("g1", "g2"), stringsAsFactors = FALSE),
        TERM2NAME = data.frame(term = "p1", name = "P1", stringsAsFactors = FALSE)))

    out <- suppressMessages(run_gsea_all(ranked, lt, output_dir = NULL, workers = 1))
    m <- out$manifest
    expect_equal(nrow(m), 1L)
    expect_equal(m$status, "empty")
    expect_true(m$evaluated)               # -> greyed "no significant results"
    expect_equal(m$n_significant, 0L)
})

# ===========================================================================
# Integration — .run_local_enrichment() return contract (KEGG only, gated)
# ===========================================================================

make_offline_fixture <- function() {
    dir <- withr::local_tempdir(.local_envir = parent.frame())
    ann <- file.path(dir, "annotation"); dir.create(ann)
    # KEGG only (no GO -> no GO.db/GOSemSim dependency, no simplify).
    writeLines(c("term\tgene",
                 paste0("path1\t", c("g1", "g2", "g3", "g4", "g5")),
                 paste0("path2\t", c("g6", "g7", "g8", "g9", "g10"))),
               file.path(ann, "KEGG_pathway2gene.tab"))
    writeLines(c("term\tname", "path1\tPathway One", "path2\tPathway Two"),
               file.path(ann, "KEGG_pathway2name.tab"))

    genes <- paste0("g", 1:20)
    # c1: g1..g10 strongly significant & up; c2: nothing significant.
    c1 <- data.frame(
        FeatureID = genes,
        log2FoldChange = c(rep(3, 10), rep(0.1, 10)),
        pvalue = c(rep(1e-8, 10), rep(0.9, 10)),
        padj   = c(rep(1e-6, 10), rep(0.99, 10)),
        stringsAsFactors = FALSE)
    c2 <- data.frame(
        FeatureID = genes, log2FoldChange = rep(0.05, 20),
        pvalue = rep(0.8, 20), padj = rep(0.99, 20), stringsAsFactors = FALSE)
    de_tables <- list(c1 = c1, c2 = c2)

    clustering_res <- list(
        objects = list(patterns = data.frame(
            feature_id = paste0("g", 1:8),
            best_pattern = c("01", "01", "10", "10", "11", "11", "00", "00"),
            stringsAsFactors = FALSE)),
        excel_order = list(partition_clusters =
            setNames(c(1, 1, 2, 2, 3, 3), paste0("g", 1:6))))

    config <- list(
        project = list(dir = dir),
        params  = list(seed = 1L),
        modes = list(rna = list(
            de = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5),
            enrichment = list(
                annotation_dir = ann, databases = "KEGG",
                pvalue_cutoff = 0.05, padj_method = "fdr", workers = 1,
                plots = list(dotplot = FALSE, ridgeplot = FALSE,
                             ridgeplot_all_genes = FALSE, shared_genes = FALSE,
                             pathway_heatmaps = FALSE),
                gsea = list(per_pathway_artifacts = FALSE)))))

    list(de_tables = de_tables, feature_ids = genes, config = config,
         clustering_res = clustering_res, out_dir = file.path(dir, "out"))
}

test_that(".run_local_enrichment returns the four additive sibling keys with correct shapes", {
    skip_if_not_installed("clusterProfiler")
    fx <- make_offline_fixture()

    res <- suppressWarnings(suppressMessages(.run_local_enrichment(
        de_tables      = fx$de_tables,
        feature_ids    = fx$feature_ids,
        enr_cfg        = fx$config$modes$rna$enrichment,
        config         = fx$config,
        out_dir        = fx$out_dir,
        clustering_res = fx$clustering_res,
        pre            = NULL)))

    # (1) existing keys still present and unchanged in kind
    expect_true(all(c("annotation", "pathway_results", "plot_files") %in% names(res)))
    expect_type(res$pathway_results, "list")
    # (1) four new additive sibling keys
    expect_true(all(c("enrichment_manifest", "enrichment_index",
                      "gsea_rankings", "pathway_membership") %in% names(res)))

    # (2) new metadata are SIBLINGS, not nested inside pathway_results
    expect_false(any(c("enrichment_manifest", "enrichment_index",
                       "gsea_rankings", "pathway_membership") %in%
                     names(res$pathway_results)))

    # manifest schema + evaluated GSEA units present (incl. any zero-significant)
    m <- res$enrichment_manifest
    expect_identical(colnames(m), colnames(.empty_enrichment_manifest()))
    expect_true(all(m$status %in% c("significant", "empty", "failed")))
    expect_true(all(m$evaluated == (m$status != "failed")))
    # successful rows carry a real count; failed rows carry NA (never a fake 0)
    expect_true(all(m$n_significant[m$status != "failed"] >= 0))
    expect_true(all(is.na(m$n_significant[m$status == "failed"])))
    expect_true(any(m$analysis == "GSEA"))

    # (6) c2 has zero DE genes -> NO ORA contrasts row for it, but GSEA still evaluates it
    ora_rows <- m[m$analysis == "ORA", ]
    expect_false("c2" %in% ora_rows$item[ora_rows$group == "contrasts"])
    expect_true("c2" %in% m$item[m$analysis == "GSEA"])

    # index rows point at real stored objects
    idx <- res$enrichment_index
    expect_identical(colnames(idx), colnames(.empty_enrichment_index()))
    for (i in seq_len(nrow(idx))) {
        obj <- res$pathway_results[[idx$container[i]]][[idx$storage_key[i]]]
        expect_false(is.null(obj))
    }

    # (8) rankings are EXACTLY build_ranked_gene_lists(de_tables) — shared, not recomputed
    expect_identical(res$gsea_rankings, build_ranked_gene_lists(fx$de_tables))

    # (9,10) membership: only significant GSEA pathways, all present in TERM2GENE
    t2g <- load_local_pathway_tables(fx$config$modes$rna$enrichment$annotation_dir,
                                     databases = "KEGG")$KEGG$TERM2GENE
    for (db in names(res$pathway_membership)) {
        expect_true(all(names(res$pathway_membership[[db]]) %in% t2g$term))
    }
})
