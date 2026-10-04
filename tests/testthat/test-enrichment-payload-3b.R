# tests/testthat/test-enrichment-payload-3b.R
#
# Stage 3B — build_enrichment_payload(): compact, integer-indexed enrichment
# block from the Stage-3A primitives. Representation transform only; no biology
# recomputed. Deterministic synthetic pathway_res (no clusterProfiler).

# A synthetic rna_pathway_res mirroring the Stage-3A output shape. `g99` is a
# membership gene deliberately OUTSIDE the expression universe (background gene).
make_pathway_res_3b <- function() {
    gsea_kegg <- data.frame(
        ID = c("p1", "p2"), Description = c("P1", "P2"), NES = c(2, -1.5),
        setSize = c(5L, 4L), enrichmentScore = c(.6, -.5), pvalue = c(1e-3, 1e-2),
        p.adjust = c(1e-2, 3e-2), qvalue = c(1e-2, 3e-2), rank = c(3L, 7L),
        leading_edge = c("tags=40%", "tags=50%"),
        core_enrichment = c("g1/g2/g3", "g6/g7"),
        contrast = "c1", database = "KEGG", ranking_method = "fc",
        padj = c(1e-2, 3e-2), pathway = c("P1", "P2"), stringsAsFactors = FALSE)
    gsea_go <- data.frame(
        ID = "GO:0001", Description = "bp one", NES = 1.8, setSize = 6L,
        enrichmentScore = .55, pvalue = 1e-3, p.adjust = 2e-2, qvalue = 2e-2,
        rank = 4L, leading_edge = "tags=33%", core_enrichment = "g1/g8",
        contrast = "c1", database = "GO_BP", ranking_method = "fc",
        padj = 2e-2, pathway = "bp one", stringsAsFactors = FALSE)
    ora_kegg <- data.frame(
        Cluster = c("up", "down"), ID = c("p1", "p2"),
        Description = c("P1", "P2"), geneID = c("g1/g2", "g4/g5"), Count = c(2L, 2L),
        p.adjust = c(1e-2, 2e-2), padj = c(1e-2, 2e-2),
        pathway = c("P1", "P2"), stringsAsFactors = FALSE)
    ora_go <- data.frame(
        Cluster = "2", ID = "GO:0001", Description = "bp one", geneID = "g1/g8",
        Count = 2L, p.adjust = 1e-2, padj = 1e-2, pathway = "bp one",
        stringsAsFactors = FALSE)

    pathway_results <- list(
        cluster_ora = list(
            KEGG_contrasts_c1_ora          = ora_kegg,
            GO_BP_partition_k_ora          = ora_go,
            GO_BP_partition_k_ora_simplify = ora_go),
        c1 = list(KEGG_gsea_fc = gsea_kegg, GO_BP_gsea_fc = gsea_go))

    row_m <- function(a, d, g, it, ev, st, n, hs, sk) data.frame(
        analysis = a, database = d, group = g, item = it, evaluated = ev,
        status = st, n_significant = n, has_simplify = hs, storage_key = sk,
        stringsAsFactors = FALSE)
    manifest <- rbind(
        row_m("GSEA","KEGG","fc","c1",TRUE,"significant",2L,FALSE,"KEGG_gsea_fc"),
        row_m("GSEA","GO_BP","fc","c1",TRUE,"significant",1L,FALSE,"GO_BP_gsea_fc"),
        row_m("GSEA","GO_MF","fc","c1",FALSE,"failed",NA_integer_,FALSE,"GO_MF_gsea_fc"),
        row_m("ORA","KEGG","contrasts","c1",TRUE,"significant",2L,FALSE,"KEGG_contrasts_c1_ora"),
        row_m("ORA","GO_BP","partition","2",TRUE,"significant",1L,TRUE,"GO_BP_partition_k_ora"),
        row_m("ORA","GO_BP","partition","5",TRUE,"empty",0L,TRUE,"GO_BP_partition_k_ora"))

    row_i <- function(a, d, g, it, ct, sk, hs, sim) data.frame(
        analysis = a, database = d, group = g, item = it, container = ct,
        storage_key = sk, has_simplify = hs, simplify_key = sim,
        stringsAsFactors = FALSE)
    index <- rbind(
        row_i("GSEA","KEGG","fc","c1","c1","KEGG_gsea_fc",FALSE,NA_character_),
        row_i("GSEA","GO_BP","fc","c1","c1","GO_BP_gsea_fc",FALSE,NA_character_),
        row_i("ORA","KEGG","contrasts","c1","cluster_ora","KEGG_contrasts_c1_ora",FALSE,NA_character_),
        row_i("ORA","GO_BP","partition","k","cluster_ora","GO_BP_partition_k_ora",TRUE,"GO_BP_partition_k_ora_simplify"))

    list(
        annotation = NULL, pathway_results = pathway_results, plot_files = list(),
        enrichment_manifest = manifest, enrichment_index = index,
        gsea_rankings = list(fc = list(c1 = setNames(c(3,2.5,2,1,0,-1,-2), paste0("g",1:7)))),
        pathway_membership = list(
            KEGG  = list(p1 = c("g1","g2","g3","g4","g5"), p2 = c("g6","g7","g8","g9","g99")),
            GO_BP = list(`GO:0001` = c("g1","g8","g2"))))
}

GU <- paste0("g", 1:10)  # measured universe (expr_norm rownames)

# recursive class collector (for "no plots" / "no matrices" invariants)
collect_classes <- function(x, acc = character(0)) {
    acc <- c(acc, class(x))
    if (is.list(x)) for (e in x) acc <- collect_classes(e, acc)
    acc
}

test_that("build_enrichment_payload returns the approved top-level schema", {
    e <- build_enrichment_payload(make_pathway_res_3b(), gene_universe = GU)
    expect_identical(names(e), c("available", "gene_index", "config", "manifest",
                                 "ora", "gsea", "gsea_leading_edge",
                                 "gsea_rankings", "pathway_membership"))
    expect_true(isTRUE(e$available))
})

test_that("gene_index is unique, deterministic, universe-first, background appended", {
    e <- build_enrichment_payload(make_pathway_res_3b(), gene_universe = GU)
    gi <- e$gene_index
    expect_equal(anyDuplicated(gi), 0L)
    expect_identical(gi[seq_along(GU)], GU)          # measured universe first, native order
    expect_identical(setdiff(gi, GU), "g99")         # background gene appended
    expect_true(match("g99", gi) == length(gi))      # after the universe
    # deterministic: same input -> identical index
    e2 <- build_enrichment_payload(make_pathway_res_3b(), gene_universe = GU)
    expect_identical(gi, e2$gene_index)
})

test_that(".encode_genes fails explicitly on an unresolved gene (no silent drop)", {
    expect_error(.encode_genes(c("g1", "gZZZ"), GU), "unresolved")
    expect_identical(.encode_genes(c("g3", "g1", "g3"), GU), c(3L, 1L, 3L))  # order + dups
})

test_that("leading-edge round-trips exactly and is order-preserving", {
    e <- build_enrichment_payload(make_pathway_res_3b(), gene_universe = GU)
    le <- e$gsea_leading_edge$KEGG$fc$c1
    expect_identical(e$gene_index[le$p1], c("g1","g2","g3"))
    expect_identical(e$gene_index[le$p2], c("g6","g7"))
    expect_identical(e$gene_index[e$gsea_leading_edge$GO_BP$fc$c1$`GO:0001`], c("g1","g8"))
})

test_that("core_enrichment is removed after encoding; leading_edge summary kept", {
    e <- build_enrichment_payload(make_pathway_res_3b(), gene_universe = GU)
    tab <- e$gsea$KEGG$fc$c1
    expect_false("core_enrichment" %in% names(tab))
    expect_true("leading_edge" %in% names(tab))
})

test_that("rankings round-trip exactly (names, order, scores) and are stored once per ranking x contrast", {
    pr <- make_pathway_res_3b()
    e <- build_enrichment_payload(pr, gene_universe = GU)
    r <- e$gsea_rankings$fc$c1
    rebuilt <- setNames(r$score, e$gene_index[r$idx])
    expect_identical(rebuilt, pr$gsea_rankings$fc$c1)
    # not duplicated by database: top level is ranking method, then contrast — no db keys
    expect_identical(names(e$gsea_rankings), "fc")
    expect_identical(names(e$gsea_rankings$fc), "c1")
    expect_false(any(c("KEGG","GO_BP") %in% names(e$gsea_rankings$fc)))
})

test_that("pathway-membership round-trips exactly and is keyed by database only", {
    pr <- make_pathway_res_3b()
    e <- build_enrichment_payload(pr, gene_universe = GU)
    expect_identical(e$gene_index[e$pathway_membership$KEGG$p1], c("g1","g2","g3","g4","g5"))
    expect_identical(e$gene_index[e$pathway_membership$KEGG$p2], c("g6","g7","g8","g9","g99"))
    expect_identical(e$gene_index[e$pathway_membership$GO_BP$`GO:0001`], c("g1","g8","g2"))
    # keyed by db -> pathway (no ranking/contrast dimension)
    expect_setequal(names(e$pathway_membership), c("KEGG","GO_BP"))
})

test_that("ORA nests correctly; GO simplify preserved; KEGG has no fabricated simplify; all rows kept", {
    e <- build_enrichment_payload(make_pathway_res_3b(), gene_universe = GU)
    # contrasts -> keyed by contrast; partition -> single df
    expect_true(is.data.frame(e$ora$KEGG$contrasts_with_direction$c1))
    expect_equal(nrow(e$ora$KEGG$contrasts_with_direction$c1), 2L)   # all rows preserved
    expect_true("geneID" %in% names(e$ora$KEGG$contrasts_with_direction$c1))
    expect_true(is.data.frame(e$ora$GO_BP$partition))
    # GO simplify present; KEGG has none
    expect_true(is.data.frame(e$ora$GO_BP$simplify$partition))
    expect_null(e$ora$KEGG$simplify)
})

test_that("GSEA nests correctly and preserves all significant rows", {
    e <- build_enrichment_payload(make_pathway_res_3b(), gene_universe = GU)
    expect_true(is.data.frame(e$gsea$KEGG$fc$c1))
    expect_equal(nrow(e$gsea$KEGG$fc$c1), 2L)
    expect_true(is.data.frame(e$gsea$GO_BP$fc$c1))
})

test_that("manifest is carried through unchanged, incl. failed and evaluated-empty states", {
    pr <- make_pathway_res_3b()
    e <- build_enrichment_payload(pr, gene_universe = GU)
    expect_identical(e$manifest, pr$enrichment_manifest)
    expect_true(any(e$manifest$status == "failed" & !e$manifest$evaluated &
                    is.na(e$manifest$n_significant)))
    expect_true(any(e$manifest$status == "empty" & e$manifest$evaluated &
                    e$manifest$n_significant == 0))
    expect_true(any(e$manifest$status == "significant" & e$manifest$n_significant > 0))
})

test_that("no plot/image objects and no duplicated expr/de/annotation are embedded", {
    e <- build_enrichment_payload(make_pathway_res_3b(), gene_universe = GU)
    cls <- collect_classes(e)
    expect_false(any(c("ggplot","gg","pheatmap","recordedplot","gtable","trellis") %in% cls))
    expect_false(any(c("expr_norm","de_stats","feature_annot") %in% names(e)))
    expect_false("matrix" %in% cls)   # no expression matrices copied in
})

test_that("build_enrichment_payload returns NULL for absent/incompatible input", {
    expect_null(build_enrichment_payload(NULL))
    # online-shape: no enrichment_index / gsea_rankings
    expect_null(build_enrichment_payload(list(pathway_results = list(a = 1))))
    # empty enrichment
    expect_null(build_enrichment_payload(list(
        pathway_results = list(), enrichment_index = .empty_enrichment_index(),
        gsea_rankings = list())))
})
