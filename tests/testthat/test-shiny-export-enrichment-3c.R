# tests/testthat/test-shiny-export-enrichment-3c.R
#
# Stage 3C — integrate build_enrichment_payload() into the RNA Shiny export +
# contract. Backward-compatible, additive `payload$enrichment`.

# --- minimal but contract-valid RNA export inputs -------------------------
make_export_fixture <- function() {
    samples <- paste0("s", 1:4)
    genes   <- paste0("g", 1:10)
    mat <- matrix(as.numeric(1:40), nrow = 10, ncol = 4,
                  dimnames = list(genes, samples))
    meta <- data.frame(sample = samples, grp = c("A","A","B","B"),
                       row.names = samples, stringsAsFactors = FALSE)
    pre <- list(expr_filt = mat, expr_work = mat, meta = meta,
                row_data = data.frame(gene_id = genes, stringsAsFactors = FALSE))
    config <- list(modes = list(rna = list(
        de = list(padj_cutoff = 0.05, linear_fc_cutoff = 1.5),
        normalization = list(method = "TMMlogCPM"),
        effects = list(color = list("grp")),
        id_columns = list(gene_id = "gene_id"))))
    list(pre = pre, inputs = list(contrasts = NULL), config = config)
}

# --- minimal Stage-3A-shaped pathway_res over genes g1..g10 (+ background g99) ---
make_pathway_res_3c <- function() {
    gsea_kegg <- data.frame(
        ID = c("p1","p2"), Description = c("P1","P2"), NES = c(2,-1.5),
        setSize = c(5L,4L), enrichmentScore = c(.6,-.5), pvalue = c(1e-3,1e-2),
        p.adjust = c(1e-2,3e-2), qvalue = c(1e-2,3e-2), rank = c(3L,7L),
        leading_edge = c("tags=40%","tags=50%"),
        core_enrichment = c("g1/g2/g3","g6/g7"),
        contrast = "c1", database = "KEGG", ranking_method = "fc",
        padj = c(1e-2,3e-2), pathway = c("P1","P2"), stringsAsFactors = FALSE)
    ora_go <- data.frame(
        Cluster = "2", ID = "GO:0001", Description = "bp one", geneID = "g1/g8",
        Count = 2L, p.adjust = 1e-2, padj = 1e-2, pathway = "bp one",
        stringsAsFactors = FALSE)
    pathway_results <- list(
        cluster_ora = list(
            GO_BP_partition_k_ora = ora_go,
            GO_BP_partition_k_ora_simplify = ora_go),
        c1 = list(KEGG_gsea_fc = gsea_kegg))
    manifest <- rbind(
        data.frame(analysis="GSEA",database="KEGG",group="fc",item="c1",evaluated=TRUE,
                   status="significant",n_significant=2L,has_simplify=FALSE,
                   storage_key="KEGG_gsea_fc",stringsAsFactors=FALSE),
        data.frame(analysis="GSEA",database="GO_MF",group="fc",item="c1",evaluated=FALSE,
                   status="failed",n_significant=NA_integer_,has_simplify=FALSE,
                   storage_key="GO_MF_gsea_fc",stringsAsFactors=FALSE),
        data.frame(analysis="ORA",database="GO_BP",group="partition",item="2",evaluated=TRUE,
                   status="significant",n_significant=1L,has_simplify=TRUE,
                   storage_key="GO_BP_partition_k_ora",stringsAsFactors=FALSE))
    index <- rbind(
        data.frame(analysis="GSEA",database="KEGG",group="fc",item="c1",container="c1",
                   storage_key="KEGG_gsea_fc",has_simplify=FALSE,simplify_key=NA_character_,
                   stringsAsFactors=FALSE),
        data.frame(analysis="ORA",database="GO_BP",group="partition",item="k",
                   container="cluster_ora",storage_key="GO_BP_partition_k_ora",
                   has_simplify=TRUE,simplify_key="GO_BP_partition_k_ora_simplify",
                   stringsAsFactors=FALSE))
    list(annotation=NULL, pathway_results=pathway_results, plot_files=list(),
         enrichment_manifest=manifest, enrichment_index=index,
         gsea_rankings=list(fc=list(c1=setNames(c(3,2.5,2,1,0,-1,-2), paste0("g",1:7)))),
         pathway_membership=list(KEGG=list(p1=c("g1","g2","g3","g4","g5"),
                                           p2=c("g6","g7","g8","g9","g99"))))
}

collect_classes_3c <- function(x, acc = character(0)) {
    acc <- c(acc, class(x))
    if (is.list(x)) for (e in x) acc <- collect_classes_3c(e, acc)
    acc
}

build_payload <- function(pathway_res = NULL) {
    fx <- make_export_fixture()
    suppressWarnings(suppressMessages(build_shiny_payload_rnaseq(
        pre = fx$pre, de_res = NULL, inputs = fx$inputs, config = fx$config,
        pathway_res = pathway_res)))
}

test_that("RNA payload includes enrichment when valid enrichment results are supplied", {
    pr <- make_pathway_res_3c()
    p <- build_payload(pr)
    expect_true("enrichment" %in% names(p))
    expect_false(is.null(p$enrichment))
    expect_identical(names(p$enrichment),
                     c("available","gene_index","config","manifest","ora","gsea",
                       "gsea_leading_edge","gsea_rankings","pathway_membership"))
    # exported block == direct assembler output (export preserves it exactly)
    direct <- build_enrichment_payload(pr, gene_universe = rownames(make_export_fixture()$pre$expr_work),
                                       config = make_export_fixture()$config)
    expect_identical(p$enrichment, direct)
})

test_that("enrichment is NULL (key present) when enrichment was not run", {
    p <- build_payload(NULL)
    expect_true("enrichment" %in% names(p))       # key present...
    expect_null(p$enrichment)                     # ...but NULL (no fake structure)
})

test_that("payload_version remains '2.0' with and without enrichment", {
    expect_equal(build_payload(make_pathway_res_3c())$payload_version, "2.0")
    expect_equal(build_payload(NULL)$payload_version, "2.0")
})

test_that("enrichment is an optional canonical key; contract accepts present and absent", {
    expect_true("enrichment" %in% get_canonical_keys())
    expect_false("enrichment" %in% get_required_keys())
    # present + valid -> no error (strict)
    p <- build_payload(make_pathway_res_3c())
    expect_error(assert_shiny_payload_contract(p, strict = TRUE, context = "rnaseq"), NA)
    # enrichment removed ENTIRELY (old payload) -> still valid, no error (strict)
    p2 <- p; p2$enrichment <- NULL
    expect_error(assert_shiny_payload_contract(p2, strict = TRUE, context = "rnaseq"), NA)
    # enrichment present but NULL (no results) -> valid
    p3 <- build_payload(NULL)
    expect_error(assert_shiny_payload_contract(p3, strict = TRUE, context = "rnaseq"), NA)
})

test_that("exported GSEA tables keep padj + pathway and drop core_enrichment", {
    p <- build_payload(make_pathway_res_3c())
    tab <- p$enrichment$gsea$KEGG$fc$c1
    expect_true(all(c("padj","pathway") %in% names(tab)))
    expect_false("core_enrichment" %in% names(tab))
})

test_that("no plot/image objects and no duplicated expr/de/annot inside enrichment", {
    p <- build_payload(make_pathway_res_3c())
    cls <- collect_classes_3c(p$enrichment)
    expect_false(any(c("ggplot","gg","pheatmap","recordedplot","gtable","trellis") %in% cls))
    expect_false("matrix" %in% cls)
    expect_false(any(c("expr_norm","de_stats","feature_annot") %in% names(p$enrichment)))
})

test_that("rankings once per ranking x contrast; membership once per database", {
    e <- build_payload(make_pathway_res_3c())$enrichment
    expect_identical(names(e$gsea_rankings), "fc")
    expect_identical(names(e$gsea_rankings$fc), "c1")
    expect_false(any(c("KEGG","GO_BP") %in% names(e$gsea_rankings$fc)))
    expect_setequal(names(e$pathway_membership), "KEGG")
})

test_that("manifest statuses survive export unchanged (significant/empty/failed semantics)", {
    pr <- make_pathway_res_3c()
    e <- build_payload(pr)$enrichment
    expect_identical(e$manifest, pr$enrichment_manifest)
    expect_true(any(e$manifest$status == "failed" & !e$manifest$evaluated &
                    is.na(e$manifest$n_significant)))
    expect_true(any(e$manifest$status == "significant"))
})

test_that("validate_enrichment_payload accepts NULL and the assembled block, rejects core_enrichment", {
    expect_true(validate_enrichment_payload(NULL))
    e <- build_payload(make_pathway_res_3c())$enrichment
    expect_true(validate_enrichment_payload(e, strict = TRUE))
    bad <- e; bad$gsea$KEGG$fc$c1$core_enrichment <- "g1/g2"
    expect_error(validate_enrichment_payload(bad, strict = TRUE), "core_enrichment")
})
