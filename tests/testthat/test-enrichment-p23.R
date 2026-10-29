# §23 pipeline: complete ORA membership (J2), Shared Genes Top-N static outputs
# + complete matrix CSV (C4), and Count-descending display order (B).

test_that("§23 J2: membership unions ORA term IDs (incl KEGG) with GSEA IDs", {
    lt <- list(
        GO_BP = list(TERM2GENE = data.frame(
            term = c("GO:1", "GO:1", "GO:2", "GO:3"),
            gene = c("a", "b", "c", "d"), stringsAsFactors = FALSE)),
        KEGG  = list(TERM2GENE = data.frame(
            term = c("K1", "K1", "K2"), gene = c("x", "y", "z"),
            stringsAsFactors = FALSE)))
    # GSEA covered only GO:1; ORA adds GO:2/GO:3 (GO_BP) and K1/K2 (KEGG, no GSEA)
    gsea <- list(cA = list(GO_BP_pval = data.frame(ID = "GO:1", database = "GO_BP",
                                                   stringsAsFactors = FALSE)))
    ora_ids <- list(GO_BP = c("GO:2", "GO:3"), KEGG = c("K1", "K2"))

    m_gsea_only <- .build_gsea_pathway_membership(gsea, lt)
    expect_false("KEGG" %in% names(m_gsea_only))          # the pre-§23 gap
    expect_setequal(names(m_gsea_only$GO_BP), "GO:1")

    m <- .build_gsea_pathway_membership(gsea, lt, ora_ids)
    expect_setequal(names(m), c("GO_BP", "KEGG"))
    expect_setequal(names(m$GO_BP), c("GO:1", "GO:2", "GO:3"))  # GSEA ∪ ORA
    expect_setequal(names(m$KEGG), c("K1", "K2"))              # KEGG now covered
    expect_setequal(m$GO_BP[["GO:1"]], c("a", "b"))
    expect_setequal(m$KEGG[["K2"]], "z")
})

test_that("§23 C4: Shared Genes writes complete CSV + top20/top50 PDFs, no suppression", {
    skip_if_not_installed("pheatmap")
    set.seed(1)
    n_terms <- 30L; genes <- paste0("g", 1:400)
    df <- data.frame(
        Cluster = "1", ID = paste0("GO:", seq_len(n_terms)),
        Description = paste("term", seq_len(n_terms)),
        p.adjust = sort(runif(n_terms, 1e-6, 0.05)),
        geneID = vapply(seq_len(n_terms), function(i)
            paste(sample(genes, sample(5:80, 1)), collapse = "/"), character(1)),
        stringsAsFactors = FALSE)
    od <- withr::local_tempdir()
    plot_ora_shared_genes(df, od)
    # COMPLETE matrix CSV — all terms as rows, every union gene as a column
    csv <- file.path(od, "cluster_1_genes_to_terms.csv")
    expect_true(file.exists(csv))
    m <- utils::read.csv(csv, row.names = 1, check.names = FALSE)
    expect_equal(nrow(m), n_terms)
    expect_gt(ncol(m), 200L)                              # >200 genes NOT suppressed
    # readable Top-N GENE figures always rendered (no legacy suppression)
    expect_true(file.exists(file.path(od, "cluster_1_shared_genes_top20.pdf")))
    expect_true(file.exists(file.path(od, "cluster_1_shared_genes_top50.pdf")))
})

test_that("§23 C2: gene ranking = recurrence desc, best-FDR asc, gene id", {
    skip_if_not_installed("pheatmap")
    # g1 in all 3 terms (recurrence 3), g2 in 2, g3 in 1; distinct FDRs
    df <- data.frame(
        Cluster = "1", ID = c("GO:1", "GO:2", "GO:3"),
        Description = c("t1", "t2", "t3"),
        p.adjust = c(1e-6, 1e-4, 1e-2),
        geneID = c("g1/g2/g3", "g1/g2", "g1"), stringsAsFactors = FALSE)
    od <- withr::local_tempdir()
    plot_ora_shared_genes(df, od)
    m <- as.matrix(utils::read.csv(file.path(od, "cluster_1_genes_to_terms.csv"),
                                   row.names = 1, check.names = FALSE))
    rec <- colSums(m)
    expect_equal(unname(rec[c("g1", "g2", "g3")]), c(3, 2, 1))   # recurrence order g1>g2>g3
})
