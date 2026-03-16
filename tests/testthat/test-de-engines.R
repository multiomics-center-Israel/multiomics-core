# Tests for differential expression engines and ORA enrichment
# Validates run_limma_proteomics() and run_ora()

# ============================================================
# run_ora() — Over-Representation Analysis
# ============================================================

test_that("run_ora returns correct columns", {
    sig_genes <- c("GENE1", "GENE2", "GENE3")
    gene_sets <- list(
        PATHWAY_A = c("GENE1", "GENE2", "GENE4", "GENE5", "GENE6",
                       "GENE7", "GENE8", "GENE9", "GENE10", "GENE11"),
        PATHWAY_B = c("GENE20", "GENE21", "GENE22", "GENE23", "GENE24",
                       "GENE25", "GENE26", "GENE27", "GENE28", "GENE29")
    )
    background <- paste0("GENE", 1:100)

    res <- run_ora(sig_genes, gene_sets, background, min_size = 5)

    expect_true(is.data.frame(res))
    expect_true(all(c("pathway", "size", "overlap", "pvalue", "padj") %in% names(res)))
})

test_that("run_ora finds enrichment when overlap is high", {
    # 8 out of 10 sig genes are in PATHWAY_A
    sig_genes <- paste0("GENE", 1:10)
    gene_sets <- list(
        PATHWAY_A = paste0("GENE", 1:15),  # 10 of 10 sig genes overlap
        PATHWAY_B = paste0("GENE", 50:65)  # 0 overlap
    )
    background <- paste0("GENE", 1:200)

    res <- run_ora(sig_genes, gene_sets, background, min_size = 5)

    # PATHWAY_A should be significant
    pa <- res[res$pathway == "PATHWAY_A", ]
    expect_true(nrow(pa) == 1)
    expect_true(pa$overlap >= 10)
    expect_true(pa$pvalue < 0.05)

    # PATHWAY_B should have 0 overlap
    pb <- res[res$pathway == "PATHWAY_B", ]
    expect_true(pb$overlap == 0)
})

test_that("run_ora returns empty data.frame when no sets pass size filter", {
    sig_genes <- c("GENE1")
    gene_sets <- list(
        TINY = c("GENE1", "GENE2")  # size 2, below min_size=10
    )
    background <- paste0("GENE", 1:50)

    res <- run_ora(sig_genes, gene_sets, background, min_size = 10)
    expect_true(is.data.frame(res))
    expect_equal(nrow(res), 0)
})

test_that("run_ora filters gene sets by max_size", {
    sig_genes <- paste0("GENE", 1:5)
    gene_sets <- list(
        BIG = paste0("GENE", 1:600),    # too large
        OK  = paste0("GENE", 1:20)      # within limits
    )
    background <- paste0("GENE", 1:1000)

    res <- run_ora(sig_genes, gene_sets, background, min_size = 5, max_size = 500)
    expect_true("OK" %in% res$pathway)
    expect_false("BIG" %in% res$pathway)
})

test_that("run_ora padj uses BH correction", {
    sig_genes <- paste0("GENE", 1:5)
    # Create multiple pathways so BH correction is meaningful
    gene_sets <- lapply(1:5, function(i) {
        paste0("GENE", ((i - 1) * 15 + 1):(i * 15))
    })
    names(gene_sets) <- paste0("PATH_", LETTERS[1:5])
    # Ensure sig genes overlap with PATH_A
    gene_sets$PATH_A <- c(paste0("GENE", 1:5), paste0("GENE", 100:110))
    background <- paste0("GENE", 1:500)

    res <- run_ora(sig_genes, gene_sets, background, min_size = 5)
    if (nrow(res) > 1) {
        # padj should be >= pvalue (BH can only increase)
        expect_true(all(res$padj >= res$pvalue - 1e-15))
    }
})

# ============================================================
# run_limma_proteomics() — Limma DE
# ============================================================

test_that("run_limma_proteomics returns expected structure", {
    skip_if_not_installed("limma")

    set.seed(42)
    n_prot <- 50
    n_samples <- 8

    # Expression matrix (no NAs — already imputed)
    expr <- matrix(rnorm(n_prot * n_samples, mean = 20, sd = 2),
                   nrow = n_prot, ncol = n_samples,
                   dimnames = list(paste0("PROT", seq_len(n_prot)),
                                   paste0("S", seq_len(n_samples))))
    # Add signal to first 5 proteins in group B
    expr[1:5, 5:8] <- expr[1:5, 5:8] + 4

    meta <- data.frame(
        SampleID = paste0("S", seq_len(n_samples)),
        Group = rep(c("A", "B"), each = 4),
        stringsAsFactors = FALSE
    )
    rownames(meta) <- meta$SampleID

    contrasts_df <- data.frame(
        Factor = "Group",
        Numerator = "B",
        Denominator = "A",
        Contrast_name = "B_vs_A",
        stringsAsFactors = FALSE
    )

    prot_tbl <- data.frame(
        Protein.IDs = rownames(expr),
        Gene.Name = paste0("GENE", seq_len(n_prot)),
        stringsAsFactors = FALSE
    )

    cfg <- list(
        modes = list(
            proteomics = list(
                effects = list(samples = "SampleID", color = "Group"),
                id_columns = list(protein_id = "Protein.IDs", protein_annot = "Gene.Name"),
                de = list(fdrtool_correction = FALSE, p_adjust_method = "BH"),
                de_table = list(id_col = "FeatureID")
            )
        )
    )

    res <- run_limma_proteomics(expr, meta, contrasts_df, prot_tbl, cfg)

    expect_true(is.list(res))
    expect_true("de_tables" %in% names(res))
    expect_true("B_vs_A" %in% names(res$de_tables))

    de_tbl <- res$de_tables$B_vs_A
    expect_true(is.data.frame(de_tbl))
    expect_true(all(c("logFC", "P.Value", "adj.P.Val") %in% names(de_tbl)))
    expect_equal(nrow(de_tbl), n_prot)
})

test_that("run_limma_proteomics detects true signal", {
    skip_if_not_installed("limma")

    set.seed(123)
    n_prot <- 30
    n_samples <- 6

    expr <- matrix(rnorm(n_prot * n_samples, mean = 20, sd = 1),
                   nrow = n_prot, ncol = n_samples,
                   dimnames = list(paste0("P", seq_len(n_prot)),
                                   paste0("S", seq_len(n_samples))))
    # Strong signal in first 3 proteins
    expr[1:3, 4:6] <- expr[1:3, 4:6] + 6

    meta <- data.frame(
        SampleID = colnames(expr),
        Group = rep(c("Ctrl", "Treat"), each = 3),
        stringsAsFactors = FALSE
    )
    rownames(meta) <- meta$SampleID

    contrasts_df <- data.frame(
        Factor = "Group",
        Numerator = "Treat",
        Denominator = "Ctrl",
        Contrast_name = "Treat_vs_Ctrl",
        stringsAsFactors = FALSE
    )

    prot_tbl <- data.frame(
        Protein.IDs = rownames(expr),
        Gene.Name = paste0("G", seq_len(n_prot)),
        stringsAsFactors = FALSE
    )

    cfg <- list(
        modes = list(
            proteomics = list(
                effects = list(samples = "SampleID", color = "Group"),
                id_columns = list(protein_id = "Protein.IDs", protein_annot = "Gene.Name"),
                de = list(fdrtool_correction = FALSE, p_adjust_method = "BH"),
                de_table = list(id_col = "FeatureID")
            )
        )
    )

    res <- run_limma_proteomics(expr, meta, contrasts_df, prot_tbl, cfg)
    de_tbl <- res$de_tables$Treat_vs_Ctrl

    # The 3 signal proteins should have positive logFC and low p-values
    signal_rows <- de_tbl[de_tbl$FeatureID %in% c("P1", "P2", "P3"), ]
    expect_true(all(signal_rows$logFC > 2))
    expect_true(all(signal_rows$adj.P.Val < 0.05))
})
