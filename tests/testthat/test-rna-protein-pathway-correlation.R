# tests/testthat/test-rna-protein-pathway-correlation.R
#
# compute_pathway_rna_protein_correlation() (06c): per-pathway Pearson and
# Spearman of RNA vs protein log2FC, kept only above min_pairs, ordered from
# most positive to most negative r. Synthetic data only.

make_pairs <- function(n, slope, prefix, contrast = "A_vs_B", seed = 1) {
    withr::with_seed(seed, {
        x <- seq(-2, 2, length.out = n)
        data.frame(gene_id = paste0("ENSG_", prefix, seq_len(n)),
                   protein_id = paste0("P_", prefix, seq_len(n)),
                   gene_symbol = paste0(prefix, seq_len(n)),
                   rna_log2FC = x,
                   protein_log2FC = slope * x + stats::rnorm(n, sd = 0.1),
                   contrast = contrast, stringsAsFactors = FALSE)
    })
}

pairs <- rbind(make_pairs(60, 1, "POS"), make_pairs(60, -1, "NEG", seed = 2),
               make_pairs(50, 1, "SMALL", seed = 3))
gene_sets <- list(hsaPOS = paste0("POS", 1:60), hsaNEG = paste0("NEG", 1:60),
                  hsaSMALL = paste0("SMALL", 1:50), hsaNONE = c("ABSENT1", "ABSENT2"))
attr(gene_sets, "descriptions") <- c(hsaPOS = "Positive set", hsaNEG = "Negative set",
                                     hsaSMALL = "Small set", hsaNONE = "Absent set")

test_that("pathways are ordered from most positive to most negative r", {
    res <- compute_pathway_rna_protein_correlation(pairs, gene_sets)
    expect_identical(res$stats$pathway_id, c("hsaPOS", "hsaNEG"))
    expect_gt(res$stats$pearson_r[1], 0.9)
    expect_lt(res$stats$pearson_r[2], -0.9)
    expect_identical(res$stats$pathway_name, c("Positive set", "Negative set"))
})

test_that("the default keeps only pathways with more than 50 pairs", {
    res <- compute_pathway_rna_protein_correlation(pairs, gene_sets)
    expect_false("hsaSMALL" %in% res$stats$pathway_id)
    expect_true(all(res$stats$n_pairs > 50))
    res50 <- compute_pathway_rna_protein_correlation(pairs, gene_sets, min_pairs = 50)
    expect_true("hsaSMALL" %in% res50$stats$pathway_id)
})

test_that("the statistics match cor.test on the member pairs", {
    res <- compute_pathway_rna_protein_correlation(pairs, gene_sets)
    pos <- pairs[startsWith(pairs$gene_symbol, "POS"), ]
    expect_equal(res$stats$pearson_r[1],
                 unname(stats::cor.test(pos$rna_log2FC, pos$protein_log2FC)$estimate))
    expect_equal(res$stats$spearman_rho[1],
                 unname(stats::cor(pos$rna_log2FC, pos$protein_log2FC, method = "spearman")))
    expect_equal(nrow(res$members[res$members$pathway_id == "hsaPOS", ]), 60L)
})

test_that("symbols match case-insensitively and NA log2FC pairs are dropped", {
    p <- pairs
    p$gene_symbol[startsWith(p$gene_symbol, "POS")] <- tolower(p$gene_symbol[startsWith(p$gene_symbol, "POS")])
    p$rna_log2FC[p$gene_symbol == "pos1"] <- NA
    res <- compute_pathway_rna_protein_correlation(p, gene_sets)
    expect_equal(res$stats$n_pairs[res$stats$pathway_id == "hsaPOS"], 59L)
})

test_that("contrasts are scored separately", {
    two <- rbind(make_pairs(60, 1, "POS"), make_pairs(60, -1, "POS", contrast = "C_vs_D", seed = 4))
    res <- compute_pathway_rna_protein_correlation(two, gene_sets["hsaPOS"])
    expect_identical(res$stats$contrast, c("A_vs_B", "C_vs_D"))
    expect_equal(sign(res$stats$pearson_r), c(1, -1))
})

test_that("no qualifying pathway or no pairs gives zero-row tables, not an error", {
    res <- compute_pathway_rna_protein_correlation(pairs, gene_sets["hsaNONE"])
    expect_equal(nrow(res$stats), 0L)
    expect_equal(nrow(res$members), 0L)
    expect_equal(nrow(compute_pathway_rna_protein_correlation(pairs[0, ], gene_sets)$stats), 0L)
})

test_that("pathway_gmt unset means the view is skipped", {
    expect_null(resolve_rna_protein_pathway_gmt(list(modes = list(multiomics = list()))))
    cfg <- list(modes = list(multiomics = list(rna_protein_correlation = list(pathway_gmt = c("a.gmt", "b.gmt")))))
    expect_error(resolve_rna_protein_pathway_gmt(cfg), "takes one GMT")
})

test_that("the report reads the written tables instead of recomputing", {
    rmd <- file.path(root_dir, "R", "domain", "multiomics", "report_template_multiomics.Rmd")
    src <- paste(readLines(rmd), collapse = "\n")
    expect_true(grepl("rna_protein_pathway_correlation.csv", src, fixed = TRUE))
    chunk <- regmatches(src, regexpr("```\\{r rna-prot-pathway-scatter[\\s\\S]*?\n```", src, perl = TRUE))
    expect_false(grepl("cor.test|cor\\(", chunk))
})
