# One join for RNA-protein pairing, and the artifacts that make it checkable.
#
# Two code paths used to pair RNA with protein independently and disagree on how
# many pairs exist, so the number a reader got depended on which table they
# opened. The difference was every gene whose Protein.Group is a
# semicolon-separated set of accessions: 06b matched feature ids through a
# long-form mapping keyed on gene_symbol and lost them, while 06_concordance.R
# merged on the mapping directly and kept them. build_rna_protein_pairs() is now
# the single implementation and both callers use it.
#
# Two further defects in the same path, both of which silently emptied the join
# rather than failing:
#
#   linearFC is a SIGNED linear ratio, so log2() of a down-regulated feature is
#   NaN. Half the proteome disappeared before the join could see it.
#
#   .mae_to_legacy() built its merge key from gene_id on the RNA side and
#   protein_id on the protein side when the mapping file carried no gene_symbol
#   column. The two sides then keyed on values that existed nowhere else, and
#   the join matched nothing at all for any non-model organism.
#
# All fixtures synthetic.

# Twelve pairs clears compute_de_concordance()'s ten-row floor. The fold changes
# are fixed rather than sampled: this file asserts that two code paths return
# the same numbers, which is only meaningful if the numbers are deterministic.
pair_fixture <- function(n = 12, semicolon_protein = TRUE) {
    genes <- sprintf("EHI_%03d", seq_len(n))
    prots <- sprintf("XP_%d", seq_len(n))
    # A Protein.Group naming several accessions -- the shape the gene_symbol
    # path used to drop.
    if (semicolon_protein) prots[1] <- "XP_1;XP_101;XP_102"

    rna_lfc  <- c(-1.5, 1.2, -0.8, 0.3, 2.1, -2.0, 0.9, -1.1, 0.2, 1.7, -0.4, 1.4)
    prot_lfc <- c(-2.0, 0.7, -1.4, -0.1, 1.3, -0.5, 1.6, -1.9, 1.0, 0.4, -0.6, 0.8)

    list(
        mapping = data.frame(gene_id = genes, protein_id = prots,
                             mapping_source = "custom_file",
                             stringsAsFactors = FALSE),
        rna_de = data.frame(feature_id = genes, logFC = rna_lfc[seq_len(n)],
                            padj = rep(0.2, n), stringsAsFactors = FALSE),
        prot_de = data.frame(feature_id = prots, logFC = prot_lfc[seq_len(n)],
                             padj = rep(0.2, n), stringsAsFactors = FALSE)
    )
}


# ---- the single join --------------------------------------------------------

test_that("a protein group naming several accessions keeps its pair", {
    fx <- pair_fixture()

    pairs <- build_rna_protein_pairs(fx$rna_de, fx$prot_de, fx$mapping)

    expect_equal(nrow(pairs), 12L)
    expect_true("XP_1;XP_101;XP_102" %in% pairs$protein_id)
    # And it carries the gene it belongs to, not a collapsed or blank id.
    expect_identical(pairs$gene_id[pairs$protein_id == "XP_1;XP_101;XP_102"],
                     "EHI_001")
})

test_that("a gene pairing with several protein groups keeps every pair", {
    # Isoform-level grouping, not duplication: collapsing it would silently pick
    # one protein per gene and change the correlation the report reports.
    mapping <- data.frame(
        gene_id    = c("EHI_001", "EHI_001", "EHI_002"),
        protein_id = c("XP_1", "XP_2", "XP_3"),
        stringsAsFactors = FALSE)
    rna_de <- data.frame(feature_id = c("EHI_001", "EHI_002"),
                         logFC = c(1.5, -0.5), padj = c(0.01, 0.4),
                         stringsAsFactors = FALSE)
    prot_de <- data.frame(feature_id = c("XP_1", "XP_2", "XP_3"),
                          logFC = c(1.1, 0.9, -0.2), padj = c(0.02, 0.03, 0.5),
                          stringsAsFactors = FALSE)

    pairs <- build_rna_protein_pairs(rna_de, prot_de, mapping)

    expect_equal(nrow(pairs), 3L)
    expect_equal(sum(pairs$gene_id == "EHI_001"), 2L)
    expect_setequal(pairs$protein_id, c("XP_1", "XP_2", "XP_3"))
})

test_that("mapping_source rides along when the mapping carries it, and not otherwise", {
    fx <- pair_fixture()

    with_source <- build_rna_protein_pairs(fx$rna_de, fx$prot_de, fx$mapping)
    expect_true("mapping_source" %in% names(with_source))
    expect_true(all(with_source$mapping_source == "custom_file"))

    bare <- fx$mapping[, c("gene_id", "protein_id")]
    expect_false("mapping_source" %in%
                     names(build_rna_protein_pairs(fx$rna_de, fx$prot_de, bare)))
})

test_that("a gene_symbol the mapping supplies survives the join and the export", {
    # The concordance table already exported gene_symbol wherever the mapping
    # file carried one. Narrowing the join to the two ids would have dropped a
    # column projects read, so it is carried through both paths.
    fx <- pair_fixture()
    fx$mapping$gene_symbol <- paste0("SYM", seq_len(nrow(fx$mapping)))

    pairs <- build_rna_protein_pairs(fx$rna_de, fx$prot_de, fx$mapping)
    expect_true("gene_symbol" %in% names(pairs))
    expect_identical(pairs$gene_symbol[pairs$gene_id == "EHI_001"], "SYM1")

    res <- suppressMessages(compute_de_concordance(
        fx$rna_de, fx$prot_de, mapping = fx$mapping, out_dir = NULL))
    expect_identical(res$gene_symbol[res$gene_id == "EHI_001"], "SYM1")

    # With no symbol in the mapping it falls back to the gene id rather than
    # to NA, which is what the report template expects to display.
    plain <- suppressMessages(compute_de_concordance(
        fx$rna_de, fx$prot_de, mapping = fx$mapping[, c("gene_id", "protein_id")],
        out_dir = NULL))
    expect_identical(plain$gene_symbol, plain$gene_id)
})

test_that("nothing to join is zero rows, not an error", {
    fx <- pair_fixture()
    empty <- fx$mapping[0, , drop = FALSE]

    expect_equal(nrow(build_rna_protein_pairs(fx$rna_de, fx$prot_de, empty)), 0L)
    expect_equal(nrow(build_rna_protein_pairs(fx$rna_de, fx$prot_de, NULL)), 0L)

    # A mapping whose ids match neither table joins nothing either.
    unrelated <- data.frame(gene_id = "NOPE_1", protein_id = "NOPE_2",
                            stringsAsFactors = FALSE)
    expect_equal(nrow(build_rna_protein_pairs(fx$rna_de, fx$prot_de, unrelated)), 0L)
})

test_that("both callers return the same pairs and the same values", {
    # The regression. These two used to join independently; whichever table a
    # reader opened decided the pair count and the correlation they quoted.
    fx <- pair_fixture()

    via_concordance <- analyze_rna_protein_concordance(
        fx$rna_de, fx$prot_de, fx$mapping, config = list(), out_dir = NULL)
    tbl <- via_concordance$concordance_table

    via_correlation <- suppressMessages(compute_de_concordance(
        fx$rna_de, fx$prot_de, mapping = fx$mapping, out_dir = NULL))

    key <- function(g, p) paste(g, p, sep = "||")
    expect_setequal(key(tbl$gene_id, tbl$protein_id),
                    key(via_correlation$gene_id, via_correlation$protein_id))

    # Same pairs AND the same numbers on them, compared pair by pair rather than
    # in whatever order each path happens to return.
    ord_a <- order(key(tbl$gene_id, tbl$protein_id))
    ord_b <- order(key(via_correlation$gene_id, via_correlation$protein_id))
    expect_equal(tbl$logFC_rna[ord_a], via_correlation$rna_log2FC[ord_b])
    expect_equal(tbl$logFC_prot[ord_a], via_correlation$protein_log2FC[ord_b])
})

test_that("the concordance path still emits the legacy column names", {
    # The report template and other projects read rna_log2FC / protein_log2FC /
    # gene_symbol. The identifiers were added beside them, not instead of them.
    fx <- pair_fixture()

    res <- suppressMessages(compute_de_concordance(
        fx$rna_de, fx$prot_de, mapping = fx$mapping, out_dir = NULL))

    expect_true(all(c("gene_id", "protein_id", "gene_symbol", "rna_log2FC",
                      "rna_padj", "protein_log2FC", "protein_padj",
                      "category", "concordant", "te_log2FC") %in% names(res)))
    expect_identical(res$gene_symbol, res$gene_id)
})


test_that("the mapping's own symbol wins over one a DE table happens to carry", {
    # merge() suffixes only the columns both frames share, so a DE table with
    # its own gene_symbol used to push the mapping's copy to gene_symbol.x
    # while leaving the other side's under the bare name -- and the bare one
    # was what got exported. Which side's symbol survives must not depend on
    # what the other side happens to hold.
    fx <- pair_fixture()
    fx$mapping$gene_symbol <- paste0("MAP", seq_len(nrow(fx$mapping)))
    fx$rna_de$gene_symbol  <- paste0("RNA", seq_len(nrow(fx$rna_de)))
    fx$prot_de$gene_symbol <- paste0("PROT", seq_len(nrow(fx$prot_de)))

    pairs <- build_rna_protein_pairs(fx$rna_de, fx$prot_de, fx$mapping)

    expect_true(all(grepl("^MAP", pairs$gene_symbol)))
    expect_identical(pairs$gene_symbol[pairs$gene_id == "EHI_001"], "MAP1")
})

test_that("limma and edgeR adjusted-p spellings reach the join", {
    # compute_de_concordance() accepted padj, adj.P.Val and FDR. Losing an
    # alias does not fail -- it marks every feature non-significant and
    # computes the significant-union result from the wrong set.
    fx <- pair_fixture()
    limma_rna <- fx$rna_de
    names(limma_rna)[names(limma_rna) == "padj"] <- "adj.P.Val"
    limma_rna$adj.P.Val <- 0.001
    edger_prot <- fx$prot_de
    names(edger_prot)[names(edger_prot) == "padj"] <- "FDR"
    edger_prot$FDR <- 0.002

    pairs <- build_rna_protein_pairs(limma_rna, edger_prot, fx$mapping)

    expect_false(any(is.na(pairs$padj_rna)))
    expect_false(any(is.na(pairs$padj_prot)))
    expect_true(all(pairs$padj_rna == 0.001))
    expect_true(all(pairs$padj_prot == 0.002))
})

test_that("an alias present on only one side still lands on that side", {
    # The hazard is asymmetric input: with padj on one table and adj.P.Val on
    # the other, neither column is shared, so neither used to pick up a side
    # suffix and both went missing.
    fx <- pair_fixture()
    limma_rna <- fx$rna_de
    names(limma_rna)[names(limma_rna) == "padj"] <- "adj.P.Val"
    limma_rna$adj.P.Val <- 0.01

    pairs <- build_rna_protein_pairs(limma_rna, fx$prot_de, fx$mapping)

    expect_true(all(pairs$padj_rna == 0.01))
    expect_true(all(pairs$padj_prot == 0.2))
})

test_that("a DESeq2-style fold-change column is found too", {
    fx <- pair_fixture()
    deseq_rna <- fx$rna_de
    names(deseq_rna)[names(deseq_rna) == "logFC"] <- "log2FoldChange"

    pairs <- build_rna_protein_pairs(deseq_rna, fx$prot_de, fx$mapping)

    expect_false(any(is.na(pairs$logFC_rna)))
    expect_equal(pairs$logFC_rna[pairs$gene_id == "EHI_001"], -1.5)
})

test_that("no usable fold change on a side is an error, not silent NAs", {
    fx <- pair_fixture()
    no_fc <- fx$rna_de[, c("feature_id", "padj")]

    expect_error(build_rna_protein_pairs(no_fc, fx$prot_de, fx$mapping),
                 "logFC")
})


# ---- the signed linear fold change ------------------------------------------

test_that("a negative linear fold change becomes a negative log2, not NaN", {
    # linearFC stores -1/2^x for a decrease. log2() of that is NaN, which took
    # every down-regulated protein out of the join without a word.
    expect_equal(signed_linear_fc_to_log2(c(4, -4)), c(2, -2))
    expect_equal(signed_linear_fc_to_log2(c(1, -1)), c(0, 0))
    expect_false(any(is.nan(signed_linear_fc_to_log2(c(-2, -8, -1.5)))))
})

test_that("values with no log2 are NA rather than an error or a warning", {
    expect_true(all(is.na(signed_linear_fc_to_log2(c(0, NA_real_)))))
    expect_silent(signed_linear_fc_to_log2(c(-3, 0, NA_real_, 3)))
    expect_equal(length(signed_linear_fc_to_log2(numeric(0))), 0L)
})


# ---- the long-form mapping both layers are keyed on -------------------------

test_that("both omics rows of the long-form mapping carry the same merge key", {
    skip_if_not_installed("MultiAssayExperiment")
    skip_if_not_installed("SummarizedExperiment")

    mapping <- data.frame(gene_id = c("EHI_001", "EHI_002"),
                          protein_id = c("XP_1", "XP_2"),
                          stringsAsFactors = FALSE)
    mini_se <- function(ids) {
        m <- matrix(0, nrow = length(ids), ncol = 2,
                    dimnames = list(ids, c("S1", "S2")))
        SummarizedExperiment::SummarizedExperiment(assays = list(counts = m))
    }
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(transcriptomics = mini_se(mapping$gene_id),
                           proteomics = mini_se(mapping$protein_id)))

    legacy <- suppressMessages(.mae_to_legacy(mae, NULL, mapping))
    gm <- legacy$gene_mapping

    rna <- gm[gm$omics == "transcriptomics", ]
    prot <- gm[gm$omics == "proteomics", ]

    # The mapping carries no gene_symbol column, so one is derived -- and the
    # derived value must be the SAME on both sides, or the join on it matches
    # nothing at all.
    expect_setequal(rna$gene_symbol, prot$gene_symbol)
    expect_identical(rna$gene_symbol[match("EHI_001", rna$feature_id)],
                     prot$gene_symbol[match("XP_1", prot$feature_id)])
    # The feature ids stay in their own ID spaces.
    expect_setequal(rna$feature_id, mapping$gene_id)
    expect_setequal(prot$feature_id, mapping$protein_id)
})


# ---- the checkable artifact -------------------------------------------------

test_that("the pairs table stacks every contrast and leads with both ids", {
    fx <- pair_fixture()
    out_dir <- withr::local_tempdir()
    one <- suppressMessages(compute_de_concordance(
        fx$rna_de, fx$prot_de, mapping = fx$mapping, out_dir = NULL))

    f <- suppressMessages(
        .write_log2fc_pairs_table(list(A_vs_B = one, C_vs_D = one), out_dir))

    expect_true(file.exists(f))
    tbl <- read.csv(f, stringsAsFactors = FALSE)
    expect_identical(names(tbl)[1:2], c("gene_id", "protein_id"))
    expect_true(all(c("rna_log2FC", "protein_log2FC", "contrast") %in% names(tbl)))
    expect_setequal(unique(tbl$contrast), c("A_vs_B", "C_vs_D"))
    expect_equal(nrow(tbl), 2L * nrow(one))

    # Recomputing the correlation from the two columns is the whole point of
    # the file: it must reproduce the figure rather than approximate it.
    a <- tbl[tbl$contrast == "A_vs_B", ]
    expect_equal(stats::cor(a$rna_log2FC, a$protein_log2FC),
                 stats::cor(one$rna_log2FC, one$protein_log2FC))
})

test_that("nothing to write leaves no file behind", {
    out_dir <- withr::local_tempdir()

    expect_null(.write_log2fc_pairs_table(list(), out_dir))
    # A contrast table without the two fold-change columns contributes nothing.
    expect_null(.write_log2fc_pairs_table(
        list(A = data.frame(gene_id = "g1", stringsAsFactors = FALSE)), out_dir))
    expect_false(file.exists(file.path(out_dir, "tables",
                                       "rna_protein_log2FC_pairs.csv")))
})

test_that("the MAE fallback label comes from the design, not a literal", {
    expect_identical(.first_contrast_name(list()), "contrast_1")
    expect_identical(
        .first_contrast_name(list(design = list(contrasts = list("S_vs_NS")))),
        "S_vs_NS")
    expect_identical(
        .first_contrast_name(list(design = list(
            contrasts = list(list(name = "Treated_vs_Control"))))),
        "Treated_vs_Control")
})


# ---- the report section switches --------------------------------------------

test_that("the shipped template leaves every report section on", {
    # show_section() treats anything but an explicit FALSE as on, so a template
    # that ships the block empty -- and every config written before it existed --
    # renders exactly as it did.
    cfg <- yaml::read_yaml(file.path(root_dir, "config", "templates",
                                     "multiomics_config.yaml"))
    sections <- (cfg$modes$multiomics$report %||% list())$sections %||% list()

    show_section <- function(name) !identical(sections[[name]], FALSE)
    expect_true(show_section("translation_efficiency"))
    expect_true(show_section("group_distance"))
    expect_true(show_section("a_section_nobody_has_named"))
})
