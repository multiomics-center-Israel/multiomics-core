test_that("resolve_de_summary_col prefers an exact bare column name", {
  cn <- c("FeatureID", "logFC", "P.Value", "adj.P.Val")
  expect_equal(
    resolve_de_summary_col(cn, bare = c("logFC", "log2FC"),
                           prefixes = c("log2FC"), contrast_label = "a_vs_b"),
    "logFC"
  )
})

test_that("resolve_de_summary_col reads our own contrast-suffixed exports", {
  # Shape of Datasets/deseq2_summary_p0.05.tsv and limma_multimp_summary_p0.05.tsv.
  rna <- c("Gene", "log2FC.S_vs_NS", "linearFC.S_vs_NS",
           "pvalue.S_vs_NS", "padj.S_vs_NS")
  expect_equal(
    resolve_de_summary_col(rna, c("log2FoldChange", "logFC"),
                           c("log2FC", "logFC"), "with_bacteria_vs_No_bacteria"),
    "log2FC.S_vs_NS"
  )

  prot <- c("FeatureID", "linearFC.imputs.SP_vs_NSP",
            "pvalue.imputs.SP_vs_NSP", "padj.imputs.SP_vs_NSP")
  expect_equal(
    resolve_de_summary_col(prot, c("adj.P.Val", "padj"),
                           c("padj.imputs", "padj"), "with_bacteria_vs_No_bacteria"),
    "padj.imputs.SP_vs_NSP"
  )
})

test_that("a single contrast resolves even when its short code differs", {
  # The single-omics runs label the contrast S_vs_NS / SP_vs_NSP while the
  # multiomics run calls it with_bacteria_vs_No_bacteria. One contrast in the
  # file means there is nothing to disambiguate.
  cn <- c("Gene", "log2FC.S_vs_NS")
  expect_equal(
    resolve_de_summary_col(cn, "logFC", "log2FC", "totally_different_label"),
    "log2FC.S_vs_NS"
  )
})

test_that("several contrasts require a matching label and abort otherwise", {
  cn <- c("Gene", "log2FC.A_vs_B", "log2FC.C_vs_D")
  expect_equal(
    resolve_de_summary_col(cn, "logFC", "log2FC", "A_vs_B"),
    "log2FC.A_vs_B"
  )
  expect_error(
    resolve_de_summary_col(cn, "logFC", "log2FC", "Z_vs_Y"),
    "several contrasts"
  )
})

test_that("resolve_de_summary_col returns NA when nothing matches", {
  expect_true(is.na(
    resolve_de_summary_col(c("Gene", "foo"), "logFC", "log2FC", "a_vs_b")
  ))
})

test_that("signed_linear_fc_to_log2 keeps down-regulated features finite", {
  # linearFC is a signed linear ratio: -2 means two-fold LOWER, so log2 = -1.
  # A plain log2() would return NaN and silently drop every such feature.
  x <- c(4, -4, 1, -1, NA, 0)
  out <- signed_linear_fc_to_log2(x)
  expect_equal(out[1], 2)
  expect_equal(out[2], -2)
  expect_equal(out[3], 0)
  expect_equal(out[4], 0)
  expect_true(is.na(out[5]))
  expect_true(is.na(out[6]))
  expect_false(any(is.nan(out[1:4])))
})

test_that("signed_linear_fc_to_log2 emits no warning for negative input", {
  expect_silent(signed_linear_fc_to_log2(c(-2, -8, 3)))
})

test_that("load_precomputed_rna_de aborts rather than returning all-NA stats", {
  # Regression: a table with no fold-change column used to load "successfully",
  # yielding row numbers as feature ids and NA statistics, which read downstream
  # as a run with no differential expression at all.
  dir <- withr::local_tempdir()
  utils::write.table(
    data.frame(Gene = c("g1", "g2"), unrelated = c(1, 2)),
    file.path(dir, "nofc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE
  )
  cfg <- list(
    project = list(dir = dir), paths = list(raw = "."),
    modes = list(rna = list(files = list(de_table = "nofc.tsv")))
  )
  expect_error(
    load_precomputed_rna_de(
      cfg, contrasts_df = data.frame(Contrast_name = "c1", stringsAsFactors = FALSE)
    ),
    "no recognisable log2 fold-change column"
  )
})
