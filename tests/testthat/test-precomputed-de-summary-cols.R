# Reading back the summary tables this pipeline itself writes.
#
# The pre-computed DE branch existed in both loaders, but neither could read our
# own Datasets/*_summary_p0.05.tsv exports: those hold every contrast in one
# table and suffix the statistics with the contrast name (log2FC.S_vs_NS,
# padj.imputs.SP_vs_NSP), while the loaders matched bare names only. Pointing a
# config at one of them "succeeded", logged a plausible feature count, and
# returned every statistic as NA -- which downstream is indistinguishable from a
# run in which nothing was differentially expressed. That is the failure these
# tests exist to prevent: not a crash, a silent nothing.
#
# All fixtures synthetic; no real export is read.

# ---- which column the resolver picks ----------------------------------------

test_that("an exact bare name wins over any suffixed column", {
    cn <- c("Gene", "log2FoldChange", "log2FC.A_vs_B")

    expect_identical(
        resolve_de_summary_col(cn, bare = "log2FoldChange",
                               prefixes = "log2FC", contrast_label = "A_vs_B"),
        "log2FoldChange")
})

test_that("bare candidates are ranked by the caller, not by the file", {
    # The caller lists its candidates in preference order. Subsetting the column
    # names instead -- cn[cn %in% bare] -- silently hands that decision to
    # whichever candidate the table happens to print first, so a table written
    # with logFC before log2FoldChange would flip a caller's stated preference.
    cn <- c("Gene", "logFC", "log2FoldChange")

    expect_identical(
        resolve_de_summary_col(cn, bare = c("log2FoldChange", "logFC"),
                               prefixes = character(0)),
        "log2FoldChange")
    expect_identical(
        resolve_de_summary_col(cn, bare = c("logFC", "log2FoldChange"),
                               prefixes = character(0)),
        "logFC")
})

test_that("our own contrast-suffixed exports resolve", {
    cn <- c("Gene", "log2FC.S_vs_NS", "log2FC.SP_vs_NSP",
            "padj.S_vs_NS", "padj.SP_vs_NSP")

    expect_identical(
        resolve_de_summary_col(cn, bare = "log2FoldChange",
                               prefixes = "log2FC", contrast_label = "SP_vs_NSP"),
        "log2FC.SP_vs_NSP")
    expect_identical(
        resolve_de_summary_col(cn, bare = "padj", prefixes = "padj",
                               contrast_label = "S_vs_NS"),
        "padj.S_vs_NS")
})

test_that("a lone contrast resolves even when its short code differs", {
    # This is what lets a single-omics export labelled S_vs_NS serve a
    # multiomics run whose contrast is called with_bacteria_vs_No_bacteria.
    cn <- c("Gene", "log2FC.S_vs_NS", "padj.S_vs_NS")

    expect_identical(
        resolve_de_summary_col(cn, bare = "log2FoldChange", prefixes = "log2FC",
                               contrast_label = "with_bacteria_vs_No_bacteria"),
        "log2FC.S_vs_NS")
})

test_that("an exact contrast match is preferred over a lone column under an earlier prefix", {
    # Checked prefix by prefix, "logFC" here holds exactly one column and would
    # be taken as the lone-contrast case -- never reaching the exact match for
    # the requested contrast sitting under "log2FC". The exact sweep runs across
    # every prefix first so the lone-column fallback means what it says.
    cn <- c("Gene", "logFC.OTHER", "log2FC.A_vs_B", "log2FC.C_vs_D")

    expect_identical(
        resolve_de_summary_col(cn, bare = "log2FoldChange",
                               prefixes = c("logFC", "log2FC"),
                               contrast_label = "A_vs_B"),
        "log2FC.A_vs_B")
})

test_that("a lone column under one prefix is not a lone contrast in the table", {
    # Three contrasts are represented here, spread across two prefixes. Deciding
    # prefix by prefix, "logFC" holds exactly one column and would be taken as
    # the single-contrast case -- returning OTHER's statistics under a label
    # that is not OTHER. The count has to be over every prefix.
    cn <- c("Gene", "logFC.OTHER", "log2FC.A_vs_B", "log2FC.C_vs_D")

    expect_error(
        resolve_de_summary_col(cn, bare = "log2FoldChange",
                               prefixes = c("logFC", "log2FC"),
                               contrast_label = "E_vs_F"),
        "OTHER")
    expect_true(is.na(
        resolve_de_summary_col(cn, bare = "log2FoldChange",
                               prefixes = c("logFC", "log2FC"))))
})

test_that("an overlapping prefix does not make one contrast look like two", {
    # "pvalue" is a prefix of "pvalue.imputs", so a single proteomics column can
    # be read twice -- as contrast S_vs_NS under the longer prefix and as
    # "imputs.S_vs_NS" under the shorter. That would count two contrasts where
    # the table holds one, and abort on a table that resolves perfectly well.
    cn <- c("Protein.Group", "pvalue.imputs.S_vs_NS", "padj.imputs.S_vs_NS")

    expect_identical(
        resolve_de_summary_col(cn, bare = "P.Value",
                               prefixes = c("pvalue.imputs", "P.Value", "pvalue"),
                               contrast_label = "some_other_label"),
        "pvalue.imputs.S_vs_NS")
})

test_that("several contrasts and no matching label aborts, naming what it found", {
    # Guessing here would attach one contrast's statistics to another's label.
    cn <- c("Gene", "log2FC.A_vs_B", "log2FC.C_vs_D")

    expect_error(
        resolve_de_summary_col(cn, bare = "log2FoldChange", prefixes = "log2FC",
                               contrast_label = "E_vs_F"),
        "A_vs_B, C_vs_D")
})

test_that("nothing matching resolves to NA rather than an error", {
    # The caller decides whether a missing column is fatal: a missing fold
    # change is, a missing padj is not -- it gets derived from the raw p.
    expect_true(is.na(
        resolve_de_summary_col(c("Gene", "unrelated"), bare = "padj",
                               prefixes = "padj", contrast_label = "A_vs_B")))
})


# ---- what the loaders do with an unreadable table ---------------------------

write_tsv_fixture <- function(df, dir, name) {
    path <- file.path(dir, name)
    utils::write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
    path
}

rna_cfg <- function(dir, file) {
    list(project = list(dir = dir), paths = list(raw = "."),
         modes = list(rna = list(files = list(de_table = file))))
}

test_that("the RNA loader aborts rather than returning all-NA statistics", {
    # Regression: this table used to load "successfully", yielding row numbers
    # as feature ids and NA for every statistic.
    dir <- withr::local_tempdir()
    write_tsv_fixture(data.frame(Gene = c("g1", "g2"), unrelated = c(1, 2)),
                      dir, "nofc.tsv")

    expect_error(
        suppressMessages(load_precomputed_rna_de(
            rna_cfg(dir, "nofc.tsv"),
            contrasts_df = data.frame(Contrast_name = "c1",
                                      stringsAsFactors = FALSE))),
        "no recognisable log2 fold-change column")
})

test_that("the RNA loader reads a suffixed export, ids and statistics alike", {
    dir <- withr::local_tempdir()
    write_tsv_fixture(
        data.frame(Gene = c("g1", "g2"),
                   log2FC.S_vs_NS = c(1.5, -2),
                   pvalue.S_vs_NS = c(0.01, 0.2),
                   padj.S_vs_NS = c(0.03, 0.4)),
        dir, "summary.tsv")

    res <- suppressMessages(load_precomputed_rna_de(
        rna_cfg(dir, "summary.tsv"),
        contrasts_df = data.frame(Contrast_name = "S_vs_NS",
                                  stringsAsFactors = FALSE)))
    tab <- res$tables[["S_vs_NS"]]

    # "Gene" is the id column our own RNA export writes; without it the loader
    # fell through to row names and the ids came back as "1", "2".
    expect_identical(as.character(tab$FeatureID), c("g1", "g2"))
    expect_equal(tab$log2FoldChange, c(1.5, -2))
    expect_false(any(is.na(tab$padj)))
})

test_that("the proteomics loader converts a linearFC-only summary", {
    # The multi-imputation export carries no logFC at all. linearFC is a signed
    # linear ratio, so log2() of it returns NaN for every down-regulated
    # protein -- about half of them, silently dropped.
    dir <- withr::local_tempdir()
    write_tsv_fixture(
        data.frame(Protein.Group = c("p1", "p2"),
                   linearFC.imputs.S_vs_NS = c(4, -4),
                   pvalue.imputs.S_vs_NS = c(0.01, 0.02),
                   padj.imputs.S_vs_NS = c(0.03, 0.04)),
        dir, "prot_summary.tsv")

    cfg <- list(project = list(dir = dir), paths = list(raw = "."),
                modes = list(proteomics = list(
                    files = list(de_table = "prot_summary.tsv"))))

    res <- suppressMessages(load_precomputed_proteomics_de(
        cfg, contrasts_df = data.frame(Contrast_name = "S_vs_NS",
                                       stringsAsFactors = FALSE)))
    tab <- res$runs_de_tables[[1]][["S_vs_NS"]]

    # 4 and -4 are a four-fold rise and a four-fold fall: +2 and -2 in log2.
    # log2(-4) would have been NaN.
    expect_equal(tab$logFC, c(2, -2))
    expect_false(any(is.nan(tab$logFC)))
})

# ---- one wide file holding several contrasts --------------------------------

test_that("the RNA loader splits one wide summary into a table per contrast", {
    # The shape our own export actually has. Pairing labels to files by
    # position fell through to the filename here -- one file, several contrasts
    # -- and the filename matches no contrast in the table, so the loader
    # aborted on precisely the round-trip it was taught to do.
    dir <- withr::local_tempdir()
    write_tsv_fixture(
        data.frame(Gene = c("g1", "g2"),
                   log2FC.A_vs_B = c(1.5, -2),
                   pvalue.A_vs_B = c(0.01, 0.2),
                   padj.A_vs_B   = c(0.03, 0.4),
                   log2FC.C_vs_D = c(0.5, -0.25),
                   pvalue.C_vs_D = c(0.04, 0.5),
                   padj.C_vs_D   = c(0.06, 0.7)),
        dir, "deseq2_summary_p0.05.tsv")

    res <- suppressMessages(load_precomputed_rna_de(
        rna_cfg(dir, "deseq2_summary_p0.05.tsv"),
        contrasts_df = data.frame(Contrast_name = c("A_vs_B", "C_vs_D"),
                                  stringsAsFactors = FALSE)))

    expect_named(res$tables, c("A_vs_B", "C_vs_D"))
    # Each table carries its own contrast's statistics, not the other's.
    expect_equal(res$tables[["A_vs_B"]]$log2FoldChange, c(1.5, -2))
    expect_equal(res$tables[["C_vs_D"]]$log2FoldChange, c(0.5, -0.25))
    expect_equal(res$tables[["A_vs_B"]]$padj, c(0.03, 0.4))
    expect_equal(res$tables[["C_vs_D"]]$padj, c(0.06, 0.7))
    expect_identical(as.character(res$tables[["A_vs_B"]]$FeatureID),
                     c("g1", "g2"))
})

test_that("the proteomics loader splits a wide limma_multimp summary", {
    # log2FC.imputs is what that export writes. It was missing from the
    # fold-change prefixes, so the log2FC. stem claimed the column and read its
    # contrast as "imputs.A_vs_B" -- survivable for one contrast, an abort for
    # two, and it preferred the signif()-rounded linearFC beside it.
    dir <- withr::local_tempdir()
    write_tsv_fixture(
        data.frame(FeatureID = c("p1", "p2"),
                   log2FC.imputs.A_vs_B   = c(2, -2),
                   pvalue.imputs.A_vs_B   = c(0.01, 0.2),
                   padj.imputs.A_vs_B     = c(0.03, 0.4),
                   log2FC.imputs.C_vs_D   = c(0.5, -0.25),
                   pvalue.imputs.C_vs_D   = c(0.04, 0.5),
                   padj.imputs.C_vs_D     = c(0.06, 0.7)),
        dir, "limma_multimp_summary.tsv")

    cfg <- list(project = list(dir = dir), paths = list(raw = "."),
                modes = list(proteomics = list(
                    files = list(de_table = "limma_multimp_summary.tsv"))))

    res <- suppressMessages(load_precomputed_proteomics_de(
        cfg, contrasts_df = data.frame(Contrast_name = c("A_vs_B", "C_vs_D"),
                                       stringsAsFactors = FALSE)))
    tabs <- res$runs_de_tables[[1]]

    expect_named(tabs, c("A_vs_B", "C_vs_D"))
    expect_equal(tabs[["A_vs_B"]]$logFC, c(2, -2))
    expect_equal(tabs[["C_vs_D"]]$logFC, c(0.5, -0.25))
    expect_equal(tabs[["A_vs_B"]]$adj.P.Val, c(0.03, 0.4))
    expect_equal(tabs[["C_vs_D"]]$adj.P.Val, c(0.06, 0.7))
})

test_that("a wide summary with no contrasts table is refused, not guessed at", {
    # mod_rnaseq_de() returns into the pre-computed branch before the
    # auto_generate_contrasts() fallback, and files$contrasts is documented
    # optional -- so contrasts_df can genuinely be NULL here. A wide file names
    # its contrasts only in column suffixes, so nothing is left to say which
    # was wanted. The error has to name the remedy: the one raised further down
    # blames the contrast naming and sends the reader to fix the wrong thing.
    dir <- withr::local_tempdir()
    write_tsv_fixture(
        data.frame(Gene = c("g1", "g2"),
                   log2FC.A_vs_B = c(1.5, -2),
                   pvalue.A_vs_B = c(0.01, 0.2),
                   log2FC.C_vs_D = c(0.5, -0.25),
                   pvalue.C_vs_D = c(0.04, 0.5)),
        dir, "deseq2_summary_p0.05.tsv")

    expect_error(
        suppressMessages(load_precomputed_rna_de(
            rna_cfg(dir, "deseq2_summary_p0.05.tsv"), contrasts_df = NULL)),
        "modes\\.rna\\.files\\.contrasts")
})

test_that("a wide summary holding one contrast still loads without a contrasts table", {
    # The refusal is about ambiguity, not about wide files. One contrast leaves
    # nothing to choose between, and that case resolved before this check
    # existed -- it must keep doing so.
    dir <- withr::local_tempdir()
    write_tsv_fixture(
        data.frame(Gene = c("g1", "g2"),
                   log2FC.S_vs_NS = c(1.5, -2),
                   pvalue.S_vs_NS = c(0.01, 0.2)),
        dir, "deseq2_summary_p0.05.tsv")

    res <- suppressMessages(load_precomputed_rna_de(
        rna_cfg(dir, "deseq2_summary_p0.05.tsv"), contrasts_df = NULL))

    expect_length(res$tables, 1L)
    expect_equal(res$tables[[1]]$log2FoldChange, c(1.5, -2))
})

test_that("one file per contrast still pairs by position", {
    # The other supported shape, unchanged: two files, two contrasts, paired in
    # order rather than split out of one table.
    dir <- withr::local_tempdir()
    write_tsv_fixture(data.frame(Gene = c("g1", "g2"),
                                 log2FoldChange = c(1.5, -2),
                                 pvalue = c(0.01, 0.2)),
                      dir, "de_A_vs_B.tsv")
    write_tsv_fixture(data.frame(Gene = c("g1", "g2"),
                                 log2FoldChange = c(0.5, -0.25),
                                 pvalue = c(0.04, 0.5)),
                      dir, "de_C_vs_D.tsv")

    cfg <- list(project = list(dir = dir), paths = list(raw = "."),
                modes = list(rna = list(files = list(
                    de_table = list("de_A_vs_B.tsv", "de_C_vs_D.tsv")))))

    res <- suppressMessages(load_precomputed_rna_de(
        cfg, contrasts_df = data.frame(Contrast_name = c("A_vs_B", "C_vs_D"),
                                       stringsAsFactors = FALSE)))

    expect_named(res$tables, c("A_vs_B", "C_vs_D"))
    expect_equal(res$tables[["A_vs_B"]]$log2FoldChange, c(1.5, -2))
    expect_equal(res$tables[["C_vs_D"]]$log2FoldChange, c(0.5, -0.25))
})

test_that("the proteomics loader aborts when no fold change can be resolved", {
    dir <- withr::local_tempdir()
    write_tsv_fixture(data.frame(Protein.Group = c("p1", "p2"),
                                 unrelated = c(1, 2)),
                      dir, "prot_nofc.tsv")

    cfg <- list(project = list(dir = dir), paths = list(raw = "."),
                modes = list(proteomics = list(
                    files = list(de_table = "prot_nofc.tsv"))))

    expect_error(
        suppressMessages(load_precomputed_proteomics_de(
            cfg, contrasts_df = data.frame(Contrast_name = "S_vs_NS",
                                           stringsAsFactors = FALSE))),
        "no recognisable fold-change column")
})
