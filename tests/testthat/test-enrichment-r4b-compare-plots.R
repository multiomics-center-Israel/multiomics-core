# tests/testthat/test-enrichment-r4b-compare-plots.R
#
# §17 R4B — Pipeline comparison plots (NOT Shiny-only). The pipeline generates and
# SAVES, over ALL significant terms, the D9 UP+DOWN diverging horizontal bar plot
# (into up_down_compare/) and the D10 comparison dot matrix (into compare_clusters/
# and compare_patterns/), reusing the existing ORA result tables (no recompute).

# Synthetic processed-ORA table (process_enrichment_table output shape).
.mk_proc_tbl <- function(cluster, count, padj, term = NULL) {
    n <- length(cluster)
    data.frame(
        Cluster = as.character(cluster),
        ID = paste0("GO:", seq_len(n)),
        Description = if (is.null(term)) paste0("term ", seq_len(n)) else term,
        geneID = paste0("g", seq_len(n)),
        Count = as.numeric(count),
        p.adjust = as.numeric(padj),
        in_term = as.numeric(count) * 10,
        Fold_enrichment = round(seq(2, 4, length.out = n), 1),
        stringsAsFactors = FALSE)
}

# Real, definitely-significant ORA result with caller-supplied Cluster labels.
.mk_ora_b <- function(labels) {
    t2g <- data.frame(term = rep(c("T1", "T2", "T3"), each = 20),
                      gene = paste0("g", 1:60), stringsAsFactors = FALSE)
    t2n <- data.frame(term = c("T1", "T2", "T3"),
                      name = c("term one", "term two", "term three"),
                      stringsAsFactors = FALSE)
    run_cluster_ora_compute(labels, TERM2GENE = t2g, TERM2NAME = t2n, type = "KEGG")
}

# ---------------------------------------------------------------------------
# plot functions — write a PDF, fail-soft, no recompute
# ---------------------------------------------------------------------------

test_that("plot_ora_diverging_bars writes a PDF from an up/down table", {
    skip_if_not_installed("ggplot2")
    df <- .mk_proc_tbl(c("up", "up", "down"), c(5, 3, 4), c(1e-3, 1e-2, 2e-2))
    f  <- file.path(withr::local_tempdir(), "diverging_bars.pdf")
    plot_ora_diverging_bars(df, f)
    expect_true(file.exists(f))
})

test_that("plot_ora_diverging_bars is fail-soft: no up/down rows -> no file, no error", {
    skip_if_not_installed("ggplot2")
    df <- .mk_proc_tbl(c("all", "all"), c(5, 3), c(1e-3, 1e-2))  # no up/down
    f  <- file.path(withr::local_tempdir(), "none.pdf")
    expect_null(plot_ora_diverging_bars(df, f))
    expect_false(file.exists(f))
})

test_that("plot_ora_diverging_bars tolerates p.adjust == 0 (no Inf blow-up)", {
    skip_if_not_installed("ggplot2")
    df <- .mk_proc_tbl(c("up", "down"), c(5, 4), c(0, 2e-2))
    f  <- file.path(withr::local_tempdir(), "z.pdf")
    expect_error(plot_ora_diverging_bars(df, f), NA)   # does not error
    expect_true(file.exists(f))
})

test_that("plot_ora_compare_dotmatrix writes a PDF from a multi-group table", {
    skip_if_not_installed("ggplot2")
    df <- .mk_proc_tbl(c("1", "1", "2", "3"), c(5, 3, 4, 2), c(1e-3, 1e-2, 2e-2, 3e-2))
    f  <- file.path(withr::local_tempdir(), "compare_dotmatrix.pdf")
    plot_ora_compare_dotmatrix(df, f)
    expect_true(file.exists(f))
})

# ---------------------------------------------------------------------------
# write_ora_unit_tree wiring — compare leaves get the R4B plots (ALL terms)
# ---------------------------------------------------------------------------

test_that("contrasts compare_up_down gets compare_dotmatrix.pdf over ALL terms (§23 G/H)", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("ggplot2")
    res <- .mk_ora_b(setNames(c(rep("up", 20), rep("down", 20)), paste0("g", 1:40)))
    skip_if(length(res) == 0, "fixture produced no ORA enrichment")

    root <- withr::local_tempdir()
    write_ora_unit_tree(res, ora_root = root, db_name = "KEGG",
                        collection = "contrasts", round = "C",
                        type = "KEGG", dotplot = TRUE)
    cmp <- file.path(root, "KEGG", "contrasts", "C", "compare_up_down")   # §23 H dir
    # §23 G: Compare Up & Down now uses the compare dot matrix (diverging retired)
    expect_true(file.exists(file.path(cmp, "compare_dotmatrix.pdf")))
    expect_false(file.exists(file.path(cmp, "diverging_bars.pdf")))
    # ALL significant terms present in the compare artifact (no Top-N on disk)
    n_all <- nrow(utils::read.csv(file.path(cmp, "results.csv"), colClasses = "character"))
    expect_equal(n_all, nrow(res[[3]]))
})

test_that("partition compare_clusters gets compare_dotmatrix.pdf; per-cluster gets enrichment bars", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("ggplot2")
    res <- .mk_ora_b(setNames(c(rep("1", 20), rep("2", 20)), paste0("g", 1:40)))
    skip_if(length(res) == 0, "fixture produced no ORA enrichment")

    root <- withr::local_tempdir()
    write_ora_unit_tree(res, ora_root = root, db_name = "KEGG",
                        collection = "partition", round = "k",
                        type = "KEGG", dotplot = TRUE)
    base <- file.path(root, "KEGG", "partition")
    expect_true(file.exists(file.path(base, "compare_clusters", "compare_dotmatrix.pdf")))
    expect_false(file.exists(file.path(base, "compare_clusters", "enrichment_bar_top20.pdf")))
    # per-group leaf: canonical summary is the enrichment bar pair (§18 D18.1)
    expect_true(file.exists(file.path(base, "cluster_1", "enrichment_bar_top20.pdf")))
    expect_true(file.exists(file.path(base, "cluster_1", "enrichment_bar_full.pdf")))
})

test_that("binary compare_patterns gets compare_dotmatrix.pdf", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("ggplot2")
    res <- .mk_ora_b(setNames(c(rep("0011", 20), rep("0101", 20)), paste0("g", 1:40)))
    skip_if(length(res) == 0, "fixture produced no ORA enrichment")

    root <- withr::local_tempdir()
    write_ora_unit_tree(res, ora_root = root, db_name = "KEGG",
                        collection = "binary_patterns", round = "best",
                        type = "KEGG", dotplot = TRUE)
    base <- file.path(root, "KEGG", "binary_patterns")
    expect_true(file.exists(file.path(base, "compare_patterns", "compare_dotmatrix.pdf")))
})

test_that("dotplot=FALSE suppresses the R4B comparison plots too", {
    skip_if_not_installed("clusterProfiler")
    res <- .mk_ora_b(setNames(c(rep("up", 20), rep("down", 20)), paste0("g", 1:40)))
    skip_if(length(res) == 0, "fixture produced no ORA enrichment")

    root <- withr::local_tempdir()
    write_ora_unit_tree(res, ora_root = root, db_name = "KEGG",
                        collection = "contrasts", round = "C",
                        type = "KEGG", dotplot = FALSE)
    cmp <- file.path(root, "KEGG", "contrasts", "C", "compare_up_down")
    expect_true(file.exists(file.path(cmp, "results.csv")))       # data still written
    expect_false(file.exists(file.path(cmp, "compare_dotmatrix.pdf")))
})
