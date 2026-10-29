# tests/testthat/test-enrichment-r3-output-tree.R
#
# §17 R3 — Output Directory Redesign (Option B, contrast-first). The physical ORA
# tree becomes ORA/<db>/contrasts/<contrast>/{up,down,compare_up_down,de}/ (§23 H),
# ORA/<db>/partition/{cluster_<x>,compare_clusters}/,
# ORA/<db>/binary_patterns/{pattern_<x>,compare_patterns}/, all_DE/<round>/.
# The per-subgroup artifacts REUSE the same result objects (Cluster filter); no
# ORA is recomputed. The retired with_direction/without_direction/clustering
# splits must never be written.

# A synthetic, definitely-significant ORA result with caller-supplied Cluster
# labels: 3 terms x 20 genes each; each label's genes map cleanly onto a term.
.mk_ora_labeled <- function(labels) {
    t2g <- data.frame(term = rep(c("T1", "T2", "T3"), each = 20),
                      gene = paste0("g", 1:60), stringsAsFactors = FALSE)
    t2n <- data.frame(term = c("T1", "T2", "T3"),
                      name = c("term one", "term two", "term three"),
                      stringsAsFactors = FALSE)
    run_cluster_ora_compute(labels, TERM2GENE = t2g, TERM2NAME = t2n, type = "KEGG")
}

# Read the Cluster column as character — read.csv would otherwise coerce a
# numeric-looking pattern label like "0011" to 11 (leading zeros lost).
.rd_cluster <- function(csv) {
    df <- utils::read.csv(csv, colClasses = "character")
    as.character(df$Cluster)
}

# ora_unit_dir sanitizes a round to a safe path component (dots/colons -> "_").
.sanitize_round <- function(x) gsub("[^a-zA-Z0-9_-]", "_", x)

# ---------------------------------------------------------------------------
# .filter_ora_result_by_cluster — pure row filter, no recompute
# ---------------------------------------------------------------------------

test_that(".filter_ora_result_by_cluster subsets tables + S4 by Cluster, no recompute", {
    skip_if_not_installed("clusterProfiler")
    res <- .mk_ora_labeled(setNames(c(rep("up", 20), rep("down", 20)), paste0("g", 1:40)))
    skip_if(length(res) == 0, "fixture produced no ORA enrichment")

    f <- .filter_ora_result_by_cluster(res, "up")
    # table (res[[3]]) filtered to exactly the original up rows
    orig_up <- res[[3]][as.character(res[[3]]$Cluster) == "up", , drop = FALSE]
    expect_true(all(as.character(f[[3]]$Cluster) == "up"))
    expect_equal(nrow(f[[3]]), nrow(orig_up))
    expect_setequal(f[[3]]$ID, orig_up$ID)          # same IDs => reused, not recomputed
    # S4 compareCluster filtered + factor levels dropped to just "up"
    expect_s4_class(f[[1]], "compareClusterResult")
    expect_true(all(as.character(f[[1]]@compareClusterResult$Cluster) == "up"))
    expect_identical(levels(f[[1]]@compareClusterResult$Cluster), "up")
    # a label with no rows => NULL S4, 0-row table (never an error)
    g <- .filter_ora_result_by_cluster(res, "nonesuch")
    expect_null(g[[1]])
    expect_equal(nrow(g[[3]]), 0L)
})

# ---------------------------------------------------------------------------
# §23 H: contrasts -> up / down / compare_up_down  (disk dirs renamed; internal
# subgroup ids stay up/down/up_down/all). §23 G: compare uses the dot matrix.
# ---------------------------------------------------------------------------

test_that("contrasts unit writes up / down / compare_up_down (§23 H disk dirs)", {
    skip_if_not_installed("clusterProfiler")
    res <- .mk_ora_labeled(setNames(c(rep("up", 20), rep("down", 20)), paste0("g", 1:40)))
    skip_if(length(res) == 0, "fixture produced no ORA enrichment")

    root <- withr::local_tempdir()
    written <- write_ora_unit_tree(res, ora_root = root, db_name = "KEGG",
                                   collection = "contrasts", round = "M_HFD.vs.NC",
                                   type = "KEGG", dotplot = FALSE)
    # the contrast dir component is sanitized (dots -> underscores)
    base <- file.path(root, "KEGG", "contrasts", .sanitize_round("M_HFD.vs.NC"))

    expect_true(file.exists(file.path(base, "up", "results.csv")))
    expect_true(file.exists(file.path(base, "down", "results.csv")))
    expect_true(file.exists(file.path(base, "compare_up_down", "results.csv")))   # §23 H
    # per-direction leaves carry ONLY their direction; compare carries BOTH. The
    # stored Cluster VALUES stay up/down (internal ids unchanged, §23 H).
    expect_true(all(.rd_cluster(file.path(base, "up", "results.csv")) == "up"))
    expect_true(all(.rd_cluster(file.path(base, "down", "results.csv")) == "down"))
    expect_setequal(unique(.rd_cluster(file.path(base, "compare_up_down", "results.csv"))),
                    c("up", "down"))
    # old disk name must NOT be created
    expect_false(dir.exists(file.path(base, "up_down_compare")))
    expect_length(written, 3L)
})

test_that("contrasts_wo_direction writes ONLY the 'de' leaf under the SAME contrast base (§23 H)", {
    skip_if_not_installed("clusterProfiler")
    res <- .mk_ora_labeled(setNames(rep("all", 40), paste0("g", 1:40)))
    skip_if(length(res) == 0, "fixture produced no ORA enrichment")

    root <- withr::local_tempdir()
    write_ora_unit_tree(res, ora_root = root, db_name = "KEGG",
                        collection = "contrasts_wo_direction", round = "M_HFD.vs.NC",
                        type = "KEGG", dotplot = FALSE)
    base <- file.path(root, "KEGG", "contrasts", .sanitize_round("M_HFD.vs.NC"))

    # §23 H: the direction-agnostic view is written to the `de/` dir; the internal
    # Cluster value stays "all" (id unchanged, only the disk dir renamed).
    expect_true(file.exists(file.path(base, "de", "results.csv")))
    expect_true(all(.rd_cluster(file.path(base, "de", "results.csv")) == "all"))
    expect_false(dir.exists(file.path(base, "all")))            # old disk name gone
    expect_false(dir.exists(file.path(base, "up")))
    expect_false(dir.exists(file.path(base, "compare_up_down")))
})

# ---------------------------------------------------------------------------
# partition / binary_patterns -> per-label + compare
# ---------------------------------------------------------------------------

test_that("partition writes cluster_<x> per cluster + compare_clusters (>=2 clusters)", {
    skip_if_not_installed("clusterProfiler")
    res <- .mk_ora_labeled(setNames(c(rep("1", 20), rep("2", 20)), paste0("g", 1:40)))
    skip_if(length(res) == 0, "fixture produced no ORA enrichment")

    root <- withr::local_tempdir()
    write_ora_unit_tree(res, ora_root = root, db_name = "KEGG",
                        collection = "partition", round = "k",
                        type = "KEGG", dotplot = FALSE)
    base <- file.path(root, "KEGG", "partition")

    expect_true(file.exists(file.path(base, "cluster_1", "results.csv")))
    expect_true(file.exists(file.path(base, "cluster_2", "results.csv")))
    expect_true(file.exists(file.path(base, "compare_clusters", "results.csv")))
    expect_true(all(.rd_cluster(file.path(base, "cluster_1", "results.csv")) == "1"))
    expect_setequal(unique(.rd_cluster(file.path(base, "compare_clusters", "results.csv"))),
                    c("1", "2"))
    # retired clustering/ wrapper gone
    expect_false(dir.exists(file.path(root, "KEGG", "clustering")))
})

test_that("partition with a single cluster writes NO compare_clusters", {
    skip_if_not_installed("clusterProfiler")
    res <- .mk_ora_labeled(setNames(rep("1", 20), paste0("g", 1:20)))
    skip_if(length(res) == 0, "fixture produced no ORA enrichment")

    root <- withr::local_tempdir()
    write_ora_unit_tree(res, ora_root = root, db_name = "KEGG",
                        collection = "partition", round = "k",
                        type = "KEGG", dotplot = FALSE)
    base <- file.path(root, "KEGG", "partition")
    expect_true(file.exists(file.path(base, "cluster_1", "results.csv")))
    expect_false(dir.exists(file.path(base, "compare_clusters")))
})

test_that("binary_patterns writes pattern_<x> per pattern + compare_patterns", {
    skip_if_not_installed("clusterProfiler")
    res <- .mk_ora_labeled(setNames(c(rep("0011", 20), rep("0101", 20)), paste0("g", 1:40)))
    skip_if(length(res) == 0, "fixture produced no ORA enrichment")

    root <- withr::local_tempdir()
    write_ora_unit_tree(res, ora_root = root, db_name = "KEGG",
                        collection = "binary_patterns", round = "best",
                        type = "KEGG", dotplot = FALSE)
    base <- file.path(root, "KEGG", "binary_patterns")

    expect_true(file.exists(file.path(base, "pattern_0011", "results.csv")))
    expect_true(file.exists(file.path(base, "pattern_0101", "results.csv")))
    expect_true(file.exists(file.path(base, "compare_patterns", "results.csv")))
    expect_true(all(.rd_cluster(file.path(base, "pattern_0011", "results.csv")) == "0011"))
})

# ---------------------------------------------------------------------------
# all_DE -> single whole-table leaf at the base
# ---------------------------------------------------------------------------

test_that("all_DE writes a single whole-table leaf at all_DE/<round> (no subgroup dirs)", {
    skip_if_not_installed("clusterProfiler")
    res <- .mk_ora_labeled(setNames(rep("all_genes", 40), paste0("g", 1:40)))
    skip_if(length(res) == 0, "fixture produced no ORA enrichment")

    root <- withr::local_tempdir()
    written <- write_ora_unit_tree(res, ora_root = root, db_name = "KEGG",
                                   collection = "all_DE", round = "any_contrast",
                                   type = "KEGG", dotplot = FALSE)
    base <- file.path(root, "KEGG", "all_DE", "any_contrast")
    expect_true(file.exists(file.path(base, "results.csv")))
    expect_length(written, 1L)
    expect_equal(normalizePath(written[1]), normalizePath(base))
})

# ---------------------------------------------------------------------------
# plots land in the new leaves (enrichment bars per subgroup, config-gated)
# ---------------------------------------------------------------------------

test_that("plots (when enabled) land in each leaf: per-group bars; compare_up_down gets the dot matrix (§23 G/H)", {
    skip_if_not_installed("clusterProfiler")
    skip_if_not_installed("ggplot2")
    res <- .mk_ora_labeled(setNames(c(rep("up", 20), rep("down", 20)), paste0("g", 1:40)))
    skip_if(length(res) == 0, "fixture produced no ORA enrichment")

    root <- withr::local_tempdir()
    write_ora_unit_tree(res, ora_root = root, db_name = "KEGG",
                        collection = "contrasts", round = "C",
                        type = "KEGG", dotplot = TRUE)
    base <- file.path(root, "KEGG", "contrasts", "C")
    # per-group leaf: canonical summary is the enrichment bar pair (§18/§23 B)
    expect_true(file.exists(file.path(base, "up", "enrichment_bar_top20.pdf")))
    expect_true(file.exists(file.path(base, "up", "enrichment_bar_full.pdf")))
    # §23 G/H: compare_up_down gets the compare DOT MATRIX (diverging retired)
    expect_true(file.exists(file.path(base, "compare_up_down", "compare_dotmatrix.pdf")))
    expect_false(file.exists(file.path(base, "compare_up_down", "diverging_bars.pdf")))
    expect_false(file.exists(file.path(base, "compare_up_down", "enrichment_bar_top20.pdf")))
})
