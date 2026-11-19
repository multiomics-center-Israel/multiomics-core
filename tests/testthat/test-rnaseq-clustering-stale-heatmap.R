# Tests that mod_rnaseq_clustering() cannot leave a previous run's
# hierarchical heatmap behind.
# test-rnaseq-clustering-stale-heatmap.R
#
# The module returns early when there is too little to cluster, and it skips
# the hierarchical step when that step is disabled. Neither path used to touch
# Hierarchical_DE_heatmap.png, so a rerun into the same output directory kept
# the earlier run's image, which the report shows when
# clustering.steps.hierarchical.show_in_report is on. Same guard as proteomics
# (test-proteomics-clustering-stale-heatmap.R).
#
# All data here is synthetic; the fixtures are in helper-rna-batch.R.
# Scope: the hierarchical heatmap only.

rsh_config <- function(hierarchical_enabled = TRUE) {
    cfg <- bc_config()
    cfg$modes$rna$clustering <- list(
        enabled   = TRUE,
        group_col = "Condition",
        steps = list(
            hierarchical    = list(enabled = hierarchical_enabled),
            partition       = list(enabled = FALSE),
            binary_patterns = list(enabled = FALSE)
        )
    )
    cfg
}

# One contrast table in which the first `n_pass` genes pass the DE cutoffs.
rsh_de <- function(n_pass = 0L) {
    ids <- bc_ids()
    pass <- seq_along(ids) <= n_pass
    list(tables = list(Treat_vs_Ctrl = data.frame(
        FeatureID      = ids,
        log2FoldChange = ifelse(pass, 3, 0),
        pvalue         = ifelse(pass, 1e-6, 1),
        padj           = ifelse(pass, 1e-4, 1),
        stringsAsFactors = FALSE
    )))
}

# A heatmap left behind by an earlier run into the same output directory.
rsh_plant_stale_heatmap <- function(out_dir) {
    hier_dir <- file.path(out_dir, "Clustering", "Hierarchical")
    dir.create(hier_dir, recursive = TRUE, showWarnings = FALSE)
    f <- file.path(hier_dir, "Hierarchical_DE_heatmap.png")
    writeLines("not a real png", f)
    f
}

test_that("an RNA run with too few DE genes clears the previous heatmap", {
    for (n_pass in c(0L, 1L)) {
        out_dir <- withr::local_tempdir()
        stale <- rsh_plant_stale_heatmap(out_dir)
        expect_true(file.exists(stale))

        res <- suppressWarnings(suppressMessages(
            mod_rnaseq_clustering(make_bc_pre(), rsh_de(n_pass), rsh_config(), out_dir)
        ))

        expect_false(file.exists(stale), info = sprintf("n_pass = %d", n_pass))
        expect_length(res$files, 0L)
    }
})

test_that("an RNA run with too few features in the matrix clears the heatmap", {
    out_dir <- withr::local_tempdir()
    stale <- rsh_plant_stale_heatmap(out_dir)
    pre <- make_bc_pre()
    pre$expr_work <- pre$expr_work[1:5, , drop = FALSE]

    suppressWarnings(suppressMessages(
        mod_rnaseq_clustering(pre, rsh_de(3L), rsh_config(), out_dir)
    ))

    expect_false(file.exists(stale))
})

test_that("an RNA run with hierarchical clustering disabled clears the heatmap", {
    out_dir <- withr::local_tempdir()
    stale <- rsh_plant_stale_heatmap(out_dir)

    # Enough DE genes to pass the early returns, but the hierarchical step is
    # off (and two groups keep partition/binary off), so nothing rewrites it.
    suppressWarnings(suppressMessages(
        mod_rnaseq_clustering(make_bc_pre(), rsh_de(5L),
                              rsh_config(hierarchical_enabled = FALSE), out_dir)
    ))

    expect_false(file.exists(stale))
})

test_that("the RNA clearing happens before every early return after the enabled check", {
    # A source-level check, because the ordering is the whole point.
    f <- file.path(root_dir, "R", "modules", "rnaseq", "04_mod_clustering.R")
    skip_if_not(file.exists(f), "RNA clustering module not found")
    lines <- readLines(f, warn = FALSE)

    i_rm      <- grep("file\\.remove\\(f_hier_hm\\)", lines)[1]
    i_returns <- grep("return\\(list\\(plots = plots, files = written", lines)
    i_flags   <- grep("flags <- clustering_run_flags\\(", lines)[1]

    expect_false(is.na(i_rm))
    # Only the "clustering disabled" return may come first.
    expect_equal(sum(i_returns < i_rm), 1L)
    expect_lt(i_rm, i_flags)
})
