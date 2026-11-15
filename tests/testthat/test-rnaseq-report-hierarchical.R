# tests/testthat/test-rnaseq-report-hierarchical.R
#
# The RNA report hides the hierarchical clustering heatmap unless
# clustering.steps.hierarchical.show_in_report is true. The flag is display
# only: the heatmap shown is the file the clustering step wrote, so a disabled
# step or a missing file skips the tab instead of failing the render.

rna_report_lines <- function() {
    f <- file.path(root_dir, "R", "domain", "rnaseq", "report_template.Rmd")
    if (!file.exists(f)) stop("R/domain/rnaseq/report_template.Rmd not found")
    readLines(f, warn = FALSE)
}

# Evaluate the template's own gate lines against a config and a Clustering dir.
rna_report_has_hier <- function(hier_cfg, clustering_root) {
    lines <- grep("^(hier_cfg_rna|hier_hm_rna|has_hier_rna) <- ", rna_report_lines(), value = TRUE)
    stopifnot(length(lines) == 3L)
    env <- new.env(parent = globalenv())
    env$rna_cfg <- list(clustering = list(enabled = TRUE, steps = list(hierarchical = hier_cfg)))
    env$clustering_root_rna <- clustering_root
    for (ln in lines) eval(parse(text = ln), envir = env)
    env$has_hier_rna
}

test_that("the RNA report shows the hierarchical heatmap only when asked, enabled and present", {
    root <- withr::local_tempdir()
    hier_dir <- file.path(root, "Hierarchical")
    dir.create(hier_dir)
    file.create(file.path(hier_dir, "Hierarchical_DE_heatmap.png"))

    # Default: flag absent or false keeps the heatmap out of the report.
    expect_false(rna_report_has_hier(list(enabled = TRUE), root))
    expect_false(rna_report_has_hier(list(enabled = TRUE, show_in_report = FALSE), root))
    expect_false(rna_report_has_hier(NULL, root))

    expect_true(rna_report_has_hier(list(enabled = TRUE, show_in_report = TRUE), root))
    # `enabled` missing resolves to TRUE, as in clustering_run_flags().
    expect_true(rna_report_has_hier(list(show_in_report = TRUE), root))

    # Step disabled: a heatmap left by an earlier run is not shown.
    expect_false(rna_report_has_hier(list(enabled = FALSE, show_in_report = TRUE), root))

    # Heatmap missing: skipped, no error.
    empty <- withr::local_tempdir()
    expect_false(rna_report_has_hier(list(enabled = TRUE, show_in_report = TRUE), empty))
})

test_that("the RNA report wires the hierarchical tab to that gate and drops the 'omitted' note", {
    src <- rna_report_lines()
    expect_true(any(grepl("{r clustering-hier-rna, eval=clustering_show_rna && has_hier_rna,", src,
                          fixed = TRUE)))
    expect_true(any(grepl("(has_hier_rna || has_part_rna || has_binary_rna || bc_clust_show_rna)",
                          src, fixed = TRUE)))
    expect_length(grep("Hierarchical clustering is omitted", src, fixed = TRUE), 0L)
})

test_that("the RNA config template ships show_in_report: false", {
    cfg <- yaml::read_yaml(file.path(root_dir, "config", "templates", "rna_config.yaml"))
    hier <- cfg$modes$rna$clustering$steps$hierarchical
    expect_true(hier$enabled)
    expect_false(hier$show_in_report)
})
