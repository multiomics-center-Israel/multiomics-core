# Tests that mod_proteomics_clustering() cannot leave a previous run's
# hierarchical heatmap behind.
# test-proteomics-clustering-stale-heatmap.R
#
# The module returns early when fewer than two DE features are present in the
# matrix, and it skips the hierarchical step when that step is disabled. Neither
# path used to touch Hierarchical_DE_heatmap.png, so a rerun into the same
# output directory kept the earlier run's image -- which the report displays,
# and which the Methods generator now reads as evidence that clustering ran.
# 01_mod_qc_pre.R already clears PCA_robust.png for exactly this reason.
#
# Scope: the hierarchical heatmap only. Partition and binary-pattern artefacts
# are deliberately not covered here.

csh_meta <- function() {
    data.frame(
        SampleID  = c("S1", "S2", "S3", "S4"),
        Condition = c("ctl", "ctl", "trt", "trt"),
        stringsAsFactors = FALSE
    )
}

csh_expr <- function(meta, features) {
    matrix(
        seq_len(length(features) * nrow(meta)),
        nrow = length(features), byrow = TRUE,
        dimnames = list(features, meta$SampleID)
    )
}

csh_pre <- function(features = c("P1", "P2", "P3")) {
    meta <- csh_meta()
    rownames(meta) <- meta$SampleID
    m <- csh_expr(meta, features)
    list(
        expr_raw        = m,
        expr_filt       = m,
        expr_work       = m,
        expr_imp_single = m,
        meta            = meta,
        row_data        = data.frame(FeatureID = features, stringsAsFactors = FALSE)
    )
}

# pass_any_contrast drives get_de_features(): 1 marks a DE feature.
csh_de <- function(features = c("P1", "P2", "P3"), n_pass = 0L) {
    list(
        method      = "limma",
        imputations = list(),
        summary_df  = data.frame(
            FeatureID         = features,
            pass_any_contrast = c(rep(1L, n_pass), rep(NA_integer_, length(features) - n_pass)),
            stringsAsFactors  = FALSE
        )
    )
}

csh_cfg <- function(hierarchical_enabled = TRUE) {
    list(modes = list(proteomics = list(
        effects    = list(samples = "SampleID", color = "Condition"),
        de_table   = list(id_col = "FeatureID"),
        de         = list(p_cutoff = 0.05, linear_fc_cutoff = 1.5),
        clustering = list(
            enabled = TRUE,
            # clustering_run_flags() calls get_n_groups_from_effects(), which
            # requires this key and aborts without it -- it is not only the
            # binary-patterns gate.
            group_col = "Condition",
            steps     = list(hierarchical = list(enabled = hierarchical_enabled))
        )
    )))
}

# A heatmap left behind by an earlier run into the same output directory.
plant_stale_heatmap <- function(out_dir) {
    hier_dir <- file.path(out_dir, "Clustering", "Hierarchical")
    dir.create(hier_dir, recursive = TRUE, showWarnings = FALSE)
    f <- file.path(hier_dir, "Hierarchical_DE_heatmap.png")
    writeLines("not a real png", f)
    f
}

test_that("a run with too few DE features clears the previous heatmap", {
    out_dir <- withr::local_tempdir()
    stale <- plant_stale_heatmap(out_dir)
    expect_true(file.exists(stale))

    # No feature passes, so get_de_features() returns nothing and the module
    # takes its <2-DE-features early return.
    suppressMessages(
        res <- mod_proteomics_clustering(csh_pre(), csh_de(n_pass = 0L),
                                         csh_cfg(), out_dir)
    )

    expect_false(file.exists(stale))
    expect_length(res$files, 0L)
})

test_that("one DE feature is still too few, and still clears it", {
    out_dir <- withr::local_tempdir()
    stale <- plant_stale_heatmap(out_dir)

    suppressMessages(
        mod_proteomics_clustering(csh_pre(), csh_de(n_pass = 1L), csh_cfg(), out_dir)
    )

    expect_false(file.exists(stale))
})

test_that("a run with hierarchical clustering disabled clears the heatmap", {
    out_dir <- withr::local_tempdir()
    stale <- plant_stale_heatmap(out_dir)

    # Enough DE features to pass the early return, but the hierarchical step is
    # off, so nothing rewrites the image.
    suppressMessages(
        mod_proteomics_clustering(csh_pre(), csh_de(n_pass = 3L),
                                  csh_cfg(hierarchical_enabled = FALSE), out_dir)
    )

    expect_false(file.exists(stale))
})

test_that("the clearing happens before the flags are resolved, not after", {
    # A source-level check, because the ordering is the whole point: the delete
    # must sit between ensure_dir() and clustering_run_flags(), so that every
    # path which can skip rewriting the image has already passed it.
    f <- c(testthat::test_path("..", "..", "R", "modules", "proteomics",
                               "04_mod_clustering.R"),
           "R/modules/proteomics/04_mod_clustering.R")
    f <- f[file.exists(f)][1]
    skip_if(is.na(f), "proteomics clustering module not found")
    lines <- readLines(f, warn = FALSE)

    i_dir   <- grep("ensure_dir\\(clustering_dir\\)", lines)[1]
    i_rm    <- grep("file\\.remove\\(f_hier_hm\\)", lines)[1]
    i_flags <- grep("flags <- clustering_run_flags\\(", lines)[1]

    expect_false(is.na(i_rm))
    expect_gt(i_rm, i_dir)
    expect_lt(i_rm, i_flags)
})
