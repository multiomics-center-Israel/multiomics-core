# tests/testthat/test-build-de-row-annotations.R
#
# Unit tests for build_de_row_annotations() — verifying all three
# column style detection branches (proteomics, RNA-seq, metabolomics).

# Source if not already loaded
if (!exists("build_de_row_annotations")) {
    source(file.path("../../R/core/09_clustering.R"))
}

# ---- Test data builders ----

make_summary_proteomics <- function() {
    data.frame(
        FeatureID             = paste0("P", 1:5),
        linearFC.imputs.TvsC  = c(2.5, -1.8, 0.3, 3.0, -2.1),
        padj.imputs.TvsC      = c(0.01, 0.03, 0.80, 0.001, 0.04),
        pass.imputs.TvsC      = c(1, 1, 0, 1, 1),
        pass_any_contrast     = c(1, 1, 0, 1, 1),
        stringsAsFactors      = FALSE
    )
}

make_summary_rnaseq <- function() {
    data.frame(
        FeatureID         = paste0("G", 1:5),
        linearFC.TvsC     = c(2.5, -1.8, 0.3, 3.0, -2.1),
        padj.TvsC         = c(0.01, 0.03, 0.80, 0.001, 0.04),
        TvsC_pass         = c("up", "down", "", "up", "down"),
        pass_any_contrast = c(1, 1, 0, 1, 1),
        stringsAsFactors  = FALSE
    )
}

make_summary_metabolomics <- function() {
    data.frame(
        feature_id        = paste0("M", 1:5),
        linearFC.TvsC     = c(2.5, -1.8, 0.3, 3.0, -2.1),
        padj.TvsC         = c(0.01, 0.03, 0.80, 0.001, 0.04),
        pass.TvsC         = c(1, 1, 0, 1, 1),
        pass_any_contrast = c(1, 1, 0, 1, 1),
        stringsAsFactors  = FALSE
    )
}


# ---- Proteomics branch ----

test_that("build_de_row_annotations detects proteomics columns", {
    df <- make_summary_proteomics()
    res <- build_de_row_annotations(df, feature_ids = paste0("P", 1:5),
                                    id_col = "FeatureID")
    expect_false(is.null(res))
    expect_true("TvsC" %in% colnames(res))
    expect_equal(res["P1", "TvsC"], "up")
    expect_equal(res["P2", "TvsC"], "down")
    expect_true(is.na(res["P3", "TvsC"]))
})


# ---- RNA-seq branch ----

test_that("build_de_row_annotations detects RNA-seq columns", {
    df <- make_summary_rnaseq()
    res <- build_de_row_annotations(df, feature_ids = paste0("G", 1:5),
                                    id_col = "FeatureID")
    expect_false(is.null(res))
    expect_true("TvsC" %in% colnames(res))
    expect_equal(res["G1", "TvsC"], "up")
    expect_equal(res["G2", "TvsC"], "down")
})


# ---- Metabolomics branch ----

test_that("build_de_row_annotations detects metabolomics columns (pass.<cn>)", {
    df <- make_summary_metabolomics()
    res <- build_de_row_annotations(df, feature_ids = paste0("M", 1:5),
                                    id_col = "feature_id")
    expect_false(is.null(res))
    expect_true("TvsC" %in% colnames(res))
    expect_equal(res["M1", "TvsC"], "up")
    expect_equal(res["M2", "TvsC"], "down")
    expect_true(is.na(res["M3", "TvsC"]))
    expect_equal(res["M4", "TvsC"], "up")
    expect_equal(res["M5", "TvsC"], "down")
})


# ---- No pass columns → NULL with warning ----

test_that("build_de_row_annotations returns NULL when no pass columns found", {
    df <- data.frame(
        feature_id    = paste0("X", 1:3),
        logFC.TvsC    = c(1, -1, 0.5),
        stringsAsFactors = FALSE
    )
    expect_warning(
        res <- build_de_row_annotations(df, feature_ids = paste0("X", 1:3),
                                        id_col = "feature_id"),
        "No pass columns"
    )
    expect_null(res)
})


# ---- Multi-contrast metabolomics ----

test_that("build_de_row_annotations handles multiple metabolomics contrasts", {
    df <- data.frame(
        feature_id    = paste0("M", 1:4),
        linearFC.AvsB = c(2.0, -1.5, 0.1, 3.0),
        pass.AvsB     = c(1, 1, 0, 1),
        linearFC.CvsD = c(-2.0, 0.5, 1.8, -0.2),
        pass.CvsD     = c(1, 0, 1, 0),
        pass_any_contrast = c(1, 1, 1, 1),
        stringsAsFactors  = FALSE
    )
    res <- build_de_row_annotations(df, feature_ids = paste0("M", 1:4),
                                    id_col = "feature_id")
    expect_false(is.null(res))
    expect_true(all(c("AvsB", "CvsD") %in% colnames(res)))
    expect_equal(res["M1", "AvsB"], "up")
    expect_equal(res["M1", "CvsD"], "down")
    expect_equal(res["M3", "CvsD"], "up")
})
