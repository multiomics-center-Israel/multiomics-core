# Tests for write_diablo_results() skipping the mixOmics "Y" outcome block.
# block.plsda() one-hot codes the outcome factor and hands it back as an extra
# block named "Y" alongside the real omics. Ten downstream consumers already
# strip it; the writer used to leak it as diablo_scores_Y.csv /
# diablo_top_features_Y.csv, which the report then rendered as an omics layer.
#
# All data below is synthetic.

make_fake_diablo_results <- function() {
    fake_scores <- function(prefix) {
        m <- matrix(
            c(1, 2, 3, 4, 5, 6),
            nrow = 3, ncol = 2,
            dimnames = list(c("S1", "S2", "S3"), c("comp1", "comp2"))
        )
        m + nchar(prefix)  # keep the blocks distinguishable
    }

    fake_top_features <- function(prefix) {
        data.frame(
            feature = paste0(prefix, "_f", 1:2),
            loading = c(0.8, -0.4),
            abs_loading = c(0.8, 0.4),
            component = c("comp1", "comp1"),
            original_name = paste0(prefix, "_name", 1:2),
            stringsAsFactors = FALSE
        )
    }

    list(
        sample_scores = list(
            omicsA = fake_scores("omicsA"),
            omicsB = fake_scores("omicsB"),
            Y = fake_scores("Y")
        ),
        top_features = list(
            omicsA = fake_top_features("omicsA"),
            omicsB = fake_top_features("omicsB"),
            Y = fake_top_features("Y")
        ),
        design = matrix(
            c(0, 1, 1, 0),
            nrow = 2,
            dimnames = list(c("omicsA", "omicsB"), c("omicsA", "omicsB"))
        )
    )
}

test_that("write_diablo_results() does not write the Y outcome block", {
    out_dir <- withr::local_tempdir()

    write_diablo_results(make_fake_diablo_results(), out_dir)

    written <- list.files(out_dir)

    expect_false("diablo_scores_Y.csv" %in% written)
    expect_false("diablo_top_features_Y.csv" %in% written)
    expect_length(grep("_Y\\.csv$", written), 0)
})

test_that("write_diablo_results() still writes every real omics block", {
    out_dir <- withr::local_tempdir()

    write_diablo_results(make_fake_diablo_results(), out_dir)

    written <- list.files(out_dir)

    expect_true(all(c(
        "diablo_scores_omicsA.csv",
        "diablo_scores_omicsB.csv",
        "diablo_top_features_omicsA.csv",
        "diablo_top_features_omicsB.csv",
        "diablo_design_matrix.csv"
    ) %in% written))

    # Content of the real blocks must be untouched by the skip
    scores_a <- read.csv(file.path(out_dir, "diablo_scores_omicsA.csv"),
                         row.names = 1)
    expect_equal(nrow(scores_a), 3)
    expect_equal(ncol(scores_a), 2)

    feats_a <- read.csv(file.path(out_dir, "diablo_top_features_omicsA.csv"),
                        stringsAsFactors = FALSE)
    expect_equal(nrow(feats_a), 2)
    expect_equal(
        colnames(feats_a),
        c("feature", "loading", "abs_loading", "component", "original_name")
    )
})

test_that("write_diablo_results() is unaffected when no Y block is present", {
    res <- make_fake_diablo_results()
    res$sample_scores$Y <- NULL
    res$top_features$Y <- NULL

    out_dir <- withr::local_tempdir()
    write_diablo_results(res, out_dir)

    expect_setequal(
        list.files(out_dir),
        c("diablo_scores_omicsA.csv", "diablo_scores_omicsB.csv",
          "diablo_top_features_omicsA.csv", "diablo_top_features_omicsB.csv",
          "diablo_design_matrix.csv")
    )
})
