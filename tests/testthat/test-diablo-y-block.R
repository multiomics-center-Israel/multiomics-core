# tests/testthat/test-diablo-y-block.R
#
# write_diablo_results() must not write the mixOmics "Y" outcome block.
# block.plsda() one-hot codes the outcome factor and hands it back as an extra
# block named "Y" alongside the real omics. Seven places downstream already
# strip it (04_integration_diablo.R:385/632/1146, 07_enrichment.R:1830,
# 10_integration_consensus.R:341/702/890), but the writer leaked it as
# diablo_scores_Y.csv / diablo_top_features_Y.csv and the report rendered it as
# an omics layer.
#
# The writer also retires those two paths on the way in, so a rerun into a
# directory from before the fix does not leave them behind for the report's
# diablo_top_features_*.csv glob to find as stale results.
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

test_that("stale Y files from an earlier run are retired", {
    # Not writing the Y block is not enough on a rerun: the report globs
    # diablo_top_features_*.csv, so a pair left by a run from before this change
    # would keep surfacing the outcome as an omics layer -- and by then as a
    # stale result that no longer tracks the data.
    out_dir <- withr::local_tempdir()
    writeLines("stale,content", file.path(out_dir, "diablo_scores_Y.csv"))
    writeLines("stale,content", file.path(out_dir, "diablo_top_features_Y.csv"))

    write_diablo_results(make_fake_diablo_results(), out_dir)

    expect_false(file.exists(file.path(out_dir, "diablo_scores_Y.csv")))
    expect_false(file.exists(file.path(out_dir, "diablo_top_features_Y.csv")))

    # The real blocks are written as usual, and nothing else was swept up.
    expect_setequal(
        list.files(out_dir),
        c("diablo_scores_omicsA.csv", "diablo_scores_omicsB.csv",
          "diablo_top_features_omicsA.csv", "diablo_top_features_omicsB.csv",
          "diablo_design_matrix.csv")
    )
})

test_that("only the two Y artifacts are retired, not other files in the directory", {
    # The removal is by exact path, so anything else a results directory holds --
    # including an unrelated file whose name contains Y -- is left alone.
    out_dir <- withr::local_tempdir()
    writeLines("stale", file.path(out_dir, "diablo_scores_Y.csv"))
    writeLines("keep",  file.path(out_dir, "diablo_scores_Yeast.csv"))
    writeLines("keep",  file.path(out_dir, "notes_Y.txt"))

    write_diablo_results(make_fake_diablo_results(), out_dir)

    expect_false(file.exists(file.path(out_dir, "diablo_scores_Y.csv")))
    expect_true(file.exists(file.path(out_dir, "diablo_scores_Yeast.csv")))
    expect_true(file.exists(file.path(out_dir, "notes_Y.txt")))
})

test_that("a block genuinely named Y-something is not caught by the skip", {
    # setdiff() matches the whole name, so only the exact "Y" block is dropped.
    res <- make_fake_diablo_results()
    res$sample_scores$Yeast <- res$sample_scores$omicsA
    res$top_features$Yeast <- res$top_features$omicsA

    out_dir <- withr::local_tempdir()
    write_diablo_results(res, out_dir)

    written <- list.files(out_dir)
    expect_true("diablo_scores_Yeast.csv" %in% written)
    expect_false("diablo_scores_Y.csv" %in% written)
})
