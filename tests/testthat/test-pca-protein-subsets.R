# tests/testthat/test-pca-protein-subsets.R
#
# Tests for the PCA panels written by R/modules/proteomics/01_mod_qc_pre.R and
# listed by the proteomics report:
#   - select_complete_case_features(): the rule behind PCA_robust.png
#   - list_pca_feature_panels(): what the dropdown shows, and in which order
#   - list_pca_subset_panels(): only PCA_subset_<name>.png counts as a subset
# All data are synthetic; the panel files are empty placeholders.

touch_pngs <- function(dir, names) {
    for (n in names) file.create(file.path(dir, n))
    invisible(dir)
}

test_that("complete-case features have no imputed value in any sample", {
    flag <- rbind(
        f0 = c(FALSE, FALSE, FALSE, FALSE),
        f1 = c(TRUE,  FALSE, FALSE, FALSE),
        f2 = c(TRUE,  TRUE,  FALSE, FALSE),
        f4 = c(TRUE,  TRUE,  TRUE,  TRUE)
    )
    expect_equal(rownames(flag)[select_complete_case_features(flag)], "f0")
})

test_that("missingness concentrated in one sample does not pass", {
    # Every protein is missing only in sample 3. Each row has a single imputed
    # value, which a per-protein allowance of one would have kept.
    flag <- matrix(FALSE, nrow = 50, ncol = 4)
    flag[, 3] <- TRUE
    expect_length(select_complete_case_features(flag), 0)
})

test_that("complete-case selection keeps everything when nothing was imputed", {
    expect_length(select_complete_case_features(matrix(FALSE, nrow = 5, ncol = 4)), 5)
})

test_that("NA flags are not counted as imputed", {
    expect_length(select_complete_case_features(rbind(a = c(FALSE, NA, FALSE, FALSE))), 1)
})

test_that("the dropdown lists all proteins first, top-N descending, complete-case last", {
    d <- tempfile("pca-panels-")
    dir.create(d)
    on.exit(unlink(d, recursive = TRUE), add = TRUE)
    touch_pngs(d, c("PCA_top1000.png", "PCA_robust.png",
                    "PCA_top2000.png", "PCA_top500.png", "PCA_PC1.vs.PC2.png"))

    panels <- list_pca_feature_panels(d)
    expect_equal(panels$key, c("all", "top2000", "top1000", "top500", "robust"))
    # "All proteins" is the existing full-matrix PCA, not a second copy of it.
    expect_equal(basename(panels$path[panels$key == "all"]), "PCA_PC1.vs.PC2.png")
    expect_equal(panels$label[panels$key == "top1000"], "Top 1,000 variable proteins")
    # The label states the selection, not an imputation-independence it cannot
    # promise once batch correction has been fitted on the full imputed matrix.
    expect_equal(panels$label[panels$key == "robust"], "Proteins observed in every sample")
    expect_false(grepl("imputed", panels$label[panels$key == "robust"], fixed = TRUE))
})

test_that("a stale PCA_all.png from an older run is never listed", {
    d <- tempfile("pca-panels-")
    dir.create(d)
    on.exit(unlink(d, recursive = TRUE), add = TRUE)
    # PCA_all.png is no longer written; one left behind must not be preferred
    # over this run's full-matrix PCA.
    touch_pngs(d, c("PCA_PC1.vs.PC2.png", "PCA_all.png", "PCA_top500.png"))

    panels <- list_pca_feature_panels(d)
    expect_equal(panels$key, c("all", "top500"))
    expect_equal(basename(panels$path[panels$key == "all"]), "PCA_PC1.vs.PC2.png")
})

test_that("an empty directory lists no panels and no subsets", {
    d <- tempfile("pca-panels-")
    dir.create(d)
    on.exit(unlink(d, recursive = TRUE), add = TRUE)
    expect_equal(nrow(list_pca_feature_panels(d)), 0)
    expect_equal(nrow(list_pca_subset_panels(d)), 0)
})

test_that("only PCA_subset_<name>.png files are listed as sample subsets", {
    d <- tempfile("pca-panels-")
    dir.create(d)
    on.exit(unlink(d, recursive = TRUE), add = TRUE)
    touch_pngs(d, c("PCA_all.png", "PCA_robust.png", "PCA_top500.png",
                    "PCA_PC1.vs.PC2.png", "PCA_PC1.vs.PC3_labeled.png",
                    "PCA_subset_treated.png", "PCA_subset_robust.png"))

    subsets <- list_pca_subset_panels(d)
    expect_equal(subsets$name, c("robust", "treated"))

    # A subset named "robust" is its own file; the dropdown panel is untouched.
    panels <- list_pca_feature_panels(d)
    expect_equal(basename(panels$path[panels$key == "robust"]), "PCA_robust.png")
})

test_that("the QC module uses the shared rule, clears old panels and namespaces subsets", {
    candidates <- c(
        testthat::test_path("..", "..", "R", "modules", "proteomics", "01_mod_qc_pre.R"),
        "R/modules/proteomics/01_mod_qc_pre.R"
    )
    f <- candidates[file.exists(candidates)][1]
    skip_if(is.na(f), "01_mod_qc_pre.R not found from the test working directory")

    src <- readLines(f, warn = FALSE)
    expect_true(any(grepl("select_complete_case_features(", src, fixed = TRUE)))
    expect_true(any(grepl('"PCA_subset_%s.png"', src, fixed = TRUE)))
    expect_false(any(grepl('sprintf("PCA_%s.png"', src, fixed = TRUE)))
    # complete-case and top-N panels, plus the stale PCA_all.png and earlier
    # subset images
    expect_gte(sum(grepl("file.remove(", src, fixed = TRUE)), 4)

    # The panel is no longer skipped on the configured method name. With
    # imputation.method "none" the flags are plain missingness and the regular
    # PCA median-imputes inside compute_pca_scores(), so the complete-case view
    # is if anything more informative there.
    expect_false(any(grepl('identical(imp_method, "none")', src, fixed = TRUE)))

    # The mask is judged on the samples the PCA runs on, and says so loudly if
    # it cannot be aligned.
    expect_true(any(grepl("pca_samples", src, fixed = TRUE)))
    expect_true(any(grepl("missing_flag_cols", src, fixed = TRUE)))

    # "All proteins" reuses the existing full-matrix PCA rather than writing a
    # second copy under another name: the full matrix is projected twice, for
    # PC1-vs-PC2 and PC1-vs-PC3, and not a third time.
    expect_equal(sum(grepl("qc_pca_scatter(pre$expr_imp_single", src, fixed = TRUE)), 2)
})


# =============================================================================
# Non-finite intensities are missing measurements
# =============================================================================

test_that("min-count filtering counts -Inf as an observation, which is why it is normalised upstream", {
    # The hazard the preprocessing fix exists for: pass_filter() asks !is.na(),
    # and -Inf is not NA. A feature observed once and zero twice would look like
    # three observations to the filter.
    expr <- matrix(c(10, -Inf, -Inf,
                     10,   11,   12), nrow = 2, byrow = TRUE,
                   dimnames = list(c("f_zeroes", "f_real"), c("S1", "S2", "S3")))
    grp <- c("A", "A", "A")

    kept_raw <- pass_filter(expr, group = grp, min_per_group = 3, min_groups = 1)
    expect_true(kept_raw[["f_zeroes"]])   # -Inf counted as measured

    expr[!is.finite(expr)] <- NA_real_
    kept_norm <- pass_filter(expr, group = grp, min_per_group = 3, min_groups = 1)
    expect_false(kept_norm[["f_zeroes"]])
    expect_true(kept_norm[["f_real"]])
})

test_that("the mask counts only Inf, -Inf and NaN, not cells that are already NA", {
    # !is.finite() is TRUE for NA as well, so counting with it would report
    # ordinary missingness as values this run converted.
    expr <- matrix(c(1,  NA,  -Inf,
                     Inf, NaN, 6), nrow = 2, byrow = TRUE,
                   dimnames = list(c("f1", "f2"), c("S1", "S2", "S3")))

    nonfinite <- is.infinite(expr) | is.nan(expr)
    expect_equal(sum(nonfinite), 3)          # -Inf, Inf, NaN
    expect_equal(sum(!is.finite(expr)), 4)   # the over-count being avoided
})

test_that("converting non-finite cells changes those cells and nothing else", {
    expr <- matrix(c(1,  NA,  -Inf,
                     Inf, NaN, 6), nrow = 2, byrow = TRUE,
                   dimnames = list(c("f1", "f2"), c("S1", "S2", "S3")))

    nonfinite <- is.infinite(expr) | is.nan(expr)
    expr[nonfinite] <- NA_real_

    # Stated in full rather than probed cell by cell: an earlier version of this
    # test asserted a 2x2 block was all NA while also asserting one cell of that
    # block was 6, and contradicted itself.
    expected <- matrix(c(1,  NA, NA,
                         NA, NA, 6), nrow = 2, byrow = TRUE,
                       dimnames = dimnames(expr))
    expect_equal(expr, expected)
})

test_that("preprocessing normalises non-finite intensities before it filters or imputes", {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "proteomics", "04_preprocess.R"),
        "R/domain/proteomics/04_preprocess.R"
    )
    f <- candidates[file.exists(candidates)][1]
    skip_if(is.na(f), "04_preprocess.R not found from the test working directory")

    # Match the assignments, not the function names: a name also appears in the
    # comments that explain the ordering, and an earlier version of this test
    # matched one of those instead of the call.
    src <- readLines(f, warn = FALSE)
    line_of <- function(pattern) {
        hits <- grep(pattern, src, fixed = TRUE)
        expect_length(hits, 1)
        hits
    }
    norm_line <- line_of("expr_raw[nonfinite] <- NA_real_")
    filt_line <- line_of("filt <- filter_proteomics_by_min_count(")
    imp_line  <- line_of("imp_res <- impute_proteomics(")

    # Order is the point: after this, is.na() means the same thing to the
    # filter, to the imputation flags and to the complete-case selection.
    expect_lt(norm_line, filt_line)
    expect_lt(norm_line, imp_line)
})


# =============================================================================
# The mask and the PCA must share a sample universe
# =============================================================================

test_that("restricting the flags to the PCA samples changes which features qualify", {
    # S4 is in the flag matrix but not in the PCA. Judged on all four samples
    # f_ok would be disqualified by a sample the comparison never shows.
    flag <- rbind(
        f_ok   = c(FALSE, FALSE, FALSE, TRUE),
        f_bad  = c(FALSE, TRUE,  FALSE, FALSE)
    )
    colnames(flag) <- c("S1", "S2", "S3", "S4")

    expect_length(select_complete_case_features(flag), 0)

    pca_samples <- c("S1", "S2", "S3")
    expect_equal(
        rownames(flag)[select_complete_case_features(flag[, pca_samples, drop = FALSE])],
        "f_ok"
    )
})
