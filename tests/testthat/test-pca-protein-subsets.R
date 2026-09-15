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
    touch_pngs(d, c("PCA_top1000.png", "PCA_robust.png", "PCA_all.png",
                    "PCA_top2000.png", "PCA_top500.png", "PCA_PC1.vs.PC2.png"))

    panels <- list_pca_feature_panels(d)
    expect_equal(panels$key, c("all", "top2000", "top1000", "top500", "robust"))
    expect_equal(basename(panels$path[panels$key == "all"]), "PCA_all.png")
    expect_equal(panels$label[panels$key == "top1000"], "Top 1,000 variable proteins")
    expect_match(panels$label[panels$key == "robust"], "measured in every sample")
})

test_that("a run without PCA_all.png falls back to the PC1-vs-PC2 plot", {
    d <- tempfile("pca-panels-")
    dir.create(d)
    on.exit(unlink(d, recursive = TRUE), add = TRUE)
    touch_pngs(d, c("PCA_PC1.vs.PC2.png", "PCA_top500.png"))

    panels <- list_pca_feature_panels(d)
    expect_equal(panels$key, c("all", "top500"))
    expect_equal(basename(panels$path[1]), "PCA_PC1.vs.PC2.png")
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
    expect_true(any(grepl('identical(imp_method, "none")', src, fixed = TRUE)))
    expect_true(any(grepl('"PCA_subset_%s.png"', src, fixed = TRUE)))
    expect_false(any(grepl('sprintf("PCA_%s.png"', src, fixed = TRUE)))
    # all, complete-case and top-N panels, plus earlier subset images
    expect_gte(sum(grepl("file.remove(", src, fixed = TRUE)), 4)
})
