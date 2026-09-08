# tests/testthat/test-pca-protein-subsets.R
#
# Tests for the PCA protein-subset panels (R/modules/proteomics/01_mod_qc_pre.R
# and the report template's selector). The "robust" panel is the one with a
# real rule behind it: keep only features that needed at most one imputed
# value, so the projection cannot be driven by the imputation draw.

select_robust_features <- function(imputed_flag, max_imputed = 1) {
    which(rowSums(imputed_flag, na.rm = TRUE) <= max_imputed)
}

test_that("robust subset keeps features with 0 or 1 imputed values", {
    flag <- rbind(
        f0 = c(FALSE, FALSE, FALSE, FALSE),
        f1 = c(TRUE,  FALSE, FALSE, FALSE),
        f2 = c(TRUE,  TRUE,  FALSE, FALSE),
        f4 = c(TRUE,  TRUE,  TRUE,  TRUE)
    )
    expect_equal(rownames(flag)[select_robust_features(flag)], c("f0", "f1"))
})

test_that("robust subset is everything when nothing was imputed", {
    expect_length(select_robust_features(matrix(FALSE, nrow = 5, ncol = 4)), 5)
})

test_that("robust subset is empty when every feature is heavily imputed", {
    expect_length(select_robust_features(matrix(TRUE, nrow = 5, ncol = 4)), 0)
})

test_that("NA flags are not counted as imputed", {
    expect_length(select_robust_features(rbind(a = c(TRUE, NA, FALSE, FALSE))), 1)
})

# --- selector ordering: all, descending top-N, robust -----------------------

order_pca_panels <- function(keys) {
    n_of <- function(k) if (grepl("^top[0-9]+$", k)) as.numeric(sub("^top", "", k)) else NA_real_
    rank <- vapply(keys, function(k) {
        if (identical(k, "all")) -Inf else if (identical(k, "robust")) Inf else -n_of(k)
    }, numeric(1))
    keys[order(rank)]
}

test_that("the selector lists all proteins first and the robust set last", {
    keys <- c("top1000", "robust", "all", "top2000", "top500")
    expect_equal(order_pca_panels(keys),
                 c("all", "top2000", "top1000", "top500", "robust"))
})

test_that("top-N panels are ordered largest to smallest", {
    expect_equal(order_pca_panels(c("top500", "top2000", "top1000")),
                 c("top2000", "top1000", "top500"))
})

test_that("a missing robust panel does not disturb the rest", {
    expect_equal(order_pca_panels(c("top1000", "all")), c("all", "top1000"))
})

test_that("the module writes both new panels and keeps the top-N ones", {
    candidates <- c(
        testthat::test_path("..", "..", "R", "modules", "proteomics", "01_mod_qc_pre.R"),
        "R/modules/proteomics/01_mod_qc_pre.R"
    )
    f <- candidates[file.exists(candidates)][1]
    skip_if(is.na(f), "01_mod_qc_pre.R not found from the test working directory")

    src <- readLines(f, warn = FALSE)
    expect_gte(length(grep('"PCA_all.png"', src, fixed = TRUE)), 1)
    expect_gte(length(grep('"PCA_robust.png"', src, fixed = TRUE)), 1)
    expect_gte(length(grep("n_top_values <- c(500, 1000, 2000)", src, fixed = TRUE)), 1)
})
