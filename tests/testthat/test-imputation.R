# Tests for proteomics imputation functions
# Validates impute_proteomics() dispatching and individual methods

# --- Helper: create expression matrix with NAs ---
make_expr_with_na <- function(nrow = 10, ncol = 6, na_frac = 0.2, seed = 42) {
    set.seed(seed)
    mat <- matrix(rnorm(nrow * ncol, mean = 20, sd = 3),
                  nrow = nrow, ncol = ncol,
                  dimnames = list(paste0("PROT", seq_len(nrow)),
                                  paste0("S", seq_len(ncol))))
    # Introduce NAs (MNAR-like: lower values more likely missing)
    n_na <- round(nrow * ncol * na_frac)
    na_idx <- sample(length(mat), n_na)
    mat[na_idx] <- NA
    mat
}

# Minimal config for imputation tests
make_imp_cfg <- function(method = "perseus") {
    list(
        modes = list(
            proteomics = list(
                imputation = list(
                    method = method,
                    width = 0.3,
                    downshift = 1.8,
                    no_repetitions = 3,
                    min_no_passed = 2,
                    multi_imputation = FALSE,
                    dep2_method = "MinDet",
                    dep2_random_seed = 1,
                    qrilc_random_seed = 1
                )
            )
        ),
        params = list(seed = 42)
    )
}

# ============================================================
# perseus_impute_with_flags()
# ============================================================

test_that("perseus_impute_with_flags fills all NAs", {
    mat <- make_expr_with_na()
    na_before <- sum(is.na(mat))
    expect_true(na_before > 0)

    res <- perseus_impute_with_flags(mat, width = 0.3, downshift = 1.8)

    expect_true(is.list(res))
    expect_true("imputed" %in% names(res))
    expect_true("imputed_flag" %in% names(res))
    expect_equal(sum(is.na(res$imputed)), 0)
    expect_equal(dim(res$imputed), dim(mat))
})

test_that("perseus_impute_with_flags flags match original NAs", {
    mat <- make_expr_with_na()
    na_positions <- is.na(mat)
    res <- perseus_impute_with_flags(mat, width = 0.3, downshift = 1.8)
    expect_equal(res$imputed_flag, na_positions)
})

test_that("perseus imputed values are downshifted from observed", {
    mat <- make_expr_with_na()
    res <- perseus_impute_with_flags(mat, width = 0.3, downshift = 1.8)

    # Imputed values should generally be lower than observed mean
    for (j in seq_len(ncol(mat))) {
        obs_mean <- mean(mat[, j], na.rm = TRUE)
        imp_vals <- res$imputed[is.na(mat[, j]), j]
        if (length(imp_vals) > 0) {
            expect_true(mean(imp_vals) < obs_mean)
        }
    }
})

# ============================================================
# impute_proteomics() dispatcher
# ============================================================

test_that("impute_proteomics method='none' returns input unchanged", {
    mat <- make_expr_with_na()
    cfg <- make_imp_cfg("none")
    res <- impute_proteomics(mat, cfg$modes$proteomics)
    expect_equal(res, mat)
})

test_that("impute_proteomics method='perseus' fills NAs", {
    mat <- make_expr_with_na()
    cfg <- make_imp_cfg("perseus")
    res <- impute_proteomics(mat, cfg$modes$proteomics)
    expect_equal(sum(is.na(res)), 0)
    expect_equal(dim(res), dim(mat))
    expect_equal(rownames(res), rownames(mat))
})

test_that("impute_proteomics method='dep2' fills NAs", {
    mat <- make_expr_with_na()
    cfg <- make_imp_cfg("dep2")
    res <- impute_proteomics(mat, cfg$modes$proteomics)
    expect_equal(sum(is.na(res)), 0)
    expect_equal(dim(res), dim(mat))
})

test_that("impute_proteomics return_flags=TRUE returns list", {
    mat <- make_expr_with_na()
    cfg <- make_imp_cfg("perseus")
    res <- impute_proteomics(mat, cfg$modes$proteomics, return_flags = TRUE)
    expect_true(is.list(res))
    expect_true("imputed" %in% names(res))
    expect_true("imputed_flag" %in% names(res))
})

test_that("impute_proteomics preserves non-NA values", {
    mat <- make_expr_with_na()
    not_na <- !is.na(mat)
    cfg <- make_imp_cfg("perseus")
    res <- impute_proteomics(mat, cfg$modes$proteomics)

    # Original non-NA values should be unchanged
    expect_equal(res[not_na], mat[not_na])
})
