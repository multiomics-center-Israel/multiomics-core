# Tests for metabolomics normalization pipeline
# Validates normalize_samples(), transform_metab(), scale_metab()

# --- Helper: create a small intensity matrix ---
make_intensity_matrix <- function(nrow = 5, ncol = 4, seed = 42) {
    set.seed(seed)
    mat <- matrix(abs(rnorm(nrow * ncol, mean = 1000, sd = 300)),
                  nrow = nrow, ncol = ncol,
                  dimnames = list(paste0("feat", seq_len(nrow)),
                                  paste0("S", seq_len(ncol))))
    mat
}

# ============================================================
# normalize_samples()
# ============================================================

test_that("normalize_samples with method='none' returns input unchanged", {
    mat <- make_intensity_matrix()
    res <- normalize_samples(mat, method = "none")
    expect_equal(res, mat)
})

test_that("normalize_samples with method='sum' preserves dimensions", {
    mat <- make_intensity_matrix()
    res <- normalize_samples(mat, method = "sum")
    expect_equal(dim(res), dim(mat))
    expect_equal(rownames(res), rownames(mat))
    expect_equal(colnames(res), colnames(mat))
})

test_that("normalize_samples 'sum' equalizes column sums", {
    mat <- make_intensity_matrix()
    res <- normalize_samples(mat, method = "sum")
    col_sums <- colSums(res)
    # All column sums should be equal (scaled to median)
    expect_true(max(col_sums) - min(col_sums) < 1e-8)
})

test_that("normalize_samples 'median' preserves dimensions and is numeric", {
    mat <- make_intensity_matrix()
    res <- normalize_samples(mat, method = "median")
    expect_equal(dim(res), dim(mat))
    expect_true(is.numeric(res))
    expect_false(any(is.na(res)))
})

test_that("normalize_samples 'pqn' preserves dimensions", {
    mat <- make_intensity_matrix()
    res <- normalize_samples(mat, method = "pqn")
    expect_equal(dim(res), dim(mat))
    expect_true(is.numeric(res))
})

test_that("normalize_samples rejects invalid method", {
    mat <- make_intensity_matrix()
    expect_error(normalize_samples(mat, method = "invalid_method"))
})

# ============================================================
# transform_metab()
# ============================================================

test_that("transform_metab 'none' returns input unchanged", {
    mat <- make_intensity_matrix()
    res <- transform_metab(mat, method = "none")
    expect_equal(res, mat)
})

test_that("transform_metab 'log2' produces correct values", {
    mat <- make_intensity_matrix()
    res <- transform_metab(mat, method = "log2", pseudocount = 1)
    expect_equal(res, log2(mat + 1))
    expect_true(all(is.finite(res)))
})

test_that("transform_metab 'log10' produces correct values", {
    mat <- make_intensity_matrix()
    res <- transform_metab(mat, method = "log10", pseudocount = 1)
    expect_equal(res, log10(mat + 1))
})

test_that("transform_metab 'log10' produces correct values", {
    mat <- make_intensity_matrix()
    res <- transform_metab(mat, method = "log10", pseudocount = 1)
    expect_equal(res, log10(mat + 1))
    expect_true(all(is.finite(res)))
})

# ============================================================
# scale_metab()
# ============================================================

test_that("scale_metab 'none' returns input unchanged", {
    mat <- make_intensity_matrix()
    res <- scale_metab(mat, method = "none")
    expect_equal(res, mat)
})

test_that("scale_metab 'auto' centers rows to mean zero", {
    mat <- make_intensity_matrix()
    res <- scale_metab(mat, method = "auto")
    row_means <- rowMeans(res)
    expect_true(all(abs(row_means) < 1e-10))
})

test_that("scale_metab 'pareto' preserves dimensions", {
    mat <- make_intensity_matrix()
    res <- scale_metab(mat, method = "pareto")
    expect_equal(dim(res), dim(mat))
    expect_true(is.numeric(res))
})

test_that("scale_metab 'range' preserves dimensions", {
    mat <- make_intensity_matrix()
    res <- scale_metab(mat, method = "range")
    expect_equal(dim(res), dim(mat))
})

# ============================================================
# apply_normalization_pipeline()
# ============================================================

test_that("apply_normalization_pipeline runs full 3-step pipeline", {
    mat <- make_intensity_matrix()
    norm_cfg <- list(
        sample_norm = "sum",
        transform   = "log2",
        scaling     = "auto",
        pseudocount = 1,
        is_ref_col  = NULL
    )
    res <- apply_normalization_pipeline(mat, norm_cfg)
    expect_true(is.list(res))
    expect_true("expr_norm" %in% names(res))
    expect_equal(dim(res$expr_norm), dim(mat))
    expect_true(is.numeric(res$expr_norm))
    expect_false(any(is.na(res$expr_norm)))
})
