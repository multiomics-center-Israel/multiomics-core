# Pin the proteomics `expr_raw` scale contract (issue #138).
#
# `expr_raw` must hold the data on its ORIGINAL input scale, not the log2 assay:
#   * scale_in = "linear" -> expr_raw is linear, expr_log2 = log2(expr_raw + 1)
#   * scale_in = "log2"   -> expr_raw == expr_log2 (both already log2)
# and `expr_raw_scale` reports which of the two applies.
#
# All fixtures are synthetic. Normalization and imputation are disabled so the
# preprocessing is fully deterministic (no RNG, no seed needed).

# --- Synthetic fixtures -------------------------------------------------------

make_prot_inputs <- function(mat, protein_id_col = "Protein.Group") {
    protein <- data.frame(
        id_placeholder = rownames(mat),
        mat,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    names(protein)[1] <- protein_id_col
    meta <- data.frame(
        SampleID  = colnames(mat),
        Condition = rep(c("Ctrl", "Treat"), length.out = ncol(mat)),
        stringsAsFactors = FALSE
    )
    list(source_type = "preprocessed", protein = protein, metadata = meta)
}

make_prot_config <- function(scale_in, is_logtransformed) {
    list(modes = list(proteomics = list(
        input         = list(format = "preprocessed"),
        engine        = "preprocessed",
        id_columns    = list(protein_id = "Protein.Group"),
        effects       = list(samples = "SampleID", color = "Condition"),
        de_table      = list(group_col = "Condition"),
        files         = list(is_logtransformed = is_logtransformed),
        scale_in      = scale_in,
        filtering     = list(min_count = 1, min_groups = 1),
        normalization = list(method = "none"),
        imputation    = list(method = "none")
    )))
}

linear_matrix <- function() {
    matrix(
        c(10, 100, 1000,
          20, 200, 2000,
          30, 300, 3000,
          40, 400, 4000),
        nrow = 3,
        dimnames = list(c("P1", "P2", "P3"), c("S1", "S2", "S3", "S4"))
    )
}

# --- Tests --------------------------------------------------------------------

test_that("scale_in = 'linear' keeps expr_raw linear and derives log2", {
    skip_if_not(exists("preprocess_proteomics"), "preprocess_proteomics not loaded")

    lin <- linear_matrix()
    inputs <- make_prot_inputs(lin)
    config <- make_prot_config(scale_in = "linear", is_logtransformed = FALSE)

    pre <- preprocess_proteomics(inputs, config)

    expect_equal(pre$expr_raw_scale, "linear")

    # expr_raw preserves the original linear intensities (no log applied).
    expect_equal(pre$expr_raw[rownames(lin), colnames(lin)], lin)

    # expr_log2 is the log2(x + 1) transform of the linear matrix.
    expect_equal(
        pre$expr_log2[rownames(lin), colnames(lin)],
        log2(lin + 1)
    )

    # The two assays share identical shape and dimnames.
    expect_identical(dim(pre$expr_raw), dim(pre$expr_log2))
    expect_identical(rownames(pre$expr_raw), rownames(pre$expr_log2))
    expect_identical(colnames(pre$expr_raw), colnames(pre$expr_log2))

    # expr_raw really is on a different scale than expr_log2 here.
    expect_false(isTRUE(all.equal(pre$expr_raw, pre$expr_log2)))
})

test_that("scale_in = 'log2' leaves expr_raw equal to the log2 assay", {
    skip_if_not(exists("preprocess_proteomics"), "preprocess_proteomics not loaded")

    log2_mat <- log2(linear_matrix() + 1)
    inputs <- make_prot_inputs(log2_mat)
    config <- make_prot_config(scale_in = "log2", is_logtransformed = TRUE)

    pre <- preprocess_proteomics(inputs, config)

    expect_equal(pre$expr_raw_scale, "log2")

    # No linear matrix exists, so expr_raw falls back to the (already log2) assay.
    expect_equal(pre$expr_raw, pre$expr_log2)

    # And no further transform was applied to the declared-log2 input.
    expect_equal(
        pre$expr_log2[rownames(log2_mat), colnames(log2_mat)],
        log2_mat
    )
})
