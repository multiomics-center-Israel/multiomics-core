# tests/testthat/test-proteomics-vsn.R
#
# normalization$method: "vsn" for proteomics. The loader hands preprocessing a
# log2 matrix with NA for missing values; normalize_proteomics_vsn() has to
# undo the log2 for vsn, keep the matrix shape and names, and leave missing
# values missing.

vsn_expr <- function(n_prot = 200, n_samp = 6, seed = 11) {
    withr::with_seed(seed, {
        lin <- matrix(stats::rlnorm(n_prot * n_samp, meanlog = 15, sdlog = 2),
                      nrow = n_prot,
                      dimnames = list(paste0("P", seq_len(n_prot)),
                                      paste0("S", seq_len(n_samp))))
        # A per-sample scale offset, the kind of difference vsn calibrates away.
        lin <- sweep(lin, 2, 2^seq(0, 1, length.out = n_samp), `*`)
        expr <- log2(lin)
        expr[sample(length(expr), 0.05 * length(expr))] <- NA_real_
        expr
    })
}

test_that("normalize_proteomics_vsn keeps shape, names and missing values", {
    skip_if_not_installed("vsn")
    expr <- vsn_expr()
    out <- normalize_proteomics_vsn(expr)

    expect_true(is.matrix(out))
    expect_identical(dim(out), dim(expr))
    expect_identical(dimnames(out), dimnames(expr))
    expect_identical(is.na(out), is.na(expr))
    expect_true(all(is.finite(out[!is.na(expr)])))
})

test_that("normalize_proteomics_vsn is deterministic", {
    # vsn only draws random numbers when subsample > 0; the helper relies on
    # the default of 0, so two calls must agree exactly without any seed.
    skip_if_not_installed("vsn")
    expr <- vsn_expr()
    expect_identical(normalize_proteomics_vsn(expr), normalize_proteomics_vsn(expr))
})

test_that("preprocess_proteomics dispatches normalization$method = 'vsn'", {
    skip_if_not_installed("vsn")
    expr <- vsn_expr(n_samp = 6)
    meta <- data.frame(SampleID = colnames(expr),
                       group = rep(c("A", "B"), each = 3),
                       stringsAsFactors = FALSE)
    row_data <- data.frame(Protein.Group = rownames(expr), stringsAsFactors = FALSE)

    cfg_for <- function(method) {
        list(
            modes = list(proteomics = list(
                id_columns = list(protein_id = "Protein.Group", sample_col = "SampleID"),
                effects = list(color = "group", samples = "SampleID"),
                filtering = list(remove_contaminants = FALSE, min_count = 1,
                                 min_groups = 1),
                normalization = list(method = method),
                imputation = list(method = "none")
            )),
            params = list(seed = 1)
        )
    }

    # The loader is bypassed: this test is about the dispatcher, not input
    # parsing. The pipeline functions are sourced globals rather than a
    # namespace, so local_mocked_bindings() has nothing to rebind; assign into
    # preprocess_proteomics()'s own environment and restore on exit (as in
    # test-gmt-kegg-class-exclusion.R).
    target <- environment(preprocess_proteomics)
    stubs <- list(
        validate_proteomics_inputs = function(inputs, cfg) invisible(TRUE),
        get_proteomics_expression_matrix = function(inputs, config) {
            list(assay_log2 = expr, row_data = row_data, col_data = meta)
        }
    )
    old <- mget(names(stubs), envir = target)
    withr::defer(for (nm in names(old)) assign(nm, old[[nm]], envir = target))
    for (nm in names(stubs)) assign(nm, stubs[[nm]], envir = target)

    res <- preprocess_proteomics(inputs = list(), config = cfg_for("vsn"))
    expected <- normalize_proteomics_vsn(res$expr_raw[rownames(res$expr_filt),
                                                      colnames(res$expr_filt)])
    expect_equal(res$expr_filt, expected)

    expect_error(preprocess_proteomics(inputs = list(), config = cfg_for("quantile")),
                 "Supported: 'none', 'median', 'vsn'", fixed = TRUE)
})
