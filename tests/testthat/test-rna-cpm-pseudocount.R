# normalization$cpm_pseudocount: TMMlogCPM as log2(TMM CPM + x).
#
# edgeR's prior.count is in reads, scaled by library size, so prior.count = 1
# is a much smaller offset than the common log2(CPM + 1). cpm_pseudocount gives
# a fixed offset in CPM units instead. These tests pin the transform, the
# unchanged default, the config contract and the one resolver every consumer
# reads. All fixtures are synthetic.

synthetic_counts <- function() {
    m <- matrix(c(0, 5, 120, 3000,
                  2, 0, 150, 2500,
                  1, 8, 90, 4000,
                  0, 3, 200, 3500),
                nrow = 4, dimnames = list(paste0("g", 1:4), paste0("S", 1:4)))
    storage.mode(m) <- "integer"
    m
}

# Drop the method/source_type attributes so only values and dimnames compare.
values_of <- function(x) {
    attr(x, "method") <- NULL
    attr(x, "source_type") <- NULL
    x
}

test_that("cpm_pseudocount gives log2(TMM CPM + x)", {
    skip_if_not_installed("edgeR")
    counts <- synthetic_counts()
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts), method = "TMM")
    expected <- log2(edgeR::cpm(dge) + 1)

    res <- normalize_counts(counts, method = "TMMlogCPM", cpm_pseudocount = 1)

    expect_equal(values_of(res), expected)
    expect_equal(attr(res, "method"), "TMMlogCPM")
    # Unexpressed genes sit at log2(0 + 1) = 0, not at a negative prior floor.
    expect_equal(res["g1", "S1"], 0)
})

test_that("without cpm_pseudocount the edgeR prior.count transform is unchanged", {
    skip_if_not_installed("edgeR")
    counts <- synthetic_counts()
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts), method = "TMM")

    res <- normalize_counts(counts, method = "TMMlogCPM", prior.count = 2)

    expect_equal(values_of(res), edgeR::cpm(dge, log = TRUE, prior.count = 2))
})

test_that("resolve_rna_log_offset reads one key and labels it", {
    pc <- resolve_rna_log_offset(list(method = "TMMlogCPM", cpm_pseudocount = 1))
    expect_identical(pc$type, "cpm_pseudocount")
    expect_identical(pc$value, 1)
    expect_match(pc$label, "log2(TMM CPM + 1)", fixed = TRUE)

    pr <- resolve_rna_log_offset(list(method = "TMMlogCPM", prior.count = 3))
    expect_identical(pr$type, "prior_count")
    expect_identical(pr$value, 3)
    expect_match(pr$label, "prior count 3 reads")

    # Neither key set: log2(TMM CPM + 1), not edgeR's prior count.
    def <- resolve_rna_log_offset(NULL)
    expect_identical(def$type, "cpm_pseudocount")
    expect_identical(def$value, 1)
    expect_match(def$label, "log2(TMM CPM + 1)", fixed = TRUE)
    expect_identical(resolve_rna_log_offset(list(method = "TMMlogCPM"))$type, "cpm_pseudocount")
})

rna_cfg_with_norm <- function(norm) {
    list(
        id_columns = list(gene_id = "gene_id"),
        files = list(counts = "counts.csv", metadata = "metadata.csv"),
        normalization = norm,
        filtering = list(group_col = "Condition"),
        de = list(linear_fc_cutoff = 1.5)
    )
}

test_that("validate_rna_config accepts cpm_pseudocount on its own", {
    expect_true(validate_rna_config(rna_cfg_with_norm(
        list(method = "TMMlogCPM", cpm_pseudocount = 1))))
})

test_that("validate_rna_config rejects cpm_pseudocount together with prior.count", {
    expect_error(
        validate_rna_config(rna_cfg_with_norm(
            list(method = "TMMlogCPM", prior.count = 1, cpm_pseudocount = 1))),
        "not both"
    )
})

test_that("validate_rna_config rejects a zero, negative or non-numeric cpm_pseudocount", {
    expect_error(validate_rna_config(rna_cfg_with_norm(
        list(method = "TMMlogCPM", cpm_pseudocount = 0))), "above 0")
    expect_error(validate_rna_config(rna_cfg_with_norm(
        list(method = "TMMlogCPM", cpm_pseudocount = -1))), "cpm_pseudocount")
    expect_error(validate_rna_config(rna_cfg_with_norm(
        list(method = "TMMlogCPM", cpm_pseudocount = "one"))), "must be a number")
})
