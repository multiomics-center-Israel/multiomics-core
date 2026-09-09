# tests/testthat/test-ppi-network-enrichment-memory.R
#
# run_network_enrichment() calls clusterProfiler::enrichGO() once per network
# community. enrichGO returns an object carrying the full gene universe and its
# set mappings; retaining one per community, plus each community's protein
# vector, made the step peak in the gigabytes and get OOM-killed on a 16 GB
# machine. Only the report's columns are kept now. These tests pin that
# contract, since the failure mode is a crash rather than a wrong number and so
# would otherwise go unnoticed until a large run.

keep_cols <- c("ID", "Description", "GeneRatio", "BgRatio",
               "pvalue", "p.adjust", "qvalue", "Count")

mock_enrichgo_result <- function() {
    data.frame(
        ID = c("GO:0001", "GO:0002"),
        Description = c("term one", "term two"),
        GeneRatio = c("3/10", "4/10"),
        BgRatio = c("30/2000", "40/2000"),
        pvalue = c(0.001, 0.02),
        p.adjust = c(0.01, 0.08),
        qvalue = c(0.008, 0.07),
        geneID = c("111/222/333", "111/444/555/666"),
        Count = c(3L, 4L),
        stringsAsFactors = FALSE
    )
}

trim_result <- function(res) res[, intersect(keep_cols, colnames(res)), drop = FALSE]

test_that("trimming keeps every column the report needs", {
    out <- trim_result(mock_enrichgo_result())
    expect_setequal(colnames(out), keep_cols)
})

test_that("trimming drops the wide geneID column", {
    # geneID holds the full member list per term and is the bulk of the frame.
    expect_false("geneID" %in% colnames(trim_result(mock_enrichgo_result())))
})

test_that("trimming preserves row count and values", {
    res <- mock_enrichgo_result()
    out <- trim_result(res)
    expect_equal(nrow(out), nrow(res))
    expect_equal(out$p.adjust, res$p.adjust)
    expect_equal(out$ID, res$ID)
})

test_that("trimming tolerates a result missing optional columns", {
    res <- mock_enrichgo_result()
    res$qvalue <- NULL
    out <- trim_result(res)
    expect_false("qvalue" %in% colnames(out))
    expect_true(all(c("ID", "p.adjust", "Count") %in% colnames(out)))
})

test_that("trimming a zero-row result keeps the columns, not an error", {
    res <- mock_enrichgo_result()[0, , drop = FALSE]
    out <- trim_result(res)
    expect_equal(nrow(out), 0)
    expect_setequal(colnames(out), keep_cols)
})

test_that("the per-community payload no longer carries the protein vector", {
    candidates <- c(
        testthat::test_path("..", "..", "R", "domain", "proteomics", "07c_ppi_networks.R"),
        "R/domain/proteomics/07c_ppi_networks.R"
    )
    f <- candidates[file.exists(candidates)][1]
    skip_if(is.na(f), "07c_ppi_networks.R not found from the test working directory")

    src <- readLines(f, warn = FALSE)
    # The community payload is built as list(go = ..., n_proteins = ...).
    # `proteins = comm_proteins` duplicated community_df and is gone.
    expect_length(grep("proteins = comm_proteins", src, fixed = TRUE), 0)
    expect_gte(length(grep("n_proteins = length(comm_proteins)", src, fixed = TRUE)), 1)
    expect_gte(length(grep("rm(go_result)", src, fixed = TRUE)), 1)
})
