# The wizard writes the RNA normalization block of new configs. Any
# prior.count it emits opts the run back into the edgeR transform (see
# resolve_rna_log_offset()), so new configs must leave it out to get the
# log2(TMM CPM + 1) default.

repo_file <- function(...) {
    candidates <- c(testthat::test_path("..", "..", ...), file.path(...))
    candidates[file.exists(candidates)][1]
}

test_that("the wizard's RNA normalization block emits no prior.count", {
    src <- readLines(repo_file("wizard.html"), warn = FALSE)

    method_line <- grep("v('rna_norm')", src, fixed = TRUE)
    expect_length(method_line, 1)
    block_starts <- grep("normalization:", src, fixed = TRUE)
    start <- max(block_starts[block_starts <= method_line])
    end <- grep("filtering:", src, fixed = TRUE)
    end <- end[end > method_line][1]
    block <- src[start:end]

    expect_true(any(grepl("method: ", block, fixed = TRUE)))
    expect_false(any(grepl("prior.count", block, fixed = TRUE)))
    expect_false(any(grepl("prior\\.count", src[grepl("yaml +=", src, fixed = TRUE)])))
})
