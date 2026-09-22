# tests/testthat/test-report-explorer-source.R
#
# The Protein Expression Explorer's source contract.
#
# The explorer shows MEASURED expression: the filtered, normalised matrix with
# missing measurements left missing. Two files carry that -- the `<sample>`
# columns of final_results.tsv, and protein_log2_filtered_unimputed.tsv. The
# imputed matrices (protein_log2_filtered_imputed_once_*.tsv and the per-run
# files under imputed_repetitions/) carry model input, which answers a
# different question and must never stand in for the measured view.
#
# The defect these tests close: the fallback was
#
#     list.files(datasets_dir, pattern = "imputed.*\\.tsv$")[1]
#
# which is not anchored, so it matched the unimputed file as readily as the
# imputed one -- and list.files() sorts, so "...filtered_imputed_once..." came
# before "...filtered_unimputed..." and won. Whenever final_results.tsv was
# absent or carried no columns matching the metadata, the explorer plotted
# imputed intensities under a measured heading and averaged them into a mean
# indistinguishable from a measured one.
#
# Source-level assertions: this code runs only inside a knitr chunk during a
# render, so the suite cannot execute it. Scoped to a named chunk's body rather
# than grepping the whole template, because the peptide explorer does its own
# file discovery and would otherwise be caught by these assertions.

repo_file <- function(...) {
    candidates <- c(testthat::test_path("..", "..", ...), file.path(...))
    candidates[file.exists(candidates)][1]
}

# The body of one labelled chunk, without its fences.
explorer_chunk <- function(label) {
    f <- repo_file("R", "domain", "proteomics", "report_template_proteomics.Rmd")
    skip_if(is.na(f) || !file.exists(f), "proteomics report template not found")
    lines <- readLines(f, warn = FALSE)

    start <- grep(sprintf("^```\\{r %s[ ,}]", label), lines)
    expect_length(start, 1L)

    rest <- lines[(start + 1L):length(lines)]
    end <- start + which(grepl("^```\\s*$", rest))[1L]
    lines[(start + 1L):(end - 1L)]
}

# Code only. The comments in this chunk necessarily name the imputed files they
# warn against, and these assertions are about what the chunk DOES.
code_only <- function(body) body[!grepl("^\\s*#", body)]


# =============================================================================
# The measured sources, in order
# =============================================================================

test_that("final_results.tsv is the primary measured source", {
    body <- code_only(explorer_chunk("protein-explorer-data"))
    txt <- paste(body, collapse = "\n")

    expect_match(txt, 'file.path(datasets_dir, "final_results.tsv")', fixed = TRUE)
    # The `<sample>` columns, selected by name against the metadata -- not the
    # `<sample>.norm` block beside them, which is model input.
    expect_match(txt, "intersect(meta_sample_ids, names(final_df))", fixed = TRUE)
})

test_that("the fallback is the unimputed matrix, named exactly", {
    body <- code_only(explorer_chunk("protein-explorer-data"))
    txt <- paste(body, collapse = "\n")

    expect_match(txt, 'file.path(datasets_dir, "protein_log2_filtered_unimputed.tsv")',
                 fixed = TRUE)
    # Primary first: the fallback must not be able to pre-empt final_results.tsv.
    i_primary  <- grep("final_results.tsv", body, fixed = TRUE)[1]
    i_fallback <- grep("protein_log2_filtered_unimputed.tsv", body, fixed = TRUE)[1]
    expect_lt(i_primary, i_fallback)
})


# =============================================================================
# No imputed matrix can reach the measured view
# =============================================================================

test_that("the explorer names its files and never matches them", {
    # A glob is how the defect got in. Naming both sources means a new file
    # dropped into the datasets directory cannot become the measured source by
    # sorting ahead of the right one.
    body <- code_only(explorer_chunk("protein-explorer-data"))

    expect_false(any(grepl("list.files(", body, fixed = TRUE)))
    expect_false(any(grepl("imputed.*\\\\.tsv", body, fixed = TRUE)))
})

test_that("no imputed matrix is referenced as a source", {
    body <- code_only(explorer_chunk("protein-explorer-data"))

    # The single-imputation QC draw, and the per-run DE matrices.
    expect_false(any(grepl("imputed_once", body, fixed = TRUE)))
    expect_false(any(grepl("imputed_repetitions", body, fixed = TRUE)))
    # The exported model-input block, which sits beside the measured columns in
    # the very file this chunk reads.
    expect_false(any(grepl(".norm", body, fixed = TRUE)))
})


# =============================================================================
# Fail closed
# =============================================================================

test_that("the explorer stays unavailable without an aligned measured matrix", {
    body <- code_only(explorer_chunk("protein-explorer-data"))
    txt <- paste(body, collapse = "\n")

    expect_match(txt, "explorer_ready <- FALSE", fixed = TRUE)
    # Both loaders leave the matrix NULL when nothing aligns to the metadata:
    # the primary only assigns inside its match guard, the fallback nulls out.
    expect_match(txt, "length(expr_cols_exp) > 0 && id_col_exp %in% names(final_df)",
                 fixed = TRUE)
    expect_match(txt, "else norm_expr_exp <- NULL", fixed = TRUE)
    # And readiness is gated on having one.
    expect_match(txt, "!is.null(norm_expr_exp) && !is.null(meta_df_exp)", fixed = TRUE)
})

test_that("the unavailable message says what is missing and what was refused", {
    txt <- paste(explorer_chunk("protein-explorer-fallback"), collapse = "\n")

    expect_match(txt, "measured expression data could not be loaded", fixed = TRUE)
    expect_match(txt, "protein_log2_filtered_unimputed.tsv", fixed = TRUE)
    # The old wording named a matrix the explorer does not show.
    expect_false(grepl("Normalized expression matrix not found", txt, fixed = TRUE))
})


# =============================================================================
# What the reader is told
# =============================================================================

test_that("the heading states that these are measured values", {
    txt <- paste(explorer_chunk("protein-explorer-heading"), collapse = "\n")

    expect_match(txt, "measured", fixed = TRUE)
    expect_match(txt, "missing measurements left missing", fixed = TRUE)
    # The reconciliation a reader needs, and the half of it that holds in every
    # mode: a fold change is a model estimate, not a difference of these means.
    expect_match(txt, "not the difference of the two means drawn here", fixed = TRUE)
})

test_that("the imputation claim is conditional on the configured method", {
    # imputation.method = "none" returns the matrix with its NAs intact and the
    # model fits that, so a report that unconditionally told the reader their
    # statistics came from an imputed matrix would be wrong for those projects.
    body <- explorer_chunk("protein-explorer-heading")
    txt <- paste(body, collapse = "\n")

    expect_match(txt, 'prot_cfg$imputation$method %||% "none"', fixed = TRUE)
    expect_match(txt, 'if (!identical(.imp_method_exp, "none")) {', fixed = TRUE)
    # The method is named rather than assumed.
    expect_match(txt, "imputes missing values (`%s`)", fixed = TRUE)

    # The unconditional half must say nothing about which matrix produced the
    # statistics: with an externally precomputed DE table, none of the
    # pipeline's own matrices did.
    unconditional <- body[!grepl("^\\s*#", body)]
    unconditional <- unconditional[seq_len(grep("^.imp_method_exp <-", unconditional)[1])]
    expect_false(any(grepl("computed on a complete matrix", unconditional, fixed = TRUE)))
})
