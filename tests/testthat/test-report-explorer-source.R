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
})

test_that("the reconciliation claims only what holds in every configuration", {
    # The fold changes may or may not equal a difference of group means:
    # de$method ttest and welch compute exactly that (05c_de_ttest.R), an
    # unblocked limma coefficient can equal it, and an externally precomputed
    # table was not computed here at all. What holds regardless is that they
    # come from the differential analysis rather than from these means.
    body <- code_only(explorer_chunk("protein-explorer-heading"))
    txt <- paste(body, collapse = "\n")

    expect_match(txt, "not recomputed from the measured-only means shown here", fixed = TRUE)
    # "may differ", never "must".
    expect_match(txt, "may therefore differ", fixed = TRUE)
})

test_that("the heading asserts no DE or imputation provenance", {
    # Every one of these was wrong in some supported configuration, and each
    # was written here before being removed. The measured-view contract does
    # not need any of them.
    body <- code_only(explorer_chunk("protein-explorer-heading"))

    expect_false(any(grepl("model estimate", body, fixed = TRUE)))
    expect_false(any(grepl("computed on a complete matrix", body, fixed = TRUE)))
    # No prose conditioned on how imputation is configured.
    expect_false(any(grepl("imputation$method", body, fixed = TRUE)))
    expect_false(any(grepl(".imp_method_exp", body, fixed = TRUE)))
})


# =============================================================================
# The model-input overlay
# =============================================================================
#
# The overlay adds ONLY the cells that were not measured and do have a
# model-input value. Everything below exists to keep that definition in one
# place and to keep an observed value out of the overlay layer entirely: the
# matrix is built in R with every non-overlay cell set to NA, so the browser is
# never in a position to decide what counts as imputed.

test_that("the overlay is paired column by column with the measured block", {
    body <- code_only(explorer_chunk("protein-explorer-data"))
    txt <- paste(body, collapse = "\n")

    # Exact <sample>.norm counterparts of the measured columns, by name.
    expect_match(txt, 'paste0(colnames(measured_expr_exp), ".norm")', fixed = TRUE)
    # Every one of them, or no overlay: a partial block is a misalignment.
    expect_match(txt, "all(norm_cols_exp %in% names(final_df))", fixed = TRUE)
    # Rows aligned by feature id, not by position.
    expect_match(txt, "model_input_exp[rownames(measured_expr_exp), , drop = FALSE]", fixed = TRUE)
    expect_match(txt, "identical(dim(model_input_exp), dim(measured_expr_exp))", fixed = TRUE)
})

test_that("an overlay cell is measured-missing and model-input-present", {
    body <- code_only(explorer_chunk("protein-explorer-data"))
    txt <- paste(body, collapse = "\n")

    expect_match(txt, "keep_ov <- is.na(meas_m) & !is.na(mod_m)", fixed = TRUE)
    # Everything else is blanked, so an observed value cannot reach the layer.
    expect_match(txt, "mod_m[!keep_ov] <- NA_real_", fixed = TRUE)
    # The overlay travels as its own key, never merged into `expr`.
    expect_match(txt, "explorer_data$imp <- imp_matrix_list", fixed = TRUE)
    expect_false(any(grepl("expr = imp_matrix_list", body, fixed = TRUE)))
})

test_that("overlay points are their own traces, never part of a box", {
    body <- code_only(explorer_chunk("protein-explorer-js"))
    txt <- paste(body, collapse = "\n")

    # A scatter trace, pushed after the box traces are complete.
    expect_match(txt, 'type: "scatter", mode: "markers"', fixed = TRUE)
    expect_match(txt, "traces.push(bpOverlayTrace(sgPts, false));", fixed = TRUE)
    expect_match(txt, "traces.push(bpOverlayTrace(mpPts, true));", fixed = TRUE)
    # The measured arrays that feed the boxes are never appended to from the
    # overlay: gValues/yVals are filled only inside the bpIsMeasured guards.
    expect_false(any(grepl("gValues.push(vals[i])", body, fixed = TRUE)))
    expect_false(any(grepl("yVals.push(geneExprData.imp", body, fixed = TRUE)))
})

test_that("the summary line still counts only measured values", {
    body <- code_only(explorer_chunk("protein-explorer-js"))
    txt <- paste(body, collapse = "\n")

    # The mean and N read `expr`, never `imp`.
    expect_match(txt, "var vals = values.filter(bpIsMeasured);", fixed = TRUE)
    expect_match(txt, "measured samples", fixed = TRUE)
    # And a group with nothing measured stays named as such, whatever the
    # overlay draws into its slot.
    expect_match(txt, "no measured values in: ", fixed = TRUE)
    expect_match(txt, "if (groups[i] === g && bpIsMeasured(values[i])) return false;", fixed = TRUE)
    expect_false(any(grepl("emptyGroups", body, fixed = TRUE) &
                     grepl("geneExprData.imp", body, fixed = TRUE)))
})

test_that("precomputed DE is detected with the pipeline's own condition", {
    body <- code_only(explorer_chunk("protein-explorer-data"))
    txt <- paste(body, collapse = "\n")

    # Same three clauses as mod_proteomics_de(); see R/modules/proteomics/02_mod_de.R.
    expect_match(txt, "prot_cfg$files$de_table", fixed = TRUE)
    expect_match(txt, "length(.de_table_files_ov) > 0", fixed = TRUE)
    expect_match(txt, "all(nzchar(unlist(.de_table_files_ov)))", fixed = TRUE)
    # And it gates the overlay.
    expect_match(txt, "!.has_precomputed_ov", fixed = TRUE)
})

test_that("every mode that cannot support an overlay is gated out", {
    body <- code_only(explorer_chunk("protein-explorer-data"))
    txt <- paste(body, collapse = "\n")

    # method "none": .norm is the measured matrix, NAs included.
    expect_match(txt, 'prot_cfg$imputation$method %||% "perseus_like"', fixed = TRUE)
    expect_match(txt, '!identical(.imp_method_ov, "none")', fixed = TRUE)
    # The #239 fallback source carries no .norm block, and none is reconstructed.
    expect_match(txt, "measured_from_final_results &&", fixed = TRUE)
    expect_match(txt, "measured_from_final_results <- TRUE", fixed = TRUE)
    # Fails closed: the flag starts FALSE and is set only inside the guards.
    expect_match(txt, "overlay_ready <- FALSE", fixed = TRUE)
    expect_equal(length(grep("overlay_ready <- TRUE", body, fixed = TRUE)), 1L)
})

test_that("the control is emitted only where an overlay exists", {
    body <- code_only(explorer_chunk("protein-explorer-ui"))
    txt <- paste(body, collapse = "\n")

    expect_match(txt, "if (isTRUE(overlay_ready)) {", fixed = TRUE)
    # The label says what the control does: it adds imputed values, it does not
    # swap the plot for the full model-input matrix.
    expect_match(txt, "Show imputed model-input values", fixed = TRUE)
    expect_false(any(grepl("disabled", body, fixed = TRUE)))
})

test_that("no imputed file discovery is reintroduced for the overlay", {
    # The overlay comes from final_results.tsv and nowhere else. This is the
    # guarantee #239 established for the measured view, extended to the second
    # layer: no glob, no imputed_once, no imputed_repetitions.
    body <- code_only(explorer_chunk("protein-explorer-data"))

    expect_false(any(grepl("list.files(", body, fixed = TRUE)))
    expect_false(any(grepl("imputed_once", body, fixed = TRUE)))
    expect_false(any(grepl("imputed_repetitions", body, fixed = TRUE)))
})
