#' Resolve the blocking factor for a differential-expression fit
#'
#' Reads `de$block_col` from a mode's config branch and returns it as a factor
#' aligned to the sample order of `meta_aligned`, or NULL when no blocking was
#' requested. Refuses degenerate blockings rather than letting the fit return a
#' meaningless estimate: proteomics passes the factor to
#' [limma::duplicateCorrelation()], RNA-seq adds it to the DESeq2 design.
#'
#' @param meta_aligned Sample metadata, already ordered to match the expression
#'   matrix columns.
#' @param mode_cfg A mode's config branch (e.g. `modes$proteomics`); only its
#'   `de$block_col` is read.
#' @param sample_col Name of the sample identifier column, used in messages.
#' @return A factor with one entry per sample, or NULL if blocking is disabled.
resolve_de_block <- function(meta_aligned, mode_cfg, sample_col) {
    block_col <- mode_cfg$de$block_col
    if (is.null(block_col) || !nzchar(block_col)) return(NULL)

    if (!block_col %in% colnames(meta_aligned)) {
        stop("de$block_col is '", block_col, "' but that column is not in the sample metadata.\n",
             "  Available columns: ", paste(colnames(meta_aligned), collapse = ", "), "\n",
             "  Set de$block_col to the column holding the repeated-measures unit ",
             "(donor, line, subject), or remove it to fit without blocking.")
    }

    block <- factor(meta_aligned[[block_col]])
    n_samples <- nrow(meta_aligned)

    if (anyNA(block)) {
        stop("de$block_col '", block_col, "' has missing values for ", sum(is.na(block)),
             " of ", n_samples, " samples.\n",
             "  Every sample needs a block label; fill them in or drop those samples ",
             "via the mode's sample_filter.")
    }
    if (nlevels(block) < 2) {
        stop("de$block_col '", block_col, "' has a single level, so there is nothing to block on.\n",
             "  Remove de$block_col, or point it at a column that varies across ", sample_col, ".")
    }
    if (nlevels(block) >= n_samples) {
        stop("de$block_col '", block_col, "' has ", nlevels(block), " levels for ", n_samples,
             " samples, so no two samples share a block.\n",
             "  Blocking needs repeated measures within a block; this looks ",
             "like a sample identifier rather than a blocking variable.")
    }

    block
}
