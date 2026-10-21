#' Describe the RNA-seq differential-expression model that actually ran
#'
#' The one reading of the \code{info} list the DE step returns, shared by the
#' HTML report, the fact sheet and the AI commentary so that all three name the
#' same model. The method, design and blocking column come from what was
#' fitted, not from \code{de$method} or \code{de$block_col} in the config.
#'
#' @param info The \code{info} list of \code{mod_rnaseq_de()}'s result, or NULL
#'   when it is not available (e.g. a report re-knit by hand); the description
#'   then says so instead of guessing.
#' @return A list with \code{method} (\code{"deseq2"}, \code{"limma"},
#'   \code{"precomputed"} or \code{"unknown"}), \code{label} (method name for
#'   tables), \code{design} (formula as text, or NA), \code{block_col} (or
#'   NULL), \code{short} (one-line label) and \code{methods} (reader-facing
#'   sentences).
describe_rna_de <- function(info) {
    method <- info$method %||% "unknown"
    design <- info$design %||% NA_character_
    block_col <- info$block_col
    paired <- !is.null(block_col) && nzchar(block_col)

    out <- switch(method,
        deseq2 = {
            mode <- info$deseq_mode %||% "default"
            methods <- sprintf(paste0(
                "Differential expression was tested with DESeq2 (Wald test, %s mode) ",
                "using the design %s."), mode, design)
            if (paired) {
                methods <- paste(methods, sprintf(paste0(
                    "The design includes %s as a fixed term, so groups are compared within ",
                    "each of its %s levels. DESeq2 has no random effect; each %s level beyond ",
                    "the first uses one residual degree of freedom."),
                    block_col, info$n_blocks %||% "its", block_col))
            }
            list(label = "DESeq2",
                 short = sprintf("DESeq2 (%s mode), design %s", mode, design),
                 methods = methods)
        },
        limma = list(
            label = "limma",
            short = sprintf("limma on preprocessed expression, design %s", design),
            methods = sprintf(paste0(
                "Differential expression was tested with limma (moderated t-test with ",
                "empirical Bayes variance shrinkage) on the preprocessed expression values, ",
                "using the design %s."), design)),
        precomputed = list(
            label = "pre-computed",
            short = "pre-computed DE tables (no model fitted by this pipeline)",
            methods = paste0(
                "Differential expression tables were supplied pre-computed; this pipeline ",
                "did not fit a model.")),
        list(
            label = "not recorded",
            short = "not recorded for this report",
            methods = "The differential-expression model was not recorded for this report.")
    )

    c(list(method = method, design = design, block_col = block_col), out)
}
