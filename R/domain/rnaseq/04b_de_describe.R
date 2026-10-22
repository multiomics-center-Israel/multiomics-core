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
#' @return A list with \code{method} (\code{"deseq2"}, \code{"limma_voom"},
#'   \code{"edger_ql"}, \code{"limma"},
#'   \code{"precomputed"} or \code{"unknown"}), \code{label} (method name for
#'   tables), \code{design} (formula as text, or NA), \code{block_col} (or
#'   NULL), \code{block_model} ("fixed", "random" or NULL), \code{short} (one-line label) and \code{methods} (reader-facing
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
        limma_voom = {
            methods <- sprintf(paste0(
                "Differential expression was tested with limma-voom: counts were TMM-normalised ",
                "(edgeR), converted to log2-CPM with voom precision weights, and fitted with ",
                "moderated t-tests (empirical Bayes variance shrinkage) using the design %s."),
                design)
            r <- info$block_correlation
            fixed <- paired && identical(info$block_model, "fixed")
            if (fixed) {
                methods <- paste(methods, sprintf(paste0(
                    "The design includes %s as a fixed term, so groups are compared within ",
                    "each of its %s levels; each %s level beyond the first uses one residual ",
                    "degree of freedom."), block_col, info$n_blocks %||% "its", block_col))
            } else if (paired) {
                methods <- paste(methods, sprintf(paste0(
                    "Samples sharing a level of %s (%s levels) were treated as correlated: ",
                    "the within-block correlation was estimated with duplicateCorrelation() ",
                    "(consensus %.3f) and included in the fit, so %s enters as a random effect ",
                    "rather than a fixed term. voom and duplicateCorrelation were each run twice."),
                    block_col, info$n_blocks %||% "its", r, block_col))
            } else if (!is.null(info$block_requested) && nzchar(info$block_requested)) {
                methods <- paste(methods, sprintf(paste0(
                    "Blocking on %s was requested, but duplicateCorrelation() returned a ",
                    "non-finite consensus, so the model was fitted without blocking."),
                    info$block_requested))
            }
            list(label = "limma-voom",
                 short = if (fixed) {
                     sprintf("limma-voom, design %s (%s as a fixed term)", design, block_col)
                 } else if (paired) {
                     sprintf("limma-voom, design %s, block %s (consensus correlation %.3f)",
                             design, block_col, r)
                 } else {
                     sprintf("limma-voom, design %s", design)
                 },
                 methods = methods)
        },
        edger_ql = {
            methods <- sprintf(paste0(
                "Differential expression was tested with edgeR (version %s): counts were ",
                "TMM-normalised, negative-binomial dispersions were estimated with estimateDisp(), ",
                "and each contrast was tested with the quasi-likelihood F-test (glmQLFit, ",
                "glmQLFTest) using the design %s."), info$edger_version %||% "not recorded", design)
            if (paired) {
                methods <- paste(methods, sprintf(paste0(
                    "The design includes %s as a fixed term, so groups are compared within ",
                    "each of its %s levels; edgeR has no random effect, and each %s level beyond ",
                    "the first uses one residual degree of freedom."),
                    block_col, info$n_blocks %||% "its", block_col))
            }
            list(label = "edgeR",
                 short = sprintf("edgeR quasi-likelihood, design %s", design),
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

    c(list(method = method, design = design, block_col = block_col,
           block_model = if (paired) info$block_model %||% "fixed" else NULL), out)
}
