#' Collect the finished RNA-seq objects that make up the run's omics bundle
#'
#' Packages what the run already computed; nothing is recalculated. Three
#' things are left out because integration does not need them and they would
#' only make the bundle large or machine-specific: the DESeq2 model object
#' (\code{de_res$dds}; the DE tables and the design record stay), the paths of
#' the rendered pathway plots (\code{pathway_res$plot_files}), and the two
#' bookkeeping fields \code{load_config()} adds to the config
#' (\code{.config_path}, \code{.config_mtime}). Annotation parts that the run
#' did not produce are omitted rather than stored as NULL.
#'
#' @param pre Preprocessed RNA object (\code{rna_pre}): raw, filtered and
#'   working matrices, the DE input, sample metadata, feature annotation, the
#'   record of what ran (\code{info}) and the contrasts.
#' @param de_res DE result (\code{rna_de_res}), with at least one non-empty
#'   table in \code{tables}.
#' @param pathway_res Enrichment result (\code{rna_pathway_res}), or NULL.
#' @param annot Gene annotation (\code{rna_annot}), or NULL when not configured.
#' @param trinotate Trinotate annotation (\code{rna_trinotate_main}), or NULL.
#' @param config Validated config of the run.
#' @return Named list of bundle parts: \code{pre}, \code{de}, \code{config},
#'   and \code{pathway}, \code{annotation}, \code{trinotate} when present.
build_rnaseq_bundle_parts <- function(pre, de_res, pathway_res = NULL, annot = NULL,
                                      trinotate = NULL, config) {
    assert_pre_contract(pre, stage = "rna")
    # Tables with no rows count as missing too: run_deseq2_de() returns
    # tables = list() when every contrast was skipped, and a header-only
    # pre-computed DE file loads as a zero-row table. Either way a bundle would
    # be advertised as complete while holding nothing to integrate.
    has_rows <- is.list(de_res) && is.list(de_res$tables) &&
        any(vapply(de_res$tables, function(t) NROW(t) > 0L, logical(1)))
    if (!has_rows) {
        stop("RNA bundle: de_res has no DE tables, so there is nothing to ",
             "integrate. Check the warnings from the DE step (a contrast whose ",
             "levels are not in the metadata is skipped); expected the output ",
             "of mod_rnaseq_de() with at least one contrast.", call. = FALSE)
    }
    assert_named_list(config, "config")

    parts <- list(
        pre        = pre,
        de         = de_res[setdiff(names(de_res), "dds")],
        pathway    = if (!is.null(pathway_res)) {
            pathway_res[setdiff(names(pathway_res), "plot_files")]
        },
        annotation = annot,
        trinotate  = trinotate,
        config     = config[!startsWith(names(config), ".")]
    )
    Filter(Negate(is.null), parts)
}
