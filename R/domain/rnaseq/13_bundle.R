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
#' @param de_res DE result (\code{rna_de_res}), with \code{tables}.
#' @param pathway_res Enrichment result (\code{rna_pathway_res}), or NULL.
#' @param annot Gene annotation (\code{rna_annot}), or NULL when not configured.
#' @param trinotate Trinotate annotation (\code{rna_trinotate_main}), or NULL.
#' @param config Validated config of the run.
#' @return Named list of bundle parts: \code{pre}, \code{de}, \code{config},
#'   and \code{pathway}, \code{annotation}, \code{trinotate} when present.
build_rnaseq_bundle_parts <- function(pre, de_res, pathway_res = NULL, annot = NULL,
                                      trinotate = NULL, config) {
    assert_pre_contract(pre, stage = "rna")
    if (!is.list(de_res) || is.null(de_res$tables)) {
        stop("RNA bundle: de_res has no DE tables; expected the output of ",
             "mod_rnaseq_de().", call. = FALSE)
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
