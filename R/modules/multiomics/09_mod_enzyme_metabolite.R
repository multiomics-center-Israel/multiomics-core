#' Module: enzyme-metabolite pairs
#'
#' Builds the pair table and writes it beside the cross-omics enrichment
#' outputs, where the report looks for it. Computation and the write are kept
#' apart in the domain layer; this wrapper is what a target calls.
#'
#' Never raises: an annotation lookup that cannot reach KEGG, or a run with
#' only one of the two layers, returns NULL and leaves no file behind.
#'
#' @param de_results Named list of DE results per omics layer.
#' @param harmonization_res Harmonization result.
#' @param config Full config object.
#' @param out_dir Cross-enrichment output directory.
#' @return List with \code{pairs} (data frame or NULL) and \code{file} (path or
#'   NULL).
mod_multiomics_enzyme_metabolite <- function(de_results, harmonization_res,
                                             config, out_dir) {
    message("\n=== Enzyme-Metabolite Pairs ===\n")

    pairs <- tryCatch(
        build_enzyme_metabolite_pairs(
            de_results = de_results,
            harmonization_res = harmonization_res,
            config = config,
            out_dir = out_dir
        ),
        error = function(e) {
            message("  Enzyme-metabolite pairs failed: ", conditionMessage(e))
            NULL
        }
    )

    # Written on every path, including the NULL one, so a rerun that produces
    # nothing clears the previous run's file rather than leaving the report
    # showing it as current.
    file <- write_enzyme_metabolite_pairs(pairs, out_dir)
    list(pairs = pairs, file = file)
}
