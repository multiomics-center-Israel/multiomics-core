#' RNA-seq clustering on the batch-corrected matrix
#'
#' Runs [mod_rnaseq_clustering()] a second time, on the batch-corrected
#' analysis matrix ([compute_rna_batch_correction()]), and keeps both sets of
#' cluster assignments side by side so that they can be compared.
#'
#' The original clustering (`rna_clustering_obj`, under `Clustering/`) is not
#' touched and stays the one behind the Excel order, the cluster-based
#' enrichment and the Shiny payload. This run uses the same DE genes, the same
#' config and the same code; only the expression values differ. Its outputs go
#' to `Batch_corrected_<method>/Clustering/`, with the same file names as the
#' original.
#'
#' Nothing is run unless `batch_corrected$ran` records that a correction ran.
#' Files this module wrote in an earlier run are removed first either way, so
#' they never outlive their correction.
#'
#' @param pre List returned by `preprocess_rna()`.
#' @param batch_corrected Output of [compute_rna_batch_correction()].
#' @param de_res DE results list from `mod_rnaseq_de()`.
#' @param clustering_res The original clustering ([mod_rnaseq_clustering()]).
#' @param config Full config list.
#' @param out_dir RNA output directory.
#' @return A list with `ran` (logical). When `ran` is `TRUE` also `method`,
#'   `clustering` (the result of [mod_rnaseq_clustering()] on the corrected
#'   matrix), `comparison` (see [compare_rna_cluster_assignments()], or `NULL`)
#'   and `files`.
mod_rnaseq_clustering_batch_corr <- function(pre, batch_corrected, de_res,
                                             clustering_res, config, out_dir) {
    stopifnot(is.character(out_dir), length(out_dir) == 1)
    datasets_dir <- create_legacy_output_dirs(out_dir, create = TRUE)$datasets
    comparison_file <- file.path(datasets_dir, "clustering_batch_comparison.tsv")

    # Only files, and only inside this module's own directories.
    own_roots <- list.files(out_dir, pattern = "^Batch_corrected_(limma|combat)$",
                            full.names = TRUE)
    stale <- c(list.files(own_roots, recursive = TRUE, full.names = TRUE), comparison_file)
    stale <- stale[file.exists(stale)]
    if (length(stale) > 0) file.remove(stale)

    if (!isTRUE(batch_corrected$ran)) {
        message("[rnaseq clustering] No batch correction ran -- no batch-corrected clustering.")
        return(list(ran = FALSE))
    }

    method <- batch_corrected$record$method
    # A local copy: only the clustering below sees the corrected values. The
    # counts used to gate the binary patterns (expr_filt) stay as they are.
    pre_corrected <- pre
    pre_corrected$expr_work <- batch_corrected$expr_work

    corrected_res <- mod_rnaseq_clustering(
        pre     = pre_corrected,
        de_res  = de_res,
        config  = config,
        out_dir = file.path(out_dir, sprintf("Batch_corrected_%s", method))
    )

    comparison <- compare_rna_cluster_assignments(
        clustering_res$objects$clusters, corrected_res$objects$clusters
    )
    files <- corrected_res$files
    if (!is.null(comparison)) {
        files <- c(files, save_tsv(comparison, datasets_dir, basename(comparison_file)))
    }

    list(
        ran        = TRUE,
        method     = method,
        clustering = corrected_res,
        comparison = comparison,
        files      = unique(files)
    )
}

#' Put two cluster assignments of the same genes side by side
#'
#' Cluster numbers are arbitrary in each run, so equal numbers do not mean the
#' same cluster; the table is meant to be cross-tabulated.
#'
#' @param uncorrected,corrected Named cluster vectors (gene -> cluster) from
#'   the clustering of the uncorrected and of the batch-corrected matrix, or
#'   `NULL` when that run assigned no clusters.
#' @return A data.frame with `gene`, `cluster_uncorrected` and
#'   `cluster_batch_corrected` (`NA` where a run did not assign the gene),
#'   ordered by gene; `NULL` when either run assigned no clusters.
#' @examples
#' compare_rna_cluster_assignments(c(g1 = 1, g2 = 2), c(g2 = 1, g1 = 1))
compare_rna_cluster_assignments <- function(uncorrected, corrected) {
    if (length(uncorrected) == 0 || length(corrected) == 0 ||
        is.null(names(uncorrected)) || is.null(names(corrected))) {
        return(NULL)
    }
    genes <- sort(union(names(uncorrected), names(corrected)))
    data.frame(
        gene                    = genes,
        cluster_uncorrected     = unname(uncorrected[genes]),
        cluster_batch_corrected = unname(corrected[genes]),
        stringsAsFactors = FALSE
    )
}
