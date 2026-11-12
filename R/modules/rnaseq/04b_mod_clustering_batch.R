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
#' The correction kept the group the DE compares. Before clustering, the
#' relation of `clustering$group_col` to that group is checked
#' ([classify_rna_clustering_group()]). When the clustering group is a batch
#' column, or cannot be told apart from batch within the kept groups, the
#' clustering of the corrected values is skipped: that is a property of the
#' design, not a failure, so the pipeline goes on and the skip is recorded.
#' When the clustering group is identifiable but was not in the correction's
#' model, the clustering runs with a warning. Errors of the clustering itself
#' are not caught.
#'
#' What happened is written to `Datasets/clustering_batch_correction.tsv`
#' ([rna_batch_clustering_record()]), which the report reads. Nothing is
#' recorded when no correction ran or clustering is disabled. Files this
#' module wrote in an earlier run, and the directories they leave empty, are
#' removed first in every case, so they never outlive their run.
#'
#' @param pre List returned by `preprocess_rna()`.
#' @param batch_corrected Output of [compute_rna_batch_correction()].
#' @param de_res DE results list from `mod_rnaseq_de()`.
#' @param clustering_res The original clustering ([mod_rnaseq_clustering()]).
#' @param config Full config list.
#' @param out_dir RNA output directory.
#' @return A list with `status` (`"not_run"`, `"completed"` or `"skipped"`),
#'   `record` (the written record, or `NULL` when not run), `clustering` (the
#'   result of [mod_rnaseq_clustering()] on the corrected matrix, or `NULL`),
#'   `comparison` (see [compare_rna_cluster_assignments()], or `NULL`) and
#'   `files` (paths written, existing on disk; `character(0)` when not run).
mod_rnaseq_clustering_batch_corr <- function(pre, batch_corrected, de_res,
                                             clustering_res, config, out_dir) {
    stopifnot(is.character(out_dir), length(out_dir) == 1)
    datasets_dir <- create_legacy_output_dirs(out_dir, create = TRUE)$datasets
    comparison_file <- file.path(datasets_dir, "clustering_batch_comparison.tsv")
    record_file <- file.path(datasets_dir, "clustering_batch_correction.tsv")
    remove_rna_batch_clustering_outputs(out_dir, c(comparison_file, record_file))

    cfg <- config$modes$rna
    cl <- cfg$clustering
    if (!isTRUE(batch_corrected$ran) || is.null(cl) || isFALSE(cl$enabled)) {
        message("[rnaseq clustering] No batch correction ran, or clustering is disabled -- ",
                "no batch-corrected clustering.")
        return(list(status = "not_run", record = NULL, clustering = NULL,
                    comparison = NULL, files = character(0)))
    }

    rec <- batch_corrected$record
    method <- rec$method
    kept_col <- rec$group_col
    batch_cols <- strsplit(rec$batch_cols, ";", fixed = TRUE)[[1]]
    clustering_col <- get_clustering_group_col(cfg, pre$meta)

    # The samples of the corrected matrix, as the correction aligned them.
    sample_col <- cfg$id_columns$sample_col %||% "SampleID"
    meta_ids <- if (sample_col %in% colnames(pre$meta)) {
        as.character(pre$meta[[sample_col]])
    } else {
        rownames(pre$meta)
    }
    rows <- match(colnames(batch_corrected$expr_work), meta_ids)
    relation <- classify_rna_clustering_group(
        kept           = pre$meta[[kept_col]][rows],
        clustering     = pre$meta[[clustering_col]][rows],
        batch          = lapply(batch_cols, function(b) pre$meta[[b]][rows]),
        kept_col       = kept_col,
        clustering_col = clustering_col,
        batch_cols     = batch_cols
    )
    finish <- function(status, reason, clustering = NULL, comparison = NULL,
                       files = character(0)) {
        record <- rna_batch_clustering_record(status, reason, method, kept_col,
                                              clustering_col, relation$relation, batch_cols)
        files <- unique(c(files, save_tsv(record, datasets_dir, basename(record_file))))
        list(status = status, record = record, clustering = clustering,
             comparison = comparison, files = files[file.exists(files)])
    }

    if (!relation$proceed) {
        skipped <- finish("skipped", relation$relation)
        message("[rnaseq clustering] ", rna_batch_clustering_note(skipped$record))
        return(skipped)
    }
    if (identical(relation$relation, "not_protected")) {
        warning("[rnaseq clustering] ", rna_batch_clustering_note(
            rna_batch_clustering_record("completed", NA_character_, method, kept_col,
                                        clustering_col, relation$relation, batch_cols)),
            call. = FALSE)
    }

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
    if (length(corrected_res$files) == 0) {
        skipped <- finish("skipped", "no_output", clustering = corrected_res)
        message("[rnaseq clustering] ", rna_batch_clustering_note(skipped$record))
        return(skipped)
    }

    comparison <- compare_rna_cluster_assignments(
        clustering_res$objects$clusters, corrected_res$objects$clusters
    )
    files <- corrected_res$files
    if (!is.null(comparison)) {
        files <- c(files, save_tsv(comparison, datasets_dir, basename(comparison_file)))
    }
    finish("completed", NA_character_, clustering = corrected_res,
           comparison = comparison, files = files)
}

#' Remove what an earlier batch-corrected clustering left behind
#'
#' Removes the files inside this module's own `Batch_corrected_<method>`
#' directories and the given files, then the directories that this leaves
#' empty. An emptied `Partition_clustering_<k>_clusters` directory of an
#' earlier run with another `k` would otherwise stay next to the new one.
#' Nothing outside those directories is touched, and a directory that still
#' holds a file is kept.
#'
#' @param out_dir RNA output directory.
#' @param files Further files of this module to remove (its tables).
#' @return `invisible(NULL)`.
remove_rna_batch_clustering_outputs <- function(out_dir, files = character(0)) {
    own_roots <- list.files(out_dir, pattern = "^Batch_corrected_(limma|combat)$",
                            full.names = TRUE)
    own_roots <- own_roots[dir.exists(own_roots)]
    stale <- c(list.files(own_roots, recursive = TRUE, full.names = TRUE, all.files = TRUE),
               files)
    stale <- stale[file.exists(stale) & !dir.exists(stale)]
    if (length(stale) > 0) file.remove(stale)

    # Deepest first, so a parent emptied by its children goes too.
    dirs <- c(list.dirs(own_roots, recursive = TRUE, full.names = TRUE))
    dirs <- unique(dirs[order(nchar(dirs), decreasing = TRUE)])
    for (d in dirs) {
        if (length(list.files(d, all.files = TRUE, no.. = TRUE)) == 0) {
            unlink(d, recursive = TRUE)
        }
    }
    invisible(NULL)
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
