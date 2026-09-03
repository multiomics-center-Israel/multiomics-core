#' Load RNA inputs
#'
#' Supports both raw count matrices and tximport objects.
#' - For raw counts: specify files$counts pointing to CSV/TSV
#' - For tximport: specify files$txi pointing to RDS file containing tximport object
#'
#' @param config list as returned by load_config()
#' @return list with one of: (counts, ...) for raw counts OR (txi, ...) for tximport
load_rna_inputs <- function(config) {
    inputs <- load_omics_inputs(config, mode = "rna")

    # Check if txi was loaded
    if (!is.null(inputs$txi)) {
        # Validate tximport structure
        if (!is_valid_tximport_structure(inputs$txi, validate_only = TRUE)) {
            stop(
                "[load_rna_inputs] File loaded as 'txi' is not a valid tximport object. ",
                "Expected list with 'counts', 'abundance', 'length' matrices.",
                call. = FALSE
            )
        }
        message("[load_rna_inputs] Loaded tximport object from RDS")
        inputs$source_type <- "tximport"
    } else if (!is.null(inputs$counts)) {
        message("[load_rna_inputs] Loaded raw count matrix")
        inputs$source_type <- "matrix"
    } else if (!is.null(inputs$preprocessed_counts)) {
        # Fallback: preprocessed_counts provided instead of raw counts
        message("[load_rna_inputs] Using preprocessed_counts as count matrix")
        inputs$counts <- inputs$preprocessed_counts
        inputs$source_type <- "matrix"
    }

    inputs
}

#' Validate RNA inputs (for raw count matrix path only)
#'
#' @param inputs List returned by \code{load_rna_inputs()}
#' @param cfg    The \code{config$modes$rna} sub-list
validate_rna_inputs <- function(inputs, cfg) {
    # This validation is only for raw count matrices
    # tximport validation is handled separately by validate_tximport()
    if (identical(inputs$source_type, "tximport") || !is.null(inputs$txi)) {
        # For tximport, metadata validation only
        sample_col <- cfg$id_columns$sample_col %||% "SampleID"
        meta <- inputs$metadata
        if (!sample_col %in% names(meta)) {
            stop("metadata missing sample column: ", sample_col)
        }
        return(invisible(TRUE))
    }

    gene_id_col <- cfg$id_columns$gene_id
    sample_col <- cfg$id_columns$sample_col %||% "SampleID"

    if (!gene_id_col %in% names(inputs$counts)) {
        stop("RNA counts missing gene id column: ", gene_id_col)
    }

    meta <- inputs$metadata
    if (!sample_col %in% names(meta)) {
        stop("metadata missing sample column: ", sample_col)
    }

    # Basic checks for counts (numeric columns only — exclude annotation columns)
    non_id_cols <- setdiff(names(inputs$counts), gene_id_col)
    is_numeric <- vapply(inputs$counts[, non_id_cols, drop = FALSE],
                         function(x) is.numeric(x) || is.integer(x),
                         logical(1))
    sample_cols <- non_id_cols[is_numeric]
    if (length(sample_cols) == 0) stop("Counts table has no numeric sample columns.")

    # Check alignment consistency (informational, strict check happens later)
    meta_samples <- unique(as.character(meta[[sample_col]]))
    count_samples <- sample_cols

    missing_in_counts <- setdiff(meta_samples, count_samples)
    if (length(missing_in_counts) > 0) {
        stop(
            "metadata contains samples not present in counts: ",
            paste(missing_in_counts, collapse = ", ")
        )
    }

    invisible(TRUE)
}


#' Import a folder of RSEM results as a tximport object
#'
#' Reads per-sample RSEM gene-level result files and summarises them into a
#' tximport-style list (\code{counts}, \code{abundance}, \code{length}) ready to
#' be saved as the \code{.rds} the pipeline consumes via
#' \code{modes.rna.files.txi}. This is a pre-pipeline helper: it does no I/O of
#' its own beyond reading the RSEM files — persist the result yourself with
#' \code{saveRDS()} so the write stays a separate, explicit step.
#'
#' @param rsem_dir Path to the directory holding the RSEM result files.
#' @param pattern Regex matching the per-sample files to import. Defaults to
#'   gene-level RSEM output (\code{"\\.genes\\.results$"}); pass an isoform
#'   pattern together with \code{tx2gene} to summarise transcripts to genes.
#' @param sample_names Optional character vector of sample IDs, one per matched
#'   file in sorted order. When \code{NULL} (default) IDs are derived by
#'   stripping \code{pattern} from each file's basename. These become the
#'   tximport column names and must match the metadata \code{SampleID} column.
#' @param tx2gene Optional two-column data.frame (transcript, gene) used to
#'   summarise isoform-level RSEM output to genes. Leave \code{NULL} for
#'   gene-level input.
#' @param fix_zero_length When \code{TRUE} (default) replaces effective lengths
#'   of 0 with 1. RSEM reports length 0 for features with no reads, which
#'   otherwise makes \code{DESeq2::DESeqDataSetFromTximport()} abort — this is
#'   the fix recommended in the tximport documentation.
#' @return A tximport list with \code{counts}, \code{abundance} and
#'   \code{length} matrices (genes x samples) plus \code{countsFromAbundance},
#'   validated against the structure the pipeline expects.
#' @examples
#' \dontrun{
#' txi <- load_rsem_as_tximport("data/rsem")
#' saveRDS(txi, "data/rsem/txi.rds") # then point modes.rna.files.txi here
#' }
load_rsem_as_tximport <- function(rsem_dir,
                                  pattern = "\\.genes\\.results$",
                                  sample_names = NULL,
                                  tx2gene = NULL,
                                  fix_zero_length = TRUE) {
    if (!requireNamespace("tximport", quietly = TRUE)) {
        stop("[load_rsem_as_tximport] Package 'tximport' is required but not installed.", call. = FALSE)
    }
    if (!dir.exists(rsem_dir)) {
        stop("[load_rsem_as_tximport] Directory does not exist: ", rsem_dir, call. = FALSE)
    }

    files <- sort(list.files(rsem_dir, pattern = pattern, full.names = TRUE))
    if (length(files) == 0) {
        stop(
            "[load_rsem_as_tximport] No files matching /", pattern, "/ found in: ", rsem_dir,
            ". Gene-level RSEM output is usually named '<sample>.genes.results'.",
            call. = FALSE
        )
    }

    if (is.null(sample_names)) {
        sample_names <- sub(pattern, "", basename(files))
    } else if (length(sample_names) != length(files)) {
        stop(
            "[load_rsem_as_tximport] 'sample_names' has ", length(sample_names),
            " entries but ", length(files), " files were matched.",
            call. = FALSE
        )
    }
    dupes <- sample_names[duplicated(sample_names)]
    if (length(dupes) > 0) {
        stop("[load_rsem_as_tximport] Sample names are not unique: ",
             paste(unique(dupes), collapse = ", "), call. = FALSE)
    }
    names(files) <- sample_names

    # Isoform-level input needs a tx2gene map; gene-level (the default) does not.
    txi <- tximport::tximport(
        files,
        type = "rsem",
        txIn = !is.null(tx2gene),
        txOut = FALSE,
        tx2gene = tx2gene
    )

    if (fix_zero_length) {
        # RSEM reports effective_length 0 for features with no reads; DESeq2's
        # tximport path rejects a zero-length offset. Setting them to 1 keeps
        # the offset matrix finite without affecting the counts.
        txi$length[txi$length == 0] <- 1
    }

    validate_tximport(txi)
    message(sprintf(
        "[load_rsem_as_tximport] Imported %d samples x %d genes from RSEM",
        ncol(txi$counts), nrow(txi$counts)
    ))
    txi
}


# load_omics_inputs and validate_contrasts_content live in R/core/01_io.R
