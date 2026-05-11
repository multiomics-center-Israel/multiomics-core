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


# load_omics_inputs and validate_contrasts_content live in R/core/01_io.R
