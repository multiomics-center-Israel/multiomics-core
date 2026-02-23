# 00a_deseq_factory.R
# DESeq2 Dataset Factory with tximport Support
#
# This module provides a unified factory for creating DESeqDataSet objects
# from either raw count matrices or tximport objects. It enforces strict
# validation and prevents silent data corruption.

# =============================================================================
# Input Type Detection and Validation
# =============================================================================

#' Detect source type from expression input
#'
#' Determines whether input is a tximport object or raw count matrix.
#' Detection priority: annotation > structure > class
#'
#' @param expr Expression input (matrix, data.frame, or tximport list)
#' @return Character: "tximport" or "matrix"
#' @keywords internal
detect_source_type <- function(expr) {
    # Priority 1: Check for explicit annotation
    annotated_type <- attr(expr, "source_type")
    if (!is.null(annotated_type)) {
        # Validate annotation against structure
        structural_type <- infer_source_type_from_structure(expr)
        if (annotated_type != structural_type) {
            stop(
                sprintf(
                    "[source_type] Annotation '%s' does not match object structure (detected: '%s'). ",
                    annotated_type, structural_type
                ),
                "This indicates a pipeline configuration error. ",
                "Verify that the input object matches the declared source_type.",
                call. = FALSE
            )
        }
        message(sprintf("[DESeq2] Input source_type: %s (annotation validated)", annotated_type))
        return(annotated_type)
    }

    # Priority 2-4: Infer from structure
    inferred_type <- infer_source_type_from_structure(expr)
    message(sprintf("[DESeq2] Input source_type: %s (structure inferred)", inferred_type))
    return(inferred_type)
}

#' Infer source type from object structure
#'
#' @param expr Expression input
#' @return Character: "tximport" or "matrix"
#' @keywords internal
infer_source_type_from_structure <- function(expr) {
    # Check for tximport signature: list with counts, abundance, length matrices
    if (is.list(expr) && !is.data.frame(expr)) {
        if (is_valid_tximport_structure(expr, validate_only = TRUE)) {
            return("tximport")
        }
    }

    # Check for matrix/data.frame (raw counts)
    if (is.matrix(expr) || is.data.frame(expr)) {
        return("matrix")
    }

    stop(
        "[source_type] Unable to determine input type. ",
        "Expected either: (1) tximport list with 'counts', 'abundance', 'length' matrices, ",
        "or (2) numeric matrix/data.frame of counts.",
        call. = FALSE
    )
}

#' Validate tximport object structure
#'
#' Checks that a tximport object has the required fields with matching dimensions.
#'
#' @param txi Suspected tximport object
#' @param validate_only If TRUE, returns FALSE instead of stopping on invalid structure
#' @return TRUE if valid, FALSE or stops if invalid
#' @keywords internal
is_valid_tximport_structure <- function(txi, validate_only = FALSE) {
    required_fields <- c("counts", "abundance", "length")

    # Check required fields exist
    missing_fields <- setdiff(required_fields, names(txi))
    if (length(missing_fields) > 0) {
        if (validate_only) return(FALSE)
        stop(
            "[tximport] Missing required fields: ", paste(missing_fields, collapse = ", "),
            ". A valid tximport object must contain 'counts', 'abundance', and 'length' matrices.",
            call. = FALSE
        )
    }

    # Check all required fields are matrices
    for (field in required_fields) {
        if (!is.matrix(txi[[field]])) {
            if (validate_only) return(FALSE)
            stop(
                sprintf("[tximport] Field '%s' must be a matrix, got %s.", field, class(txi[[field]])[1]),
                call. = FALSE
            )
        }
    }

    # Check dimensions match
    dims <- lapply(txi[required_fields], dim)
    if (!all(sapply(dims, function(d) identical(d, dims[[1]])))) {
        if (validate_only) return(FALSE)
        dim_str <- sapply(required_fields, function(f) paste(dim(txi[[f]]), collapse = "x"))
        stop(
            "[tximport] Dimension mismatch between matrices. ",
            paste(required_fields, dim_str, sep = ": ", collapse = ", "),
            call. = FALSE
        )
    }

    # Check column names match (sample IDs must be consistent)
    colnames_list <- lapply(txi[required_fields], colnames)
    if (!all(sapply(colnames_list, function(cn) identical(cn, colnames_list[[1]])))) {
        if (validate_only) return(FALSE)
        stop(
            "[tximport] Column names (sample IDs) must be identical across 'counts', 'abundance', and 'length' matrices.",
            call. = FALSE
        )
    }

    # Check rownames match (gene IDs must be consistent)
    rownames_list <- lapply(txi[required_fields], rownames)
    if (!all(sapply(rownames_list, function(rn) identical(rn, rownames_list[[1]])))) {
        if (validate_only) return(FALSE)
        stop(
            "[tximport] Row names (gene IDs) must be identical across 'counts', 'abundance', and 'length' matrices.",
            call. = FALSE
        )
    }

    TRUE
}

#' Validate tximport object (strict mode)
#'
#' Performs full validation of a tximport object. Terminates on failure.
#'
#' @param txi tximport object
#' @return TRUE invisibly if valid
#' @keywords internal
validate_tximport <- function(txi) {
    is_valid_tximport_structure(txi, validate_only = FALSE)

    # Log countsFromAbundance for reproducibility (optional field)
    cfa <- txi$countsFromAbundance
    if (!is.null(cfa)) {
        message(sprintf("[tximport] countsFromAbundance: %s", cfa))
    } else {
        message("[tximport] countsFromAbundance: not set (using tximport defaults)")
    }

    invisible(TRUE)
}

#' Validate raw count matrix (strict integer enforcement)
#'
#' Ensures counts are strictly integer. Non-integer values terminate with error.
#'
#' @param counts Count matrix
#' @return TRUE invisibly if valid
#' @keywords internal
validate_count_matrix <- function(counts) {
    if (!is.matrix(counts) && !is.data.frame(counts)) {
        stop("[counts] Input must be a matrix or data.frame.", call. = FALSE)
    }

    counts_mat <- as.matrix(counts)

    # Check for negative values
    if (any(counts_mat < 0, na.rm = TRUE)) {
        stop("[counts] Count matrix contains negative values, which is invalid.", call. = FALSE)
    }

    # Strict integer check
    if (!is.integer(counts_mat)) {
        # Check if values are integer-like (whole numbers stored as numeric)
        if (is.numeric(counts_mat)) {
            non_integer_mask <- (counts_mat != floor(counts_mat)) & !is.na(counts_mat)
            if (any(non_integer_mask)) {
                # Find examples of non-integer values
                non_int_idx <- which(non_integer_mask, arr.ind = TRUE)[1:min(3, sum(non_integer_mask)), , drop = FALSE]
                examples <- sapply(1:nrow(non_int_idx), function(i) {
                    sprintf("  [%d,%d] = %g", non_int_idx[i, 1], non_int_idx[i, 2],
                            counts_mat[non_int_idx[i, 1], non_int_idx[i, 2]])
                })

                stop(
                    "[counts] Non-integer counts detected in matrix-type input.\n",
                    "Examples:\n", paste(examples, collapse = "\n"), "\n\n",
                    "This suggests tximport-derived data that should use the tximport path.\n",
                    "If this is raw count data, ensure integer values. Do NOT round or coerce.",
                    call. = FALSE
                )
            }
        }
    }

    # Validate row/column names
    if (is.null(rownames(counts_mat)) || any(rownames(counts_mat) == "")) {
        stop("[counts] Count matrix must have non-empty rownames (gene IDs).", call. = FALSE)
    }
    if (is.null(colnames(counts_mat)) || any(colnames(counts_mat) == "")) {
        stop("[counts] Count matrix must have non-empty colnames (sample IDs).", call. = FALSE)
    }

    invisible(TRUE)
}

# =============================================================================
# Sample Alignment
# =============================================================================

#' Align samples between expression data and metadata (strict mode by default)
#'
#' @param expr_samples Character vector of sample IDs from expression data
#' @param meta Data frame with sample metadata
#' @param sample_col Column name containing sample IDs in metadata
#' @param lenient If TRUE, subset to intersection with warning; if FALSE (default), fail on mismatch
#' @return Data frame: metadata subset and reordered to match expr_samples
#' @keywords internal
align_samples_strict <- function(expr_samples, meta, sample_col, lenient = TRUE) {
    if (!sample_col %in% colnames(meta)) {
        stop(sprintf("[alignment] Metadata missing sample column: '%s'", sample_col), call. = FALSE)
    }

    meta_samples <- as.character(meta[[sample_col]])

    # Find discrepancies
    in_expr_not_meta <- setdiff(expr_samples, meta_samples)
    in_meta_not_expr <- setdiff(meta_samples, expr_samples)

    has_mismatch <- length(in_expr_not_meta) > 0 || length(in_meta_not_expr) > 0

    if (has_mismatch) {
        # Build detailed error message
        msg_parts <- c("[alignment] Sample mismatch detected:")

        if (length(in_expr_not_meta) > 0) {
            msg_parts <- c(msg_parts, sprintf(
                "  - In expression data but not metadata (%d): %s",
                length(in_expr_not_meta),
                paste(head(in_expr_not_meta, 5), collapse = ", "),
                if (length(in_expr_not_meta) > 5) "..." else ""
            ))
        }

        if (length(in_meta_not_expr) > 0) {
            msg_parts <- c(msg_parts, sprintf(
                "  - In metadata but not expression data (%d): %s",
                length(in_meta_not_expr),
                paste(head(in_meta_not_expr, 5), collapse = ", "),
                if (length(in_meta_not_expr) > 5) "..." else ""
            ))
        }

        if (!lenient) {
            stop(
                paste(msg_parts, collapse = "\n"), "\n\n",
                "To allow partial overlap, set sample_alignment: lenient in config.",
                call. = FALSE
            )
        } else {
            warning(paste(msg_parts, collapse = "\n"))
        }
    }

    # Compute intersection
    common_samples <- intersect(expr_samples, meta_samples)
    if (length(common_samples) == 0) {
        stop("[alignment] No samples in common between expression data and metadata.", call. = FALSE)
    }

    # Subset and reorder metadata to match expression sample order
    # (preserving the order from expr_samples for downstream consistency)
    meta_aligned <- meta[match(common_samples, meta_samples), , drop = FALSE]
    rownames(meta_aligned) <- common_samples

    message(sprintf("[alignment] %d samples matched between expression data and metadata", length(common_samples)))

    meta_aligned
}

# =============================================================================
# DESeqDataSet Factory
# =============================================================================
#' Create DESeqDataSet from polymorphic input
#'
#' Factory function that creates a DESeqDataSet using the appropriate constructor
#' based on input type (raw counts or tximport). Enforces strict validation and
#' preserves tximport length offsets.
#'
#' @param expr Expression input: either a count matrix or tximport object
#' @param meta Data frame with sample metadata
#' @param design Formula for the DESeq2 design (e.g., ~ condition)
#' @param sample_col Column name in meta containing sample IDs
#' @param lenient_alignment If TRUE, allow sample subsetting with warning; if FALSE (default), fail on mismatch
#' @return DESeqDataSet object
#' @export
create_deseq_dataset <- function(expr, meta, design, sample_col, lenient_alignment = FALSE) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required.", call. = FALSE)
  }
  
  # Detect and validate source type
  source_type <- detect_source_type(expr)
  
  if (source_type == "tximport") {
    # Validate tximport structure
    validate_tximport(expr)
    
    # Get sample IDs from tximport
    expr_samples <- colnames(expr$counts)
    
    # Align samples (strict by default)
    meta_aligned <- align_samples_strict(expr_samples, meta, sample_col, lenient = lenient_alignment)
    
    # Subset tximport to aligned samples if needed
    aligned_samples <- rownames(meta_aligned)
    if (!identical(expr_samples, aligned_samples)) {
      expr <- subset_tximport(expr, samples = aligned_samples)
    }
    
    # --- OPTION A: Fix zero/negative lengths for DESeq2 robustness ---
    if (any(expr$length <= 0)) {
      n_fix <- sum(expr$length <= 0)
      message(sprintf("[tximport] Fixing %d cells with length <= 0 (setting to 1) to prevent DESeq2 crash.", n_fix))
      expr$length[expr$length <= 0] <- 1
    }
    # -----------------------------------------------------------------
    
    # Create DESeqDataSet via tximport path
    message("[DESeq2] Constructing dataset via DESeqDataSetFromTximport")
    dds <- DESeq2::DESeqDataSetFromTximport(
      txi = expr,
      colData = as.data.frame(meta_aligned),
      design = design
    )
    
  } else {
    # source_type == "matrix"
    # Validate count matrix (strict integer enforcement)
    validate_count_matrix(expr)
    
    # Convert to matrix if data.frame
    counts_mat <- as.matrix(expr)
    storage.mode(counts_mat) <- "integer"
    
    # Get sample IDs
    expr_samples <- colnames(counts_mat)
    
    # Align samples (strict by default)
    meta_aligned <- align_samples_strict(expr_samples, meta, sample_col, lenient = lenient_alignment)
    
    # Subset counts to aligned samples if needed
    aligned_samples <- rownames(meta_aligned)
    if (!identical(expr_samples, aligned_samples)) {
      counts_mat <- counts_mat[, aligned_samples, drop = FALSE]
    }
    
    # Create DESeqDataSet via matrix path
    message("[DESeq2] Constructing dataset via DESeqDataSetFromMatrix")
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts_mat,
      colData = as.data.frame(meta_aligned),
      design = design
    )
  }
  
  # Attach source_type to DDS for downstream reference
  S4Vectors::metadata(dds)$source_type <- source_type
  
  return(dds)
}

# =============================================================================
# tximport Subsetting Utilities
# =============================================================================

#' Subset tximport object by samples
#'
#' Subsets all three matrices (counts, abundance, length) together.
#'
#' @param txi tximport object
#' @param samples Character vector of sample IDs to keep
#' @return Subsetted tximport object
#' @keywords internal
subset_tximport <- function(txi, samples) {
    available_samples <- colnames(txi$counts)
    missing <- setdiff(samples, available_samples)
    if (length(missing) > 0) {
        stop(
            "[tximport] Cannot subset: samples not found in tximport object: ",
            paste(head(missing, 5), collapse = ", "),
            call. = FALSE
        )
    }

    # Subset all three matrices together (invariant: always subset together)
    txi$counts <- txi$counts[, samples, drop = FALSE]
    txi$abundance <- txi$abundance[, samples, drop = FALSE]
    txi$length <- txi$length[, samples, drop = FALSE]
    
    # Subset inferential replicates if they exist
    if (!is.null(txi$infReps)) {
      txi$infReps <- lapply(txi$infReps, function(m) m[, samples, drop = FALSE])
    }
    
    # Preserve other attributes
    txi
}

#' Subset tximport object by genes
#'
#' Subsets all three matrices (counts, abundance, length) by gene IDs.
#' Used for gene filtering.
#'
#' @param txi tximport object
#' @param genes Character vector of gene IDs to keep, OR logical vector (keep_vec)
#' @return Subsetted tximport object
#' @keywords internal
subset_tximport_genes <- function(txi, genes) {
  if (is.logical(genes)) {
    # genes is a keep_vec
    if (length(genes) != nrow(txi$counts)) {
      stop(
        sprintf("[tximport] keep_vec length (%d) != number of genes (%d)",
                length(genes), nrow(txi$counts)),
        call. = FALSE
      )
    }
    keep_idx <- genes
  } else {
    # genes is a character vector of gene IDs
    available_genes <- rownames(txi$counts)
    missing <- setdiff(genes, available_genes)
    if (length(missing) > 0) {
      warning(
        sprintf("[tximport] %d genes not found in tximport object (ignored)", length(missing))
      )
    }
    keep_idx <- available_genes %in% genes
  }
  
  # Subset the three core matrices
  txi$counts    <- txi$counts[keep_idx, , drop = FALSE]
  txi$abundance <- txi$abundance[keep_idx, , drop = FALSE]
  txi$length    <- txi$length[keep_idx, , drop = FALSE]
  
  # Handle inferential replicates if they exist
  if (!is.null(txi$infReps)) {
    txi$infReps <- lapply(txi$infReps, function(m) m[keep_idx, , drop = FALSE])
  }
  
  txi
}

# =============================================================================
# Annotation Utilities
# =============================================================================

#' Annotate expression object with source_type
#'
#' Sets the source_type attribute on an expression object.
#'
#' @param expr Expression object (matrix or tximport)
#' @param source_type Character: "matrix" or "tximport"
#' @return expr with source_type attribute set
#' @export
annotate_source_type <- function(expr, source_type) {
    if (!source_type %in% c("matrix", "tximport")) {
        stop("source_type must be 'matrix' or 'tximport'", call. = FALSE)
    }
    attr(expr, "source_type") <- source_type
    expr
}

#' Get source_type from expression object or DESeqDataSet
#'
#' @param obj Expression object or DESeqDataSet
#' @return Character: "matrix" or "tximport", or NULL if not annotated
#' @export
get_source_type <- function(obj) {
    # Check for attribute (on raw expression objects)
    st <- attr(obj, "source_type")
    if (!is.null(st)) return(st)

    # Check for DESeqDataSet metadata
    if (inherits(obj, "DESeqDataSet")) {
        st <- S4Vectors::metadata(obj)$source_type
        if (!is.null(st)) return(st)
    }

    NULL
}
