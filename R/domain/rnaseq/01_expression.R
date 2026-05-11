#' Normalize RNA-seq counts
#'
#' Supports both raw count matrices and tximport objects. For tximport inputs,
#' VST is performed on a DESeqDataSet constructed via DESeqDataSetFromTximport,
#' preserving length offset information.
#'
#' Options:
#' - "TMMlogCPM" – edgeR TMM normalization + logCPM (matrix input only)
#' - "VST"       – DESeq2 variance stabilizing transform (matrix or tximport)
#'
#' @param counts Expression input: matrix/data.frame of raw counts (genes x samples)
#'   OR tximport object with 'counts', 'abundance', 'length' matrices.
#' @param meta   data.frame of sample metadata (required for VST)
#' @param method one of c("TMMlogCPM","VST")
#' @param prior.count numeric added before log in logCPM (default 1)
#' @param sample_col Column name in meta containing sample IDs (default "SampleID")
#' @return numeric matrix; attr(., "method") indicates method used
normalize_counts <- function(counts, meta = NULL, method = c("TMMlogCPM", "VST"),
                             prior.count = 1, sample_col = "SampleID", filter_zero_count = TRUE) {
    method <- match.arg(method)
   message("THE FILTER ZERO COUNT IS:  ", filter_zero_count )
    # Detect source type
    source_type <- detect_source_type(counts)

    if (method == "TMMlogCPM") {
        # TMMlogCPM only supports matrix input
        if (source_type == "tximport") {
            stop(
                "[normalize_counts] TMMlogCPM requires raw count matrix input. ",
                "For tximport data, use method = 'VST' instead.",
                call. = FALSE
            )
        }

        if (!requireNamespace("edgeR", quietly = TRUE)) {
            stop("Package 'edgeR' is required.")
        }

        counts_mat <- as.matrix(counts)
        dge <- edgeR::DGEList(counts = counts_mat)
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        mat <- edgeR::cpm(dge, log = TRUE, prior.count = prior.count)
        attr(mat, "method") <- "TMMlogCPM"
        attr(mat, "source_type") <- "matrix"
        return(mat)
    }

    # VST branch - supports both matrix and tximport
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Package 'DESeq2' is required.")
    }
    if (is.null(meta)) {
        stop("For method = 'VST', 'meta' must be provided.")
    }

    # Ensure sample_col exists in meta
    if (!sample_col %in% colnames(meta)) {
        if (!is.null(rownames(meta)) && !any(rownames(meta) == "")) {
            meta[[sample_col]] <- rownames(meta)
        } else {
            stop(sprintf(
                "Sample column '%s' not found in metadata and rownames not usable",
                sample_col
            ))
        }
    }

    # Create DESeqDataSet using factory (handles both matrix and tximport)
    # This ensures tximport data uses DESeqDataSetFromTximport with length offsets
    dds <- create_deseq_dataset(
        expr = counts,
        meta = meta,
        design = ~1,
        sample_col = sample_col,
        lenient_alignment = FALSE
    )

    # Filter zero-count genes
    if (isTRUE(filter_zero_count)) {
      keep_nonzero <- rowSums(DESeq2::counts(dds)) > 0
      n_dropped <- sum(!keep_nonzero)
      if (n_dropped > 0) {
        message(sprintf("Dropped %d zero-count features before VST", n_dropped))
      }
      dds <- dds[keep_nonzero, , drop = FALSE]
      if (nrow(dds) == 0) {
        stop("No rows with nonzero counts available for VST.")
      }
    }
    # Store original counts for potential fallback
    original_counts <- if (source_type == "matrix") as.matrix(counts) else NULL

    mat <- tryCatch(
        {
              vt <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
              SummarizedExperiment::assay(vt)
            
        },
        error = function(e1) {
            # Fallback to TMMlogCPM only if source is matrix (not tximport)
            if (source_type == "matrix" && !is.null(original_counts)) {
                message("[VST] Fallback to TMMlogCPM due to: ", conditionMessage(e1))
                normalize_counts(original_counts, meta, "TMMlogCPM", prior.count, sample_col)
            } else {
                stop(
                    "[VST] Failed for tximport input: ", conditionMessage(e1), "\n",
                    "Fallback to TMMlogCPM is not available for tximport data.",
                    call. = FALSE
                )
            }
        }
    )

    if (is.null(attr(mat, "method"))) attr(mat, "method") <- "VST"
    attr(mat, "source_type") <- source_type
    mat
}

# compute CPM
compute_cpm <- function(counts) {
    counts <- as.matrix(counts)
    lib <- colSums(counts, na.rm = TRUE)
    if (any(lib == 0 | !is.finite(lib))) stop("Zero/invalid library size for CPM.")
    sweep(counts, 2, lib, "/") * 1e6
}

# compute TPM if gene_length exist
compute_tpm <- function(counts, gene_lengths_bp) {
    stopifnot(!is.null(rownames(counts)))
    stopifnot(all(c("gene", "length") %in% colnames(gene_lengths_bp)))

    common <- intersect(rownames(counts), gene_lengths_bp$gene)
    if (length(common) == 0) stop("No overlapping genes for TPM.")

    counts <- as.matrix(counts[common, , drop = FALSE])
    gl_kb <- gene_lengths_bp$length[match(common, gene_lengths_bp$gene)] / 1000
    if (any(!is.finite(gl_kb) | gl_kb <= 0)) stop("Invalid gene lengths.")

    rpk <- counts / gl_kb
    sf <- colSums(rpk, na.rm = TRUE)
    if (any(sf == 0 | !is.finite(sf))) stop("Zero/invalid TPM scaling factor.")

    tpm <- sweep(rpk, 2, sf, "/") * 1e6
    rownames(tpm) <- common
    tpm
}
