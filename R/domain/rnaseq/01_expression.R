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

#' Linear TMM-normalized CPM
#'
#' The linear counterpart of the TMMlogCPM branch of \code{normalize_counts()}:
#' the same \code{DGEList} and TMM factors, so the same effective library sizes,
#' but no log and no prior count. Unlike \code{compute_cpm()}, library sizes
#' are TMM-scaled.
#'
#' @param counts Raw count matrix (genes x samples), the same matrix that
#'   \code{normalize_counts()} received.
#' @return Numeric matrix of CPM values with the dimnames of \code{counts}.
#' @examples
#' counts <- matrix(c(10, 0, 5, 20, 3, 8), nrow = 3,
#'                  dimnames = list(c("G1", "G2", "G3"), c("S1", "S2")))
#' compute_tmm_cpm(counts)
compute_tmm_cpm <- function(counts) {
    dge <- edgeR::DGEList(counts = as.matrix(counts))
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    edgeR::cpm(dge, log = FALSE)
}

#' Resolve the transformation used for RNA-seq QC plots
#'
#' The one reading of \code{qc$cpm_pseudocount}, shared by the QC module, the
#' fact sheet, the report and the pipeline summary. It is decided by the method
#' that actually produced \code{expr_work} (\code{pre$info$expr_work_method},
#' which is "TMMlogCPM" also when VST fell back to it), not by the configured
#' method. Only TMMlogCPM has linear TMM CPM behind it; VST and preprocessed
#' input are plotted as they are.
#'
#' @param info The \code{info} list returned by \code{preprocess_rna()}.
#' @param qc_cfg The \code{modes$rna$qc} config list (may be NULL).
#' @return A list with \code{applied} (logical), \code{pseudocount} (numeric,
#'   NA when not applied) and \code{label} (reader-facing description).
#' @examples
#' resolve_rna_qc_transform(list(expr_work_method = "TMMlogCPM"), list(cpm_pseudocount = 1))
resolve_rna_qc_transform <- function(info, qc_cfg = NULL) {
    method <- info$expr_work_method %||% NA_character_
    if (identical(method, "TMMlogCPM")) {
        pseudocount <- as.numeric(qc_cfg$cpm_pseudocount %||% 1)
        return(list(applied = TRUE, pseudocount = pseudocount,
                    label = sprintf("log2(TMM CPM + %s)", format(pseudocount))))
    }
    label <- if (identical(method, "VST")) {
        "the VST values of the analysis matrix"
    } else if (identical(method, "preprocessed")) {
        "the preprocessed input values"
    } else {
        "the analysis matrix"
    }
    list(applied = FALSE, pseudocount = NA_real_, label = label)
}

#' Expression matrix for RNA-seq QC plots
#'
#' \code{log2(TMM CPM + qc$cpm_pseudocount)} when \code{expr_work} is
#' TMMlogCPM, otherwise \code{expr_work} itself. edgeR's \code{prior.count} is
#' in reads and, at typical depths, a far smaller offset than 1 CPM, so genes
#' with a handful of reads dominate the variance of \code{expr_work}; a CPM
#' offset damps them for sample-level QC without changing the matrix that
#' clustering, the post-DE heatmaps and the exported normalized table use.
#'
#' @param pre The list returned by \code{preprocess_rna()}.
#' @param qc_cfg The \code{modes$rna$qc} config list (may be NULL).
#' @return Numeric matrix with the dimnames of \code{pre$expr_work}.
rna_qc_expr <- function(pre, qc_cfg = NULL) {
    tr <- resolve_rna_qc_transform(pre$info, qc_cfg)
    if (!tr$applied) return(pre$expr_work)
    # expr_filt is the count matrix normalize_counts() turned into expr_work,
    # so the same genes, samples and TMM factors are used.
    if (!identical(rownames(pre$expr_filt), rownames(pre$expr_work)) ||
        !identical(colnames(pre$expr_filt), colnames(pre$expr_work))) {
        stop("[rnaseq QC] expr_filt and expr_work have different genes or samples, ",
             "so the QC matrix cannot be built from the counts behind expr_work. ",
             "preprocess_rna() should return them aligned.", call. = FALSE)
    }
    log2(compute_tmm_cpm(pre$expr_filt) + tr$pseudocount)
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
