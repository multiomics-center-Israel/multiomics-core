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

#' Resolve the normalization method that actually produced expr_work
#'
#' The one place downstream code asks "how was this run normalized?". It reads
#' the method recorded on the matrix by normalize_counts() (or set to
#' "preprocessed" by preprocess_rna()), not \code{normalization$method} from the
#' config: the config is what was requested, while the attribute is what ran.
#' The two differ when the config omits the key, and when VST on a count matrix
#' falls back to TMMlogCPM.
#'
#' @param pre List returned by preprocess_rna().
#' @return Character scalar, e.g. "VST", "TMMlogCPM", "preprocessed", or
#'   "unknown" when neither the matrix nor \code{pre$info} records a method.
resolve_rna_norm_method <- function(pre) {
    m <- attr(pre$expr_work, "method") %||% pre$info$expr_work_method
    if (is.null(m) || length(m) != 1 || is.na(m) || !nzchar(m)) return("unknown")
    as.character(m)
}

#' Locate this run's normalized-expression export and the method it records
#'
#' write_rnaseq_datasets_legacy() names the file after the method that ran
#' (\code{rna_norm_<method>.tsv}) and removes any other \code{rna_norm_*.tsv},
#' so a single match is this run's matrix. Zero or several matches return NULL
#' rather than guessing which file is current.
#'
#' @param datasets_dir Path to the run's \code{Datasets/} directory.
#' @return \code{list(path, method)}, or NULL when there is not exactly one
#'   \code{rna_norm_*.tsv} file.
find_rna_norm_file <- function(datasets_dir) {
    hits <- list.files(datasets_dir, pattern = "^rna_norm_.+\\.tsv$", full.names = TRUE)
    if (length(hits) != 1) return(NULL)
    list(path = hits, method = sub("^rna_norm_(.+)\\.tsv$", "\\1", basename(hits)))
}

#' Reader-facing wording for an RNA normalization method
#'
#' Shared by the report, pipeline summary, PowerPoint and commentary so they
#' describe the same method the same way. Unrecognised methods get neutral
#' wording instead of borrowing the TMMlogCPM description.
#'
#' @param method Method string, normally from resolve_rna_norm_method() or
#'   find_rna_norm_file().
#' @param prior_count Prior count used for logCPM; only quoted for TMMlogCPM.
#' @return List with \code{label} (short name), \code{name} (name plus tool),
#'   \code{value_label} (axis / column wording), \code{sentence} (one sentence
#'   for Methods-style prose) and \code{tags} (short chips).
describe_rna_norm_method <- function(method, prior_count = 1) {
    method <- if (is.null(method) || length(method) != 1 || is.na(method) ||
                  !nzchar(method)) "unknown" else as.character(method)
    switch(method,
        TMMlogCPM = list(
            label       = "TMMlogCPM",
            name        = "TMM + log2 CPM (edgeR)",
            value_label = "log2 CPM (TMM)",
            sentence    = sprintf(paste0(
                "Expression values were normalized with edgeR TMM (Trimmed Mean of M-values) ",
                "scaling factors followed by log2 counts-per-million (prior count = %s)."),
                prior_count),
            tags        = c("TMM scaling", "log2 CPM", sprintf("prior.count = %s", prior_count))
        ),
        VST = list(
            label       = "VST",
            name        = "VST (DESeq2)",
            value_label = "VST (log2 scale)",
            sentence    = paste0(
                "Expression values were transformed with the DESeq2 variance-stabilizing ",
                "transformation (VST, blind to the experimental design)."),
            tags        = c("VST", "DESeq2", "blind = TRUE")
        ),
        preprocessed = list(
            label       = "preprocessed",
            name        = "none (preprocessed input used as supplied)",
            value_label = "Expression (as supplied)",
            sentence    = paste0(
                "The input was supplied already normalized; the pipeline used the values ",
                "as provided and applied no further normalization."),
            tags        = c("preprocessed input", "no normalization")
        ),
        list(
            label       = method,
            name        = method,
            value_label = "Normalized expression",
            sentence    = if (identical(method, "unknown")) {
                "The normalization method could not be determined from the run outputs."
            } else {
                sprintf("Expression values were normalized using %s.", method)
            },
            tags        = method
        )
    )
}
