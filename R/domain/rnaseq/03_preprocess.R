# 03_preprocess.R
#
# Supports both raw count matrices and tximport objects.
# Input type is detected automatically and handled appropriately.

preprocess_rna <- function(inputs, config, gene_lengths = NULL, verbose = FALSE) {
    cfg <- config$modes$rna
    sample_col <- cfg$id_columns$sample_col %||% "SampleID"

    # Determine alignment mode from config (strict by default)
    lenient_alignment <- identical(cfg$sample_alignment, "lenient")

    # =========================================================================
    # Detect input type: tximport vs raw counts
    # =========================================================================
    has_txi <- !is.null(inputs$txi) && is_valid_tximport_structure(inputs$txi, validate_only = TRUE)
    has_counts <- !is.null(inputs$counts)

    if (has_txi) {
        source_type <- "tximport"
        message("[preprocess_rna] Input detected as tximport object")

        # Validate tximport structure (will stop on failure)
        validate_tximport(inputs$txi)

        # Extract expression data from tximport
        txi <- inputs$txi
        counts <- txi$counts
        abundance <- txi$abundance  # TPM for filtering
        message(sprintf("[txi_counts] %d", dim(counts)[1]))
        # Gene IDs are rownames
        gene_ids <- rownames(counts)
        row_data <- data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)

    } else if (has_counts) {
        source_type <- "matrix"
        message("[preprocess_rna] Input detected as raw count matrix")

        # Validate using existing function
        validate_rna_inputs(inputs, cfg)

        gene_id_col <- cfg$id_columns$gene_id
        if (is.null(gene_id_col) || !nzchar(gene_id_col)) {
            stop(
                "[preprocess_rna] 'id_columns$gene_id' must be specified in config ",
                "when using raw count matrix input. ",
                "This column identifies which column in the counts table contains gene IDs.",
                call. = FALSE
            )
        }
        if (!gene_id_col %in% names(inputs$counts)) {
            stop(
                "[preprocess_rna] Gene ID column '", gene_id_col, "' not found in counts table. ",
                "Available columns: ", paste(names(inputs$counts), collapse = ", "),
                call. = FALSE
            )
        }
        # Identify annotation columns vs sample columns:
        # Keep only columns that are numeric (counts); exclude gene_id and any other
        # annotation columns (e.g. gene_position, transcript_id).
        all_cols <- names(inputs$counts)
        non_id_cols <- setdiff(all_cols, gene_id_col)
        is_numeric <- vapply(inputs$counts[, non_id_cols, drop = FALSE],
                             function(x) is.numeric(x) || is.integer(x),
                             logical(1))
        anno_cols <- non_id_cols[!is_numeric]
        sample_cols <- non_id_cols[is_numeric]

        if (length(sample_cols) == 0) stop("Counts table has no numeric sample columns.")

        row_data <- inputs$counts[, c(gene_id_col, anno_cols), drop = FALSE]

        counts <- as.matrix(inputs$counts[, sample_cols, drop = FALSE])
        storage.mode(counts) <- "numeric"

        # Aggregate duplicate gene IDs (multiple transcripts per gene) by summing counts
        gene_ids <- inputs$counts[[gene_id_col]]
        if (anyDuplicated(gene_ids)) {
            if (verbose) message("Aggregating ", length(gene_ids), " transcript rows to ",
                                 length(unique(gene_ids)), " unique gene IDs (summing counts).")
            counts <- rowsum(counts, group = gene_ids, reorder = FALSE)
            # Keep first annotation row per gene for row_data
            first_idx <- !duplicated(gene_ids)
            row_data <- row_data[first_idx, , drop = FALSE]
            gene_ids <- gene_ids[first_idx]
        }
        rownames(counts) <- gene_ids

        # Validate strict integer counts
        validate_count_matrix(counts)

        txi <- NULL
        abundance <- NULL

    } else {
        stop(
            "[preprocess_rna] No valid expression input found. ",
            "Expected either 'txi' (tximport object) or 'counts' (data.frame) in inputs.",
            call. = FALSE
        )
    }

    # =========================================================================
    # Optional sample name mapping (applies to both input types)
    # =========================================================================
    map_from <- cfg$id_columns$map_from
    map_to <- cfg$id_columns$map_to
    if (!is.null(inputs$sample_map) && !is.null(map_from) && !is.null(map_to) &&
        all(c(map_from, map_to) %in% names(inputs$sample_map))) {
        m <- setNames(inputs$sample_map[[map_to]], inputs$sample_map[[map_from]])
        mapped <- m[colnames(counts)]
        if (sum(!is.na(mapped)) > 0) {
            old_names <- colnames(counts)
            colnames(counts)[!is.na(mapped)] <- unname(mapped[!is.na(mapped)])
            message(sprintf(
                "[preprocess_rna] Mapped %d sample names",
                sum(!is.na(mapped))
            ))

            # Also update txi matrices if present
            if (!is.null(txi)) {
                colnames(txi$counts) <- colnames(counts)
                colnames(txi$abundance) <- colnames(counts)
                colnames(txi$length) <- colnames(counts)
                abundance <- txi$abundance
            }
        }
    }

    # =========================================================================
    # Sample alignment (strict by default)
    # =========================================================================
    meta <- inputs$metadata
    expr_samples <- colnames(counts)

    meta2 <- align_samples_strict(
        expr_samples = expr_samples,
        meta = meta,
        sample_col = sample_col,
        lenient = lenient_alignment
    )

    # Subset expression data to aligned samples
    aligned_samples <- rownames(meta2)
    if (!identical(expr_samples, aligned_samples)) {
        counts <- counts[, aligned_samples, drop = FALSE]
        if (!is.null(txi)) {
            txi <- subset_tximport(txi, samples = aligned_samples)
            abundance <- txi$abundance
        }
    }

    # =========================================================================
    # Exclude specific samples (if configured)
    # =========================================================================
    excl <- cfg$exclude_samples
    if (!is.null(excl) && length(excl) > 0) {
        excl <- as.character(unlist(excl))
        n_before <- nrow(meta2)
        meta2 <- meta2[!meta2[[sample_col]] %in% excl, , drop = FALSE]
        message(sprintf("[preprocess_rna] Excluded %d sample(s): %s",
                        n_before - nrow(meta2), paste(excl, collapse = ", ")))
        counts <- counts[, meta2[[sample_col]], drop = FALSE]
        if (!is.null(txi)) {
            txi <- subset_tximport(txi, samples = meta2[[sample_col]])
            abundance <- txi$abundance
        }
    }

    # =========================================================================
    # Sample filtering (rule-based, if configured)
    # =========================================================================
    rules <- get_sample_filter_rules(config, mode = "rna")
    if (!is.null(rules)) {
        if (exists("apply_sample_filter")) {
            filtered <- apply_sample_filter(
                sample_col = sample_col,
                meta = meta2,
                expr = counts,
                rules = rules,
                mode = "rna",
                strict_cols = FALSE
            )
            meta2 <- filtered$meta
            counts <- filtered$expr

            # Also subset txi if present
            if (!is.null(txi)) {
                txi <- subset_tximport(txi, samples = colnames(counts))
                abundance <- txi$abundance
            }
        }
    }

    # =========================================================================
    # Gene filtering
    # =========================================================================
    group_col <- cfg$filtering$group_col %||% "Line"
    if (!group_col %in% names(meta2)) {
        stop("metadata missing filter group column: ", group_col)
    }

    # Determine what to use for filtering threshold evaluation
    if (source_type == "tximport") {
        # For tximport: prefer txi$abundance (TPM) for filtering
        norm_for_filter <- abundance
        norm_mode <- "TPM (tximport)"
        message("[Filtering] Threshold evaluation using txi$abundance (TPM)")
    } else {
        # For raw counts: compute CPM or TPM if gene lengths available
        use_tpm <- !is.null(gene_lengths) &&
            all(c("gene", "length") %in% colnames(gene_lengths)) &&
            length(intersect(rownames(counts), gene_lengths$gene)) > 0

        norm_for_filter <- if (use_tpm) {
            compute_tpm(counts, gene_lengths)
        } else {
            compute_cpm(counts)
        }
        norm_mode <- if (use_tpm) "TPM" else "CPM"
        message(sprintf("[Filtering] Threshold evaluation using %s", norm_mode))
    }

    filter_mode <- cfg$filtering$mode %||% "adaptive"
    filt_flag <- cfg$filtering$enable %||% TRUE

    if (!filt_flag) {
        # Filtering disabled entirely
        message("Filtering disabled by config.")
        keep_vec <- rep(TRUE, nrow(norm_for_filter))
        thr <- 0

    } else if (filter_mode == "deseq2_only") {
        # Only remove all-zero genes; let DESeq2 handle the rest
        message("Filtering mode: deseq2_only \u2014 removing all-zero genes only.")
        keep_vec <- rowSums(counts) > 0
        thr <- NA

    } else if (filter_mode == "fixed") {
        # Classic fixed CPM threshold
        fixed_thr <- as.numeric(cfg$filtering$fixed_threshold %||% 1.0)
        message(sprintf("Filtering mode: fixed CPM threshold = %.2f", fixed_thr))
        fr <- filter_features_optimized(
            norm_mat    = norm_for_filter,
            meta        = meta2,
            sample_col  = sample_col,
            group_col   = group_col,
            threshold   = fixed_thr
        )
        keep_vec <- fr$keep_vec
        thr <- fixed_thr

    } else {
        # Default: adaptive KDE
        plot_path <- file.path(config$paths$out, "rnaseq", "filtering_threshold_qc.png")
        fr <- run_auto_filter_pipeline(
            cpm_mat     = norm_for_filter,
            meta        = meta2,
            sample_col  = sample_col,
            group_col   = group_col,
            output_plot = plot_path
        )
        keep_vec <- fr$keep_vec
        thr <- fr$used_threshold
    }

    if (sum(keep_vec) == 0) {
        stop("Filtering removed all features.")
    }

    # Apply gene filter to counts
    counts_filt <- counts[keep_vec, , drop = FALSE]
    row_data_filt <- row_data[keep_vec, , drop = FALSE]

    # Apply gene filter to txi (all three matrices together - invariant)
    txi_filt <- NULL
    if (!is.null(txi)) {
        txi_filt <- subset_tximport_genes(txi, genes = keep_vec)
    }

    # =========================================================================
    # Normalized expression for QC/visualization (expr_work)
    # =========================================================================
    ncfg <- cfg$normalization

    # Determine input for normalization
    norm_input <- if (!is.null(txi_filt)) {
        annotate_source_type(txi_filt, "tximport")
    } else {
        annotate_source_type(counts_filt, "matrix")
    }

    expr_work <- normalize_counts(
        counts      = norm_input,
        meta        = meta2,
        method      = ncfg$method %||% "TMMlogCPM",
        prior.count = as.numeric(ncfg$prior.count %||% 1),
        sample_col  = sample_col
    )
    
    # match norm expr and filt expr in case filtering step wasn't apply
    if (!filt_flag) counts_filt = counts_filt[rownames(expr_work), ]
    
    # =========================================================================
    # Build return object
    # =========================================================================
    result <- list(
        # Expression data (counts-like, NOT necessarily integer for tximport)
        expr_raw = counts,
        expr_filt = counts_filt,
        expr_work = expr_work,

        # For DE: pass the appropriate object based on source type
        # This is what run_deseq2_de() should receive
        de_input = if (!is.null(txi_filt)) {
            annotate_source_type(txi_filt, "tximport")
        } else {
            annotate_source_type(counts_filt, "matrix")
        },

        # Metadata
        row_data = row_data_filt,
        meta = meta2,
        qc = list(),

        # Pipeline info
        info = list(
            source_type = source_type,
            filter_norm = norm_mode,
            filter_threshold = thr,
            filter_group = group_col,
            expr_work_method = attr(expr_work, "method"),
            n_genes_raw = nrow(counts),
            n_genes_filt = nrow(counts_filt),
            n_samples = ncol(counts_filt)
        ),

        contrasts = inputs$contrasts
    )

    # Attach source_type to result for downstream checks
    attr(result, "source_type") <- source_type

    result
}

get_sample_filter_rules <- function(config, mode) {
    cfg <- config$modes[[mode]]
    if (is.null(cfg)) {
        return(NULL)
    }
    sf <- cfg$sample_filter
    if (is.null(sf) || !isTRUE(sf$enabled)) {
        return(NULL)
    }
    sf$rules %||% NULL
}
