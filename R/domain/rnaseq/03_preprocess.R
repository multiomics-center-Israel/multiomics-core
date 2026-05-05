# 03_preprocess.R
#
preprocess_rna <- function(inputs, config, gene_lengths = NULL, verbose = FALSE) {
  cfg <- config$modes$rna
  sample_col <- cfg$id_columns$sample_col %||% "SampleID"
  lenient_alignment <- identical(cfg$sample_alignment, "lenient")
  
  # =========================================================================
  # Detect input type: preprocessed vs tximport vs raw counts
  # =========================================================================
  is_preprocessed <- identical(cfg$input$format, "preprocessed")
  has_txi <- !is.null(inputs$txi) && is_valid_tximport_structure(inputs$txi, validate_only = TRUE)
  has_counts <- !is.null(inputs$counts)
  
  if (is_preprocessed && has_counts) {
    source_type <- "preprocessed"
    if(verbose) message("[preprocess_rna] Input detected as preprocessed expression matrix")
    
    gene_id_col <- cfg$id_columns$gene_id
    if (is.null(gene_id_col) || !nzchar(gene_id_col)) {
      stop("[preprocess_rna] 'id_columns$gene_id' must be specified in config for preprocessed input.", call. = FALSE)
    }
    
    all_cols <- names(inputs$counts)
    non_id_cols <- setdiff(all_cols, gene_id_col)
    is_numeric <- vapply(inputs$counts[, non_id_cols, drop = FALSE], 
                         function(x) is.numeric(x) || is.integer(x), logical(1))
    
    sample_cols <- non_id_cols[is_numeric]
    anno_cols <- non_id_cols[!is_numeric]
    
    if (length(sample_cols) == 0) stop("Preprocessed table has no numeric sample columns.")
    
    row_data <- inputs$counts[, c(gene_id_col, anno_cols), drop = FALSE]
    counts <- as.matrix(inputs$counts[, sample_cols, drop = FALSE])
    storage.mode(counts) <- "numeric"
    
    # Handle log-transform (from your oz_multi logic)
    if (!isTRUE(cfg$files$is_logtransformed)) {
      message("[preprocess_rna] Applying log2(x + 1) to preprocessed counts")
      counts <- log2(counts + 1)
    }
    
    # Aggregating duplicate gene IDs (averaging for preprocessed, not summing)
    gene_ids <- inputs$counts[[gene_id_col]]
    if (anyDuplicated(gene_ids)) {
      counts_sum <- rowsum(counts, group = gene_ids, reorder = FALSE)
      group_sizes <- as.numeric(table(gene_ids)[rownames(counts_sum)])
      counts <- counts_sum / group_sizes
      row_data <- row_data[!duplicated(gene_ids), , drop = FALSE]
      gene_ids <- gene_ids[!duplicated(gene_ids)]
    }
    rownames(counts) <- gene_ids
    txi <- NULL; abundance <- NULL

  } else if (has_txi) {
    source_type <- "tximport"
    if(verbose) message("[preprocess_rna] Input detected as tximport object")
    validate_tximport(inputs$txi)
    txi <- inputs$txi; counts <- txi$counts; abundance <- txi$abundance
    gene_ids <- rownames(counts)
    row_data <- data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)

  } else if (has_counts) {
    source_type <- "matrix"
    if(verbose) message("[preprocess_rna] Input detected as raw count matrix")
    validate_rna_inputs(inputs, cfg)
    
    gene_id_col <- cfg$id_columns$gene_id
    all_cols <- names(inputs$counts)
    sample_cols <- all_cols[vapply(inputs$counts, is.numeric, logical(1))]
    sample_cols <- setdiff(sample_cols, gene_id_col)
    
    counts <- as.matrix(inputs$counts[, sample_cols, drop = FALSE])
    gene_ids <- inputs$counts[[gene_id_col]]
    if (anyDuplicated(gene_ids)) {
      counts <- rowsum(counts, group = gene_ids, reorder = FALSE)
      gene_ids <- rownames(counts)
    }
    row_data <- data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)
    rownames(counts) <- gene_ids
    validate_count_matrix(counts)
    txi <- NULL; abundance <- NULL
  }
  
  # =========================================================================
  # Optional sample name mapping (applies to both input types)
  # =========================================================================
  map_from <- cfg$id_columns$map_from
  map_to   <- cfg$id_columns$map_to
  if (!is.null(inputs$sample_map) && !is.null(map_from) && !is.null(map_to)) {
    counts <- apply_sample_map_to_colnames(counts, inputs$sample_map, map_from, map_to)

    # Also update txi matrices if present
    if (!is.null(txi)) {
      colnames(txi$counts)    <- colnames(counts)
      colnames(txi$abundance) <- colnames(counts)
      colnames(txi$length)    <- colnames(counts)
      abundance <- txi$abundance
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
  # Sample filtering (rule-based, if configured)
  # =========================================================================
  rules <- get_sample_filter_rules(config, mode = "rna")
  if (!is.null(rules) && exists("apply_sample_filter")) {
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
      
      if (!is.null(txi)) {
        txi <- subset_tximport(txi, samples = colnames(counts))
        abundance <- txi$abundance
      }
  }
  
  # =========================================================================
  # Gene filtering
  # =========================================================================
  filter_mode <- cfg$filtering$mode %||% "adaptive"
  group_col <- cfg$filtering$group_col %||% cfg$de_table$group_col %||% cfg$effects$color[1]
  

  if (filter_mode == "none" || source_type == "preprocessed") {
    message("Filtering mode: none/preprocessed — keeping all features.")
    keep_vec <- rep(TRUE, nrow(counts))
    fr <- list(keep_vec = keep_vec, used_threshold = NA)
    norm_mode <- if (source_type == "preprocessed") "preprocessed" else "none"
  } else {
    # Determine evaluation matrix for filtering
    if (source_type == "tximport") {
      norm_for_filter <- abundance
      norm_mode <- "TPM (tximport)"
    } else {
      use_tpm <- !is.null(gene_lengths) && all(c("gene", "length") %in% colnames(gene_lengths))
      norm_for_filter <- if (use_tpm) compute_tpm(counts, gene_lengths) else compute_cpm(counts)
      norm_mode <- if (use_tpm) "TPM" else "CPM"
    }
    
    if (filter_mode == "deseq2_only") {
      keep_vec <- rowSums(counts) > 0
      fr <- list(keep_vec = keep_vec, used_threshold = NA)
    } else if (filter_mode == "fixed") {
      fixed_thr <- as.numeric(cfg$filtering$fixed_threshold %||% 1.0)
      fr <- filter_features_optimized(norm_for_filter, meta2, sample_col, group_col, threshold = fixed_thr)
      fr$used_threshold <- fixed_thr
    } else {
      plot_path <- file.path(config$paths$out, "rnaseq", "filtering_threshold_qc.png")
      fr <- run_auto_filter_pipeline(norm_for_filter, meta2, sample_col, group_col, output_plot = plot_path)
    }
  }

  # FIX: Extract keep_vec and threshold from result list
  keep_vec <- fr$keep_vec
  thr <- fr$used_threshold
  if (sum(keep_vec) == 0) stop("Filtering removed all features.")

  # Apply gene filter
  counts_filt <- counts[keep_vec, , drop = FALSE]
  row_data_filt <- row_data[keep_vec, , drop = FALSE]
  txi_filt <- if (!is.null(txi)) subset_tximport_genes(txi, genes = keep_vec) else NULL

  # =========================================================================
  # Normalized expression for QC/visualization (expr_work)
  # =========================================================================
  if (source_type == "preprocessed") {
      expr_work <- counts_filt
      attr(expr_work, "method") <- "preprocessed"
      message("[Normalization] Skipped — using preprocessed values directly")
  } else {
      norm_input <- if (!is.null(txi_filt)) {
          annotate_source_type(txi_filt, "tximport")
      } else {
          annotate_source_type(counts_filt, "matrix")
      }

      expr_work <- normalize_counts(
          counts            = norm_input,
          meta              = meta2,
          method            = cfg$normalization$method %||% "TMMlogCPM",
          prior.count       = as.numeric(cfg$normalization$prior.count %||% 1),
          sample_col        = sample_col,
          filter_zero_count = cfg$de$filter_zero_count
      )
  }

  # FIX: Ensure row consistency if filtering was skipped or modified
  if (!identical(rownames(counts_filt), rownames(expr_work))) {
    counts_filt <- counts_filt[rownames(expr_work), , drop = FALSE]
    row_data_filt <- row_data_filt[rownames(expr_work), , drop = FALSE]
  }

  # =========================================================================
  # Build return object
  # =========================================================================
  result <- list(
    expr_raw = counts,
    expr_filt = counts_filt,
    expr_work = expr_work,
    
    # de_input for downstream analysis
    de_input = if (source_type == "preprocessed") {
        annotate_source_type(counts_filt, "preprocessed")
    } else if (!is.null(txi_filt)) {
        annotate_source_type(txi_filt, "tximport")
    } else {
        annotate_source_type(counts_filt, "matrix")
    },
    
    row_data = row_data_filt,
    meta = meta2,
    
    info = list(
      source_type = source_type,
      filter_norm = norm_mode,
      filter_threshold = thr,
      expr_work_method = attr(expr_work, "method"),
      n_genes_raw = nrow(counts),
      n_genes_filt = nrow(counts_filt),
      n_samples = ncol(counts_filt)
    ),
    contrasts = inputs$contrasts
  )
  
  attr(result, "source_type") <- source_type
  return(result)
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
