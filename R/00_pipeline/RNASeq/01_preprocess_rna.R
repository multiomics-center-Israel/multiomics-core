# 01_preprocess_rna
preprocess_rna <- function(inputs, config, gene_lengths = NULL, verbose = FALSE) {
  cfg <- config$modes$rna
  
  gene_id_col <- cfg$id_columns$gene_id
  sample_col  <- cfg$id_columns$sample_col %||% "SampleID"
  
  if (!gene_id_col %in% names(inputs$counts)) {
    stop("RNA counts missing gene id column: ", gene_id_col)
  }
  
  row_data <- inputs$counts[, gene_id_col, drop = FALSE]
  sample_cols <- setdiff(names(inputs$counts), gene_id_col)
  
  counts <- as.matrix(inputs$counts[, sample_cols, drop = FALSE])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- inputs$counts[[gene_id_col]]
  
  # optional mapping
  map_from <- cfg$id_columns$map_from
  map_to   <- cfg$id_columns$map_to
  if (!is.null(inputs$sample_map) && all(c(map_from, map_to) %in% names(inputs$sample_map))) {
    m <- setNames(inputs$sample_map[[map_to]], inputs$sample_map[[map_from]])
    mapped <- m[colnames(counts)]
    if (sum(!is.na(mapped)) > 0) colnames(counts)[!is.na(mapped)] <- unname(mapped[!is.na(mapped)])
  }
  
  meta <- inputs$metadata
  if (!sample_col %in% names(meta)) {
    stop("metadata missing sample column: ", sample_col)
  }
  
  # --- rule: metadata samples must be subset of counts samples ---
  meta_samples   <- unique(as.character(meta[[sample_col]]))
  count_samples  <- colnames(counts)
  
  missing_in_counts <- setdiff(meta_samples, count_samples)
  if (length(missing_in_counts) > 0) {
    stop(
      "metadata contains samples not present in counts (not allowed): ",
      paste(missing_in_counts, collapse = ", ")
    )
  }
  
  # --- allowed: counts has extra samples not in metadata -> warn + subset counts ---
  extra_in_counts <- setdiff(count_samples, meta_samples)
  if (length(extra_in_counts) > 0) {
    warning(sprintf(
      "Counts contains %d samples but metadata contains %d. Dropping %d count samples not in metadata.",
      length(count_samples), length(meta_samples), length(extra_in_counts)
    ))
    counts <- counts[, meta_samples, drop = FALSE]
    count_samples <- colnames(counts)
  }
  
  # --- now reorder metadata to match counts ---
  meta2 <- meta[match(count_samples, meta[[sample_col]]), , drop = FALSE]
  
  # safety: should never happen now
  if (anyNA(meta2[[sample_col]])) {
    stop("Internal error: failed to align metadata to counts after subset.")
  }
  
  rownames(meta2) <- meta2[[sample_col]]
  
  rules <- get_sample_filter_rules(config, mode = "rna")
  if (!is.null(rules)) {
    filtered <- apply_sample_filter(
      meta  = meta2,
      expr  = counts,
      rules = rules,
      mode  = "rna",
      strict_cols = FALSE
    )
    meta2  <- filtered$meta
    counts <- filtered$expr
  }
  
  # ---- dynamic filter (by Line) ----
  group_col <- cfg$filtering$group_col %||% "Line"
  if (!group_col %in% names(meta2)) stop("metadata missing filter group column: ", group_col)
  
  thr <- as.numeric(cfg$filtering$threshold %||% cfg$filter_tpm %||% 3)
  
  use_tpm <- !is.null(gene_lengths) &&
    all(c("gene", "length") %in% colnames(gene_lengths)) &&
    length(intersect(rownames(counts), gene_lengths$gene)) > 0
  
  norm_for_filter <- if (use_tpm) compute_tpm(counts, gene_lengths) else compute_cpm(counts)
  norm_mode <- if (use_tpm) "TPM" else "CPM"
  
  fr <- filter_features_dynamic(
    norm_mat   = norm_for_filter,
    meta       = meta2,
    sample_col = sample_col,
    group_col  = group_col,
    threshold  = thr
  )
  
  if (sum(fr$keep_vec) == 0) stop("Filtering removed all features. Check threshold/grouping.")
  
  counts_filt   <- counts[fr$keep_vec, , drop = FALSE]
  row_data_filt <- row_data[fr$keep_vec, , drop = FALSE]
  
  # ---- expr_work for QC/DE ----
  ncfg <- cfg$normalization
  expr_work <- normalize_counts(
    counts      = counts_filt,
    meta        = meta2,
    method      = ncfg$method %||% "TMMlogCPM",
    prior.count = as.numeric(ncfg$prior.count %||% 1)
  )
  
  if (verbose) {
    message(sprintf(
      "RNA filter=%s thr=%s by=%s kept=%d/%d; expr_work=%s",
      norm_mode, thr, group_col, sum(fr$keep_vec), length(fr$keep_vec),
      attr(expr_work, "method") %||% (ncfg$method %||% "TMMlogCPM")
    ))
  }
  
  list(
    expr_raw  = counts,
    expr_filt = counts_filt,
    expr_work = expr_work,
    row_data  = row_data_filt,
    meta      = meta2,
    qc        = list(),
    info      = list(
      filter_norm = norm_mode,
      filter_threshold = thr,
      filter_group = group_col,
      expr_work_method = attr(expr_work, "method")
    ),
    contrasts = inputs$contrasts
  )
}
