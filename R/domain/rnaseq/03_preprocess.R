# 03_preprocess.R
preprocess_rna <- function(inputs, config, gene_lengths = NULL, verbose = FALSE) {
    cfg <- config$modes$rna
    # Validate
    validate_rna_inputs(inputs, cfg)

    gene_id_col <- cfg$id_columns$gene_id
    sample_col <- cfg$id_columns$sample_col %||% "SampleID"

    row_data <- inputs$counts[, gene_id_col, drop = FALSE]
    sample_cols <- setdiff(names(inputs$counts), gene_id_col)

    counts <- as.matrix(inputs$counts[, sample_cols, drop = FALSE])
    storage.mode(counts) <- "numeric"
    rownames(counts) <- inputs$counts[[gene_id_col]]

    # optional mapping
    map_from <- cfg$id_columns$map_from
    map_to <- cfg$id_columns$map_to
    if (!is.null(inputs$sample_map) && !is.null(map_from) && !is.null(map_to) &&
        all(c(map_from, map_to) %in% names(inputs$sample_map))) {
        m <- setNames(inputs$sample_map[[map_to]], inputs$sample_map[[map_from]])
        mapped <- m[colnames(counts)]
        if (sum(!is.na(mapped)) > 0) colnames(counts)[!is.na(mapped)] <- unname(mapped[!is.na(mapped)])
    }

    meta <- inputs$metadata
    meta_samples <- unique(as.character(meta[[sample_col]]))
    count_samples <- colnames(counts)

    # Extra samples drop
    extra_in_counts <- setdiff(count_samples, meta_samples)
    if (length(extra_in_counts) > 0) {
        counts <- counts[, meta_samples, drop = FALSE]
        count_samples <- colnames(counts)
    }

    meta2 <- meta[match(count_samples, meta[[sample_col]]), , drop = FALSE]
    rownames(meta2) <- meta2[[sample_col]]

    # Sample filtering
    rules <- get_sample_filter_rules(config, mode = "rna")
    if (!is.null(rules)) {
        if (exists("apply_sample_filter")) {
            filtered <- apply_sample_filter(
                sample_col = sample_col, meta = meta2, expr = counts, rules = rules, mode = "rna", strict_cols = FALSE
            )
            meta2 <- filtered$meta
            counts <- filtered$expr
        }
    }

    # Dynamic feature filtering
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

    if (sum(fr$keep_vec) == 0) stop("Filtering removed all features.")

    counts_filt <- counts[fr$keep_vec, , drop = FALSE]
    row_data_filt <- row_data[fr$keep_vec, , drop = FALSE]

    # expr_work (for QC/DE)
    ncfg <- cfg$normalization
    expr_work <- normalize_counts(
        counts      = counts_filt,
        meta        = meta2,
        method      = ncfg$method %||% "TMMlogCPM",
        prior.count = as.numeric(ncfg$prior.count %||% 1)
    )

    list(
        expr_raw = counts,
        expr_filt = counts_filt,
        expr_work = expr_work,
        row_data = row_data_filt,
        meta = meta2,
        qc = list(),
        info = list(
            filter_norm = norm_mode,
            filter_threshold = thr,
            filter_group = group_col,
            expr_work_method = attr(expr_work, "method")
        ),
        contrasts = inputs$contrasts
    )
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
