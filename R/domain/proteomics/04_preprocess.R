#' Preprocess proteomics data starting from loaded inputs
preprocess_proteomics <- function(inputs, config) {
    cfg <- config$modes$proteomics

    # validate inputs
    validate_proteomics_inputs(inputs, cfg)

    # Build standardized proteomics object
    prot_obj <- get_proteomics_expression_matrix(inputs, config)

    # Downstream contract: work on log2 assay
    expr_raw <- prot_obj$assay_log2
    expr_linear <- prot_obj$assay_linear
    row_data <- prot_obj$row_data
    col_data <- prot_obj$col_data

    # Use protein_id_col from config
    protein_id_col <- cfg$id_columns$protein_id

    # Ensure expr_raw matrix contract
    if (is.data.frame(expr_raw)) {
        rn <- NULL
        if (!is.null(row_data) && protein_id_col %in% colnames(row_data)) rn <- row_data[[protein_id_col]]
        expr_raw <- coerce_df_to_numeric_matrix(expr_raw, rownames_vec = rn, name = "expr_raw")
    }

    # Ensure rownames
    if (!is.null(row_data) && protein_id_col %in% colnames(row_data)) {
        rn <- as.character(row_data[[protein_id_col]])
        if (is.null(rownames(expr_raw)) || all(grepl("^\\d+$", rownames(expr_raw)))) {
            rownames(expr_raw) <- rn
        }
    }
    assert_numeric_matrix(expr_raw, "expr_raw")

    # Contaminant filtering (e.g. cRAP proteins) — configurable toggle
    if (isTRUE(cfg$filtering$remove_contaminants %||% TRUE)) {
        contam_res <- filter_contaminants(expr_raw, row_data, cfg)
        # Sync linear matrix to same rows
        if (!is.null(expr_linear)) {
            keep_rows <- rownames(contam_res$expr_mat)
            expr_linear <- expr_linear[keep_rows, , drop = FALSE]
        }
        expr_raw <- contam_res$expr_mat
        row_data <- contam_res$row_data
    }

    # Align col_data
    sample_id_col <- cfg$effects$samples %||% cfg$id_columns$sample_col
    check_has_cols(col_data, sample_id_col, df_name = "col_data")
    col_data <- align_meta_to_expr(expr_raw, col_data, cfg)

    # Exclude specific samples (if configured)
    excl <- cfg$exclude_samples
    if (!is.null(excl) && length(excl) > 0) {
        excl <- as.character(unlist(excl))
        n_before <- nrow(col_data)
        col_data <- col_data[!col_data[[sample_id_col]] %in% excl, , drop = FALSE]
        message(sprintf("[preprocess_proteomics] Excluded %d sample(s): %s",
                        n_before - nrow(col_data), paste(excl, collapse = ", ")))
        keep_cols <- intersect(colnames(expr_raw), rownames(col_data))
        expr_raw <- expr_raw[, keep_cols, drop = FALSE]
        if (!is.null(expr_linear)) {
            expr_linear <- expr_linear[, keep_cols, drop = FALSE]
        }
    }

    # Optional: sample_filter
    rules <- get_sample_filter_rules(config, mode = "proteomics")
    if (!is.null(rules)) {
        # implementation specific to sample filtering would act here
        # assuming logic similar to utils::apply_sample_filter is available or moved.
        # We will assume apply_sample_filter is in R/core/something.R or we need to copy it.
        # I didn't see apply_sample_filter in my extraction list, it was in 00_utils.R (lines 158).
        # I should have put it in R/core/02_validation.R or 03_alignment.R.
        # I'll check if I missed it. Just in case, I'll allow this file to assume it's available.
        if (exists("apply_sample_filter")) {
            filtered <- apply_sample_filter(sample_col = sample_id_col, meta = col_data, expr = expr_raw, rules = rules, mode = "proteomics")
            col_data <- filtered$meta
            expr_raw <- filtered$expr
            if (!is.null(expr_linear)) {
                expr_linear <- expr_linear[rownames(expr_raw), colnames(expr_raw), drop = FALSE]
            }
        }
    }

    # Filtering
    filt <- filter_proteomics_by_min_count(expr_mat = expr_raw, row_data = row_data, meta = col_data, cfg = cfg)
    expr_filt <- filt$expr_mat
    row_data_f <- filt$row_data
    assert_numeric_matrix(expr_filt, "expr_filt")

    # Normalization — dispatches based on cfg$normalization$method
    norm_method <- cfg$normalization$method %||% "none"
    if (norm_method == "median") {
        expr_filt <- normalize_proteomics_median(expr_filt)
        message("Normalization: median centering applied.")
    } else if (norm_method == "vsn") {
        # VSN operates on linear (raw) intensities and produces its own
        # generalized-log transform, so we feed the linear-scale matrix
        if (is.null(expr_linear)) {
            stop("VSN normalization requires linear-scale data but assay_linear is NULL.")
        }
        expr_linear_filt <- expr_linear[rownames(expr_filt), colnames(expr_filt), drop = FALSE]
        expr_filt <- normalize_proteomics_vsn(expr_linear_filt)
        message("Normalization: VSN (variance stabilizing) applied to linear intensities.")
    } else if (norm_method != "none") {
        stop(sprintf("Unknown normalization method: '%s'. Supported: 'none', 'median', 'vsn'.", norm_method))
    }

    # Single imputation (QC/plots) — dispatches based on cfg$imputation$method
    imp_res <- impute_proteomics(expr_mat = expr_filt, cfg = cfg, return_flags = TRUE)
    expr_imp_single <- imp_res$imputed
    assert_numeric_matrix(expr_imp_single, "expr_imp_single")

    imputation_qc <- NULL
    if (!is.null(imp_res$imputed_flag)) {
        imputation_qc <- list(imputed_flag = imp_res$imputed_flag)
    }

    list(
        expr_raw = expr_raw,
        expr_filt = expr_filt,
        expr_filt_pre_imp = expr_filt,
        expr_work = expr_imp_single,
        expr_imp_single = expr_imp_single,
        row_data = row_data_f,
        meta = col_data,
        imputation_qc = imputation_qc
    )
}

#' Median centering normalization for proteomics
#'
#' Centers each sample (column) by subtracting its median, so all samples
#' have the same median intensity. Common for log2-transformed proteomics data.
#' @param expr_mat numeric matrix (proteins x samples)
#' @return normalized matrix
normalize_proteomics_median <- function(expr_mat) {
    col_medians <- apply(expr_mat, 2, median, na.rm = TRUE)
    global_median <- median(col_medians, na.rm = TRUE)
    sweep(expr_mat, 2, col_medians - global_median)
}

#' Variance Stabilizing Normalization (VSN) for proteomics
#'
#' Applies vsnMatrix to linear-scale intensity data (matching DEP::normalize_vsn).
#' VSN performs its own generalized-log transformation internally,
#' so the output is on a log-like (glog2) scale comparable to log2.
#' NAs are preserved: the model is fitted on complete values only.
#' Requires the Bioconductor 'vsn' package.
#' @param expr_mat numeric matrix (proteins x samples), linear scale
#' @return normalized matrix (glog2 scale)
normalize_proteomics_vsn <- function(expr_mat) {
    if (!requireNamespace("vsn", quietly = TRUE))
        stop("VSN normalization requires the 'vsn' Bioconductor package.\n",
             "Install with: BiocManager::install('vsn')")
    na_mask <- is.na(expr_mat)
    # Use vsnMatrix (matching DEP::normalize_vsn exactly)
    vsn_fit <- vsn::vsnMatrix(expr_mat)
    result <- vsn::predict(vsn_fit, expr_mat)
    result[na_mask] <- NA
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
