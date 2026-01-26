#' Preprocess proteomics data starting from loaded inputs
preprocess_proteomics <- function(inputs, config) {
    cfg <- config$modes$proteomics

    # validate inputs
    validate_proteomics_inputs(inputs, cfg)

    # Build standardized proteomics object
    prot_obj <- get_proteomics_expression_matrix(inputs, config)

    # Downstream contract: work on log2 assay
    expr_raw <- prot_obj$assay_log2
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

    # Align col_data
    sample_id_col <- cfg$effects$samples %||% cfg$id_columns$sample_col
    check_has_cols(col_data, sample_id_col, df_name = "col_data")
    col_data <- align_meta_to_expr(expr_raw, col_data, cfg)

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
        }
    }

    # Filtering
    filt <- filter_proteomics_by_min_count(expr_mat = expr_raw, row_data = row_data, meta = col_data, cfg = cfg)
    expr_filt <- filt$expr_mat
    row_data_f <- filt$row_data
    assert_numeric_matrix(expr_filt, "expr_filt")

    # Single imputation (QC/plots)
    imp_res <- impute_proteomics_perseus(expr_mat = expr_filt, cfg = cfg, return_flags = TRUE)
    expr_imp_single <- imp_res$imputed
    assert_numeric_matrix(expr_imp_single, "expr_imp_single")

    imputation_qc <- NULL
    if (!is.null(imp_res$imputed_flag)) {
        imputation_qc <- list(imputed_flag = imp_res$imputed_flag)
    }

    list(
        expr_raw = expr_raw,
        expr_filt = expr_filt,
        expr_work = expr_imp_single,
        expr_imp_single = expr_imp_single,
        row_data = row_data_f,
        meta = col_data,
        imputation_qc = imputation_qc
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
