# R/domain/lipidomics/02_preprocess.R
#
# Preprocessing assembly for lipidomics.
# Dispatches to format-specific parsers, builds the standard internal
# representation, and applies normalization.
#
# Normalization reuses the metabolomics normalization functions:
#   normalize_samples, transform_metab, scale_metab, apply_normalization_pipeline
#
# Reuses: assert_numeric_matrix, assert_meta_contract, align_meta_to_matrix,
#         build_minimal_meta


#' Preprocess lipidomics data
#'
#' @param inputs  List from load_lipidomics_inputs().
#' @param config  Full pipeline config.
#' @return list matching the pre-processing contract:
#'   expr_raw, expr_filt, expr_work, meta, row_data, info,
#'   normalization_eval (NULL if evaluate_methods not set)
preprocess_lipidomics <- function(inputs, config) {
    cfg <- config$modes$lipidomics
    format <- inputs$format %||% cfg$input$format %||% "lipidsearch"

    # ---- 1. Parse input into expr_raw / row_data / sample_ids ----
    parsed <- switch(format,
        lipidsearch     = parse_lipidsearch(inputs$data, cfg, inputs$metadata),
        translated_csv  = parse_translated_csv(inputs$data, cfg, inputs$metadata),
        processed_wide  = parse_lipid_processed_wide(inputs$data, cfg, inputs$metadata),
        stop("Unsupported lipidomics input format: ", format)
    )

    expr_raw <- parsed$expr_raw
    row_data <- parsed$row_data
    sample_ids <- colnames(expr_raw)

    assert_numeric_matrix(expr_raw, "lipid_expr_raw")

    # ---- 2. Build / align metadata ----
    meta <- inputs$metadata
    sample_col <- cfg$effects$samples %||% "sample_id"

    if (is.null(meta)) {
        meta <- build_minimal_meta(sample_ids)
        sample_col <- "sample_id"

        if (is.null(cfg$effects$color)) cfg$effects$color <- "sample_type"
        if (is.null(cfg$effects$samples)) cfg$effects$samples <- "sample_id"
    } else {
        check_has_cols(meta, sample_col, df_name = "metadata")
        meta <- align_meta_to_matrix(sample_ids, meta, sample_col)
    }

    assert_meta_contract(meta, sample_col)

    # ---- 3. Optional sample filtering (QC/blank/pool removal) ----
    rules <- get_sample_filter_rules_lipid(cfg)
    if (!is.null(rules)) {
        keep <- apply_sample_filter_lipid(sample_ids, meta, rules, sample_col)
        expr_raw <- expr_raw[, keep, drop = FALSE]
        meta <- meta[meta[[sample_col]] %in% keep, , drop = FALSE]
        sample_ids <- keep
    }

    # ---- 4. Basic feature filtering: remove all-NA/all-zero rows ----
    row_valid <- apply(expr_raw, 1, function(x) {
        any(!is.na(x) & x != 0)
    })
    expr_filt <- expr_raw[row_valid, , drop = FALSE]
    row_data <- row_data[row_valid, , drop = FALSE]

    # ---- 4b. Ensure lipid class annotation ----
    if (!"lipid_class" %in% colnames(row_data)) {
        if ("Name" %in% colnames(row_data)) {
            row_data$lipid_class <- parse_lipid_class(row_data$Name)
        } else {
            row_data$lipid_class <- parse_lipid_class(row_data$feature_id)
        }
    }

    # ---- 4c. Parse chain info if not present ----
    if (!"total_carbons" %in% colnames(row_data)) {
        name_src <- row_data$Name %||% row_data$feature_id
        chain_info <- parse_lipid_chains(name_src)
        row_data$total_carbons      <- chain_info$total_carbons
        row_data$total_double_bonds <- chain_info$total_double_bonds
    }

    # ---- 4d. Parse bond type if not present ----
    if (!"bond_type" %in% colnames(row_data)) {
        bt_src <- row_data$Name %||% row_data$feature_id
        row_data$bond_type <- detect_lipid_bond_type(bt_src)
    }

    # Also ensure a Name column exists for downstream display
    if (!"Name" %in% colnames(row_data)) {
        row_data$Name <- row_data$feature_id
    }

    norm_cfg_raw <- cfg$normalization %||% list()

    # ---- 5. Handle missing values ----
    na_count <- sum(is.na(expr_filt))
    na_policy <- norm_cfg_raw$na_policy %||% "keep"
    na_policy <- tolower(na_policy)

    expr_for_norm <- expr_filt
    if (na_policy == "zero") {
        if (na_count > 0) {
            message(sprintf("lipidomics: replacing %d NA values with 0 (na_policy='zero').",
                            na_count))
        }
        expr_for_norm[is.na(expr_for_norm)] <- 0
    } else if (na_policy == "lod") {
        lod_frac <- norm_cfg_raw$na_lod_fraction %||% 0.2
        n_replaced <- 0
        for (i in seq_len(nrow(expr_for_norm))) {
            row_vals <- expr_for_norm[i, ]
            is_missing <- is.na(row_vals) | row_vals == 0
            if (any(is_missing)) {
                pos_vals <- row_vals[!is_missing & row_vals > 0]
                fill_val <- if (length(pos_vals) > 0) min(pos_vals) * lod_frac else NA_real_
                expr_for_norm[i, is_missing] <- fill_val
                n_replaced <- n_replaced + sum(is_missing)
            }
        }
        if (n_replaced > 0) {
            message(sprintf("lipidomics: imputed %d zero/NA values with %.1f%% of min per feature (na_policy='lod').",
                            n_replaced, lod_frac * 100))
        }
    } else if (na_policy == "min_half") {
        n_replaced <- 0
        for (i in seq_len(nrow(expr_for_norm))) {
            row_vals <- expr_for_norm[i, ]
            is_missing <- is.na(row_vals) | row_vals == 0
            if (any(is_missing)) {
                pos_vals <- row_vals[!is_missing & row_vals > 0]
                fill_val <- if (length(pos_vals) > 0) min(pos_vals) / 2 else NA_real_
                expr_for_norm[i, is_missing] <- fill_val
                n_replaced <- n_replaced + sum(is_missing)
            }
        }
        if (n_replaced > 0) {
            message(sprintf("lipidomics: imputed %d zero/NA values with min/2 per feature (na_policy='min_half').",
                            n_replaced))
        }
    } else if (na_count > 0) {
        message(sprintf("lipidomics: keeping %d NA values as-is (na_policy='keep').",
                        na_count))
    }

    # ---- 6. Normalization (reuses metabolomics normalization functions) ----
    norm_cfg <- norm_cfg_raw
    input_already_normalized <- isTRUE(norm_cfg$input_already_normalized)

    norm_cfg$sample_norm <- norm_cfg$sample_norm %||% "none"
    norm_cfg$transform   <- norm_cfg$transform   %||% "none"
    norm_cfg$scaling     <- norm_cfg$scaling     %||% "none"

    # Pre-scaling matrix for logFC computation
    norm_cfg_no_scale <- norm_cfg
    norm_cfg_no_scale$scaling <- "none"
    pre_scale_result <- apply_normalization_pipeline(expr_for_norm, norm_cfg_no_scale, row_data)
    expr_log <- pre_scale_result$expr_norm

    # Full pipeline (with scaling)
    norm_result <- apply_normalization_pipeline(expr_for_norm, norm_cfg, row_data)
    expr_work <- norm_result$expr_norm

    assert_numeric_matrix(expr_work, "lipid_expr_work")

    # ---- 7. Optional: evaluate alternative methods ----
    norm_eval <- NULL
    eval_methods <- norm_cfg$evaluate_methods
    if (!is.null(eval_methods) && length(eval_methods) > 0) {
        if (input_already_normalized) {
            for (j in seq_along(eval_methods)) {
                if (is.null(eval_methods[[j]]$sample_norm)) {
                    eval_methods[[j]]$sample_norm <- "none"
                }
            }
        }
        norm_eval <- evaluate_normalization_methods(expr_for_norm, eval_methods, row_data)
    }

    # ---- 8. Missingness summary ----
    miss_per_sample <- colSums(is.na(expr_filt)) / nrow(expr_filt)
    miss_per_feature <- rowSums(is.na(expr_filt)) / ncol(expr_filt)
    miss_summary <- list(
        per_sample  = miss_per_sample,
        per_feature = miss_per_feature,
        total_na    = na_count,
        total_cells = prod(dim(expr_filt)),
        pct_missing = round(100 * na_count / prod(dim(expr_filt)), 2)
    )

    # ---- 9. Lipid class summary ----
    class_summary <- table(row_data$lipid_class)
    message(sprintf("lipidomics: %d features across %d lipid classes",
                    nrow(expr_filt), length(class_summary)))

    # ---- return pre-contract ----
    list(
        expr_raw   = expr_raw,
        expr_filt  = expr_filt,
        expr_log   = expr_log,
        expr_work  = expr_work,
        meta       = meta,
        row_data   = row_data,
        info       = list(
            mode     = "lipidomics",
            format   = format,
            input_already_normalized = input_already_normalized,
            na_policy = na_policy,
            normalization = norm_result$applied,
            n_features_raw  = nrow(expr_raw),
            n_features_filt = nrow(expr_filt),
            n_samples       = ncol(expr_work),
            n_lipid_classes = length(class_summary),
            lipid_class_counts = as.list(class_summary),
            missingness     = miss_summary
        ),
        normalization_eval = norm_eval
    )
}


# ---- helpers ----------------------------------------------------------------

#' Extract sample filter rules for lipidomics
get_sample_filter_rules_lipid <- function(cfg) {
    sf <- cfg$sample_filter
    if (is.null(sf) || !isTRUE(sf$enabled)) return(NULL)
    sf$rules %||% NULL
}


#' Apply sample filtering for lipidomics
#'
#' Supports: exclude_blanks, exclude_qc, exclude_pools, exclude_samples.
#' @return Character vector of sample IDs to keep.
apply_sample_filter_lipid <- function(sample_ids, meta, rules, sample_col) {
    keep <- rep(TRUE, length(sample_ids))

    if (isTRUE(rules$exclude_blanks) && "is_blank" %in% colnames(meta)) {
        keep <- keep & !meta$is_blank
    }

    if (isTRUE(rules$exclude_qc) && "is_QC" %in% colnames(meta)) {
        keep <- keep & !meta$is_QC
    }

    # Exclude pools by name pattern
    if (isTRUE(rules$exclude_pools)) {
        is_pool <- grepl("^Pool", sample_ids, ignore.case = TRUE)
        keep <- keep & !is_pool
    }

    # Explicit exclude list
    if (!is.null(rules$exclude_samples)) {
        keep <- keep & !(sample_ids %in% rules$exclude_samples)
    }

    sample_ids[keep]
}
