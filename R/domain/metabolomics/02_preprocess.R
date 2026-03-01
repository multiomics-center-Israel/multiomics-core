# R/domain/metabolomics/02_preprocess.R
#
# Preprocessing assembly for metabolomics.
# Dispatches to format-specific parsers, builds the standard internal
# representation, and applies normalization.
#
# Reuses: assert_numeric_matrix, assert_meta_contract, align_meta_to_expr,
#         apply_normalization_pipeline, evaluate_normalization_methods,
#         build_minimal_meta, parse_cd_raw, parse_processed_wide


#' Preprocess metabolomics data
#'
#' @param inputs  List from load_metabolomics_inputs().
#' @param config  Full pipeline config.
#' @return list matching the pre-processing contract:
#'   expr_raw, expr_filt, expr_work, meta, row_data, info,
#'   normalization_eval (NULL if evaluate_methods not set)
preprocess_metabolomics <- function(inputs, config) {
    cfg <- config$modes$metabolomics
    format <- inputs$format %||% cfg$input$format %||% "cd_raw"

    # ---- 1. Parse input into expr_raw / row_data / sample_ids ----
    parsed <- switch(format,
        cd_raw = parse_cd_raw(inputs$data, cfg),
        processed_wide = parse_processed_wide(inputs$data, cfg, inputs$metadata),
        stop("Unsupported metabolomics input format: ", format)
    )

    expr_raw <- parsed$expr_raw
    row_data <- parsed$row_data
    sample_ids <- colnames(expr_raw)

    assert_numeric_matrix(expr_raw, "metab_expr_raw")

    # ---- 2. Build / align metadata ----
    meta <- inputs$metadata
    sample_col <- cfg$effects$samples %||% "sample_id"

    if (is.null(meta)) {
        # No external metadata provided — build minimal
        meta <- build_minimal_meta(sample_ids)
        sample_col <- "sample_id"

        # Update effects for downstream QC
        if (is.null(cfg$effects$color)) cfg$effects$color <- "sample_type"
        if (is.null(cfg$effects$samples)) cfg$effects$samples <- "sample_id"
    } else {
        check_has_cols(meta, sample_col, df_name = "metadata")
        meta <- align_meta_to_matrix(sample_ids, meta, sample_col)
    }

    assert_meta_contract(meta, sample_col)

    # ---- 3. Optional sample filtering (QC/blank removal) ----
    rules <- get_sample_filter_rules_metab(cfg)
    if (!is.null(rules)) {
        keep <- apply_sample_filter_metab(sample_ids, meta, rules, sample_col)
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

    # ---- 4b. Annotate feature names from HMDB lookup (if missing) ----
    row_data <- annotate_hmdb_names(row_data, config)

    norm_cfg_raw <- cfg$normalization %||% list()

    # ---- 5. Handle missing values ----
    na_count <- sum(is.na(expr_filt))
    na_policy <- norm_cfg_raw$na_policy %||% "keep"
    na_policy <- tolower(na_policy)

    expr_for_norm <- expr_filt
    if (na_policy == "zero") {
        if (na_count > 0) {
            message(sprintf("metabolomics: replacing %d NA values with 0 (na_policy='zero').", na_count))
        }
        expr_for_norm[is.na(expr_for_norm)] <- 0
    } else if (na_policy == "lod") {
        # Replace zeros and NAs with a fraction of the minimum positive value per feature
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
            message(sprintf("metabolomics: imputed %d zero/NA values with %.1f%% of min per feature (na_policy='lod').",
                            n_replaced, lod_frac * 100))
        }
    } else if (na_policy == "min_half") {
        # Replace zeros and NAs with half the minimum positive value per feature
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
            message(sprintf("metabolomics: imputed %d zero/NA values with min/2 per feature (na_policy='min_half').", n_replaced))
        }
    } else if (na_count > 0) {
        message(sprintf("metabolomics: keeping %d NA values as-is (na_policy='keep').", na_count))
    }

    # ---- 6. Normalization ----
    norm_cfg <- norm_cfg_raw
    input_already_normalized <- isTRUE(norm_cfg$input_already_normalized)

    norm_cfg$sample_norm <- norm_cfg$sample_norm %||% "none"
    norm_cfg$transform   <- norm_cfg$transform   %||% "none"
    norm_cfg$scaling     <- norm_cfg$scaling     %||% "none"

    # Run sample_norm + transform first (pre-scaling matrix for logFC)
    norm_cfg_no_scale <- norm_cfg
    norm_cfg_no_scale$scaling <- "none"
    pre_scale_result <- apply_normalization_pipeline(expr_for_norm, norm_cfg_no_scale, row_data)
    expr_log <- pre_scale_result$expr_norm

    # Full pipeline (with scaling) for statistical tests
    norm_result <- apply_normalization_pipeline(expr_for_norm, norm_cfg, row_data)
    expr_work <- norm_result$expr_norm

    assert_numeric_matrix(expr_work, "metab_expr_work")

    # ---- 7. Optional: evaluate alternative methods ----
    norm_eval <- NULL
    eval_methods <- norm_cfg$evaluate_methods
    if (!is.null(eval_methods) && length(eval_methods) > 0) {
        # If input_already_normalized, filter out methods that apply sample_norm
        # unless explicitly requested
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

    # ---- return pre-contract ----
    list(
        expr_raw   = expr_raw,
        expr_filt  = expr_filt,
        expr_log   = expr_log,
        expr_work  = expr_work,
        meta       = meta,
        row_data   = row_data,
        info       = list(
            mode     = "metabolomics",
            format   = format,
            input_already_normalized = input_already_normalized,
            na_policy = na_policy,
            normalization = norm_result$applied,
            n_features_raw  = nrow(expr_raw),
            n_features_filt = nrow(expr_filt),
            n_samples       = ncol(expr_work),
            missingness     = miss_summary
        ),
        normalization_eval = norm_eval,
        sample_map = parsed$sample_map %||% NULL
    )
}


# ---- helpers ----------------------------------------------------------------

#' Extract sample filter rules for metabolomics
get_sample_filter_rules_metab <- function(cfg) {
    sf <- cfg$sample_filter
    if (is.null(sf) || !isTRUE(sf$enabled)) return(NULL)
    sf$rules %||% NULL
}


#' Apply simple sample filtering
#'
#' Supports rules like exclude_blanks=TRUE and exclude_qc=TRUE.
#' @return Character vector of sample IDs to keep.
apply_sample_filter_metab <- function(sample_ids, meta, rules, sample_col) {
    keep <- rep(TRUE, length(sample_ids))

    if (isTRUE(rules$exclude_blanks) && "is_blank" %in% colnames(meta)) {
        keep <- keep & !meta$is_blank
    }

    if (isTRUE(rules$exclude_qc) && "is_QC" %in% colnames(meta)) {
        keep <- keep & !meta$is_QC
    }

    # Explicit exclude list
    if (!is.null(rules$exclude_samples)) {
        keep <- keep & !(sample_ids %in% rules$exclude_samples)
    }

    sample_ids[keep]
}


#' Annotate row_data with metabolite names from HMDB lookup table
#'
#' If row_data already has a populated Name column, this is a no-op.
#' Otherwise, looks up feature IDs in the bundled HMDB compound names table
#' at {project$dir}/data/hmdb_compound_names.tsv.
#' For non-HMDB feature IDs (already human-readable names), uses the ID itself.
#'
#' @param row_data  data.frame with at least a feature_id column.
#' @param config    Full pipeline config (for project$dir path).
#' @return row_data with a Name column populated.
annotate_hmdb_names <- function(row_data, config) {
    # Skip if Name column already exists and is mostly populated
    if ("Name" %in% colnames(row_data)) {
        n_populated <- sum(!is.na(row_data$Name) & nzchar(trimws(row_data$Name)))
        if (n_populated > nrow(row_data) * 0.5) return(row_data)
    }

    feat_ids <- as.character(row_data$feature_id)
    is_hmdb <- grepl("^HMDB[0-9]+$", feat_ids)

    # For non-HMDB IDs, the feature ID is itself a name
    names_out <- ifelse(is_hmdb, NA_character_, feat_ids)

    # Look for bundled lookup table
    proj_dir <- config$project$dir %||% "."
    lookup_path <- file.path(proj_dir, "data", "hmdb_compound_names.tsv")

    if (file.exists(lookup_path)) {
        lookup <- utils::read.delim(lookup_path, stringsAsFactors = FALSE)
        if (all(c("HMDB", "Name") %in% colnames(lookup))) {
            lut <- stats::setNames(lookup$Name, lookup$HMDB)
            hmdb_ids <- feat_ids[is_hmdb]
            matched <- lut[hmdb_ids]
            names_out[is_hmdb] <- ifelse(is.na(matched), hmdb_ids, matched)
            n_annotated <- sum(!is.na(matched))
            n_hmdb <- sum(is_hmdb)
            message(sprintf("metabolomics: annotated %d/%d HMDB features with compound names.",
                            n_annotated, n_hmdb))
        }
    } else {
        # No lookup available — use HMDB IDs as-is
        names_out[is_hmdb] <- feat_ids[is_hmdb]
    }

    row_data$Name <- names_out
    row_data
}
