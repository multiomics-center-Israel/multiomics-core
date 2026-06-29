# R/domain/lipidomics/02_preprocess.R
#
# Preprocessing assembly for lipidomics.
# Dispatches to format-specific parsers, builds the standard internal
# representation, and applies normalization.
#
# Normalization reuses the metabolomics primitives:
#   normalize_samples, transform_metab, scale_metab
# composed here by apply_lipid_normalization() (this layer's only consumer of
# the sample_norm -> transform -> scaling chain; metabolomics drives the same
# primitives through its own module path).
#
# Reuses: assert_numeric_matrix, assert_meta_contract, align_meta_to_matrix,
#         build_minimal_meta


#' Apply the full normalization pipeline: sample_norm -> transform -> scaling
#'
#' Thin lipidomics-side wrapper that orchestrates the metabolomics
#' normalization primitives (normalize_samples / transform_metab /
#' scale_metab). Lives here because lipidomics is its only caller.
#'
#' @param mat       Numeric matrix (features x samples).
#' @param norm_cfg  normalization config section.
#' @param row_data  data.frame (for IS normalization only).
#' @param groups    Character/factor vector of group labels per sample
#'                  (required for eigenms; ignored by other methods).
#' @param meta      Sample metadata data.frame (for bio_factor normalization).
#' @return list(expr_norm, applied) where applied records what was done.
apply_lipid_normalization <- function(mat, norm_cfg, row_data = NULL,
                                      groups = NULL, meta = NULL) {
    sample_norm <- norm_cfg$sample_norm %||% "none"
    transform   <- norm_cfg$transform   %||% "none"
    scaling     <- norm_cfg$scaling     %||% "none"
    pseudocount <- norm_cfg$pseudocount %||% 1
    is_ref_col     <- norm_cfg$is_ref_col            %||% NULL
    bio_factor_col <- norm_cfg$biological_factor_col %||% NULL

    mat <- normalize_samples(mat, method = sample_norm,
                             ref_col = is_ref_col, row_data = row_data,
                             groups = groups, meta = meta,
                             bio_factor_col = bio_factor_col)
    mat <- transform_metab(mat, method = transform, pseudocount = pseudocount)
    mat <- scale_metab(mat, method = scaling)

    list(
        expr_norm = mat,
        applied = list(
            sample_norm = sample_norm,
            transform   = transform,
            scaling     = scaling,
            pseudocount = if (transform != "none") pseudocount else NA
        )
    )
}


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

    # ---- 2.5 Pool sub-matrix for QC ----
    # Captured upstream in parse_lipidsearch() before exclude_patterns strips
    # them. Pull through unchanged here; downstream Pool CV QC uses it even
    # when pools are also excluded from the analysis matrix.
    pool_matrix <- parsed$pool_matrix

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

    # ---- 4a. Group-aware missingness filter ----
    # Drop features whose missing fraction (NA or zero) exceeds
    # max_group_missing_frac in EVERY biological group. Retains features
    # with robust signal in at least one group.
    feat_cfg               <- cfg$feature_filter %||% list()
    feat_filter_enabled    <- isTRUE(feat_cfg$enabled %||% TRUE)
    max_group_missing_frac <- feat_cfg$max_group_missing_frac %||% 0.30

    n_features_pre_group_filter <- nrow(expr_filt)
    n_dropped_group_filter      <- 0L
    group_filter_applied        <- FALSE
    group_filter_skip_reason    <- NA_character_

    if (feat_filter_enabled) {
        group_col <- cfg$de$condition_column %||% cfg$effects$color
        if (is.null(group_col) || !group_col %in% colnames(meta)) {
            group_filter_skip_reason <- sprintf(
                "group column '%s' not found on metadata", group_col %||% "(unset)"
            )
            message("lipidomics: group-missingness filter skipped (",
                    group_filter_skip_reason, ").")
        } else {
            keep_feat <- filter_features_by_group_missingness(
                expr_filt, meta, sample_col, group_col, max_group_missing_frac
            )
            n_dropped_group_filter <- sum(!keep_feat)
            expr_filt <- expr_filt[keep_feat, , drop = FALSE]
            row_data  <- row_data[keep_feat, , drop = FALSE]
            group_filter_applied <- TRUE
            message(sprintf(
                "lipidomics: group-missingness filter (%.0f%% in every group of '%s') dropped %d/%d features.",
                max_group_missing_frac * 100, group_col,
                n_dropped_group_filter, n_features_pre_group_filter
            ))
        }
    }

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

    # Extract group labels for methods that need them (e.g. eigenms)
    group_col <- cfg$de$condition_column %||% cfg$effects$color %||% "sample_type"
    norm_groups <- if (group_col %in% colnames(meta)) as.character(meta[[group_col]]) else NULL

    # Pre-scaling matrix for logFC computation
    norm_cfg_no_scale <- norm_cfg
    norm_cfg_no_scale$scaling <- "none"
    pre_scale_result <- apply_lipid_normalization(expr_for_norm, norm_cfg_no_scale, row_data,
                                                  groups = norm_groups, meta = meta)
    expr_log <- pre_scale_result$expr_norm

    # Full pipeline (with scaling)
    norm_result <- apply_lipid_normalization(expr_for_norm, norm_cfg, row_data,
                                             groups = norm_groups, meta = meta)
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
        pool_matrix = pool_matrix,
        info       = list(
            mode     = "lipidomics",
            format   = format,
            input_already_normalized = input_already_normalized,
            na_policy = na_policy,
            normalization = norm_result$applied,
            n_features_raw  = nrow(expr_raw),
            n_features_filt = nrow(expr_filt),
            n_samples       = ncol(expr_work),
            n_pools         = if (is.null(pool_matrix)) 0L else ncol(pool_matrix),
            n_lipid_classes = length(class_summary),
            lipid_class_counts = as.list(class_summary),
            missingness     = miss_summary,
            feature_filter  = list(
                enabled                = feat_filter_enabled,
                max_group_missing_frac = max_group_missing_frac,
                applied                = group_filter_applied,
                skip_reason            = group_filter_skip_reason,
                n_before               = n_features_pre_group_filter,
                n_dropped              = n_dropped_group_filter
            )
        ),
        normalization_eval = norm_eval
    )
}


# ---- helpers ----------------------------------------------------------------

#' Filter features by group-aware missingness
#'
#' Returns a logical keep-vector. A feature is kept if at least one biological
#' group has missing fraction (NA or zero) <= `max_frac`. A feature is dropped
#' only when missingness exceeds `max_frac` in EVERY group. Samples whose
#' group label is NA are ignored. Groups with zero non-NA samples after that
#' contribute no information and are skipped.
#'
#' @param expr        Numeric matrix, features x samples.
#' @param meta        Sample metadata data.frame (must contain sample_col, group_col).
#' @param sample_col  Column in `meta` whose values match colnames(expr).
#' @param group_col   Column in `meta` giving the biological group.
#' @param max_frac    Maximum tolerated missingness in any group (default 0.30).
#' @return Logical vector of length nrow(expr): TRUE = keep.
filter_features_by_group_missingness <- function(expr, meta, sample_col,
                                                  group_col, max_frac = 0.30) {
    if (nrow(expr) == 0L) return(logical(0))

    samp_to_grp <- setNames(as.character(meta[[group_col]]),
                            as.character(meta[[sample_col]]))
    grp_per_col <- samp_to_grp[colnames(expr)]
    valid_cols  <- !is.na(grp_per_col)
    if (!any(valid_cols)) {
        warning("filter_features_by_group_missingness: no samples have a valid group label; keeping all features.")
        return(rep(TRUE, nrow(expr)))
    }

    expr_v <- expr[, valid_cols, drop = FALSE]
    grp    <- grp_per_col[valid_cols]
    groups <- unique(grp)

    is_miss <- is.na(expr_v) | expr_v == 0

    # Per-group missingness fraction, features x groups
    miss_frac <- vapply(
        groups,
        function(g) {
            cols_g <- which(grp == g)
            rowMeans(is_miss[, cols_g, drop = FALSE])
        },
        numeric(nrow(expr_v))
    )
    if (!is.matrix(miss_frac)) {
        miss_frac <- matrix(miss_frac, nrow = nrow(expr_v))
    }

    # Keep if at least one group is below the threshold
    apply(miss_frac, 1, function(x) any(x <= max_frac))
}


#' Extract sample filter rules for lipidomics
get_sample_filter_rules_lipid <- function(cfg) {
    sf <- cfg$sample_filter
    if (is.null(sf) || !isTRUE(sf$enabled)) return(NULL)
    # The web wizard writes exclusions (e.g. exclude_samples) directly under
    # sample_filter rather than nesting them under sample_filter$rules; fall
    # back to sf itself so those direct rules are honoured.
    sf$rules %||% sf
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


#' Compare alternative normalization pipelines on a single matrix
#'
#' For each method spec (sample_norm/transform/scaling), applies
#' apply_lipid_normalization() and reports median RSD and PC1 variance —
#' the same diagnostics used by mod_met_norm_comparison.
#'
#' @param expr_for_norm Numeric matrix (features x samples), pre-normalization.
#' @param methods       List of method specs, each a list with optional
#'                      sample_norm/transform/scaling/pseudocount/na_policy.
#' @param row_data      data.frame of feature metadata (for is/bio_factor norms).
#' @return data.frame: label, sample_norm, transform, scaling, median_rsd, pc1_var.
#'   The per-method full-pipeline matrices are attached as `attr(df, "matrices")`
#'   (named list keyed by label) for downstream per-method visualisation. The
#'   attribute survives RDS serialisation but is dropped by TSV writers.
evaluate_normalization_methods <- function(expr_for_norm, methods, row_data = NULL) {
    if (!length(methods)) return(NULL)

    run_norm <- function(spec, label) {
        tryCatch(
            apply_lipid_normalization(expr_for_norm, spec, row_data)$expr_norm,
            error = function(e) {
                message("evaluate_normalization_methods: ", label, " failed: ",
                        conditionMessage(e))
                NULL
            }
        )
    }

    method_mats <- list()

    rows <- lapply(methods, function(spec) {
        sn <- spec$sample_norm %||% "none"
        tr <- spec$transform   %||% "none"
        sc <- spec$scaling     %||% "none"
        label <- paste(sn, tr, sc, sep = "+")

        # RSD on pre-scaling matrix: after auto/pareto scaling, feature means
        # collapse to ~0 and SD/|mean| explodes. RSD is meaningful only when
        # features retain their intensity scale.
        spec_no_scale <- spec
        spec_no_scale$scaling <- "none"
        mat_for_rsd <- run_norm(spec_no_scale, paste0(label, " (RSD)"))

        median_rsd <- if (is.null(mat_for_rsd) || ncol(mat_for_rsd) < 2) {
            NA_real_
        } else {
            feat_means <- rowMeans(mat_for_rsd, na.rm = TRUE)
            feat_sds   <- apply(mat_for_rsd, 1, stats::sd, na.rm = TRUE)
            rsd        <- feat_sds / abs(feat_means)
            stats::median(rsd[is.finite(rsd)], na.rm = TRUE)
        }

        # PC1 variance on the fully-pipelined matrix (the one actually used
        # downstream) — captures whether scaling spreads variance off PC1.
        mat_full <- if (identical(sc, "none")) mat_for_rsd else run_norm(spec, label)
        pc1_var <- if (is.null(mat_full) || ncol(mat_full) < 2) {
            NA_real_
        } else {
            tryCatch(
                compute_pca_scores(mat_full, pcs = 1L,
                                   center = TRUE, scale = FALSE)$var_expl[1],
                error = function(e) NA_real_
            )
        }

        if (!is.null(mat_full)) method_mats[[label]] <<- mat_full

        data.frame(label = label, sample_norm = sn, transform = tr, scaling = sc,
                   median_rsd = round(median_rsd, 4),
                   pc1_var = round(pc1_var, 4),
                   stringsAsFactors = FALSE)
    })

    out <- do.call(rbind, rows)
    attr(out, "matrices") <- method_mats
    out
}
