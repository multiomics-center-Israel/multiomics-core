#' Select features for multi-omics integration
#'
#' Implements multiple strategies to identify the most informative features
#' for integration analysis:
#' 1. Differential expression-based (top DE features per omics)
#' 2. Variance-based (high-variance features)
#' 3. Combined (intersection or union of methods)
#'
#' @param mae MultiAssayExperiment object
#' @param de_results Named list of DE results per omics (optional)
#' @param config Full config object
#' @return List with selected feature IDs per omics and selection metadata
select_features_for_integration <- function(mae, de_results = NULL, config) {

    cfg <- config$modes$multiomics$feature_selection %||% list()
    method <- cfg$method %||% "variance"  # "de", "variance", "combined"
    default_top_n <- cfg$top_n %||% 500
    per_omics_n <- cfg$per_omics %||% list()
    min_var_percentile <- cfg$min_var_percentile %||% 50

    omics <- names(MultiAssayExperiment::experiments(mae))
    selected_features <- list()
    selection_stats <- list()

    for (om in omics) {
        expr <- SummarizedExperiment::assay(mae[[om]], "expr")

        # Resolve per-omics top_n (fall back to global default)
        top_n <- per_omics_n[[om]] %||% default_top_n

        # Normalize DE result to a flat data frame (feature_id, logFC, padj)
        # before passing to select_features_by_de(). Upstream DE objects are
        # structured lists (RNA: list(dds, tables); prot/metab: list(summary_df,
        # ...)), not flat data frames, so normalization is always required.
        de_flat <- NULL
        if (method %in% c("de", "combined") && !is.null(de_results[[om]])) {
            de_flat <- tryCatch(
                normalize_de_for_concordance(de_results[[om]], omics_type = om),
                error = function(e) {
                    warning(
                        "Could not normalize DE results for ", om, ": ", e$message,
                        ". Falling back to variance-based selection."
                    )
                    NULL
                }
            )
        }

        if (method == "de" && !is.null(de_flat)) {
            # Select top DE features
            sel <- select_features_by_de(expr, de_flat, top_n)
            selected_features[[om]] <- sel$features
            selection_stats[[om]] <- sel$stats

        } else if (method == "de" && is.null(de_flat)) {
            # DE normalization failed or no DE results -- fall back to variance
            message(sprintf("%s: no usable DE results; falling back to variance selection", om))
            sel <- select_features_by_variance(expr, top_n, min_var_percentile)
            selected_features[[om]] <- sel$features
            selection_stats[[om]] <- sel$stats

        } else if (method == "variance") {
            # Select high-variance features
            sel <- select_features_by_variance(expr, top_n, min_var_percentile)
            selected_features[[om]] <- sel$features
            selection_stats[[om]] <- sel$stats

        } else if (method == "combined") {
            # Combine DE and variance
            sel_de <- if (!is.null(de_flat)) {
                select_features_by_de(expr, de_flat, top_n)$features
            } else character(0)

            sel_var <- select_features_by_variance(expr, top_n, min_var_percentile)$features

            # Union or intersection
            combine_mode <- cfg$combine_mode %||% "union"
            if (combine_mode == "union") {
                sel_combined <- union(sel_de, sel_var)
            } else {
                sel_combined <- intersect(sel_de, sel_var)
            }

            selected_features[[om]] <- sel_combined
            selection_stats[[om]] <- list(
                method = "combined",
                n_de = length(sel_de),
                n_var = length(sel_var),
                n_final = length(sel_combined)
            )

        } else {
            stop("Unknown feature selection method: ", method)
        }

        message(sprintf(
            "%s feature selection (top_n=%d): %d / %d features (%.1f%%)",
            om, top_n, length(selected_features[[om]]), nrow(expr),
            100 * length(selected_features[[om]]) / nrow(expr)
        ))
    }

    list(
        selected_features = selected_features,
        selection_stats = selection_stats,
        method = method,
        config = cfg
    )
}


#' Select features by differential expression
select_features_by_de <- function(expr, de_results, top_n) {

    # Assume de_results is a data frame with columns: feature_id, padj, logFC
    if (!"feature_id" %in% names(de_results)) {
        de_results$feature_id <- rownames(de_results)
    }

    # Filter to significant features
    sig <- de_results[!is.na(de_results$padj) & de_results$padj < 0.05, ]

    if (nrow(sig) == 0) {
        warning("No significant DE features found; using all features")
        sig <- de_results
    }

    # Rank by absolute logFC
    if ("logFC" %in% names(sig)) {
        sig <- sig[order(abs(sig$logFC), decreasing = TRUE), ]
    } else if ("log2FoldChange" %in% names(sig)) {
        sig <- sig[order(abs(sig$log2FoldChange), decreasing = TRUE), ]
    }

    # Take top N
    n_select <- min(top_n, nrow(sig))
    selected <- sig$feature_id[seq_len(n_select)]

    # Filter to features present in expression matrix
    selected <- intersect(selected, rownames(expr))

    list(
        features = selected,
        stats = list(
            method = "de",
            n_sig = nrow(sig),
            n_selected = length(selected),
            mean_padj = mean(sig$padj[seq_len(n_select)], na.rm = TRUE)
        )
    )
}


#' Select features by variance
select_features_by_variance <- function(expr, top_n, min_var_percentile = 50) {

    # Compute row variances
    row_vars <- apply(expr, 1, var, na.rm = TRUE)

    # Filter to top percentile
    var_threshold <- quantile(row_vars, probs = min_var_percentile / 100, na.rm = TRUE)
    high_var_idx <- which(row_vars >= var_threshold)

    # Rank by variance
    ranked_idx <- high_var_idx[order(row_vars[high_var_idx], decreasing = TRUE)]

    # Take top N
    n_select <- min(top_n, length(ranked_idx))
    selected <- rownames(expr)[ranked_idx[seq_len(n_select)]]

    list(
        features = selected,
        stats = list(
            method = "variance",
            n_high_var = length(high_var_idx),
            n_selected = n_select,
            var_threshold = var_threshold,
            median_var_selected = median(row_vars[ranked_idx[seq_len(n_select)]], na.rm = TRUE)
        )
    )
}


#' Subset MAE to selected features
#'
#' @param mae MultiAssayExperiment object
#' @param selected_features Named list of feature IDs per omics
#' @return Subsetted MAE
subset_mae_by_features <- function(mae, selected_features) {

    omics <- names(MultiAssayExperiment::experiments(mae))

    for (om in omics) {
        if (om %in% names(selected_features)) {
            feat <- selected_features[[om]]
            # Subset to selected features
            mae[[om]] <- mae[[om]][feat, ]
        }
    }

    message(sprintf(
        "MAE subsetted to selected features: %s",
        paste(
            sprintf("%s=%d", omics, sapply(mae@ExperimentList, nrow)),
            collapse = ", "
        )
    ))

    mae
}
