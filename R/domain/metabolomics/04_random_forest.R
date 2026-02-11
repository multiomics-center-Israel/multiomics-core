# R/domain/metabolomics/04_random_forest.R
#
# Feature importance via random forest classification.
# Uses ranger (preferred) or randomForest as fallback.
#
# Operates on the standard pre-processing contract (expr_work, meta).
# Reuses: assert_numeric_matrix, %||%


# ---- public entry point -----------------------------------------------------

#' Run random forest feature importance
#'
#' @param pre    List from preprocess_metabolomics() (pre contract).
#' @param config Full pipeline config.
#' @return list(importance_df, model, method) or NULL if disabled/unavailable.
run_metabolomics_rf <- function(pre, config) {
    cfg <- config$modes$metabolomics
    rf_cfg <- cfg$rf %||% list()

    if (!isTRUE(rf_cfg$run_rf)) {
        message("metabolomics RF: disabled in config — skipping")
        return(NULL)
    }

    condition_col <- rf_cfg$condition_column %||%
                     cfg$de$condition_column %||%
                     cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"

    mat  <- pre$expr_work
    meta <- pre$meta
    assert_numeric_matrix(mat, "metab_expr_work")

    # Align and build response
    meta <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- factor(meta[[condition_col]])

    # Prepare X: samples (rows) x features (cols)
    X <- t(mat)

    # Impute NAs with column medians
    for (j in seq_len(ncol(X))) {
        nas <- is.na(X[, j])
        if (any(nas)) X[nas, j] <- stats::median(X[!nas, j], na.rm = TRUE)
    }

    n_trees    <- rf_cfg$n_trees    %||% 500
    importance <- rf_cfg$importance  %||% "permutation"
    seed       <- rf_cfg$seed       %||% 1234

    # Try ranger first, then randomForest
    use_ranger <- requireNamespace("ranger", quietly = TRUE)

    if (use_ranger) {
        message("metabolomics RF: using ranger (", n_trees, " trees, ",
                importance, " importance)")
        set.seed(seed)
        rf_fit <- ranger::ranger(
            x = X, y = condition,
            num.trees = n_trees, importance = importance, seed = seed
        )
        imp <- ranger::importance(rf_fit)
        backend <- "ranger"
    } else if (requireNamespace("randomForest", quietly = TRUE)) {
        message("metabolomics RF: ranger unavailable, using randomForest")
        set.seed(seed)
        # randomForest needs safe column names
        orig_names <- colnames(X)
        safe_names <- paste0("F", seq_len(ncol(X)))
        colnames(X) <- safe_names
        name_map <- stats::setNames(orig_names, safe_names)

        rf_fit <- randomForest::randomForest(
            x = X, y = condition, ntree = n_trees, importance = TRUE
        )
        imp <- randomForest::importance(rf_fit)[, "MeanDecreaseAccuracy"]
        names(imp) <- name_map[names(imp)]
        backend <- "randomForest"
    } else {
        message("metabolomics RF: no RF package available — skipping")
        return(NULL)
    }

    # Build importance table
    imp_df <- data.frame(
        feature_id = names(imp),
        importance = unname(imp),
        stringsAsFactors = FALSE
    )
    imp_df <- imp_df[order(imp_df$importance, decreasing = TRUE), ]
    rownames(imp_df) <- NULL

    message("metabolomics RF complete: ", nrow(imp_df), " features scored")

    list(
        importance_df = imp_df,
        model         = rf_fit,
        method        = backend
    )
}


# ---- plot helper ------------------------------------------------------------

#' Random forest importance bar plot
#'
#' @param imp_df data.frame with feature_id and importance columns.
#' @param top_n  Number of top features to display.
#' @return ggplot object.
plot_rf_importance <- function(imp_df, top_n = 20) {
    top_df <- utils::head(imp_df, top_n)
    top_df$feature_id <- factor(top_df$feature_id,
                                 levels = rev(top_df$feature_id))

    ggplot2::ggplot(top_df, ggplot2::aes(x = feature_id, y = importance)) +
        ggplot2::geom_col(fill = "steelblue") +
        ggplot2::coord_flip() +
        ggplot2::labs(
            title = paste0("Random Forest — Top ", nrow(top_df), " Features"),
            x = NULL, y = "Importance"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13)
        )
}
