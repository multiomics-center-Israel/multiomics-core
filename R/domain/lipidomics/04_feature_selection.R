# R/domain/lipidomics/04_feature_selection.R
#
# Feature selection for lipidomics: RF importance + PLS-DA VIP.
# Mirrors metabolomics feature selection but reads from config$modes$lipidomics.


# ==== RANDOM FOREST ==========================================================

#' Run random forest feature importance for lipidomics
#'
#' @param pre    List from preprocess_lipidomics() (pre contract).
#' @param config Full pipeline config.
#' @return list(importance_df, model, method) or NULL if disabled/unavailable.
run_lipidomics_rf <- function(pre, config) {
    cfg <- config$modes$lipidomics
    rf_cfg <- cfg$rf %||% list()

    if (!isTRUE(rf_cfg$run_rf)) {
        message("lipidomics RF: disabled in config — skipping")
        return(NULL)
    }

    condition_col <- rf_cfg$condition_column %||%
                     cfg$de$condition_column %||%
                     cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"

    mat  <- pre$expr_work
    meta <- pre$meta
    assert_numeric_matrix(mat, "lipid_expr_work")

    meta <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- factor(meta[[condition_col]])

    X <- t(mat)
    for (j in seq_len(ncol(X))) {
        nas <- is.na(X[, j])
        if (any(nas)) X[nas, j] <- stats::median(X[!nas, j], na.rm = TRUE)
    }

    n_trees    <- rf_cfg$n_trees    %||% 500
    importance <- rf_cfg$importance  %||% "permutation"
    seed       <- rf_cfg$seed       %||% 1234

    use_ranger <- requireNamespace("ranger", quietly = TRUE)

    if (use_ranger) {
        message("lipidomics RF: using ranger (", n_trees, " trees)")
        set.seed(seed)
        rf_fit <- ranger::ranger(
            x = X, y = condition,
            num.trees = n_trees, importance = importance, seed = seed
        )
        imp <- ranger::importance(rf_fit)
        backend <- "ranger"
    } else if (requireNamespace("randomForest", quietly = TRUE)) {
        message("lipidomics RF: using randomForest")
        set.seed(seed)
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
        message("lipidomics RF: no RF package available — skipping")
        return(NULL)
    }

    imp_df <- data.frame(
        feature_id = names(imp),
        importance = unname(imp),
        stringsAsFactors = FALSE
    )
    imp_df <- imp_df[order(imp_df$importance, decreasing = TRUE), ]
    rownames(imp_df) <- NULL

    list(importance_df = imp_df, model = rf_fit, method = backend)
}


# ==== PLS-DA ==================================================================

#' Run PLS-DA analysis with VIP scores for lipidomics
#'
#' @param pre    List from preprocess_lipidomics() (pre contract).
#' @param config Full pipeline config.
#' @return list(model, vip_scores, explained_variance, vip_df) or NULL.
run_lipidomics_plsda <- function(pre, config) {
    cfg <- config$modes$lipidomics
    plsda_cfg <- cfg$plsda %||% list()

    if (!isTRUE(plsda_cfg$run_plsda)) {
        message("lipidomics PLS-DA: disabled in config — skipping")
        return(NULL)
    }

    if (!requireNamespace("mixOmics", quietly = TRUE)) {
        message("lipidomics PLS-DA: mixOmics not available — skipping")
        return(NULL)
    }

    condition_col <- plsda_cfg$condition_column %||%
                     cfg$de$condition_column %||%
                     cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"

    mat  <- pre$expr_work
    meta <- pre$meta
    assert_numeric_matrix(mat, "lipid_expr_work")

    meta <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- factor(meta[[condition_col]])

    X <- t(mat)
    for (j in seq_len(ncol(X))) {
        nas <- is.na(X[, j])
        if (any(nas)) X[nas, j] <- stats::median(X[!nas, j], na.rm = TRUE)
    }

    col_vars <- apply(X, 2, stats::var, na.rm = TRUE)
    keep <- col_vars > 0
    if (sum(keep) < ncol(X)) {
        message("lipidomics PLS-DA: removed ", sum(!keep), " zero-variance features")
        X <- X[, keep, drop = FALSE]
    }

    ncomp <- min(2, ncol(X), nrow(X) - 1)
    message("lipidomics PLS-DA: fitting ", nrow(X), " samples x ",
            ncol(X), " features, ncomp=", ncomp)

    plsda_model <- mixOmics::plsda(X, condition, ncomp = ncomp)
    expl_var    <- plsda_model$prop_expl_var$X
    vip_scores  <- mixOmics::vip(plsda_model)

    vip_df <- data.frame(
        feature_id = rownames(vip_scores),
        VIP_comp1  = vip_scores[, 1],
        stringsAsFactors = FALSE
    )
    if (ncol(vip_scores) >= 2) vip_df$VIP_comp2 <- vip_scores[, 2]

    if (!is.null(pre$row_data) && "Name" %in% colnames(pre$row_data)) {
        name_map <- stats::setNames(
            as.character(pre$row_data$Name),
            as.character(pre$row_data$feature_id)
        )
        vip_df$Name <- name_map[vip_df$feature_id]
    }

    vip_df <- vip_df[order(vip_df$VIP_comp1, decreasing = TRUE), ]
    rownames(vip_df) <- NULL

    list(
        model              = plsda_model,
        vip_scores         = vip_scores,
        explained_variance = expl_var,
        vip_df             = vip_df,
        X                  = X,
        condition          = condition
    )
}
