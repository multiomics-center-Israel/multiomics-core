# R/domain/metabolomics/04_feature_selection.R
#
# Feature selection for metabolomics:
#   1. Random forest importance (ranger / randomForest)
#   2. PLS-DA + VIP scores (mixOmics)
#
# Operates on the standard pre-processing contract (expr_work, meta).
# Reuses: assert_numeric_matrix, %||%


# ==== RANDOM FOREST ==========================================================

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


# ---- RF plot helper ---------------------------------------------------------

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
            title = paste0("Random Forest \u2014 Top ", nrow(top_df), " Features"),
            x = NULL, y = "Importance"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13)
        )
}


# ==== PLS-DA ==================================================================

# ---- public entry point -----------------------------------------------------

#' Run PLS-DA analysis with VIP scores
#'
#' @param pre    List from preprocess_metabolomics() (pre contract).
#' @param config Full pipeline config.
#' @return list(model, vip_scores, explained_variance, vip_df) or NULL.
run_metabolomics_plsda <- function(pre, config) {
    cfg <- config$modes$metabolomics
    plsda_cfg <- cfg$plsda %||% list()

    if (!isTRUE(plsda_cfg$run_plsda)) {
        message("metabolomics PLS-DA: disabled in config — skipping")
        return(NULL)
    }

    if (!requireNamespace("mixOmics", quietly = TRUE)) {
        message("metabolomics PLS-DA: mixOmics not available — skipping")
        return(NULL)
    }

    condition_col <- plsda_cfg$condition_column %||%
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

    # Remove zero-variance features
    col_vars <- apply(X, 2, stats::var, na.rm = TRUE)
    keep <- col_vars > 0
    if (sum(keep) < ncol(X)) {
        message("metabolomics PLS-DA: removed ", sum(!keep),
                " zero-variance features")
        X <- X[, keep, drop = FALSE]
    }

    # Fit PLS-DA
    ncomp <- min(2, ncol(X), nrow(X) - 1)
    message("metabolomics PLS-DA: fitting ", nrow(X), " samples x ",
            ncol(X), " features, ncomp=", ncomp)

    plsda_model <- mixOmics::plsda(X, condition, ncomp = ncomp)
    expl_var    <- plsda_model$prop_expl_var$X
    vip_scores  <- mixOmics::vip(plsda_model)

    # Build VIP table
    vip_df <- data.frame(
        feature_id = rownames(vip_scores),
        VIP_comp1  = vip_scores[, 1],
        stringsAsFactors = FALSE
    )
    if (ncol(vip_scores) >= 2) vip_df$VIP_comp2 <- vip_scores[, 2]

    # Annotate with molecule names from row_data
    if (!is.null(pre$row_data) && "Name" %in% colnames(pre$row_data)) {
        name_map <- stats::setNames(
            as.character(pre$row_data$Name),
            as.character(pre$row_data$feature_id)
        )
        vip_df$Name <- name_map[vip_df$feature_id]
    }

    vip_df <- vip_df[order(vip_df$VIP_comp1, decreasing = TRUE), ]
    rownames(vip_df) <- NULL

    message("metabolomics PLS-DA complete: variance explained = ",
            paste(round(expl_var * 100, 1), "%", collapse = ", "))

    list(
        model              = plsda_model,
        vip_scores         = vip_scores,
        explained_variance = expl_var,
        vip_df             = vip_df,
        X                  = X,
        condition          = condition
    )
}


# ---- PLS-DA plot helpers ----------------------------------------------------

#' PLS-DA scores plot
#'
#' @param plsda_res Result from run_metabolomics_plsda().
#' @param colors    Optional character vector of colors.
#' @return ggplot object.
plot_plsda_scores <- function(plsda_res, colors = NULL) {
    model     <- plsda_res$model
    condition <- plsda_res$condition
    expl_var  <- plsda_res$explained_variance

    scores  <- model$variates$X
    pct_var <- round(expl_var * 100, 1)

    df <- data.frame(
        Comp1     = scores[, 1],
        Comp2     = scores[, 2],
        condition = condition,
        stringsAsFactors = FALSE
    )

    default_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                        "#FF7F00", "#A65628")
    colors <- colors %||% default_colors
    n_groups <- nlevels(condition)
    colors <- colors[seq_len(n_groups)]

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Comp1, y = Comp2,
                                           color = condition)) +
        ggplot2::geom_point(size = 3.5, alpha = 0.85) +
        ggplot2::scale_color_manual(values = colors)

    # 95% confidence ellipses (need >= 3 per group)
    grp_counts <- table(condition)
    if (all(grp_counts >= 3)) {
        p <- p + ggplot2::stat_ellipse(type = "norm", level = 0.95,
                                        linewidth = 0.8, show.legend = FALSE)
    }

    p +
        ggplot2::labs(
            title = "PLS-DA Scores Plot",
            x = paste0("Component 1 (", pct_var[1], "%)"),
            y = paste0("Component 2 (", pct_var[2], "%)"),
            color = ""
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title     = ggplot2::element_text(hjust = 0.5, face = "bold",
                                                    size = 14),
            axis.title     = ggplot2::element_text(size = 12),
            legend.text    = ggplot2::element_text(size = 11),
            legend.position = "bottom"
        )
}


#' PLS-DA VIP scores lollipop plot
#'
#' @param plsda_res Result from run_metabolomics_plsda().
#' @param top_n     Number of top features to display.
#' @param colors    Optional character vector of colors.
#' @return ggplot object.
plot_plsda_vip <- function(plsda_res, top_n = 15, colors = NULL) {
    vip_scores <- plsda_res$vip_scores
    X          <- plsda_res$X
    condition  <- plsda_res$condition
    row_data   <- plsda_res$vip_df

    # Top features by comp1 VIP
    vip_vec    <- vip_scores[, 1]
    vip_sorted <- sort(vip_vec, decreasing = TRUE)
    top_feats  <- names(vip_sorted)[seq_len(min(top_n, length(vip_sorted)))]

    # Determine which group has highest mean
    groups <- levels(condition)
    high_group <- vapply(top_feats, function(feat) {
        means <- tapply(X[, feat], condition, mean, na.rm = TRUE)
        groups[which.max(means)]
    }, character(1))

    # Display names from vip_df if Name column exists
    display_names <- top_feats
    if ("Name" %in% colnames(row_data)) {
        name_map <- stats::setNames(row_data$Name, row_data$feature_id)
        display_names <- ifelse(
            is.na(name_map[top_feats]) | name_map[top_feats] == "",
            top_feats, name_map[top_feats]
        )
    }

    df <- data.frame(
        feature = factor(display_names, levels = rev(display_names)),
        VIP     = vip_sorted[top_feats],
        high_in = high_group,
        stringsAsFactors = FALSE
    )

    default_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                        "#FF7F00", "#A65628")
    colors <- colors %||% default_colors
    n_groups <- length(groups)
    group_colors <- stats::setNames(colors[seq_len(n_groups)], groups)

    ggplot2::ggplot(df, ggplot2::aes(x = VIP, y = feature, color = high_in)) +
        ggplot2::geom_segment(
            ggplot2::aes(x = 0, xend = VIP, y = feature, yend = feature),
            color = "grey70", linewidth = 0.5
        ) +
        ggplot2::geom_point(size = 4) +
        ggplot2::scale_color_manual(
            values = group_colors,
            labels = paste("High in", names(group_colors))
        ) +
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                             color = "grey40", linewidth = 0.4) +
        ggplot2::labs(
            title = paste0("PLS-DA VIP Scores \u2014 Top ", nrow(df), " Features"),
            x = "VIP Score", y = NULL, color = ""
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title     = ggplot2::element_text(hjust = 0.5, face = "bold",
                                                    size = 14),
            axis.title     = ggplot2::element_text(size = 12),
            axis.text.y    = ggplot2::element_text(size = 9),
            legend.text    = ggplot2::element_text(size = 11),
            legend.position = "bottom"
        )
}
