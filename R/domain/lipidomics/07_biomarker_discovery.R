# R/domain/lipidomics/07_biomarker_discovery.R
#
# Biomarker discovery analysis for lipidomics (LipidOne approach).
#   1. Per-feature biomarker metrics (t-test, AUC, Cohen's d, power, fold change)
#   2. Biomarker table visualisation
#   3. Per-feature ROC curves
#   4. Multivariate ROC via Random Forest


# ==== BIOMARKER METRICS =====================================================

#' Compute biomarker metrics for each lipid feature
#'
#' For every feature in the expression matrix, computes a t-test p-value,
#' AUC (via pROC or manual Wilcoxon-based), Cohen's d, approximate
#' statistical power, and log2 fold change between conditions.
#'
#' @param pre    List from preprocess_lipidomics() (pre contract).
#' @param config Full pipeline config.
#' @return data.frame with columns: feature_id, p_value, AUC, cohens_d,
#'   power_pct, log2FC, is_candidate.  Sorted by p_value ascending;
#'   top 20 rows are flagged as candidate biomarkers.
compute_biomarker_metrics <- function(pre, config) {
    cfg <- config$modes$lipidomics

    condition_col <- cfg$de$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col    <- cfg$effects$samples %||% "sample_id"

    mat  <- pre$expr_work
    meta <- pre$meta

    meta <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- factor(meta[[condition_col]])

    lvls <- levels(condition)
    if (length(lvls) != 2) {
        message("biomarker metrics: condition has ", length(lvls),
                " levels — need exactly 2, skipping")
        return(NULL)
    }

    idx1 <- which(condition == lvls[1])
    idx2 <- which(condition == lvls[2])
    n1   <- length(idx1)
    n2   <- length(idx2)

    has_pROC <- requireNamespace("pROC", quietly = TRUE)
    has_pwr  <- requireNamespace("pwr",  quietly = TRUE)

    # Fold change must come from the pre-scaling log matrix; expr_work may be
    # centered/unit-variance scaled (e.g. scaling: auto), which would turn the
    # group-mean difference into a standardized effect size, not a log2FC.
    mat_log <- pre$expr_log %||% mat

    results <- lapply(seq_len(nrow(mat)), function(i) {
        x1 <- mat[i, idx1]
        x2 <- mat[i, idx2]

        x1 <- x1[!is.na(x1)]
        x2 <- x2[!is.na(x2)]

        if (length(x1) < 2 || length(x2) < 2) {
            return(data.frame(
                feature_id = rownames(mat)[i],
                p_value    = NA_real_,
                AUC        = NA_real_,
                cohens_d   = NA_real_,
                power_pct  = NA_real_,
                log2FC     = NA_real_,
                stringsAsFactors = FALSE
            ))
        }

        # --- p-value (Welch t-test) ---
        tt <- tryCatch(
            stats::t.test(x1, x2),
            error = function(e) NULL
        )
        pval <- if (!is.null(tt)) tt$p.value else NA_real_

        # --- AUC ---
        auc_val <- tryCatch({
            if (has_pROC) {
                resp   <- c(rep(0L, length(x1)), rep(1L, length(x2)))
                pred   <- c(x1, x2)
                roc_obj <- pROC::roc(resp, pred, quiet = TRUE)
                as.numeric(pROC::auc(roc_obj))
            } else {
                # Manual via Wilcoxon U-statistic
                wt <- stats::wilcox.test(x2, x1)
                U  <- wt$statistic
                as.numeric(U / (length(x1) * length(x2)))
            }
        }, error = function(e) NA_real_)

        # --- Cohen's d ---
        m1 <- mean(x1)
        m2 <- mean(x2)
        sd1 <- stats::sd(x1)
        sd2 <- stats::sd(x2)
        pooled_sd <- sqrt(((length(x1) - 1) * sd1^2 + (length(x2) - 1) * sd2^2) /
                          (length(x1) + length(x2) - 2))
        d <- if (pooled_sd > 0) (m1 - m2) / pooled_sd else NA_real_

        # --- Power (%) ---
        power_val <- tryCatch({
            if (has_pwr && !is.na(d) && abs(d) > 0) {
                pw <- pwr::pwr.t2n.test(
                    n1 = length(x1), n2 = length(x2),
                    d = abs(d), sig.level = 0.05,
                    alternative = "two.sided"
                )
                pw$power * 100
            } else if (!is.na(d) && abs(d) > 0) {
                # Approximate power using normal approximation
                n_harm <- 2 / (1 / length(x1) + 1 / length(x2))
                se <- sqrt(2 / n_harm)
                ncp <- abs(d) / se
                crit <- stats::qnorm(1 - 0.05 / 2)
                pw_approx <- stats::pnorm(ncp - crit) +
                             stats::pnorm(-ncp - crit)
                pw_approx * 100
            } else {
                NA_real_
            }
        }, error = function(e) NA_real_)

        # --- log2 fold change (from pre-scaling log matrix) ---
        l1 <- mat_log[i, idx1]
        l2 <- mat_log[i, idx2]
        l1 <- l1[!is.na(l1)]
        l2 <- l2[!is.na(l2)]
        log2fc <- if (length(l1) && length(l2)) mean(l1) - mean(l2) else NA_real_

        data.frame(
            feature_id = rownames(mat)[i],
            p_value    = pval,
            AUC        = auc_val,
            cohens_d   = round(d, 4),
            power_pct  = round(power_val, 2),
            log2FC     = round(log2fc, 4),
            stringsAsFactors = FALSE
        )
    })

    result_df <- do.call(rbind, results)
    result_df <- result_df[order(result_df$p_value, na.last = TRUE), ]
    rownames(result_df) <- NULL

    # Add full display names
    result_df$display_name <- format_lipid_display_name(result_df$feature_id)

    # Flag top 20 as candidate biomarkers
    result_df$is_candidate <- seq_len(nrow(result_df)) <= 20

    result_df
}


# ==== BIOMARKER TABLE PLOT ==================================================

#' Plot a publication-quality biomarker table
#'
#' Renders the top N biomarker candidates as a heatmap-style table using
#' ggplot2::geom_tile + geom_text.
#'
#' @param biomarker_df data.frame from compute_biomarker_metrics().
#' @param top_n        Number of candidates to display (default 20).
#' @return ggplot object or NULL.
plot_biomarker_table <- function(biomarker_df, top_n = 20) {
    if (is.null(biomarker_df) || nrow(biomarker_df) == 0) return(NULL)

    top_df <- utils::head(biomarker_df, top_n)
    top_df$rank <- seq_len(nrow(top_df))

    # Columns to display
    display_cols <- c("p_value", "AUC", "cohens_d", "power_pct", "log2FC")

    # Use display_name if available, otherwise feature_id
    label_col <- if ("display_name" %in% colnames(top_df)) "display_name" else "feature_id"

    # Build long-form data for tiles
    long <- do.call(rbind, lapply(display_cols, function(col) {
        data.frame(
            lipid      = top_df[[label_col]],
            rank       = top_df$rank,
            metric     = col,
            value      = top_df[[col]],
            stringsAsFactors = FALSE
        )
    }))

    # Normalise values to [0, 1] per metric for fill colour
    long$norm_val <- ave(long$value, long$metric, FUN = function(v) {
        rng <- range(v, na.rm = TRUE)
        if (diff(rng) == 0) rep(0.5, length(v))
        else (v - rng[1]) / diff(rng)
    })

    # For p_value, invert so small p = intense colour
    long$norm_val[long$metric == "p_value"] <-
        1 - long$norm_val[long$metric == "p_value"]

    # Format display text
    long$label <- ifelse(is.na(long$value), "NA",
        ifelse(long$metric == "p_value", formatC(long$value, format = "e", digits = 2),
        ifelse(long$metric == "power_pct", paste0(round(long$value, 1), "%"),
               round(long$value, 3))))

    long$lipid  <- factor(long$lipid, levels = rev(top_df[[label_col]]))
    long$metric <- factor(long$metric, levels = display_cols)

    ggplot2::ggplot(long, ggplot2::aes(x = metric, y = lipid,
                                        fill = norm_val)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.5) +
        ggplot2::geom_text(ggplot2::aes(label = label), size = 3) +
        ggplot2::scale_fill_gradient(low = "white", high = "#2166AC",
                                      na.value = "grey90", guide = "none") +
        ggplot2::scale_x_discrete(position = "top") +
        ggplot2::labs(
            title = paste0("Top ", top_n, " Biomarker Candidates"),
            x = NULL, y = NULL
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title    = ggplot2::element_text(face = "bold", size = 13,
                                                   hjust = 0.5),
            axis.text.x   = ggplot2::element_text(face = "bold", size = 10),
            axis.text.y   = ggplot2::element_text(size = 9),
            panel.grid     = ggplot2::element_blank()
        )
}


# ==== PER-FEATURE ROC CURVES ================================================

#' Compute ROC curves for top lipid biomarkers
#'
#' @param pre    List from preprocess_lipidomics() (pre contract).
#' @param config Full pipeline config.
#' @param top_n  Number of top features (by p-value) to include.
#' @return list(roc_data, auc_values) or NULL.
compute_roc_curves <- function(pre, config, top_n = 10) {
    cfg <- config$modes$lipidomics

    condition_col <- cfg$de$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col    <- cfg$effects$samples %||% "sample_id"

    mat  <- pre$expr_work
    meta <- pre$meta

    meta <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- factor(meta[[condition_col]])

    lvls <- levels(condition)
    if (length(lvls) != 2) {
        message("ROC curves: need exactly 2 condition levels — skipping")
        return(NULL)
    }

    # Binary response: 0 for first level, 1 for second
    resp <- as.integer(condition == lvls[2])

    has_pROC <- requireNamespace("pROC", quietly = TRUE)

    # Select top features from biomarker metrics (recompute p-values quickly)
    idx1 <- which(condition == lvls[1])
    idx2 <- which(condition == lvls[2])

    pvals <- vapply(seq_len(nrow(mat)), function(i) {
        tryCatch(
            stats::t.test(mat[i, idx1], mat[i, idx2])$p.value,
            error = function(e) NA_real_
        )
    }, numeric(1))

    top_idx <- utils::head(order(pvals, na.last = TRUE), top_n)
    top_features <- rownames(mat)[top_idx]

    roc_list <- list()
    auc_vals <- stats::setNames(numeric(length(top_features)), top_features)

    for (feat in top_features) {
        pred <- mat[feat, ]

        if (has_pROC) {
            roc_obj <- tryCatch(
                pROC::roc(resp, pred, quiet = TRUE),
                error = function(e) NULL
            )
            if (is.null(roc_obj)) next

            auc_vals[feat] <- as.numeric(pROC::auc(roc_obj))
            roc_list[[feat]] <- data.frame(
                feature_id = feat,
                fpr        = 1 - roc_obj$specificities,
                tpr        = roc_obj$sensitivities,
                stringsAsFactors = FALSE
            )
        } else {
            # Manual ROC via sorted thresholds
            ord <- order(pred)
            sorted_pred <- pred[ord]
            sorted_resp <- resp[ord]

            thresholds <- c(-Inf, sorted_pred, Inf)
            n_pos <- sum(resp == 1)
            n_neg <- sum(resp == 0)

            fpr_vec <- numeric(length(thresholds))
            tpr_vec <- numeric(length(thresholds))

            for (ti in seq_along(thresholds)) {
                predicted_pos <- pred >= thresholds[ti]
                tpr_vec[ti] <- sum(predicted_pos & resp == 1) / max(n_pos, 1)
                fpr_vec[ti] <- sum(predicted_pos & resp == 0) / max(n_neg, 1)
            }

            # Sort by fpr for proper curve
            ord2 <- order(fpr_vec, tpr_vec)
            fpr_vec <- fpr_vec[ord2]
            tpr_vec <- tpr_vec[ord2]

            # AUC via trapezoidal rule
            auc_val <- sum(diff(fpr_vec) * (tpr_vec[-1] + tpr_vec[-length(tpr_vec)]) / 2)
            auc_vals[feat] <- abs(auc_val)

            roc_list[[feat]] <- data.frame(
                feature_id = feat,
                fpr        = fpr_vec,
                tpr        = tpr_vec,
                stringsAsFactors = FALSE
            )
        }
    }

    if (length(roc_list) == 0) return(NULL)

    roc_data <- do.call(rbind, roc_list)
    rownames(roc_data) <- NULL

    list(roc_data = roc_data, auc_values = auc_vals)
}


#' Plot ROC curves for top lipid biomarkers
#'
#' @param roc_result List from compute_roc_curves().
#' @param top_n     Number of curves to plot (default 10).
#' @return ggplot object or NULL.
plot_roc_curves <- function(roc_result, top_n = 10) {
    if (is.null(roc_result)) return(NULL)

    roc_data   <- roc_result$roc_data
    auc_values <- roc_result$auc_values

    # Keep only top_n features
    top_feats <- names(utils::head(sort(auc_values, decreasing = TRUE), top_n))
    roc_data  <- roc_data[roc_data$feature_id %in% top_feats, ]

    # Build legend labels with full names and AUC
    display_names <- format_lipid_display_name(top_feats)
    legend_labels <- paste0(display_names, " (AUC=",
                             round(auc_values[top_feats], 3), ")")
    names(legend_labels) <- top_feats

    roc_data$label <- legend_labels[roc_data$feature_id]
    roc_data$label <- factor(roc_data$label, levels = legend_labels)

    ggplot2::ggplot(roc_data, ggplot2::aes(x = fpr, y = tpr,
                                            color = label)) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::geom_abline(slope = 1, intercept = 0,
                              linetype = "dashed", color = "grey50") +
        ggplot2::scale_color_brewer(palette = "Set3") +
        ggplot2::labs(
            title = "ROC Curves - Top Lipid Biomarkers",
            x = "False Positive Rate (1 - Specificity)",
            y = "True Positive Rate (Sensitivity)",
            color = "Feature"
        ) +
        ggplot2::coord_equal() +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title     = ggplot2::element_text(face = "bold", size = 13,
                                                    hjust = 0.5),
            legend.position = "right",
            legend.text     = ggplot2::element_text(size = 8)
        )
}


# ==== MULTIVARIATE ROC ======================================================

#' Compute multivariate ROC using Random Forest
#'
#' Fits a Random Forest on all features, extracts OOB predicted class
#' probabilities, and builds a single combined ROC curve.
#'
#' @param pre    List from preprocess_lipidomics() (pre contract).
#' @param config Full pipeline config.
#' @return list(roc_data, auc, model) or NULL.
compute_multivariate_roc <- function(pre, config) {
    cfg <- config$modes$lipidomics

    condition_col <- cfg$de$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col    <- cfg$effects$samples %||% "sample_id"

    mat  <- pre$expr_work
    meta <- pre$meta

    meta <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- factor(meta[[condition_col]])

    lvls <- levels(condition)
    if (length(lvls) != 2) {
        message("multivariate ROC: need exactly 2 condition levels — skipping")
        return(NULL)
    }

    resp <- as.integer(condition == lvls[2])

    X <- t(mat)
    # Impute NAs with column medians
    for (j in seq_len(ncol(X))) {
        nas <- is.na(X[, j])
        if (any(nas)) X[nas, j] <- stats::median(X[!nas, j], na.rm = TRUE)
    }

    # Remove zero-variance columns
    col_vars <- apply(X, 2, stats::var, na.rm = TRUE)
    keep <- col_vars > 0
    if (sum(keep) < 2) {
        message("multivariate ROC: fewer than 2 variable features — skipping")
        return(NULL)
    }
    X <- X[, keep, drop = FALSE]

    seed <- cfg$rf$seed %||% 1234
    n_trees <- cfg$rf$n_trees %||% 500

    # Fit RF and get OOB probabilities
    use_ranger <- requireNamespace("ranger", quietly = TRUE)

    if (use_ranger) {
        set.seed(seed)
        rf_fit <- tryCatch(
            ranger::ranger(
                x = X, y = factor(resp),
                num.trees = n_trees,
                probability = TRUE,
                seed = seed
            ),
            error = function(e) NULL
        )
        if (is.null(rf_fit)) {
            message("multivariate ROC: ranger fitting failed — skipping")
            return(NULL)
        }
        prob_pos <- rf_fit$predictions[, "1"]
    } else if (requireNamespace("randomForest", quietly = TRUE)) {
        set.seed(seed)
        orig_names <- colnames(X)
        safe_names <- paste0("F", seq_len(ncol(X)))
        colnames(X) <- safe_names

        rf_fit <- tryCatch(
            randomForest::randomForest(
                x = X, y = factor(resp),
                ntree = n_trees, importance = FALSE
            ),
            error = function(e) NULL
        )
        if (is.null(rf_fit)) {
            message("multivariate ROC: randomForest fitting failed — skipping")
            return(NULL)
        }
        # OOB votes as proportions
        votes <- rf_fit$votes
        prob_pos <- votes[, "1"]
        colnames(X) <- orig_names
    } else {
        message("multivariate ROC: no RF package available — skipping")
        return(NULL)
    }

    # Build ROC from OOB probabilities
    has_pROC <- requireNamespace("pROC", quietly = TRUE)

    if (has_pROC) {
        roc_obj <- tryCatch(
            pROC::roc(resp, prob_pos, quiet = TRUE, levels = c(0, 1),
                       direction = "<"),
            error = function(e) NULL
        )
        if (is.null(roc_obj)) return(NULL)

        auc_val <- as.numeric(pROC::auc(roc_obj))
        # Ensure AUC >= 0.5 (accounts for label/direction mismatches)
        if (auc_val < 0.5) auc_val <- 1 - auc_val
        roc_data <- data.frame(
            fpr = 1 - roc_obj$specificities,
            tpr = roc_obj$sensitivities,
            stringsAsFactors = FALSE
        )
    } else {
        # Manual ROC
        thresholds <- sort(unique(c(-Inf, prob_pos, Inf)), decreasing = TRUE)
        n_pos <- sum(resp == 1)
        n_neg <- sum(resp == 0)

        fpr_vec <- numeric(length(thresholds))
        tpr_vec <- numeric(length(thresholds))

        for (ti in seq_along(thresholds)) {
            predicted_pos <- prob_pos >= thresholds[ti]
            tpr_vec[ti] <- sum(predicted_pos & resp == 1) / max(n_pos, 1)
            fpr_vec[ti] <- sum(predicted_pos & resp == 0) / max(n_neg, 1)
        }

        ord <- order(fpr_vec, tpr_vec)
        fpr_vec <- fpr_vec[ord]
        tpr_vec <- tpr_vec[ord]

        auc_val <- abs(sum(diff(fpr_vec) *
                           (tpr_vec[-1] + tpr_vec[-length(tpr_vec)]) / 2))

        # Ensure AUC >= 0.5
        if (auc_val < 0.5) auc_val <- 1 - auc_val

        roc_data <- data.frame(fpr = fpr_vec, tpr = tpr_vec,
                                stringsAsFactors = FALSE)
    }

    list(roc_data = roc_data, auc = auc_val, model = rf_fit)
}


#' Plot multivariate ROC curve
#'
#' @param mv_roc_result List from compute_multivariate_roc().
#' @return ggplot object or NULL.
plot_multivariate_roc <- function(mv_roc_result) {
    if (is.null(mv_roc_result)) return(NULL)

    roc_data <- mv_roc_result$roc_data
    auc_val  <- mv_roc_result$auc

    auc_label <- paste0("AUC = ", round(auc_val, 3))

    ggplot2::ggplot(roc_data, ggplot2::aes(x = fpr, y = tpr)) +
        ggplot2::geom_line(color = "#D6604D", linewidth = 1.2) +
        ggplot2::geom_abline(slope = 1, intercept = 0,
                              linetype = "dashed", color = "grey50") +
        ggplot2::annotate("text", x = 0.6, y = 0.2, label = auc_label,
                           size = 5, fontface = "bold", color = "#D6604D") +
        ggplot2::coord_equal() +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::labs(
            title = "Multivariate ROC (Random Forest)",
            x = "False Positive Rate (1 - Specificity)",
            y = "True Positive Rate (Sensitivity)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13,
                                                hjust = 0.5)
        )
}
