# R/domain/lipidomics/08_qc_enhancements.R
#
# PCA and PLS-DA plot enhancements for lipidomics QC.
# Adds loading plots, scree plots, cross-validation diagnostics,
# confusion matrices, lipid class correlation, and PERMANOVA.
#
# Optional dependencies: mixOmics, vegan (loaded via tryCatch).


# ==== PCA ENHANCEMENTS =======================================================

#' PCA loading plot for top contributing features
#'
#' Plots loading vectors on the selected PCs and labels the top N features
#' with the highest absolute loading.
#'
#' @param pca_result  A prcomp (or equivalent) PCA result object.
#' @param pcs         Integer vector of length 2: which PCs to plot (default c(1,2)).
#' @param top_n       Number of top-loading features to label (default 20).
#' @param title       Plot title (auto-generated if NULL).
#' @return ggplot object.
plot_pca_loading <- function(pca_result, pcs = c(1, 2), top_n = 20,
                             title = NULL) {
    stopifnot(length(pcs) == 2)

    loadings <- pca_result$rotation
    if (is.null(loadings)) {
        warning("plot_pca_loading: no rotation matrix found in PCA result")
        return(NULL)
    }

    pc_labels <- paste0("PC", pcs)
    if (!all(pc_labels %in% colnames(loadings))) {
        warning("plot_pca_loading: requested PCs not found in rotation matrix")
        return(NULL)
    }

    var_expl <- summary(pca_result)$importance[2, pcs] * 100

    df <- data.frame(
        feature = rownames(loadings),
        x       = loadings[, pc_labels[1]],
        y       = loadings[, pc_labels[2]],
        stringsAsFactors = FALSE
    )
    df$abs_load <- sqrt(df$x^2 + df$y^2)

    top_n <- min(top_n, nrow(df))
    top_idx <- order(df$abs_load, decreasing = TRUE)[seq_len(top_n)]
    df$label <- ""
    df$label[top_idx] <- df$feature[top_idx]
    df$is_top <- FALSE
    df$is_top[top_idx] <- TRUE

    title <- title %||% paste0("PCA Loadings: ", pc_labels[1], " vs ", pc_labels[2])

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point(
            ggplot2::aes(color = is_top),
            size = 1.5, alpha = 0.6
        ) +
        ggplot2::scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick"),
                                     guide = "none") +
        ggplot2::geom_segment(
            data = df[df$is_top, ],
            ggplot2::aes(xend = 0, yend = 0),
            arrow = grid::arrow(length = grid::unit(0.15, "cm")),
            color = "firebrick", alpha = 0.5
        ) +
        ggplot2::geom_text(
            ggplot2::aes(label = label),
            size = 2.8, hjust = -0.1, vjust = -0.3, check_overlap = TRUE
        ) +
        ggplot2::labs(
            title = title,
            x = sprintf("%s (%.1f%%)", pc_labels[1], var_expl[1]),
            y = sprintf("%s (%.1f%%)", pc_labels[2], var_expl[2])
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5)
        )

    p
}


#' PCA scree plot
#'
#' Bar plot of per-PC variance explained with cumulative variance line overlay.
#'
#' @param pca_result  A prcomp (or equivalent) PCA result object.
#' @param max_pcs     Maximum number of PCs to display (default 10).
#' @param title       Plot title.
#' @return ggplot object.
plot_pca_scree <- function(pca_result, max_pcs = 10,
                           title = "PCA Scree Plot") {
    importance <- summary(pca_result)$importance
    var_pct <- importance[2, ] * 100
    cum_pct <- importance[3, ] * 100

    n_show <- min(max_pcs, length(var_pct))
    df <- data.frame(
        PC         = factor(seq_len(n_show), levels = seq_len(n_show)),
        variance   = var_pct[seq_len(n_show)],
        cumulative = cum_pct[seq_len(n_show)]
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = PC, y = variance)) +
        ggplot2::geom_col(fill = "steelblue", width = 0.6, alpha = 0.85) +
        ggplot2::geom_line(
            ggplot2::aes(x = as.numeric(PC), y = cumulative),
            color = "firebrick", linewidth = 0.8, group = 1
        ) +
        ggplot2::geom_point(
            ggplot2::aes(x = as.numeric(PC), y = cumulative),
            color = "firebrick", size = 2
        ) +
        ggplot2::scale_y_continuous(
            name = "Variance Explained (%)",
            sec.axis = ggplot2::sec_axis(~ ., name = "Cumulative Variance (%)")
        ) +
        ggplot2::labs(title = title, x = "Principal Component") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5)
        )

    p
}


# ==== PLS-DA ENHANCEMENTS ====================================================

#' PLS-DA cross-validation
#'
#' Runs cross-validation on a PLS-DA model using mixOmics::perf().
#' Falls back gracefully if mixOmics is unavailable.
#'
#' @param X           Numeric matrix (samples x features).
#' @param condition   Factor vector of group labels.
#' @param ncomp_max   Maximum number of components to evaluate (default 5).
#' @param n_folds     Number of CV folds; NULL = auto (LOO if n<10, else 5-fold).
#' @return list(accuracy, R2Y, Q2Y, optimal_ncomp, perf_result) or NULL.
compute_plsda_cv <- function(X, condition, ncomp_max = 5, n_folds = NULL) {
    if (!requireNamespace("mixOmics", quietly = TRUE)) {
        message("compute_plsda_cv: mixOmics not available -- skipping")
        return(NULL)
    }

    n <- nrow(X)
    ncomp_max <- min(ncomp_max, ncol(X), n - 1)
    if (ncomp_max < 1) {
        message("compute_plsda_cv: insufficient dimensions -- skipping")
        return(NULL)
    }

    # Determine fold strategy
    validation <- "Mfold"
    folds <- n_folds %||% if (n < 10) n else 5L
    if (is.null(n_folds) && n < 10) validation <- "loo"

    tryCatch({
        model <- mixOmics::plsda(X, condition, ncomp = ncomp_max)
        perf_res <- if (validation == "loo") {
            mixOmics::perf(model, validation = "loo", progressBar = FALSE)
        } else {
            mixOmics::perf(model, validation = "Mfold", folds = folds,
                           nrepeat = 10, progressBar = FALSE)
        }

        # Extract overall error rates (BER)
        err <- perf_res$error.rate$BER
        accuracy <- 1 - err[, "max.dist"]

        # Optimal ncomp based on max.dist BER
        optimal_ncomp <- perf_res$choice.ncomp["BER", "max.dist"]

        list(
            accuracy      = accuracy,
            optimal_ncomp = optimal_ncomp,
            perf_result   = perf_res
        )
    }, error = function(e) {
        warning("compute_plsda_cv failed: ", conditionMessage(e))
        NULL
    })
}


#' Plot PLS-DA cross-validation results
#'
#' Line plot showing accuracy (1 - BER) across components.
#'
#' @param cv_result  List from compute_plsda_cv().
#' @param title      Plot title.
#' @return ggplot object or NULL.
plot_plsda_cv <- function(cv_result, title = "PLS-DA Cross-Validation") {
    if (is.null(cv_result)) return(NULL)

    accuracy <- cv_result$accuracy
    ncomp_vals <- seq_along(accuracy)

    df <- data.frame(
        ncomp    = ncomp_vals,
        Accuracy = accuracy,
        stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = ncomp, y = Accuracy)) +
        ggplot2::geom_line(color = "steelblue", linewidth = 0.9) +
        ggplot2::geom_point(color = "steelblue", size = 2.5) +
        ggplot2::scale_x_continuous(breaks = ncomp_vals) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(
            title = title,
            x = "Number of Components",
            y = "Accuracy (1 - BER)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5)
        )

    if (!is.null(cv_result$optimal_ncomp)) {
        p <- p + ggplot2::geom_vline(
            xintercept = cv_result$optimal_ncomp,
            linetype = "dashed", color = "grey50"
        )
    }

    p
}


#' Compute PLS-DA confusion matrix
#'
#' Derives a confusion matrix from PLS-DA predictions (max.dist).
#'
#' @param plsda_res  List from run_lipidomics_plsda() (must contain model, X, condition).
#' @return list(confusion_matrix, accuracy, per_class_accuracy) or NULL.
compute_plsda_confusion <- function(plsda_res) {
    if (is.null(plsda_res)) return(NULL)

    if (!requireNamespace("mixOmics", quietly = TRUE)) {
        message("compute_plsda_confusion: mixOmics not available -- skipping")
        return(NULL)
    }

    tryCatch({
        model     <- plsda_res$model
        condition <- plsda_res$condition

        # Predict using training data (resubstitution).
        # mixOmics doesn't export `predict`; it ships predict.mixo_pls as an S3
        # method. Pull it from the namespace explicitly so dispatch works even
        # when the package is only loaded via ::, not library().
        predict_method <- get("predict.mixo_pls", envir = asNamespace("mixOmics"))
        pred <- predict_method(model, plsda_res$X)
        ncomp <- model$ncomp[1]
        predicted <- pred$class$max.dist[, ncomp]

        cm <- table(True = condition, Predicted = factor(predicted, levels = levels(condition)))

        overall_acc <- sum(diag(cm)) / sum(cm)
        per_class <- diag(cm) / rowSums(cm)
        per_class[is.nan(per_class)] <- 0

        list(
            confusion_matrix   = cm,
            accuracy           = overall_acc,
            per_class_accuracy = per_class
        )
    }, error = function(e) {
        warning("compute_plsda_confusion failed: ", conditionMessage(e))
        NULL
    })
}


#' Plot PLS-DA confusion matrix as heatmap
#'
#' @param confusion_result  List from compute_plsda_confusion().
#' @param title             Plot title.
#' @return ggplot object or NULL.
plot_plsda_confusion <- function(confusion_result,
                                  title = "PLS-DA Confusion Matrix") {
    if (is.null(confusion_result)) return(NULL)

    cm <- confusion_result$confusion_matrix
    df <- as.data.frame(cm, stringsAsFactors = FALSE)
    colnames(df) <- c("True", "Predicted", "Count")

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Predicted, y = True, fill = Count)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.8) +
        ggplot2::geom_text(ggplot2::aes(label = Count), size = 5, fontface = "bold") +
        ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
        ggplot2::labs(
            title = title,
            subtitle = sprintf("Overall Accuracy: %.1f%%", confusion_result$accuracy * 100),
            x = "Predicted Class", y = "True Class"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title    = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5),
            panel.grid    = ggplot2::element_blank()
        )

    p
}


# ==== LIPID CLASS CORRELATION ================================================

#' Lipid class correlation heatmap
#'
#' Aggregates intensities per lipid class and computes a correlation matrix,
#' then renders a heatmap.
#'
#' @param pre     List from preprocess_lipidomics() (pre contract).
#' @param config  Full pipeline config.
#' @param method  Correlation method: "pearson", "spearman", or "kendall".
#' @return ggplot object or NULL.
plot_lipid_class_correlation <- function(pre, config, method = "pearson") {
    row_data <- pre$row_data
    mat <- pre$expr_filt

    if (is.null(row_data) || !"lipid_class" %in% colnames(row_data)) {
        warning("plot_lipid_class_correlation: lipid_class column not found")
        return(NULL)
    }

    classes <- as.character(row_data$lipid_class)
    unique_classes <- sort(unique(classes))

    if (length(unique_classes) < 2) {
        message("plot_lipid_class_correlation: fewer than 2 classes -- skipping")
        return(NULL)
    }

    # Aggregate: sum per lipid class per sample
    class_mat <- matrix(0, nrow = length(unique_classes), ncol = ncol(mat),
                         dimnames = list(unique_classes, colnames(mat)))
    for (cl in unique_classes) {
        idx <- which(classes == cl)
        if (length(idx) == 1) {
            class_mat[cl, ] <- mat[idx, ]
        } else {
            class_mat[cl, ] <- colSums(mat[idx, , drop = FALSE], na.rm = TRUE)
        }
    }

    cor_mat <- stats::cor(t(class_mat), method = method, use = "pairwise.complete.obs")

    # Convert to long format for ggplot
    df <- expand.grid(
        Class1 = rownames(cor_mat),
        Class2 = colnames(cor_mat),
        stringsAsFactors = FALSE
    )
    df$correlation <- as.vector(cor_mat)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Class1, y = Class2, fill = correlation)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.3) +
        ggplot2::geom_text(
            ggplot2::aes(label = sprintf("%.2f", correlation)),
            size = 2.5
        ) +
        ggplot2::scale_fill_gradient2(
            low = "#2166AC", mid = "white", high = "#B2182B",
            midpoint = 0, limits = c(-1, 1),
            name = method
        ) +
        ggplot2::labs(
            title = paste0("Lipid Class Correlation (", method, ")"),
            x = NULL, y = NULL
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title    = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
            axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
            axis.text.y   = ggplot2::element_text(size = 9),
            panel.grid     = ggplot2::element_blank()
        )

    p
}


# ==== PERMANOVA ==============================================================

#' Run PERMANOVA on lipidomics expression matrix
#'
#' Tests multivariate group differences using vegan::adonis2 with Euclidean
#' distance. Falls back gracefully if vegan is unavailable.
#'
#' @param pre     List from preprocess_lipidomics() (pre contract).
#' @param config  Full pipeline config.
#' @param n_perm  Number of permutations (default 999).
#' @return list(F_statistic, R2, p_value, full_result) or NULL.
run_permanova <- function(pre, config, n_perm = 999) {
    if (!requireNamespace("vegan", quietly = TRUE)) {
        message("run_permanova: vegan not available -- skipping")
        return(NULL)
    }

    cfg <- config$modes$lipidomics
    condition_col <- cfg$de$condition_column %||% cfg$effects$color %||% "sample_type"
    sample_col <- cfg$effects$samples %||% "sample_id"

    mat  <- pre$expr_work
    meta <- pre$meta

    meta <- meta[match(colnames(mat), meta[[sample_col]]), , drop = FALSE]
    condition <- factor(meta[[condition_col]])

    X <- t(mat)
    # Impute NAs with column medians
    for (j in seq_len(ncol(X))) {
        nas <- is.na(X[, j])
        if (any(nas)) X[nas, j] <- stats::median(X[!nas, j], na.rm = TRUE)
    }

    tryCatch({
        dist_mat <- stats::dist(X, method = "euclidean")
        perm_df <- data.frame(condition = condition)

        result <- vegan::adonis2(
            dist_mat ~ condition,
            data      = perm_df,
            permutations = n_perm,
            method    = "euclidean"
        )

        list(
            F_statistic = result$F[1],
            R2          = result$R2[1],
            p_value     = result$`Pr(>F)`[1],
            full_result = result
        )
    }, error = function(e) {
        warning("run_permanova failed: ", conditionMessage(e))
        NULL
    })
}
