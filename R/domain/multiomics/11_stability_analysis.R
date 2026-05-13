# =============================================================================
# Integration Stability Analysis
# =============================================================================
# Assess robustness and stability of integration results through resampling
#
# Key analyses:
# - Bootstrap resampling to assess feature selection stability
# - Leave-one-out cross-validation
# - K-fold cross-validation
# - Factor/loading confidence intervals
# - Cluster stability assessment
# =============================================================================

# -----------------------------------------------------------------------------
# Main Orchestrator Function
# -----------------------------------------------------------------------------

#' Run Stability Analysis
#'
#' Comprehensive stability assessment of integration results
#'
#' @param mae MultiAssayExperiment object
#' @param feature_data Feature-selected data for integration
#' @param integration_results Results from integration methods
#' @param config Pipeline configuration list
#' @param out_dir Output directory
#' @return List containing stability analysis results
run_stability_analysis <- function(mae, feature_data, integration_results, config, out_dir) {
    message("=== Running Stability Analysis ===")

    # Get config settings
    stability_config <- config$modes$multiomics$stability %||% list()

    if (!(stability_config$run_stability %||% TRUE)) {
        message("Stability analysis disabled in config. Skipping.")
        return(NULL)
    }

    # Set up output directories
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    plots_dir <- file.path(out_dir, "plots")
    dir.create(plots_dir, showWarnings = FALSE)

    # Get available methods
    available_methods <- names(integration_results)[!sapply(integration_results, is.null)]
    message("Methods to analyze: ", paste(available_methods, collapse = ", "))

    # 1. Bootstrap feature selection stability
    feature_stability <- NULL
    if (stability_config$bootstrap_features %||% TRUE) {
        n_boot <- stability_config$bootstrap_n %||% 100
        feature_stability <- run_bootstrap_feature_stability(
            feature_data = feature_data,
            integration_results = integration_results,
            n_bootstrap = n_boot,
            config = config,
            output_dir = out_dir
        )
    }

    # 2. Sample-level stability (leave-one-out)
    sample_stability <- NULL
    if (stability_config$loo_analysis %||% FALSE) {
        sample_stability <- run_leave_one_out_analysis(
            feature_data = feature_data,
            integration_results = integration_results,
            config = config,
            output_dir = out_dir
        )
    }

    # 3. K-fold cross-validation stability
    kfold_stability <- NULL
    if (stability_config$kfold_analysis %||% TRUE) {
        n_folds <- stability_config$n_folds %||% 5
        kfold_stability <- run_kfold_stability(
            feature_data = feature_data,
            integration_results = integration_results,
            n_folds = n_folds,
            config = config,
            output_dir = out_dir
        )
    }

    # 4. Factor loading confidence intervals
    loading_ci <- NULL
    if (stability_config$loading_ci %||% TRUE) {
        n_boot <- stability_config$bootstrap_n %||% 100
        loading_ci <- compute_loading_confidence_intervals(
            feature_data = feature_data,
            integration_results = integration_results,
            n_bootstrap = n_boot,
            config = config,
            output_dir = out_dir
        )
    }

    # 5. Cluster stability
    cluster_stability <- NULL
    if (stability_config$cluster_stability %||% TRUE) {
        cluster_stability <- assess_cluster_stability(
            feature_data = feature_data,
            integration_results = integration_results,
            n_bootstrap = stability_config$bootstrap_n %||% 100,
            config = config,
            output_dir = out_dir
        )
    }

    # Generate visualizations
    figures <- generate_stability_plots(
        feature_stability = feature_stability,
        sample_stability = sample_stability,
        kfold_stability = kfold_stability,
        loading_ci = loading_ci,
        cluster_stability = cluster_stability,
        config = config,
        plots_dir = plots_dir
    )

    # Compile results
    results <- list(
        feature_stability = feature_stability,
        sample_stability = sample_stability,
        kfold_stability = kfold_stability,
        loading_ci = loading_ci,
        cluster_stability = cluster_stability,
        figures = figures,
        summary = create_stability_summary(feature_stability, kfold_stability, cluster_stability)
    )

    # Save outputs
    save_stability_outputs(results, out_dir)

    message("Stability analysis complete!")
    return(results)
}

# -----------------------------------------------------------------------------
# Bootstrap Feature Selection Stability
# -----------------------------------------------------------------------------

#' Run Bootstrap Feature Selection Stability Analysis
#'
#' @param feature_data Feature-selected data
#' @param integration_results Integration results
#' @param n_bootstrap Number of bootstrap iterations
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return Feature stability results
run_bootstrap_feature_stability <- function(feature_data, integration_results,
                                             n_bootstrap = 100, config, output_dir) {
    message("Running bootstrap feature selection stability...")
    message("  Bootstrap iterations: ", n_bootstrap)

    stability_results <- list()

    # MOFA stability
    if (!is.null(integration_results$mofa)) {
        stability_results$mofa <- bootstrap_mofa_features(
            feature_data = feature_data,
            mofa_results = integration_results$mofa,
            n_bootstrap = n_bootstrap,
            config = config
        )
    }

    # DIABLO stability
    if (!is.null(integration_results$diablo)) {
        stability_results$diablo <- bootstrap_diablo_features(
            feature_data = feature_data,
            diablo_results = integration_results$diablo,
            n_bootstrap = n_bootstrap,
            config = config
        )
    }

    # Combine results
    if (length(stability_results) == 0) {
        message("  No methods available for bootstrap analysis.")
        return(NULL)
    }

    # Create combined stability data frame
    all_stability <- lapply(names(stability_results), function(method) {
        df <- stability_results[[method]]
        df$method <- method
        df
    })

    combined_df <- do.call(rbind, all_stability)

    # Save results
    write.csv(
        combined_df,
        file.path(output_dir, "feature_selection_stability.csv"),
        row.names = FALSE
    )

    return(list(
        per_method = stability_results,
        combined = combined_df
    ))
}

#' Bootstrap MOFA Feature Stability
bootstrap_mofa_features <- function(feature_data, mofa_results, n_bootstrap, config) {
    message("  Bootstrapping MOFA features...")

    # Use pre-extracted weights (survives targets serialization)
    original_weights <- mofa_results$weights
    if (is.null(original_weights) && !is.null(mofa_results$model)) {
        original_weights <- tryCatch(MOFA2::get_weights(mofa_results$model), error = function(e) NULL)
    }

    if (is.null(original_weights)) return(NULL)

    # Get sample indices
    samples <- colnames(feature_data$matrices[[1]])
    n_samples <- length(samples)

    # Track feature selection frequency
    all_features <- unique(unlist(lapply(original_weights, rownames)))
    feature_counts <- setNames(rep(0, length(all_features)), all_features)

    # Bootstrap iterations
    pb <- if (interactive()) txtProgressBar(min = 0, max = n_bootstrap, style = 3) else NULL

    for (b in 1:n_bootstrap) {
        if (!is.null(pb)) setTxtProgressBar(pb, b)

        # Bootstrap sample indices
        boot_idx <- sample(1:n_samples, replace = TRUE)
        boot_samples <- samples[boot_idx]

        # Subset data
        boot_data <- lapply(feature_data$matrices, function(mat) {
            mat[, boot_idx, drop = FALSE]
        })

        # Run MOFA on bootstrap sample
        tryCatch({
            # Create MOFA object
            mofa_object <- create_mofa_for_bootstrap(boot_data, config)

            if (!is.null(mofa_object)) {
                # Get top features from this bootstrap
                boot_weights <- MOFA2::get_weights(mofa_object)

                for (view in names(boot_weights)) {
                    w <- boot_weights[[view]]
                    # Top features by max absolute weight
                    top_features <- names(sort(apply(abs(w), 1, max), decreasing = TRUE))[1:min(100, nrow(w))]
                    feature_counts[top_features] <- feature_counts[top_features] + 1
                }
            }
        }, error = function(e) {
            # Skip failed bootstrap iterations
        })
    }

    if (!is.null(pb)) close(pb)

    # Calculate stability scores
    stability_df <- data.frame(
        feature = names(feature_counts),
        selection_count = as.numeric(feature_counts),
        stability = as.numeric(feature_counts) / n_bootstrap,
        stringsAsFactors = FALSE
    )

    stability_df <- stability_df[order(-stability_df$stability), ]

    return(stability_df)
}

#' Bootstrap DIABLO Feature Stability
bootstrap_diablo_features <- function(feature_data, diablo_results, n_bootstrap, config) {
    message("  Bootstrapping DIABLO features...")

    # Use pre-extracted loadings (survives targets serialization)
    original_loadings <- diablo_results$loadings %||% diablo_results$model$loadings

    if (is.null(original_loadings)) return(NULL)

    # Get sample indices
    samples <- colnames(feature_data$matrices[[1]])
    n_samples <- length(samples)

    # Track feature selection frequency
    all_features <- unique(unlist(lapply(original_loadings, rownames)))
    feature_counts <- setNames(rep(0, length(all_features)), all_features)

    # Get condition labels
    condition_col <- config$design$condition_column %||% "condition"
    Y <- feature_data$metadata[[condition_col]]

    # Bootstrap iterations
    pb <- if (interactive()) txtProgressBar(min = 0, max = n_bootstrap, style = 3) else NULL

    for (b in 1:n_bootstrap) {
        if (!is.null(pb)) setTxtProgressBar(pb, b)

        # Bootstrap sample indices (stratified by condition)
        boot_idx <- stratified_bootstrap_indices(Y)
        boot_samples <- samples[boot_idx]

        # Subset data
        boot_data <- lapply(feature_data$matrices, function(mat) {
            mat[, boot_idx, drop = FALSE]
        })
        boot_Y <- Y[boot_idx]

        # Run DIABLO on bootstrap sample
        tryCatch({
            diablo_object <- run_diablo_for_bootstrap(boot_data, boot_Y, config)

            if (!is.null(diablo_object)) {
                # Get selected features
                boot_loadings <- diablo_object$loadings

                for (view in names(boot_loadings)) {
                    l <- boot_loadings[[view]]
                    # Features with non-zero loadings
                    selected <- rownames(l)[rowSums(abs(l)) > 0]
                    feature_counts[selected] <- feature_counts[selected] + 1
                }
            }
        }, error = function(e) {
            # Skip failed bootstrap iterations
        })
    }

    if (!is.null(pb)) close(pb)

    # Calculate stability scores
    stability_df <- data.frame(
        feature = names(feature_counts),
        selection_count = as.numeric(feature_counts),
        stability = as.numeric(feature_counts) / n_bootstrap,
        stringsAsFactors = FALSE
    )

    stability_df <- stability_df[order(-stability_df$stability), ]

    return(stability_df)
}

#' Stratified Bootstrap Indices
stratified_bootstrap_indices <- function(Y) {
    unique_levels <- unique(Y)
    boot_idx <- c()

    for (level in unique_levels) {
        level_idx <- which(Y == level)
        boot_level_idx <- sample(level_idx, replace = TRUE)
        boot_idx <- c(boot_idx, boot_level_idx)
    }

    return(boot_idx)
}

#' Create MOFA Object for Bootstrap
create_mofa_for_bootstrap <- function(boot_data, config) {
    tryCatch({
        # Create MOFA object
        mofa_object <- MOFA2::create_mofa(boot_data)

        # Set options (simpler for bootstrap)
        data_opts <- MOFA2::get_default_data_options(mofa_object)
        model_opts <- MOFA2::get_default_model_options(mofa_object)

        mofa_config <- config$modes$multiomics$mofa %||% list()
        model_opts$num_factors <- mofa_config$n_factors %||% 5

        train_opts <- MOFA2::get_default_training_options(mofa_object)
        train_opts$maxiter <- 500
        train_opts$verbose <- FALSE
        train_opts$seed <- sample(1:10000, 1)

        mofa_object <- MOFA2::prepare_mofa(
            mofa_object,
            data_options = data_opts,
            model_options = model_opts,
            training_options = train_opts
        )

        # Run MOFA
        mofa_trained <- MOFA2::run_mofa(mofa_object, save_data = FALSE)

        return(mofa_trained)
    }, error = function(e) {
        return(NULL)
    })
}

#' Run DIABLO for Bootstrap
run_diablo_for_bootstrap <- function(boot_data, boot_Y, config) {
    tryCatch({
        # Prepare data as list
        X <- lapply(boot_data, function(mat) t(mat))

        # Simpler design for bootstrap
        design <- matrix(0.1, nrow = length(X), ncol = length(X))
        diag(design) <- 0

        # Keep variables setting
        keepX <- lapply(X, function(x) rep(min(50, ncol(x)), 2))

        # Run DIABLO
        result <- mixOmics::block.splsda(
            X = X,
            Y = boot_Y,
            ncomp = 2,
            keepX = keepX,
            design = design,
            max.iter = 100
        )

        return(result)
    }, error = function(e) {
        return(NULL)
    })
}

# -----------------------------------------------------------------------------
# Leave-One-Out Analysis
# -----------------------------------------------------------------------------

#' Run Leave-One-Out Stability Analysis
#'
#' @param feature_data Feature-selected data
#' @param integration_results Integration results
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return LOO stability results
run_leave_one_out_analysis <- function(feature_data, integration_results, config, output_dir) {
    message("Running leave-one-out analysis...")

    # Get samples
    samples <- colnames(feature_data$matrices[[1]])
    n_samples <- length(samples)

    message("  Processing ", n_samples, " LOO iterations...")

    # Track cluster assignments
    loo_clusters <- matrix(NA, nrow = n_samples, ncol = n_samples - 1)
    rownames(loo_clusters) <- samples

    # For each sample left out
    pb <- if (interactive()) txtProgressBar(min = 0, max = n_samples, style = 3) else NULL

    for (i in 1:n_samples) {
        if (!is.null(pb)) setTxtProgressBar(pb, i)

        # Leave out sample i
        loo_idx <- setdiff(1:n_samples, i)
        loo_data <- lapply(feature_data$matrices, function(mat) {
            mat[, loo_idx, drop = FALSE]
        })

        # Run quick integration
        tryCatch({
            # Simple PCA-based clustering for LOO
            combined <- do.call(rbind, lapply(loo_data, function(mat) {
                scale(t(mat))
            }))
            combined <- combined[, colSums(is.na(combined)) == 0]

            pca <- prcomp(combined, scale. = FALSE)
            factors <- pca$x[, 1:min(5, ncol(pca$x))]

            clusters <- kmeans(factors, centers = 2, nstart = 10)$cluster
            loo_clusters[i, ] <- clusters

        }, error = function(e) {
            # Skip failed iterations
        })
    }

    if (!is.null(pb)) close(pb)

    # Calculate stability metrics
    # How often are pairs of samples clustered together?
    co_cluster <- matrix(0, n_samples, n_samples, dimnames = list(samples, samples))

    for (i in 1:n_samples) {
        if (all(is.na(loo_clusters[i, ]))) next

        loo_samples <- samples[-i]
        clust <- loo_clusters[i, ]

        for (j in 1:(length(loo_samples) - 1)) {
            for (k in (j + 1):length(loo_samples)) {
                if (!is.na(clust[j]) && !is.na(clust[k]) && clust[j] == clust[k]) {
                    s1 <- loo_samples[j]
                    s2 <- loo_samples[k]
                    co_cluster[s1, s2] <- co_cluster[s1, s2] + 1
                    co_cluster[s2, s1] <- co_cluster[s2, s1] + 1
                }
            }
        }
    }

    co_cluster <- co_cluster / (n_samples - 1)
    diag(co_cluster) <- 1

    # Save results
    write.csv(as.data.frame(co_cluster), file.path(output_dir, "loo_cocluster_matrix.csv"))

    return(list(
        loo_clusters = loo_clusters,
        cocluster_matrix = co_cluster
    ))
}

# -----------------------------------------------------------------------------
# K-Fold Cross-Validation Stability
# -----------------------------------------------------------------------------

#' Run K-Fold Stability Analysis
#'
#' @param feature_data Feature-selected data
#' @param integration_results Integration results
#' @param n_folds Number of folds
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return K-fold stability results
run_kfold_stability <- function(feature_data, integration_results, n_folds = 5,
                                 config, output_dir) {
    message("Running ", n_folds, "-fold cross-validation stability...")

    # Get samples
    samples <- colnames(feature_data$matrices[[1]])
    n_samples <- length(samples)

    # Create folds
    fold_assignment <- create_balanced_folds(n_samples, n_folds)

    # Track feature importance across folds
    fold_importance <- list()

    for (fold in 1:n_folds) {
        message("  Processing fold ", fold, "/", n_folds)

        # Training set (all except this fold)
        train_idx <- which(fold_assignment != fold)
        train_data <- lapply(feature_data$matrices, function(mat) {
            mat[, train_idx, drop = FALSE]
        })

        # Run integration on training fold
        tryCatch({
            # Quick MOFA for stability assessment
            mofa_result <- create_mofa_for_bootstrap(train_data, config)

            if (!is.null(mofa_result)) {
                weights <- MOFA2::get_weights(mofa_result)

                fold_importance[[paste0("fold_", fold)]] <- lapply(weights, function(w) {
                    importance <- apply(abs(w), 1, max)
                    data.frame(
                        feature = names(importance),
                        importance = importance,
                        stringsAsFactors = FALSE
                    )
                })
            }
        }, error = function(e) {
            # Skip failed folds
        })
    }

    if (length(fold_importance) == 0) {
        message("  No successful folds.")
        return(NULL)
    }

    # Calculate consistency across folds
    # Get all features
    all_features <- unique(unlist(lapply(fold_importance, function(fi) {
        unlist(lapply(fi, function(x) x$feature))
    })))

    # Create importance matrix
    importance_matrix <- matrix(
        NA, nrow = length(all_features), ncol = length(fold_importance),
        dimnames = list(all_features, names(fold_importance))
    )

    for (fold in names(fold_importance)) {
        fi <- fold_importance[[fold]]
        combined_fi <- do.call(rbind, fi)

        for (i in 1:nrow(combined_fi)) {
            feature <- combined_fi$feature[i]
            importance_matrix[feature, fold] <- max(importance_matrix[feature, fold],
                                                     combined_fi$importance[i], na.rm = TRUE)
        }
    }

    # Calculate stability metrics
    stability_df <- data.frame(
        feature = rownames(importance_matrix),
        mean_importance = rowMeans(importance_matrix, na.rm = TRUE),
        sd_importance = apply(importance_matrix, 1, sd, na.rm = TRUE),
        cv_importance = apply(importance_matrix, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)),
        n_folds_present = rowSums(!is.na(importance_matrix)),
        stringsAsFactors = FALSE
    )

    stability_df <- stability_df[order(-stability_df$mean_importance), ]

    # Save results
    write.csv(stability_df, file.path(output_dir, "kfold_feature_stability.csv"), row.names = FALSE)

    return(list(
        fold_importance = fold_importance,
        importance_matrix = importance_matrix,
        stability = stability_df
    ))
}

#' Create Balanced Folds
create_balanced_folds <- function(n_samples, n_folds) {
    # Shuffle and assign to folds
    shuffled <- sample(1:n_samples)
    fold_assignment <- rep(1:n_folds, length.out = n_samples)
    fold_assignment <- fold_assignment[order(shuffled)]
    return(fold_assignment)
}

# -----------------------------------------------------------------------------
# Loading Confidence Intervals
# -----------------------------------------------------------------------------

#' Compute Loading Confidence Intervals
#'
#' @param feature_data Feature-selected data
#' @param integration_results Integration results
#' @param n_bootstrap Number of bootstrap iterations
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return Loading CI results
compute_loading_confidence_intervals <- function(feature_data, integration_results,
                                                   n_bootstrap = 100, config, output_dir) {
    message("Computing loading confidence intervals...")

    if (is.null(integration_results$mofa)) {
        message("  MOFA results not available for loading CI.")
        return(NULL)
    }

    # Get original loadings
    # Use pre-extracted weights
    original_weights <- integration_results$mofa$weights
    if (is.null(original_weights) && !is.null(integration_results$mofa$model)) {
        original_weights <- tryCatch(MOFA2::get_weights(integration_results$mofa$model), error = function(e) NULL)
    }

    if (is.null(original_weights)) return(NULL)

    # Focus on top features for CI estimation
    stability_config <- config$modes$multiomics$stability %||% list()
    top_n <- stability_config$ci_top_features %||% 50

    # Get samples
    samples <- colnames(feature_data$matrices[[1]])
    n_samples <- length(samples)

    # Bootstrap loadings
    boot_loadings <- list()

    pb <- if (interactive()) txtProgressBar(min = 0, max = n_bootstrap, style = 3) else NULL

    for (b in 1:n_bootstrap) {
        if (!is.null(pb)) setTxtProgressBar(pb, b)

        # Bootstrap sample indices
        boot_idx <- sample(1:n_samples, replace = TRUE)

        # Subset data
        boot_data <- lapply(feature_data$matrices, function(mat) {
            mat[, boot_idx, drop = FALSE]
        })

        # Run MOFA
        tryCatch({
            mofa_object <- create_mofa_for_bootstrap(boot_data, config)

            if (!is.null(mofa_object)) {
                boot_loadings[[b]] <- MOFA2::get_weights(mofa_object)
            }
        }, error = function(e) {
            # Skip failed iterations
        })
    }

    if (!is.null(pb)) close(pb)

    if (length(boot_loadings) < 10) {
        message("  Too few successful bootstrap iterations.")
        return(NULL)
    }

    # Calculate CIs for each view
    loading_cis <- list()

    for (view in names(original_weights)) {
        orig_w <- original_weights[[view]]

        # Get top features by importance
        importance <- apply(abs(orig_w), 1, max)
        top_features <- names(sort(importance, decreasing = TRUE))[1:min(top_n, length(importance))]

        # Collect bootstrap values for top features
        boot_values <- lapply(top_features, function(feature) {
            values <- sapply(boot_loadings, function(bl) {
                if (!is.null(bl) && view %in% names(bl) && feature %in% rownames(bl[[view]])) {
                    max(abs(bl[[view]][feature, ]))
                } else {
                    NA
                }
            })
            values[!is.na(values)]
        })
        names(boot_values) <- top_features

        # Calculate CIs
        ci_level <- stability_config$ci_level %||% 0.95
        alpha <- 1 - ci_level

        ci_df <- data.frame(
            feature = top_features,
            view = view,
            original = importance[top_features],
            boot_mean = sapply(boot_values, mean, na.rm = TRUE),
            boot_sd = sapply(boot_values, sd, na.rm = TRUE),
            ci_lower = sapply(boot_values, function(x) quantile(x, alpha / 2, na.rm = TRUE)),
            ci_upper = sapply(boot_values, function(x) quantile(x, 1 - alpha / 2, na.rm = TRUE)),
            n_boot = sapply(boot_values, length),
            stringsAsFactors = FALSE
        )

        loading_cis[[view]] <- ci_df
    }

    # Combine all views
    combined_ci <- do.call(rbind, loading_cis)

    # Save results
    write.csv(
        combined_ci,
        file.path(output_dir, "loading_confidence_intervals.csv"),
        row.names = FALSE
    )

    return(list(
        per_view = loading_cis,
        combined = combined_ci
    ))
}

# -----------------------------------------------------------------------------
# Cluster Stability Assessment
# -----------------------------------------------------------------------------

#' Assess Cluster Stability
#'
#' @param feature_data Feature-selected data
#' @param integration_results Integration results
#' @param n_bootstrap Number of bootstrap iterations
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return Cluster stability results
assess_cluster_stability <- function(feature_data, integration_results,
                                      n_bootstrap = 100, config, output_dir) {
    message("Assessing cluster stability...")

    # Get samples
    samples <- colnames(feature_data$matrices[[1]])
    n_samples <- length(samples)

    # Get original clusters
    original_clusters <- NULL

    if (!is.null(integration_results$mofa)) {
        tryCatch({
            # Use pre-extracted factors
            fdf <- integration_results$mofa$factors
            factors <- NULL
            if (!is.null(fdf)) {
                if (is.data.frame(fdf) && "sample_id" %in% colnames(fdf)) {
                    rn <- fdf$sample_id
                    factor_cols <- grep("^Factor|^V\\d", colnames(fdf), value = TRUE)
                    if (length(factor_cols) == 0) factor_cols <- setdiff(colnames(fdf), c("sample_id", "condition", "group"))
                    factors <- as.matrix(fdf[, factor_cols, drop = FALSE])
                    rownames(factors) <- rn
                } else if (is.matrix(fdf)) {
                    factors <- fdf
                }
            }
            # Fallback to model
            if (is.null(factors) && !is.null(integration_results$mofa$model)) {
                factors <- MOFA2::get_factors(integration_results$mofa$model)[[1]]
            }
            if (!is.null(factors) && nrow(factors) >= 5) {
                k <- 2
                original_clusters <- kmeans(factors, centers = k, nstart = 25)$cluster
                names(original_clusters) <- rownames(factors)
            }
        }, error = function(e) NULL)
    }

    if (is.null(original_clusters)) {
        message("  Could not obtain original clusters.")
        return(NULL)
    }

    # Bootstrap cluster assignments
    k <- length(unique(original_clusters))

    # Co-clustering matrix
    co_cluster <- matrix(0, n_samples, n_samples, dimnames = list(samples, samples))

    pb <- if (interactive()) txtProgressBar(min = 0, max = n_bootstrap, style = 3) else NULL

    for (b in 1:n_bootstrap) {
        if (!is.null(pb)) setTxtProgressBar(pb, b)

        # Bootstrap sample indices
        boot_idx <- sample(1:n_samples, replace = TRUE)
        boot_samples <- samples[boot_idx]

        # Subset data
        boot_data <- lapply(feature_data$matrices, function(mat) {
            mat[, boot_idx, drop = FALSE]
        })

        # Run quick clustering
        tryCatch({
            # Combined PCA
            combined <- do.call(rbind, lapply(boot_data, function(mat) t(mat)))
            combined <- combined[, colSums(is.na(combined)) == 0]

            if (ncol(combined) >= 5) {
                pca <- prcomp(combined, scale. = TRUE)
                factors <- pca$x[, 1:min(5, ncol(pca$x))]

                boot_clusters <- kmeans(factors, centers = k, nstart = 10)$cluster

                # Update co-clustering matrix
                for (i in 1:(n_samples - 1)) {
                    for (j in (i + 1):n_samples) {
                        s1 <- boot_samples[i]
                        s2 <- boot_samples[j]
                        if (boot_clusters[i] == boot_clusters[j]) {
                            co_cluster[s1, s2] <- co_cluster[s1, s2] + 1
                            co_cluster[s2, s1] <- co_cluster[s2, s1] + 1
                        }
                    }
                }
            }
        }, error = function(e) {
            # Skip failed iterations
        })
    }

    if (!is.null(pb)) close(pb)

    # Normalize co-clustering matrix
    co_cluster <- co_cluster / n_bootstrap
    diag(co_cluster) <- 1

    # Calculate per-sample stability (Jaccard index with original)
    sample_stability <- sapply(samples, function(s) {
        orig_same <- names(original_clusters)[original_clusters == original_clusters[s]]
        co_probs <- co_cluster[s, orig_same]
        mean(co_probs, na.rm = TRUE)
    })

    stability_df <- data.frame(
        sample = samples,
        original_cluster = original_clusters[samples],
        stability = sample_stability,
        stringsAsFactors = FALSE
    )

    # Cluster-level stability
    cluster_stability <- tapply(stability_df$stability, stability_df$original_cluster, mean)

    # Save results
    write.csv(stability_df, file.path(output_dir, "sample_cluster_stability.csv"), row.names = FALSE)
    write.csv(as.data.frame(co_cluster), file.path(output_dir, "cocluster_matrix.csv"))

    return(list(
        cocluster_matrix = co_cluster,
        sample_stability = stability_df,
        cluster_stability = cluster_stability,
        original_clusters = original_clusters
    ))
}

# -----------------------------------------------------------------------------
# Visualization Functions
# -----------------------------------------------------------------------------

#' Generate Stability Analysis Plots
generate_stability_plots <- function(feature_stability, sample_stability, kfold_stability,
                                      loading_ci, cluster_stability, config, plots_dir) {
    message("Generating stability analysis plots...")

    figures <- list()

    # 1. Feature stability barplot
    if (!is.null(feature_stability$combined)) {
        figures$feature_stability <- plot_feature_stability(
            stability_df = feature_stability$combined,
            plots_dir = plots_dir
        )
    }

    # 2. K-fold stability
    if (!is.null(kfold_stability)) {
        figures$kfold_stability <- plot_kfold_stability(
            kfold_results = kfold_stability,
            plots_dir = plots_dir
        )
    }

    # 3. Loading confidence intervals
    if (!is.null(loading_ci)) {
        figures$loading_ci <- plot_loading_ci(
            loading_ci = loading_ci,
            plots_dir = plots_dir
        )
    }

    # 4. Cluster stability heatmap
    if (!is.null(cluster_stability)) {
        figures$cluster_stability <- plot_cluster_stability(
            cluster_results = cluster_stability,
            plots_dir = plots_dir
        )
    }

    return(figures)
}

#' Plot Feature Stability
plot_feature_stability <- function(stability_df, plots_dir) {
    message("  Creating feature stability plot...")

    # Top 50 features by stability
    top_df <- head(stability_df[order(-stability_df$stability), ], 50)
    top_df$feature <- factor(top_df$feature, levels = rev(top_df$feature))

    p <- ggplot2::ggplot(top_df, ggplot2::aes(x = feature, y = stability, fill = method)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::coord_flip() +
        ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
        ggplot2::labs(
            title = "Feature Selection Stability",
            subtitle = "Proportion of bootstrap iterations where feature was selected",
            x = "Feature",
            y = "Stability (selection frequency)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))

    fig_path <- file.path(plots_dir, "feature_stability.png")
    ggplot2::ggsave(fig_path, p, width = 10, height = 12)

    return(fig_path)
}

#' Plot K-Fold Stability
plot_kfold_stability <- function(kfold_results, plots_dir) {
    message("  Creating k-fold stability plot...")

    stability_df <- kfold_results$stability

    # Top 50 features
    top_df <- head(stability_df, 50)
    top_df$feature <- factor(top_df$feature, levels = rev(top_df$feature))

    p <- ggplot2::ggplot(top_df, ggplot2::aes(x = feature, y = mean_importance)) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_importance - sd_importance,
                      ymax = mean_importance + sd_importance),
                  width = 0.3) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            title = "K-Fold Cross-Validation Feature Importance",
            subtitle = "Mean importance across folds with standard deviation",
            x = "Feature",
            y = "Mean Importance"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))

    fig_path <- file.path(plots_dir, "kfold_stability.png")
    ggplot2::ggsave(fig_path, p, width = 10, height = 12)

    return(fig_path)
}

#' Plot Loading Confidence Intervals
plot_loading_ci <- function(loading_ci, plots_dir) {
    message("  Creating loading CI plot...")

    combined_ci <- loading_ci$combined

    # Top 30 by original importance
    top_ci <- head(combined_ci[order(-combined_ci$original), ], 30)
    top_ci$feature <- factor(top_ci$feature, levels = rev(top_ci$feature))

    p <- ggplot2::ggplot(top_ci, ggplot2::aes(x = feature, y = original, color = view)) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), width = 0.3) +
        ggplot2::coord_flip() +
        ggplot2::labs(
            title = "Loading Confidence Intervals (Bootstrap)",
            x = "Feature",
            y = "Loading (max absolute)",
            color = "View"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))

    fig_path <- file.path(plots_dir, "loading_confidence_intervals.png")
    ggplot2::ggsave(fig_path, p, width = 10, height = 10)

    return(fig_path)
}

#' Plot Cluster Stability
plot_cluster_stability <- function(cluster_results, plots_dir) {
    message("  Creating cluster stability heatmap...")

    co_cluster <- cluster_results$cocluster_matrix
    original_clusters <- cluster_results$original_clusters

    # Order by cluster and stability
    sample_order <- names(sort(original_clusters))
    co_cluster_ordered <- co_cluster[sample_order, sample_order]

    # Create annotation
    cluster_anno <- data.frame(
        Cluster = factor(original_clusters[sample_order])
    )

    # Create heatmap using ComplexHeatmap if available
    fig_path <- file.path(plots_dir, "cluster_stability_heatmap.png")

    if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
        requireNamespace("circlize", quietly = TRUE)) {
        ht <- ComplexHeatmap::Heatmap(
            co_cluster_ordered,
            name = "Co-clustering\nProbability",
            col = circlize::colorRamp2(c(0, 0.5, 1), c("white", "orange", "red")),
            show_row_names = FALSE,
            show_column_names = FALSE,
            column_title = "Bootstrap Cluster Stability",
            top_annotation = ComplexHeatmap::HeatmapAnnotation(
                Cluster = cluster_anno$Cluster,
                col = list(Cluster = c("1" = "blue", "2" = "red", "3" = "green")[1:length(unique(cluster_anno$Cluster))])
            ),
            left_annotation = ComplexHeatmap::rowAnnotation(
                Cluster = cluster_anno$Cluster,
                col = list(Cluster = c("1" = "blue", "2" = "red", "3" = "green")[1:length(unique(cluster_anno$Cluster))])
            ),
          rect_gp = grid::gpar(col = NA)
        )

        png(fig_path, width = 10, height = 10, units = "in", res = 150)
        ComplexHeatmap::draw(ht)
        dev.off()
    } else {
        # Fallback to base R heatmap
        png(fig_path, width = 10, height = 10, units = "in", res = 150)
        heatmap(co_cluster_ordered, main = "Bootstrap Cluster Stability",
                col = colorRampPalette(c("white", "orange", "red"))(100))
        dev.off()
    }

    return(fig_path)
}

# -----------------------------------------------------------------------------
# Output Functions
# -----------------------------------------------------------------------------

#' Create Stability Summary
create_stability_summary <- function(feature_stability, kfold_stability, cluster_stability) {
    summary <- list()

    if (!is.null(feature_stability$combined)) {
        # Features with >80% stability
        stable_features <- sum(feature_stability$combined$stability > 0.8)
        summary$n_highly_stable_features <- stable_features
        summary$top_stable_features <- head(feature_stability$combined$feature[feature_stability$combined$stability > 0.8], 10)
    }

    if (!is.null(kfold_stability)) {
        # Features present in all folds
        summary$n_consistent_features <- sum(kfold_stability$stability$n_folds_present == max(kfold_stability$stability$n_folds_present))
    }

    if (!is.null(cluster_stability)) {
        summary$mean_cluster_stability <- mean(cluster_stability$sample_stability$stability)
        summary$per_cluster_stability <- cluster_stability$cluster_stability
    }

    return(summary)
}

#' Save Stability Outputs
save_stability_outputs <- function(results, output_dir) {
    message("Saving stability analysis outputs...")

    saveRDS(results, file.path(output_dir, "stability_results.rds"))

    message("  Outputs saved to: ", output_dir)
}
