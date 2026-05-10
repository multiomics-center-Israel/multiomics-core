# =============================================================================
# Integration Method Consensus Analysis
# =============================================================================
# Compare and synthesize results from multiple integration methods
# (MOFA, DIABLO, SNF)
#
# Key analyses:
# - Sample clustering comparison across methods (ARI, NMI)
# - Feature importance ranking comparison
# - Identification of robust patterns (consistent across methods)
# - Identification of method-specific patterns
# - Meta-integration of signals from all methods
# =============================================================================

# -----------------------------------------------------------------------------
# Main Orchestrator Function
# -----------------------------------------------------------------------------

#' Run Integration Consensus Analysis
#'
#' Compare and synthesize results from multiple integration methods
#'
#' @param integration_results List containing results from MOFA, DIABLO, SNF
#' @param mae MultiAssayExperiment object
#' @param config Pipeline configuration list
#' @param out_dir Output directory
#' @return List containing consensus analysis results
run_integration_consensus <- function(integration_results, mae, config, out_dir) {
    message("=== Running Integration Consensus Analysis ===")

    # Get config settings
    consensus_config <- config$modes$multiomics$consensus %||% list()

    if (!(consensus_config$compare_methods %||% TRUE)) {
        message("Consensus analysis disabled in config. Skipping.")
        return(NULL)
    }

    # Check which methods have results
    available_methods <- names(integration_results)[!sapply(integration_results, is.null)]

    if (length(available_methods) < 2) {
        message("Need at least 2 integration methods for consensus analysis.")
        message("Available methods: ", paste(available_methods, collapse = ", "))
        return(NULL)
    }

    message("Comparing methods: ", paste(available_methods, collapse = ", "))

    # Set up output directories
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    plots_dir <- file.path(out_dir, "plots")
    dir.create(plots_dir, showWarnings = FALSE)

    # Convert MAE to legacy format for compatibility
    mae_data <- .mae_to_legacy_consensus(mae, config)

    # 1. Sample clustering comparison
    sample_comparison <- compare_sample_clustering(
        integration_results = integration_results,
        mae_data = mae_data,
        output_dir = out_dir
    )

    # 2. Feature importance comparison
    feature_comparison <- compare_feature_importance(
        integration_results = integration_results,
        output_dir = out_dir
    )

    # Resolve generic feature IDs (GENE_N, feature_N) to display names
    feature_name_map <- build_feature_name_map(mae)
    if (!is.null(feature_comparison) && length(feature_name_map) > 0) {
        feature_comparison <- resolve_consensus_feature_names(
            feature_comparison, feature_name_map, out_dir
        )
    }

    # 3. Identify robust patterns
    robust_patterns <- identify_robust_patterns(
        sample_comparison = sample_comparison,
        feature_comparison = feature_comparison,
        config = config,
        output_dir = out_dir
    )

    # 4. Identify method-specific patterns
    specific_patterns <- identify_method_specific_patterns(
        feature_comparison = feature_comparison,
        config = config,
        output_dir = out_dir
    )

    # 5. Meta-integration
    meta_integration <- NULL
    if (consensus_config$run_meta_integration %||% TRUE) {
        meta_integration <- run_meta_integration(
            integration_results = integration_results,
            mae_data = mae_data,
            config = config,
            output_dir = out_dir
        )
    }

    # Generate visualizations
    figures <- generate_consensus_plots(
        sample_comparison = sample_comparison,
        feature_comparison = feature_comparison,
        robust_patterns = robust_patterns,
        specific_patterns = specific_patterns,
        meta_integration = meta_integration,
        config = config,
        plots_dir = plots_dir
    )

    # Cross-method bucket enrichment (MOFA factor x DIABLO component story).
    # Reads the on-disk MOFA/DIABLO top-feature CSVs that the integration modules
    # already wrote, so it depends only on out_dir's parent layout.
    bucket_results <- tryCatch({
        mo_root <- dirname(dirname(out_dir))  # consensus/consensus -> multiomics
        if (!dir.exists(file.path(mo_root, "integration"))) mo_root <- dirname(out_dir)
        run_cross_method_buckets(mo_root = mo_root, config = config,
                                 out_dir = file.path(out_dir, "buckets"))
    }, error = function(e) {
        message("  Cross-method bucket analysis failed: ", e$message)
        NULL
    })

    # Compile results
    results <- list(
        available_methods = available_methods,
        sample_comparison = sample_comparison,
        feature_comparison = feature_comparison,
        robust_patterns = robust_patterns,
        method_specific = specific_patterns,
        meta_integration = meta_integration,
        bucket_results = bucket_results,
        figures = figures,
        summary = create_consensus_summary(sample_comparison, feature_comparison, robust_patterns)
    )

    # Save outputs
    save_consensus_outputs(results, out_dir)

    message("Integration consensus analysis complete!")
    return(results)
}

# -----------------------------------------------------------------------------
# MAE to Legacy Adapter
# -----------------------------------------------------------------------------

#' Convert MAE to legacy format for consensus analysis
#' @param mae MultiAssayExperiment object
#' @return Legacy-format list
.mae_to_legacy_consensus <- function(mae, config = NULL) {
    metadata <- as.data.frame(MultiAssayExperiment::colData(mae))
    common_samples <- colnames(mae)

    condition_col <- config$design$condition_column %||% "condition"

    list(
        mae = mae,
        metadata = metadata,
        common_samples = common_samples,
        config = list(design = list(condition_column = condition_col))
    )
}

# -----------------------------------------------------------------------------
# Sample Clustering Comparison
# -----------------------------------------------------------------------------

#' Compare Sample Clustering Across Methods
#'
#' @param integration_results Integration results list
#' @param mae_data MAE data object
#' @param output_dir Output directory
#' @return List with clustering comparison metrics
compare_sample_clustering <- function(integration_results, mae_data, output_dir) {
    message("Comparing sample clustering across methods...")

    # Check for aricode package
    if (!requireNamespace("aricode", quietly = TRUE)) {
        message("WARNING: aricode package not available. Using base R implementation.")
    }

    # Get sample clusters from each method
    sample_clusters <- list()

    # MOFA clusters (from factor scores)
    if (!is.null(integration_results$mofa)) {
        sample_clusters$mofa <- get_mofa_clusters(integration_results$mofa)
    }

    # DIABLO clusters (from component scores)
    if (!is.null(integration_results$diablo)) {
        sample_clusters$diablo <- get_diablo_clusters(integration_results$diablo)
    }

    # SNF clusters (from fused network)
    if (!is.null(integration_results$snf)) {
        sample_clusters$snf <- get_snf_clusters(integration_results$snf)
    }

    if (length(sample_clusters) < 2) {
        message("  Not enough methods with clustering results.")
        return(NULL)
    }

    # Align samples across methods
    common_samples <- Reduce(intersect, lapply(sample_clusters, function(x) names(x)))
    message("  Common samples across methods: ", length(common_samples))

    if (length(common_samples) < 5) {
        message("  Too few common samples for comparison.")
        return(NULL)
    }

    # Align clusters to common samples
    aligned_clusters <- lapply(sample_clusters, function(x) x[common_samples])

    # Compute pairwise agreement metrics
    methods <- names(aligned_clusters)
    n_methods <- length(methods)

    ari_matrix <- matrix(NA, n_methods, n_methods, dimnames = list(methods, methods))
    nmi_matrix <- matrix(NA, n_methods, n_methods, dimnames = list(methods, methods))

    for (i in 1:(n_methods - 1)) {
        for (j in (i + 1):n_methods) {
            clust1 <- aligned_clusters[[i]]
            clust2 <- aligned_clusters[[j]]

            # Adjusted Rand Index
            ari_matrix[i, j] <- compute_ari(clust1, clust2)
            ari_matrix[j, i] <- ari_matrix[i, j]

            # Normalized Mutual Information
            nmi_matrix[i, j] <- compute_nmi(clust1, clust2)
            nmi_matrix[j, i] <- nmi_matrix[i, j]
        }
    }

    diag(ari_matrix) <- 1
    diag(nmi_matrix) <- 1

    # Compare with ground truth (condition labels)
    ground_truth <- NULL
    condition_ari <- NULL
    condition_nmi <- NULL

    if (!is.null(mae_data$metadata)) {
        condition_col <- mae_data$config$design$condition_column %||% "condition"
        if (condition_col %in% colnames(mae_data$metadata)) {
            ground_truth <- mae_data$metadata[[condition_col]][match(common_samples, rownames(mae_data$metadata))]

            if (!all(is.na(ground_truth))) {
                condition_ari <- sapply(aligned_clusters, function(clust) {
                    compute_ari(clust, ground_truth)
                })

                condition_nmi <- sapply(aligned_clusters, function(clust) {
                    compute_nmi(clust, ground_truth)
                })
            }
        }
    }

    # Compute consensus clustering
    consensus_clusters <- compute_consensus_clustering(aligned_clusters)

    # Save results
    write.csv(as.data.frame(ari_matrix), file.path(output_dir, "method_ari_matrix.csv"))
    write.csv(as.data.frame(nmi_matrix), file.path(output_dir, "method_nmi_matrix.csv"))

    cluster_df <- data.frame(
        sample_id = common_samples,
        as.data.frame(aligned_clusters),
        consensus = consensus_clusters,
        stringsAsFactors = FALSE
    )
    write.csv(cluster_df, file.path(output_dir, "sample_clusters_comparison.csv"), row.names = FALSE)

    return(list(
        clusters = aligned_clusters,
        consensus = consensus_clusters,
        ari_matrix = ari_matrix,
        nmi_matrix = nmi_matrix,
        condition_ari = condition_ari,
        condition_nmi = condition_nmi,
        common_samples = common_samples
    ))
}

#' Get MOFA Clusters
get_mofa_clusters <- function(mofa_results) {
    tryCatch({
        # Use pre-extracted factors (survives targets serialization)
        factors <- NULL
        if (!is.null(mofa_results$factors)) {
            fdf <- mofa_results$factors
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

        # Fallback: try extracting from model
        if (is.null(factors) && !is.null(mofa_results$model)) {
            factors <- MOFA2::get_factors(mofa_results$model)[[1]]
        }

        if (is.null(factors) || nrow(factors) < 5) return(NULL)

        # Cluster based on first few factors
        n_factors <- min(5, ncol(factors))
        factor_mat <- factors[, 1:n_factors, drop = FALSE]

        best_k <- choose_optimal_k(factor_mat)
        clusters <- kmeans(factor_mat, centers = best_k, nstart = 25)$cluster

        setNames(clusters, rownames(factors))
    }, error = function(e) {
        message("  Could not extract MOFA clusters: ", e$message)
        NULL
    })
}

#' Get DIABLO Clusters
get_diablo_clusters <- function(diablo_results) {
    tryCatch({
        # Use pre-extracted sample_scores (survives targets serialization)
        variates <- diablo_results$sample_scores %||% diablo_results$model$variates
        if (is.null(variates)) return(NULL)

        # Exclude Y component if present
        variates <- variates[setdiff(names(variates), "Y")]
        if (length(variates) == 0) return(NULL)

        # Combine across views
        combined <- do.call(cbind, variates)

        # Cluster
        n_comp <- min(5, ncol(combined))
        factor_mat <- combined[, 1:n_comp, drop = FALSE]

        best_k <- choose_optimal_k(factor_mat)
        clusters <- kmeans(factor_mat, centers = best_k, nstart = 25)$cluster

        setNames(clusters, rownames(combined))
    }, error = function(e) {
        message("  Could not extract DIABLO clusters: ", e$message)
        NULL
    })
}

#' Get SNF Clusters
get_snf_clusters <- function(snf_results) {
    if (is.null(snf_results)) return(NULL)

    tryCatch({
        # SNF already provides clusters
        if (!is.null(snf_results$clusters)) {
            return(snf_results$clusters)
        }

        # Otherwise cluster from affinity matrix
        W <- snf_results$fused_network

        # Spectral clustering
        best_k <- snf_results$optimal_k %||% choose_optimal_k_spectral(W)
        clusters <- SNFtool::spectralClustering(W, best_k)

        setNames(clusters, rownames(W))
    }, error = function(e) {
        message("  Could not extract SNF clusters: ", e$message)
        NULL
    })
}

#' Choose Optimal K for Clustering
choose_optimal_k <- function(data_matrix, max_k = 10) {
    n <- nrow(data_matrix)
    max_k <- min(max_k, n - 1)

    if (max_k < 2) return(2)

    # Silhouette method
    sil_scores <- sapply(2:max_k, function(k) {
        tryCatch({
            clust <- kmeans(data_matrix, centers = k, nstart = 10)
            sil <- cluster::silhouette(clust$cluster, dist(data_matrix))
            mean(sil[, 3])
        }, error = function(e) NA)
    })

    optimal_k <- which.max(sil_scores) + 1
    if (is.na(optimal_k) || length(optimal_k) == 0) optimal_k <- 2

    return(optimal_k)
}

#' Choose Optimal K for Spectral Clustering
choose_optimal_k_spectral <- function(W, max_k = 10) {
    n <- nrow(W)
    max_k <- min(max_k, n - 1)

    if (max_k < 2) return(2)

    # Eigengap heuristic
    L <- diag(rowSums(W)) - W
    D_inv_sqrt <- diag(1 / sqrt(rowSums(W) + 1e-10))
    L_norm <- D_inv_sqrt %*% L %*% D_inv_sqrt

    eigenvalues <- sort(eigen(L_norm, only.values = TRUE)$values)
    gaps <- diff(eigenvalues[1:min(max_k + 1, length(eigenvalues))])

    optimal_k <- which.max(gaps) + 1
    if (is.na(optimal_k) || optimal_k < 2) optimal_k <- 2

    return(optimal_k)
}

#' Compute Adjusted Rand Index
compute_ari <- function(clust1, clust2) {
    if (requireNamespace("aricode", quietly = TRUE)) {
        return(aricode::ARI(clust1, clust2))
    }

    # Base R implementation
    tab <- table(clust1, clust2)
    n <- length(clust1)

    sum_comb <- sum(choose(tab, 2))
    sum_rows <- sum(choose(rowSums(tab), 2))
    sum_cols <- sum(choose(colSums(tab), 2))

    expected <- (sum_rows * sum_cols) / choose(n, 2)
    max_index <- (sum_rows + sum_cols) / 2

    ari <- (sum_comb - expected) / (max_index - expected)
    return(ari)
}

#' Compute Normalized Mutual Information
compute_nmi <- function(clust1, clust2) {
    if (requireNamespace("aricode", quietly = TRUE)) {
        return(aricode::NMI(clust1, clust2))
    }

    # Base R implementation
    tab <- table(clust1, clust2)
    n <- sum(tab)

    # Marginal probabilities
    p_i <- rowSums(tab) / n
    p_j <- colSums(tab) / n
    p_ij <- tab / n

    # Mutual information
    mi <- 0
    for (i in 1:nrow(tab)) {
        for (j in 1:ncol(tab)) {
            if (p_ij[i, j] > 0) {
                mi <- mi + p_ij[i, j] * log(p_ij[i, j] / (p_i[i] * p_j[j]))
            }
        }
    }

    # Entropies
    h_i <- -sum(p_i[p_i > 0] * log(p_i[p_i > 0]))
    h_j <- -sum(p_j[p_j > 0] * log(p_j[p_j > 0]))

    nmi <- 2 * mi / (h_i + h_j)
    return(nmi)
}

#' Compute Consensus Clustering
compute_consensus_clustering <- function(cluster_list) {
    # Create co-clustering matrix
    samples <- names(cluster_list[[1]])
    n <- length(samples)

    co_cluster <- matrix(0, n, n, dimnames = list(samples, samples))

    for (clusters in cluster_list) {
        for (i in 1:(n - 1)) {
            for (j in (i + 1):n) {
                if (clusters[i] == clusters[j]) {
                    co_cluster[i, j] <- co_cluster[i, j] + 1
                    co_cluster[j, i] <- co_cluster[j, i] + 1
                }
            }
        }
    }

    co_cluster <- co_cluster / length(cluster_list)
    diag(co_cluster) <- 1

    # Cluster the co-clustering matrix
    hc <- hclust(as.dist(1 - co_cluster), method = "ward.D2")
    k <- choose_optimal_k(co_cluster)
    consensus <- cutree(hc, k = k)

    return(consensus)
}

# -----------------------------------------------------------------------------
# Feature Importance Comparison
# -----------------------------------------------------------------------------

#' Compare Feature Importance Across Methods
#'
#' @param integration_results Integration results list
#' @param output_dir Output directory
#' @return List with feature importance comparison
compare_feature_importance <- function(integration_results, output_dir) {
    message("Comparing feature importance across methods...")

    # Get feature importance from each method
    feature_importance <- list()

    # MOFA feature weights
    if (!is.null(integration_results$mofa)) {
        feature_importance$mofa <- get_mofa_feature_importance(integration_results$mofa)
    }

    # DIABLO loadings
    if (!is.null(integration_results$diablo)) {
        feature_importance$diablo <- get_diablo_feature_importance(integration_results$diablo)
    }

    # SNF doesn't directly provide feature importance
    # Can derive from contribution to similarity

    if (length(feature_importance) < 2) {
        message("  Not enough methods with feature importance results.")
        return(NULL)
    }

    # Find common features
    all_features <- unique(unlist(lapply(feature_importance, function(x) x$feature)))

    # Create importance matrix
    importance_matrix <- matrix(
        NA,
        nrow = length(all_features),
        ncol = length(feature_importance),
        dimnames = list(all_features, names(feature_importance))
    )

    for (method in names(feature_importance)) {
        fi <- feature_importance[[method]]
        importance_matrix[match(fi$feature, all_features), method] <- fi$importance
    }

    # Compute rank correlations between methods
    methods <- colnames(importance_matrix)
    n_methods <- length(methods)

    rank_cor_matrix <- matrix(NA, n_methods, n_methods, dimnames = list(methods, methods))

    for (i in 1:(n_methods - 1)) {
        for (j in (i + 1):n_methods) {
            valid <- !is.na(importance_matrix[, i]) & !is.na(importance_matrix[, j])
            if (sum(valid) >= 10) {
                cor_test <- cor.test(
                    importance_matrix[valid, i],
                    importance_matrix[valid, j],
                    method = "spearman"
                )
                rank_cor_matrix[i, j] <- cor_test$estimate
                rank_cor_matrix[j, i] <- rank_cor_matrix[i, j]
            }
        }
    }
    diag(rank_cor_matrix) <- 1

    # Compute combined importance score
    # Rank-based aggregation
    rank_matrix <- apply(importance_matrix, 2, function(x) {
        r <- rank(-x, na.last = "keep", ties.method = "average")
        r / max(r, na.rm = TRUE)
    })

    combined_importance <- rowMeans(rank_matrix, na.rm = TRUE)

    importance_df <- data.frame(
        feature = all_features,
        importance_matrix,
        combined_rank = combined_importance,
        n_methods = rowSums(!is.na(importance_matrix)),
        stringsAsFactors = FALSE
    )

    importance_df <- importance_df[order(importance_df$combined_rank), ]

    # Save results
    write.csv(importance_df, file.path(output_dir, "feature_importance_comparison.csv"), row.names = FALSE)
    write.csv(as.data.frame(rank_cor_matrix), file.path(output_dir, "feature_rank_correlation.csv"))

    return(list(
        per_method = feature_importance,
        importance_matrix = importance_matrix,
        rank_correlation = rank_cor_matrix,
        combined = importance_df
    ))
}

#' Resolve generic feature IDs in consensus results to display names
#'
#' Replaces GENE_N / feature_N with actual gene symbols or metabolite names
#' in the importance matrix, combined table, and re-saves the CSV.
resolve_consensus_feature_names <- function(feature_comparison, feature_name_map, output_dir) {
    # Map feature names
    old_features <- rownames(feature_comparison$importance_matrix)
    new_features <- ifelse(
        old_features %in% names(feature_name_map) &
            !is.na(feature_name_map[old_features]) &
            feature_name_map[old_features] != "",
        unname(feature_name_map[old_features]),
        old_features
    )
    new_features <- make.unique(new_features, sep = "_")

    # Update importance matrix rownames
    rownames(feature_comparison$importance_matrix) <- new_features

    # Update combined table
    if (!is.null(feature_comparison$combined)) {
        idx <- match(feature_comparison$combined$feature, old_features)
        feature_comparison$combined$feature <- ifelse(
            !is.na(idx), new_features[idx], feature_comparison$combined$feature
        )
        # Re-save CSV with resolved names
        write.csv(feature_comparison$combined,
                  file.path(output_dir, "feature_importance_comparison.csv"),
                  row.names = FALSE)
    }

    # Update per-method feature names
    for (m in names(feature_comparison$per_method)) {
        fi <- feature_comparison$per_method[[m]]
        if (!is.null(fi) && "feature" %in% colnames(fi)) {
            mapped <- feature_name_map[fi$feature]
            fi$feature <- ifelse(!is.na(mapped) & mapped != "",
                                  unname(mapped), fi$feature)
            feature_comparison$per_method[[m]] <- fi
        }
    }

    feature_comparison
}

#' Get MOFA Feature Importance
get_mofa_feature_importance <- function(mofa_results) {
    tryCatch({
        # Use pre-extracted weights (survives targets serialization)
        weights <- mofa_results$weights
        if (is.null(weights) && !is.null(mofa_results$model)) {
            weights <- MOFA2::get_weights(mofa_results$model)
        }
        if (is.null(weights)) return(NULL)

        # Combine across views and factors
        all_weights <- lapply(names(weights), function(view) {
            w <- weights[[view]]
            # Aggregate across factors (max absolute weight)
            importance <- apply(abs(w), 1, max)

            data.frame(
                feature = names(importance),
                view = view,
                importance = importance,
                stringsAsFactors = FALSE
            )
        })

        combined <- do.call(rbind, all_weights)

        # Take max across views for features in multiple views
        aggregate(importance ~ feature, data = combined, FUN = max)

    }, error = function(e) {
        message("  Could not extract MOFA feature importance: ", e$message)
        NULL
    })
}

#' Get DIABLO Feature Importance
get_diablo_feature_importance <- function(diablo_results) {
    tryCatch({
        # Use pre-extracted loadings (survives targets serialization)
        loadings <- diablo_results$loadings %||% diablo_results$model$loadings
        if (is.null(loadings)) return(NULL)

        # Remove Y component (response variable, not an omics view)
        loadings <- loadings[setdiff(names(loadings), "Y")]
        if (length(loadings) == 0) return(NULL)

        # Combine across views and components
        all_loadings <- lapply(names(loadings), function(view) {
            l <- loadings[[view]]
            # Aggregate across components
            importance <- apply(abs(l), 1, max)

            data.frame(
                feature = names(importance),
                view = view,
                importance = importance,
                stringsAsFactors = FALSE
            )
        })

        combined <- do.call(rbind, all_loadings)

        aggregate(importance ~ feature, data = combined, FUN = max)

    }, error = function(e) {
        message("  Could not extract DIABLO feature importance: ", e$message)
        NULL
    })
}

# -----------------------------------------------------------------------------
# Robust Pattern Identification
# -----------------------------------------------------------------------------

#' Identify Robust Patterns
#'
#' @param sample_comparison Sample comparison results
#' @param feature_comparison Feature comparison results
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return List with robust patterns
identify_robust_patterns <- function(sample_comparison, feature_comparison, config, output_dir) {
    message("Identifying robust patterns across methods...")

    consensus_config <- config$modes$multiomics$consensus %||% list()
    robust_threshold <- consensus_config$robust_threshold %||% 0.7

    # Robust features: high importance in multiple methods
    robust_features <- NULL
    if (!is.null(feature_comparison$combined)) {
        n_methods <- ncol(feature_comparison$importance_matrix)
        min_methods <- ceiling(n_methods * robust_threshold)

        fc <- feature_comparison$combined
        robust_features <- fc[fc$n_methods >= min_methods & fc$combined_rank <= 0.2, ]
        robust_features <- robust_features[order(robust_features$combined_rank), ]

        message("  Robust features (top 20% in >=",
                round(robust_threshold * 100), "% of methods): ", nrow(robust_features))
    }

    # Robust sample groupings: consensus clusters with high agreement
    robust_clusters <- NULL
    if (!is.null(sample_comparison$ari_matrix)) {
        mean_ari <- mean(sample_comparison$ari_matrix[upper.tri(sample_comparison$ari_matrix)], na.rm = TRUE)
        mean_nmi <- mean(sample_comparison$nmi_matrix[upper.tri(sample_comparison$nmi_matrix)], na.rm = TRUE)

        message("  Mean ARI between methods: ", round(mean_ari, 3))
        message("  Mean NMI between methods: ", round(mean_nmi, 3))

        robust_clusters <- list(
            mean_ari = mean_ari,
            mean_nmi = mean_nmi,
            consensus = sample_comparison$consensus,
            agreement = mean_ari >= robust_threshold
        )
    }

    # Save results
    if (!is.null(robust_features)) {
        write.csv(robust_features, file.path(output_dir, "robust_features.csv"), row.names = FALSE)
    }

    return(list(
        features = robust_features,
        clusters = robust_clusters
    ))
}

# -----------------------------------------------------------------------------
# Method-Specific Pattern Identification
# -----------------------------------------------------------------------------

#' Identify Method-Specific Patterns
#'
#' @param feature_comparison Feature comparison results
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return List with method-specific patterns
identify_method_specific_patterns <- function(feature_comparison, config, output_dir) {
    message("Identifying method-specific patterns...")

    if (is.null(feature_comparison$importance_matrix)) {
        return(NULL)
    }

    methods <- colnames(feature_comparison$importance_matrix)
    specific_features <- list()

    for (method in methods) {
        # Features unique to this method (high importance only in this method)
        method_importance <- feature_comparison$importance_matrix[, method]
        other_importance <- feature_comparison$importance_matrix[, setdiff(methods, method), drop = FALSE]

        # Threshold: top 10% in this method
        method_threshold <- quantile(method_importance, 0.9, na.rm = TRUE)

        # Max importance in other methods
        max_other <- apply(other_importance, 1, max, na.rm = TRUE)
        other_threshold <- quantile(max_other, 0.5, na.rm = TRUE)

        specific_idx <- method_importance >= method_threshold &
            (is.na(max_other) | max_other <= other_threshold)

        specific_features[[method]] <- data.frame(
            feature = rownames(feature_comparison$importance_matrix)[specific_idx],
            importance = method_importance[specific_idx],
            stringsAsFactors = FALSE
        )

        message("  ", method, "-specific features: ", sum(specific_idx, na.rm = TRUE))
    }

    # Save results
    for (method in names(specific_features)) {
        if (nrow(specific_features[[method]]) > 0) {
            write.csv(
                specific_features[[method]],
                file.path(output_dir, paste0(method, "_specific_features.csv")),
                row.names = FALSE
            )
        }
    }

    return(specific_features)
}

# -----------------------------------------------------------------------------
# Meta-Integration
# -----------------------------------------------------------------------------

#' Run Meta-Integration Analysis
#'
#' @param integration_results Integration results list
#' @param mae_data MAE data object
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return Meta-integration results
run_meta_integration <- function(integration_results, mae_data, config, output_dir) {
    message("Running meta-integration...")

    # Combine latent factors from all methods
    latent_factors <- list()

    # MOFA factors (use pre-extracted data)
    if (!is.null(integration_results$mofa)) {
        tryCatch({
            fdf <- integration_results$mofa$factors
            if (!is.null(fdf)) {
                if (is.data.frame(fdf) && "sample_id" %in% colnames(fdf)) {
                    rn <- fdf$sample_id
                    factor_cols <- grep("^Factor|^V\\d", colnames(fdf), value = TRUE)
                    if (length(factor_cols) == 0) factor_cols <- setdiff(colnames(fdf), c("sample_id", "condition", "group"))
                    factors <- as.matrix(fdf[, factor_cols, drop = FALSE])
                    rownames(factors) <- rn
                } else if (is.matrix(fdf)) {
                    factors <- fdf
                } else {
                    factors <- NULL
                }
                if (!is.null(factors)) latent_factors$mofa <- factors
            }
        }, error = function(e) NULL)
    }

    # DIABLO variates (use pre-extracted data)
    if (!is.null(integration_results$diablo)) {
        tryCatch({
            variates <- integration_results$diablo$sample_scores %||%
                         integration_results$diablo$model$variates
            if (!is.null(variates)) {
                variates <- variates[setdiff(names(variates), "Y")]
                combined <- do.call(cbind, variates)
                latent_factors$diablo <- combined
            }
        }, error = function(e) NULL)
    }

    if (length(latent_factors) < 2) {
        message("  Not enough latent factors for meta-integration.")
        return(NULL)
    }

    # Align samples
    common_samples <- Reduce(intersect, lapply(latent_factors, rownames))

    if (length(common_samples) < 5) {
        message("  Too few common samples for meta-integration.")
        return(NULL)
    }

    aligned_factors <- lapply(latent_factors, function(x) x[common_samples, , drop = FALSE])

    # Concatenate all factors
    meta_matrix <- do.call(cbind, aligned_factors)

    # Perform meta-level dimensionality reduction (PCA on combined factors)
    meta_pca <- prcomp(meta_matrix, scale. = TRUE, center = TRUE)

    # Extract meta-factors
    meta_factors <- meta_pca$x[, 1:min(5, ncol(meta_pca$x)), drop = FALSE]

    # Cluster on meta-factors — prefer config n_clusters if specified
    config_k <- config$modes$multiomics$integration$snf$n_clusters
    if (!is.null(config_k) && is.numeric(config_k) && config_k >= 2) {
        best_k <- as.integer(config_k)
        message("  Using config n_clusters = ", best_k, " for meta-integration")
    } else {
        best_k <- choose_optimal_k(meta_factors)
        message("  Optimal k = ", best_k, " (silhouette method)")
    }
    meta_clusters <- kmeans(meta_factors, centers = best_k, nstart = 25)$cluster

    # Association with condition
    condition_assoc <- NULL
    if (!is.null(mae_data$metadata)) {
        condition_col <- mae_data$config$design$condition_column %||% "condition"
        if (condition_col %in% colnames(mae_data$metadata)) {
            condition <- mae_data$metadata[[condition_col]][match(common_samples, rownames(mae_data$metadata))]

            # Test each meta-factor for association
            condition_assoc <- sapply(1:ncol(meta_factors), function(i) {
                anova_test <- aov(meta_factors[, i] ~ factor(condition))
                summary(anova_test)[[1]]["Pr(>F)"][1, 1]
            })
        }
    }

    # Save results
    meta_df <- data.frame(
        sample_id = common_samples,
        meta_factors,
        meta_cluster = meta_clusters,
        stringsAsFactors = FALSE
    )
    write.csv(meta_df, file.path(output_dir, "meta_integration_results.csv"), row.names = FALSE)

    return(list(
        factors = meta_factors,
        clusters = meta_clusters,
        pca = meta_pca,
        condition_association = condition_assoc
    ))
}

# -----------------------------------------------------------------------------
# Visualization Functions
# -----------------------------------------------------------------------------

#' Generate Consensus Analysis Plots
generate_consensus_plots <- function(sample_comparison, feature_comparison, robust_patterns,
                                      specific_patterns, meta_integration, config, plots_dir) {
    message("Generating consensus analysis plots...")

    figures <- list()

    # 1. Method agreement heatmap (ARI/NMI)
    if (!is.null(sample_comparison$ari_matrix)) {
        figures$method_agreement <- plot_method_agreement(
            sample_comparison = sample_comparison,
            plots_dir = plots_dir
        )
    }

    # 2. Feature importance correlation
    if (!is.null(feature_comparison$rank_correlation)) {
        figures$feature_correlation <- plot_feature_correlation(
            feature_comparison = feature_comparison,
            plots_dir = plots_dir
        )
    }

    # 3. Robust features
    if (!is.null(robust_patterns$features) && nrow(robust_patterns$features) > 0) {
        figures$robust_features <- plot_robust_features(
            robust_features = robust_patterns$features,
            feature_comparison = feature_comparison,
            plots_dir = plots_dir
        )
    }

    # 4. Meta-integration
    if (!is.null(meta_integration)) {
        figures$meta_pca <- plot_meta_integration(
            meta_integration = meta_integration,
            plots_dir = plots_dir
        )
    }

    return(figures)
}

#' Plot Method Agreement Heatmap
plot_method_agreement <- function(sample_comparison, plots_dir) {
    message("  Creating method agreement heatmap...")

    # Combine ARI and NMI into one plot
    ari_df <- reshape2::melt(sample_comparison$ari_matrix)
    ari_df$metric <- "ARI"

    nmi_df <- reshape2::melt(sample_comparison$nmi_matrix)
    nmi_df$metric <- "NMI"

    combined_df <- rbind(ari_df, nmi_df)
    colnames(combined_df)[1:2] <- c("Method1", "Method2")

    p <- ggplot2::ggplot(combined_df, ggplot2::aes(x = Method1, y = Method2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label = round(value, 2)), color = "black", size = 3) +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
        ggplot2::facet_wrap(~metric) +
        ggplot2::labs(
            title = "Sample Clustering Agreement Between Methods",
            fill = "Score"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    fig_path <- file.path(plots_dir, "method_agreement_heatmap.png")
    ggplot2::ggsave(fig_path, p, width = 12, height = 5)

    return(fig_path)
}

#' Plot Feature Importance Correlation
plot_feature_correlation <- function(feature_comparison, plots_dir) {
    message("  Creating feature correlation heatmap...")

    cor_df <- reshape2::melt(feature_comparison$rank_correlation)
    colnames(cor_df) <- c("Method1", "Method2", "Correlation")

    p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Method1, y = Method2, fill = Correlation)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 2)), color = "black", size = 4) +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        ggplot2::labs(
            title = "Feature Importance Rank Correlation Between Methods",
            subtitle = "Spearman correlation of feature rankings",
            fill = "Correlation"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    fig_path <- file.path(plots_dir, "feature_rank_correlation.png")
    ggplot2::ggsave(fig_path, p, width = 8, height = 6)

    return(fig_path)
}

#' Plot Robust Features
plot_robust_features <- function(robust_features, feature_comparison, plots_dir) {
    message("  Creating robust features plot...")

    # Top 30 robust features
    top_features <- head(robust_features, 30)

    # Get importance across methods
    methods <- colnames(feature_comparison$importance_matrix)

    # Rank features within each method (1 = most important) for comparable scale
    rank_matrix <- apply(feature_comparison$importance_matrix, 2, function(x) {
        r <- rep(NA_real_, length(x))
        valid <- !is.na(x)
        r[valid] <- rank(-x[valid], ties.method = "average")
        r
    })
    rownames(rank_matrix) <- rownames(feature_comparison$importance_matrix)
    n_total <- colSums(!is.na(feature_comparison$importance_matrix))

    plot_data <- lapply(methods, function(m) {
        idx <- match(top_features$feature, rownames(rank_matrix))
        # Express rank as percentile (0% = top, 100% = bottom)
        pct <- rank_matrix[idx, m] / n_total[m] * 100
        data.frame(
            feature = top_features$feature,
            method = m,
            rank_pct = pct,
            stringsAsFactors = FALSE
        )
    })

    plot_df <- do.call(rbind, plot_data)
    plot_df$feature <- factor(plot_df$feature, levels = rev(top_features$feature))

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = method, y = feature, fill = rank_pct)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c(option = "plasma", na.value = "grey90",
                                       direction = -1,
                                       labels = function(x) paste0(round(x), "%")) +
        ggplot2::labs(
            title = "Top Robust Features Across Methods",
            subtitle = "Features consistently important in multiple methods",
            x = "Method",
            y = "Feature",
            fill = "Rank\n(percentile)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

    fig_path <- file.path(plots_dir, "robust_features.png")
    ggplot2::ggsave(fig_path, p, width = 10, height = 12)

    return(fig_path)
}

#' Plot Meta-Integration Results
plot_meta_integration <- function(meta_integration, plots_dir) {
    message("  Creating meta-integration plot...")

    factors <- as.data.frame(meta_integration$factors)
    factors$cluster <- factor(meta_integration$clusters)
    factors$sample_id <- rownames(meta_integration$factors)

    # Variance explained
    var_exp <- summary(meta_integration$pca)$importance[2, 1:2] * 100

    base_aes <- ggplot2::aes(x = PC1, y = PC2, color = cluster)
    base_labs <- ggplot2::labs(
        title = "Meta-Integration: Combined Latent Space",
        subtitle = "PCA of concatenated latent factors from all methods",
        x = paste0("Meta-PC1 (", round(var_exp[1], 1), "%)"),
        y = paste0("Meta-PC2 (", round(var_exp[2], 1), "%)"),
        color = "Cluster"
    )

    # Unlabeled version
    p <- ggplot2::ggplot(factors, base_aes) +
        ggplot2::geom_point(size = 3, alpha = 0.7) +
        ggplot2::stat_ellipse(level = 0.95) +
        base_labs +
        ggplot2::theme_minimal()

    fig_path <- file.path(plots_dir, "meta_integration_pca.png")
    ggplot2::ggsave(fig_path, p, width = 8, height = 6)

    # Labeled version (with sample names)
    p_labeled <- ggplot2::ggplot(factors, base_aes) +
        ggplot2::geom_point(size = 3, alpha = 0.7) +
        ggplot2::stat_ellipse(level = 0.95) +
        base_labs +
        ggplot2::theme_minimal()

    if (requireNamespace("ggrepel", quietly = TRUE)) {
        p_labeled <- p_labeled + ggrepel::geom_text_repel(
            ggplot2::aes(label = sample_id),
            size = 2.8, max.overlaps = 20, show.legend = FALSE
        )
    } else {
        p_labeled <- p_labeled + ggplot2::geom_text(
            ggplot2::aes(label = sample_id),
            size = 2.8, vjust = -0.8, show.legend = FALSE
        )
    }

    fig_path_labeled <- file.path(plots_dir, "meta_integration_pca_labeled.png")
    ggplot2::ggsave(fig_path_labeled, p_labeled, width = 8, height = 6)

    return(fig_path)
}

# -----------------------------------------------------------------------------
# Output Functions
# -----------------------------------------------------------------------------

#' Create Consensus Summary
create_consensus_summary <- function(sample_comparison, feature_comparison, robust_patterns) {
    summary <- list()

    if (!is.null(sample_comparison)) {
        summary$mean_method_ari <- mean(sample_comparison$ari_matrix[upper.tri(sample_comparison$ari_matrix)], na.rm = TRUE)
        summary$mean_method_nmi <- mean(sample_comparison$nmi_matrix[upper.tri(sample_comparison$nmi_matrix)], na.rm = TRUE)
        summary$n_consensus_clusters <- length(unique(sample_comparison$consensus))
    }

    if (!is.null(feature_comparison)) {
        summary$mean_feature_correlation <- mean(feature_comparison$rank_correlation[upper.tri(feature_comparison$rank_correlation)], na.rm = TRUE)
    }

    if (!is.null(robust_patterns$features)) {
        summary$n_robust_features <- nrow(robust_patterns$features)
        summary$top_robust_features <- head(robust_patterns$features$feature, 10)
    }

    return(summary)
}

#' Save Consensus Outputs
save_consensus_outputs <- function(results, output_dir) {
    message("Saving consensus analysis outputs...")

    saveRDS(results, file.path(output_dir, "consensus_results.rds"))

    message("  Outputs saved to: ", output_dir)
}
