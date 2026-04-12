# =============================================================================
# Foundational Cross-Omics Correlations and Concordance Analysis
# =============================================================================
# This script runs BEFORE advanced integration methods (MOFA, DIABLO, SNF)
# to establish foundational understanding of cross-omics relationships.
#
# Analyses:
#   1. Cross-omics feature correlations (pairwise, partial)
#   2. Sample-level concordance across omics
#   3. Pathway overlap analysis
#   4. Cross-omics co-expression modules
# =============================================================================

#' Main function: Run foundational cross-omics analysis
#'
#' @param mae MultiAssayExperiment object with aligned samples
#' @param de_results Named list of DE results per omics (optional)
#' @param gene_protein_mapping Gene mapping table for cross-omics matching
#' @param config Pipeline configuration
#' @param out_dir Output directory for results and plots
#' @return List containing all foundational analysis results
run_foundational_analysis <- function(mae, de_results = NULL,
                                       gene_protein_mapping = NULL,
                                       config, out_dir = NULL) {
    message("=== Running Foundational Cross-Omics Analysis ===")
    message("This runs BEFORE integration methods to understand basic relationships")

    # Convert MAE to legacy format
    mae_data <- .mae_to_legacy_foundational(mae, de_results, gene_protein_mapping)
    harmonized <- mae_data$harmonized_omics
    omics_names <- names(harmonized)

    # Create output directories if specified
    if (!is.null(out_dir)) {
        dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
    }

    # Check we have at least 2 omics
    if (length(omics_names) < 2) {
        message("Need at least 2 omics for cross-omics analysis. Skipping.")
        return(NULL)
    }

    # Get foundational config with defaults
    fc <- get_foundational_config(config)

    results <- list()

    # 1. Cross-omics feature correlations
    results$feature_correlations <- tryCatch(
        compute_crossomics_feature_correlations(harmonized, mae_data$gene_mapping, fc, out_dir),
        error = function(e) {
            message("WARNING: Feature correlation analysis failed: ", conditionMessage(e))
            NULL
        }
    )

    # 2. Sample-level concordance
    results$sample_concordance <- tryCatch(
        compute_sample_concordance(harmonized, mae_data$metadata, fc, config, out_dir),
        error = function(e) {
            message("WARNING: Sample concordance analysis failed: ", conditionMessage(e))
            NULL
        }
    )

    # 3. Pathway overlap analysis
    results$pathway_overlap <- tryCatch(
        analyze_pathway_overlap(harmonized, mae_data$gene_mapping, fc, out_dir),
        error = function(e) {
            message("WARNING: Pathway overlap analysis failed: ", conditionMessage(e))
            NULL
        }
    )

    # 4. Cross-omics module detection (if enabled)
    if (fc$module_detection) {
        results$crossomics_modules <- tryCatch(
            find_crossomics_modules(harmonized, mae_data$gene_mapping, fc, out_dir),
            error = function(e) {
                message("WARNING: Module detection failed: ", conditionMessage(e))
                NULL
            }
        )
    }

    # Generate summary
    results$summary <- tryCatch(
        summarize_foundational_results(results, out_dir),
        error = function(e) {
            message("WARNING: Summary generation failed: ", conditionMessage(e))
            list(status = "partial", error = conditionMessage(e))
        }
    )

    message("=== Foundational Analysis Complete ===")
    return(results)
}

#' Get foundational analysis configuration with defaults
#' @param config Full config object
#' @return List of foundational config parameters
get_foundational_config <- function(config) {
    fc <- config$modes$multiomics$foundational %||% list()

    list(
        run_foundational = fc$run_foundational %||% TRUE,
        correlation_method = fc$correlation_method %||% "spearman",
        partial_correlations = fc$partial_correlations %||% TRUE,
        confounder_columns = fc$confounder_columns %||% NULL,
        min_correlation = fc$min_correlation %||% 0.3,
        fdr_threshold = fc$fdr_threshold %||% 0.05,
        module_detection = fc$module_detection %||% TRUE,
        min_common_features = fc$min_common_features %||% 50,
        top_variable_features = fc$top_variable_features %||% 2000
    )
}


# =============================================================================
# 1. Cross-Omics Feature Correlations
# =============================================================================

#' Compute cross-omics feature correlations
#' @param harmonized List of harmonized omics data
#' @param gene_mapping Gene mapping table for cross-omics matching
#' @param fc Foundational config
#' @param out_dir Output directory
compute_crossomics_feature_correlations <- function(harmonized, gene_mapping, fc, out_dir = NULL) {
    message("Computing cross-omics feature correlations...")

    omics_names <- names(harmonized)
    results <- list()

    # Generate all pairwise combinations
    pairs <- utils::combn(omics_names, 2, simplify = FALSE)

    for (pair in pairs) {
        omics1 <- pair[1]
        omics2 <- pair[2]
        pair_name <- paste(omics1, omics2, sep = "_vs_")

        message("Processing: ", omics1, " vs ", omics2)

        # Get matrices and find common features
        mat1 <- harmonized[[omics1]]$normalized_matrix
        mat2 <- harmonized[[omics2]]$normalized_matrix

        # Get common samples
        common_samples <- intersect(colnames(mat1), colnames(mat2))
        if (length(common_samples) < 5) {
            message("  Insufficient common samples (", length(common_samples), "). Skipping.")
            next
        }

        mat1 <- mat1[, common_samples, drop = FALSE]
        mat2 <- mat2[, common_samples, drop = FALSE]

        # Map features to common identifiers (gene symbols)
        mapping_result <- map_features_to_common_ids(
            mat1, mat2, omics1, omics2, gene_mapping
        )

        if (is.null(mapping_result) || mapping_result$n_common < fc$min_common_features) {
            message("  Insufficient common features. Skipping.")
            next
        }

        # Compute correlations
        pair_result <- compute_pairwise_correlations(
            mapping_result$mat1_matched,
            mapping_result$mat2_matched,
            mapping_result$common_ids,
            fc
        )

        pair_result$omics1 <- omics1
        pair_result$omics2 <- omics2
        pair_result$n_samples <- length(common_samples)

        # Compute partial correlations if confounders specified
        if (fc$partial_correlations && !is.null(fc$confounder_columns)) {
            metadata <- harmonized[[omics1]]$metadata %||% NULL
            if (!is.null(metadata)) {
                pair_result$partial <- compute_partial_correlations(
                    mapping_result$mat1_matched,
                    mapping_result$mat2_matched,
                    mapping_result$common_ids,
                    metadata[common_samples, , drop = FALSE],
                    fc
                )
            }
        }

        # Build correlation network
        pair_result$network <- build_correlation_network(pair_result, fc)

        # Save results
        if (!is.null(out_dir)) {
            save_pairwise_correlation_results(pair_result, pair_name, out_dir)
        }

        results[[pair_name]] <- pair_result
    }

    # Check if we have any results
    pair_names <- setdiff(names(results), "summary")
    if (length(pair_names) == 0) {
        message("No cross-omics correlations could be computed.")
        return(list(
            summary = data.frame(
                pair = character(),
                omics1 = character(),
                omics2 = character(),
                n_features = integer(),
                n_samples = integer(),
                mean_correlation = numeric(),
                median_correlation = numeric(),
                n_significant = integer(),
                pct_positive = numeric(),
                stringsAsFactors = FALSE
            )
        ))
    }

    # Create summary
    results$summary <- summarize_correlation_results(results, out_dir)

    # Generate visualizations
    if (!is.null(out_dir)) {
        plot_correlation_overview(results, out_dir)
    }

    return(results)
}

#' Map features from two omics to common identifiers
map_features_to_common_ids <- function(mat1, mat2, omics1, omics2, gene_mapping) {
    # Always try direct matching by rownames first — works when MAE was
    # pre-harmonized (features renamed to GENE_1, GENE_2, ...)
    common_ids <- intersect(rownames(mat1), rownames(mat2))
    if (length(common_ids) > 0) {
        message("  Direct row-name matching: ", length(common_ids), " common features")
        return(list(
            mat1_matched = mat1[common_ids, , drop = FALSE],
            mat2_matched = mat2[common_ids, , drop = FALSE],
            common_ids = common_ids,
            n_common = length(common_ids)
        ))
    }

    # Fall back to gene mapping if no direct matches
    if (is.null(gene_mapping)) return(NULL)

    # Use gene mapping
    map1 <- gene_mapping[gene_mapping$omics == omics1, ]
    map2 <- gene_mapping[gene_mapping$omics == omics2, ]

    # Find common gene symbols
    common_symbols <- intersect(map1$gene_symbol, map2$gene_symbol)
    common_symbols <- common_symbols[!is.na(common_symbols) & common_symbols != ""]

    if (length(common_symbols) == 0) return(NULL)

    # Get feature IDs for each omics
    features1 <- map1$feature_id[match(common_symbols, map1$gene_symbol)]
    features2 <- map2$feature_id[match(common_symbols, map2$gene_symbol)]

    # Filter to features present in matrices
    valid_idx <- features1 %in% rownames(mat1) & features2 %in% rownames(mat2)

    if (sum(valid_idx) == 0) return(NULL)

    common_symbols <- common_symbols[valid_idx]
    features1 <- features1[valid_idx]
    features2 <- features2[valid_idx]

    list(
        mat1_matched = mat1[features1, , drop = FALSE],
        mat2_matched = mat2[features2, , drop = FALSE],
        common_ids = common_symbols,
        n_common = length(common_symbols)
    )
}

#' Compute pairwise correlations between matched features
compute_pairwise_correlations <- function(mat1, mat2, feature_ids, fc) {
    n_features <- length(feature_ids)

    # Initialize results
    correlations <- numeric(n_features)
    pvalues <- numeric(n_features)
    names(correlations) <- feature_ids
    names(pvalues) <- feature_ids

    # Compute correlations
    for (i in seq_len(n_features)) {
        x <- mat1[i, ]
        y <- mat2[i, ]

        # Skip if no variance
        if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
            correlations[i] <- NA
            pvalues[i] <- NA
            next
        }

        # Correlation test
        ct <- tryCatch({
            cor.test(x, y, method = fc$correlation_method, use = "pairwise.complete.obs")
        }, error = function(e) NULL)

        if (!is.null(ct)) {
            correlations[i] <- ct$estimate
            pvalues[i] <- ct$p.value
        } else {
            correlations[i] <- NA
            pvalues[i] <- NA
        }
    }

    # Adjust p-values
    padj <- p.adjust(pvalues, method = "BH")

    # Create results data frame
    cor_df <- data.frame(
        feature_id = feature_ids,
        correlation = correlations,
        pvalue = pvalues,
        padj = padj,
        significant = padj < fc$fdr_threshold & abs(correlations) >= fc$min_correlation,
        stringsAsFactors = FALSE
    )

    # Summary statistics
    valid_cors <- correlations[!is.na(correlations)]
    summary_stats <- list(
        n_features = n_features,
        n_valid = length(valid_cors),
        mean_cor = mean(valid_cors),
        median_cor = median(valid_cors),
        sd_cor = sd(valid_cors),
        n_positive = sum(valid_cors > 0),
        n_negative = sum(valid_cors < 0),
        n_significant = sum(cor_df$significant, na.rm = TRUE),
        pct_positive = 100 * mean(valid_cors > 0)
    )

    message("  Computed ", length(valid_cors), " correlations: ",
            "mean=", round(summary_stats$mean_cor, 3),
            ", median=", round(summary_stats$median_cor, 3),
            ", ", summary_stats$n_significant, " significant")

    list(
        correlations = correlations,
        pvalues = pvalues,
        padj = padj,
        cor_df = cor_df,
        summary = summary_stats
    )
}

#' Compute partial correlations controlling for confounders
compute_partial_correlations <- function(mat1, mat2, feature_ids, metadata, fc) {
    message("  Computing partial correlations controlling for: ",
            paste(fc$confounder_columns, collapse = ", "))

    # Check if ppcor is available
    if (!requireNamespace("ppcor", quietly = TRUE)) {
        message("  ppcor package not available. Skipping partial correlations.")
        return(NULL)
    }

    # Get confounder matrix
    confounder_cols <- intersect(fc$confounder_columns, colnames(metadata))
    if (length(confounder_cols) == 0) {
        message("  No confounder columns found in metadata. Skipping.")
        return(NULL)
    }

    # Prepare confounder matrix (numeric encoding)
    Z <- as.data.frame(metadata[, confounder_cols, drop = FALSE])
    for (col in colnames(Z)) {
        if (is.character(Z[[col]]) || is.factor(Z[[col]])) {
            Z[[col]] <- as.numeric(as.factor(Z[[col]]))
        }
    }
    Z <- as.matrix(Z)

    n_features <- length(feature_ids)
    partial_cors <- numeric(n_features)
    partial_pvals <- numeric(n_features)
    names(partial_cors) <- feature_ids
    names(partial_pvals) <- feature_ids

    for (i in seq_len(n_features)) {
        x <- mat1[i, ]
        y <- mat2[i, ]

        # Skip if no variance
        if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
            partial_cors[i] <- NA
            partial_pvals[i] <- NA
            next
        }

        # Partial correlation
        pt <- tryCatch({
            data_mat <- cbind(x, y, Z)
            complete_idx <- complete.cases(data_mat)
            if (sum(complete_idx) < 10) return(NULL)

            ppcor::pcor.test(x[complete_idx], y[complete_idx], Z[complete_idx, , drop = FALSE],
                             method = fc$correlation_method)
        }, error = function(e) NULL)

        if (!is.null(pt)) {
            partial_cors[i] <- pt$estimate
            partial_pvals[i] <- pt$p.value
        } else {
            partial_cors[i] <- NA
            partial_pvals[i] <- NA
        }
    }

    partial_padj <- p.adjust(partial_pvals, method = "BH")

    partial_df <- data.frame(
        feature_id = feature_ids,
        partial_correlation = partial_cors,
        partial_pvalue = partial_pvals,
        partial_padj = partial_padj,
        stringsAsFactors = FALSE
    )

    valid_pcors <- partial_cors[!is.na(partial_cors)]
    message("  Partial correlations: mean=", round(mean(valid_pcors), 3),
            ", median=", round(median(valid_pcors), 3))

    list(
        partial_correlations = partial_cors,
        partial_df = partial_df,
        confounders = confounder_cols
    )
}

#' Build correlation network from significant correlations
build_correlation_network <- function(cor_results, fc) {
    cor_df <- cor_results$cor_df

    # Filter to significant correlations
    sig_df <- cor_df[cor_df$significant == TRUE & !is.na(cor_df$significant), ]

    if (nrow(sig_df) == 0) {
        message("  No significant correlations for network. Skipping.")
        return(NULL)
    }

    message("  Building network from ", nrow(sig_df), " significant correlations")

    # Check if igraph is available
    if (!requireNamespace("igraph", quietly = TRUE)) {
        message("  igraph package not available. Returning edge list only.")
        return(list(edges = sig_df, graph = NULL))
    }

    # Create edge list
    edges <- data.frame(
        from = paste0("omics1_", sig_df$feature_id),
        to = paste0("omics2_", sig_df$feature_id),
        weight = abs(sig_df$correlation),
        correlation = sig_df$correlation,
        stringsAsFactors = FALSE
    )

    g <- igraph::graph_from_data_frame(edges, directed = FALSE)

    # Network metrics
    metrics <- list(
        n_nodes = igraph::vcount(g),
        n_edges = igraph::ecount(g),
        density = igraph::edge_density(g),
        mean_weight = mean(edges$weight)
    )

    list(
        edges = sig_df,
        graph = g,
        metrics = metrics
    )
}

#' Save pairwise correlation results
save_pairwise_correlation_results <- function(pair_result, pair_name, out_dir) {
    # Save correlation table
    write.csv(pair_result$cor_df,
              file.path(out_dir, "tables", paste0("crossomics_correlations_", pair_name, ".csv")),
              row.names = FALSE)

    # Save partial correlations if available
    if (!is.null(pair_result$partial)) {
        write.csv(pair_result$partial$partial_df,
                  file.path(out_dir, "tables", paste0("partial_correlations_", pair_name, ".csv")),
                  row.names = FALSE)
    }

    # Save network edges if available
    if (!is.null(pair_result$network) && !is.null(pair_result$network$edges)) {
        write.csv(pair_result$network$edges,
                  file.path(out_dir, "tables", paste0("correlation_network_edges_", pair_name, ".csv")),
                  row.names = FALSE)
    }
}

#' Summarize correlation results across all pairs
summarize_correlation_results <- function(results, out_dir = NULL) {
    # Filter to actual pair results (not summary)
    pair_names <- setdiff(names(results), "summary")

    summary_rows <- lapply(pair_names, function(pn) {
        pr <- results[[pn]]
        data.frame(
            pair = pn,
            omics1 = pr$omics1,
            omics2 = pr$omics2,
            n_features = pr$summary$n_features,
            n_samples = pr$n_samples,
            mean_correlation = round(pr$summary$mean_cor, 4),
            median_correlation = round(pr$summary$median_cor, 4),
            n_significant = pr$summary$n_significant,
            pct_positive = round(pr$summary$pct_positive, 1),
            stringsAsFactors = FALSE
        )
    })

    summary_df <- do.call(rbind, summary_rows)

    if (!is.null(out_dir)) {
        write.csv(summary_df,
                  file.path(out_dir, "tables", "crossomics_correlation_summary.csv"),
                  row.names = FALSE)
    }

    summary_df
}

#' Plot correlation overview
plot_correlation_overview <- function(results, out_dir) {
    pair_names <- setdiff(names(results), "summary")
    if (length(pair_names) == 0) return(NULL)

    # Combine all correlations for overview
    all_cors <- lapply(pair_names, function(pn) {
        pr <- results[[pn]]
        data.frame(
            pair = pn,
            correlation = pr$correlations[!is.na(pr$correlations)],
            stringsAsFactors = FALSE
        )
    })
    all_cors_df <- do.call(rbind, all_cors)

    # Density plot by pair
    p1 <- ggplot2::ggplot(all_cors_df, ggplot2::aes(x = correlation, fill = pair)) +
        ggplot2::geom_density(alpha = 0.5) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "Cross-Omics Correlation Distributions",
            subtitle = "Feature-level correlations between omics layers",
            x = "Correlation",
            y = "Density",
            fill = "Omics Pair"
        ) +
        ggplot2::theme(legend.position = "bottom")

    ggplot2::ggsave(file.path(out_dir, "plots", "crossomics_correlation_density.png"),
                    p1, width = 10, height = 6)

    # Box plot by pair
    p2 <- ggplot2::ggplot(all_cors_df, ggplot2::aes(x = pair, y = correlation, fill = pair)) +
        ggplot2::geom_boxplot(alpha = 0.7) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "Cross-Omics Correlation Summary",
            x = "Omics Pair",
            y = "Correlation"
        ) +
        ggplot2::theme(legend.position = "none",
                       axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    ggplot2::ggsave(file.path(out_dir, "plots", "crossomics_correlation_boxplot.png"),
                    p2, width = 8, height = 6)

    list(density_plot = p1, boxplot = p2)
}


# =============================================================================
# 2. Sample-Level Concordance Analysis
# =============================================================================

#' Compute sample-level concordance across omics
#' @param harmonized List of harmonized omics data
#' @param metadata Sample metadata
#' @param fc Foundational config
#' @param config Full config
#' @param out_dir Output directory
compute_sample_concordance <- function(harmonized, metadata, fc, config, out_dir = NULL) {
    message("Computing sample-level concordance...")

    omics_names <- names(harmonized)
    results <- list()

    # Get common samples across all omics
    all_samples <- lapply(harmonized, function(h) colnames(h$normalized_matrix))
    common_samples <- Reduce(intersect, all_samples)

    if (length(common_samples) < 5) {
        message("Insufficient common samples (", length(common_samples), "). Skipping.")
        return(NULL)
    }

    message("  Analyzing ", length(common_samples), " samples across ",
            length(omics_names), " omics")

    # 1. Per-sample rank correlation across omics
    results$sample_rank_cors <- tryCatch(
        compute_sample_rank_correlations(harmonized, common_samples, fc),
        error = function(e) {
            message("  WARNING: Sample rank correlation failed: ", conditionMessage(e))
            NULL
        }
    )

    # 2. Sample clustering consistency (NMI, ARI)
    results$clustering_consistency <- tryCatch(
        compute_clustering_consistency(harmonized, common_samples, metadata, fc, config),
        error = function(e) {
            message("  WARNING: Clustering consistency failed: ", conditionMessage(e))
            NULL
        }
    )

    # 3. Within-condition vs between-condition similarity
    if (!is.null(metadata)) {
        condition_col <- config$design$condition_column %||% "condition"
        if (condition_col %in% colnames(metadata)) {
            results$condition_similarity <- tryCatch(
                compute_condition_similarity(harmonized, common_samples, metadata, condition_col, fc),
                error = function(e) {
                    message("  WARNING: Condition similarity failed: ", conditionMessage(e))
                    NULL
                }
            )
        }
    }

    # 4. Consensus sample clustering
    results$consensus_clustering <- tryCatch(
        compute_foundational_consensus_clustering(harmonized, common_samples, fc),
        error = function(e) {
            message("  WARNING: Consensus clustering failed: ", conditionMessage(e))
            NULL
        }
    )

    # 5. Identify discordant samples
    results$discordant_samples <- tryCatch(
        identify_discordant_samples(results, common_samples, fc),
        error = function(e) {
            message("  WARNING: Discordant sample identification failed: ", conditionMessage(e))
            NULL
        }
    )

    # Save and visualize
    if (!is.null(out_dir)) {
        tryCatch(
            save_sample_concordance_results(results, out_dir),
            error = function(e) message("  WARNING: Saving concordance results failed: ", conditionMessage(e))
        )
        tryCatch(
            plot_sample_concordance(results, metadata, out_dir),
            error = function(e) message("  WARNING: Plotting concordance failed: ", conditionMessage(e))
        )
    }

    return(results)
}

#' Compute per-sample rank correlations
compute_sample_rank_correlations <- function(harmonized, common_samples, fc) {
    omics_names <- names(harmonized)

    # For each sample, correlate its expression profile across omics
    omics_matrices <- lapply(harmonized, function(h) {
        mat <- h$normalized_matrix[, common_samples, drop = FALSE]
        vars <- apply(mat, 1, var, na.rm = TRUE)
        top_idx <- order(vars, decreasing = TRUE)[1:min(fc$top_variable_features, nrow(mat))]
        mat[top_idx, , drop = FALSE]
    })

    # Pairwise sample correlations per omics
    pairs <- utils::combn(omics_names, 2, simplify = FALSE)

    n_perm <- 999
    sample_cors <- lapply(pairs, function(pair) {
        omics1 <- pair[1]
        omics2 <- pair[2]

        mat1 <- omics_matrices[[omics1]]
        mat2 <- omics_matrices[[omics2]]

        # Correlate sample-sample distance matrices
        dist1 <- as.matrix(dist(t(mat1)))
        dist2 <- as.matrix(dist(t(mat2)))

        obs_r <- cor(as.vector(dist1), as.vector(dist2), method = fc$correlation_method)

        # Permutation test for significance
        n_s <- ncol(mat1)
        perm_r <- vapply(seq_len(n_perm), function(i) {
            idx <- sample(n_s)
            cor(as.vector(dist1), as.vector(dist2[idx, idx]), method = fc$correlation_method)
        }, numeric(1))
        p_value <- (sum(perm_r >= obs_r) + 1) / (n_perm + 1)

        list(
            pair = paste(pair, collapse = "_vs_"),
            mantel_correlation = obs_r,
            mantel_pvalue = p_value
        )
    })

    # Summary
    mantel_cors <- sapply(sample_cors, function(x) x$mantel_correlation)
    mantel_pvals <- sapply(sample_cors, function(x) x$mantel_pvalue)
    names(mantel_cors) <- sapply(sample_cors, function(x) x$pair)
    names(mantel_pvals) <- names(mantel_cors)

    message("  Sample distance correlations (Mantel): ",
            paste(names(mantel_cors), "=", round(mantel_cors, 3),
                  " (p=", format(mantel_pvals, digits = 3), ")", collapse = ", "))

    list(
        mantel_correlations = mantel_cors,
        mantel_pvalues = mantel_pvals,
        details = sample_cors,
        distance_matrices = omics_matrices
    )
}

#' Compute clustering consistency across omics (NMI, ARI)
compute_clustering_consistency <- function(harmonized, common_samples, metadata, fc, config) {
    message("  Computing clustering consistency...")

    omics_names <- names(harmonized)

    # Cluster samples in each omics
    cluster_results <- lapply(omics_names, function(omics_name) {
        mat <- harmonized[[omics_name]]$normalized_matrix[, common_samples, drop = FALSE]

        # PCA
        vars <- apply(mat, 1, var, na.rm = TRUE)
        top_idx <- order(vars, decreasing = TRUE)[1:min(1000, nrow(mat))]
        mat_sub <- mat[top_idx, ]
        mat_sub <- mat_sub[complete.cases(mat_sub), ]

        if (nrow(mat_sub) < 10) return(NULL)

        pca <- prcomp(t(mat_sub), scale. = TRUE, center = TRUE)

        # K-means clustering
        k <- 3
        if (!is.null(metadata)) {
            condition_col <- config$design$condition_column %||% "condition"
            if (condition_col %in% colnames(metadata)) {
                k <- length(unique(metadata[[condition_col]]))
            }
        }
        k <- min(k, floor(length(common_samples) / 2))
        k <- max(k, 2)

        km <- kmeans(pca$x[, 1:min(5, ncol(pca$x))], centers = k, nstart = 25)

        list(
            clusters = km$cluster,
            pca = pca,
            k = k
        )
    })
    names(cluster_results) <- omics_names

    # Remove NULL results
    cluster_results <- cluster_results[!sapply(cluster_results, is.null)]

    if (length(cluster_results) < 2) {
        message("  Insufficient omics for clustering comparison")
        return(NULL)
    }

    # Compute pairwise ARI and NMI
    pairs <- utils::combn(names(cluster_results), 2, simplify = FALSE)

    ari_nmi <- lapply(pairs, function(pair) {
        c1 <- cluster_results[[pair[1]]]$clusters
        c2 <- cluster_results[[pair[2]]]$clusters

        ari <- compute_ari(c1, c2)
        nmi <- compute_nmi(c1, c2)

        list(
            pair = paste(pair, collapse = "_vs_"),
            ari = ari,
            nmi = nmi
        )
    })

    # Summary
    ari_values <- sapply(ari_nmi, function(x) x$ari)
    nmi_values <- sapply(ari_nmi, function(x) x$nmi)
    names(ari_values) <- names(nmi_values) <- sapply(ari_nmi, function(x) x$pair)

    message("  Clustering consistency (ARI): ",
            paste(names(ari_values), "=", round(ari_values, 3), collapse = ", "))

    list(
        ari = ari_values,
        nmi = nmi_values,
        cluster_results = cluster_results,
        details = ari_nmi
    )
}

#' Compute Adjusted Rand Index
compute_ari <- function(c1, c2) {
    if (requireNamespace("aricode", quietly = TRUE)) {
        return(aricode::ARI(c1, c2))
    }

    # Simple implementation
    n <- length(c1)
    tab <- table(c1, c2)

    sum_comb_a <- sum(choose(rowSums(tab), 2))
    sum_comb_b <- sum(choose(colSums(tab), 2))
    sum_comb <- sum(choose(tab, 2))

    expected <- (sum_comb_a * sum_comb_b) / choose(n, 2)
    max_index <- 0.5 * (sum_comb_a + sum_comb_b)

    if (max_index == expected) return(1)

    (sum_comb - expected) / (max_index - expected)
}

#' Compute Normalized Mutual Information
compute_nmi <- function(c1, c2) {
    if (requireNamespace("aricode", quietly = TRUE)) {
        return(aricode::NMI(c1, c2))
    }

    # Simple implementation
    tab <- table(c1, c2)
    n <- sum(tab)

    # Joint entropy
    p_joint <- tab / n
    p_joint <- p_joint[p_joint > 0]
    h_joint <- -sum(p_joint * log(p_joint))

    # Marginal entropies
    p1 <- table(c1) / n
    p2 <- table(c2) / n
    h1 <- -sum(p1 * log(p1))
    h2 <- -sum(p2 * log(p2))

    # Mutual information
    mi <- h1 + h2 - h_joint

    # Normalized
    if (h1 + h2 == 0) return(0)
    2 * mi / (h1 + h2)
}

#' Compute within-condition vs between-condition similarity
compute_condition_similarity <- function(harmonized, common_samples, metadata,
                                          condition_col, fc) {
    message("  Computing within vs between-condition similarity...")

    omics_names <- names(harmonized)
    conditions <- metadata[[condition_col]][match(common_samples, rownames(metadata))]

    results <- lapply(omics_names, function(omics_name) {
        mat <- harmonized[[omics_name]]$normalized_matrix[, common_samples, drop = FALSE]

        # Compute sample correlation matrix
        cor_mat <- cor(mat, use = "pairwise.complete.obs")

        # Separate within and between condition correlations
        within_cors <- c()
        between_cors <- c()

        for (i in 1:(length(common_samples) - 1)) {
            for (j in (i + 1):length(common_samples)) {
                r <- cor_mat[i, j]
                if (is.na(r)) next

                if (conditions[i] == conditions[j]) {
                    within_cors <- c(within_cors, r)
                } else {
                    between_cors <- c(between_cors, r)
                }
            }
        }

        list(
            omics = omics_name,
            within_mean = mean(within_cors, na.rm = TRUE),
            within_sd = sd(within_cors, na.rm = TRUE),
            between_mean = mean(between_cors, na.rm = TRUE),
            between_sd = sd(between_cors, na.rm = TRUE),
            separation = mean(within_cors, na.rm = TRUE) - mean(between_cors, na.rm = TRUE)
        )
    })
    names(results) <- omics_names

    # Summary data frame
    summary_df <- do.call(rbind, lapply(results, function(r) {
        data.frame(
            omics = r$omics,
            within_condition_cor = round(r$within_mean, 4),
            between_condition_cor = round(r$between_mean, 4),
            separation = round(r$separation, 4),
            stringsAsFactors = FALSE
        )
    }))

    message("  Condition separation: ",
            paste(summary_df$omics, "=", round(summary_df$separation, 3), collapse = ", "))

    list(
        per_omics = results,
        summary = summary_df
    )
}

#' Compute consensus clustering across omics
compute_foundational_consensus_clustering <- function(harmonized, common_samples, fc) {
    message("  Computing consensus clustering...")

    # Check if ConsensusClusterPlus is available
    if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
        message("  ConsensusClusterPlus not available. Skipping.")
        return(NULL)
    }

    omics_names <- names(harmonized)

    # Get PCA embeddings for each omics
    embeddings <- lapply(omics_names, function(omics_name) {
        mat <- harmonized[[omics_name]]$normalized_matrix[, common_samples, drop = FALSE]
        vars <- apply(mat, 1, var, na.rm = TRUE)
        top_idx <- order(vars, decreasing = TRUE)[1:min(1000, nrow(mat))]
        mat_sub <- mat[top_idx, ]
        mat_sub <- mat_sub[complete.cases(mat_sub), ]

        if (nrow(mat_sub) < 10) return(NULL)

        pca <- prcomp(t(mat_sub), scale. = TRUE, center = TRUE)
        pca$x[, 1:min(10, ncol(pca$x))]
    })
    names(embeddings) <- omics_names
    embeddings <- embeddings[!sapply(embeddings, is.null)]

    if (length(embeddings) < 2) return(NULL)

    # Combine embeddings
    combined <- do.call(cbind, embeddings)
    combined <- scale(combined)

    # Run consensus clustering
    result <- tryCatch({
        ConsensusClusterPlus::ConsensusClusterPlus(
            t(combined),
            maxK = min(6, floor(length(common_samples) / 3)),
            reps = 100,
            pItem = 0.8,
            pFeature = 1,
            clusterAlg = "hc",
            distance = "euclidean",
            seed = 42,
            plot = NULL,
            verbose = FALSE
        )
    }, error = function(e) {
        message("  Consensus clustering failed: ", e$message)
        NULL
    })

    if (is.null(result)) return(NULL)

    # Get best K
    k_best <- 2
    for (k in 3:min(6, length(result))) {
        if (!is.null(result[[k]])) {
            if (mean(result[[k]]$consensusMatrix) > 0.7) {
                k_best <- k
            }
        }
    }

    list(
        k_best = k_best,
        clusters = result[[k_best]]$consensusClass,
        consensus_matrix = result[[k_best]]$consensusMatrix
    )
}

#' Identify discordant samples
identify_discordant_samples <- function(results, common_samples, fc) {
    message("  Identifying discordant samples...")

    if (is.null(results$clustering_consistency)) {
        return(NULL)
    }

    cluster_results <- results$clustering_consistency$cluster_results
    if (length(cluster_results) < 2) return(NULL)

    omics_names <- names(cluster_results)

    # For each sample, check if it clusters consistently
    discordance_scores <- sapply(common_samples, function(s) {
        clusters <- sapply(cluster_results, function(cr) cr$clusters[s])

        # Count how many omics pairs disagree
        pairs <- utils::combn(omics_names, 2, simplify = FALSE)
        n_disagree <- sum(sapply(pairs, function(pair) {
            clusters[pair[1]] != clusters[pair[2]]
        }))

        n_disagree / length(pairs)
    })

    # Identify highly discordant samples (>50% disagreement)
    discordant_samples <- names(discordance_scores)[discordance_scores > 0.5]

    discordance_df <- data.frame(
        sample_id = common_samples,
        discordance_score = discordance_scores,
        is_discordant = discordance_scores > 0.5,
        stringsAsFactors = FALSE
    )

    message("  Found ", length(discordant_samples), " discordant samples")

    list(
        discordance_df = discordance_df,
        discordant_samples = discordant_samples,
        n_discordant = length(discordant_samples)
    )
}

#' Save sample concordance results
save_sample_concordance_results <- function(results, out_dir) {
    # Save clustering consistency
    if (!is.null(results$clustering_consistency)) {
        consistency_df <- data.frame(
            pair = names(results$clustering_consistency$ari),
            ARI = results$clustering_consistency$ari,
            NMI = results$clustering_consistency$nmi,
            stringsAsFactors = FALSE
        )
        write.csv(consistency_df,
                  file.path(out_dir, "tables", "sample_clustering_consistency.csv"),
                  row.names = FALSE)
    }

    # Save condition similarity
    if (!is.null(results$condition_similarity)) {
        write.csv(results$condition_similarity$summary,
                  file.path(out_dir, "tables", "sample_condition_similarity.csv"),
                  row.names = FALSE)
    }

    # Save discordant samples
    if (!is.null(results$discordant_samples)) {
        write.csv(results$discordant_samples$discordance_df,
                  file.path(out_dir, "tables", "sample_discordance.csv"),
                  row.names = FALSE)
    }

    # Overall summary
    summary_df <- data.frame(
        metric = c("mantel_correlation", "mean_ari", "mean_nmi", "n_discordant"),
        value = c(
            if (!is.null(results$sample_rank_cors)) mean(results$sample_rank_cors$mantel_correlations, na.rm = TRUE) else NA,
            if (!is.null(results$clustering_consistency)) mean(results$clustering_consistency$ari, na.rm = TRUE) else NA,
            if (!is.null(results$clustering_consistency)) mean(results$clustering_consistency$nmi, na.rm = TRUE) else NA,
            results$discordant_samples$n_discordant %||% NA
        ),
        stringsAsFactors = FALSE
    )
    write.csv(summary_df,
              file.path(out_dir, "tables", "sample_concordance_summary.csv"),
              row.names = FALSE)
}

#' Plot sample concordance
plot_sample_concordance <- function(results, metadata, out_dir) {
    # Plot 0: Mantel test — per-omics sample correlation heatmaps with
    # dendrograms, plus a summary panel with test statistics
    sc <- results$sample_rank_cors
    if (!is.null(sc) && !is.null(sc$mantel_correlations)) {
        mantel_r <- sc$mantel_correlations
        mantel_p <- sc$mantel_pvalues

        # Build symmetric matrix of Mantel correlations
        pair_names <- names(mantel_r)
        omics_names <- unique(unlist(strsplit(pair_names, "_vs_")))
        n_omics <- length(omics_names)
        mat <- matrix(1, nrow = n_omics, ncol = n_omics,
                       dimnames = list(omics_names, omics_names))
        pmat <- matrix(0, nrow = n_omics, ncol = n_omics,
                        dimnames = list(omics_names, omics_names))

        for (i in seq_along(pair_names)) {
            parts <- strsplit(pair_names[i], "_vs_")[[1]]
            mat[parts[1], parts[2]] <- mantel_r[i]
            mat[parts[2], parts[1]] <- mantel_r[i]
            pmat[parts[1], parts[2]] <- mantel_p[i]
            pmat[parts[2], parts[1]] <- mantel_p[i]
        }

        dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

        # --- Per-omics sample distance/correlation heatmaps with dendrograms ---
        dist_mats <- sc$distance_matrices
        if (!is.null(dist_mats) && length(dist_mats) > 0) {
            for (om_name in names(dist_mats)) {
                om_mat <- dist_mats[[om_name]]
                # Compute sample-sample correlation matrix
                cor_mat <- cor(om_mat, method = "spearman")
                # Hierarchical clustering for dendrogram
                hc <- hclust(as.dist(1 - cor_mat), method = "ward.D2")

                out_file <- file.path(out_dir, "plots",
                                       paste0("mantel_corr_heatmap_", om_name, ".png"))

                if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
                    col_fun <- circlize::colorRamp2(
                        c(-1, 0, 1), c("steelblue", "white", "firebrick")
                    )
                    ht <- ComplexHeatmap::Heatmap(
                        cor_mat,
                        name = "Spearman r",
                        col = col_fun,
                        cluster_rows = hc,
                        cluster_columns = hc,
                        show_row_dend = TRUE,
                        show_column_dend = TRUE,
                        row_dend_width = grid::unit(25, "mm"),
                        column_dend_height = grid::unit(25, "mm"),
                        row_names_gp = grid::gpar(fontsize = 8),
                        column_names_gp = grid::gpar(fontsize = 8),
                        column_title = paste0(om_name, " — Sample Correlation Matrix"),
                        column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                            grid::grid.text(
                                sprintf("%.2f", cor_mat[i, j]),
                                x, y, gp = grid::gpar(fontsize = 6)
                            )
                        }
                    )
                    n_samples <- ncol(cor_mat)
                    plot_size <- max(6, 3 + n_samples * 0.35)
                    png(out_file, width = plot_size, height = plot_size,
                        units = "in", res = 150)
                    ComplexHeatmap::draw(ht)
                    dev.off()
                } else {
                    # Fallback: pheatmap
                    png(out_file, width = 8, height = 8, units = "in", res = 150)
                    pheatmap::pheatmap(
                        cor_mat,
                        color = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
                        clustering_method = "ward.D2",
                        display_numbers = TRUE,
                        number_format = "%.2f",
                        fontsize_number = 6,
                        fontsize_row = 8, fontsize_col = 8,
                        main = paste0(om_name, " — Sample Correlation Matrix")
                    )
                    dev.off()
                }
                message("  Saved ", om_name, " correlation heatmap: ", out_file)
            }
        }

        # --- Mantel test summary heatmap (existing tile plot, enhanced) ---
        df_tile <- expand.grid(Omics1 = omics_names, Omics2 = omics_names,
                               stringsAsFactors = FALSE)
        df_tile$r <- mapply(function(o1, o2) mat[o1, o2], df_tile$Omics1, df_tile$Omics2)
        df_tile$p <- mapply(function(o1, o2) pmat[o1, o2], df_tile$Omics1, df_tile$Omics2)
        df_tile$label <- ifelse(
            df_tile$Omics1 == df_tile$Omics2, "1.000",
            sprintf("%.3f\n(p=%s)", df_tile$r,
                    ifelse(df_tile$p < 0.001, "<0.001", format(round(df_tile$p, 3), nsmall = 3)))
        )

        # Build test results text for annotation
        sig_labels <- ifelse(mantel_p < 0.001, "***",
                      ifelse(mantel_p < 0.01, "**",
                      ifelse(mantel_p < 0.05, "*", "ns")))
        results_text <- paste(
            pair_names, ": r =", sprintf("%.3f", mantel_r),
            ", p =", sprintf("%.4f", mantel_p), sig_labels,
            collapse = "\n"
        )

        p0 <- ggplot2::ggplot(df_tile, ggplot2::aes(x = Omics1, y = Omics2, fill = r)) +
            ggplot2::geom_tile(color = "white", linewidth = 1) +
            ggplot2::geom_text(ggplot2::aes(label = label), size = 4, color = "black") +
            ggplot2::scale_fill_gradient2(
                low = "steelblue", mid = "white", high = "firebrick",
                midpoint = 0, limits = c(-1, 1), name = "Mantel r"
            ) +
            ggplot2::theme_minimal(base_size = 14) +
            ggplot2::labs(
                title = "Sample Distance Concordance (Mantel Test)",
                subtitle = "Correlation between sample distance matrices (999 permutations)",
                x = NULL, y = NULL,
                caption = paste0("Mantel test results:\n", results_text,
                                 "\nSignificance: *** p<0.001, ** p<0.01, * p<0.05")
            ) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                panel.grid = ggplot2::element_blank(),
                plot.caption = ggplot2::element_text(hjust = 0, size = 10,
                                                      face = "italic", lineheight = 1.3)
            ) +
            ggplot2::coord_fixed()

        ggplot2::ggsave(file.path(out_dir, "plots", "mantel_test_heatmap.png"),
                        p0, width = 8, height = 8, dpi = 150)
        message("  Saved Mantel test heatmap: ", file.path(out_dir, "plots", "mantel_test_heatmap.png"))

        # Save Mantel results as CSV
        mantel_df <- data.frame(
            pair = pair_names,
            mantel_r = mantel_r,
            mantel_p = mantel_p,
            stringsAsFactors = FALSE
        )
        write.csv(mantel_df, file.path(out_dir, "tables", "mantel_test_results.csv"),
                  row.names = FALSE)
    }

    # Plot 1: Clustering consistency bar plot
    if (!is.null(results$clustering_consistency)) {
        df <- data.frame(
            pair = names(results$clustering_consistency$ari),
            ARI = results$clustering_consistency$ari,
            NMI = results$clustering_consistency$nmi
        )
        df_long <- tidyr::pivot_longer(df, cols = c("ARI", "NMI"),
                                        names_to = "metric", values_to = "value")

        p1 <- ggplot2::ggplot(df_long, ggplot2::aes(x = pair, y = value, fill = metric)) +
            ggplot2::geom_bar(stat = "identity", position = "dodge") +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = "Sample Clustering Consistency Across Omics",
                subtitle = "Higher values indicate more consistent sample clustering",
                x = "Omics Pair",
                y = "Score",
                fill = "Metric"
            ) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
            ggplot2::ylim(0, 1)

        ggplot2::ggsave(file.path(out_dir, "plots", "sample_clustering_consistency.png"),
                        p1, width = 8, height = 6)
    }

    # Plot 2: Condition separation
    if (!is.null(results$condition_similarity)) {
        df <- results$condition_similarity$summary

        p2 <- ggplot2::ggplot(df, ggplot2::aes(x = omics, y = separation, fill = omics)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = "Within vs Between-Condition Similarity",
                subtitle = "Positive = samples cluster by condition within omics",
                x = "Omics",
                y = "Separation (within - between)"
            ) +
            ggplot2::theme(legend.position = "none")

        ggplot2::ggsave(file.path(out_dir, "plots", "sample_condition_separation.png"),
                        p2, width = 6, height = 5)
    }

    # Plot 3: Discordance distribution
    if (!is.null(results$discordant_samples)) {
        df <- results$discordant_samples$discordance_df

        p3 <- ggplot2::ggplot(df, ggplot2::aes(x = discordance_score)) +
            ggplot2::geom_histogram(bins = 20, fill = "steelblue", color = "white") +
            ggplot2::geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = "Sample Discordance Distribution",
                subtitle = "Samples right of red line cluster differently across omics",
                x = "Discordance Score",
                y = "Count"
            )

        ggplot2::ggsave(file.path(out_dir, "plots", "sample_discordance_distribution.png"),
                        p3, width = 8, height = 5)
    }
}


# =============================================================================
# 3. Pathway Overlap Analysis
# =============================================================================

#' Analyze pathway overlap across omics
#' @param harmonized List of harmonized omics data
#' @param gene_mapping Gene mapping table
#' @param fc Foundational config
#' @param out_dir Output directory
analyze_pathway_overlap <- function(harmonized, gene_mapping, fc, out_dir = NULL) {
    message("Analyzing pathway overlap across omics...")

    omics_names <- names(harmonized)

    # Check if enrichment results are available
    enrichment_available <- sapply(harmonized, function(h) {
        !is.null(h$pathway_results) || !is.null(h$enrichment_results)
    })

    if (!any(enrichment_available)) {
        message("  No pre-computed pathway results available.")
        message("  Running basic enrichment on DE/DA results...")
        pathway_results <- run_basic_enrichment(harmonized, gene_mapping, fc)
    } else {
        pathway_results <- lapply(names(enrichment_available)[enrichment_available], function(on) {
            h <- harmonized[[on]]
            list(
                omics = on,
                pathways = h$pathway_results %||% h$enrichment_results
            )
        })
        names(pathway_results) <- names(enrichment_available)[enrichment_available]
    }

    if (length(pathway_results) < 2) {
        message("  Insufficient omics with pathway results for overlap analysis")
        return(NULL)
    }

    results <- list()

    # 1. Pathway overlap matrix (Jaccard similarity)
    results$overlap_matrix <- compute_pathway_overlap_matrix(pathway_results, fc)

    # 2. Shared vs omics-specific pathways
    results$shared_specific <- identify_shared_specific_pathways(pathway_results, fc)

    # 3. Pathway direction concordance
    results$direction_concordance <- compute_pathway_direction_concordance(pathway_results, fc)

    # Save and visualize
    if (!is.null(out_dir)) {
        save_pathway_overlap_results(results, out_dir)
        plot_pathway_overlap(results, out_dir)
    }

    return(results)
}

#' Run basic enrichment for omics with DE/DA results
run_basic_enrichment <- function(harmonized, gene_mapping, fc) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        message("  clusterProfiler not available. Skipping enrichment.")
        return(list())
    }

    omics_names <- names(harmonized)

    results <- lapply(omics_names, function(omics_name) {
        h <- harmonized[[omics_name]]
        de_table <- h$de_table %||% h$da_table

        if (is.null(de_table)) return(NULL)

        # Get significant genes
        padj_col <- intersect(c("padj", "adj.P.Val", "FDR"), colnames(de_table))[1]
        if (is.na(padj_col)) return(NULL)

        sig_genes <- de_table$gene_symbol[de_table[[padj_col]] < 0.05]
        sig_genes <- sig_genes[!is.na(sig_genes) & sig_genes != ""]

        if (length(sig_genes) < 10) return(NULL)

        # Run GO enrichment
        enrich_result <- tryCatch({
            clusterProfiler::enrichGO(
                gene = sig_genes,
                OrgDb = "org.Hs.eg.db",
                keyType = "SYMBOL",
                ont = "BP",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.1
            )
        }, error = function(e) NULL)

        if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) == 0) {
            return(NULL)
        }

        enrich_df <- as.data.frame(enrich_result)

        list(
            omics = omics_name,
            pathways = enrich_df,
            sig_genes = sig_genes
        )
    })
    names(results) <- omics_names
    results <- results[!sapply(results, is.null)]

    results
}

#' Compute pathway overlap matrix (Jaccard similarity)
compute_pathway_overlap_matrix <- function(pathway_results, fc) {
    omics_names <- names(pathway_results)
    n_omics <- length(omics_names)

    # Get significant pathways per omics
    sig_pathways <- lapply(pathway_results, function(pr) {
        df <- pr$pathways
        padj_col <- intersect(c("p.adjust", "qvalue", "padj"), colnames(df))[1]
        if (is.na(padj_col)) return(character(0))

        term_col <- intersect(c("ID", "Description", "pathway"), colnames(df))[1]
        if (is.na(term_col)) return(character(0))

        df[[term_col]][df[[padj_col]] < fc$fdr_threshold]
    })

    # Compute Jaccard similarity matrix
    jaccard_mat <- matrix(0, n_omics, n_omics)
    rownames(jaccard_mat) <- colnames(jaccard_mat) <- omics_names

    for (i in 1:n_omics) {
        for (j in 1:n_omics) {
            a <- sig_pathways[[i]]
            b <- sig_pathways[[j]]

            if (length(a) == 0 || length(b) == 0) {
                jaccard_mat[i, j] <- 0
            } else {
                intersection <- length(intersect(a, b))
                union_size <- length(union(a, b))
                jaccard_mat[i, j] <- intersection / union_size
            }
        }
    }

    # Hypergeometric test for overlap significance
    overlap_tests <- list()
    pairs <- utils::combn(omics_names, 2, simplify = FALSE)

    all_pathways <- unique(unlist(sig_pathways))
    universe_size <- length(all_pathways)

    for (pair in pairs) {
        a <- sig_pathways[[pair[1]]]
        b <- sig_pathways[[pair[2]]]

        if (length(a) > 0 && length(b) > 0 && universe_size > 0) {
            overlap <- length(intersect(a, b))

            pval <- phyper(overlap - 1, length(a), universe_size - length(a),
                           length(b), lower.tail = FALSE)

            overlap_tests[[paste(pair, collapse = "_vs_")]] <- list(
                overlap = overlap,
                size1 = length(a),
                size2 = length(b),
                jaccard = jaccard_mat[pair[1], pair[2]],
                pvalue = pval
            )
        }
    }

    list(
        jaccard_matrix = jaccard_mat,
        significant_pathways = sig_pathways,
        overlap_tests = overlap_tests
    )
}

#' Identify shared vs omics-specific pathways
identify_shared_specific_pathways <- function(pathway_results, fc) {
    omics_names <- names(pathway_results)

    # Get significant pathways per omics
    sig_pathways <- lapply(pathway_results, function(pr) {
        df <- pr$pathways
        padj_col <- intersect(c("p.adjust", "qvalue", "padj"), colnames(df))[1]
        term_col <- intersect(c("ID", "Description", "pathway"), colnames(df))[1]

        if (is.na(padj_col) || is.na(term_col)) return(character(0))
        df[[term_col]][df[[padj_col]] < fc$fdr_threshold]
    })

    # All pathways across all omics
    all_pathways <- unique(unlist(sig_pathways))

    if (length(all_pathways) == 0) {
        return(NULL)
    }

    # Classify each pathway
    pathway_classification <- data.frame(
        pathway = all_pathways,
        stringsAsFactors = FALSE
    )

    for (omics in omics_names) {
        pathway_classification[[omics]] <- all_pathways %in% sig_pathways[[omics]]
    }

    # Number of omics enriched
    pathway_classification$n_omics <- rowSums(pathway_classification[, omics_names, drop = FALSE])

    # Classification
    pathway_classification$category <- ifelse(
        pathway_classification$n_omics == length(omics_names), "shared_all",
        ifelse(pathway_classification$n_omics > 1, "shared_some", "omics_specific")
    )

    # Summary
    summary_df <- data.frame(
        category = c("shared_all", "shared_some", "omics_specific"),
        count = c(
            sum(pathway_classification$category == "shared_all"),
            sum(pathway_classification$category == "shared_some"),
            sum(pathway_classification$category == "omics_specific")
        ),
        stringsAsFactors = FALSE
    )

    message("  Pathway classification: ",
            sum(pathway_classification$category == "shared_all"), " shared by all, ",
            sum(pathway_classification$category == "omics_specific"), " omics-specific")

    list(
        classification = pathway_classification,
        summary = summary_df,
        shared_pathways = all_pathways[pathway_classification$category != "omics_specific"]
    )
}

#' Compute pathway direction concordance
compute_pathway_direction_concordance <- function(pathway_results, fc) {
    omics_names <- names(pathway_results)

    # Check if we have direction information
    has_direction <- sapply(pathway_results, function(pr) {
        "NES" %in% colnames(pr$pathways) || "direction" %in% colnames(pr$pathways)
    })

    if (!any(has_direction)) {
        message("  No direction information available for pathway concordance")
        return(NULL)
    }

    # Get pathways with direction
    pathway_directions <- lapply(pathway_results, function(pr) {
        df <- pr$pathways

        if ("NES" %in% colnames(df)) {
            term_col <- intersect(c("ID", "Description", "pathway"), colnames(df))[1]
            data.frame(
                pathway = df[[term_col]],
                direction = sign(df$NES),
                stringsAsFactors = FALSE
            )
        } else if ("direction" %in% colnames(df)) {
            term_col <- intersect(c("ID", "Description", "pathway"), colnames(df))[1]
            data.frame(
                pathway = df[[term_col]],
                direction = ifelse(df$direction == "up", 1, -1),
                stringsAsFactors = FALSE
            )
        } else {
            NULL
        }
    })
    pathway_directions <- pathway_directions[!sapply(pathway_directions, is.null)]

    if (length(pathway_directions) < 2) {
        return(NULL)
    }

    # Find common pathways and check direction concordance
    common_pathways <- Reduce(intersect, lapply(pathway_directions, function(x) x$pathway))

    if (length(common_pathways) == 0) {
        return(NULL)
    }

    # Build direction matrix
    direction_mat <- matrix(NA, length(common_pathways), length(pathway_directions))
    rownames(direction_mat) <- common_pathways
    colnames(direction_mat) <- names(pathway_directions)

    for (omics in names(pathway_directions)) {
        pd <- pathway_directions[[omics]]
        idx <- match(common_pathways, pd$pathway)
        direction_mat[, omics] <- pd$direction[idx]
    }

    # Concordance: same direction across all omics
    concordance <- apply(direction_mat, 1, function(row) {
        row <- row[!is.na(row)]
        if (length(row) < 2) return(NA)
        all(row == row[1])
    })

    concordance_df <- data.frame(
        pathway = common_pathways,
        concordant = concordance,
        stringsAsFactors = FALSE
    )

    pct_concordant <- 100 * mean(concordance, na.rm = TRUE)
    message("  Pathway direction concordance: ", round(pct_concordant, 1), "%")

    list(
        direction_matrix = direction_mat,
        concordance = concordance_df,
        pct_concordant = pct_concordant
    )
}

#' Save pathway overlap results
save_pathway_overlap_results <- function(results, out_dir) {
    if (!is.null(results$overlap_matrix)) {
        # Save Jaccard matrix
        jaccard_df <- as.data.frame(results$overlap_matrix$jaccard_matrix)
        jaccard_df$omics <- rownames(jaccard_df)
        write.csv(jaccard_df,
                  file.path(out_dir, "tables", "pathway_overlap_jaccard.csv"),
                  row.names = FALSE)

        # Save overlap tests
        if (length(results$overlap_matrix$overlap_tests) > 0) {
            tests_df <- do.call(rbind, lapply(names(results$overlap_matrix$overlap_tests), function(n) {
                t <- results$overlap_matrix$overlap_tests[[n]]
                data.frame(
                    pair = n,
                    overlap = t$overlap,
                    size1 = t$size1,
                    size2 = t$size2,
                    jaccard = round(t$jaccard, 4),
                    pvalue = t$pvalue,
                    stringsAsFactors = FALSE
                )
            }))
            write.csv(tests_df,
                      file.path(out_dir, "tables", "pathway_overlap_tests.csv"),
                      row.names = FALSE)
        }
    }

    if (!is.null(results$shared_specific)) {
        write.csv(results$shared_specific$classification,
                  file.path(out_dir, "tables", "pathway_shared_specific.csv"),
                  row.names = FALSE)
    }

    if (!is.null(results$direction_concordance)) {
        write.csv(results$direction_concordance$concordance,
                  file.path(out_dir, "tables", "pathway_direction_concordance.csv"),
                  row.names = FALSE)
    }
}

#' Plot pathway overlap
plot_pathway_overlap <- function(results, out_dir) {
    # Plot 1: Jaccard similarity heatmap
    if (!is.null(results$overlap_matrix)) {
        jaccard_mat <- results$overlap_matrix$jaccard_matrix

        jaccard_df <- as.data.frame(as.table(jaccard_mat))
        colnames(jaccard_df) <- c("Omics1", "Omics2", "Jaccard")

        p1 <- ggplot2::ggplot(jaccard_df, ggplot2::aes(x = Omics1, y = Omics2, fill = Jaccard)) +
            ggplot2::geom_tile() +
            ggplot2::geom_text(ggplot2::aes(label = round(Jaccard, 2)), color = "white", size = 4) +
            ggplot2::scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = "Pathway Overlap Across Omics",
                subtitle = "Jaccard similarity of enriched pathways",
                x = "", y = "",
                fill = "Jaccard\nSimilarity"
            ) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

        ggplot2::ggsave(file.path(out_dir, "plots", "pathway_overlap_heatmap.png"),
                        p1, width = 7, height = 6)
    }

    # Plot 2: Shared vs specific pie chart
    if (!is.null(results$shared_specific)) {
        summary_df <- results$shared_specific$summary

        p3 <- ggplot2::ggplot(summary_df, ggplot2::aes(x = "", y = count, fill = category)) +
            ggplot2::geom_bar(stat = "identity", width = 1) +
            ggplot2::coord_polar("y", start = 0) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = "Pathway Classification",
                fill = "Category"
            ) +
            ggplot2::scale_fill_manual(values = c("shared_all" = "#2166ac",
                                                   "shared_some" = "#67a9cf",
                                                   "omics_specific" = "#d1e5f0"))

        ggplot2::ggsave(file.path(out_dir, "plots", "pathway_classification_pie.png"),
                        p3, width = 7, height = 6)
    }
}


# =============================================================================
# 4. Cross-Omics Module Detection
# =============================================================================

#' Find cross-omics co-expression modules
#' @param harmonized List of harmonized omics data
#' @param gene_mapping Gene mapping table
#' @param fc Foundational config
#' @param out_dir Output directory
find_crossomics_modules <- function(harmonized, gene_mapping, fc, out_dir = NULL) {
    message("Finding cross-omics co-expression modules...")

    # Check if WGCNA is available
    if (!requireNamespace("WGCNA", quietly = TRUE)) {
        message("  WGCNA not available. Skipping module detection.")
        return(NULL)
    }

    omics_names <- names(harmonized)

    # Get common samples
    all_samples <- lapply(harmonized, function(h) colnames(h$normalized_matrix))
    common_samples <- Reduce(intersect, all_samples)

    if (length(common_samples) < 10) {
        message("  Insufficient common samples for module detection")
        return(NULL)
    }

    # Map features to common identifiers across omics
    combined_data <- list()
    feature_info <- data.frame()

    for (omics in omics_names) {
        mat <- harmonized[[omics]]$normalized_matrix[, common_samples, drop = FALSE]

        # Use top variable features
        vars <- apply(mat, 1, var, na.rm = TRUE)
        top_idx <- order(vars, decreasing = TRUE)[1:min(fc$top_variable_features, nrow(mat))]
        mat_sub <- mat[top_idx, ]

        # Prefix feature names
        rownames(mat_sub) <- paste0(omics, "_", rownames(mat_sub))

        combined_data[[omics]] <- mat_sub

        # Track feature info
        feature_info <- rbind(feature_info, data.frame(
            feature_id = rownames(mat_sub),
            omics = omics,
            original_id = gsub(paste0("^", omics, "_"), "", rownames(mat_sub)),
            stringsAsFactors = FALSE
        ))
    }

    # Combine matrices
    combined_mat <- do.call(rbind, combined_data)
    combined_mat <- combined_mat[complete.cases(combined_mat), ]

    if (nrow(combined_mat) < 100) {
        message("  Insufficient features after filtering for module detection")
        return(NULL)
    }

    message("  Running WGCNA on ", nrow(combined_mat), " features across ",
            length(common_samples), " samples")

    # Transpose for WGCNA (samples in rows)
    expr_data <- t(combined_mat)

    # Pick soft threshold
    powers <- c(1:10, seq(12, 20, 2))
    sft <- tryCatch({
        WGCNA::pickSoftThreshold(
            expr_data,
            powerVector = powers,
            verbose = 0
        )
    }, error = function(e) NULL)

    if (is.null(sft)) {
        message("  Soft threshold selection failed")
        return(NULL)
    }

    # Use recommended power or default
    soft_power <- sft$powerEstimate
    if (is.na(soft_power)) {
        soft_power <- 6
        message("  Using default soft power: ", soft_power)
    } else {
        message("  Selected soft power: ", soft_power)
    }

    # Build network and detect modules
    net <- tryCatch({
        WGCNA::blockwiseModules(
            expr_data,
            power = soft_power,
            TOMType = "unsigned",
            minModuleSize = 30,
            mergeCutHeight = 0.25,
            numericLabels = TRUE,
            verbose = 0
        )
    }, error = function(e) {
        message("  Module detection failed: ", e$message)
        NULL
    })

    if (is.null(net)) return(NULL)

    # Process results
    module_labels <- net$colors
    module_colors <- WGCNA::labels2colors(module_labels)

    # Create module membership table
    module_df <- data.frame(
        feature_id = colnames(expr_data),
        module_number = module_labels,
        module_color = module_colors,
        stringsAsFactors = FALSE
    )

    # Add omics information
    module_df <- merge(module_df, feature_info, by = "feature_id", all.x = TRUE)

    # Module summary
    module_summary <- table(module_df$module_number, module_df$omics)

    # Calculate omics composition per module
    module_composition <- as.data.frame.matrix(module_summary)
    module_composition$module <- rownames(module_composition)
    module_composition$total <- rowSums(module_composition[, omics_names])

    # Identify cross-omics modules
    n_omics_per_module <- apply(module_composition[, omics_names], 1, function(x) sum(x > 0))
    crossomics_modules <- names(n_omics_per_module)[n_omics_per_module > 1]

    message("  Detected ", length(unique(module_labels)), " modules, ",
            length(crossomics_modules), " contain features from multiple omics")

    # Module eigengenes
    ME <- net$MEs
    colnames(ME) <- gsub("^ME", "module_", colnames(ME))

    # Save results
    if (!is.null(out_dir)) {
        write.csv(module_df,
                  file.path(out_dir, "tables", "crossomics_modules.csv"),
                  row.names = FALSE)
        write.csv(module_composition,
                  file.path(out_dir, "tables", "crossomics_module_composition.csv"),
                  row.names = FALSE)

        ME_df <- as.data.frame(ME)
        ME_df$sample_id <- rownames(ME_df)
        write.csv(ME_df,
                  file.path(out_dir, "tables", "crossomics_module_eigengenes.csv"),
                  row.names = FALSE)

        # Plot
        plot_module_composition(module_composition, omics_names, out_dir)
    }

    list(
        module_assignments = module_df,
        module_composition = module_composition,
        module_eigengenes = ME,
        n_modules = length(unique(module_labels)),
        n_crossomics = length(crossomics_modules),
        soft_power = soft_power
    )
}

#' Plot module composition
plot_module_composition <- function(module_composition, omics_names, out_dir) {
    # Prepare data for stacked bar plot
    mc_long <- tidyr::pivot_longer(
        module_composition,
        cols = tidyr::all_of(omics_names),
        names_to = "omics",
        values_to = "count"
    )

    # Filter out grey module (unassigned)
    mc_long <- mc_long[mc_long$module != "0", ]

    if (nrow(mc_long) == 0) return(NULL)

    p <- ggplot2::ggplot(mc_long, ggplot2::aes(x = module, y = count, fill = omics)) +
        ggplot2::geom_bar(stat = "identity", position = "stack") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "Cross-Omics Module Composition",
            subtitle = "Features per omics in each detected module",
            x = "Module",
            y = "Number of Features",
            fill = "Omics"
        ) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    ggplot2::ggsave(file.path(out_dir, "plots", "crossomics_module_composition.png"),
                    p, width = 10, height = 6)
}


# =============================================================================
# Summary Function
# =============================================================================

#' Summarize foundational analysis results
summarize_foundational_results <- function(results, out_dir = NULL) {
    message("Generating foundational analysis summary...")

    summary_list <- list()

    # Feature correlations summary
    if (!is.null(results$feature_correlations)) {
        cor_summary <- results$feature_correlations$summary
        if (!is.null(cor_summary) && nrow(cor_summary) > 0) {
            summary_list$correlations <- list(
                n_pairs_analyzed = nrow(cor_summary),
                mean_correlation = mean(cor_summary$mean_correlation, na.rm = TRUE),
                total_significant = sum(cor_summary$n_significant, na.rm = TRUE)
            )
        }
    }

    # Sample concordance summary
    if (!is.null(results$sample_concordance)) {
        sc <- results$sample_concordance
        summary_list$sample_concordance <- list(
            mean_ari = if (!is.null(sc$clustering_consistency)) mean(sc$clustering_consistency$ari, na.rm = TRUE) else NA,
            mean_nmi = if (!is.null(sc$clustering_consistency)) mean(sc$clustering_consistency$nmi, na.rm = TRUE) else NA,
            n_discordant = sc$discordant_samples$n_discordant %||% NA
        )
    }

    # Pathway overlap summary
    if (!is.null(results$pathway_overlap)) {
        po <- results$pathway_overlap
        if (!is.null(po$shared_specific)) {
            summary_list$pathway_overlap <- list(
                n_shared_all = sum(po$shared_specific$classification$category == "shared_all"),
                n_omics_specific = sum(po$shared_specific$classification$category == "omics_specific"),
                mean_jaccard = mean(po$overlap_matrix$jaccard_matrix[upper.tri(po$overlap_matrix$jaccard_matrix)])
            )
        }
    }

    # Cross-omics modules summary
    if (!is.null(results$crossomics_modules)) {
        cm <- results$crossomics_modules
        summary_list$modules <- list(
            n_modules = cm$n_modules,
            n_crossomics = cm$n_crossomics
        )
    }

    # Save summary
    if (!is.null(out_dir) && length(summary_list) > 0) {
        flat_summary <- do.call(rbind, lapply(names(summary_list), function(cat) {
            s <- summary_list[[cat]]
            data.frame(
                category = cat,
                metric = names(s),
                value = unlist(s),
                stringsAsFactors = FALSE
            )
        }))

        write.csv(flat_summary,
                  file.path(out_dir, "tables", "foundational_analysis_summary.csv"),
                  row.names = FALSE)
    }

    message("=== Foundational Analysis Summary ===")
    for (cat in names(summary_list)) {
        message("  ", cat, ": ", paste(names(summary_list[[cat]]), "=",
                                        round(unlist(summary_list[[cat]]), 3),
                                        collapse = ", "))
    }

    summary_list
}


# =============================================================================
# Internal: MAE to Legacy Adapter
# =============================================================================

#' Convert MAE to legacy mae_data format (internal)
#'
#' @param mae MultiAssayExperiment object
#' @param de_results Named list of DE results per omics
#' @param gene_protein_mapping Gene-protein mapping table
#' @return List in legacy mae_data format
.mae_to_legacy_foundational <- function(mae, de_results = NULL, gene_protein_mapping = NULL) {

    harmonized_omics <- lapply(names(mae@ExperimentList), function(nm) {
        exp_data <- mae@ExperimentList[[nm]]
        de <- if (!is.null(de_results) && nm %in% names(de_results)) de_results[[nm]] else NULL

        list(
            normalized_matrix = as.matrix(SummarizedExperiment::assay(exp_data)),
            de_table = de,
            da_table = de,
            feature_annotation = as.data.frame(SummarizedExperiment::rowData(exp_data)),
            metadata = as.data.frame(SummarizedExperiment::colData(mae))
        )
    })
    names(harmonized_omics) <- names(mae@ExperimentList)

    # Convert gene_protein_mapping (gene_id, protein_id) to the format expected
    # by foundational analysis (omics, feature_id, gene_symbol)
    gene_mapping <- NULL
    if (!is.null(gene_protein_mapping) && nrow(gene_protein_mapping) > 0) {
        rna_rows <- data.frame(
            omics = "transcriptomics",
            feature_id = gene_protein_mapping$gene_id,
            gene_symbol = if ("gene_symbol" %in% colnames(gene_protein_mapping)) {
                gene_protein_mapping$gene_symbol
            } else {
                gene_protein_mapping$gene_id
            },
            stringsAsFactors = FALSE
        )
        prot_rows <- data.frame(
            omics = "proteomics",
            feature_id = gene_protein_mapping$protein_id,
            gene_symbol = if ("gene_symbol" %in% colnames(gene_protein_mapping)) {
                gene_protein_mapping$gene_symbol
            } else {
                gene_protein_mapping$protein_id
            },
            stringsAsFactors = FALSE
        )
        gene_mapping <- rbind(rna_rows, prot_rows)
    }

    list(
        mae = mae,
        harmonized_omics = harmonized_omics,
        metadata = as.data.frame(SummarizedExperiment::colData(mae)),
        common_samples = colnames(mae),
        gene_mapping = gene_mapping
    )
}
