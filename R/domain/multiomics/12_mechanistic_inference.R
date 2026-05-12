# =============================================================================
# Mechanistic Inference and Regulatory Network Analysis
# =============================================================================
# This script provides deep mechanistic insights from multi-omics data:
#   1. RNA-Protein regulatory analysis (translation efficiency, protein stability)
#   2. Transcription factor activity inference
#   3. Gene regulatory network inference
#   4. Causal mediation analysis (RNA -> Protein -> Phenotype)
# =============================================================================

#' Main function: Run mechanistic analysis
#'
#' @param mae MultiAssayExperiment object with aligned samples
#' @param de_results Named list of DE results per omics (optional)
#' @param gene_protein_mapping Gene-protein mapping table
#' @param foundational_results Results from foundational analysis (optional)
#' @param config Pipeline configuration
#' @param out_dir Output directory for results
#' @return List containing mechanistic analysis results
run_mechanistic_analysis <- function(mae, de_results = NULL,
                                      gene_protein_mapping = NULL,
                                      foundational_results = NULL,
                                      config, out_dir = NULL) {

    message("=== Running Mechanistic Inference Analysis ===")

    # Convert MAE to legacy format
    mae_data <- .mae_to_legacy_mechanistic(mae, de_results, gene_protein_mapping)

    # Get mechanistic config with defaults
    mc <- get_mechanistic_config(config)

    if (!mc$run_mechanistic) {
        message("Mechanistic analysis disabled in config. Skipping.")
        return(NULL)
    }

    # Create output directories
    if (!is.null(out_dir)) {
        dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
    }

    results <- list()
    harmonized <- mae_data$harmonized_omics

    # 1. RNA-Protein regulatory analysis
    if ("transcriptomics" %in% names(harmonized) &&
        "proteomics" %in% names(harmonized)) {
        results$rna_protein <- tryCatch({
            analyze_rna_protein_regulation(
                harmonized, mae_data$gene_mapping, mc, config, out_dir
            )
        }, error = function(e) {
            message("RNA-Protein regulatory analysis failed: ", e$message)
            NULL
        })
    }

    # 2. Transcription factor activity inference (RNA-seq)
    if ("transcriptomics" %in% names(harmonized) && mc$tf_activity) {
        results$tf_activity <- tryCatch({
            infer_tf_activity(
                harmonized$transcriptomics, mc, config, out_dir
            )
        }, error = function(e) {
            message("TF activity inference failed: ", e$message)
            NULL
        })
    }

    # 2b. Pathway activity inference from proteomics (PROGENy via decoupleR)
    if ("proteomics" %in% names(harmonized) && mc$tf_activity) {
        results$pathway_activity <- tryCatch({
            infer_pathway_activity_proteomics(
                harmonized$proteomics, mc, config, out_dir
            )
        }, error = function(e) {
            message("Pathway activity inference from proteomics failed: ", e$message)
            NULL
        })
    }

    # 3. Regulatory / co-expression network inference
    #    Uses transcriptomics if available, falls back to proteomics
    if (mc$infer_regulatory_network) {
        results$regulatory_network <- tryCatch({
            infer_regulatory_network(
                harmonized, mae_data$gene_mapping, mc, config, out_dir
            )
        }, error = function(e) {
            message("Regulatory network inference failed: ", e$message)
            NULL
        })
    }

    # 3b. Protein-metabolite correlation network (works without transcriptomics)
    if ("proteomics" %in% names(harmonized) &&
        "metabolomics" %in% names(harmonized)) {
        results$protein_metabolite_network <- tryCatch({
            analyze_protein_metabolite_correlations(
                harmonized, de_results, mc, config, out_dir
            )
        }, error = function(e) {
            message("Protein-metabolite correlation failed: ", e$message)
            NULL
        })
    }

    # 4. Causal mediation analysis (RNA -> Protein -> Phenotype)
    if (mc$run_mediation && !is.null(mae_data$metadata) &&
        "transcriptomics" %in% names(harmonized) &&
        "proteomics" %in% names(harmonized)) {
        results$mediation <- tryCatch({
            run_mediation_analysis(
                harmonized, mae_data$metadata, mae_data$gene_mapping, mc, config, out_dir
            )
        }, error = function(e) {
            message("Mediation analysis failed: ", e$message)
            NULL
        })
    }

    # 4b. Protein -> Metabolite -> Phenotype mediation
    #     Works WITHOUT transcriptomics — needs proteomics + metabolomics + outcome
    if (mc$run_mediation && !is.null(mae_data$metadata) &&
        "proteomics" %in% names(harmonized) &&
        "metabolomics" %in% names(harmonized)) {
        results$prot_metab_mediation <- tryCatch({
            run_protein_metabolite_mediation(
                harmonized, mae_data$metadata, mc, config, out_dir
            )
        }, error = function(e) {
            message("Protein-metabolite mediation failed: ", e$message)
            NULL
        })
    }

    # 5. Multi-layer mediation (RNA -> Protein -> Metabolite -> Phenotype)
    if (mc$run_mediation && !is.null(mae_data$metadata) &&
        "transcriptomics" %in% names(harmonized) &&
        "metabolomics" %in% names(harmonized)) {
        results$multilayer_mediation <- tryCatch({
            run_multilayer_mediation(
                harmonized, mae_data$metadata, mc, config, out_dir
            )
        }, error = function(e) {
            message("Multi-layer mediation failed: ", e$message)
            NULL
        })
    }

    # 6. COSMOS causal network (if cosmosR is available)
    if ((mc$run_cosmos %||% TRUE) && requireNamespace("cosmosR", quietly = TRUE)) {
        results$cosmos <- tryCatch({
            run_cosmos_causal_network(
                mae = mae_data$mae,
                de_results = de_results,
                config = config,
                out_dir = out_dir
            )
        }, error = function(e) {
            message("COSMOS causal network failed: ", e$message)
            NULL
        })
    }

    # Generate summary
    results$summary <- summarize_mechanistic_results(results, out_dir)

    message("=== Mechanistic Analysis Complete ===")
    return(results)
}

#' Get mechanistic analysis configuration with defaults
#' @param config Full config object
#' @return List of mechanistic config parameters
get_mechanistic_config <- function(config) {
    mc <- config$modes$multiomics$mechanistic %||% list()

    list(
        run_mechanistic = mc$run_mechanistic %||% TRUE,
        tf_activity = mc$tf_activity %||% TRUE,
        tf_database = mc$tf_database %||% "dorothea",
        infer_regulatory_network = mc$infer_regulatory_network %||% TRUE,
        network_method = mc$network_method %||% "genie3",
        run_mediation = mc$run_mediation %||% TRUE,
        run_cosmos = mc$run_cosmos %||% TRUE,
        n_top_regulators = mc$n_top_regulators %||% 50,
        fdr_threshold = mc$fdr_threshold %||% 0.05,
        min_targets = mc$min_targets %||% 10
    )
}


# =============================================================================
# 1. RNA-Protein Regulatory Analysis
# =============================================================================

#' Analyze RNA-Protein regulatory relationships
#' @param harmonized List of harmonized omics data
#' @param gene_mapping Gene mapping table
#' @param mc Mechanistic config
#' @param config Full config
#' @param out_dir Output directory
analyze_rna_protein_regulation <- function(harmonized, gene_mapping, mc, config, out_dir = NULL) {
    message("Analyzing RNA-Protein regulatory relationships...")

    rna_mat <- harmonized$transcriptomics$normalized_matrix
    prot_mat <- harmonized$proteomics$normalized_matrix

    # Get common samples
    common_samples <- intersect(colnames(rna_mat), colnames(prot_mat))
    if (length(common_samples) < 5) {
        message("  Insufficient common samples. Skipping.")
        return(NULL)
    }

    rna_mat <- rna_mat[, common_samples, drop = FALSE]
    prot_mat <- prot_mat[, common_samples, drop = FALSE]

    # Map to common gene symbols
    mapping_result <- map_rna_protein_features(rna_mat, prot_mat, gene_mapping)
    if (is.null(mapping_result)) {
        message("  Could not map RNA and protein features. Skipping.")
        return(NULL)
    }

    rna_matched <- mapping_result$rna
    prot_matched <- mapping_result$protein
    gene_symbols <- mapping_result$gene_symbols

    message("  Analyzing ", length(gene_symbols), " matched genes")

    results <- list()

    # 1. Per-gene RNA-Protein correlation analysis
    results$correlations <- compute_rna_protein_correlations(
        rna_matched, prot_matched, gene_symbols
    )

    # 2. Translation efficiency inference
    results$translation_efficiency <- compute_translation_efficiency(
        rna_matched, prot_matched, gene_symbols
    )

    # 3. Post-transcriptional regulation detection
    results$post_transcriptional <- detect_post_transcriptional_regulation(
        rna_matched, prot_matched, gene_symbols, mc
    )

    # 4. Protein stability inference
    results$protein_stability <- infer_protein_stability(
        rna_matched, prot_matched, gene_symbols
    )

    # Save results
    if (!is.null(out_dir)) {
        save_rna_protein_results(results, out_dir)
        plot_rna_protein_regulation(results, out_dir)
    }

    return(results)
}

#' Map RNA and protein features to common identifiers
map_rna_protein_features <- function(rna_mat, prot_mat, gene_mapping) {
    # Always try direct row-name matching first — works when MAE was
    # pre-harmonized (features renamed to GENE_1, GENE_2, ...)
    common <- intersect(rownames(rna_mat), rownames(prot_mat))
    if (length(common) >= 50) {
        message("  Direct row-name matching: ", length(common), " common features")
        return(list(
            rna = rna_mat[common, , drop = FALSE],
            protein = prot_mat[common, , drop = FALSE],
            gene_symbols = common
        ))
    }

    if (is.null(gene_mapping)) return(NULL)

    # Use gene mapping
    rna_map <- gene_mapping[gene_mapping$omics == "transcriptomics", ]
    prot_map <- gene_mapping[gene_mapping$omics == "proteomics", ]

    common_symbols <- intersect(rna_map$gene_symbol, prot_map$gene_symbol)
    common_symbols <- common_symbols[!is.na(common_symbols) & common_symbols != ""]

    if (length(common_symbols) < 50) return(NULL)

    # Get feature IDs
    rna_features <- rna_map$feature_id[match(common_symbols, rna_map$gene_symbol)]
    prot_features <- prot_map$feature_id[match(common_symbols, prot_map$gene_symbol)]

    # Filter to present features
    valid <- rna_features %in% rownames(rna_mat) & prot_features %in% rownames(prot_mat)
    common_symbols <- common_symbols[valid]
    rna_features <- rna_features[valid]
    prot_features <- prot_features[valid]

    if (length(common_symbols) < 50) return(NULL)

    list(
        rna = rna_mat[rna_features, , drop = FALSE],
        protein = prot_mat[prot_features, , drop = FALSE],
        gene_symbols = common_symbols
    )
}

#' Compute RNA-Protein correlations per gene
compute_rna_protein_correlations <- function(rna_mat, prot_mat, gene_symbols) {
    n_genes <- length(gene_symbols)

    correlations <- numeric(n_genes)
    pvalues <- numeric(n_genes)
    names(correlations) <- names(pvalues) <- gene_symbols

    for (i in seq_len(n_genes)) {
        rna_expr <- rna_mat[i, ]
        prot_expr <- prot_mat[i, ]

        if (sd(rna_expr, na.rm = TRUE) == 0 || sd(prot_expr, na.rm = TRUE) == 0) {
            correlations[i] <- NA
            pvalues[i] <- NA
            next
        }

        ct <- tryCatch({
            cor.test(rna_expr, prot_expr, method = "spearman", use = "pairwise.complete.obs")
        }, error = function(e) NULL)

        if (!is.null(ct)) {
            correlations[i] <- ct$estimate
            pvalues[i] <- ct$p.value
        } else {
            correlations[i] <- NA
            pvalues[i] <- NA
        }
    }

    padj <- p.adjust(pvalues, method = "BH")

    cor_df <- data.frame(
        gene_symbol = gene_symbols,
        rna_protein_cor = correlations,
        pvalue = pvalues,
        padj = padj,
        stringsAsFactors = FALSE
    )

    # Summary
    valid_cors <- correlations[!is.na(correlations)]
    summary_stats <- list(
        n_genes = n_genes,
        n_valid = length(valid_cors),
        mean_cor = mean(valid_cors),
        median_cor = median(valid_cors),
        n_positive = sum(valid_cors > 0),
        n_strong_positive = sum(valid_cors > 0.5),
        n_negative = sum(valid_cors < 0),
        pct_positive = 100 * mean(valid_cors > 0)
    )

    message("  RNA-Protein correlations: mean=", round(summary_stats$mean_cor, 3),
            ", ", summary_stats$n_strong_positive, " with r>0.5")

    list(
        correlations = correlations,
        cor_df = cor_df,
        summary = summary_stats
    )
}

#' Compute translation efficiency scores
compute_translation_efficiency <- function(rna_mat, prot_mat, gene_symbols) {
    # Translation efficiency = log2(protein) - log2(RNA)
    # Higher values indicate more efficient translation

    n_genes <- length(gene_symbols)

    # Compute mean expression
    rna_mean <- rowMeans(rna_mat, na.rm = TRUE)
    prot_mean <- rowMeans(prot_mat, na.rm = TRUE)

    # Add small value to avoid log(0)
    rna_mean <- rna_mean + 1
    prot_mean <- prot_mean + 1

    # Translation efficiency (mean across samples)
    te_scores <- log2(prot_mean) - log2(rna_mean)
    names(te_scores) <- gene_symbols

    # Per-sample TE
    te_matrix <- log2(prot_mat + 1) - log2(rna_mat + 1)

    # TE variability across samples
    te_cv <- apply(te_matrix, 1, function(x) {
        sd(x, na.rm = TRUE) / abs(mean(x, na.rm = TRUE))
    })

    te_df <- data.frame(
        gene_symbol = gene_symbols,
        translation_efficiency = te_scores,
        te_cv = te_cv,
        rna_mean = log2(rna_mean),
        protein_mean = log2(prot_mean),
        stringsAsFactors = FALSE
    )

    # Classify genes by TE
    te_df$te_class <- cut(te_df$translation_efficiency,
                          breaks = c(-Inf, -1, 1, Inf),
                          labels = c("low_TE", "moderate_TE", "high_TE"))

    # Summary
    summary_stats <- list(
        mean_te = mean(te_scores, na.rm = TRUE),
        sd_te = sd(te_scores, na.rm = TRUE),
        n_high_te = sum(te_df$te_class == "high_TE", na.rm = TRUE),
        n_low_te = sum(te_df$te_class == "low_TE", na.rm = TRUE)
    )

    message("  Translation efficiency: mean=", round(summary_stats$mean_te, 3),
            ", ", summary_stats$n_high_te, " high TE, ",
            summary_stats$n_low_te, " low TE")

    list(
        te_scores = te_scores,
        te_matrix = te_matrix,
        te_df = te_df,
        summary = summary_stats
    )
}

#' Detect post-transcriptional regulation
detect_post_transcriptional_regulation <- function(rna_mat, prot_mat, gene_symbols, mc) {
    n_genes <- length(gene_symbols)

    # RNA variability
    rna_var <- apply(rna_mat, 1, var, na.rm = TRUE)
    rna_cv <- apply(rna_mat, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))

    # Protein variability
    prot_var <- apply(prot_mat, 1, var, na.rm = TRUE)
    prot_cv <- apply(prot_mat, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))

    # RNA-protein correlation
    cors <- sapply(seq_len(n_genes), function(i) {
        cor(rna_mat[i, ], prot_mat[i, ], use = "pairwise.complete.obs")
    })

    # Residual protein variance not explained by RNA
    residual_var <- sapply(seq_len(n_genes), function(i) {
        fit <- tryCatch({
            lm(prot_mat[i, ] ~ rna_mat[i, ])
        }, error = function(e) NULL)

        if (!is.null(fit)) {
            var(residuals(fit))
        } else {
            NA
        }
    })

    # Fraction of protein variance explained by RNA
    r_squared <- cors^2
    unexplained_var <- 1 - r_squared

    # Post-transcriptional regulation score
    ptr_score <- unexplained_var * prot_cv / (prot_cv + 0.1)

    ptr_df <- data.frame(
        gene_symbol = gene_symbols,
        rna_cv = rna_cv,
        protein_cv = prot_cv,
        rna_protein_cor = cors,
        r_squared = r_squared,
        unexplained_variance = unexplained_var,
        residual_variance = residual_var,
        ptr_score = ptr_score,
        stringsAsFactors = FALSE
    )

    # Identify likely PTR-regulated genes
    ptr_df$likely_ptr <- ptr_df$rna_protein_cor < 0.3 &
        ptr_df$protein_cv > median(ptr_df$protein_cv, na.rm = TRUE)

    n_ptr <- sum(ptr_df$likely_ptr, na.rm = TRUE)
    message("  Post-transcriptional regulation: ", n_ptr,
            " genes with likely PTR")

    list(
        ptr_df = ptr_df,
        n_ptr = n_ptr,
        ptr_genes = gene_symbols[ptr_df$likely_ptr & !is.na(ptr_df$likely_ptr)]
    )
}

#' Infer protein stability
infer_protein_stability <- function(rna_mat, prot_mat, gene_symbols) {
    # Fit global RNA-protein relationship
    rna_mean <- rowMeans(rna_mat, na.rm = TRUE)
    prot_mean <- rowMeans(prot_mat, na.rm = TRUE)

    global_fit <- lm(log2(prot_mean + 1) ~ log2(rna_mean + 1))

    # Expected protein from global fit
    expected_prot <- predict(global_fit, newdata = data.frame(`log2(rna_mean + 1)` = log2(rna_mean + 1)))

    # Deviation from expected (stability score)
    stability_score <- log2(prot_mean + 1) - expected_prot

    stability_df <- data.frame(
        gene_symbol = gene_symbols,
        rna_mean = rna_mean,
        protein_mean = prot_mean,
        expected_protein = 2^expected_prot - 1,
        stability_score = stability_score,
        stringsAsFactors = FALSE
    )

    # Classify
    stability_df$stability_class <- cut(
        stability_df$stability_score,
        breaks = c(-Inf, -0.5, 0.5, Inf),
        labels = c("unstable", "normal", "stable")
    )

    # Summary
    summary_stats <- list(
        n_stable = sum(stability_df$stability_class == "stable", na.rm = TRUE),
        n_unstable = sum(stability_df$stability_class == "unstable", na.rm = TRUE),
        mean_stability = mean(stability_score, na.rm = TRUE),
        global_r2 = summary(global_fit)$r.squared
    )

    message("  Protein stability: ", summary_stats$n_stable, " stable, ",
            summary_stats$n_unstable, " unstable (global R2=",
            round(summary_stats$global_r2, 3), ")")

    list(
        stability_df = stability_df,
        global_fit = global_fit,
        summary = summary_stats
    )
}

#' Save RNA-protein regulation results
save_rna_protein_results <- function(results, out_dir) {
    if (!is.null(results$correlations)) {
        write.csv(results$correlations$cor_df,
                  file.path(out_dir, "tables", "mech_rna_protein_correlations.csv"),
                  row.names = FALSE)
    }

    if (!is.null(results$translation_efficiency)) {
        write.csv(results$translation_efficiency$te_df,
                  file.path(out_dir, "tables", "mech_translation_efficiency.csv"),
                  row.names = FALSE)
    }

    if (!is.null(results$post_transcriptional)) {
        write.csv(results$post_transcriptional$ptr_df,
                  file.path(out_dir, "tables", "mech_post_transcriptional_regulation.csv"),
                  row.names = FALSE)
    }

    if (!is.null(results$protein_stability)) {
        write.csv(results$protein_stability$stability_df,
                  file.path(out_dir, "tables", "mech_protein_stability.csv"),
                  row.names = FALSE)
    }
}

#' Plot RNA-protein regulation
plot_rna_protein_regulation <- function(results, out_dir) {
    # Plot 1: RNA-Protein correlation distribution
    if (!is.null(results$correlations)) {
        cor_df <- results$correlations$cor_df

        p1 <- ggplot2::ggplot(cor_df, ggplot2::aes(x = rna_protein_cor)) +
            ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
            ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
            ggplot2::geom_vline(xintercept = median(cor_df$rna_protein_cor, na.rm = TRUE),
                                linetype = "dashed", color = "darkgreen") +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = "RNA-Protein Correlation Distribution",
                subtitle = paste0("Median = ", round(median(cor_df$rna_protein_cor, na.rm = TRUE), 3)),
                x = "Spearman Correlation",
                y = "Count"
            )

        ggplot2::ggsave(file.path(out_dir, "plots", "mech_rna_protein_correlation_dist.png"),
                        p1, width = 8, height = 5)
    }

    # Plot 2: Translation efficiency
    if (!is.null(results$translation_efficiency)) {
        te_df <- results$translation_efficiency$te_df

        te_r <- cor(te_df$rna_mean, te_df$protein_mean, use = "pairwise.complete.obs")

        p2 <- ggplot2::ggplot(te_df, ggplot2::aes(x = rna_mean, y = protein_mean)) +
            ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
            ggplot2::geom_point(alpha = 0.6, size = 1.5, color = "gray40") +
            ggplot2::geom_smooth(method = "lm", se = TRUE, color = "black",
                                 linetype = "dotted", linewidth = 0.5) +
            ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = "RNA vs Protein Expression (Translation Efficiency)",
                subtitle = paste0("Pearson r = ", round(te_r, 3), ", n = ", nrow(te_df), " genes"),
                x = "log2(RNA)",
                y = "log2(Protein)"
            )

        ggplot2::ggsave(file.path(out_dir, "plots", "mech_translation_efficiency_scatter.png"),
                        p2, width = 8, height = 6)

        # TE distribution
        p2b <- ggplot2::ggplot(te_df, ggplot2::aes(x = translation_efficiency)) +
            ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white") +
            ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = "Translation Efficiency Distribution",
                x = "Translation Efficiency (log2 Protein - log2 RNA)",
                y = "Count"
            )

        ggplot2::ggsave(file.path(out_dir, "plots", "mech_translation_efficiency_dist.png"),
                        p2b, width = 8, height = 5)
    }

    # Plot 3: Protein stability
    if (!is.null(results$protein_stability)) {
        stab_df <- results$protein_stability$stability_df

        p3 <- ggplot2::ggplot(stab_df, ggplot2::aes(x = stability_score, fill = stability_class)) +
            ggplot2::geom_histogram(bins = 50, color = "white") +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = "Protein Stability Score Distribution",
                subtitle = "Deviation from expected protein level given RNA",
                x = "Stability Score",
                y = "Count",
                fill = "Stability"
            ) +
            ggplot2::scale_fill_manual(values = c("unstable" = "#d73027",
                                                   "normal" = "#fee090",
                                                   "stable" = "#4575b4"))

        ggplot2::ggsave(file.path(out_dir, "plots", "mech_protein_stability_dist.png"),
                        p3, width = 8, height = 5)
    }
}


# =============================================================================
# 2. Transcription Factor Activity Inference
# =============================================================================

#' Infer transcription factor activity
#' @param transcriptomics Transcriptomics data
#' @param mc Mechanistic config
#' @param config Full config
#' @param out_dir Output directory
infer_tf_activity <- function(transcriptomics, mc, config, out_dir = NULL) {
    message("Inferring transcription factor activity...")

    # Check organism support — dorothea regulons only cover human and mouse
    organism <- config$global$organism %||% "human"
    supported_organisms <- c("human", "homo_sapiens", "mouse", "mus_musculus")
    if (!tolower(organism) %in% supported_organisms) {
        message("  TF activity inference not available for '", organism,
                "' (dorothea regulons only cover human and mouse). Skipping.")
        return(NULL)
    }

    # Check if decoupleR or dorothea is available
    has_decoupler <- requireNamespace("decoupleR", quietly = TRUE)
    has_dorothea <- requireNamespace("dorothea", quietly = TRUE)

    if (!has_decoupler && !has_dorothea) {
        message("  Neither decoupleR nor dorothea available. Skipping TF activity.")
        return(NULL)
    }

    expr_mat <- transcriptomics$normalized_matrix

    # Get TF-target regulons
    regulons <- get_tf_regulons(mc$tf_database)
    if (is.null(regulons)) {
        message("  Could not load TF regulons. Skipping.")
        return(NULL)
    }

    message("  Using ", length(unique(regulons$tf)), " TFs with ",
            nrow(regulons), " TF-target pairs")

    # Infer TF activity
    if (has_decoupler) {
        tf_activity <- infer_tf_activity_decoupler(expr_mat, regulons)
    } else {
        tf_activity <- infer_tf_activity_dorothea(expr_mat, regulons)
    }

    if (is.null(tf_activity)) return(NULL)

    # Test TF activity associations with conditions
    if (!is.null(transcriptomics$metadata)) {
        tf_activity$associations <- test_tf_associations(
            tf_activity$activity_matrix,
            transcriptomics$metadata,
            config
        )
    }

    # Save results
    if (!is.null(out_dir)) {
        save_tf_activity_results(tf_activity, out_dir)
        plot_tf_activity(tf_activity, out_dir)
    }

    return(tf_activity)
}

#' Get TF-target regulons
get_tf_regulons <- function(database) {
    if (database == "dorothea") {
        if (!requireNamespace("dorothea", quietly = TRUE)) return(NULL)

        # Get high-confidence regulons (A and B)
        regulons <- tryCatch({
            data("dorothea_hs", package = "dorothea", envir = environment())
            dorothea_hs <- get("dorothea_hs", envir = environment())
            dorothea_hs[dorothea_hs$confidence %in% c("A", "B"), ]
        }, error = function(e) NULL)

        if (!is.null(regulons)) {
            regulons <- data.frame(
                tf = regulons$tf,
                target = regulons$target,
                mor = regulons$mor,
                confidence = regulons$confidence,
                stringsAsFactors = FALSE
            )
        }

        return(regulons)
    }

    return(NULL)
}

#' Infer TF activity using decoupleR
infer_tf_activity_decoupler <- function(expr_mat, regulons) {
    message("  Using decoupleR for TF activity inference...")

    activity <- tryCatch({
        decoupleR::run_ulm(
            mat = expr_mat,
            net = regulons,
            .source = "tf",
            .target = "target",
            .mor = "mor",
            minsize = 5
        )
    }, error = function(e) {
        message("  decoupleR::run_ulm failed: ", e$message)
        NULL
    })

    if (is.null(activity)) return(NULL)

    # Pivot to wide format
    activity_mat <- tryCatch({
        tidyr::pivot_wider(activity,
                           id_cols = "source",
                           names_from = "condition",
                           values_from = "score")
    }, error = function(e) NULL)

    if (is.null(activity_mat)) return(NULL)

    # Convert to matrix
    tf_names <- activity_mat$source
    activity_mat <- as.matrix(activity_mat[, -1])
    rownames(activity_mat) <- tf_names

    list(
        activity_matrix = activity_mat,
        activity_long = activity,
        method = "decoupleR_ulm"
    )
}

#' Infer TF activity using dorothea directly
infer_tf_activity_dorothea <- function(expr_mat, regulons) {
    message("  Using dorothea for TF activity inference...")

    tfs <- unique(regulons$tf)
    n_tfs <- length(tfs)
    n_samples <- ncol(expr_mat)

    activity_mat <- matrix(NA, n_tfs, n_samples)
    rownames(activity_mat) <- tfs
    colnames(activity_mat) <- colnames(expr_mat)

    for (i in seq_along(tfs)) {
        tf <- tfs[i]
        targets <- regulons[regulons$tf == tf, ]

        # Get target expression
        target_genes <- intersect(targets$target, rownames(expr_mat))
        if (length(target_genes) < 5) next

        target_expr <- expr_mat[target_genes, , drop = FALSE]
        target_mor <- targets$mor[match(target_genes, targets$target)]

        # Weight by mode of regulation
        weighted_expr <- sweep(target_expr, 1, target_mor, `*`)

        # TF activity = mean of weighted target expression
        activity_mat[i, ] <- colMeans(weighted_expr, na.rm = TRUE)
    }

    # Scale activity scores
    activity_mat <- t(scale(t(activity_mat)))

    # Remove TFs with no valid activity
    valid_tfs <- !apply(is.na(activity_mat), 1, all)
    activity_mat <- activity_mat[valid_tfs, , drop = FALSE]

    list(
        activity_matrix = activity_mat,
        method = "dorothea_weighted_mean"
    )
}

#' Test TF activity associations with conditions
test_tf_associations <- function(activity_mat, metadata, config) {
    condition_col <- config$design$condition_column %||% "condition"

    if (!condition_col %in% colnames(metadata)) return(NULL)

    # Get conditions for samples
    samples <- intersect(colnames(activity_mat), rownames(metadata))
    if (length(samples) < 5) return(NULL)

    activity_mat <- activity_mat[, samples, drop = FALSE]
    conditions <- metadata[samples, condition_col]

    # Test each TF
    tfs <- rownames(activity_mat)

    results <- lapply(tfs, function(tf) {
        activity <- activity_mat[tf, ]

        # ANOVA or t-test
        if (length(unique(conditions)) == 2) {
            test <- tryCatch({
                t.test(activity ~ conditions)
            }, error = function(e) NULL)

            if (!is.null(test)) {
                return(data.frame(
                    tf = tf,
                    test = "t.test",
                    statistic = test$statistic,
                    pvalue = test$p.value,
                    stringsAsFactors = FALSE
                ))
            }
        } else {
            test <- tryCatch({
                summary(aov(activity ~ conditions))[[1]]
            }, error = function(e) NULL)

            if (!is.null(test)) {
                return(data.frame(
                    tf = tf,
                    test = "anova",
                    statistic = test[1, "F value"],
                    pvalue = test[1, "Pr(>F)"],
                    stringsAsFactors = FALSE
                ))
            }
        }

        NULL
    })

    results_df <- do.call(rbind, results)
    if (is.null(results_df)) return(NULL)

    results_df$padj <- p.adjust(results_df$pvalue, method = "BH")
    results_df <- results_df[order(results_df$pvalue), ]

    n_sig <- sum(results_df$padj < 0.05, na.rm = TRUE)
    message("  TF associations: ", n_sig, " significantly associated with condition")

    results_df
}

#' Save TF activity results
save_tf_activity_results <- function(tf_activity, out_dir) {
    # Save activity matrix
    activity_df <- as.data.frame(tf_activity$activity_matrix)
    activity_df$tf <- rownames(activity_df)
    write.csv(activity_df,
              file.path(out_dir, "tables", "mech_tf_activity_scores.csv"),
              row.names = FALSE)

    # Save associations if available
    if (!is.null(tf_activity$associations)) {
        write.csv(tf_activity$associations,
                  file.path(out_dir, "tables", "mech_tf_condition_associations.csv"),
                  row.names = FALSE)
    }
}

#' Plot TF activity
plot_tf_activity <- function(tf_activity, out_dir) {
    activity_mat <- tf_activity$activity_matrix

    # Heatmap of top variable TFs
    tf_var <- apply(activity_mat, 1, var, na.rm = TRUE)
    top_tfs <- names(sort(tf_var, decreasing = TRUE))[1:min(50, length(tf_var))]

    mat_plot <- activity_mat[top_tfs, , drop = FALSE]

    # Use pheatmap if available
    if (requireNamespace("pheatmap", quietly = TRUE)) {
        grDevices::png(file.path(out_dir, "plots", "mech_tf_activity_heatmap.png"),
                       width = 12, height = 10, units = "in", res = 300)

        tryCatch({
            pheatmap::pheatmap(
                mat_plot,
                scale = "row",
                clustering_method = "ward.D2",
                show_colnames = ncol(mat_plot) < 50,
                main = "Transcription Factor Activity",
                color = grDevices::colorRampPalette(c("blue", "white", "red"))(100),
                border_color = NA
            )
        }, error = function(e) {
            message("  TF heatmap failed: ", e$message)
        })

        grDevices::dev.off()
    }

    # Bar plot of most significant TFs
    if (!is.null(tf_activity$associations)) {
        assoc <- tf_activity$associations
        assoc_sig <- assoc[assoc$padj < 0.1, ]

        if (nrow(assoc_sig) > 0) {
            assoc_sig <- utils::head(assoc_sig[order(assoc_sig$pvalue), ], 20)

            p <- ggplot2::ggplot(assoc_sig, ggplot2::aes(x = reorder(tf, -log10(pvalue)),
                                                          y = -log10(pvalue))) +
                ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
                ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
                ggplot2::coord_flip() +
                ggplot2::theme_minimal() +
                ggplot2::labs(
                    title = "TF Activity Association with Condition",
                    x = "Transcription Factor",
                    y = "-log10(p-value)"
                )

            ggplot2::ggsave(file.path(out_dir, "plots", "mech_tf_association_barplot.png"),
                            p, width = 8, height = 6)
        }
    }
}


# =============================================================================
# 3. Gene Regulatory Network Inference
# =============================================================================

#' Infer gene regulatory network
#' @param harmonized List of harmonized omics data
#' @param gene_mapping Gene mapping table
#' @param mc Mechanistic config
#' @param config Full config
#' @param out_dir Output directory
infer_regulatory_network <- function(harmonized, gene_mapping, mc, config, out_dir = NULL) {
    message("Inferring regulatory / co-expression network...")

    # Use transcriptomics if available, fall back to proteomics
    if ("transcriptomics" %in% names(harmonized)) {
        expr_mat <- harmonized$transcriptomics$normalized_matrix
        network_type <- "gene"
        message("  Using transcriptomics for gene regulatory network")
    } else if ("proteomics" %in% names(harmonized)) {
        expr_mat <- harmonized$proteomics$normalized_matrix
        network_type <- "protein"
        message("  Using proteomics for protein co-expression network")
    } else {
        message("  No transcriptomics or proteomics data. Skipping.")
        return(NULL)
    }

    # Select top variable genes
    gene_vars <- apply(expr_mat, 1, var, na.rm = TRUE)
    top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(500, nrow(expr_mat))]
    expr_sub <- expr_mat[top_genes, ]

    message("  Using ", length(top_genes), " genes for network inference")

    # Run network inference
    if (mc$network_method == "genie3") {
        network <- infer_network_genie3(expr_sub, mc)
    } else if (mc$network_method == "aracne") {
        network <- infer_network_aracne(expr_sub, mc)
    } else {
        network <- infer_network_correlation(expr_sub, mc)
    }

    if (is.null(network)) return(NULL)

    network$network_type <- network_type

    # Identify hub regulators / hub proteins
    network$hubs <- identify_hub_regulators(network$edges, mc)

    # Save results
    if (!is.null(out_dir)) {
        write.csv(network$edges,
                  file.path(out_dir, "tables", "mech_regulatory_network_edges.csv"),
                  row.names = FALSE)
        write.csv(network$hubs,
                  file.path(out_dir, "tables", "mech_hub_regulators.csv"),
                  row.names = FALSE)

        # Build gene symbol map from feature annotation
        gene_symbol_map <- NULL
        omics_used <- if ("transcriptomics" %in% names(harmonized)) "transcriptomics" else "proteomics"
        anno <- harmonized[[omics_used]]$feature_annotation
        if (!is.null(anno)) {
            feat_ids <- rownames(anno)
            for (col in c("gene_symbol", "Genes", "Gene.Symbol", "symbol")) {
                if (col %in% colnames(anno)) {
                    syms <- as.character(anno[[col]])
                    valid <- !is.na(syms) & syms != "" & syms != feat_ids
                    if (sum(valid) > 0) {
                        gene_symbol_map <- setNames(syms, feat_ids)
                        break
                    }
                }
            }
        }

        # Fallback: try DE table gene symbol map (proteomics UniProt -> Genes)
        if (is.null(gene_symbol_map) || sum(!is.na(gene_symbol_map)) == 0) {
            de_map <- tryCatch(
                build_gene_symbol_map_from_de(config),
                error = function(e) NULL
            )
            if (!is.null(de_map) && length(de_map) > 0) {
                gene_symbol_map <- de_map
            }
        }

        # Build omics source map: feature_id -> omics layer name
        omics_source_map <- setNames(
            rep(omics_used, length(top_genes)),
            top_genes
        )

        # Visualizations
        plot_regulatory_network(network, out_dir,
                                gene_symbol_map = gene_symbol_map,
                                omics_source_map = omics_source_map)
    }

    return(network)
}

#' Infer network using GENIE3
infer_network_genie3 <- function(expr_mat, mc) {
    if (!requireNamespace("GENIE3", quietly = TRUE)) {
        message("  GENIE3 not available. Using correlation-based network.")
        return(infer_network_correlation(expr_mat, mc))
    }

    message("  Running GENIE3...")

    weight_mat <- tryCatch({
        GENIE3::GENIE3(
            exprMatrix = expr_mat,
            nCores = 1,
            verbose = FALSE
        )
    }, error = function(e) {
        message("  GENIE3 failed: ", e$message)
        NULL
    })

    if (is.null(weight_mat)) {
        return(infer_network_correlation(expr_mat, mc))
    }

    # Convert to edge list
    link_list <- GENIE3::getLinkList(weight_mat)

    # Filter to top edges
    n_edges <- min(5000, nrow(link_list))
    edges <- link_list[1:n_edges, ]
    colnames(edges) <- c("regulator", "target", "weight")

    message("  GENIE3 network: ", nrow(edges), " edges")

    list(
        edges = edges,
        method = "GENIE3"
    )
}

#' Infer network using ARACNe
infer_network_aracne <- function(expr_mat, mc) {
    if (!requireNamespace("minet", quietly = TRUE)) {
        message("  minet not available. Using correlation-based network.")
        return(infer_network_correlation(expr_mat, mc))
    }

    message("  Running ARACNe (via minet)...")

    # Compute mutual information
    mi_mat <- tryCatch({
        minet::build.mim(t(expr_mat), estimator = "spearman")
    }, error = function(e) {
        message("  MI estimation failed: ", e$message)
        NULL
    })

    if (is.null(mi_mat)) {
        return(infer_network_correlation(expr_mat, mc))
    }

    # ARACNe network
    aracne_mat <- tryCatch({
        minet::aracne(mi_mat)
    }, error = function(e) NULL)

    if (is.null(aracne_mat)) {
        return(infer_network_correlation(expr_mat, mc))
    }

    # Convert to edge list
    edges <- mat_to_edge_list(aracne_mat, "MI")

    message("  ARACNe network: ", nrow(edges), " edges")

    list(
        edges = edges,
        method = "ARACNe"
    )
}

#' Infer network using correlation
infer_network_correlation <- function(expr_mat, mc) {
    message("  Using correlation-based network...")

    # Compute correlation matrix
    cor_mat <- cor(t(expr_mat), use = "pairwise.complete.obs", method = "spearman")

    # Convert to edge list (filter by threshold)
    edges <- mat_to_edge_list(cor_mat, "correlation", threshold = 0.5)

    message("  Correlation network: ", nrow(edges), " edges")

    list(
        edges = edges,
        method = "correlation"
    )
}

#' Convert matrix to edge list
mat_to_edge_list <- function(mat, weight_name, threshold = 0) {
    genes <- rownames(mat)
    n <- length(genes)

    edges <- data.frame()

    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            w <- mat[i, j]
            if (!is.na(w) && abs(w) > threshold) {
                edges <- rbind(edges, data.frame(
                    regulator = genes[i],
                    target = genes[j],
                    weight = w,
                    stringsAsFactors = FALSE
                ))
            }
        }
    }

    colnames(edges)[3] <- weight_name
    edges <- edges[order(abs(edges[[3]]), decreasing = TRUE), ]

    edges
}

#' Identify hub regulators
identify_hub_regulators <- function(edges, mc) {
    # Count outgoing edges per regulator
    reg_counts <- table(edges$regulator)
    reg_counts <- sort(reg_counts, decreasing = TRUE)

    hub_df <- data.frame(
        gene = names(reg_counts),
        n_targets = as.numeric(reg_counts),
        stringsAsFactors = FALSE
    )

    # Mark top regulators
    hub_df$is_hub <- hub_df$n_targets >= quantile(hub_df$n_targets, 0.9)

    message("  Identified ", sum(hub_df$is_hub), " hub regulators")

    hub_df
}

#' Plot regulatory network
#' @param network Network object with edges, hubs, method
#' @param out_dir Output directory
#' @param gene_symbol_map Named character vector: feature_id -> gene_symbol (optional)
#' @param omics_source_map Named character vector: feature_id -> omics layer (optional)
plot_regulatory_network <- function(network, out_dir, gene_symbol_map = NULL,
                                    omics_source_map = NULL) {
    edges <- network$edges

    # Use igraph for visualization
    if (!requireNamespace("igraph", quietly = TRUE)) return(NULL)

    # Create graph from top edges
    top_edges <- utils::head(edges, 500)

    g <- igraph::graph_from_data_frame(
        top_edges[, c("regulator", "target")],
        directed = TRUE
    )

    # Identify hubs and compute degree
    if (!is.null(network$hubs)) {
        hub_genes <- network$hubs$gene[network$hubs$is_hub]
        igraph::V(g)$is_hub <- igraph::V(g)$name %in% hub_genes
    } else {
        igraph::V(g)$is_hub <- FALSE
    }
    node_deg <- igraph::degree(g, mode = "all")
    is_hub <- igraph::V(g)$is_hub

    # Resolve gene symbols for nodes
    node_names <- igraph::V(g)$name
    gene_symbols <- if (!is.null(gene_symbol_map)) {
        mapped <- gene_symbol_map[node_names]
        ifelse(!is.na(mapped) & mapped != "" & mapped != node_names,
               mapped, NA_character_)
    } else {
        rep(NA_character_, length(node_names))
    }
    has_gene_symbols <- sum(!is.na(gene_symbols)) > length(gene_symbols) * 0.1

    # Compute layout once so both static plots use identical node positions
    set.seed(42)
    net_layout <- igraph::layout_with_fr(g)

    # Helper: generate static network PNG
    .plot_static_network <- function(out_path, use_gene_symbols = FALSE, title_suffix = "") {
        grDevices::png(out_path, width = 12, height = 12, units = "in", res = 300)
        tryCatch({
            v_colors <- ifelse(is_hub, "#d73027", "#4575b4")
            deg_scaled <- 3 + (node_deg / max(node_deg, 1)) * 10

            # Choose labels
            if (use_gene_symbols) {
                display <- ifelse(!is.na(gene_symbols), gene_symbols, node_names)
            } else {
                display <- node_names
            }

            max_labels <- 60
            n_hub_labels <- sum(is_hub)
            non_hub_order <- order(node_deg[!is_hub], decreasing = TRUE)
            n_extra <- min(max_labels - n_hub_labels, length(non_hub_order))
            top_non_hub_names <- node_names[!is_hub][non_hub_order[seq_len(n_extra)]]
            show_label <- is_hub | (node_names %in% top_non_hub_names)
            v_labels <- ifelse(show_label, display, NA)
            v_label_cex <- ifelse(is_hub, 0.8, 0.55)

            igraph::plot.igraph(
                g,
                vertex.size = deg_scaled,
                vertex.color = v_colors,
                vertex.label = v_labels,
                vertex.label.cex = v_label_cex,
                vertex.label.dist = 0.5,
                vertex.label.color = "black",
                edge.arrow.size = 0.3,
                edge.color = "gray70",
                layout = net_layout,
                main = paste0("Regulatory Network (", network$method, ")", title_suffix)
            )
            legend("bottomleft",
                   legend = c("Hub node (high degree)", "Non-hub node",
                              "Larger node = more connections"),
                   col = c("#d73027", "#4575b4", NA),
                   pt.bg = c("#d73027", "#4575b4", NA),
                   pch = c(21, 21, NA),
                   pt.cex = c(2.5, 1.5, NA),
                   cex = 0.9, bty = "n", title = "Legend")
        }, error = function(e) {
            message("  Network plot failed: ", e$message)
        })
        grDevices::dev.off()
    }

    # Static PNG with feature IDs
    .plot_static_network(file.path(out_dir, "plots", "mech_regulatory_network.png"))

    # Static PNG with gene symbols (only if symbols are available)
    if (has_gene_symbols) {
        .plot_static_network(
            file.path(out_dir, "plots", "mech_regulatory_network_genesymbol.png"),
            use_gene_symbols = TRUE,
            title_suffix = " - Gene Symbols"
        )
        message("  Gene symbol network plot saved")
    }

    # --- Interactive visNetwork widget ---
    tryCatch({
        if (!requireNamespace("visNetwork", quietly = TRUE)) {
            message("  visNetwork not available — skipping interactive network")
            return(invisible(NULL))
        }

        # Build tooltip with gene symbol and omics source if available
        tooltip_lines <- paste0("<b>", node_names, "</b>")
        if (has_gene_symbols) {
            tooltip_lines <- paste0(
                tooltip_lines, "<br>Gene Symbol: ",
                ifelse(!is.na(gene_symbols), gene_symbols, "N/A")
            )
        }
        # Add omics source information
        if (!is.null(omics_source_map)) {
            omics_labels <- omics_source_map[node_names]
            omics_labels[is.na(omics_labels)] <- network$network_type %||% "unknown"
            tooltip_lines <- paste0(tooltip_lines, "<br>Omics: ", omics_labels)
        } else if (!is.null(network$network_type)) {
            tooltip_lines <- paste0(tooltip_lines, "<br>Omics: ", network$network_type)
        }
        tooltip_lines <- paste0(
            tooltip_lines, "<br>Degree: ", node_deg, "<br>",
            ifelse(is_hub, "Hub regulator", "Non-hub")
        )

        # Use the same layout as the static plot (scaled for visNetwork coords)
        layout_scale <- 200
        vis_nodes <- data.frame(
            id = node_names,
            label = node_names,
            title = tooltip_lines,
            color = ifelse(is_hub, "#d73027", "#4575b4"),
            size = 10 + (node_deg / max(node_deg, 1)) * 30,
            font.size = ifelse(is_hub, 16, 11),
            borderWidth = ifelse(is_hub, 3, 1),
            x = net_layout[, 1] * layout_scale,
            y = -net_layout[, 2] * layout_scale,
            stringsAsFactors = FALSE
        )

        vis_edges <- data.frame(
            from = top_edges$regulator,
            to = top_edges$target,
            color = "rgba(150,150,150,0.4)",
            arrows = "to",
            width = 0.5,
            stringsAsFactors = FALSE
        )

        net <- visNetwork::visNetwork(
            vis_nodes, vis_edges,
            main = paste0("Regulatory Network (", network$method, ") - Interactive"),
            width = "100%", height = "700px"
        ) |>
            visNetwork::visPhysics(enabled = FALSE) |>
            visNetwork::visInteraction(
                hover = TRUE,
                tooltipDelay = 100,
                navigationButtons = TRUE,
                zoomView = TRUE
            ) |>
            visNetwork::visLegend(
                addNodes = list(
                    list(label = "Hub (high degree)", shape = "dot",
                         color = "#d73027", size = 20),
                    list(label = "Non-hub", shape = "dot",
                         color = "#4575b4", size = 10)
                ),
                useGroups = FALSE,
                position = "left",
                width = 0.15
            )

        # Save as RDS for embedding in report
        rds_path <- file.path(out_dir, "plots", "mech_regulatory_network_interactive.rds")
        saveRDS(net, rds_path)
        message("  Interactive network widget saved: ", rds_path)

    }, error = function(e) {
        message("  Interactive network creation failed: ", e$message)
    })
}


# =============================================================================
# 4. Causal Mediation Analysis
# =============================================================================

#' Run mediation analysis: RNA -> Protein -> Phenotype
#' @param harmonized List of harmonized omics data
#' @param metadata Sample metadata
#' @param gene_mapping Gene mapping table
#' @param mc Mechanistic config
#' @param config Full config
#' @param out_dir Output directory
run_mediation_analysis <- function(harmonized, metadata, gene_mapping, mc, config, out_dir = NULL) {
    message("Running mediation analysis (RNA -> Protein -> Phenotype)...")

    # Check requirements
    if (!"transcriptomics" %in% names(harmonized) ||
        !"proteomics" %in% names(harmonized)) {
        message("  Both transcriptomics and proteomics required. Skipping.")
        return(NULL)
    }

    if (!requireNamespace("mediation", quietly = TRUE)) {
        message("  mediation package not available. Skipping.")
        return(NULL)
    }

    # Get outcome variable
    outcome_col <- config$design$outcome_column %||% config$design$condition_column
    if (is.null(outcome_col) || !outcome_col %in% colnames(metadata)) {
        message("  No outcome variable specified. Skipping mediation.")
        return(NULL)
    }

    # Map features
    rna_mat <- harmonized$transcriptomics$normalized_matrix
    prot_mat <- harmonized$proteomics$normalized_matrix

    mapping_result <- map_rna_protein_features(rna_mat, prot_mat, gene_mapping)
    if (is.null(mapping_result)) return(NULL)

    # Get common samples with outcome
    common_samples <- Reduce(intersect, list(
        colnames(rna_mat),
        colnames(prot_mat),
        rownames(metadata)
    ))

    if (length(common_samples) < 20) {
        message("  Insufficient samples for mediation (need >= 20, have ",
                length(common_samples), "). Skipping.")
        return(NULL)
    }

    rna_matched <- mapping_result$rna[, common_samples, drop = FALSE]
    prot_matched <- mapping_result$protein[, common_samples, drop = FALSE]
    gene_symbols <- mapping_result$gene_symbols
    outcome <- metadata[common_samples, outcome_col]

    # Convert outcome to numeric if factor
    if (is.factor(outcome) || is.character(outcome)) {
        outcome <- as.numeric(as.factor(outcome)) - 1
    }

    message("  Testing ", length(gene_symbols), " genes for mediation")

    # Run mediation for each gene (limit to significant genes)
    rna_assoc <- apply(rna_matched, 1, function(x) {
        tryCatch({
            cor.test(x, outcome, method = "spearman")$p.value
        }, error = function(e) 1)
    })

    # Focus on genes associated with outcome
    sig_genes <- gene_symbols[rna_assoc < 0.1]
    sig_genes <- utils::head(sig_genes, 100)

    if (length(sig_genes) < 5) {
        message("  Too few genes associated with outcome. Skipping.")
        return(NULL)
    }

    message("  Running mediation for ", length(sig_genes), " candidate genes")

    # Run mediation
    mediation_results <- lapply(seq_along(sig_genes), function(i) {
        gene <- sig_genes[i]
        gene_idx <- which(gene_symbols == gene)

        rna_expr <- rna_matched[gene_idx, ]
        prot_expr <- prot_matched[gene_idx, ]

        run_single_mediation(rna_expr, prot_expr, outcome, gene)
    })

    # Combine results
    mediation_df <- do.call(rbind, mediation_results)
    mediation_df <- mediation_df[!is.na(mediation_df$prop_mediated), ]

    if (nrow(mediation_df) == 0) {
        message("  No successful mediation analyses.")
        return(NULL)
    }

    # Adjust p-values
    mediation_df$acme_padj <- p.adjust(mediation_df$acme_pval, method = "BH")

    # Sort by proportion mediated
    mediation_df <- mediation_df[order(mediation_df$prop_mediated, decreasing = TRUE), ]

    # Summary
    n_mediated <- sum(mediation_df$prop_mediated > 0.1 &
                        mediation_df$acme_padj < 0.1, na.rm = TRUE)
    message("  Found ", n_mediated, " genes with significant mediation")

    # Save results
    if (!is.null(out_dir)) {
        write.csv(mediation_df,
                  file.path(out_dir, "tables", "mech_mediation_analysis_results.csv"),
                  row.names = FALSE)
        plot_mediation_results(mediation_df, out_dir)
    }

    list(
        mediation_df = mediation_df,
        n_mediated = n_mediated
    )
}

#' Run single gene mediation analysis with sensitivity
run_single_mediation <- function(rna_expr, prot_expr, outcome, gene_name) {
    tryCatch({
        # Create data frame
        df <- data.frame(
            rna = as.numeric(rna_expr),
            protein = as.numeric(prot_expr),
            outcome = as.numeric(outcome)
        )
        df <- df[complete.cases(df), ]

        if (nrow(df) < 20) {
            return(data.frame(gene = gene_name, prop_mediated = NA, acme_pval = NA,
                              ade_pval = NA, rho_at_acme_zero = NA,
                              stringsAsFactors = FALSE))
        }

        # Mediator model: Protein ~ RNA
        med_model <- lm(protein ~ rna, data = df)

        # Outcome model: Outcome ~ RNA + Protein
        out_model <- lm(outcome ~ rna + protein, data = df)

        # Run mediation
        med_result <- mediation::mediate(
            med_model, out_model,
            treat = "rna",
            mediator = "protein",
            boot = FALSE,
            sims = 100
        )

        # Sensitivity analysis: find rho at which ACME = 0
        # medsens probes the sequential ignorability assumption
        # by varying the correlation between mediator and outcome residuals
        rho_at_zero <- NA_real_
        tryCatch({
            sens <- mediation::medsens(
                med_result,
                rho.by = 0.05,
                effect.type = "indirect",
                sims = 100
            )
            # Extract rho value where ACME crosses zero
            # sens$rho gives the rho grid, sens$d0 gives ACME at each rho
            if (!is.null(sens$rho) && !is.null(sens$d0)) {
                # Find sign changes in ACME across rho values
                signs <- sign(sens$d0)
                sign_changes <- which(diff(signs) != 0)
                if (length(sign_changes) > 0) {
                    # Interpolate to find rho where ACME = 0
                    idx <- sign_changes[1]
                    rho_at_zero <- sens$rho[idx] +
                        (sens$rho[idx + 1] - sens$rho[idx]) *
                        abs(sens$d0[idx]) / (abs(sens$d0[idx]) + abs(sens$d0[idx + 1]))
                } else {
                    # ACME never crosses zero — highly robust
                    rho_at_zero <- Inf
                }
            }
        }, error = function(e) {
            # Sensitivity analysis can fail with singular models
            NULL
        })

        data.frame(
            gene = gene_name,
            acme = med_result$d0,
            acme_pval = med_result$d0.p,
            ade = med_result$z0,
            ade_pval = med_result$z0.p,
            total_effect = med_result$tau.coef,
            prop_mediated = med_result$n0,
            rho_at_acme_zero = rho_at_zero,
            stringsAsFactors = FALSE
        )
    }, error = function(e) {
        data.frame(gene = gene_name, prop_mediated = NA, acme_pval = NA,
                   rho_at_acme_zero = NA, stringsAsFactors = FALSE)
    })
}

#' Plot mediation results
plot_mediation_results <- function(mediation_df, out_dir) {
    # Filter to significant
    med_sig <- mediation_df[mediation_df$acme_padj < 0.1 & !is.na(mediation_df$acme_padj), ]

    if (nrow(med_sig) == 0) {
        med_sig <- utils::head(mediation_df[!is.na(mediation_df$prop_mediated), ], 20)
    }

    if (nrow(med_sig) == 0) return(NULL)

    # Bar plot of proportion mediated
    med_sig <- utils::head(med_sig[order(med_sig$prop_mediated, decreasing = TRUE), ], 20)

    p1 <- ggplot2::ggplot(med_sig, ggplot2::aes(x = reorder(gene, prop_mediated),
                                                 y = prop_mediated)) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        ggplot2::geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
        ggplot2::coord_flip() +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "Mediation Analysis: RNA -> Protein -> Phenotype",
            subtitle = "Proportion of RNA effect mediated by protein",
            x = "Gene",
            y = "Proportion Mediated"
        ) +
        ggplot2::ylim(0, 1)

    ggplot2::ggsave(file.path(out_dir, "plots", "mech_mediation_proportion_barplot.png"),
                    p1, width = 8, height = 6)

    # Effect decomposition
    if (all(c("acme", "ade") %in% colnames(med_sig))) {
        med_long <- tidyr::pivot_longer(med_sig,
                                         cols = c("acme", "ade"),
                                         names_to = "effect_type",
                                         values_to = "effect_size")

        p2 <- ggplot2::ggplot(med_long, ggplot2::aes(x = reorder(gene, effect_size),
                                                       y = effect_size,
                                                       fill = effect_type)) +
            ggplot2::geom_bar(stat = "identity", position = "dodge") +
            ggplot2::coord_flip() +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = "Mediation Effect Decomposition",
                subtitle = "ACME = mediated effect, ADE = direct effect",
                x = "Gene",
                y = "Effect Size",
                fill = "Effect Type"
            ) +
            ggplot2::scale_fill_manual(values = c("acme" = "#4575b4", "ade" = "#d73027"),
                                        labels = c("acme" = "Mediated (ACME)",
                                                  "ade" = "Direct (ADE)"))

        ggplot2::ggsave(file.path(out_dir, "plots", "mech_mediation_effect_decomposition.png"),
                        p2, width = 8, height = 6)
    }
}


# =============================================================================
# 5. Multi-Layer Mediation (RNA -> Protein -> Metabolite -> Phenotype)
# =============================================================================

#' Run multi-layer mediation: RNA -> Protein -> Metabolite -> Phenotype
#'
#' Tests cascading mediation chains across three omics layers.
#' Each step is a separate mediation analysis:
#'   Step 1: RNA -> Protein -> Phenotype (already done in standard mediation)
#'   Step 2: Protein -> Metabolite -> Phenotype
#'
#' NOTE: Each additional step compounds the untestable sequential ignorability
#' assumption. Results should be treated as hypothesis-generating, not as
#' proof of causality.
#'
#' @param harmonized List of harmonized omics data
#' @param metadata Sample metadata
#' @param mc Mechanistic config
#' @param config Full config
#' @param out_dir Output directory
#' @return List with multilayer mediation results, or NULL
run_multilayer_mediation <- function(harmonized, metadata, mc, config,
                                     out_dir = NULL) {
    message("Running multi-layer mediation (RNA -> Protein -> Metabolite -> Phenotype)...")

    required <- c("transcriptomics", "proteomics", "metabolomics")
    present <- intersect(required, names(harmonized))
    if (length(present) < 3) {
        message("  Requires all three omics layers. Present: ",
                paste(present, collapse = ", "), ". Skipping.")
        return(NULL)
    }

    if (!requireNamespace("mediation", quietly = TRUE)) {
        message("  mediation package not available. Skipping.")
        return(NULL)
    }

    outcome_col <- config$design$outcome_column %||% config$design$condition_column
    if (is.null(outcome_col) || !outcome_col %in% colnames(metadata)) {
        message("  No outcome variable. Skipping multi-layer mediation.")
        return(NULL)
    }

    rna_mat <- harmonized$transcriptomics$normalized_matrix
    prot_mat <- harmonized$proteomics$normalized_matrix
    metab_mat <- harmonized$metabolomics$normalized_matrix

    # Find common samples across all three omics + metadata
    common_samples <- Reduce(intersect, list(
        colnames(rna_mat), colnames(prot_mat), colnames(metab_mat),
        rownames(metadata)
    ))

    if (length(common_samples) < 20) {
        message("  Insufficient common samples (", length(common_samples),
                "). Need >= 20. Skipping.")
        return(NULL)
    }

    rna_mat <- rna_mat[, common_samples, drop = FALSE]
    prot_mat <- prot_mat[, common_samples, drop = FALSE]
    metab_mat <- metab_mat[, common_samples, drop = FALSE]
    outcome <- metadata[common_samples, outcome_col]

    if (is.factor(outcome) || is.character(outcome)) {
        outcome <- as.numeric(as.factor(outcome)) - 1
    }

    # For each protein, find the most correlated metabolite
    # This is a data-driven pairing since we lack explicit protein-metabolite mapping
    n_prot <- min(nrow(prot_mat), 100)
    # Rank proteins by variance for efficiency
    prot_var <- apply(prot_mat, 1, var, na.rm = TRUE)
    top_prot_idx <- order(prot_var, decreasing = TRUE)[seq_len(n_prot)]

    message("  Testing ", n_prot, " proteins x top metabolite correlations")

    results_list <- list()
    for (i in seq_along(top_prot_idx)) {
        pi <- top_prot_idx[i]
        prot_id <- rownames(prot_mat)[pi]
        prot_vals <- as.numeric(prot_mat[pi, ])

        # Find best-correlated metabolite
        metab_cors <- apply(metab_mat, 1, function(m) {
            tryCatch(cor(prot_vals, m, use = "pairwise.complete.obs"),
                     error = function(e) 0)
        })
        best_metab_idx <- which.max(abs(metab_cors))
        if (abs(metab_cors[best_metab_idx]) < 0.2) next  # skip weak pairs

        metab_id <- rownames(metab_mat)[best_metab_idx]
        metab_vals <- as.numeric(metab_mat[best_metab_idx, ])

        # Mediation: Protein -> Metabolite -> Phenotype
        res <- tryCatch({
            df <- data.frame(
                protein = prot_vals,
                metabolite = metab_vals,
                outcome = outcome
            )
            df <- df[complete.cases(df), ]
            if (nrow(df) < 20) return(NULL)

            med_mod <- lm(metabolite ~ protein, data = df)
            out_mod <- lm(outcome ~ protein + metabolite, data = df)

            med_res <- mediation::mediate(
                med_mod, out_mod,
                treat = "protein", mediator = "metabolite",
                boot = FALSE, sims = 100
            )

            data.frame(
                protein = prot_id,
                metabolite = metab_id,
                prot_metab_cor = metab_cors[best_metab_idx],
                acme = med_res$d0,
                acme_pval = med_res$d0.p,
                ade = med_res$z0,
                ade_pval = med_res$z0.p,
                total_effect = med_res$tau.coef,
                prop_mediated = med_res$n0,
                stringsAsFactors = FALSE
            )
        }, error = function(e) NULL)

        if (!is.null(res)) results_list[[length(results_list) + 1]] <- res
    }

    if (length(results_list) == 0) {
        message("  No successful multi-layer mediation results.")
        return(NULL)
    }

    multilayer_df <- do.call(rbind, results_list)
    multilayer_df$acme_padj <- p.adjust(multilayer_df$acme_pval, method = "BH")
    multilayer_df <- multilayer_df[order(abs(multilayer_df$prop_mediated),
                                         decreasing = TRUE), ]

    n_sig <- sum(multilayer_df$acme_padj < 0.1 &
                     abs(multilayer_df$prop_mediated) > 0.1, na.rm = TRUE)
    message("  Found ", n_sig, " significant Protein->Metabolite mediation paths")

    if (!is.null(out_dir)) {
        write.csv(multilayer_df,
                  file.path(out_dir, "tables", "mech_multilayer_mediation_results.csv"),
                  row.names = FALSE)
        plot_multilayer_mediation(multilayer_df, out_dir)
    }

    list(multilayer_df = multilayer_df, n_significant = n_sig)
}


#' Plot multi-layer mediation results
plot_multilayer_mediation <- function(multilayer_df, out_dir) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)

    # Show top results
    top <- utils::head(multilayer_df[!is.na(multilayer_df$prop_mediated), ], 20)
    if (nrow(top) == 0) return(NULL)

    top$label <- paste0(top$protein, " -> ", top$metabolite)
    top$label <- ifelse(nchar(top$label) > 40,
                        paste0(substr(top$label, 1, 37), "..."), top$label)

    p <- ggplot2::ggplot(top, ggplot2::aes(
            x = reorder(label, abs(prop_mediated)),
            y = prop_mediated)) +
        ggplot2::geom_col(ggplot2::aes(fill = prop_mediated > 0)) +
        ggplot2::scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
                                    guide = "none") +
        ggplot2::geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed",
                            color = "grey50") +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid = ggplot2::element_blank()) +
        ggplot2::labs(
            title = "Multi-Layer Mediation: Protein -> Metabolite -> Phenotype",
            subtitle = "Proportion of protein effect mediated by metabolite (hypothesis-generating)",
            x = NULL,
            y = "Proportion Mediated"
        )

    ggplot2::ggsave(
        file.path(out_dir, "plots", "mech_multilayer_mediation_barplot.png"),
        p, width = 10, height = 7, dpi = 300
    )
}


# =============================================================================
# 5b. Protein -> Metabolite -> Phenotype Mediation (no RNA-seq required)
# =============================================================================

#' Run Protein -> Metabolite -> Phenotype mediation
#'
#' Self-contained mediation analysis that works with proteomics + metabolomics
#' only. For each protein, finds the best-correlated metabolite and tests
#' whether the metabolite mediates the protein's effect on the phenotype.
#'
#' @param harmonized Named list with proteomics and metabolomics data
#' @param metadata Sample metadata
#' @param mc Mechanistic config
#' @param config Full pipeline config
#' @param out_dir Output directory
#' @return List with mediation results, or NULL
run_protein_metabolite_mediation <- function(harmonized, metadata, mc, config,
                                              out_dir = NULL) {
    message("Running Protein -> Metabolite -> Phenotype mediation...")

    if (!requireNamespace("mediation", quietly = TRUE)) {
        message("  mediation package not available. Skipping.")
        return(NULL)
    }

    outcome_col <- config$design$outcome_column %||% config$design$condition_column
    if (is.null(outcome_col) || !outcome_col %in% colnames(metadata)) {
        message("  No outcome variable. Skipping.")
        return(NULL)
    }

    prot_mat <- harmonized$proteomics$normalized_matrix %||%
        harmonized$proteomics$expr_work
    metab_mat <- harmonized$metabolomics$normalized_matrix %||%
        harmonized$metabolomics$expr_work

    if (is.null(prot_mat) || is.null(metab_mat)) {
        message("  Missing expression matrices. Skipping.")
        return(NULL)
    }

    # Align samples
    common_samples <- Reduce(intersect, list(
        colnames(prot_mat), colnames(metab_mat), rownames(metadata)
    ))

    if (length(common_samples) < 6) {
        message("  Insufficient common samples (", length(common_samples),
                "). Skipping.")
        return(NULL)
    }

    prot_mat <- prot_mat[, common_samples, drop = FALSE]
    metab_mat <- metab_mat[, common_samples, drop = FALSE]
    outcome <- metadata[common_samples, outcome_col]

    if (is.factor(outcome) || is.character(outcome)) {
        outcome <- as.numeric(as.factor(outcome)) - 1
    }

    # Rank proteins by association with outcome
    prot_assoc_p <- apply(prot_mat, 1, function(x) {
        tryCatch(cor.test(x, outcome, method = "spearman")$p.value,
                 error = function(e) 1)
    })
    # Use top outcome-associated proteins (or all if few)
    candidate_prots <- names(sort(prot_assoc_p))[seq_len(min(100, nrow(prot_mat)))]

    message("  Testing ", length(candidate_prots), " proteins for mediation via metabolites")

    results_list <- list()
    for (prot_id in candidate_prots) {
        prot_vals <- as.numeric(prot_mat[prot_id, ])

        # Find best-correlated metabolite
        metab_cors <- apply(metab_mat, 1, function(m) {
            tryCatch(cor(prot_vals, m, use = "pairwise.complete.obs"),
                     error = function(e) 0)
        })
        best_idx <- which.max(abs(metab_cors))
        if (abs(metab_cors[best_idx]) < 0.3) next

        metab_id <- rownames(metab_mat)[best_idx]
        metab_vals <- as.numeric(metab_mat[best_idx, ])

        res <- tryCatch({
            df <- data.frame(
                protein = prot_vals,
                metabolite = metab_vals,
                outcome = outcome
            )
            df <- df[complete.cases(df), ]
            if (nrow(df) < 6) return(NULL)

            med_mod <- lm(metabolite ~ protein, data = df)
            out_mod <- lm(outcome ~ protein + metabolite, data = df)

            med_res <- mediation::mediate(
                med_mod, out_mod,
                treat = "protein", mediator = "metabolite",
                boot = FALSE, sims = 100
            )

            # Sensitivity analysis
            sens <- tryCatch({
                s <- mediation::medsens(med_res, rho.by = 0.05, sims = 100)
                s$rho.at.sign
            }, error = function(e) NA_real_)

            data.frame(
                protein = prot_id,
                metabolite = metab_id,
                prot_metab_cor = metab_cors[best_idx],
                prot_outcome_p = prot_assoc_p[prot_id],
                acme = med_res$d0,
                acme_pval = med_res$d0.p,
                ade = med_res$z0,
                ade_pval = med_res$z0.p,
                total_effect = med_res$tau.coef,
                prop_mediated = med_res$n0,
                rho_at_acme_zero = sens,
                stringsAsFactors = FALSE
            )
        }, error = function(e) NULL)

        if (!is.null(res)) results_list[[length(results_list) + 1]] <- res
    }

    if (length(results_list) == 0) {
        message("  No successful Protein->Metabolite mediation results.")
        return(NULL)
    }

    med_df <- do.call(rbind, results_list)
    med_df$acme_padj <- p.adjust(med_df$acme_pval, method = "BH")
    med_df <- med_df[order(abs(med_df$prop_mediated), decreasing = TRUE), ]

    n_sig <- sum(med_df$acme_padj < 0.1 &
                     abs(med_df$prop_mediated) > 0.1, na.rm = TRUE)
    message("  Found ", n_sig, " significant Protein->Metabolite->Phenotype paths")

    if (!is.null(out_dir)) {
        tables_dir <- file.path(out_dir, "tables")
        plots_dir <- file.path(out_dir, "plots")
        dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

        write.csv(med_df,
                  file.path(tables_dir, "mech_prot_metab_mediation_results.csv"),
                  row.names = FALSE)

        # Plot
        if (requireNamespace("ggplot2", quietly = TRUE) && nrow(med_df) > 0) {
            top <- utils::head(med_df[!is.na(med_df$prop_mediated), ], 20)
            if (nrow(top) > 0) {
                top$label <- paste0(top$protein, " -> ", top$metabolite)
                top$label <- ifelse(nchar(top$label) > 40,
                                    paste0(substr(top$label, 1, 37), "..."),
                                    top$label)

                p <- ggplot2::ggplot(top, ggplot2::aes(
                        x = reorder(label, abs(prop_mediated)),
                        y = prop_mediated)) +
                    ggplot2::geom_col(ggplot2::aes(fill = prop_mediated > 0)) +
                    ggplot2::scale_fill_manual(
                        values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
                        guide = "none") +
                    ggplot2::geom_hline(yintercept = c(-0.1, 0.1),
                                        linetype = "dashed", color = "grey50") +
                    ggplot2::coord_flip() +
                    ggplot2::theme_bw() +
                    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
                    ggplot2::labs(
                        title = "Protein -> Metabolite -> Phenotype Mediation",
                        subtitle = "Proportion of protein effect on phenotype mediated through metabolite",
                        x = NULL,
                        y = "Proportion Mediated"
                    )

                ggplot2::ggsave(
                    file.path(plots_dir, "mech_prot_metab_mediation_barplot.png"),
                    p, width = 10, height = 7, dpi = 300
                )
            }
        }
    }

    list(mediation_df = med_df, n_significant = n_sig)
}


# =============================================================================
# 6b. Pathway Activity Inference from Proteomics (PROGENy via decoupleR)
# =============================================================================

#' Infer pathway activity from proteomics data using PROGENy regulons
#'
#' Uses decoupleR's weighted mean method with PROGENy pathway-gene weights
#' derived from perturbation experiments. Works on protein-level expression.
#'
#' @param proteomics_data List with normalized_matrix slot
#' @param mc Mechanistic config
#' @param config Full pipeline config
#' @param out_dir Output directory
#' @return List with pathway activity matrix, associations, and plots
infer_pathway_activity_proteomics <- function(proteomics_data, mc, config,
                                               out_dir = NULL) {
    message("Inferring pathway activity from proteomics (PROGENy)...")

    # Check organism support — PROGENy only covers human and mouse
    organism <- config$global$organism %||% "human"
    supported_organisms <- c("human", "homo_sapiens", "mouse", "mus_musculus")
    if (!tolower(organism) %in% supported_organisms) {
        message("  PROGENy pathway activity not available for '", organism,
                "' (only human and mouse are supported). Skipping.")
        return(NULL)
    }

    if (!requireNamespace("decoupleR", quietly = TRUE)) {
        message("  decoupleR not available. Skipping pathway activity.")
        return(NULL)
    }

    expr_mat <- proteomics_data$normalized_matrix %||% proteomics_data$expr_work
    if (is.null(expr_mat) || nrow(expr_mat) < 10) {
        message("  Insufficient protein features. Skipping.")
        return(NULL)
    }

    # Get PROGENy pathway-gene weights
    progeny_net <- tryCatch({
        decoupleR::get_progeny(organism = "human", top = 500)
    }, error = function(e) {
        tryCatch({
            decoupleR::get_progeny(organism = "mouse", top = 500)
        }, error = function(e2) NULL)
    })

    if (is.null(progeny_net) || nrow(progeny_net) == 0) {
        message("  Could not load PROGENy regulons. Skipping.")
        return(NULL)
    }

    # Run decoupleR weighted mean
    activity <- tryCatch({
        decoupleR::run_wmean(
            mat = expr_mat,
            net = progeny_net,
            .source = "source",
            .target = "target",
            .mor = "weight",
            times = 100,
            minsize = 5
        )
    }, error = function(e) {
        message("  decoupleR::run_wmean failed: ", e$message)
        NULL
    })

    if (is.null(activity) || nrow(activity) == 0) {
        message("  No pathway activities inferred. Skipping.")
        return(NULL)
    }

    # Filter to consensus scores
    act_consensus <- activity[activity$statistic == "norm_wmean", ]
    if (nrow(act_consensus) == 0) {
        act_consensus <- activity[activity$statistic == "wmean", ]
    }

    # Pivot to matrix (pathways x samples)
    pathways <- unique(act_consensus$source)
    samples <- unique(act_consensus$condition)
    act_mat <- matrix(NA, nrow = length(pathways), ncol = length(samples),
                      dimnames = list(pathways, samples))
    for (i in seq_len(nrow(act_consensus))) {
        act_mat[act_consensus$source[i], act_consensus$condition[i]] <-
            act_consensus$score[i]
    }

    message("  Inferred activity for ", length(pathways), " pathways across ",
            length(samples), " samples")

    # Test pathway-condition associations
    condition_col <- config$design$condition_column
    metadata <- proteomics_data$metadata
    associations <- NULL
    if (!is.null(condition_col) && !is.null(metadata) &&
        condition_col %in% colnames(metadata)) {
        condition <- metadata[samples, condition_col]
        if (!is.null(condition) && length(unique(na.omit(condition))) >= 2) {
            assoc_list <- lapply(pathways, function(pw) {
                tryCatch({
                    score <- act_mat[pw, ]
                    grp <- as.factor(condition)
                    tt <- t.test(score ~ grp)
                    data.frame(
                        pathway = pw,
                        mean_diff = diff(tapply(score, grp, mean, na.rm = TRUE)),
                        pvalue = tt$p.value,
                        stringsAsFactors = FALSE
                    )
                }, error = function(e) NULL)
            })
            associations <- do.call(rbind, assoc_list[!sapply(assoc_list, is.null)])
            if (!is.null(associations) && nrow(associations) > 0) {
                associations$padj <- p.adjust(associations$pvalue, method = "BH")
                associations <- associations[order(associations$pvalue), ]
            }
        }
    }

    # Save and plot
    if (!is.null(out_dir)) {
        tables_dir <- file.path(out_dir, "tables")
        plots_dir <- file.path(out_dir, "plots")
        dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

        write.csv(as.data.frame(act_mat),
                  file.path(tables_dir, "mech_pathway_activity_matrix.csv"))

        if (!is.null(associations)) {
            write.csv(associations,
                      file.path(tables_dir, "mech_pathway_activity_associations.csv"),
                      row.names = FALSE)
        }

        # Heatmap
        if (requireNamespace("pheatmap", quietly = TRUE)) {
            tryCatch({
                clean_mat <- act_mat[!apply(is.na(act_mat), 1, all), , drop = FALSE]
                if (nrow(clean_mat) >= 2) {
                    ann_col <- NULL
                    if (!is.null(condition_col) && !is.null(metadata)) {
                        ann_col <- data.frame(
                            row.names = samples,
                            Condition = metadata[samples, condition_col]
                        )
                    }

                    grDevices::png(file.path(plots_dir,
                                  "mech_pathway_activity_heatmap.png"),
                                  width = 10, height = 8, units = "in", res = 300)
                    pheatmap::pheatmap(
                        clean_mat,
                        annotation_col = ann_col,
                        scale = "row",
                        cluster_rows = nrow(clean_mat) > 1,
                        cluster_cols = ncol(clean_mat) > 1,
                        main = "PROGENy Pathway Activity (from Proteomics)",
                        fontsize_row = 8,
                        fontsize_col = 8
                    )
                    grDevices::dev.off()
                }
            }, error = function(e) {
                message("  Pathway activity heatmap failed: ", e$message)
                try(grDevices::dev.off(), silent = TRUE)
            })
        }

        # Barplot of significant pathways
        if (!is.null(associations) && requireNamespace("ggplot2", quietly = TRUE)) {
            assoc_plot <- associations[!is.na(associations$pvalue), ]
            if (nrow(assoc_plot) > 0) {
                assoc_plot$sig <- ifelse(assoc_plot$padj < 0.05, "padj < 0.05",
                                  ifelse(assoc_plot$padj < 0.1, "padj < 0.1", "NS"))

                p <- ggplot2::ggplot(assoc_plot, ggplot2::aes(
                        x = reorder(pathway, -log10(pvalue)),
                        y = -log10(pvalue),
                        fill = sig)) +
                    ggplot2::geom_col() +
                    ggplot2::scale_fill_manual(values = c(
                        "padj < 0.05" = "#E64B35",
                        "padj < 0.1" = "#F39B7F",
                        "NS" = "grey70"
                    )) +
                    ggplot2::coord_flip() +
                    ggplot2::theme_bw() +
                    ggplot2::labs(
                        title = "Pathway Activity Differences (PROGENy from Proteomics)",
                        x = NULL,
                        y = "-log10(p-value)",
                        fill = "Significance"
                    )

                ggplot2::ggsave(
                    file.path(plots_dir, "mech_pathway_activity_barplot.png"),
                    p, width = 10, height = 6, dpi = 300
                )
            }
        }
    }

    list(
        activity_matrix = act_mat,
        associations = associations,
        n_pathways = length(pathways),
        method = "PROGENy_wmean"
    )
}


# =============================================================================
# Summary Function
# =============================================================================

#' Summarize mechanistic analysis results
summarize_mechanistic_results <- function(results, out_dir = NULL) {
    message("Generating mechanistic analysis summary...")

    summary_list <- list()

    # RNA-Protein regulation
    if (!is.null(results$rna_protein)) {
        rp <- results$rna_protein
        summary_list$rna_protein <- list(
            mean_correlation = rp$correlations$summary$mean_cor,
            n_high_te = rp$translation_efficiency$summary$n_high_te,
            n_low_te = rp$translation_efficiency$summary$n_low_te,
            n_ptr_genes = rp$post_transcriptional$n_ptr
        )
    }

    # TF activity
    if (!is.null(results$tf_activity)) {
        tf <- results$tf_activity
        n_sig_tfs <- ifelse(!is.null(tf$associations),
                            sum(tf$associations$padj < 0.05, na.rm = TRUE), NA)
        summary_list$tf_activity <- list(
            n_tfs = nrow(tf$activity_matrix),
            n_significant = n_sig_tfs
        )
    }

    # Regulatory network
    if (!is.null(results$regulatory_network)) {
        rn <- results$regulatory_network
        summary_list$regulatory_network <- list(
            method = rn$method,
            n_edges = nrow(rn$edges),
            n_hubs = sum(rn$hubs$is_hub, na.rm = TRUE)
        )
    }

    # Mediation
    if (!is.null(results$mediation)) {
        med <- results$mediation
        summary_list$mediation <- list(
            n_genes_tested = nrow(med$mediation_df),
            n_mediated = med$n_mediated
        )
    }

    # Multi-layer mediation
    if (!is.null(results$multilayer_mediation)) {
        ml <- results$multilayer_mediation
        summary_list$multilayer_mediation <- list(
            n_paths_tested = nrow(ml$multilayer_df),
            n_significant = ml$n_significant
        )
    }

    # Protein-metabolite mediation
    if (!is.null(results$prot_metab_mediation)) {
        pm <- results$prot_metab_mediation
        summary_list$prot_metab_mediation <- list(
            n_paths_tested = nrow(pm$mediation_df),
            n_significant = pm$n_significant
        )
    }

    # Pathway activity
    if (!is.null(results$pathway_activity)) {
        pa <- results$pathway_activity
        n_sig_pw <- ifelse(!is.null(pa$associations),
                           sum(pa$associations$padj < 0.05, na.rm = TRUE), NA)
        summary_list$pathway_activity <- list(
            n_pathways = pa$n_pathways,
            n_significant = n_sig_pw,
            method = pa$method
        )
    }

    # Protein-metabolite correlation network
    if (!is.null(results$protein_metabolite_network)) {
        pmn <- results$protein_metabolite_network
        summary_list$protein_metabolite_network <- list(
            n_proteins = pmn$n_proteins,
            n_metabolites = pmn$n_metabolites,
            n_significant_pairs = pmn$n_significant
        )
    }

    # COSMOS
    if (!is.null(results$cosmos)) {
        cosmos <- results$cosmos
        summary_list$cosmos <- list(
            n_nodes = cosmos$n_nodes %||% NA,
            n_edges = cosmos$n_edges %||% NA
        )
    }

    # Save summary
    if (!is.null(out_dir) && length(summary_list) > 0) {
        flat_summary <- do.call(rbind, lapply(names(summary_list), function(cat) {
            s <- summary_list[[cat]]
            data.frame(
                category = cat,
                metric = names(s),
                value = as.character(unlist(s)),
                stringsAsFactors = FALSE
            )
        }))

        write.csv(flat_summary,
                  file.path(out_dir, "tables", "mechanistic_analysis_summary.csv"),
                  row.names = FALSE)
    }

    message("=== Mechanistic Analysis Summary ===")
    for (cat in names(summary_list)) {
        message("  ", cat, ": ", paste(names(summary_list[[cat]]), "=",
                                        unlist(summary_list[[cat]]),
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
.mae_to_legacy_mechanistic <- function(mae, de_results = NULL, gene_protein_mapping = NULL) {

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

    list(
        mae = mae,
        harmonized_omics = harmonized_omics,
        metadata = as.data.frame(SummarizedExperiment::colData(mae)),
        common_samples = colnames(mae),
        gene_mapping = gene_protein_mapping
    )
}


# =============================================================================
# Protein-Metabolite Correlation Analysis
# Works with proteomics + metabolomics (no transcriptomics needed)
# =============================================================================

# =============================================================================
# Correlation Pair Annotation (enzyme/TF/gene context)
# =============================================================================

#' Map UniProt IDs to gene symbols and EC numbers via OrgDb
#' @param uniprot_ids Character vector of UniProt accessions
#' @param organism Organism name (e.g. "Rattus norvegicus")
#' @return Data frame with columns: uniprot_id, gene_symbol, ec_number, entrez_id
annotate_proteins_orgdb <- function(uniprot_ids, organism) {
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) return(NULL)

    org_info <- get_organism_info(organism)
    if (!isTRUE(org_info$supported) || is.na(org_info$orgdb)) return(NULL)

    orgdb_pkg <- org_info$orgdb
    if (!requireNamespace(orgdb_pkg, quietly = TRUE)) return(NULL)

    orgdb <- tryCatch(
        getExportedValue(orgdb_pkg, orgdb_pkg),
        error = function(e) NULL
    )
    if (is.null(orgdb)) return(NULL)

    avail_cols <- AnnotationDbi::columns(orgdb)
    want_cols <- intersect(c("SYMBOL", "ENZYME", "ENTREZID", "GENENAME"), avail_cols)
    if (!"SYMBOL" %in% want_cols) return(NULL)

    # Filter to keys that exist in OrgDb
    valid_keys <- tryCatch(
        AnnotationDbi::keys(orgdb, keytype = "UNIPROT"),
        error = function(e) character(0)
    )
    usable <- unique(uniprot_ids[uniprot_ids %in% valid_keys])
    if (length(usable) == 0) return(NULL)

    raw <- tryCatch(
        AnnotationDbi::select(orgdb, keys = usable, columns = want_cols,
                              keytype = "UNIPROT"),
        error = function(e) NULL
    )
    if (is.null(raw) || nrow(raw) == 0) return(NULL)

    # Collapse multiple EC numbers per UniProt ID
    result <- stats::aggregate(
        raw[, setdiff(colnames(raw), "UNIPROT"), drop = FALSE],
        by = list(uniprot_id = raw$UNIPROT),
        FUN = function(x) {
            vals <- unique(x[!is.na(x) & x != ""])
            if (length(vals) == 0) NA_character_ else paste(vals, collapse = "/")
        }
    )

    colnames(result) <- gsub("^SYMBOL$", "gene_symbol", colnames(result))
    colnames(result) <- gsub("^ENZYME$", "ec_number", colnames(result))
    colnames(result) <- gsub("^ENTREZID$", "entrez_id", colnames(result))
    colnames(result) <- gsub("^GENENAME$", "gene_name", colnames(result))

    result
}


#' Check enzyme-metabolite substrate/product links via KEGG
#'
#' For each EC number, queries KEGG for associated compounds and cross-references
#' against metabolite names resolved to KEGG compound IDs.
#'
#' @param ec_numbers Character vector of EC numbers (one per protein row, may be NA)
#' @param metabolite_names Character vector of metabolite names (one per row)
#' @param kegg_code KEGG organism code (e.g. "rno")
#' @return Character vector: "Yes (EC:x.x.x.x)" where a link exists, else NA
match_enzyme_metabolite_kegg <- function(ec_numbers, metabolite_names, kegg_code) {

    n <- length(ec_numbers)
    result <- rep(NA_character_, n)

    if (!requireNamespace("KEGGREST", quietly = TRUE)) return(result)

    unique_ecs <- unique(ec_numbers[!is.na(ec_numbers) & ec_numbers != ""])
    unique_metabs <- unique(metabolite_names[!is.na(metabolite_names)])

    if (length(unique_ecs) == 0 || length(unique_metabs) == 0) return(result)

    # Step 1: Map metabolite names -> KEGG compound IDs
    message("  Mapping metabolite names to KEGG compound IDs...")
    metab_to_cpd <- list()
    for (mname in unique_metabs) {
        cpds <- tryCatch({
            Sys.sleep(0.15)
            hits <- KEGGREST::keggFind("compound", mname)
            if (length(hits) > 0) names(hits) else character(0)
        }, error = function(e) character(0))
        # Strip "cpd:" prefix
        cpds <- sub("^cpd:", "", cpds)
        if (length(cpds) > 0) metab_to_cpd[[mname]] <- cpds
    }

    if (length(metab_to_cpd) == 0) {
        message("  No metabolites matched KEGG compounds")
        return(result)
    }
    message("  Matched ", length(metab_to_cpd), "/", length(unique_metabs),
            " metabolites to KEGG compounds")

    # Step 2: For each unique EC, get linked compound IDs
    message("  Querying KEGG enzyme-compound links...")
    ec_to_cpd <- list()
    # Handle slash-separated multi-EC entries (e.g. "1.2.3.4/2.3.4.5")
    all_single_ecs <- unique(unlist(strsplit(unique_ecs, "/")))
    for (ec in all_single_ecs) {
        cpds <- tryCatch({
            Sys.sleep(0.15)
            links <- KEGGREST::keggLink("compound", paste0("ec:", ec))
            if (length(links) > 0) sub("^cpd:", "", links) else character(0)
        }, error = function(e) character(0))
        if (length(cpds) > 0) ec_to_cpd[[ec]] <- cpds
    }

    if (length(ec_to_cpd) == 0) {
        message("  No EC-compound links found in KEGG")
        return(result)
    }
    message("  Found compound links for ", length(ec_to_cpd), "/",
            length(all_single_ecs), " EC numbers")

    # Step 3: Cross-reference
    for (i in seq_len(n)) {
        ec_val <- ec_numbers[i]
        met_val <- metabolite_names[i]
        if (is.na(ec_val) || ec_val == "" || is.na(met_val)) next

        met_cpds <- metab_to_cpd[[met_val]]
        if (is.null(met_cpds) || length(met_cpds) == 0) next

        # Check each EC in this protein (may be multi-EC)
        ecs <- unlist(strsplit(ec_val, "/"))
        for (ec in ecs) {
            ec_cpds <- ec_to_cpd[[ec]]
            if (!is.null(ec_cpds) && any(met_cpds %in% ec_cpds)) {
                result[i] <- paste0("Yes (EC:", ec, ")")
                break
            }
        }
    }

    result
}


#' Check whether gene symbols are known transcription factors (dorothea)
#' @param gene_symbols Character vector of gene symbols
#' @return Logical vector, TRUE if the gene is a known TF
check_tf_status <- function(gene_symbols) {
    result <- rep(FALSE, length(gene_symbols))
    regulons <- tryCatch(get_tf_regulons("dorothea"), error = function(e) NULL)
    if (is.null(regulons)) return(result)

    tf_names <- toupper(unique(regulons$tf))
    result <- toupper(gene_symbols) %in% tf_names
    result[is.na(gene_symbols)] <- NA
    result
}


#' Annotate protein-metabolite correlation pairs with biological context
#'
#' Adds gene symbol, EC number, enzyme-metabolite link, and TF status columns.
#'
#' @param sig_pairs Data frame of significant correlation pairs
#' @param config Pipeline config (needs config$global$organism)
#' @return sig_pairs with added annotation columns
annotate_correlation_pairs <- function(sig_pairs, config) {

    if (nrow(sig_pairs) == 0) {
        sig_pairs$gene_symbol <- character(0)
        sig_pairs$ec_number <- character(0)
        sig_pairs$enzyme_metabolite_link <- character(0)
        sig_pairs$is_transcription_factor <- logical(0)
        return(sig_pairs)
    }

    organism <- config$global$organism %||%
                config$organism %||%
                config$annotation$organism %||%
                "Rattus norvegicus"

    message("  Annotating protein-metabolite pairs (organism: ", organism, ")...")

    # 1. OrgDb annotation: UniProt -> gene symbol, EC number
    prot_annot <- tryCatch(
        annotate_proteins_orgdb(unique(sig_pairs$protein), organism),
        error = function(e) {
            message("  OrgDb annotation failed: ", e$message)
            NULL
        }
    )

    if (!is.null(prot_annot)) {
        sig_pairs <- merge(sig_pairs, prot_annot, by.x = "protein",
                           by.y = "uniprot_id", all.x = TRUE, sort = FALSE)
        # Restore original sort order
        sig_pairs <- sig_pairs[order(-sig_pairs$abs_cor), , drop = FALSE]
    } else {
        sig_pairs$gene_symbol <- NA_character_
        sig_pairs$ec_number <- NA_character_
    }

    # Ensure columns exist
    if (!"gene_symbol" %in% colnames(sig_pairs))
        sig_pairs$gene_symbol <- NA_character_
    if (!"ec_number" %in% colnames(sig_pairs))
        sig_pairs$ec_number <- NA_character_

    # 2. Enzyme-metabolite link via KEGG
    org_info <- get_organism_info(organism)
    sig_pairs$enzyme_metabolite_link <- tryCatch(
        match_enzyme_metabolite_kegg(
            sig_pairs$ec_number, sig_pairs$metabolite, org_info$kegg
        ),
        error = function(e) {
            message("  KEGG enzyme-metabolite lookup failed: ", e$message)
            rep(NA_character_, nrow(sig_pairs))
        }
    )

    # 3. Transcription factor status via dorothea
    sig_pairs$is_transcription_factor <- tryCatch(
        check_tf_status(sig_pairs$gene_symbol),
        error = function(e) {
            message("  TF status check failed: ", e$message)
            rep(NA, nrow(sig_pairs))
        }
    )

    n_ec <- sum(!is.na(sig_pairs$ec_number) & sig_pairs$ec_number != "", na.rm = TRUE)
    n_link <- sum(!is.na(sig_pairs$enzyme_metabolite_link), na.rm = TRUE)
    n_tf <- sum(sig_pairs$is_transcription_factor, na.rm = TRUE)
    message("  Annotation complete: ", n_ec, " enzymes, ", n_link,
            " enzyme-metabolite links, ", n_tf, " transcription factors")

    sig_pairs
}


#' Analyze correlations between DE proteins and metabolites
#'
#' Computes pairwise Spearman correlations between significant proteins and
#' metabolites, identifies significant protein-metabolite pairs, and generates
#' a correlation network plot.
#'
#' @param harmonized Named list with proteomics and metabolomics data
#' @param de_results Named list of DE results per omics
#' @param mc Mechanistic config
#' @param config Full pipeline config
#' @param out_dir Output directory
#' @return List with correlation matrix, significant pairs, and plots
analyze_protein_metabolite_correlations <- function(harmonized, de_results,
                                                     mc, config, out_dir = NULL) {

    message("Running protein-metabolite correlation analysis...")

    prot <- harmonized$proteomics
    metab <- harmonized$metabolomics

    prot_mat <- prot$normalized_matrix %||% prot$expr_work
    metab_mat <- metab$normalized_matrix %||% metab$expr_work

    if (is.null(prot_mat) || is.null(metab_mat)) {
        message("  Missing expression matrices. Skipping.")
        return(NULL)
    }

    # Align samples
    common_samples <- intersect(colnames(prot_mat), colnames(metab_mat))
    if (length(common_samples) < 4) {
        message("  Too few common samples (", length(common_samples), "). Skipping.")
        return(NULL)
    }

    prot_mat <- prot_mat[, common_samples, drop = FALSE]
    metab_mat <- metab_mat[, common_samples, drop = FALSE]

    # Select features: use DE features if available, otherwise top variable
    prot_ids <- select_features_for_correlation(prot_mat, de_results, "proteomics",
                                                 max_n = 100)
    metab_ids <- select_features_for_correlation(metab_mat, de_results, "metabolomics",
                                                  max_n = 50)

    message("  Selected ", length(prot_ids), " proteins (from ",
            nrow(prot_mat), " total) and ", length(metab_ids),
            " metabolites (from ", nrow(metab_mat), " total)")
    message("  Selection method: DE-significant features if available (pass_any filter), ",
            "otherwise top variable features (by variance)")

    if (length(prot_ids) < 2 || length(metab_ids) < 2) {
        message("  Too few features for correlation analysis. Skipping.")
        return(NULL)
    }

    prot_sub <- prot_mat[prot_ids, , drop = FALSE]
    metab_sub <- metab_mat[metab_ids, , drop = FALSE]

    message("  Computing correlations: ", length(prot_ids), " proteins x ",
            length(metab_ids), " metabolites")

    # Compute cross-omics correlations
    cor_mat <- cor(t(prot_sub), t(metab_sub), method = "spearman",
                   use = "pairwise.complete.obs")

    # Compute p-values
    n <- length(common_samples)
    t_stat <- cor_mat * sqrt((n - 2) / (1 - cor_mat^2))
    pval_mat <- 2 * pt(-abs(t_stat), df = n - 2)

    # Flatten to pairs
    pairs_df <- expand.grid(
        protein = rownames(cor_mat),
        metabolite = colnames(cor_mat),
        stringsAsFactors = FALSE
    )
    pairs_df$correlation <- as.vector(cor_mat)
    pairs_df$pvalue <- as.vector(pval_mat)
    pairs_df$padj <- p.adjust(pairs_df$pvalue, method = "BH")
    pairs_df$abs_cor <- abs(pairs_df$correlation)

    # Filter significant pairs
    sig_pairs <- pairs_df[!is.na(pairs_df$padj) &
                          pairs_df$padj < 0.2 &
                          pairs_df$abs_cor > 0.7, , drop = FALSE]
    sig_pairs <- sig_pairs[order(-sig_pairs$abs_cor), ]

    message("  Found ", nrow(sig_pairs), " significant protein-metabolite pairs",
            " (padj < 0.2, |r| > 0.7)")

    # Annotate significant pairs with biological context
    sig_pairs <- tryCatch({
        annotate_correlation_pairs(sig_pairs, config)
    }, error = function(e) {
        message("  Pair annotation failed: ", e$message)
        sig_pairs
    })

    # Save results
    if (!is.null(out_dir)) {
        tables_dir <- file.path(out_dir, "tables")
        plots_dir <- file.path(out_dir, "plots")
        dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

        write.csv(sig_pairs,
                  file.path(tables_dir, "protein_metabolite_significant_pairs.csv"),
                  row.names = FALSE)
        write.csv(pairs_df,
                  file.path(tables_dir, "protein_metabolite_all_correlations.csv"),
                  row.names = FALSE)

        # Correlation heatmap of top pairs
        if (nrow(sig_pairs) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
            top_prots <- unique(sig_pairs$protein[1:min(30, nrow(sig_pairs))])
            top_metabs <- unique(sig_pairs$metabolite[1:min(20, nrow(sig_pairs))])
            sub_cor <- cor_mat[intersect(rownames(cor_mat), top_prots),
                               intersect(colnames(cor_mat), top_metabs), drop = FALSE]

            if (nrow(sub_cor) > 1 && ncol(sub_cor) > 1) {
                plot_df <- expand.grid(
                    protein = rownames(sub_cor),
                    metabolite = colnames(sub_cor),
                    stringsAsFactors = FALSE
                )
                plot_df$correlation <- as.vector(sub_cor)

                p <- ggplot2::ggplot(plot_df,
                    ggplot2::aes(x = metabolite, y = protein, fill = correlation)) +
                    ggplot2::geom_tile() +
                    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                                  high = "#B2182B", midpoint = 0,
                                                  limits = c(-1, 1)) +
                    ggplot2::theme_minimal() +
                    ggplot2::theme(
                        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
                        axis.text.y = ggplot2::element_text(size = 7)
                    ) +
                    ggplot2::labs(title = "Protein-Metabolite Correlations",
                                 x = "Metabolite", y = "Protein",
                                 fill = "Spearman r")

                ggplot2::ggsave(file.path(plots_dir, "protein_metabolite_correlation_heatmap.png"),
                                p, width = 10, height = 8, dpi = 150)
                message("  Saved protein-metabolite correlation heatmap")
            }
        }
    }

    list(
        correlation_matrix = cor_mat,
        pvalue_matrix = pval_mat,
        significant_pairs = sig_pairs,
        all_pairs = pairs_df,
        n_proteins = length(prot_ids),
        n_metabolites = length(metab_ids),
        n_significant = nrow(sig_pairs)
    )
}


#' Select features for correlation analysis
#' Uses DE features if available, otherwise top variable features
select_features_for_correlation <- function(expr_mat, de_results, omics_type,
                                             max_n = 100) {
    # Try DE features first
    de <- de_results[[omics_type]]
    if (!is.null(de)) {
        summary_df <- de$summary_df %||% de$de_table
        if (!is.null(summary_df)) {
            pass_col <- grep("pass_any", colnames(summary_df), value = TRUE)[1]
            id_col <- grep("feature_id|FeatureID|Protein.Group", colnames(summary_df), value = TRUE)[1]
            if (!is.na(pass_col) && !is.na(id_col)) {
                sig_ids <- summary_df[[id_col]][summary_df[[pass_col]] == 1]
                sig_ids <- sig_ids[sig_ids %in% rownames(expr_mat)]
                if (length(sig_ids) >= 2) {
                    message("  ", omics_type, ": using ", length(sig_ids), " DE features")
                    return(head(sig_ids, max_n))
                }
            }
        }
    }

    # Fallback: top variable features
    vars <- apply(expr_mat, 1, var, na.rm = TRUE)
    top_ids <- names(sort(vars, decreasing = TRUE))[1:min(max_n, nrow(expr_mat))]
    message("  ", omics_type, ": using top ", length(top_ids), " variable features")
    top_ids
}
