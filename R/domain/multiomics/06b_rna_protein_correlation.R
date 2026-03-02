#' RNA-Protein Correlation Analysis
#'
#' Computes expression correlations and log2FC concordance between RNA and
#' protein data for genes that map between the two omics layers.
#'
#' @description
#' This module provides:
#'   1. Per-gene expression correlation (abundance vs abundance)
#'   2. Differential concordance (log2FC vs log2FC)
#'   3. Translation efficiency analysis
#'
#' @name rna_protein_correlation
NULL


#' Run RNA-protein correlation analysis
#'
#' @param mae MultiAssayExperiment object with aligned samples
#' @param de_results Named list of DE results per omics (requires transcriptomics and proteomics)
#' @param gene_protein_mapping Data frame mapping genes to proteins (must have columns:
#'   omics, feature_id, gene_symbol)
#' @param config Full config object
#' @param out_dir Output directory for results
#' @return List with: correlations (per-gene), summary, de_concordance
run_rna_protein_correlation <- function(mae, de_results = NULL,
                                         gene_protein_mapping = NULL,
                                         config, out_dir = NULL) {

    message("Running RNA-Protein correlation analysis...")

    # Convert MAE to legacy format for easier access
    mae_data <- .mae_to_legacy(mae, de_results, gene_protein_mapping)

    # Create output directories if specified
    if (!is.null(out_dir)) {
        dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
    }

    # =========================================================================
    # Part 1: Expression Correlation (Abundance vs Abundance per gene)
    # =========================================================================

    rna_mat <- mae_data$harmonized_omics$transcriptomics$normalized_matrix
    prot_mat <- mae_data$harmonized_omics$proteomics$normalized_matrix
    gene_mapping <- mae_data$gene_mapping

    cor_df <- NULL
    summary_stats <- list(
        n_genes = 0,
        mean_cor = NA,
        median_cor = NA,
        pct_positive = NA,
        n_significant = 0
    )

    if (!is.null(gene_mapping) && !is.null(rna_mat) && !is.null(prot_mat)) {
        # Subset mapping for each omic
        rna_genes <- gene_mapping[gene_mapping$omics == "transcriptomics", ]
        prot_genes <- gene_mapping[gene_mapping$omics == "proteomics", ]

        # Find common gene symbols
        common_genes <- intersect(rna_genes$gene_symbol, prot_genes$gene_symbol)
        common_genes <- common_genes[!is.na(common_genes)]

        if (length(common_genes) >= 10) {
            # Map to feature IDs
            rna_features <- rna_genes$feature_id[match(common_genes, rna_genes$gene_symbol)]
            prot_features <- prot_genes$feature_id[match(common_genes, prot_genes$gene_symbol)]

            # Align samples
            common_samples <- intersect(colnames(rna_mat), colnames(prot_mat))

            if (length(common_samples) > 2) {
                message("  Computing per-gene correlations for ", length(common_genes), " genes...")

                # Compute per-gene Pearson correlation
                correlations <- numeric(length(common_genes))
                names(correlations) <- common_genes

                for (i in seq_along(common_genes)) {
                    rna_expr <- rna_mat[rna_features[i], common_samples]
                    prot_expr <- prot_mat[prot_features[i], common_samples]
                    if (sd(rna_expr, na.rm = TRUE) > 0 && sd(prot_expr, na.rm = TRUE) > 0) {
                        correlations[i] <- cor(rna_expr, prot_expr, use = "pairwise.complete.obs")
                    } else {
                        correlations[i] <- NA
                    }
                }
                correlations <- correlations[!is.na(correlations)]

                if (length(correlations) > 0) {
                    cor_df <- data.frame(
                        gene_symbol = names(correlations),
                        correlation = correlations,
                        stringsAsFactors = FALSE
                    )

                    # Save results
                    if (!is.null(out_dir)) {
                        write.csv(cor_df,
                                  file.path(out_dir, "tables", "rna_protein_correlations.csv"),
                                  row.names = FALSE)

                        # Plot histogram
                        p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = correlation)) +
                            ggplot2::geom_histogram(bins = 50, fill = "steelblue",
                                                    color = "white", alpha = 0.7) +
                            ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
                            ggplot2::theme_minimal() +
                            ggplot2::labs(
                                title = "RNA-Protein Expression Correlation Distribution",
                                x = "Pearson Correlation",
                                y = "Count"
                            )
                        ggplot2::ggsave(file.path(out_dir, "plots", "rna_protein_concordance.png"),
                                        p, width = 8, height = 6, dpi = 150)
                    }

                    # Update summary
                    summary_stats$n_genes <- length(correlations)
                    summary_stats$mean_cor <- mean(correlations)
                    summary_stats$median_cor <- median(correlations)
                    summary_stats$pct_positive <- 100 * mean(correlations > 0)
                    summary_stats$n_significant <- sum(abs(correlations) > 0.5)

                    message("  Expression correlation: mean r = ", round(summary_stats$mean_cor, 3),
                            ", ", summary_stats$n_significant, " genes with |r| > 0.5")
                }
            }
        }
    } else {
        message("  Skipping expression correlation: missing gene mapping or matrices")
    }

    # =========================================================================
    # Part 2: Differential Concordance (Log2FC vs Log2FC)
    # =========================================================================

    de_concordance_df <- NULL

    rna_obj <- mae_data$harmonized_omics$transcriptomics
    prot_obj <- mae_data$harmonized_omics$proteomics

    rna_de <- rna_obj$de_table %||% rna_obj$da_table
    prot_de <- prot_obj$da_table %||% prot_obj$de_table

    if (!is.null(rna_de) && !is.null(prot_de) && !is.null(gene_mapping)) {
        de_concordance_df <- compute_de_concordance(
            rna_de = rna_de,
            prot_de = prot_de,
            gene_mapping = gene_mapping,
            out_dir = out_dir
        )
    } else {
        message("  Skipping DE concordance: missing DE tables or gene mapping")
    }

    message("RNA-Protein correlation analysis complete")

    list(
        rna_protein = list(
            correlations = cor_df,
            summary = summary_stats
        ),
        de_concordance = de_concordance_df
    )
}


#' Compute differential expression concordance between RNA and protein
#'
#' @param rna_de RNA DE table (must have feature_id and log2FC columns)
#' @param prot_de Protein DE/DA table (must have feature_id and log2FC columns)
#' @param gene_mapping Gene mapping table
#' @param out_dir Output directory (optional)
#' @return Data frame with merged DE results and concordance metrics
compute_de_concordance <- function(rna_de, prot_de, gene_mapping, out_dir = NULL) {

    # Helper to get padj column name
    get_padj_col <- function(df) {
        if ("adj.P.Val" %in% colnames(df)) return("adj.P.Val")
        if ("padj" %in% colnames(df)) return("padj")
        if ("FDR" %in% colnames(df)) return("FDR")
        return(NULL)
    }

    # Prepare RNA
    rna_map_sub <- gene_mapping[gene_mapping$omics == "transcriptomics", ]

    if (!"log2FC" %in% colnames(rna_de)) {
        message("  RNA DE table missing log2FC column")
        return(NULL)
    }

    rna_de$gene_symbol <- rna_map_sub$gene_symbol[match(rna_de$feature_id, rna_map_sub$feature_id)]

    padj_col <- get_padj_col(rna_de)
    cols_to_keep <- c("gene_symbol", "log2FC")
    if (!is.null(padj_col)) cols_to_keep <- c(cols_to_keep, padj_col)

    rna_de_clean <- rna_de[!is.na(rna_de$gene_symbol), cols_to_keep, drop = FALSE]

    if (!is.null(padj_col)) {
        colnames(rna_de_clean) <- c("gene_symbol", "rna_log2FC", "rna_padj")
    } else {
        colnames(rna_de_clean) <- c("gene_symbol", "rna_log2FC")
        rna_de_clean$rna_padj <- NA
    }

    # Prepare Protein
    prot_map_sub <- gene_mapping[gene_mapping$omics == "proteomics", ]

    if (!"log2FC" %in% colnames(prot_de)) {
        message("  Protein DE table missing log2FC column")
        return(NULL)
    }

    prot_de$gene_symbol <- prot_map_sub$gene_symbol[match(prot_de$feature_id, prot_map_sub$feature_id)]

    padj_col <- get_padj_col(prot_de)
    cols_to_keep <- c("gene_symbol", "log2FC")
    if (!is.null(padj_col)) cols_to_keep <- c(cols_to_keep, padj_col)

    prot_de_clean <- prot_de[!is.na(prot_de$gene_symbol), cols_to_keep, drop = FALSE]

    if (!is.null(padj_col)) {
        colnames(prot_de_clean) <- c("gene_symbol", "protein_log2FC", "protein_padj")
    } else {
        colnames(prot_de_clean) <- c("gene_symbol", "protein_log2FC")
        prot_de_clean$protein_padj <- NA
    }

    # Merge
    de_merged <- merge(rna_de_clean, prot_de_clean, by = "gene_symbol")

    if (nrow(de_merged) < 10) {
        message("  Fewer than 10 genes merged, skipping DE concordance")
        return(NULL)
    }

    message("  Computing DE concordance for ", nrow(de_merged), " genes...")

    # Define significance categories (FDR < 0.05)
    rna_p <- de_merged$rna_padj
    rna_p[is.na(rna_p)] <- 1

    prot_p <- de_merged$protein_padj
    prot_p[is.na(prot_p)] <- 1

    is_sig_rna <- rna_p < 0.05
    is_sig_prot <- prot_p < 0.05

    de_merged$category <- "Non-sig"
    de_merged$category[is_sig_rna & !is_sig_prot] <- "Sig RNA (Gold)"
    de_merged$category[!is_sig_rna & is_sig_prot] <- "Sig Protein (Purple)"
    de_merged$category[is_sig_rna & is_sig_prot] <- "Sig Both (Red)"

    de_merged$category <- factor(de_merged$category,
        levels = c("Non-sig", "Sig RNA (Gold)", "Sig Protein (Purple)", "Sig Both (Red)")
    )

    de_merged$concordant <- sign(de_merged$rna_log2FC) == sign(de_merged$protein_log2FC)

    # Compute correlations
    is_sig_union <- is_sig_rna | is_sig_prot

    cor_all <- cor(de_merged$rna_log2FC, de_merged$protein_log2FC, use = "complete.obs")

    cor_sig <- NA
    if (sum(is_sig_union) > 3) {
        cor_sig <- cor(de_merged$rna_log2FC[is_sig_union],
                       de_merged$protein_log2FC[is_sig_union],
                       use = "complete.obs")
    }

    # Translation efficiency
    de_merged$te_log2FC <- de_merged$protein_log2FC - de_merged$rna_log2FC

    message("  DE concordance: r = ", round(cor_all, 3),
            ", ", sum(de_merged$concordant), "/", nrow(de_merged), " concordant")

    # Save results and plots
    if (!is.null(out_dir)) {
        write.csv(de_merged,
                  file.path(out_dir, "tables", "rna_protein_de_concordance.csv"),
                  row.names = FALSE)

        write.csv(de_merged[, c("gene_symbol", "rna_log2FC", "protein_log2FC", "te_log2FC", "category")],
                  file.path(out_dir, "tables", "translation_efficiency.csv"),
                  row.names = FALSE)

        # Scatter plot
        custom_colors <- c(
            "Non-sig" = "gray80",
            "Sig RNA (Gold)" = "gold3",
            "Sig Protein (Purple)" = "purple",
            "Sig Both (Red)" = "red"
        )

        subtitle_text <- paste0("All genes: r = ", round(cor_all, 3))
        if (!is.na(cor_sig)) {
            subtitle_text <- paste0(subtitle_text, " | Sig. Union (FDR<0.05): r = ", round(cor_sig, 3))
        }

        p1 <- ggplot2::ggplot(de_merged, ggplot2::aes(x = rna_log2FC, y = protein_log2FC)) +
            ggplot2::geom_vline(xintercept = 0, color = "gray90") +
            ggplot2::geom_hline(yintercept = 0, color = "gray90") +
            ggplot2::geom_point(ggplot2::aes(color = category), alpha = 0.6, size = 1.5) +
            ggplot2::geom_smooth(method = "lm", color = "black", linetype = "dashed",
                                 se = FALSE, linewidth = 0.5) +
            ggplot2::scale_color_manual(values = custom_colors) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = "Differential Concordance: RNA vs Protein",
                subtitle = subtitle_text,
                x = "RNA log2 Fold Change",
                y = "Protein log2 Fold Change",
                color = "Significance"
            ) +
            ggplot2::theme(legend.position = "bottom")

        ggplot2::ggsave(file.path(out_dir, "plots", "rna_protein_de_scatter.png"),
                        p1, width = 8, height = 7, dpi = 300)

        # TE histogram
        p2 <- ggplot2::ggplot(de_merged, ggplot2::aes(x = te_log2FC)) +
            ggplot2::geom_histogram(bins = 40, fill = "darkcyan", color = "white", alpha = 0.8) +
            ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = "Distribution of Translation Efficiency Changes",
                subtitle = "log2(TE) = log2(Protein FC) - log2(RNA FC)",
                x = "TE log2 Fold Change",
                y = "Count"
            )
        ggplot2::ggsave(file.path(out_dir, "plots", "translation_efficiency_hist.png"),
                        p2, width = 8, height = 6)

        # TE scatter
        te_cor <- cor(de_merged$rna_log2FC, de_merged$protein_log2FC, use = "complete.obs")
        te_subtitle <- paste0(
            "r = ", round(te_cor, 3),
            " | Red = protein > RNA (high TE), Blue = protein < RNA (low TE)"
        )

        p3 <- ggplot2::ggplot(de_merged, ggplot2::aes(x = rna_log2FC, y = protein_log2FC, color = te_log2FC)) +
            ggplot2::geom_hline(yintercept = 0, color = "gray90") +
            ggplot2::geom_vline(xintercept = 0, color = "gray90") +
            ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
            ggplot2::geom_point(alpha = 0.7, size = 1.5) +
            ggplot2::scale_color_gradient2(
                low = "blue", mid = "grey90", high = "red",
                midpoint = 0, name = "log2(TE)\n(Prot FC - RNA FC)"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
                title = "Translation Efficiency: RNA vs Protein log2FC",
                subtitle = te_subtitle,
                x = "RNA log2 Fold Change",
                y = "Protein log2 Fold Change"
            )
        ggplot2::ggsave(file.path(out_dir, "plots", "translation_efficiency_scatter.png"),
                        p3, width = 8, height = 7, dpi = 300)
    }

    de_merged
}


#' Convert MAE to legacy mae_data format (internal)
#'
#' Adapter function to minimize changes when porting from legacy code.
#'
#' @param mae MultiAssayExperiment object
#' @param de_results Named list of DE results per omics
#' @param gene_protein_mapping Gene-protein mapping table
#' @return List in legacy mae_data format
.mae_to_legacy <- function(mae, de_results = NULL, gene_protein_mapping = NULL) {

    harmonized_omics <- lapply(names(mae@ExperimentList), function(nm) {
        exp_data <- mae@ExperimentList[[nm]]
        de <- if (!is.null(de_results) && nm %in% names(de_results)) de_results[[nm]] else NULL

        list(
            normalized_matrix = as.matrix(SummarizedExperiment::assay(exp_data)),
            de_table = de,
            da_table = de,
            feature_annotation = as.data.frame(SummarizedExperiment::rowData(exp_data))
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
