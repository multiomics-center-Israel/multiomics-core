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

    if (!is.null(rna_mat) && !is.null(prot_mat)) {
        # Check if features are already harmonized (same row names in RNA and protein)
        harmonized <- identical(sort(rownames(rna_mat)), sort(rownames(prot_mat)))

        if (harmonized) {
            # MAE was pre-harmonized: features are already 1:1 aligned
            common_features <- intersect(rownames(rna_mat), rownames(prot_mat))
            common_samples <- intersect(colnames(rna_mat), colnames(prot_mat))

            if (length(common_features) >= 10 && length(common_samples) > 2) {
                message("  Computing per-gene correlations for ", length(common_features),
                        " harmonized gene-protein pairs...")

                correlations <- numeric(length(common_features))
                names(correlations) <- common_features

                for (i in seq_along(common_features)) {
                    feat <- common_features[i]
                    rna_expr <- rna_mat[feat, common_samples]
                    prot_expr <- prot_mat[feat, common_samples]
                    rna_sd <- sd(rna_expr, na.rm = TRUE)
                    prot_sd <- sd(prot_expr, na.rm = TRUE)
                    if (!is.na(rna_sd) && rna_sd > 0 && !is.na(prot_sd) && prot_sd > 0) {
                        correlations[i] <- cor(rna_expr, prot_expr, use = "pairwise.complete.obs")
                    } else {
                        correlations[i] <- NA
                    }
                }
                correlations <- correlations[!is.na(correlations)]
            }
        } else if (!is.null(gene_mapping)) {
            # Use gene mapping to find matched features
            rna_genes <- gene_mapping[gene_mapping$omics == "transcriptomics", ]
            prot_genes <- gene_mapping[gene_mapping$omics == "proteomics", ]

            common_genes <- intersect(rna_genes$gene_symbol, prot_genes$gene_symbol)
            common_genes <- common_genes[!is.na(common_genes)]

            if (length(common_genes) >= 10) {
                rna_features <- rna_genes$feature_id[match(common_genes, rna_genes$gene_symbol)]
                prot_features <- prot_genes$feature_id[match(common_genes, prot_genes$gene_symbol)]

                # Filter to features actually present in matrices
                valid <- rna_features %in% rownames(rna_mat) & prot_features %in% rownames(prot_mat)
                common_genes <- common_genes[valid]
                rna_features <- rna_features[valid]
                prot_features <- prot_features[valid]

                common_samples <- intersect(colnames(rna_mat), colnames(prot_mat))

                if (length(common_genes) >= 10 && length(common_samples) > 2) {
                    message("  Computing per-gene correlations for ", length(common_genes), " genes...")

                    correlations <- numeric(length(common_genes))
                    names(correlations) <- common_genes

                    for (i in seq_along(common_genes)) {
                        rna_expr <- rna_mat[rna_features[i], common_samples]
                        prot_expr <- prot_mat[prot_features[i], common_samples]
                        rna_sd <- sd(rna_expr, na.rm = TRUE)
                        prot_sd <- sd(prot_expr, na.rm = TRUE)
                        if (!is.na(rna_sd) && rna_sd > 0 && !is.na(prot_sd) && prot_sd > 0) {
                            correlations[i] <- cor(rna_expr, prot_expr, use = "pairwise.complete.obs")
                        } else {
                            correlations[i] <- NA
                        }
                    }
                    correlations <- correlations[!is.na(correlations)]
                }
            }
        } else {
            message("  Skipping expression correlation: features not harmonized and no gene mapping")
        }

        # Save results and plots if correlations were computed
        if (exists("correlations") && length(correlations) > 0) {
            cor_df <- data.frame(
                gene_symbol = names(correlations),
                correlation = correlations,
                stringsAsFactors = FALSE
            )

            if (!is.null(out_dir)) {
                write.csv(cor_df,
                          file.path(out_dir, "tables", "rna_protein_correlations.csv"),
                          row.names = FALSE)

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

            summary_stats$n_genes <- length(correlations)
            summary_stats$mean_cor <- mean(correlations)
            summary_stats$median_cor <- median(correlations)
            summary_stats$pct_positive <- 100 * mean(correlations > 0)
            summary_stats$n_significant <- sum(abs(correlations) > 0.5)

            message("  Expression correlation: mean r = ", round(summary_stats$mean_cor, 3),
                    ", ", summary_stats$n_significant, " genes with |r| > 0.5")
        }
    } else {
        message("  Skipping expression correlation: missing RNA or protein matrices")
    }

    # =========================================================================
    # Part 2: Differential Concordance (Log2FC vs Log2FC) — per contrast
    # =========================================================================

    de_concordance_list <- list()

    # Load all precomputed DE tables (one per contrast)
    rna_de_tables <- .load_all_precomputed_de(config, "rna")
    prot_de_tables <- .load_all_precomputed_de(config, "proteomics")

    # Fall back to extracted single-table if precomputed not available
    if (length(rna_de_tables) == 0) {
        rna_obj <- mae_data$harmonized_omics$transcriptomics
        rna_de <- rna_obj$de_table %||% rna_obj$da_table
        if (!is.null(rna_de) && !all(grepl("^\\d+$", head(rna_de$feature_id, 20)))) {
            rna_de_tables <- stats::setNames(list(rna_de), .first_contrast_name(config))
        }
    }
    if (length(prot_de_tables) == 0) {
        prot_obj <- mae_data$harmonized_omics$proteomics
        prot_de <- prot_obj$da_table %||% prot_obj$de_table
        if (!is.null(prot_de) && !all(grepl("^\\d+$", head(prot_de$feature_id, 20)))) {
            prot_de_tables <- stats::setNames(list(prot_de), .first_contrast_name(config))
        }
    }

    # Match RNA and protein DE tables by position (contrast order from config)
    n_contrasts <- min(length(rna_de_tables), length(prot_de_tables))

    if (n_contrasts > 0 && !is.null(gene_mapping)) {
        contrast_names <- names(rna_de_tables)
        if (is.null(contrast_names)) contrast_names <- paste0("contrast_", seq_len(n_contrasts))

        for (ci in seq_len(n_contrasts)) {
            cname <- contrast_names[ci]
            message("  DE concordance for contrast: ", cname)

            # Create the per-contrast directories only after compute_de_concordance()
            # has something to write. Creating them up front left empty
            # plots/ and tables/ behind whenever the join produced nothing, and
            # the report then rendered a bare heading over an empty directory.
            contrast_out_dir <- if (is.null(out_dir)) NULL else
                file.path(out_dir, "per_contrast", cname)

            res <- compute_de_concordance(
                rna_de = rna_de_tables[[ci]],
                prot_de = prot_de_tables[[ci]],
                gene_mapping = gene_mapping,
                out_dir = contrast_out_dir,
                contrast_label = cname
            )

            if (!is.null(res)) {
                de_concordance_list[[cname]] <- res
            }
        }

        # Also save a combined summary plot to the main plots dir
        if (length(de_concordance_list) > 0 && !is.null(out_dir)) {
            .save_combined_de_scatter(de_concordance_list, out_dir)
        }
    } else {
        message("  Skipping DE concordance: missing DE tables or gene mapping")
    }

    message("RNA-Protein correlation analysis complete")

    list(
        rna_protein = list(
            correlations = cor_df,
            summary = summary_stats
        ),
        de_concordance = if (length(de_concordance_list) > 0) de_concordance_list else NULL
    )
}


#' First configured contrast name, for labelling fallback DE tables
#'
#' The precomputed-DE path names each table after its contrast. The MAE fallback
#' used to hard-code \code{"contrast_1"}, which surfaced verbatim as a report
#' heading and as a directory name. Prefer the real contrast label.
#'
#' @param config Full config object.
#' @return Character scalar contrast name; \code{"contrast_1"} only when no
#'   contrast can be resolved from the config.
.first_contrast_name <- function(config) {
    cand <- config$design$contrasts
    if (length(cand) > 0 && nzchar(cand[[1]])) {
        return(gsub("\\s*-\\s*", "_vs_", trimws(as.character(cand[[1]]))))
    }
    "contrast_1"
}


#' Load all precomputed DE tables from config (one per contrast)
#'
#' @param config Full config object
#' @param omics_type "rna" or "proteomics"
#' @return Named list of normalized DE data frames (one per contrast)
.load_all_precomputed_de <- function(config, omics_type) {
    data_dir <- file.path(config$project$dir, config$paths$raw)

    if (omics_type == "proteomics") {
        de_files <- config$modes$proteomics$files$de_table
        contrasts <- config$modes$proteomics$de$contrasts
    } else if (omics_type == "rna") {
        de_files <- config$modes$rna$files$de_table
        contrasts <- config$modes$rna$de$contrasts
    } else {
        return(list())
    }

    if (is.null(de_files) || length(de_files) == 0) return(list())

    result <- list()
    for (i in seq_along(de_files)) {
        de_path <- file.path(data_dir, de_files[[i]])
        if (!file.exists(de_path)) next

        df <- read.csv(de_path, stringsAsFactors = FALSE, check.names = FALSE)
        df <- .normalize_de_file_columns(df)

        if (is.null(df) || !"log2FC" %in% colnames(df)) next

        # Name by contrast if available, otherwise by file
        label <- if (!is.null(contrasts) && i <= length(contrasts)) {
            contrasts[[i]]
        } else {
            tools::file_path_sans_ext(basename(de_files[[i]]))
        }

        message("  Loaded precomputed DE: ", nrow(df), " features from ",
                basename(de_path), " (", label, ")")
        result[[label]] <- df
    }

    result
}


#' Normalize columns in a precomputed DE file
.normalize_de_file_columns <- function(df) {
    # feature ID
    if (!"feature_id" %in% colnames(df)) {
        id_candidates <- c("ID", "FeatureID", "gene_id", "protein_id", "Row.names")
        found_id <- intersect(id_candidates, colnames(df))
        if (length(found_id) > 0) {
            df$feature_id <- df[[found_id[1]]]
        } else {
            df$feature_id <- rownames(df)
        }
    }

    # log2FC
    if (!"log2FC" %in% colnames(df)) {
        fc_candidates <- c("log2FoldChange", "logFC", "log2fc")
        found_fc <- intersect(fc_candidates, colnames(df))
        if (length(found_fc) > 0) df$log2FC <- df[[found_fc[1]]]
    }

    # padj
    if (!"padj" %in% colnames(df)) {
        padj_candidates <- c("adj.P.Val", "FDR", "p.adjust", "qvalue")
        found_padj <- intersect(padj_candidates, colnames(df))
        if (length(found_padj) > 0) df$padj <- df[[found_padj[1]]]
    }

    df
}


#' Save a combined DE scatter plot showing all contrasts as facets
.save_combined_de_scatter <- function(de_concordance_list, out_dir) {
    # Combine all concordance data frames with a contrast column
    combined <- do.call(rbind, lapply(names(de_concordance_list), function(cname) {
        df <- de_concordance_list[[cname]]
        df$contrast <- cname
        df
    }))

    if (nrow(combined) == 0) return(invisible(NULL))

    # Compute per-contrast correlations for subtitles
    cor_labels <- sapply(split(combined, combined$contrast), function(d) {
        r <- cor(d$rna_log2FC, d$protein_log2FC, use = "complete.obs")
        n_conc <- sum(d$concordant, na.rm = TRUE)
        sprintf("r=%.3f, %d/%d concordant", r, n_conc, nrow(d))
    })

    # Add correlation info to contrast labels
    combined$contrast_label <- paste0(combined$contrast, "\n(",
                                       cor_labels[combined$contrast], ")")

    custom_colors <- c(
        "Non-sig" = "gray80",
        "Sig RNA (Gold)" = "gold3",
        "Sig Protein (Purple)" = "purple",
        "Sig Both (Red)" = "red"
    )

    p <- ggplot2::ggplot(combined, ggplot2::aes(x = rna_log2FC, y = protein_log2FC)) +
        ggplot2::geom_vline(xintercept = 0, color = "gray90") +
        ggplot2::geom_hline(yintercept = 0, color = "gray90") +
        ggplot2::geom_point(ggplot2::aes(color = category), alpha = 0.6, size = 1.5) +
        ggplot2::geom_smooth(method = "lm", color = "black", linetype = "dashed",
                             se = FALSE, linewidth = 0.5) +
        ggplot2::scale_color_manual(values = custom_colors) +
        ggplot2::facet_wrap(~ contrast_label, scales = "free") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "RNA-Protein Differential Concordance by Contrast",
            x = "RNA log2 Fold Change",
            y = "Protein log2 Fold Change",
            color = "Significance"
        ) +
        ggplot2::theme(legend.position = "bottom")

    ggplot2::ggsave(file.path(out_dir, "plots", "rna_protein_de_scatter.png"),
                    p, width = 14, height = 7, dpi = 300)

    # Combined TE histogram
    p_te <- ggplot2::ggplot(combined, ggplot2::aes(x = te_log2FC)) +
        ggplot2::geom_histogram(bins = 40, fill = "darkcyan", color = "white", alpha = 0.8) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::facet_wrap(~ contrast, scales = "free_y") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "Translation Efficiency Changes by Contrast",
            subtitle = "log2(TE) = log2(Protein FC) - log2(RNA FC)",
            x = "TE log2 Fold Change",
            y = "Count"
        )
    ggplot2::ggsave(file.path(out_dir, "plots", "translation_efficiency_hist.png"),
                    p_te, width = 14, height = 6)

    # Combined TE scatter
    p_te2 <- ggplot2::ggplot(combined, ggplot2::aes(x = rna_log2FC, y = protein_log2FC,
                                                     color = te_log2FC)) +
        ggplot2::geom_hline(yintercept = 0, color = "gray90") +
        ggplot2::geom_vline(xintercept = 0, color = "gray90") +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        ggplot2::geom_point(alpha = 0.7, size = 1.5) +
        ggplot2::scale_color_gradient2(
            low = "blue", mid = "grey90", high = "red",
            midpoint = 0, name = "log2(TE)"
        ) +
        ggplot2::facet_wrap(~ contrast, scales = "free") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "Translation Efficiency by Contrast",
            x = "RNA log2 Fold Change",
            y = "Protein log2 Fold Change"
        )
    ggplot2::ggsave(file.path(out_dir, "plots", "translation_efficiency_scatter.png"),
                    p_te2, width = 14, height = 7, dpi = 300)
}


#' Compute differential expression concordance between RNA and protein
#'
#' @param rna_de RNA DE table (must have feature_id and log2FC columns)
#' @param prot_de Protein DE/DA table (must have feature_id and log2FC columns)
#' @param gene_mapping Gene mapping table
#' @param out_dir Output directory (optional)
#' @return Data frame with merged DE results and concordance metrics
compute_de_concordance <- function(rna_de, prot_de, gene_mapping, out_dir = NULL,
                                   contrast_label = NULL) {

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

    # Annotate transcription factor status (dorothea)
    de_merged <- tryCatch({
        regulons <- get_tf_regulons("dorothea")
        if (!is.null(regulons)) {
            tf_names <- toupper(unique(regulons$tf))
            de_merged$is_transcription_factor <- toupper(de_merged$gene_symbol) %in% tf_names
            n_tf <- sum(de_merged$is_transcription_factor, na.rm = TRUE)
            message("  ", n_tf, "/", nrow(de_merged),
                    " RNA-protein concordance genes are known TFs")
        }
        de_merged
    }, error = function(e) {
        message("  TF annotation skipped: ", e$message)
        de_merged
    })

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
        # Directories are created here, not by the caller: reaching this point
        # means the join produced rows, so the directories will not be left empty.
        dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
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
                title = paste0("Differential Concordance: RNA vs Protein",
                               if (!is.null(contrast_label)) paste0(" (", contrast_label, ")") else ""),
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
                title = paste0("Distribution of Translation Efficiency Changes",
                               if (!is.null(contrast_label)) paste0(" (", contrast_label, ")") else ""),
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


#' Convert a signed linear fold change to log2
#'
#' The proteomics summary carries \code{linearFC} as a signed linear ratio:
#' \code{2^log2FC} when the log2 change is non-negative, and \code{-1 / 2^log2FC}
#' when it is negative (see \code{build_provenance_notes()} in
#' \code{R/core/05_export_excel.R}). A plain \code{log2()} therefore returns NaN
#' for every down-regulated feature, which silently dropped roughly half the
#' proteome from the RNA-protein correlation.
#'
#' @param x Numeric vector of signed linear fold changes.
#' @return Numeric vector of log2 fold changes, NA where \code{x} is NA or zero.
.signed_linear_fc_to_log2 <- function(x) {
    x <- as.numeric(x)
    out <- rep(NA_real_, length(x))
    ok <- !is.na(x) & x != 0
    out[ok] <- ifelse(x[ok] > 0, log2(x[ok]), -log2(abs(x[ok])))
    out
}


#' Extract a simple DE data frame from complex DE result objects
#'
#' Handles various DE result formats:
#' - RNA-seq: list with $tables (list of per-contrast data frames)
#' - Proteomics: list with $summary_df
#' - Simple data frame: passed through
#'
#' @param de_obj DE result object from single-omics pipeline
#' @return Data frame with at least feature_id, log2FC columns, or NULL
.extract_de_table <- function(de_obj) {
    if (is.null(de_obj)) return(NULL)

    # Already a data frame
    if (is.data.frame(de_obj)) {
        return(.normalize_de_columns(de_obj))
    }

    # RNA-seq format: $tables is a list of per-contrast data frames
    if (is.list(de_obj) && "tables" %in% names(de_obj) && length(de_obj$tables) > 0) {
        # Use the first contrast table
        tbl <- de_obj$tables[[1]]
        return(.normalize_de_columns(tbl))
    }

    # Proteomics format: $summary_df with contrast-specific columns
    if (is.list(de_obj) && "summary_df" %in% names(de_obj)) {
        df <- de_obj$summary_df
        # Find log2FC/linearFC columns from the first contrast
        fc_cols <- grep("^linearFC\\.", colnames(df), value = TRUE)
        padj_cols <- grep("^padj\\.", colnames(df), value = TRUE)
        id_col <- if ("FeatureID" %in% colnames(df)) "FeatureID" else if ("ID" %in% colnames(df)) "ID" else NULL

        if (length(fc_cols) > 0 && !is.null(id_col)) {
            out <- data.frame(
                feature_id = df[[id_col]],
                log2FC = .signed_linear_fc_to_log2(df[[fc_cols[1]]]),
                stringsAsFactors = FALSE
            )
            if (length(padj_cols) > 0) out$padj <- df[[padj_cols[1]]]
            return(out)
        }
    }

    # Metabolomics format: $de_tables or $summary_df
    if (is.list(de_obj) && "de_tables" %in% names(de_obj) && length(de_obj$de_tables) > 0) {
        tbl <- de_obj$de_tables[[1]]
        if (is.data.frame(tbl)) return(.normalize_de_columns(tbl))
    }

    NULL
}


#' Normalize DE column names to standard: feature_id, log2FC, padj
.normalize_de_columns <- function(df) {
    # Normalize feature ID column
    if (!"feature_id" %in% colnames(df)) {
        id_candidates <- c("FeatureID", "ID", "gene_id", "protein_id", "Row.names")
        found <- intersect(id_candidates, colnames(df))
        if (length(found) > 0) {
            df$feature_id <- df[[found[1]]]
        } else {
            df$feature_id <- rownames(df)
        }
    }

    # Normalize log2FC column
    if (!"log2FC" %in% colnames(df)) {
        fc_candidates <- c("log2FoldChange", "logFC", "log2fc")
        found <- intersect(fc_candidates, colnames(df))
        if (length(found) > 0) {
            df$log2FC <- df[[found[1]]]
        }
    }

    # Normalize padj column
    if (!"padj" %in% colnames(df)) {
        padj_candidates <- c("adj.P.Val", "FDR", "p.adjust", "qvalue")
        found <- intersect(padj_candidates, colnames(df))
        if (length(found) > 0) {
            df$padj <- df[[found[1]]]
        }
    }

    df
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
        raw_de <- if (!is.null(de_results) && nm %in% names(de_results)) de_results[[nm]] else NULL
        de_table <- .extract_de_table(raw_de)

        list(
            normalized_matrix = as.matrix(SummarizedExperiment::assay(exp_data)),
            de_table = de_table,
            da_table = de_table,
            feature_annotation = as.data.frame(SummarizedExperiment::rowData(exp_data))
        )
    })
    names(harmonized_omics) <- names(mae@ExperimentList)

    # Convert gene_protein_mapping (gene_id, protein_id) to the format expected
    # by run_rna_protein_correlation (omics, feature_id, gene_symbol)
    gene_mapping <- NULL
    if (!is.null(gene_protein_mapping) && nrow(gene_protein_mapping) > 0) {
        # gene_symbol is the JOIN KEY that compute_de_concordance() merges the two
        # layers on, so both sides must carry the same value per pair. When the
        # mapping file has no gene_symbol column (every non-model organism), the
        # RNA side used to key on gene_id and the protein side on protein_id, so
        # the merge matched nothing and DE concordance silently returned NULL.
        join_key <- if ("gene_symbol" %in% colnames(gene_protein_mapping)) {
            gene_protein_mapping$gene_symbol
        } else {
            gene_protein_mapping$gene_id
        }
        rna_rows <- data.frame(
            omics = "transcriptomics",
            feature_id = gene_protein_mapping$gene_id,
            gene_symbol = join_key,
            stringsAsFactors = FALSE
        )
        prot_rows <- data.frame(
            omics = "proteomics",
            feature_id = gene_protein_mapping$protein_id,
            gene_symbol = join_key,
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
