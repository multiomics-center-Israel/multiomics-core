#' MultiGSEA Correlation Plots
#'
#' Functions for visualizing cross-omics enrichment correlations
#' and generating pathway overlay plots using pathview.
#'
#' @name multigsea_plots
NULL

#' Run MultiGSEA Correlation Analysis
#'
#' Generates scatter plots comparing enrichment scores between pairs of omics.
#'
#' @param enrichment_results List containing enrichment results from run_multiomics_enrichment.
#' @param config Pipeline configuration list.
#' @param out_dir Output directory for plots and tables.
#' @return A list of ggplot objects.
#' @export
run_multigsea_plots <- function(enrichment_results, config, out_dir = NULL) {
    message("=== Running MultiGSEA Correlation Analysis ===")

    if (is.null(enrichment_results) || is.null(enrichment_results$per_omics)) {
        message("No enrichment results available for MultiGSEA.")
        return(NULL)
    }

    mg_config <- config$modes$multiomics$enrichment$multigsea %||% list()
    if (!(mg_config$run_multigsea %||% TRUE)) {
        message("MultiGSEA analysis disabled in config.")
        return(NULL)
    }

    p_thresh <- mg_config$pvalue_threshold %||% 0.05
    corr_method <- mg_config$correlation_method %||% "pearson"

    # Create output directory
    if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # Extract results per omic
    per_omics <- enrichment_results$per_omics
    omics_names <- names(per_omics)

    if (length(omics_names) < 2) {
        message("Need at least 2 omics with enrichment results for MultiGSEA.")
        return(NULL)
    }

    # Identify pairs
    pairs <- utils::combn(omics_names, 2, simplify = FALSE)
    plots <- list()

    for (pair in pairs) {
        omic1 <- pair[1]
        omic2 <- pair[2]

        res1 <- per_omics[[omic1]]
        res2 <- per_omics[[omic2]]

        if (is.null(res1) || is.null(res2)) next

        # Standardize 'term' column helper
        get_term_col <- function(df) {
            if ("term" %in% colnames(df)) {
                return(df$term)
            }
            if ("ID" %in% colnames(df)) {
                return(df$ID)
            }
            if ("Description" %in% colnames(df)) {
                return(df$Description)
            }
            return(rownames(df))
        }

        res1$term <- get_term_col(res1)
        res2$term <- get_term_col(res2)

        # Build ID → pathway name lookup
        id_to_name <- character(0)
        for (df_tmp in list(res1, res2)) {
            if ("pathway" %in% colnames(df_tmp) && "ID" %in% colnames(df_tmp)) {
                nms <- setNames(df_tmp$pathway, df_tmp$ID)
                id_to_name <- c(id_to_name, nms[!names(nms) %in% names(id_to_name)])
            }
        }

        # Union of terms
        common_terms <- union(res1$term, res2$term)

        if (length(common_terms) < 3) {
            message("Too few terms (union) between ", omic1, " and ", omic2)
            next
        }

        # Align data
        df1 <- res1[match(common_terms, res1$term), ]
        df2 <- res2[match(common_terms, res2$term), ]

        # Calculate -log10(padj) scores
        get_score <- function(df) {
            # Robustly find p-adj column
            padj_col <- NULL
            for (col in c("padj", "p.adjust", "adj.P.Val", "FDR", "qvalue", "pvalue")) {
                if (col %in% colnames(df)) {
                    padj_col <- col
                    break
                }
            }

            if (is.null(padj_col)) {
                warning("No p-value column found in enrichment results")
                return(rep(0, nrow(df)))
            }

            padj <- df[[padj_col]]
            padj[is.na(padj)] <- 1

            # Handle zero p-values
            non_zeros <- padj[padj > 0 & padj < 1]
            if (length(non_zeros) > 0) {
                min_nz <- min(non_zeros, na.rm = TRUE)
            } else {
                min_nz <- 1e-10
            }

            padj[padj == 0] <- min_nz / 10
            -log10(padj)
        }

        score1 <- get_score(df1)
        score2 <- get_score(df2)

        # Extract count (setSize) and compute fold enrichment from GeneRatio
        get_count <- function(df) {
            if ("setSize" %in% colnames(df)) return(df$setSize)
            if ("Count" %in% colnames(df)) return(df$Count)
            return(rep(NA_real_, nrow(df)))
        }

        parse_gene_ratio <- function(df) {
            if ("Fold Enrichment" %in% colnames(df)) return(df[["Fold Enrichment"]])
            if ("fold_enrichment" %in% colnames(df)) return(df$fold_enrichment)
            if ("GeneRatio" %in% colnames(df)) {
                gr <- as.character(df$GeneRatio)
                parts <- strsplit(gr, "/")
                ratio <- vapply(parts, function(p) {
                    if (length(p) == 2) as.numeric(p[1]) / as.numeric(p[2])
                    else NA_real_
                }, numeric(1))
                return(ratio)
            }
            return(rep(NA_real_, nrow(df)))
        }

        count1 <- get_count(df1)
        count2 <- get_count(df2)
        fold1 <- parse_gene_ratio(df1)
        fold2 <- parse_gene_ratio(df2)

        plot_df <- data.frame(
            term = common_terms,
            x = score1,
            y = score2,
            count1 = count1,
            count2 = count2,
            fold1 = fold1,
            fold2 = fold2,
            stringsAsFactors = FALSE
        )

        # Impute NAs (terms present in one omic but not the other)
        plot_df$count1[is.na(plot_df$count1)] <- 0
        plot_df$count2[is.na(plot_df$count2)] <- 0
        plot_df$fold1[is.na(plot_df$fold1)] <- 0
        plot_df$fold2[is.na(plot_df$fold2)] <- 0

        plot_df$avg_count <- (plot_df$count1 + plot_df$count2) / 2
        plot_df$avg_fold <- (plot_df$fold1 + plot_df$fold2) / 2

        # Resolve term to readable name, truncate long names
        resolve_term <- function(term) {
            # Look up pathway name from ID
            if (term %in% names(id_to_name) && nzchar(id_to_name[[term]])) {
                t <- id_to_name[[term]]
            } else {
                # Fallback: strip GO/KEGG ID prefixes
                t <- sub("^GO:\\d+~", "", term)
                t <- sub("^[a-z]{2,3}\\d{5}\\s*", "", t)
                if (nchar(t) == 0) t <- term
            }
            if (nchar(t) > 50) t <- paste0(substr(t, 1, 47), "...")
            t
        }
        plot_df$label <- vapply(plot_df$term, resolve_term, character(1))

        # Calculate correlation
        cor_res <- cor.test(plot_df$x, plot_df$y, method = corr_method)
        cor_val <- round(cor_res$estimate, 3)
        p_val <- signif(cor_res$p.value, 3)

        omic1_label <- gsub("_", " ", tools::toTitleCase(omic1))
        omic2_label <- gsub("_", " ", tools::toTitleCase(omic2))

        # Determine whether we have meaningful count/fold data
        has_count <- any(plot_df$avg_count > 0, na.rm = TRUE)
        has_fold <- any(plot_df$avg_fold > 0, na.rm = TRUE)

        # Create plot matching DAVID scatter style
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y))

        if (has_count && has_fold) {
            p <- p + ggplot2::geom_point(
                ggplot2::aes(size = avg_count, color = avg_fold), alpha = 0.7
            ) +
            ggplot2::scale_color_viridis_c(name = "Avg Fold Enrichment") +
            ggplot2::scale_size_continuous(name = "Avg Count")
        } else if (has_count) {
            p <- p + ggplot2::geom_point(
                ggplot2::aes(size = avg_count), color = "steelblue", alpha = 0.7
            ) +
            ggplot2::scale_size_continuous(name = "Avg Count")
        } else if (has_fold) {
            p <- p + ggplot2::geom_point(
                ggplot2::aes(color = avg_fold), size = 3, alpha = 0.7
            ) +
            ggplot2::scale_color_viridis_c(name = "Avg Fold Enrichment")
        } else {
            p <- p + ggplot2::geom_point(color = "steelblue", size = 3, alpha = 0.7)
        }

        p <- p +
            ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
            ggplot2::labs(
                title = paste0(omic1_label, " vs ", omic2_label),
                subtitle = paste0("Total Unique Terms: ", nrow(plot_df),
                                  "  |  ", corr_method, " r = ", cor_val, ", p = ", p_val),
                x = paste0(omic1_label, " [-log10(FDR)]"),
                y = paste0(omic2_label, " [-log10(FDR)]")
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
                axis.title = ggplot2::element_text(face = "bold"),
                plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                plot.subtitle = ggplot2::element_text(hjust = 0.5)
            )

        # Label terms with ggrepel — top terms from each axis + any with both
        if (requireNamespace("ggrepel", quietly = TRUE)) {
            cut_score <- -log10(p_thresh)
            both_sig <- plot_df[plot_df$x > cut_score & plot_df$y > cut_score, ]
            top_x <- head(plot_df[order(plot_df$x, decreasing = TRUE), ], 10)
            top_y <- head(plot_df[order(plot_df$y, decreasing = TRUE), ], 5)
            label_df <- unique(rbind(both_sig, top_x, top_y))
            label_df <- label_df[label_df$x > 0 | label_df$y > 0, ]

            p <- p + ggrepel::geom_text_repel(
                data = label_df,
                ggplot2::aes(label = label),
                size = 3,
                max.overlaps = Inf,
                box.padding = 0.5,
                force = 2,
                min.segment.length = 0
            )
        }

        # Save plot and data
        if (!is.null(out_dir)) {
            filename <- paste0("multigsea_", omic1, "_vs_", omic2)
            ggplot2::ggsave(
                file.path(out_dir, paste0(filename, ".png")),
                plot = p, width = 8, height = 8, dpi = 300
            )
            ggplot2::ggsave(
                file.path(out_dir, paste0(filename, ".pdf")),
                plot = p, width = 8, height = 8
            )
            write.csv(plot_df, file.path(out_dir, paste0(filename, ".csv")),
                      row.names = FALSE)
        }

        plots[[paste0(omic1, "_vs_", omic2)]] <- p
    }

    # Generate combined 2x2 plot
    if (length(plots) > 0 && !is.null(out_dir)) {
        combined <- plot_multigsea_combined(plots, per_omics, out_dir)
        if (!is.null(combined)) {
            plots[["combined"]] <- combined
        }
    }

    message("MultiGSEA plots generated: ", length(plots))
    return(plots)
}


#' Combined MultiGSEA Enrichment Plot
#'
#' Creates a 2x2 grid combining pairwise scatter plots with a summary dot plot.
#'
#' @param pairwise_plots Named list of pairwise ggplot objects.
#' @param per_omics Per-omics enrichment results list.
#' @param out_dir Output directory.
#' @return The combined ggplot/patchwork object, or NULL on failure.
plot_multigsea_combined <- function(pairwise_plots, per_omics, out_dir) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
        message("Package 'patchwork' not available. Skipping combined plot.")
        return(NULL)
    }

    # Build summary dot plot of top co-significant pathways across omics
    omics_names <- names(per_omics)
    summary_rows <- list()

    for (om in omics_names) {
        res <- per_omics[[om]]
        if (is.null(res) || !is.data.frame(res)) next

        # Find term column
        term_col <- NULL
        for (tc in c("term", "ID", "Description")) {
            if (tc %in% colnames(res)) { term_col <- tc; break }
        }
        if (is.null(term_col)) next

        # Find p-value column
        padj_col <- NULL
        for (pc in c("padj", "p.adjust", "adj.P.Val", "FDR", "qvalue", "pvalue")) {
            if (pc %in% colnames(res)) { padj_col <- pc; break }
        }
        if (is.null(padj_col)) next

        # Resolve pathway name: prefer 'pathway' column over raw ID
        term_ids <- res[[term_col]]
        if ("pathway" %in% colnames(res)) {
            term_names <- res$pathway
        } else {
            term_names <- term_ids
        }

        df_tmp <- data.frame(
            term = term_ids,
            pathway_name = term_names,
            padj = res[[padj_col]],
            omic = om,
            stringsAsFactors = FALSE
        )
        df_tmp <- df_tmp[!is.na(df_tmp$padj), ]
        summary_rows[[length(summary_rows) + 1]] <- df_tmp
    }

    if (length(summary_rows) == 0) {
        message("No data for summary panel.")
        return(NULL)
    }

    summary_df <- do.call(rbind, summary_rows)

    # Find terms significant in at least 2 omics
    sig_terms <- summary_df[summary_df$padj < 0.05, ]
    term_counts <- table(sig_terms$term)
    co_sig <- names(term_counts[term_counts >= 2])

    if (length(co_sig) == 0) {
        # Fall back to top terms by lowest p-value across any omic
        top_terms <- unique(summary_df$term[order(summary_df$padj)])[1:min(15, length(unique(summary_df$term)))]
        plot_data <- summary_df[summary_df$term %in% top_terms, ]
        panel_title <- "Top Enriched Pathways"
    } else {
        plot_data <- summary_df[summary_df$term %in% co_sig, ]
        # Keep top 15 by mean -log10(padj) across omics
        mean_scores <- tapply(-log10(pmax(plot_data$padj, 1e-300)), plot_data$term, mean, na.rm = TRUE)
        top_terms <- names(sort(mean_scores, decreasing = TRUE))[1:min(15, length(mean_scores))]
        plot_data <- plot_data[plot_data$term %in% top_terms, ]
        panel_title <- "Co-Significant Pathways"
    }

    plot_data$neg_log10_padj <- -log10(pmax(plot_data$padj, 1e-300))
    # Use pathway name, truncate long names
    plot_data$term_short <- ifelse(
        !is.na(plot_data$pathway_name) & nzchar(plot_data$pathway_name),
        plot_data$pathway_name,
        plot_data$term
    )
    plot_data$term_short <- ifelse(
        nchar(plot_data$term_short) > 40,
        paste0(substr(plot_data$term_short, 1, 37), "..."),
        plot_data$term_short
    )
    plot_data$omic_label <- gsub("_", " ", tools::toTitleCase(plot_data$omic))

    summary_panel <- ggplot2::ggplot(plot_data,
        ggplot2::aes(x = omic_label, y = stats::reorder(term_short, neg_log10_padj),
                     size = neg_log10_padj, color = neg_log10_padj)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_viridis_c(name = "-log10(padj)") +
        ggplot2::scale_size_continuous(name = "-log10(padj)", range = c(2, 8)) +
        ggplot2::labs(title = panel_title, x = NULL, y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            axis.title = ggplot2::element_text(face = "bold")
        )

    # Assemble 2x2 grid: up to 3 pairwise + summary panel
    scatter_plots <- pairwise_plots[!names(pairwise_plots) %in% "combined"]
    panels <- scatter_plots[1:min(3, length(scatter_plots))]
    panels[[length(panels) + 1]] <- summary_panel

    combined <- patchwork::wrap_plots(panels, ncol = 2)

    ggplot2::ggsave(
        file.path(out_dir, "multigsea_combined_enrichment.png"),
        plot = combined, width = 16, height = 16, dpi = 300
    )
    message("Saved combined MultiGSEA plot.")

    return(combined)
}


#' Run Pathview Visualization for Agreed Pathways
#'
#' Generates KEGG pathway overlay plots for pathways significant in multiple omics.
#'
#' @param enrichment_results List containing enrichment results.
#' @param mae_data MultiAssayExperiment object (or list with harmonized_omics).
#' @param config Pipeline configuration list.
#' @param out_dir Output directory for pathview plots.
#' @return List of generated plot paths.
#' @export
run_multigsea_pathview <- function(enrichment_results, mae_data, config, out_dir = NULL) {
    message("=== Running MultiGSEA Pathview Visualization ===")

    if (!requireNamespace("pathview", quietly = TRUE)) {
        message("Package 'pathview' not installed. Skipping.")
        return(NULL)
    }

    mg_config <- config$modes$multiomics$enrichment$multigsea %||% list()
    if (!(mg_config$run_pathview %||% TRUE)) {
        message("Pathview analysis disabled in config.")
        return(NULL)
    }

    if (is.null(enrichment_results) || is.null(enrichment_results$per_omics)) {
        message("No enrichment results available.")
        return(NULL)
    }

    # Create output directory
    if (is.null(out_dir)) {
        out_dir <- tempdir()
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Identify Agreed KEGG Pathways
    per_omics <- enrichment_results$per_omics
    omics_names <- names(per_omics)

    kegg_pathways <- list()

    for (omic in omics_names) {
        res <- per_omics[[omic]]
        if (!is.null(res) && nrow(res) > 0) {
            # Check for KEGG results
            if ("type" %in% colnames(res)) {
                kegg_res <- res[res$type == "KEGG", ]
            } else {
                # Heuristic: IDs start with "hsa" or "map" or numeric
                kegg_res <- res[grep("^hsa|^map|^[0-9]+$", res$ID), ]
            }

            if (nrow(kegg_res) > 0) {
                kegg_pathways[[omic]] <- kegg_res$ID
            }
        }
    }

    if (length(kegg_pathways) < 2) {
        message("Less than 2 omics have KEGG results. Skipping Pathview.")
        return(NULL)
    }

    # Find common pathways
    all_kegg <- unlist(kegg_pathways)
    if (length(all_kegg) == 0) {
        message("No KEGG pathways found.")
        return(NULL)
    }

    pathway_counts <- table(all_kegg)
    common_pathways <- names(pathway_counts)[pathway_counts >= 2]

    if (length(common_pathways) == 0) {
        message("No agreed KEGG pathways found between omics.")
        return(NULL)
    }

    message("Found ", length(common_pathways), " agreed KEGG pathways.")

    # Prepare Data for Pathview
    # Helper to get fold changes
    get_logfc <- function(omic_name, id_col = "entrez_id") {
        if (!omic_name %in% names(mae_data$harmonized_omics)) {
            return(NULL)
        }

        dat <- mae_data$harmonized_omics[[omic_name]]
        de <- dat$de_table %||% dat$da_table
        anno <- dat$feature_annotation

        if (is.null(de) || is.null(anno)) {
            return(NULL)
        }

        merged <- merge(de, anno, by = "feature_id")

        if (!id_col %in% colnames(merged)) {
            return(NULL)
        }

        # Get LogFC column
        fc_col <- grep("logFC|log2FoldChange", colnames(merged), ignore.case = TRUE, value = TRUE)[1]
        if (is.na(fc_col)) {
            return(NULL)
        }

        # Create named vector (take mean for duplicates)
        vec <- tapply(merged[[fc_col]], merged[[id_col]], mean, na.rm = TRUE)
        return(vec)
    }

    # Gene data (Transcriptomics + Proteomics)
    gene_data <- NULL
    rna_fc <- get_logfc("transcriptomics", "entrez_id")
    prot_fc <- get_logfc("proteomics", "entrez_id")

    if (!is.null(rna_fc) && !is.null(prot_fc)) {
        all_genes <- unique(c(names(rna_fc), names(prot_fc)))
        gene_data <- matrix(NA,
            nrow = length(all_genes), ncol = 2,
            dimnames = list(all_genes, c("Transcriptomics", "Proteomics"))
        )
        idx_rna <- match(names(rna_fc), all_genes)
        gene_data[idx_rna, 1] <- rna_fc
        idx_prot <- match(names(prot_fc), all_genes)
        gene_data[idx_prot, 2] <- prot_fc
    } else if (!is.null(rna_fc)) {
        gene_data <- rna_fc
    } else if (!is.null(prot_fc)) {
        gene_data <- prot_fc
    }

    # Metabolomics Data (CPD Data)
    cpd_data <- NULL
    met_fc <- get_logfc("metabolomics", "kegg_id")
    if (!is.null(met_fc)) {
        cpd_data <- met_fc
    }

    if (is.null(gene_data) && is.null(cpd_data)) {
        message("No valid Entrez/KEGG IDs found in data for Pathview.")
        return(NULL)
    }

    # Run Pathview
    generated_plots <- list()

    cwd <- getwd()
    setwd(out_dir)
    on.exit(setwd(cwd))

    # Determine organism code
    organism <- config$global$organism %||% "human"
    if (organism == "c_elegans") {
        species_code <- "cel"
    } else {
        species_code <- mg_config$organism_code %||% "hsa"
    }

    for (pid in common_pathways) {
        # Clean ID
        clean_pid <- sub("^[a-z]+", "", pid)

        tryCatch({
            pv_out <- pathview::pathview(
                gene.data = gene_data,
                cpd.data = cpd_data,
                pathway.id = clean_pid,
                species = species_code,
                out.suffix = "multiomics",
                temp.file = TRUE,
                kegg.dir = out_dir,
                keys.align = "y",
                match.data = TRUE,
                multi.state = !is.null(dim(gene_data)) && ncol(gene_data) > 1,
                same.layer = FALSE
            )

            outfile <- paste0(species_code, clean_pid, ".multiomics.png")
            if (file.exists(outfile)) {
                generated_plots[[pid]] <- file.path(out_dir, outfile)
                message("Generated Pathview: ", outfile)
            }
        }, error = function(e) {
            message("Pathview failed for ", pid, ": ", e$message)
        })
    }

    return(generated_plots)
}
