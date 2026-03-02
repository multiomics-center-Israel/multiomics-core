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

        res1 <- per_omics[[omic1]]$results
        res2 <- per_omics[[omic2]]$results

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

        plot_df <- data.frame(
            term = common_terms,
            x = score1,
            y = score2,
            stringsAsFactors = FALSE
        )

        # Determine significance status for coloring
        cut_score <- -log10(p_thresh)

        plot_df$status <- "Not Sig"
        plot_df$status[plot_df$x > cut_score & plot_df$y > cut_score] <- "Both Sig"
        plot_df$status[plot_df$x > cut_score & plot_df$y <= cut_score] <- paste0(omic1, " Sig")
        plot_df$status[plot_df$x <= cut_score & plot_df$y > cut_score] <- paste0(omic2, " Sig")

        # Calculate correlation
        cor_res <- cor.test(plot_df$x, plot_df$y, method = corr_method)
        cor_val <- round(cor_res$estimate, 3)
        p_val <- signif(cor_res$p.value, 3)

        # Create plot
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, color = status, label = term)) +
            ggplot2::geom_point(alpha = 0.7, size = 2) +
            ggplot2::geom_hline(yintercept = cut_score, linetype = "dashed", color = "gray50") +
            ggplot2::geom_vline(xintercept = cut_score, linetype = "dashed", color = "gray50") +
            ggplot2::labs(
                title = paste0("MultiGSEA: ", omic1, " vs ", omic2),
                subtitle = paste0(corr_method, " cor = ", cor_val, ", p = ", p_val),
                x = paste0(omic1, " -log10(FDR)"),
                y = paste0(omic2, " -log10(FDR)")
            ) +
            ggplot2::theme_minimal() +
            ggplot2::scale_color_manual(values = c(
                "Both Sig" = "red",
                "Not Sig" = "grey",
                setNames("steelblue", paste0(omic1, " Sig")),
                setNames("darkgreen", paste0(omic2, " Sig"))
            ))

        # Add labels for top points
        if (requireNamespace("ggrepel", quietly = TRUE)) {
            label_df <- plot_df[plot_df$status == "Both Sig", ]
            if (nrow(label_df) > 0) {
                label_df$sum_score <- label_df$x + label_df$y
                label_df <- label_df[order(label_df$sum_score, decreasing = TRUE), ]
                label_df <- head(label_df, 10)

                p <- p + ggrepel::geom_text_repel(
                    data = label_df,
                    ggplot2::aes(label = term),
                    size = 3,
                    max.overlaps = 10
                )
            }
        }

        # Save plot and data
        if (!is.null(out_dir)) {
            filename <- paste0("multigsea_", omic1, "_vs_", omic2)
            ggplot2::ggsave(
                file.path(out_dir, paste0(filename, ".png")),
                plot = p, width = 8, height = 8, dpi = 300
            )
            write.csv(plot_df, file.path(out_dir, paste0(filename, ".csv")),
                      row.names = FALSE)
        }

        plots[[paste0(omic1, "_vs_", omic2)]] <- p
    }

    message("MultiGSEA plots generated: ", length(plots))
    return(plots)
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
        res <- per_omics[[omic]]$results
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
