#' RNA-seq Pathway / Enrichment Analysis
#'
#' Runs fGSEA and/or ORA per contrast using GO / KEGG / custom GMT gene sets.
#' Adapted for multiomics-core DE result format (FeatureID, log2FoldChange, pvalue, padj).

# ==============================================================================
# GO TERM NAME LOOKUP
# ==============================================================================

#' Look up GO term names from GO IDs
#'
#' @param go_ids Character vector of GO IDs (e.g., "GO:0000001")
#' @return Named vector with GO IDs as names and term names as values
lookup_go_term_names <- function(go_ids) {
    if (length(go_ids) == 0) return(character(0))

    # Initialize with IDs as fallback
    term_names <- setNames(go_ids, go_ids)

    # Try GO.db package first
    if (requireNamespace("GO.db", quietly = TRUE)) {
        tryCatch({
            # Filter to valid GO IDs present in GO.db (some may be obsolete)
            valid_keys <- AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID")
            valid_ids <- go_ids[go_ids %in% valid_keys]
            if (length(valid_ids) > 0) {
                go_terms <- AnnotationDbi::select(
                    GO.db::GO.db,
                    keys = valid_ids,
                    columns = "TERM",
                    keytype = "GOID"
                )
                # Update names where we got results
                matched <- match(go_ids, go_terms$GOID)
                has_term <- !is.na(matched) & !is.na(go_terms$TERM[matched])
                term_names[has_term] <- go_terms$TERM[matched[has_term]]
            }
        }, error = function(e) {
            message("Could not look up GO terms: ", e$message)
        })
    }

    term_names
}

#' Add pathway names to fGSEA/ORA results
#'
#' @param pathway_df Data frame with pathway analysis results (has 'pathway' column)
#' @param database Character string indicating the database ("GO", "KEGG", etc.)
#' @return Data frame with added 'pathway_name' column
add_pathway_names <- function(pathway_df, database, gene_sets = NULL) {
    if (is.null(pathway_df) || nrow(pathway_df) == 0) return(pathway_df)
    if (!"pathway" %in% colnames(pathway_df)) return(pathway_df)

    pathway_ids <- pathway_df$pathway

    if (database == "GO" || grepl("^GO", database, ignore.case = TRUE)) {
        # Look up GO term names
        names_vec <- lookup_go_term_names(pathway_ids)
        pathway_df$pathway_name <- unname(names_vec[pathway_ids])
    } else if (database == "KEGG" || grepl("KEGG", database, ignore.case = TRUE)) {
        # For KEGG, the pathway names are usually already descriptive
        # but we'll keep the ID for now as KEGG lookup is more complex
        pathway_df$pathway_name <- pathway_ids
    } else {
        # For custom GMT, use descriptions attribute if available
        descriptions <- if (!is.null(gene_sets)) attr(gene_sets, "descriptions") else NULL
        if (!is.null(descriptions)) {
            pathway_df$pathway_name <- unname(descriptions[pathway_ids])
        } else {
            pathway_df$pathway_name <- pathway_ids
        }
        # For entries where description is still a GO ID, look up the term name
        is_go_id <- grepl("^GO:[0-9]+$", pathway_df$pathway_name)
        if (any(is_go_id)) {
            go_names <- lookup_go_term_names(pathway_df$pathway_name[is_go_id])
            pathway_df$pathway_name[is_go_id] <- unname(go_names)
        }
    }

    # Reorder columns to put pathway_name right after pathway
    cols <- colnames(pathway_df)
    pathway_idx <- which(cols == "pathway")
    new_order <- c(
        cols[1:pathway_idx],
        "pathway_name",
        setdiff(cols[(pathway_idx + 1):length(cols)], "pathway_name")
    )
    pathway_df <- pathway_df[, new_order, drop = FALSE]

    pathway_df
}

# ==============================================================================
# PATHWAY ANALYSIS (fGSEA)
# ==============================================================================

#' Run pathway analysis on DE results
#'
#' @param de_tables Named list of DE result data frames (from run_deseq2_de()$tables).
#'   Each table has columns: FeatureID, log2FoldChange, pvalue, padj, stat, etc.
#' @param gene_sets Named list of gene set collections from load_gene_sets()
#' @param annotation Gene annotation data frame (gene_id, symbol, entrez_id)
#' @param method "fgsea", "ora", or "both"
#' @param min_size Minimum gene set size for fGSEA
#' @param max_size Maximum gene set size for fGSEA
#' @return Named list (by contrast) of named lists (by db+method) of result data frames
#' @export
run_pathway_analysis <- function(de_tables,
                                  gene_sets,
                                  annotation = NULL,
                                  method = "fgsea",
                                  min_size = 10,
                                  max_size = 500) {

    if (length(gene_sets) == 0) {
        message("No gene sets available. Skipping pathway analysis.")
        return(list())
    }

    pathway_results <- list()

    for (contrast_name in names(de_tables)) {
        res <- de_tables[[contrast_name]]
        if (is.null(res) || nrow(res) == 0) next

        message("Running pathway analysis for: ", contrast_name)
        contrast_results <- list()

        for (db_name in names(gene_sets)) {
            gs <- gene_sets[[db_name]]
            if (length(gs) == 0) next

            message("  Database: ", db_name, " (", length(gs), " gene sets)")

            tryCatch({
                # ---- fGSEA ----
                if (method %in% c("fgsea", "both")) {

                    # Build ranked gene list from DE table
                    # Prefer stat column (Wald statistic); fallback to sign(lfc)*-log10(p)
                    if ("stat" %in% colnames(res)) {
                        ranks <- setNames(res$stat, res$FeatureID)
                    } else {
                        ranks <- setNames(
                            sign(res$log2FoldChange) * -log10(res$pvalue + 1e-300),
                            res$FeatureID
                        )
                    }

                    ranks <- ranks[!is.na(ranks)]
                    ranks <- sort(ranks, decreasing = TRUE)

                    fgsea_res <- fgsea::fgsea(
                        pathways = gs,
                        stats = ranks,
                        minSize = min_size,
                        maxSize = max_size,
                        nPermSimple = 10000
                    )

                    fgsea_df <- as.data.frame(fgsea_res)

                    if (nrow(fgsea_df) == 0) {
                        message("    fgsea: no gene set overlap — skipping ", db_name)
                    } else {
                        fgsea_df$contrast <- contrast_name
                        fgsea_df$database <- db_name
                        fgsea_df$method <- "fgsea"

                        # Collapse leading edge to comma-separated string
                        fgsea_df$leadingEdge <- sapply(fgsea_df$leadingEdge, paste, collapse = ",")

                        # Add pathway names (GO term names, etc.)
                        fgsea_df <- add_pathway_names(fgsea_df, db_name, gs)

                        contrast_results[[paste0(db_name, "_fgsea")]] <- fgsea_df

                        n_sig <- sum(fgsea_df$padj < 0.05, na.rm = TRUE)
                        message("    fgsea: ", n_sig, " significant pathways (padj < 0.05)")
                    }
                }

                # ---- ORA ----
                if (method %in% c("ora", "both")) {

                    # Identify significant up/down genes
                    de_cfg_padj <- 0.05
                    de_cfg_lfc  <- log2(1.5)
                    sig_up   <- res$FeatureID[!is.na(res$padj) &
                                              res$padj < de_cfg_padj &
                                              res$log2FoldChange > de_cfg_lfc]
                    sig_down <- res$FeatureID[!is.na(res$padj) &
                                              res$padj < de_cfg_padj &
                                              res$log2FoldChange < -de_cfg_lfc]
                    background <- res$FeatureID

                    if (length(sig_up) >= 5) {
                        ora_up <- run_ora(sig_up, gs, background, min_size, max_size)
                        if (nrow(ora_up) > 0) {
                            ora_up$contrast <- contrast_name
                            ora_up$database <- db_name
                            ora_up$method <- "ora"
                            ora_up$direction <- "up"
                            ora_up <- add_pathway_names(ora_up, db_name, gs)
                            contrast_results[[paste0(db_name, "_ora_up")]] <- ora_up
                        }
                    }

                    if (length(sig_down) >= 5) {
                        ora_down <- run_ora(sig_down, gs, background, min_size, max_size)
                        if (nrow(ora_down) > 0) {
                            ora_down$contrast <- contrast_name
                            ora_down$database <- db_name
                            ora_down$method <- "ora"
                            ora_down$direction <- "down"
                            ora_down <- add_pathway_names(ora_down, db_name, gs)
                            contrast_results[[paste0(db_name, "_ora_down")]] <- ora_down
                        }
                    }
                }
            }, error = function(e) {
                warning("Pathway analysis failed for ", db_name, ": ", e$message)
            })
        }

        pathway_results[[contrast_name]] <- contrast_results
    }

    pathway_results
}

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

#' Save pathway results to CSV files
#'
#' @param pathway_results Output of run_pathway_analysis()
#' @param output_dir Directory to write CSV files
#' @return Character vector of output file paths
#' @export
save_pathway_results <- function(pathway_results, output_dir) {

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_files <- c()

    for (contrast_name in names(pathway_results)) {
        contrast_res <- pathway_results[[contrast_name]]

        for (analysis_name in names(contrast_res)) {
            res_df <- contrast_res[[analysis_name]]
            if (is.null(res_df) || nrow(res_df) == 0) next

            clean_contrast <- gsub("[^a-zA-Z0-9_-]", "_", contrast_name)
            clean_analysis <- gsub("[^a-zA-Z0-9_-]", "_", analysis_name)

            output_file <- file.path(output_dir,
                                     paste0("pathway_", clean_contrast, "_", clean_analysis, ".csv"))
            write.csv(res_df, output_file, row.names = FALSE)
            output_files <- c(output_files, output_file)
        }
    }

    message("Saved ", length(output_files), " pathway result file(s) to ", output_dir)
    output_files
}

# ==============================================================================
# PATHWAY PLOTS
# ==============================================================================

#' Generate pathway visualization dotplots
#'
#' @param pathway_results Output of run_pathway_analysis()
#' @param output_dir Directory for PNG files
#' @return Named list of plot file paths
#' @export
generate_pathway_plots <- function(pathway_results, output_dir) {

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    plot_files <- list()

    for (contrast_name in names(pathway_results)) {
        contrast_res <- pathway_results[[contrast_name]]
        clean_contrast <- gsub("[^a-zA-Z0-9_-]", "_", contrast_name)

        for (analysis_name in names(contrast_res)) {
            res_df <- contrast_res[[analysis_name]]
            if (is.null(res_df) || nrow(res_df) == 0) next

            clean_analysis <- gsub("[^a-zA-Z0-9_-]", "_", analysis_name)

            if ("padj" %in% colnames(res_df)) {
                top_pathways <- head(res_df[order(res_df$padj), ], 20)
            } else if ("pvalue" %in% colnames(res_df)) {
                top_pathways <- head(res_df[order(res_df$pvalue), ], 20)
            } else {
                next
            }

            if (nrow(top_pathways) < 3) next

            # Use pathway_name if available, otherwise fall back to pathway ID
            label_col <- if ("pathway_name" %in% colnames(top_pathways)) "pathway_name" else "pathway"
            # Truncate long names for display (max 60 characters)
            top_pathways$pathway_label <- substr(top_pathways[[label_col]], 1, 60)

            if ("NES" %in% colnames(top_pathways)) {
                p <- ggplot2::ggplot(top_pathways,
                                     ggplot2::aes(x = NES, y = reorder(pathway_label, NES))) +
                    ggplot2::geom_point(ggplot2::aes(size = size, color = -log10(padj))) +
                    ggplot2::scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
                    ggplot2::labs(
                        title = paste("Top Pathways:", contrast_name),
                        subtitle = analysis_name,
                        x = "Normalized Enrichment Score",
                        y = "",
                        size = "Gene Set Size"
                    ) +
                    ggplot2::theme_minimal() +
                    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
            } else {
                p <- ggplot2::ggplot(top_pathways,
                                     ggplot2::aes(x = -log10(pvalue),
                                                  y = reorder(pathway_label, -log10(pvalue)))) +
                    ggplot2::geom_point(ggplot2::aes(size = overlap, color = -log10(padj))) +
                    ggplot2::scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
                    ggplot2::labs(
                        title = paste("Top Pathways:", contrast_name),
                        subtitle = analysis_name,
                        x = "-log10(p-value)",
                        y = "",
                        size = "Overlap"
                    ) +
                    ggplot2::theme_minimal() +
                    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
            }

            plot_file <- file.path(output_dir,
                                   paste0("pathway_", clean_contrast, "_", clean_analysis, ".png"))
            ggplot2::ggsave(plot_file, p, width = 10, height = 8)
            plot_files[[paste0(clean_contrast, "_", clean_analysis)]] <- plot_file
        }
    }

    message("Generated ", length(plot_files), " pathway plots in ", output_dir)
    plot_files
}

# ==============================================================================
# PATHWAY VOLCANO DATA
# ==============================================================================

#' Build pathway-colored volcano data
#'
#' Creates a data frame tagging each gene with its enriched pathway memberships,
#' where genes are colored by pathway membership on a volcano plot.
#'
#' @param de_table DE table with FeatureID, log2FoldChange, pvalue, padj columns
#' @param pathway_results Pathway results from run_pathway_analysis()
#' @return data.frame: FeatureID, log2FC, neg_log10_pvalue, padj, pathways
build_pathway_volcano_data <- function(de_table, pathway_results) {
    if (is.null(de_table) || nrow(de_table) == 0) return(NULL)
    if (is.null(pathway_results)) return(NULL)

    # Cap -log10(p) at the 99.5th percentile to remove extreme outliers
    # that would compress the bulk of the data to the bottom of the plot
    neg_log10_p <- -log10(de_table$pvalue + 1e-300)
    cap <- quantile(neg_log10_p[is.finite(neg_log10_p)], 0.995, na.rm = TRUE)
    cap <- max(cap, 10)  # minimum cap of 10 so small datasets aren't over-truncated
    neg_log10_p <- pmin(neg_log10_p, cap)

    volcano_df <- data.frame(
        FeatureID = de_table$FeatureID,
        log2FC = de_table$log2FoldChange,
        neg_log10_pvalue = neg_log10_p,
        padj = de_table$padj,
        stringsAsFactors = FALSE
    )

    # Collect enriched pathways and their member genes
    pathway_memberships <- character(nrow(volcano_df))

    # Flatten: collect all data.frames from the nested structure
    all_result_dfs <- list()
    for (top_name in names(pathway_results)) {
        top_elem <- pathway_results[[top_name]]
        if (is.data.frame(top_elem)) {
            all_result_dfs <- c(all_result_dfs, list(top_elem))
        } else if (is.list(top_elem)) {
            for (sub_name in names(top_elem)) {
                sub_elem <- top_elem[[sub_name]]
                if (is.data.frame(sub_elem)) {
                    all_result_dfs <- c(all_result_dfs, list(sub_elem))
                }
            }
        }
    }

    enriched_sets <- list()
    for (res in all_result_dfs) {
        if ("leadingEdge" %in% colnames(res)) {
            sig <- res[!is.na(res$padj) & res$padj < 0.05, , drop = FALSE]
            for (r in seq_len(nrow(sig))) {
                # Prefer human-readable pathway_name over GO ID
                pw_label <- if ("pathway_name" %in% colnames(sig) && nzchar(sig$pathway_name[r] %||% "")) {
                    sig$pathway_name[r]
                } else {
                    sig$pathway[r]
                }
                genes <- trimws(strsplit(as.character(sig$leadingEdge[[r]]), ",")[[1]])
                if (length(genes) > 0) enriched_sets[[pw_label]] <- genes
            }
        } else if ("geneID" %in% colnames(res)) {
            sig <- res[!is.na(res$p.adjust) & res$p.adjust < 0.05, , drop = FALSE]
            for (r in seq_len(nrow(sig))) {
                pw_label <- sig$Description[r] %||% sig$ID[r]
                genes <- unlist(strsplit(as.character(sig$geneID[r]), "/"))
                if (length(genes) > 0) enriched_sets[[pw_label]] <- genes
            }
        }
    }

    if (length(enriched_sets) == 0) return(volcano_df)

    # Tag each gene with its pathway memberships (top 5 per gene)
    for (i in seq_len(nrow(volcano_df))) {
        fid <- volcano_df$FeatureID[i]
        member_of <- names(enriched_sets)[vapply(enriched_sets, function(gs) fid %in% gs, logical(1))]
        if (length(member_of) > 5) member_of <- member_of[1:5]
        pathway_memberships[i] <- paste(member_of, collapse = "; ")
    }

    volcano_df$pathways <- pathway_memberships
    volcano_df
}
