#' Cross-omics pathway enrichment analysis
#'
#' Combines enrichment results from multiple omics layers to identify
#' pathways that are consistently dysregulated across modalities.
#'
#' @param enrichment_results Named list of enrichment results per omics
#' @param config Full config object
#' @param out_dir Output directory for plots
#' @return List with: combined_pathways, meta_analysis, plots
analyze_cross_omics_enrichment <- function(enrichment_results, config, out_dir = NULL) {

    if (length(enrichment_results) < 2) {
        message("Cross-omics enrichment requires ≥2 omics layers with enrichment results")
        return(NULL)
    }

    message("Analyzing cross-omics pathway enrichment...")

    omics <- names(enrichment_results)

    # Extract pathway-level results from each omics
    pathway_tables <- list()
    for (om in omics) {
        enrich_res <- enrichment_results[[om]]

        # Handle different enrichment result structures
        if (!is.null(enrich_res$enrichment_df)) {
            pathway_tables[[om]] <- enrich_res$enrichment_df
        } else if (is.data.frame(enrich_res)) {
            pathway_tables[[om]] <- enrich_res
        } else {
            warning("Cannot extract pathway table from ", om, " enrichment results")
            next
        }
    }

    if (length(pathway_tables) < 2) {
        warning("Insufficient pathway tables for cross-omics enrichment")
        return(NULL)
    }

    # Find common pathways across omics
    all_pathways <- lapply(pathway_tables, function(df) {
        if ("pathway" %in% names(df)) df$pathway
        else if ("ID" %in% names(df)) df$ID
        else if ("Description" %in% names(df)) df$Description
        else rownames(df)
    })

    common_pathways <- Reduce(intersect, all_pathways)

    if (length(common_pathways) == 0) {
        message("No common pathways found across omics layers")
        return(NULL)
    }

    message(sprintf("  Found %d common pathways across %d omics layers",
                    length(common_pathways), length(omics)))

    # Merge pathway p-values for meta-analysis
    merged_pathways <- merge_pathway_pvalues(pathway_tables, common_pathways, omics)

    # Combine p-values using Fisher's method
    meta_results <- fisher_combined_pvalues(merged_pathways)

    # Sort by combined p-value
    meta_results <- meta_results[order(meta_results$combined_pval), ]

    # Generate plots
    plots <- list()
    if (!is.null(out_dir) && nrow(meta_results) > 0) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        # Heatmap of pathway significance across omics
        plots$pathway_heatmap <- file.path(out_dir, "cross_omics_pathway_heatmap.png")
        png(plots$pathway_heatmap, width = 1000, height = 800, res = 120)
        tryCatch({
            plot_cross_omics_pathway_heatmap(meta_results, omics)
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Heatmap failed:", e$message), cex = 1.2)
        })
        dev.off()

        message("  Cross-omics enrichment plots saved to: ", out_dir)
    }

    list(
        common_pathways = common_pathways,
        meta_analysis = meta_results,
        pathway_tables = pathway_tables,
        plots = plots
    )
}


#' Merge pathway p-values from multiple omics
merge_pathway_pvalues <- function(pathway_tables, common_pathways, omics) {

    merged <- data.frame(pathway = common_pathways, stringsAsFactors = FALSE)

    for (om in omics) {
        df <- pathway_tables[[om]]

        # Identify pathway column
        pathway_col <- if ("pathway" %in% names(df)) "pathway"
                       else if ("ID" %in% names(df)) "ID"
                       else if ("Description" %in% names(df)) "Description"
                       else NULL

        if (is.null(pathway_col)) {
            warning("Cannot identify pathway column in ", om, " enrichment table")
            next
        }

        # Identify p-value column
        pval_col <- if ("pvalue" %in% names(df)) "pvalue"
                    else if ("pval" %in% names(df)) "pval"
                    else if ("p.adjust" %in% names(df)) "p.adjust"
                    else if ("padj" %in% names(df)) "padj"
                    else NULL

        if (is.null(pval_col)) {
            warning("Cannot identify p-value column in ", om, " enrichment table")
            next
        }

        # Subset to common pathways
        df_sub <- df[df[[pathway_col]] %in% common_pathways, c(pathway_col, pval_col)]
        colnames(df_sub) <- c("pathway", paste0("pval_", om))

        # Merge
        merged <- merge(merged, df_sub, by = "pathway", all.x = TRUE)
    }

    merged
}


#' Combine p-values using Fisher's method
fisher_combined_pvalues <- function(merged_pathways) {

    pval_cols <- grep("^pval_", names(merged_pathways), value = TRUE)

    if (length(pval_cols) < 2) {
        stop("Need at least 2 p-value columns for Fisher's method")
    }

    pval_matrix <- as.matrix(merged_pathways[, pval_cols])

    # Fisher's method: -2 * sum(log(p_i)) ~ chi-squared(2k)
    combined_pvals <- apply(pval_matrix, 1, function(pvals) {
        pvals <- pvals[!is.na(pvals) & pvals > 0]
        if (length(pvals) == 0) return(NA)

        chi_stat <- -2 * sum(log(pvals))
        df <- 2 * length(pvals)
        pchisq(chi_stat, df = df, lower.tail = FALSE)
    })

    merged_pathways$combined_pval <- combined_pvals
    merged_pathways$combined_padj <- p.adjust(combined_pvals, method = "BH")

    merged_pathways
}


#' Plot cross-omics pathway heatmap
plot_cross_omics_pathway_heatmap <- function(meta_results, omics, top_n = 30) {

    # Select top N pathways by combined p-value
    top_pathways <- meta_results[seq_len(min(top_n, nrow(meta_results))), ]

    pval_cols <- grep("^pval_", names(top_pathways), value = TRUE)
    pval_matrix <- as.matrix(top_pathways[, pval_cols])
    rownames(pval_matrix) <- top_pathways$pathway

    # Transform to -log10(p)
    log_pval_matrix <- -log10(pval_matrix + 1e-300)
    colnames(log_pval_matrix) <- gsub("^pval_", "", colnames(log_pval_matrix))

    # Heatmap
    if (requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(log_pval_matrix,
                           cluster_rows = TRUE, cluster_cols = FALSE,
                           main = "Cross-Omics Pathway Enrichment",
                           color = colorRampPalette(c("white", "orange", "red"))(50),
                           fontsize_row = 8, fontsize_col = 10,
                           angle_col = 45)
    } else {
        # Fallback: simple heatmap
        heatmap(log_pval_matrix, scale = "none", Colv = NA,
                main = "Cross-Omics Pathway Enrichment",
                col = colorRampPalette(c("white", "orange", "red"))(50))
    }
}


#' Write cross-omics enrichment results
write_cross_omics_enrichment <- function(enrichment_res, out_dir) {

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    if (!is.null(enrichment_res$meta_analysis)) {
        write.csv(enrichment_res$meta_analysis,
                  file.path(out_dir, "cross_omics_pathways_meta_analysis.csv"),
                  row.names = FALSE)
    }

    message("Cross-omics enrichment results written to: ", out_dir)
    invisible(NULL)
}
