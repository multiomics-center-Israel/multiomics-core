#' Normalize pipeline DE results to a flat data frame for concordance
#'
#' Handles various DE result structures from RNA, proteomics, and metabolomics
#' pipelines. Returns a data frame with columns: feature_id, logFC, padj.
#'
#' @param de_res DE result object from the pipeline
#' @param omics_type Character label for error messages
#' @return data.frame with feature_id, logFC, padj columns
normalize_de_for_concordance <- function(de_res, omics_type = "unknown") {

    # Already a flat data frame with expected columns
    if (is.data.frame(de_res)) {
        if (!"feature_id" %in% names(de_res)) de_res$feature_id <- rownames(de_res)
        return(de_res)
    }

    # RNA-seq: list(dds, tables) where tables is named list of per-contrast DFs
    if (!is.null(de_res$tables) && is.list(de_res$tables)) {
        tbl <- de_res$tables[[1]]

        lfc_col <- intersect(c("log2FoldChange", "logFC"), names(tbl))[1]
        padj_col <- intersect(c("padj", "adj.P.Val", "FDR"), names(tbl))[1]
        id_col <- intersect(c("FeatureID", "gene_id", "feature_id"), names(tbl))[1]

        id_vals <- if (is.na(id_col)) rownames(tbl) else tbl[[id_col]]

        if (is.na(lfc_col)) stop("Cannot find logFC column in ", omics_type, " DE results")

        out <- data.frame(
            feature_id = id_vals,
            logFC = tbl[[lfc_col]],
            padj = if (!is.na(padj_col)) tbl[[padj_col]] else NA_real_,
            stringsAsFactors = FALSE
        )
        rownames(out) <- out$feature_id
        return(out)
    }

    # Proteomics/Metabolomics: list(summary_df, ...)
    if (!is.null(de_res$summary_df) && is.data.frame(de_res$summary_df)) {
        df <- de_res$summary_df

        id_col <- intersect(c("FeatureID", "feature_id", "protein_id"), names(df))[1]
        id_vals <- if (is.na(id_col)) rownames(df) else df[[id_col]]

        # Wide-format contrast columns (proteomics naming convention)
        lfc_cols <- grep("^(linearFC\\.imputs\\.|logFC\\.|log2FoldChange\\.)", names(df), value = TRUE)
        padj_cols <- grep("^(padj\\.imputs\\.|padj\\.|adj\\.P\\.Val\\.)", names(df), value = TRUE)

        if (length(lfc_cols) == 0) {
            stop("Cannot find logFC columns in ", omics_type, " summary_df. Available: ",
                 paste(head(names(df), 20), collapse = ", "))
        }

        lfc_vals <- df[[lfc_cols[1]]]
        # Convert linearFC (signed fold change) to log2FC
        if (grepl("^linearFC", lfc_cols[1])) {
            lfc_vals <- ifelse(lfc_vals >= 0, log2(pmax(lfc_vals, 1e-10)),
                               -log2(pmax(abs(lfc_vals), 1e-10)))
        }

        padj_vals <- if (length(padj_cols) > 0) df[[padj_cols[1]]] else NA_real_

        out <- data.frame(
            feature_id = id_vals,
            logFC = lfc_vals,
            padj = padj_vals,
            stringsAsFactors = FALSE
        )
        rownames(out) <- out$feature_id
        return(out)
    }

    stop("Unrecognized DE result format for ", omics_type,
         ". Expected list with $tables or $summary_df, or a data.frame.")
}


#' Analyze concordance between omics layers
#'
#' Compares differential expression/abundance patterns across omics to identify
#' concordant and discordant changes. Useful for validating multi-omics findings.
#'
#' @param de_results Named list of DE results per omics (must have logFC, padj)
#' @param gene_protein_mapping Mapping between gene and protein IDs
#' @param config Full config object
#' @param out_dir Output directory for plots
#' @return List with: concordance_table, plots, stats
analyze_multiomics_concordance <- function(de_results, gene_protein_mapping = NULL,
                                            config, out_dir = NULL) {

    if (length(de_results) < 2) {
        stop("Concordance analysis requires at least 2 omics layers with DE results")
    }

    message("Analyzing multi-omics concordance...")

    # Normalize all DE results to flat data frames
    de_flat <- lapply(names(de_results), function(om) {
        tryCatch(
            normalize_de_for_concordance(de_results[[om]], omics_type = om),
            error = function(e) {
                warning("Could not normalize DE results for ", om, ": ", e$message)
                NULL
            }
        )
    })
    names(de_flat) <- names(de_results)
    de_flat <- Filter(Negate(is.null), de_flat)

    if (length(de_flat) < 2) {
        stop("Concordance analysis requires at least 2 normalized omics DE results")
    }

    # For RNA + Proteomics, use gene-protein mapping
    if ("transcriptomics" %in% names(de_flat) &&
        "proteomics" %in% names(de_flat) &&
        !is.null(gene_protein_mapping)) {

        conc <- analyze_rna_protein_concordance(
            rna_de = de_flat$transcriptomics,
            prot_de = de_flat$proteomics,
            mapping = gene_protein_mapping,
            config = config,
            out_dir = out_dir
        )

        return(conc)
    }

    # General pairwise concordance for any omics combination
    omics_pairs <- combn(names(de_flat), 2, simplify = FALSE)
    concordance_list <- list()

    for (pair in omics_pairs) {
        om1 <- pair[1]
        om2 <- pair[2]

        conc <- compute_pairwise_concordance(
            de1 = de_flat[[om1]],
            de2 = de_flat[[om2]],
            name1 = om1,
            name2 = om2,
            config = config
        )

        concordance_list[[paste(om1, om2, sep = "_vs_")]] <- conc
    }

    # Generate summary plots
    plots <- list()
    if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        for (pair_name in names(concordance_list)) {
            conc <- concordance_list[[pair_name]]

            # Scatter plot of logFC values
            plot_file <- file.path(out_dir, paste0("concordance_", pair_name, ".png"))
            plots[[pair_name]] <- plot_file

            png(plot_file, width = 800, height = 800, res = 120)
            tryCatch({
                plot_concordance_scatter(conc, pair_name)
            }, error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2)
            })
            dev.off()
        }

        message("  Concordance plots saved to: ", out_dir)
    }

    list(
        concordance = concordance_list,
        plots = plots
    )
}


#' Analyze RNA-protein concordance using gene-protein mapping
analyze_rna_protein_concordance <- function(rna_de, prot_de, mapping, config, out_dir = NULL) {

    message("  Analyzing RNA-protein concordance...")

    # Ensure feature_id column
    if (!"feature_id" %in% names(rna_de)) rna_de$feature_id <- rownames(rna_de)
    if (!"feature_id" %in% names(prot_de)) prot_de$feature_id <- rownames(prot_de)

    # Merge DE results via mapping
    rna_mapped <- merge(mapping, rna_de, by.x = "gene_id", by.y = "feature_id", all = FALSE)
    merged <- merge(rna_mapped, prot_de, by.x = "protein_id", by.y = "feature_id",
                    suffixes = c("_rna", "_prot"))

    if (nrow(merged) == 0) {
        warning("No overlapping features found between RNA and protein DE results")
        return(list(concordance_table = data.frame(), stats = NULL))
    }

    # Standardize column names
    logfc_rna_col <- grep("^log.*FC.*_rna$|^logFC_rna$", names(merged), value = TRUE)[1]
    logfc_prot_col <- grep("^log.*FC.*_prot$|^logFC_prot$", names(merged), value = TRUE)[1]
    padj_rna_col <- grep("padj_rna", names(merged), value = TRUE)[1]
    padj_prot_col <- grep("padj_prot", names(merged), value = TRUE)[1]

    if (is.na(logfc_rna_col) || is.na(logfc_prot_col)) {
        stop("Cannot find logFC columns in merged DE results")
    }

    # Extract logFC values
    merged$logFC_rna <- merged[[logfc_rna_col]]
    merged$logFC_prot <- merged[[logfc_prot_col]]
    merged$padj_rna <- if (!is.na(padj_rna_col)) merged[[padj_rna_col]] else NA
    merged$padj_prot <- if (!is.na(padj_prot_col)) merged[[padj_prot_col]] else NA

    # Classify concordance
    p_cutoff <- 0.05
    merged$sig_rna <- !is.na(merged$padj_rna) & merged$padj_rna < p_cutoff
    merged$sig_prot <- !is.na(merged$padj_prot) & merged$padj_prot < p_cutoff

    merged$concordance <- classify_concordance(
        merged$logFC_rna, merged$logFC_prot,
        merged$sig_rna, merged$sig_prot
    )

    # Compute Pearson correlation
    cor_test <- cor.test(merged$logFC_rna, merged$logFC_prot, method = "pearson")
    cor_val <- cor_test$estimate
    cor_pval <- cor_test$p.value

    # Stats
    stats <- list(
        n_total = nrow(merged),
        n_both_sig = sum(merged$sig_rna & merged$sig_prot, na.rm = TRUE),
        n_concordant = sum(merged$concordance == "concordant", na.rm = TRUE),
        n_discordant = sum(merged$concordance == "discordant", na.rm = TRUE),
        correlation = cor_val,
        correlation_pval = cor_pval
    )

    message(sprintf(
        "  RNA-protein concordance: %d pairs, r=%.3f (p=%.2e), %d concordant, %d discordant",
        stats$n_total, stats$correlation, stats$correlation_pval,
        stats$n_concordant, stats$n_discordant
    ))

    # Generate plot
    plots <- list()
    if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        plots$scatter <- file.path(out_dir, "concordance_rna_protein.png")
        png(plots$scatter, width = 900, height = 800, res = 120)
        tryCatch({
            plot_rna_protein_concordance(merged, stats)
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Plot failed:", e$message), cex = 1.2)
        })
        dev.off()
    }

    list(
        concordance_table = merged,
        stats = stats,
        plots = plots
    )
}


#' Classify concordance status for paired features
classify_concordance <- function(logfc1, logfc2, sig1, sig2) {
    concordance <- rep("not_sig", length(logfc1))

    both_sig <- sig1 & sig2
    same_direction <- sign(logfc1) == sign(logfc2)

    concordance[both_sig & same_direction] <- "concordant"
    concordance[both_sig & !same_direction] <- "discordant"
    concordance[sig1 & !sig2] <- "rna_only"
    concordance[!sig1 & sig2] <- "protein_only"

    concordance
}


#' Plot RNA-protein concordance scatter
plot_rna_protein_concordance <- function(merged, stats) {
    par(mar = c(5, 5, 4, 2))

    # Color by concordance
    colors <- c(
        "concordant" = "#377EB8",
        "discordant" = "#E41A1C",
        "rna_only" = "#FFA500",
        "protein_only" = "#984EA3",
        "not_sig" = "#CCCCCC"
    )

    plot(merged$logFC_rna, merged$logFC_prot,
         col = colors[merged$concordance],
         pch = 16, cex = 0.8,
         xlab = "RNA log2 Fold Change",
         ylab = "Protein log2 Fold Change",
         main = sprintf("RNA-Protein Concordance (r=%.3f, p=%.2e)",
                        stats$correlation, stats$correlation_pval))

    abline(h = 0, v = 0, lty = 2, col = "gray")
    abline(a = 0, b = 1, lty = 2, col = "black")

    legend("topleft", legend = names(colors), col = colors, pch = 16, cex = 0.9, bty = "n")

    text(par("usr")[1] * 0.95, par("usr")[4] * 0.95,
         sprintf("n=%d\nConcordant=%d\nDiscordant=%d",
                 stats$n_total, stats$n_concordant, stats$n_discordant),
         pos = 4, cex = 0.9)
}


#' Compute pairwise concordance for general omics pairs
compute_pairwise_concordance <- function(de1, de2, name1, name2, config) {

    # Simple correlation-based concordance (no mapping needed)
    # Requires overlapping feature IDs

    if (!"feature_id" %in% names(de1)) de1$feature_id <- rownames(de1)
    if (!"feature_id" %in% names(de2)) de2$feature_id <- rownames(de2)

    common_features <- intersect(de1$feature_id, de2$feature_id)

    if (length(common_features) == 0) {
        warning(sprintf("No common features between %s and %s", name1, name2))
        return(NULL)
    }

    merged <- merge(de1, de2, by = "feature_id", suffixes = c("_1", "_2"))

    # Extract logFC
    logfc1 <- merged[[grep("^log.*FC.*_1$|^logFC_1$", names(merged), value = TRUE)[1]]]
    logfc2 <- merged[[grep("^log.*FC.*_2$|^logFC_2$", names(merged), value = TRUE)[1]]]

    cor_test <- cor.test(logfc1, logfc2, method = "pearson")

    list(
        merged = merged,
        correlation = cor_test$estimate,
        correlation_pval = cor_test$p.value,
        n_features = nrow(merged)
    )
}


#' Generic concordance scatter plot
plot_concordance_scatter <- function(conc, pair_name) {
    merged <- conc$merged
    logfc1 <- merged[[grep("^log.*FC.*_1$|^logFC_1$", names(merged), value = TRUE)[1]]]
    logfc2 <- merged[[grep("^log.*FC.*_2$|^logFC_2$", names(merged), value = TRUE)[1]]]

    plot(logfc1, logfc2, pch = 16, col = "#377EB8", cex = 0.8,
         xlab = "logFC (omics 1)", ylab = "logFC (omics 2)",
         main = sprintf("%s\nr=%.3f, p=%.2e",
                        gsub("_", " ", pair_name),
                        conc$correlation, conc$correlation_pval))

    abline(h = 0, v = 0, lty = 2, col = "gray")
    abline(a = 0, b = 1, lty = 2, col = "black")
}
