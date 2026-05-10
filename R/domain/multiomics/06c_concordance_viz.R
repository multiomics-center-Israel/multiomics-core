#' Concordance visualization plots for RNA-protein comparisons
#'
#' Generates lollipop, pie, stacked bar, volcano overlay, Venn,
#' correlation scatter, and mHG enrichment test plots from
#' concordance analysis results.
#'
#' @param conc_df Concordance data frame with logFC_rna, logFC_prot,
#'   padj_rna, padj_prot, concordance columns
#' @param out_dir Output directory for PNG files
#' @return List of saved file paths (invisible)
plot_concordance_visualization <- function(conc_df, out_dir) {

    if (is.null(conc_df) || nrow(conc_df) == 0) return(invisible(NULL))

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saved <- character(0)

    # Derive flags
    conc_df$sig_prot <- !is.na(conc_df$padj_prot) & conc_df$padj_prot < 0.05
    conc_df$sig_rna  <- !is.na(conc_df$padj_rna) & conc_df$padj_rna < 0.05
    conc_df$prot_trend <- ifelse(!conc_df$sig_prot, "Not sig",
                          ifelse(conc_df$logFC_prot > 0, "Up", "Down"))
    conc_df$rna_change <- ifelse(!conc_df$sig_rna, "Not sig",
                          ifelse(conc_df$logFC_rna > 0, "Up", "Down"))

    # Gene label: prefer gene_symbol, then gene_id, then protein_id
    gene_col <- intersect(c("gene_symbol", "gene_id", "protein_id"), colnames(conc_df))[1]
    conc_df$gene_label <- if (!is.na(gene_col)) conc_df[[gene_col]] else rownames(conc_df)

    # --- Lollipop plots ---
    saved <- c(saved, .plot_lollipops(conc_df, out_dir))

    # --- Pie charts ---
    saved <- c(saved, .plot_pie_charts(conc_df, out_dir))

    # --- Stacked bar ---
    saved <- c(saved, .plot_stacked_bar(conc_df, out_dir))

    # --- Volcano overlay ---
    saved <- c(saved, .plot_volcano_overlay(conc_df, out_dir))

    # --- Venn diagram ---
    saved <- c(saved, .plot_venn(conc_df, out_dir))

    # --- Correlation scatter ---
    saved <- c(saved, .plot_correlation_scatter(conc_df, out_dir))

    # --- mHG test ---
    saved <- c(saved, .plot_mhg_test(conc_df, out_dir))

    message("  Concordance visualization: ", length(saved), " plots saved to ", out_dir)
    invisible(saved)
}


# --- Internal helpers ---

.plot_lollipops <- function(conc_df, out_dir) {
    saved <- character(0)
    sig_prot <- conc_df[conc_df$sig_prot, ]
    if (nrow(sig_prot) < 2) return(saved)

    # Cap at 100 for readability
    if (nrow(sig_prot) > 100) {
        sig_prot <- sig_prot[order(abs(sig_prot$logFC_prot), decreasing = TRUE), ]
        sig_prot <- head(sig_prot, 100)
    }

    for (subset_info in list(
        list(name = "all", data = sig_prot, title = "All Significant Proteins"),
        list(name = "up", data = sig_prot[sig_prot$logFC_prot > 0, ], title = "Upregulated Proteins"),
        list(name = "down", data = sig_prot[sig_prot$logFC_prot < 0, ], title = "Downregulated Proteins")
    )) {
        df <- subset_info$data
        if (nrow(df) < 2) next

        df$gene_f <- forcats::fct_reorder(df$gene_label, df$logFC_prot)

        p_top <- ggplot2::ggplot(df, ggplot2::aes(x = gene_f, y = logFC_prot, color = prot_trend)) +
            ggplot2::geom_segment(ggplot2::aes(xend = gene_f, y = 0, yend = logFC_prot), linewidth = 0.5) +
            ggplot2::geom_point(size = 1.5) +
            ggplot2::scale_color_manual(values = c("Up" = "springgreen4", "Down" = "violetred", "Not sig" = "gray")) +
            ggplot2::labs(y = "Protein log2FC", title = paste(subset_info$title, "— Protein")) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(), legend.position = "none",
                           panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())

        p_bot <- ggplot2::ggplot(df, ggplot2::aes(x = gene_f, y = logFC_rna, color = rna_change)) +
            ggplot2::geom_segment(ggplot2::aes(xend = gene_f, y = 0, yend = logFC_rna), linewidth = 0.5) +
            ggplot2::geom_point(size = 1.5) +
            ggplot2::scale_color_manual(values = c("Up" = "lightgreen", "Down" = "lightcoral", "Not sig" = "gray")) +
            ggplot2::labs(y = "RNA log2FC", title = "Corresponding RNA") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 4),
                           legend.position = "none",
                           panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())

        p <- p_top / p_bot
        fname <- file.path(out_dir, paste0("lollipop_", subset_info$name, ".png"))
        ggplot2::ggsave(fname, p, width = 14, height = 8, dpi = 200)
        saved <- c(saved, fname)
    }
    saved
}


.plot_pie_charts <- function(conc_df, out_dir) {
    saved <- character(0)
    for (direction in c("up", "down")) {
        subset <- if (direction == "up") {
            conc_df[conc_df$sig_prot & conc_df$logFC_prot > 0, ]
        } else {
            conc_df[conc_df$sig_prot & conc_df$logFC_prot < 0, ]
        }
        if (nrow(subset) < 2) next

        counts <- data.frame(
            RNA = c("Up", "Down", "Not sig"),
            n = c(sum(subset$rna_change == "Up"), sum(subset$rna_change == "Down"),
                  sum(subset$rna_change == "Not sig")),
            stringsAsFactors = FALSE
        )
        counts$pct <- round(100 * counts$n / sum(counts$n), 1)
        counts$label <- paste0(counts$RNA, "\n", counts$n, " (", counts$pct, "%)")
        counts$RNA <- factor(counts$RNA, levels = c("Up", "Down", "Not sig"))

        p <- ggplot2::ggplot(counts, ggplot2::aes(x = 2, y = n, fill = RNA)) +
            ggplot2::geom_col(width = 1) +
            ggplot2::coord_polar("y") +
            ggplot2::xlim(0.5, 2.5) +
            ggplot2::scale_fill_manual(values = c("Up" = "lightgreen", "Down" = "lightcoral", "Not sig" = "gray80")) +
            ggplot2::geom_text(ggplot2::aes(label = label), position = ggplot2::position_stack(vjust = 0.5), size = 3) +
            ggplot2::labs(title = paste0("RNA status of ", direction, "regulated proteins (n=", nrow(subset), ")")) +
            ggplot2::theme_void() +
            ggplot2::theme(legend.position = "none")

        fname <- file.path(out_dir, paste0("pie_", direction, ".png"))
        ggplot2::ggsave(fname, p, width = 6, height = 6, dpi = 200)
        saved <- c(saved, fname)
    }
    saved
}


.plot_stacked_bar <- function(conc_df, out_dir) {
    if (!"concordance" %in% colnames(conc_df)) return(character(0))

    counts <- as.data.frame(table(conc_df$concordance), stringsAsFactors = FALSE)
    colnames(counts) <- c("Category", "Count")
    counts$Category <- factor(counts$Category,
        levels = c("concordant", "discordant", "protein_only", "rna_only", "not_sig"))

    cols <- c("concordant" = "#2ca02c", "discordant" = "#d62728",
              "protein_only" = "#ff7f0e", "rna_only" = "#1f77b4", "not_sig" = "gray80")

    p <- ggplot2::ggplot(counts, ggplot2::aes(x = 1, y = Count, fill = Category)) +
        ggplot2::geom_col(width = 0.6) +
        ggplot2::scale_fill_manual(values = cols,
            labels = paste0(names(cols), " (", counts$Count[match(names(cols), counts$Category)], ")")) +
        ggplot2::labs(title = "RNA-Protein Concordance Distribution", y = "Number of gene-protein pairs") +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank())

    fname <- file.path(out_dir, "stacked_bar.png")
    ggplot2::ggsave(fname, p, width = 6, height = 8, dpi = 200)
    fname
}


.plot_volcano_overlay <- function(conc_df, out_dir) {
    saved <- character(0)
    conc_df$neg_log10_padj <- -log10(pmax(conc_df$padj_prot, 1e-300))

    # By protein significance
    p1 <- ggplot2::ggplot(conc_df, ggplot2::aes(x = logFC_prot, y = neg_log10_padj, color = prot_trend)) +
        ggplot2::geom_point(alpha = 0.5, size = 1) +
        ggplot2::scale_color_manual(values = c("Up" = "springgreen4", "Down" = "violetred", "Not sig" = "gray80")) +
        ggplot2::labs(title = "Protein Volcano — colored by protein significance",
                      x = "Protein log2FC", y = "-log10(padj)") +
        ggplot2::theme_minimal()
    fname1 <- file.path(out_dir, "volcano_by_prot.png")
    ggplot2::ggsave(fname1, p1, width = 8, height = 6, dpi = 200)
    saved <- c(saved, fname1)

    # By RNA change
    p2 <- ggplot2::ggplot(conc_df, ggplot2::aes(x = logFC_prot, y = neg_log10_padj, color = rna_change)) +
        ggplot2::geom_point(alpha = 0.5, size = 1) +
        ggplot2::scale_color_manual(values = c("Up" = "lightgreen", "Down" = "lightcoral", "Not sig" = "gray80")) +
        ggplot2::labs(title = "Protein Volcano — colored by RNA change status",
                      x = "Protein log2FC", y = "-log10(padj)") +
        ggplot2::theme_minimal()
    fname2 <- file.path(out_dir, "volcano_by_rna.png")
    ggplot2::ggsave(fname2, p2, width = 8, height = 6, dpi = 200)
    saved <- c(saved, fname2)

    saved
}


.plot_venn <- function(conc_df, out_dir) {
    sig_rna_n  <- sum(conc_df$sig_rna, na.rm = TRUE)
    sig_prot_n <- sum(conc_df$sig_prot, na.rm = TRUE)
    both_n     <- sum(conc_df$sig_rna & conc_df$sig_prot, na.rm = TRUE)
    rna_only_n <- sig_rna_n - both_n
    prot_only_n <- sig_prot_n - both_n

    df <- data.frame(
        Category = c("RNA only", "Both", "Protein only"),
        Count = c(rna_only_n, both_n, prot_only_n),
        stringsAsFactors = FALSE
    )
    df$Category <- factor(df$Category, levels = c("RNA only", "Both", "Protein only"))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Category, y = Count, fill = Category)) +
        ggplot2::geom_col() +
        ggplot2::geom_text(ggplot2::aes(label = Count), vjust = -0.5, size = 5) +
        ggplot2::scale_fill_manual(values = c("RNA only" = "#1f77b4", "Both" = "#2ca02c", "Protein only" = "#ff7f0e")) +
        ggplot2::labs(title = "Overlap: Significant RNA vs Protein Features",
                      y = "Number of features") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")

    fname <- file.path(out_dir, "venn.png")
    ggplot2::ggsave(fname, p, width = 6, height = 5, dpi = 200)
    fname
}


.plot_correlation_scatter <- function(conc_df, out_dir) {
    valid <- is.finite(conc_df$logFC_rna) & is.finite(conc_df$logFC_prot)
    df <- conc_df[valid, ]
    if (nrow(df) < 10) return(character(0))

    r <- cor(df$logFC_rna, df$logFC_prot, use = "pairwise.complete.obs")

    conc_col <- if ("concordance" %in% colnames(df)) "concordance" else NULL
    if (!is.null(conc_col)) {
        p <- ggplot2::ggplot(df, ggplot2::aes(x = logFC_rna, y = logFC_prot, color = .data[[conc_col]])) +
            ggplot2::geom_point(alpha = 0.4, size = 1) +
            ggplot2::scale_color_manual(values = c("concordant" = "#2ca02c", "discordant" = "#d62728",
                "rna_only" = "#1f77b4", "protein_only" = "#ff7f0e", "not_sig" = "gray80"))
    } else {
        p <- ggplot2::ggplot(df, ggplot2::aes(x = logFC_rna, y = logFC_prot)) +
            ggplot2::geom_point(alpha = 0.4, size = 1, color = "gray40")
    }

    p <- p +
        ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
        ggplot2::geom_hline(yintercept = 0, color = "gray80") +
        ggplot2::geom_vline(xintercept = 0, color = "gray80") +
        ggplot2::annotate("text", x = min(df$logFC_rna), y = max(df$logFC_prot),
                          label = sprintf("r = %.3f", r), hjust = 0, size = 5) +
        ggplot2::labs(title = "RNA vs Protein log2FC Correlation",
                      x = "RNA log2FC", y = "Protein log2FC") +
        ggplot2::theme_minimal()

    fname <- file.path(out_dir, "correlation_scatter.png")
    ggplot2::ggsave(fname, p, width = 8, height = 7, dpi = 200)
    fname
}


.plot_mhg_test <- function(conc_df, out_dir) {
    if (!requireNamespace("mHG", quietly = TRUE)) return(character(0))

    sig_prot <- conc_df[conc_df$sig_prot, ]
    if (nrow(sig_prot) < 10) return(character(0))

    results <- data.frame(
        Direction = character(0), mHG_stat = numeric(0),
        p_value = numeric(0), n_proteins = integer(0), n_rna_sig = integer(0),
        stringsAsFactors = FALSE
    )

    for (dir in c("Up", "Down")) {
        subset <- if (dir == "Up") {
            sig_prot[sig_prot$logFC_prot > 0, ]
        } else {
            sig_prot[sig_prot$logFC_prot < 0, ]
        }
        if (nrow(subset) < 5) next

        subset <- subset[order(abs(subset$logFC_prot), decreasing = TRUE), ]
        binary <- as.integer(subset$sig_rna)

        res <- tryCatch(mHG::mHG.test(binary), error = function(e) NULL)
        if (!is.null(res)) {
            results <- rbind(results, data.frame(
                Direction = paste0(dir, "regulated proteins"),
                mHG_stat = res$statistic,
                p_value = res$p.value,
                n_proteins = nrow(subset),
                n_rna_sig = sum(binary),
                stringsAsFactors = FALSE
            ))
        }
    }

    if (nrow(results) == 0) return(character(0))

    # Save as simple text plot
    fname <- file.path(out_dir, "mhg_test.png")
    grDevices::png(fname, width = 600, height = 300, res = 120)
    gridExtra::grid.table(results, rows = NULL)
    grDevices::dev.off()
    fname
}
