#' Whole-matrix sample figures on the DE matrix (proteomics)
#'
#' For each protein set -- the top N by variance among fully measured proteins
#' (\code{report$sample_heatmap_top_n}, default 500, 1000, 2000) and all fully
#' measured proteins -- writes:
#' \itemize{
#'   \item a sample x sample Pearson correlation heatmap on the protein-centred
#'     matrix (pheatmap clustering of the correlation matrix, Euclidean /
#'     complete; colour scale fixed at -1..1);
#'   \item a samples x proteins heatmap of row z-scores, Euclidean / complete on
#'     both axes, colours clipped at +/-3.
#' }
#' When the DE model fitted \code{de$covariates}, also a PCA of the DE matrix
#' before and after removing them, on the fully measured proteins.
#'
#' The matrix is [resolve_proteomics_display_matrix()]: the DE matrix, with
#' the covariates removed when they were fitted. Figures go to
#' \code{Diagnostic_plots/DE_matrix/} as PNG and PDF, with
#' \code{figures_index.tsv} listing each figure and the facts its caption needs.
#'
#' @param display Output of [resolve_proteomics_display_matrix()] (the
#'   prot_display_matrix target, shared with the clustering heatmap).
#' @param config Full pipeline config.
#' @param out_dir Proteomics output directory.
#' @return list(files, index) -- written paths and the index data.frame (NULL
#'   when the figures are switched off).
mod_proteomics_sample_figures <- function(display, config, out_dir) {
    cfg <- config$modes$proteomics
    rep_cfg <- cfg$report %||% list()
    fig_dir <- file.path(create_legacy_output_dirs(out_dir)$diagnostic_plots, "DE_matrix")
    ensure_dir(fig_dir)
    # Everything here is rebuilt each run; an earlier run's set or PCA must not
    # be left behind for the report to pick up.
    old <- list.files(fig_dir, full.names = TRUE)
    if (length(old) > 0) file.remove(old)

    if (isFALSE(rep_cfg$show_sample_heatmaps)) {
        message("Sample heatmaps on the DE matrix: off (report$show_sample_heatmaps: false).")
        return(list(files = character(0), index = NULL))
    }

    annot <- build_proteomics_sample_annotation(display$meta, cfg)
    # pheatmap gets the columns in its drawing order; the PCA keeps `annot`.
    annot_hm <- annotation_for_pheatmap(annot, cfg)
    annot_colors <- build_sample_annotation_colors(annot_hm)
    n_samples <- ncol(display$mat)
    rdbu <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100)

    fully <- select_top_variable_measured(display$mat, display$observed, n = NULL)
    sets <- list()
    for (n in sort(as.integer(unlist(rep_cfg$sample_heatmap_top_n %||% c(500, 1000, 2000))))) {
        if (n < length(fully)) {
            sets[[paste0("top", n)]] <- select_top_variable_measured(display$mat, display$observed, n)
        } else {
            message(sprintf("  Sample heatmaps: top %d skipped, only %d proteins are fully measured.",
                            n, length(fully)))
        }
    }
    sets$fully_measured <- fully

    files <- character(0)
    index <- list()
    for (nm in names(sets)) {
        ids <- sets[[nm]]
        if (length(ids) < 3) {
            message(sprintf("  Sample heatmaps: '%s' skipped, %d proteins.", nm, length(ids)))
            next
        }
        x <- display$mat[ids, , drop = FALSE]
        set_label <- if (nm == "fully_measured") {
            sprintf("all %d fully measured proteins", length(ids))
        } else {
            sprintf("top %d fully measured proteins by variance", length(ids))
        }

        ph_cor <- plot_sample_correlation_heatmap(
            x, method = "pearson", annotation_col = annot_hm,
            main = sprintf("Sample correlation, %s\n(protein-centred, Pearson)", set_label),
            colors = rdbu, fontsize = 16, base_fontsize = 16,
            center_features = TRUE, limits = c(-1, 1), annotation_colors = annot_colors
        )
        f_cor <- file.path(fig_dir, sprintf("sample_correlation_heatmap_%s.png", nm))
        files <- c(files, save_heatmap_png_pdf(ph_cor, f_cor,
                                               width_in = 0.3 * n_samples + 10,
                                               height_in = 0.3 * n_samples + 7))

        ph_z <- plot_zscore_heatmap(
            x, annotation_col = annot_hm, annotation_colors = annot_colors,
            cluster_rows = TRUE, cluster_cols = TRUE,
            distance = "euclidean", linkage = "complete", zlim = 3, fontsize = 16,
            title = sprintf("%s, row z-score, Euclidean / complete", set_label)
        )
        f_z <- file.path(fig_dir, sprintf("samples_protein_heatmap_%s.png", nm))
        files <- c(files, save_heatmap_png_pdf(ph_z, f_z,
                                               width_in = 0.35 * n_samples + 8,
                                               height_in = 10))

        index[[length(index) + 1]] <- data.frame(
            figure = c("correlation", "zscore"), set = nm, n_proteins = length(ids),
            png = basename(c(f_cor, f_z)), matrix = display$label,
            stringsAsFactors = FALSE)
    }

    if (!is.null(display$covariates)) {
        pca <- write_covariate_pca_panels(display, annot, cfg, fully, fig_dir)
        files <- c(files, pca$files)
        if (!is.null(pca$index)) index[[length(index) + 1]] <- pca$index
    }

    index <- if (length(index) > 0) do.call(rbind, index) else NULL
    if (!is.null(index)) {
        f_idx <- file.path(fig_dir, "figures_index.tsv")
        utils::write.table(index, f_idx, sep = "\t", quote = FALSE, row.names = FALSE)
        files <- c(files, f_idx)
    }
    list(files = files, index = index)
}

#' PCA of the DE matrix before and after covariate removal
#'
#' Reuses \code{qc_pca_scatter()} (the pipeline's PCA) on the same proteins
#' for both panels, coloured by the first sample annotation column and shaped
#' by the second. A third panel colours the "before" scores by each covariate.
#'
#' @param display Output of [resolve_proteomics_display_matrix()] with
#'   covariates.
#' @param annot Output of [build_proteomics_sample_annotation()].
#' @param cfg The \code{modes$proteomics} config branch.
#' @param ids Proteins to use (the fully measured set).
#' @param fig_dir Output directory.
#' @return list(files, index).
write_covariate_pca_panels <- function(display, annot, cfg, ids, fig_dir) {
    if (length(ids) < 3) {
        message("  Covariate PCA skipped: fewer than 3 fully measured proteins.")
        return(list(files = character(0), index = NULL))
    }
    sample_col <- cfg$effects$samples %||% "SampleID"
    meta_pca <- display$meta
    cfg_pca <- cfg
    if (!is.null(annot) && ncol(annot) > 0) {
        for (col in names(annot)) {
            v <- annot[[col]][match(as.character(meta_pca[[sample_col]]), rownames(annot))]
            # Fill and shape are discrete scales here.
            meta_pca[[col]] <- if (is.numeric(v)) factor(v) else v
        }
        cfg_pca$effects$color <- names(annot)[1]
        cfg_pca$effects$shape <- if (ncol(annot) > 1) names(annot)[2] else NULL
    }

    label_layer <- function() {
        if (requireNamespace("ggrepel", quietly = TRUE)) {
            ggrepel::geom_text_repel(ggplot2::aes(label = .data[["sample"]]), size = 3.2,
                                     max.overlaps = Inf, show.legend = FALSE,
                                     colour = "grey25")
        } else {
            ggplot2::geom_text(ggplot2::aes(label = .data[["sample"]]), size = 3.2,
                               vjust = -0.8, check_overlap = TRUE, show.legend = FALSE,
                               colour = "grey25")
        }
    }
    # Scores and variance are read off qc_pca_scatter()'s plot before any layer
    # is added; adding to a ggplot builds a new object without those attributes.
    # Same colours as the heatmap annotation bars, so a group looks the same in
    # every figure.
    # Points are filled shapes (21, 24, ...) with a dark grey outline, so the
    # palest fill (0h) stays visible on white. The plot is drawn from
    # qc_pca_scatter()'s scores and variance; its own point layer is not used.
    ann_cols <- build_sample_annotation_colors(annot)
    fill_col <- cfg_pca$effects$color
    shape_col <- cfg_pca$effects$shape
    filled_shapes <- c(21, 24, 22, 23, 25)
    panel <- function(mat, title, subtitle = NULL) {
        p <- qc_pca_scatter(mat[ids, , drop = FALSE], meta_pca, cfg_pca, pcs = c(1, 2))
        sc <- attr(p, "scores")
        ve <- attr(p, "var_expl")
        if (is.null(shape_col)) {
            g <- ggplot2::ggplot(sc, ggplot2::aes(x = .data[["PC1"]], y = .data[["PC2"]],
                                                  fill = .data[[fill_col]])) +
                ggplot2::geom_point(size = 4.5, shape = 21, colour = "grey25", stroke = 0.7)
        } else {
            g <- ggplot2::ggplot(sc, ggplot2::aes(x = .data[["PC1"]], y = .data[["PC2"]],
                                                  fill = .data[[fill_col]],
                                                  shape = .data[[shape_col]])) +
                ggplot2::geom_point(size = 4.5, colour = "grey25", stroke = 0.7) +
                ggplot2::scale_shape_manual(values = filled_shapes, name = shape_col)
        }
        g <- g + label_layer()
        fc <- ann_cols[[fill_col]]
        g <- g + (if (is.null(fc)) ggplot2::scale_fill_discrete(name = fill_col)
                  else ggplot2::scale_fill_manual(values = fc, name = fill_col)) +
            # The fill legend needs a fillable key shape to show its colours.
            ggplot2::guides(fill = ggplot2::guide_legend(
                override.aes = list(shape = 21, colour = "grey25"), order = 1),
                shape = if (is.null(shape_col)) "none" else ggplot2::guide_legend(
                    override.aes = list(fill = "grey80", colour = "grey25"), order = 2)) +
            ggplot2::labs(title = title, subtitle = subtitle,
                          x = sprintf("PC1: %.2f%% of variance", 100 * ve[1]),
                          y = sprintf("PC2: %.2f%% of variance", 100 * ve[2])) +
            ggplot2::theme_minimal() +
            theme_print_text(16)
        list(plot = g, scores = sc, var_expl = ve)
    }
    before <- panel(display$before, "Before: DE matrix")
    # Percent variance is relative to what is left after removal, so PC1 can
    # explain a smaller share once a dominant axis is gone.
    after <- panel(display$mat, "After: covariate removed",
                   "% variance relative to the variation left after removal")

    scores <- before$scores
    axis_lab <- sprintf("PC%d: %.2f%% of variance", 1:2, 100 * before$var_expl[1:2])
    cov_panels <- lapply(colnames(display$covariates), function(cv) {
        scores$.cov <- meta_pca[[cv]][match(scores$sample, as.character(meta_pca[[sample_col]]))]
        ggplot2::ggplot(scores, ggplot2::aes(x = .data[["PC1"]], y = .data[["PC2"]],
                                             fill = .data[[".cov"]])) +
            ggplot2::geom_point(size = 4.5, shape = 21, colour = "grey25", stroke = 0.7) +
            label_layer() +
            ggplot2::scale_fill_gradient(low = "#DADAEB", high = "#3F007D", name = cv) +
            ggplot2::labs(title = paste0("Before, coloured by ", cv),
                          x = axis_lab[1], y = axis_lab[2]) +
            ggplot2::theme_minimal() +
            theme_print_text(16)
    })

    # Before and after side by side, each panel large enough to read; the
    # covariate-coloured view is its own figure.
    combined <- patchwork::wrap_plots(list(before$plot, after$plot), ncol = 2)
    f_png <- file.path(fig_dir, "PCA_covariate_before_after.png")
    f_pdf <- sub("\\.png$", ".pdf", f_png)
    ggplot2::ggsave(f_png, combined, width = 20, height = 9, dpi = 300)
    ggplot2::ggsave(f_pdf, combined, width = 20, height = 9)
    f_cov <- character(0)
    if (length(cov_panels) > 0) {
        f_cov <- file.path(fig_dir, "PCA_covariate_colouring.png")
        cov_fig <- patchwork::wrap_plots(cov_panels, ncol = length(cov_panels))
        ggplot2::ggsave(f_cov, cov_fig, width = 10 * length(cov_panels), height = 9, dpi = 300)
        ggplot2::ggsave(sub("\\.png$", ".pdf", f_cov), cov_fig,
                        width = 10 * length(cov_panels), height = 9)
        f_cov <- c(f_cov, sub("\\.png$", ".pdf", f_cov))
    }

    index <- data.frame(figure = "pca_before_after", set = "fully_measured",
                        n_proteins = length(ids), png = basename(f_png),
                        matrix = display$label, stringsAsFactors = FALSE)
    list(files = c(f_png, f_pdf, f_cov), index = index)
}
