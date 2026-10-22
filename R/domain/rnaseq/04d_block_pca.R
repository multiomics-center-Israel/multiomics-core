#' PCA of RNA-seq expression before and after removing the blocking factor
#'
#' Shows a reader what pairing on \code{de$block_col} does. The "after" matrix
#' is [limma::removeBatchEffect()] with the block as batch and the contrast
#' factor protected: it removes each gene's estimated block effect from a
#' \code{~ factor + block} fit and keeps the within-block group differences
#' unchanged, the quantities a \code{~ block + factor} fit compares. Only for
#' balanced blocks is that the same as subtracting each block's mean. It is a picture only; no statistic is
#' computed on it. Both panels use the same genes, the \code{n_top} most
#' variable before correction, so the gene choice does not favour the factor.
#'
#' @param expr Numeric matrix (genes x samples) of log-scale expression,
#'   e.g. \code{pre$expr_work}.
#' @param meta Sample metadata with one row per column of \code{expr}.
#' @param sample_col Metadata column holding the sample IDs in \code{colnames(expr)}.
#' @param block_col Metadata column of the blocking factor (the one the DE fit used).
#' @param group_col Metadata column of the contrast factor.
#' @param n_top Number of most variable genes used in both PCAs.
#' @return A list with \code{coords} (data.frame: panel, sample, PC1, PC2,
#'   group, block), \code{var_explained} (data.frame: panel, PC1, PC2, in
#'   percent), \code{n_top} (genes used) and \code{block_sizes} (samples per
#'   block level).
build_rna_block_pca <- function(expr, meta, sample_col, block_col, group_col, n_top = 500) {
    for (col in c(sample_col, block_col, group_col)) {
        if (!col %in% colnames(meta)) {
            stop("build_rna_block_pca: column '", col, "' is not in the sample metadata.")
        }
    }
    idx <- match(colnames(expr), as.character(meta[[sample_col]]))
    if (anyNA(idx)) {
        stop("build_rna_block_pca: ", sum(is.na(idx)), " expression column(s) have no row in ",
             "the metadata under '", sample_col, "', e.g. ", colnames(expr)[is.na(idx)][1], ".")
    }
    meta <- meta[idx, , drop = FALSE]
    block <- factor(meta[[block_col]])
    group <- factor(meta[[group_col]])

    corrected <- limma::removeBatchEffect(expr, batch = block,
                                          design = stats::model.matrix(~ group))
    vars <- apply(expr, 1, stats::var)
    top <- order(vars, decreasing = TRUE)[seq_len(min(n_top, sum(is.finite(vars) & vars > 0)))]

    panels <- c(before = "Before", after = paste0("After removing ", block_col))
    run_pca <- function(m, key) {
        p <- stats::prcomp(t(m[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
        ve <- 100 * p$sdev^2 / sum(p$sdev^2)
        list(coords = data.frame(panel = panels[[key]], sample = colnames(m),
                                 PC1 = p$x[, 1], PC2 = p$x[, 2],
                                 group = group, block = block, row.names = NULL),
             ve = data.frame(panel = panels[[key]], PC1 = ve[1], PC2 = ve[2]))
    }
    b <- run_pca(expr, "before")
    a <- run_pca(corrected, "after")
    coords <- rbind(b$coords, a$coords)
    coords$panel <- factor(coords$panel, levels = unname(panels))

    list(coords = coords,
         var_explained = rbind(b$ve, a$ve),
         n_top = length(top),
         block_sizes = table(block))
}

#' Plot the before/after-block PCA built by build_rna_block_pca()
#'
#' Grey lines join the samples of each block level, so the reader can follow a
#' pair from one panel to the other.
#'
#' @param pca Result of [build_rna_block_pca()].
#' @param block_col Name of the blocking factor, for the legend.
#' @param group_col Name of the contrast factor, for the legend.
#' @return A ggplot object.
plot_rna_block_pca <- function(pca, block_col, group_col) {
    df <- pca$coords
    ve <- pca$var_explained
    labels <- stats::setNames(
        sprintf("%s\nPC1 %.0f%%, PC2 %.0f%% of variance", ve$panel, ve$PC1, ve$PC2), ve$panel)
    use_shape <- nlevels(df$block) <= 6

    p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2)) +
        ggplot2::geom_line(ggplot2::aes(group = block), colour = "grey60", linewidth = 0.5) +
        ggplot2::facet_wrap(~ panel, scales = "free", labeller = ggplot2::as_labeller(labels)) +
        ggplot2::labs(colour = group_col, shape = block_col,
                      subtitle = sprintf("Same %d most variable genes in both panels; grey lines join samples of the same %s.",
                                         pca$n_top, block_col)) +
        ggplot2::theme_bw(base_size = 14)
    if (use_shape) {
        p + ggplot2::geom_point(ggplot2::aes(colour = group, shape = block), size = 4)
    } else {
        p + ggplot2::geom_point(ggplot2::aes(colour = group), size = 4)
    }
}
