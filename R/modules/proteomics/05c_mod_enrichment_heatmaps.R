#' Enrichment summary heatmaps across contrasts (proteomics)
#'
#' One heatmap per enrichment collection (database) and analysis family --
#' fgsea (NES), ORA (up and down lists together, signed -log10 FDR), and on the
#' pooled ranking fgsea and cameraPR (signed -log10 FDR) -- with every contrast
#' as a column and the gene sets significant in at least one contrast as rows.
#' It reads the per-contrast tables the enrichment steps already wrote; nothing
#' is recomputed. STRING functional enrichment is not included: it is run once
#' on the DE proteins of all contrasts together, so it has no per-contrast
#' columns.
#'
#' Writes \code{Enrichment/contrast_heatmaps/<source>_<database>_<analysis>.png}
#' and \code{.pdf}, and \code{heatmaps_index.tsv} with one row per collection
#' (including collections with no significant set, which get no figure).
#'
#' @param pathway_res Result of \code{mod_proteomics_pathway()} (NULL when it
#'   did not run).
#' @param pooled_res Result of \code{mod_proteomics_pathway_pooled()} (NULL
#'   when it did not run).
#' @param config Full pipeline config.
#' @param out_dir Proteomics output directory.
#' @return list(files, index).
mod_proteomics_enrichment_heatmaps <- function(pathway_res, pooled_res, config, out_dir) {
    rep_cfg <- config$modes$proteomics$report %||% list()
    heat_dir <- file.path(out_dir, "Enrichment", "contrast_heatmaps")
    ensure_dir(heat_dir)
    old <- list.files(heat_dir, full.names = TRUE)
    if (length(old) > 0) file.remove(old)

    if (isFALSE(rep_cfg$show_enrichment_heatmaps)) {
        message("Enrichment summary heatmaps: off (report$show_enrichment_heatmaps: false).")
        return(list(files = character(0), index = NULL))
    }
    fdr <- as.numeric(rep_cfg$enrichment_heatmap_fdr %||% 0.05)
    n_max <- as.integer(rep_cfg$enrichment_heatmap_max_terms %||% 40)

    # The results objects say whether each step ran this time; the directories
    # alone could hold an earlier run's tables.
    sources <- list()
    if (length(pathway_res) > 0) {
        sources$pathway <- list(dir = file.path(out_dir, "Enrichment"),
                                contrasts = names(pathway_res),
                                label = "Enrichment/pathway_<contrast>_<database>_<method>.csv")
    }
    if (!is.null(pooled_res) && length(pooled_res$results) > 0) {
        sources$pooled <- list(dir = proteomics_pooled_enrichment_dir(out_dir),
                               contrasts = names(pooled_res$results),
                               label = "Enrichment/pooled_ranking/pathway_<contrast>_<database>_<method>.csv")
    }

    families <- list(GSEA = "GSEA", ORA = c("ORA up", "ORA down"), cameraPR = "cameraPR")
    files <- character(0)
    index <- list()
    for (src in names(sources)) {
        s <- sources[[src]]
        res <- collect_enrichment_results(s$dir)
        if (nrow(res) == 0) next
        for (db in sort(unique(res$database))) {
            for (fam in names(families)) {
                sub <- res[res$database == db & res$analysis %in% families[[fam]], , drop = FALSE]
                if (nrow(sub) == 0) next
                built <- build_enrichment_contrast_matrix(sub, s$contrasts, fdr, n_max)
                png_name <- NA_character_
                if (!is.null(built)) {
                    value_txt <- switch(built$kind,
                        NES = "NES",
                        signed_log10 = "-log10(FDR), signed by direction",
                        "-log10(FDR)")
                    caption <- if (built$n_left_out > 0) {
                        sprintf("Top %d of %d sets with FDR < %s; %d not shown.",
                                built$n_shown, built$n_sig, format(fdr), built$n_left_out)
                    } else {
                        sprintf("All %d sets with FDR < %s.", built$n_sig, format(fdr))
                    }
                    title <- sprintf("%s, %s%s: %s", db, fam,
                                     if (src == "pooled") " (pooled ranking)" else "", value_txt)
                    p <- plot_enrichment_matrix_heatmap(built, title = title, caption = caption)
                    safe <- gsub("[^A-Za-z0-9_-]", "_", paste(src, db, fam, sep = "_"))
                    png_name <- paste0(safe, ".png")
                    f_png <- file.path(heat_dir, png_name)
                    w <- 9 + 0.8 * ncol(built$value)
                    h <- 4.5 + 0.32 * nrow(built$value)
                    ggplot2::ggsave(f_png, p, width = w, height = h, dpi = 300, limitsize = FALSE)
                    ggplot2::ggsave(sub("\\.png$", ".pdf", f_png), p, width = w, height = h,
                                    limitsize = FALSE)
                    files <- c(files, f_png, sub("\\.png$", ".pdf", f_png))
                }
                index[[length(index) + 1]] <- data.frame(
                    source = src, database = db, analysis = fam,
                    value = if (is.null(built)) NA_character_ else built$kind,
                    n_significant = if (is.null(built)) 0L else built$n_sig,
                    n_shown = if (is.null(built)) 0L else built$n_shown,
                    n_left_out = if (is.null(built)) 0L else built$n_left_out,
                    fdr_cutoff = fdr, png = png_name, source_files = s$label,
                    stringsAsFactors = FALSE)
            }
        }
    }

    index <- if (length(index) > 0) do.call(rbind, index) else NULL
    if (!is.null(index)) {
        f_idx <- file.path(heat_dir, "heatmaps_index.tsv")
        utils::write.table(index, f_idx, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
        files <- c(files, f_idx)
    }
    list(files = files, index = index)
}
