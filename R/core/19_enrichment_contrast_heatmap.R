# Enrichment across contrasts: one gene set per row, one contrast per column.
# Reads the per-contrast pathway_*.csv files the enrichment step writes; it does
# not rerun or re-threshold the enrichment itself.

#' Read every per-contrast enrichment result in a directory
#'
#' @param enrich_dir Directory holding \code{pathway_<contrast>_<db>_<method>.csv}.
#' @return Data frame with \code{pathway}, \code{pathway_name}, \code{contrast},
#'   \code{database}, \code{analysis} ("GSEA", "ORA up" or "ORA down"),
#'   \code{NES} (NA for ORA) and \code{padj}; zero rows when nothing is found.
collect_enrichment_results <- function(enrich_dir) {
  empty <- data.frame(pathway = character(), pathway_name = character(),
                      contrast = character(), database = character(),
                      analysis = character(), NES = numeric(), padj = numeric(),
                      stringsAsFactors = FALSE)
  if (!dir.exists(enrich_dir)) return(empty)
  files <- list.files(enrich_dir, pattern = "^pathway_.+\\.csv$", full.names = TRUE)
  files <- files[!grepl("^pathway_volcano_data_", basename(files))]
  rows <- lapply(files, function(f) {
    d <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    need <- c("pathway", "padj", "contrast", "database", "method")
    if (is.null(d) || nrow(d) == 0 || !all(need %in% names(d))) return(NULL)
    analysis <- if (identical(d$method[1], "fgsea")) {
      "GSEA"
    } else if ("direction" %in% names(d)) {
      paste("ORA", d$direction)
    } else {
      return(NULL)
    }
    data.frame(
      pathway = as.character(d$pathway),
      pathway_name = as.character(if ("pathway_name" %in% names(d)) d$pathway_name else d$pathway),
      contrast = as.character(d$contrast),
      database = as.character(d$database),
      analysis = analysis,
      NES = if ("NES" %in% names(d)) as.numeric(d$NES) else NA_real_,
      padj = as.numeric(d$padj),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(empty)
  do.call(rbind, rows)
}

#' Pick the gene sets for the three cross-contrast heatmaps
#'
#' For one database and analysis: \code{shared} are gene sets significant in
#' more than half of the contrasts, \code{unique} those significant in only one
#' or two, and \code{top} the best overall. Each list is ordered by the
#' smallest padj across contrasts and cut to \code{n}.
#'
#' @param res Rows of \code{\link{collect_enrichment_results}} for one database
#'   and analysis.
#' @param contrasts Contrast names, in display order.
#' @param padj_cut Significance cutoff on padj.
#' @param n Most gene sets per list.
#' @return List with character vectors \code{shared}, \code{unique}, \code{top}.
#' @examples
#' res <- data.frame(pathway = c("p1", "p1", "p2"), contrast = c("a", "b", "a"),
#'                   padj = c(0.01, 0.02, 0.001))
#' select_enrichment_heatmap_sets(res, c("a", "b"), 0.05)
select_enrichment_heatmap_sets <- function(res, contrasts, padj_cut = 0.05, n = 30) {
  res <- res[res$contrast %in% contrasts & !is.na(res$padj), , drop = FALSE]
  sig <- res[res$padj <= padj_cut, , drop = FALSE]
  if (nrow(sig) == 0) return(list(shared = character(), unique = character(), top = character()))
  n_sig <- tapply(sig$contrast, sig$pathway, function(x) length(unique(x)))
  best <- tapply(sig$padj, sig$pathway, min)
  ordered <- names(best)[order(best, names(best))]
  pick <- function(keep) utils::head(ordered[keep[ordered]], n)
  list(
    shared = pick(n_sig > length(contrasts) / 2),
    unique = pick(n_sig <= 2),
    top = utils::head(ordered, n)
  )
}

#' Heatmap of one enrichment statistic across contrasts
#'
#' GSEA cells show NES (diverging scale), ORA cells -log10(padj). Cells
#' significant at \code{padj_cut} carry a dot; gene sets not tested in a
#' contrast are left grey.
#'
#' @param res Rows of \code{\link{collect_enrichment_results}} for one database
#'   and analysis.
#' @param pathways Gene set IDs to show, top row first.
#' @param contrasts Contrast names, in display order.
#' @param padj_cut Significance cutoff on padj.
#' @param title Plot title.
#' @return ggplot object, or NULL when \code{pathways} is empty.
plot_enrichment_contrast_heatmap <- function(res, pathways, contrasts, padj_cut = 0.05,
                                             title = NULL) {
  if (length(pathways) == 0) return(NULL)
  is_gsea <- identical(unique(res$analysis), "GSEA")
  grid <- expand.grid(pathway = pathways, contrast = contrasts, stringsAsFactors = FALSE)
  key <- paste(res$pathway, res$contrast)
  hit <- match(paste(grid$pathway, grid$contrast), key)
  grid$padj <- res$padj[hit]
  grid$value <- if (is_gsea) res$NES[hit] else -log10(pmax(res$padj[hit], 1e-300))
  names_by_id <- res$pathway_name[match(pathways, res$pathway)]
  labels <- ifelse(is.na(names_by_id) | !nzchar(names_by_id), pathways, names_by_id)
  labels <- make.unique(substr(labels, 1, 60))
  grid$label <- factor(labels[match(grid$pathway, pathways)], levels = rev(labels))
  grid$contrast <- factor(gsub("_", " ", grid$contrast), levels = gsub("_", " ", contrasts))
  grid$sig <- !is.na(grid$padj) & grid$padj <= padj_cut

  fill_scale <- if (is_gsea) {
    lim <- max(abs(grid$value), na.rm = TRUE)
    ggplot2::scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                                  midpoint = 0, limits = c(-lim, lim),
                                  na.value = "grey90", name = "NES")
  } else {
    ggplot2::scale_fill_gradient(low = "white", high = "#984EA3",
                                 na.value = "grey90", name = "-log10(padj)")
  }
  ggplot2::ggplot(grid, ggplot2::aes(x = contrast, y = label, fill = value)) +
    ggplot2::geom_tile(colour = "grey85") +
    ggplot2::geom_point(data = grid[grid$sig, , drop = FALSE], size = 1.2,
                        colour = "grey20", show.legend = FALSE) +
    fill_scale +
    ggplot2::labs(title = title, x = NULL, y = NULL,
                  caption = sprintf("Dot: padj <= %s. Grey: not tested in that contrast.", padj_cut)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0),
                   plot.caption = ggplot2::element_text(hjust = 0))
}

#' Print the cross-contrast enrichment heatmaps as nested report tabs
#'
#' One tab per database, one sub-tab per analysis (GSEA, ORA up, ORA down),
#' each holding the shared, unique and top heatmaps. Meant for a knitr chunk
#' with \code{results='asis'}.
#'
#' @param res Output of \code{\link{collect_enrichment_results}}.
#' @param contrasts Contrast names, in display order.
#' @param gsea_cut,ora_cut padj cutoffs for GSEA and ORA, as the enrichment
#'   step used them.
#' @param n Most gene sets per heatmap.
#' @return Invisibly, the number of heatmaps printed.
print_enrichment_contrast_heatmaps <- function(res, contrasts, gsea_cut = 0.05,
                                               ora_cut = 0.05, n = 30) {
  kinds <- c(shared = "Shared by most contrasts",
             unique = "In one or two contrasts",
             top = sprintf("Top %d", n))
  n_printed <- 0L
  for (db in sort(unique(res$database))) {
    cat(sprintf("\n### %s {.tabset}\n\n", db))
    for (an in intersect(c("GSEA", "ORA up", "ORA down"), res$analysis[res$database == db])) {
      cat(sprintf("\n#### %s\n\n", an))
      sub <- res[res$database == db & res$analysis == an, , drop = FALSE]
      cut <- if (an == "GSEA") gsea_cut else ora_cut
      sets <- select_enrichment_heatmap_sets(sub, contrasts, cut, n)
      for (k in names(kinds)) {
        p <- plot_enrichment_contrast_heatmap(sub, sets[[k]], contrasts, cut,
                                              title = paste0(db, ", ", an, ": ", kinds[[k]]))
        if (is.null(p)) {
          cat(sprintf("\n%s: no gene set qualifies.\n\n", kinds[[k]]))
          next
        }
        print(p)
        cat("\n\n")
        n_printed <- n_printed + 1L
      }
    }
  }
  invisible(n_printed)
}
