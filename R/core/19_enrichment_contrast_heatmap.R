# Enrichment across contrasts: one gene set per row, one contrast per column.
# Reads the per-contrast pathway_*.csv files the enrichment step writes; it does
# not rerun or re-threshold the enrichment itself.

#' Read every per-contrast enrichment result in a directory
#'
#' @param enrich_dir Directory holding \code{pathway_<contrast>_<db>_<method>.csv}.
#' @return Data frame with \code{pathway}, \code{pathway_name}, \code{contrast},
#'   \code{database}, \code{analysis} ("GSEA", "ORA up", "ORA down" or
#'   "cameraPR"), \code{NES} (NA except for GSEA), \code{padj} (cameraPR's FDR
#'   column) and \code{direction} ("up"/"down" for ORA and cameraPR, NA for
#'   GSEA); zero rows when nothing is found.
collect_enrichment_results <- function(enrich_dir) {
  empty <- data.frame(pathway = character(), pathway_name = character(),
                      contrast = character(), database = character(),
                      analysis = character(), NES = numeric(), padj = numeric(),
                      direction = character(), stringsAsFactors = FALSE)
  if (!dir.exists(enrich_dir)) return(empty)
  files <- list.files(enrich_dir, pattern = "^pathway_.+\\.csv$", full.names = TRUE)
  files <- files[!grepl("^pathway_volcano_data_", basename(files))]
  rows <- lapply(files, function(f) {
    d <- tryCatch(utils::read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    # cameraPR tables carry FDR and Direction (Up/Down) instead of padj/direction.
    if (!is.null(d) && nrow(d) > 0 && identical(d$method[1], "camerapr") &&
        all(c("FDR", "Direction") %in% names(d))) {
      d$padj <- d$FDR
      d$direction <- tolower(d$Direction)
    }
    need <- c("pathway", "padj", "contrast", "database", "method")
    if (is.null(d) || nrow(d) == 0 || !all(need %in% names(d))) return(NULL)
    analysis <- if (identical(d$method[1], "fgsea")) {
      "GSEA"
    } else if (identical(d$method[1], "camerapr")) {
      "cameraPR"
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
      direction = if ("direction" %in% names(d)) as.character(d$direction) else NA_character_,
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

#' Gene set x contrast matrix for one enrichment collection
#'
#' Rows are the gene sets with padj below \code{padj_cut} in at least one
#' contrast, ordered by their smallest padj and cut to \code{n_max}. Cells hold
#' NES for GSEA; for ORA and cameraPR, -log10(padj) signed by direction
#' (positive for the up list, negative for the down list) when the table
#' carries a direction, unsigned otherwise. Where a set was tested in both
#' directions in one contrast, the direction with the smaller padj is shown.
#' Cells of a contrast in which the set was not tested stay NA.
#'
#' @param res Rows of \code{\link{collect_enrichment_results}} for one
#'   database and one analysis family: GSEA rows, ORA up and down rows
#'   together, or cameraPR rows.
#' @param contrasts Contrast names, in column order.
#' @param padj_cut Significance cutoff (strict: padj < padj_cut).
#' @param n_max Most rows to keep.
#' @return NULL when no set passes; otherwise a list with matrices
#'   \code{value} and \code{padj} (sets x contrasts), \code{labels} (display
#'   names), \code{kind} ("NES", "signed_log10" or "log10"), \code{n_sig}
#'   (sets passing in any contrast), \code{n_shown}, \code{n_left_out} and
#'   \code{padj_cut}.
#' @examples
#' res <- data.frame(pathway = c("p1", "p1"), pathway_name = "P one",
#'                   contrast = c("a", "b"), analysis = "GSEA", NES = c(2, -1),
#'                   padj = c(0.01, 0.3), direction = NA)
#' build_enrichment_contrast_matrix(res, c("a", "b", "c"))$value
build_enrichment_contrast_matrix <- function(res, contrasts, padj_cut = 0.05, n_max = 40) {
  res <- res[res$contrast %in% contrasts & !is.na(res$padj), , drop = FALSE]
  if (nrow(res) == 0) return(NULL)
  is_gsea <- all(res$analysis == "GSEA")
  has_dir <- !is_gsea && "direction" %in% names(res) && any(!is.na(res$direction))
  kind <- if (is_gsea) "NES" else if (has_dir) "signed_log10" else "log10"

  # One row per set and contrast: the smaller padj when both directions exist.
  res <- res[order(res$pathway, res$contrast, res$padj), , drop = FALSE]
  res <- res[!duplicated(paste(res$pathway, res$contrast, sep = "\r")), , drop = FALSE]

  sig <- res[res$padj < padj_cut, , drop = FALSE]
  if (nrow(sig) == 0) return(NULL)
  best <- tapply(sig$padj, sig$pathway, min)
  ordered <- names(best)[order(best, names(best))]
  shown <- utils::head(ordered, n_max)

  key <- match(paste(rep(shown, times = length(contrasts)),
                     rep(contrasts, each = length(shown)), sep = "\r"),
               paste(res$pathway, res$contrast, sep = "\r"))
  padj <- matrix(res$padj[key], nrow = length(shown),
                 dimnames = list(shown, contrasts))
  value <- if (kind == "NES") {
    matrix(res$NES[key], nrow = length(shown), dimnames = list(shown, contrasts))
  } else {
    sgn <- if (kind == "signed_log10") ifelse(res$direction[key] %in% "down", -1, 1) else 1
    matrix(sgn * -log10(pmax(res$padj[key], 1e-300)), nrow = length(shown),
           dimnames = list(shown, contrasts))
  }
  names_by_id <- res$pathway_name[match(shown, res$pathway)]
  labels <- ifelse(is.na(names_by_id) | !nzchar(names_by_id), shown, names_by_id)

  list(value = value, padj = padj, labels = labels, kind = kind,
       n_sig = length(ordered), n_shown = length(shown),
       n_left_out = length(ordered) - length(shown), padj_cut = padj_cut)
}

#' Draw a matrix from \code{\link{build_enrichment_contrast_matrix}}
#'
#' Diverging scale centred at 0 for NES and signed -log10(padj), sequential for
#' unsigned -log10(padj). A dot marks padj below the cutoff; untested cells are
#' left blank. All text is \code{base_size} pt at the size the figure is saved.
#'
#' @param built Output of \code{\link{build_enrichment_contrast_matrix}}.
#' @param title Plot title.
#' @param caption Optional caption under the plot.
#' @param base_size Font size in pt.
#' @return ggplot object.
plot_enrichment_matrix_heatmap <- function(built, title = NULL, caption = NULL,
                                           base_size = 16) {
  sets <- rownames(built$value)
  contrasts <- colnames(built$value)
  labels <- make.unique(substr(built$labels, 1, 55))
  grid <- data.frame(
    set      = rep(sets, times = length(contrasts)),
    contrast = rep(contrasts, each = length(sets)),
    value    = as.vector(built$value),
    padj     = as.vector(built$padj),
    stringsAsFactors = FALSE
  )
  grid$label <- factor(labels[match(grid$set, sets)], levels = rev(labels))
  grid$contrast <- factor(gsub("_", " ", grid$contrast), levels = gsub("_", " ", contrasts))
  grid$sig <- !is.na(grid$padj) & grid$padj < built$padj_cut

  lim <- max(abs(grid$value), na.rm = TRUE)
  fill_scale <- switch(built$kind,
    NES = ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                        midpoint = 0, limits = c(-lim, lim),
                                        na.value = "white", name = "NES"),
    signed_log10 = ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                                 midpoint = 0, limits = c(-lim, lim),
                                                 na.value = "white",
                                                 name = "signed\n-log10(FDR)"),
    ggplot2::scale_fill_gradient(low = "white", high = "#54278F", limits = c(0, lim),
                                 na.value = "white", name = "-log10(FDR)")
  )
  ggplot2::ggplot(grid, ggplot2::aes(x = contrast, y = label, fill = value)) +
    ggplot2::geom_tile(colour = "grey85") +
    ggplot2::geom_point(data = grid[grid$sig, , drop = FALSE], size = 1.6,
                        colour = "grey20", show.legend = FALSE) +
    fill_scale +
    ggplot2::labs(title = title, caption = caption, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    theme_print_text(base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = base_size),
                   panel.grid = ggplot2::element_blank())
}
