#' Name the changed genes on a rendered KEGG map
#'
#' A KEGG box often stands for several genes -- one EC number, or a family --
#' and pathview prints only one name on it: the EC number, or the first gene of
#' the node. The box is coloured by whichever member changed, so a reader sees
#' "Plod1" coloured when the gene that moved was another member of the box.
#' These helpers write the changed genes' symbols just above each such box.


#' Which boxes carry a changed gene, and what to write above them
#'
#' Pure: no files, no drawing.
#'
#' @param node_data Gene-node table pathview returns as `plot.data.gene`, with
#'   `all.mapped` (comma-separated Entrez ids of the input genes on the node),
#'   `labels`, `x`, `y` (box centre, map pixels), `width` and `height`.
#' @param de_entrez Entrez ids of the genes that cleared the node rule.
#' @param symbols Named character vector, Entrez id -> gene symbol.
#' @param skip_if_shown Logical: leave a box alone when the name pathview
#'   already printed on it is exactly the changed gene. TRUE for gene-symbol
#'   maps; FALSE for EC-number maps, whose printed label is never a symbol.
#' @return Data frame with `x`, `y`, `width`, `height` and `label` (changed
#'   genes' symbols, sorted, joined by "/"), one row per box to annotate.
de_gene_node_labels <- function(node_data, de_entrez, symbols,
                                skip_if_shown = TRUE) {
    empty <- data.frame(x = numeric(0), y = numeric(0), width = numeric(0),
                        height = numeric(0), label = character(0),
                        stringsAsFactors = FALSE)
    need <- c("all.mapped", "x", "y", "width", "height")
    if (is.null(node_data) || !is.data.frame(node_data) || nrow(node_data) == 0 ||
        !all(need %in% names(node_data)) || length(de_entrez) == 0) {
        return(empty)
    }
    de_entrez <- as.character(de_entrez)
    mapped <- strsplit(as.character(node_data$all.mapped), ",", fixed = TRUE)
    label <- vapply(mapped, function(ids) {
        ids <- intersect(trimws(ids), de_entrez)
        if (length(ids) == 0) return(NA_character_)
        sym <- unname(symbols[ids])
        sym <- ifelse(is.na(sym) | !nzchar(sym), ids, sym)
        paste(sort(unique(sym)), collapse = "/")
    }, character(1))

    keep <- !is.na(label)
    if (skip_if_shown && "labels" %in% names(node_data)) {
        shown <- as.character(node_data$labels)
        keep <- keep & !(!is.na(shown) & !is.na(label) & label == shown)
    }
    if (!any(keep)) return(empty)
    out <- data.frame(x = as.numeric(node_data$x[keep]),
                      y = as.numeric(node_data$y[keep]),
                      width = as.numeric(node_data$width[keep]),
                      height = as.numeric(node_data$height[keep]),
                      label = label[keep], stringsAsFactors = FALSE)
    # One box can appear on several nodes of the same reaction; one label each.
    out[!duplicated(out[, c("x", "y")]), , drop = FALSE]
}


#' Write gene labels above boxes on a pathview PNG, in place
#'
#' pathview draws KEGG native maps at the KGML's own pixel size, so the node
#' coordinates it returns address the PNG directly.
#'
#' @param png_file Rendered pathview PNG.
#' @param labels Output of \code{de_gene_node_labels()}.
#' @return \code{png_file}, invisibly; the file is unchanged when there is
#'   nothing to write or the redraw fails.
annotate_pathview_png <- function(png_file, labels) {
    if (is.null(labels) || nrow(labels) == 0 || !file.exists(png_file)) {
        return(invisible(png_file))
    }
    img <- png::readPNG(png_file)
    h <- dim(img)[1]
    w <- dim(img)[2]
    tmp <- tempfile(fileext = ".png")
    before <- grDevices::dev.cur()
    ok <- tryCatch({
        # res = 72 makes one point one pixel, so the font size below is in
        # map pixels, like the KEGG boxes.
        grDevices::png(tmp, width = w, height = h, units = "px", res = 72)
        grid::grid.newpage()
        grid::grid.raster(img, interpolate = FALSE)
        grid::pushViewport(grid::viewport(xscale = c(0, w), yscale = c(h, 0)))
        fs <- 10
        for (i in seq_len(nrow(labels))) {
            lab <- labels$label[i]
            x <- labels$x[i]
            # Just above the box; below it when the box sits at the top edge.
            y <- labels$y[i] - labels$height[i] / 2 - fs / 2 - 2
            if (y < fs) y <- labels$y[i] + labels$height[i] / 2 + fs / 2 + 2
            tw <- grid::convertWidth(grid::stringWidth(lab), "points",
                                     valueOnly = TRUE) * fs / 12
            grid::grid.rect(x = x, y = y, width = tw + 4, height = fs + 2,
                            default.units = "native",
                            gp = grid::gpar(fill = "white", col = "grey30",
                                            lwd = 0.6))
            grid::grid.text(lab, x = x, y = y, default.units = "native",
                            gp = grid::gpar(fontsize = fs, col = "black"))
        }
        grid::popViewport()
        TRUE
    }, error = function(e) {
        message("    Pathview: gene labels not added to ", basename(png_file),
                " (", conditionMessage(e), ")")
        FALSE
    })
    if (!identical(grDevices::dev.cur(), before)) {
        tryCatch(grDevices::dev.off(), error = function(e) NULL)
    }
    if (isTRUE(ok) && file.exists(tmp) && file.size(tmp) > 0) {
        file.copy(tmp, png_file, overwrite = TRUE)
    }
    unlink(tmp)
    invisible(png_file)
}


#' Name the changed genes on one rendered map
#'
#' Glue for the renderers: never fails the render, only skips the labels.
#'
#' @param pv_out What \code{pathview::pathview()} returned for the map.
#' @param png_file The PNG it wrote.
#' @param de_entrez Entrez ids of the genes that cleared the node rule.
#' @param org_db OrgDb object for Entrez -> symbol, or NULL.
#' @param skip_if_shown Passed to \code{de_gene_node_labels()}.
#' @return \code{png_file}, invisibly.
label_de_genes_on_pathview <- function(pv_out, png_file, de_entrez, org_db,
                                       skip_if_shown = TRUE) {
    tryCatch({
        nodes <- pv_out$plot.data.gene
        if (is.null(nodes) || length(de_entrez) == 0) return(invisible(png_file))
        ids <- unique(as.character(de_entrez))
        symbols <- if (!is.null(org_db)) {
            suppressMessages(AnnotationDbi::mapIds(org_db, keys = ids,
                                                   column = "SYMBOL",
                                                   keytype = "ENTREZID",
                                                   multiVals = "first"))
        } else {
            stats::setNames(rep(NA_character_, length(ids)), ids)
        }
        labels <- de_gene_node_labels(nodes, ids, symbols, skip_if_shown)
        annotate_pathview_png(png_file, labels)
    }, error = function(e) {
        message("    Pathview: gene labels skipped for ", basename(png_file),
                " (", conditionMessage(e), ")")
        invisible(png_file)
    })
}


#' Entrez ids of the features that cleared the node rule
#'
#' @param de_table Standardized DE table for one contrast.
#' @param id_map Feature -> Entrez table from \code{map_feature_ids_to_entrez()}.
#' @return Character vector of Entrez ids; empty when nothing qualifies.
changed_feature_entrez <- function(de_table, id_map) {
    df <- filter_changed_features(de_table)
    if (is.null(df) || nrow(df) == 0 || is.null(id_map) || nrow(id_map) == 0) {
        return(character(0))
    }
    unique(as.character(merge(df, id_map, by = "feature_id")$ENTREZID))
}
