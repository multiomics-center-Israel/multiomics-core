#' GO hierarchy (DAG) helpers
#'
#' Build and draw the part of the Gene Ontology graph that sits above a set of
#' "seed" terms: the seeds, every ancestor reachable through is_a and part_of
#' links, and the links among them. The seeds come from an enrichment result
#' the caller already has; nothing here tests anything. Omics-agnostic: the
#' caller decides which terms are seeds and what their values mean.

# ==============================================================================
# GO.db ACCESS
# ==============================================================================

#' Parent links of one GO ontology, from GO.db
#'
#' Reads \code{GO<ONT>PARENTS} and keeps the is_a and part_of links, the two
#' relations along which GO annotations are inherited ("true path"). regulates
#' and its signed variants are left out: a term regulating a process is not a
#' kind or a part of that process, and drawing it as an ancestor would read as
#' one. The artificial root "all" is dropped.
#'
#' @param ontology One of "BP", "MF", "CC".
#' @return data.frame with character columns \code{child}, \code{parent} and
#'   \code{relation} ("is_a" or "part_of"), one row per link.
go_parent_table <- function(ontology = c("BP", "MF", "CC")) {
    ontology <- match.arg(ontology)
    if (!requireNamespace("GO.db", quietly = TRUE)) {
        stop("go_parent_table() needs the Bioconductor package GO.db, which is ",
             "not installed. Run renv::restore() to install the locked version.",
             call. = FALSE)
    }
    bimap <- switch(ontology,
                    BP = GO.db::GOBPPARENTS,
                    MF = GO.db::GOMFPARENTS,
                    CC = GO.db::GOCCPARENTS)
    lst <- AnnotationDbi::as.list(bimap)
    parents <- data.frame(
        child    = rep(names(lst), lengths(lst)),
        parent   = unname(unlist(lst, use.names = FALSE)),
        relation = unname(unlist(lapply(lst, names), use.names = FALSE)),
        stringsAsFactors = FALSE
    )
    normalize_go_relations(parents)
}

#' Normalise GO relation names and keep the inherited ones
#'
#' GO.db spells is_a as "isa". Rows whose relation is neither is_a nor part_of,
#' and rows whose parent is the artificial root "all", are dropped.
#'
#' @param parents data.frame with \code{child}, \code{parent}, \code{relation}.
#' @return The same columns, filtered, with relation in c("is_a", "part_of").
normalize_go_relations <- function(parents) {
    rel <- tolower(as.character(parents$relation))
    rel[rel %in% c("isa", "is-a", "is_a")] <- "is_a"
    rel[rel %in% c("part_of", "part-of", "partof")] <- "part_of"
    parents$relation <- rel
    keep <- parents$relation %in% c("is_a", "part_of") &
        !is.na(parents$parent) & parents$parent != "all"
    out <- unique(parents[keep, c("child", "parent", "relation"), drop = FALSE])
    rownames(out) <- NULL
    out
}

#' Names and ontologies of GO terms, from GO.db
#'
#' @param go_ids Character vector of GO identifiers.
#' @return data.frame with \code{go_id}, \code{term} and \code{ontology}
#'   (NA for an identifier GO.db does not know, e.g. an obsolete term), in the
#'   order of \code{unique(go_ids)}.
go_term_info <- function(go_ids) {
    go_ids <- unique(as.character(go_ids))
    out <- data.frame(go_id = go_ids, term = NA_character_,
                      ontology = NA_character_, stringsAsFactors = FALSE)
    if (length(go_ids) == 0 || !requireNamespace("GO.db", quietly = TRUE)) return(out)
    known <- intersect(go_ids, AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID"))
    if (length(known) == 0) return(out)
    info <- AnnotationDbi::select(GO.db::GO.db, keys = known,
                                  columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
    idx <- match(out$go_id, info$GOID)
    out$term <- info$TERM[idx]
    out$ontology <- info$ONTOLOGY[idx]
    out
}

# ==============================================================================
# GRAPH CONSTRUCTION
# ==============================================================================

#' All ancestors of a set of GO terms, the terms included
#'
#' Walks the parent table upwards from the seeds until no new term appears.
#'
#' @param seed_ids Character vector of GO identifiers.
#' @param parents Parent table as returned by go_parent_table().
#' @return Character vector: the seeds followed by their ancestors, each once.
go_dag_closure <- function(seed_ids, parents) {
    nodes <- unique(as.character(seed_ids))
    frontier <- nodes
    while (length(frontier) > 0) {
        up <- unique(parents$parent[parents$child %in% frontier])
        frontier <- setdiff(up, nodes)
        nodes <- c(nodes, frontier)
    }
    nodes
}

#' Build the GO sub-graph above a set of seed terms
#'
#' The node set is the seeds plus all their is_a/part_of ancestors; the edges
#' are the parent links with both ends in that set. When the node set would
#' exceed \code{max_nodes}, seeds are dropped from the end of \code{seeds}
#' (callers order them by adjusted p-value) until it fits, and the number
#' dropped is recorded, so a figure never silently shows fewer terms than were
#' significant.
#'
#' @param seeds data.frame with at least \code{go_id}; any further columns
#'   (e.g. \code{padj}, \code{NES}, \code{size}) are carried onto the seed
#'   nodes. Rows are in priority order.
#' @param parents Parent table as returned by go_parent_table().
#' @param term_info Optional data.frame with \code{go_id} and \code{term}; when
#'   NULL, names come from GO.db via go_term_info(). A seed GO.db does not
#'   name keeps the \code{term} column of \code{seeds}, if any.
#' @param max_nodes Largest node count to draw.
#' @return list(nodes, edges, n_seeds_requested, n_seeds_shown) where
#'   \code{nodes} has \code{go_id}, \code{term}, \code{is_seed} and the seed
#'   columns (NA on context nodes), and \code{edges} has \code{child},
#'   \code{parent}, \code{relation}.
build_go_dag <- function(seeds, parents, term_info = NULL, max_nodes = 120) {
    seeds <- seeds[!duplicated(seeds$go_id), , drop = FALSE]
    n_requested <- nrow(seeds)
    node_ids <- go_dag_closure(seeds$go_id, parents)
    while (length(node_ids) > max_nodes && nrow(seeds) > 1) {
        seeds <- seeds[-nrow(seeds), , drop = FALSE]
        node_ids <- go_dag_closure(seeds$go_id, parents)
    }

    edges <- parents[parents$child %in% node_ids & parents$parent %in% node_ids,
                     c("child", "parent", "relation"), drop = FALSE]
    edges <- unique(edges)
    rownames(edges) <- NULL

    if (is.null(term_info)) term_info <- go_term_info(node_ids)
    nodes <- data.frame(go_id = node_ids, stringsAsFactors = FALSE)
    nodes$term <- term_info$term[match(nodes$go_id, term_info$go_id)]
    if ("term" %in% colnames(seeds)) {
        miss <- is.na(nodes$term)
        nodes$term[miss] <- seeds$term[match(nodes$go_id[miss], seeds$go_id)]
    }
    nodes$term[is.na(nodes$term)] <- nodes$go_id[is.na(nodes$term)]
    nodes$is_seed <- nodes$go_id %in% seeds$go_id
    extra <- setdiff(colnames(seeds), c("go_id", "term"))
    for (col in extra) nodes[[col]] <- seeds[[col]][match(nodes$go_id, seeds$go_id)]

    list(nodes = nodes, edges = edges,
         n_seeds_requested = n_requested, n_seeds_shown = nrow(seeds))
}

# ==============================================================================
# PLOTTING
# ==============================================================================

#' Bin adjusted p-values for the node fill
#'
#' GOrilla-style bins on a fixed scale, so one colour means the same range in
#' every figure. Breaks at or above \code{fdr} are dropped and \code{fdr}
#' closes the last bin, so no bin lies outside what a seed can have.
#'
#' @param padj Numeric vector of adjusted p-values.
#' @param fdr The cutoff the seeds were selected with.
#' @return Ordered factor; levels run from the weakest bin to the strongest.
go_padj_bins <- function(padj, fdr = 0.05) {
    inner <- c(1e-9, 1e-7, 1e-5, 1e-3, 1e-2)
    inner <- inner[inner < fdr]
    breaks <- c(0, inner, fdr)
    fmt <- function(x) format(x, scientific = x < 1e-2, digits = 2)
    labels <- vapply(seq_len(length(breaks) - 1), function(i) {
        if (i == 1) paste0("< ", fmt(breaks[2]))
        else paste0(fmt(breaks[i]), " to ", fmt(breaks[i + 1]))
    }, character(1))
    bins <- cut(padj, breaks = breaks, labels = labels, right = FALSE,
                include.lowest = TRUE)
    factor(bins, levels = rev(labels), ordered = TRUE)
}

#' Fill colours of the adjusted p bins
#'
#' One single-hue scale shared by the static and the interactive figure, so a
#' colour means the same range in both.
#'
#' @param bin_levels Levels of go_padj_bins(), weakest bin first.
#' @return Named character vector of colours, one per level.
go_padj_palette <- function(bin_levels) {
    blues <- c("#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c")
    stats::setNames(utils::tail(blues, length(bin_levels)), bin_levels)
}

#' Shorten a label for display
#'
#' @param x Character vector.
#' @param width Maximum number of characters kept, ellipsis included.
#' @return Character vector.
truncate_go_label <- function(x, width = 40) {
    x <- as.character(x)
    long <- !is.na(x) & nchar(x) > width
    x[long] <- paste0(substr(x[long], 1, width - 3), "...")
    x
}

#' Draw a GO sub-graph
#'
#' Sugiyama (layered) layout with general terms at the top. Seed nodes are
#' filled by binned adjusted p-value on a single-hue scale and labelled;
#' context nodes (ancestors that are not seeds) are hollow grey and unlabelled.
#' is_a links are solid, part_of links dashed.
#'
#' @param dag Output of build_go_dag(); seed nodes need a \code{padj} column.
#' @param fdr The cutoff the seeds were selected with (sets the bins).
#' @param title,subtitle Plot title and subtitle.
#' @param label_width Maximum characters of a seed label.
#' @return A ggplot object (ggraph).
plot_go_dag <- function(dag, fdr = 0.05, title = NULL, subtitle = NULL,
                        label_width = 40) {
    nodes <- dag$nodes
    nodes$padj_bin <- go_padj_bins(ifelse(nodes$is_seed, nodes$padj, NA_real_), fdr)
    nodes$node_type <- factor(ifelse(nodes$is_seed, "seed", "context"),
                              levels = c("seed", "context"))
    nodes$label <- ifelse(nodes$is_seed, truncate_go_label(nodes$term, label_width), "")
    edges <- dag$edges

    # Edges run parent -> child so the layered layout puts general terms on top.
    graph <- igraph::graph_from_data_frame(
        d = data.frame(from = edges$parent, to = edges$child,
                       relation = edges$relation, stringsAsFactors = FALSE),
        vertices = nodes[, c("go_id", setdiff(colnames(nodes), "go_id"))],
        directed = TRUE)

    bin_levels <- levels(nodes$padj_bin)
    fill_values <- go_padj_palette(bin_levels)

    ggraph::ggraph(graph, layout = "sugiyama") +
        ggraph::geom_edge_link(ggplot2::aes(edge_linetype = relation),
                               edge_colour = "grey70", edge_width = 0.3) +
        ggraph::geom_node_point(ggplot2::aes(fill = padj_bin, colour = node_type,
                                             size = node_type),
                                shape = 21, stroke = 0.4) +
        ggraph::geom_node_text(ggplot2::aes(label = label), repel = TRUE,
                               size = 2.6, colour = "grey15",
                               max.overlaps = Inf, seed = 1) +
        ggraph::scale_edge_linetype_manual(
            values = c(is_a = "solid", part_of = "dashed"),
            name = "GO relation") +
        ggplot2::scale_fill_manual(values = fill_values, limits = bin_levels,
                                   drop = FALSE,
                                   na.value = "white", na.translate = FALSE,
                                   name = "Adjusted p (seed terms)") +
        ggplot2::scale_colour_manual(
            values = c(seed = "grey30", context = "grey70"),
            limits = c("seed", "context"),
            labels = c(seed = "Seed term",
                       context = "Context term (not significant or not tested)"),
            drop = FALSE, name = NULL) +
        ggplot2::scale_size_manual(values = c(seed = 3.2, context = 1.8),
                                   limits = c("seed", "context"),
                                   guide = "none") +
        ggplot2::guides(
            fill = ggplot2::guide_legend(order = 1,
                override.aes = list(size = 3.2, colour = "grey30")),
            colour = ggplot2::guide_legend(order = 2,
                override.aes = list(size = c(3.2, 1.8),
                                    fill = c(fill_values[[length(fill_values)]], "white")))) +
        ggplot2::labs(title = title, subtitle = subtitle) +
        ggplot2::theme_void(base_size = 10) +
        ggplot2::theme(legend.position = "right",
                       plot.title = ggplot2::element_text(size = 11),
                       plot.subtitle = ggplot2::element_text(size = 9, colour = "grey30"),
                       plot.margin = ggplot2::margin(8, 8, 8, 8))
}

#' Interactive GO sub-graph
#'
#' The same graph as plot_go_dag(), as a visNetwork widget that can be zoomed,
#' panned and searched, for graphs too dense to read at page width. Layered
#' top-down with general terms on top; every term is labelled in full (context
#' terms in grey), and hovering a node shows its GO id and, for a seed, NES,
#' adjusted p and set size. Fill, outline and line styles follow plot_go_dag().
#'
#' @param dag Output of build_go_dag(); seed nodes need a \code{padj} column.
#' @param fdr The cutoff the seeds were selected with (sets the bins).
#' @param title,subtitle Shown above the graph.
#' @param footer Optional text shown below the graph, e.g. the figure legend.
#' @param label_width Characters per line of a wrapped node label.
#' @return A visNetwork htmlwidget.
go_dag_widget <- function(dag, fdr = 0.05, title = NULL, subtitle = NULL,
                          footer = NULL, label_width = 28) {
    if (!requireNamespace("visNetwork", quietly = TRUE)) {
        stop("go_dag_widget() needs the package visNetwork, which is not ",
             "installed. Run renv::restore() to install the locked version.",
             call. = FALSE)
    }
    nodes <- dag$nodes
    is_seed <- nodes$is_seed
    bins <- go_padj_bins(ifelse(is_seed, nodes$padj, NA_real_), fdr)
    fill_values <- go_padj_palette(levels(bins))
    fill <- ifelse(is_seed, unname(fill_values[as.character(bins)]), "#ffffff")
    fill[is.na(fill)] <- "#ffffff"
    esc <- function(x) htmltools::htmlEscape(as.character(x))
    num <- function(col, digits) {
        v <- if (col %in% colnames(nodes)) nodes[[col]] else rep(NA_real_, nrow(nodes))
        ifelse(is.na(v), "NA", formatC(as.numeric(v), format = "g", digits = digits))
    }
    tip <- paste0("<b>", esc(nodes$term), "</b><br>", esc(nodes$go_id))
    tip <- ifelse(is_seed,
                  paste0(tip, "<br>NES ", num("NES", 3), ", adjusted p ",
                         num("padj", 2), ", size ", num("size", 4)),
                  paste0(tip, "<br>Context term (not significant or not tested)"))
    wrap <- vapply(as.character(nodes$term), function(x)
        paste(strwrap(x, width = label_width), collapse = "\n"), character(1))

    vn_nodes <- data.frame(
        id = nodes$go_id, label = unname(wrap), title = tip, shape = "dot",
        size = ifelse(is_seed, 14, 8),
        color.background = fill,
        color.border = ifelse(is_seed, "#4d4d4d", "#b3b3b3"),
        borderWidth = ifelse(is_seed, 1.5, 1),
        font.color = ifelse(is_seed, "#262626", "#8c8c8c"),
        font.size = ifelse(is_seed, 16, 13),
        stringsAsFactors = FALSE)
    edges <- dag$edges
    # parent -> child, so the hierarchical layout puts general terms on top.
    vn_edges <- data.frame(
        from = edges$parent, to = edges$child,
        dashes = edges$relation == "part_of",
        title = ifelse(edges$relation == "part_of", "part_of", "is_a"),
        color = "#b3b3b3", width = 1, stringsAsFactors = FALSE)

    legend_nodes <- data.frame(
        label = c(paste("Adjusted p", names(fill_values)), "Context term"),
        shape = "dot", size = 10,
        color.background = c(unname(fill_values), "#ffffff"),
        color.border = c(rep("#4d4d4d", length(fill_values)), "#b3b3b3"),
        stringsAsFactors = FALSE)
    legend_edges <- data.frame(label = c("is_a", "part_of"), color = "#8c8c8c",
                               dashes = c(FALSE, TRUE), arrows = "",
                               stringsAsFactors = FALSE)

    w <- visNetwork::visNetwork(vn_nodes, vn_edges, width = "100%", height = "85vh",
                                main = list(text = title %||% "",
                                            style = "font-family:sans-serif;font-size:16px;"),
                                submain = list(text = subtitle %||% "",
                                               style = "font-family:sans-serif;font-size:13px;color:#4d4d4d;"),
                                footer = list(text = footer %||% "",
                                              style = "font-family:sans-serif;font-size:12px;color:#4d4d4d;"))
    w <- visNetwork::visHierarchicalLayout(w, direction = "UD", sortMethod = "directed",
                                           levelSeparation = 150, nodeSpacing = 220,
                                           treeSpacing = 220, blockShifting = TRUE,
                                           edgeMinimization = TRUE, parentCentralization = TRUE)
    w <- visNetwork::visPhysics(w, enabled = FALSE)
    w <- visNetwork::visEdges(w, smooth = FALSE, arrows = list(to = list(enabled = FALSE)))
    w <- visNetwork::visInteraction(w, navigationButtons = TRUE, hover = TRUE,
                                    zoomView = TRUE, dragView = TRUE,
                                    tooltipDelay = 100)
    w <- visNetwork::visOptions(w, highlightNearest = list(enabled = TRUE, degree = 1,
                                                           hover = FALSE),
                                nodesIdSelection = list(enabled = TRUE,
                                                        main = "Find a GO term"))
    visNetwork::visLegend(w, useGroups = FALSE, addNodes = legend_nodes,
                          addEdges = legend_edges, position = "right", width = 0.18,
                          zoom = FALSE)
}

#' Save an interactive GO sub-graph as one HTML file
#'
#' Self-contained when pandoc is available, so the file opens on its own; else
#' the JavaScript goes to a shared \code{lib/} folder beside it.
#'
#' @param widget Output of go_dag_widget().
#' @param file Target .html path.
#' @param title Browser tab title.
#' @return \code{file}, invisibly.
save_go_dag_widget <- function(widget, file, title = "GO hierarchy") {
    file <- normalizePath(file, winslash = "/", mustWork = FALSE)
    if (rmarkdown::pandoc_available()) {
        htmlwidgets::saveWidget(widget, file, selfcontained = TRUE, title = title)
    } else {
        htmlwidgets::saveWidget(widget, file, selfcontained = FALSE, title = title,
                                libdir = "lib")
    }
    invisible(file)
}

#' Figure size for a GO sub-graph
#'
#' @param n_nodes Number of nodes drawn.
#' @return Named numeric vector c(width, height) in inches.
go_dag_figure_size <- function(n_nodes) {
    c(width  = min(20, max(8, 6 + 0.12 * n_nodes)),
      height = min(14, max(6, 4.5 + 0.06 * n_nodes)))
}
