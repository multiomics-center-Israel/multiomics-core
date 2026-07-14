# R/domain/metabolomics/08_metabolite_network.R
#
# Metabolite network construction from a pinned KEGG reaction-pair reference.
# Builds interactive HTML networks of DE metabolites connected by known
# biochemical reactions (substrate-product pairs from KEGG). The reaction pairs
# are read from a precomputed, checksum-pinned reference (see 08b) — this module
# never queries KEGG.
#
# Dependencies: igraph, plotly, htmlwidgets


# ==== REACTION-PAIR FILTERING ================================================

#' Filter a reaction-pair reference to undirected edges among a set of compounds
#'
#' Pure — no I/O and no KEGG access. Keeps reference rows where BOTH the
#' substrate and product are in `kegg_ids` (self-pairs excluded), collapses them
#' to undirected compound pairs, and aggregates the contributing reactions so
#' provenance survives the edge collapse.
#'
#' @param reaction_pairs data.frame with columns `reaction_id`, `substrate_id`,
#'   `product_id` (the validated reference); extra columns are ignored.
#' @param kegg_ids Character vector of KEGG compound IDs (C#####) to keep.
#' @return data.frame, one row per undirected pair: `from`, `to` (KEGG IDs,
#'   `from` <= `to`), `reaction_ids` (";"-joined, sorted unique reaction IDs) and
#'   `n_reactions` (count of unique reactions). Ordered by `from`, `to`.
filter_reaction_pairs_to_features <- function(reaction_pairs, kegg_ids) {
    empty <- data.frame(from = character(0), to = character(0),
                        reaction_ids = character(0), n_reactions = integer(0),
                        stringsAsFactors = FALSE)
    if (is.null(reaction_pairs) || nrow(reaction_pairs) == 0) return(empty)

    kegg_ids <- unique(kegg_ids)
    keep <- reaction_pairs$substrate_id %in% kegg_ids &
            reaction_pairs$product_id %in% kegg_ids &
            reaction_pairs$substrate_id != reaction_pairs$product_id
    rp <- reaction_pairs[keep, , drop = FALSE]
    if (nrow(rp) == 0) return(empty)

    # Undirected: order the two endpoints so A-B and B-A collapse to one edge.
    a <- pmin(rp$substrate_id, rp$product_id)
    b <- pmax(rp$substrate_id, rp$product_id)
    key <- paste(a, b, sep = "\r")

    rid_by_pair <- split(rp$reaction_id, key)
    keys <- sort(names(rid_by_pair))
    parts <- strsplit(keys, "\r", fixed = TRUE)
    out <- data.frame(
        from = vapply(parts, `[`, character(1), 1L),
        to   = vapply(parts, `[`, character(1), 2L),
        # reaction_ids: unique + deterministically sorted before joining;
        # n_reactions: count of unique reactions (edge provenance, corr. #5).
        reaction_ids = vapply(keys, function(k)
            paste(sort(unique(rid_by_pair[[k]])), collapse = ";"), character(1)),
        n_reactions = vapply(keys, function(k)
            length(unique(rid_by_pair[[k]])), integer(1)),
        stringsAsFactors = FALSE, row.names = NULL
    )
    out[order(out$from, out$to), , drop = FALSE]
}


# ==== NETWORK CONSTRUCTION ====================================================

#' Build a DE metabolite network with KEGG reaction pair edges
#'
#' @param de_res           data.frame with columns: feature_id, logFC, P.Value.
#' @param feature_annotations data.frame with columns: feature_id, KEGG, plus others.
#' @param reaction_pairs   Validated KEGG reaction-pair reference (from
#'   `read_kegg_reaction_pairs()`): columns reaction_id, substrate_id,
#'   product_id, ... Edges are derived from this — no KEGG queries are made.
#' @param p_cutoff         P-value cutoff for DE significance (default 0.05).
#' @param remove_isolated  Drop nodes with no edges (default TRUE).
#' @return list with:
#'   - graph: igraph object (edge attrs: reaction_ids, n_reactions)
#'   - nodes: data.frame (name, kegg_id, logFC, pvalue, neg_log10p, direction, color, size)
#'   - edges: data.frame (from, to, reaction_ids, n_reactions, from_name, to_name)
#'   - n_total_de: total DE metabolites
#'   - n_with_kegg: DE metabolites with valid KEGG IDs
#'   - n_connected: DE metabolites with at least one edge
build_de_metabolite_network <- function(de_res, feature_annotations, reaction_pairs,
                                        p_cutoff = 0.05, remove_isolated = TRUE) {
    if (!requireNamespace("igraph", quietly = TRUE))
        stop("igraph package required")

    # Filter to significant DE metabolites
    de_sig <- de_res[de_res$P.Value < p_cutoff, ]
    message("Building network for ", nrow(de_sig), " DE metabolites (p < ", p_cutoff, ")")

    # Merge with annotations to get KEGG IDs
    merged <- merge(de_sig, feature_annotations[, c("feature_id", "KEGG")],
                    by = "feature_id", all.x = FALSE)

    # Keep only valid KEGG compound IDs (C-numbers)
    merged$has_kegg <- grepl("^C[0-9]+$", merged$KEGG)
    kegg_mets <- merged[merged$has_kegg, ]
    message("  ", nrow(kegg_mets), " metabolites with valid KEGG IDs (of ",
            nrow(de_sig), " total DE)")

    if (nrow(kegg_mets) == 0) {
        warning("No DE metabolites with valid KEGG IDs — cannot build network")
        return(NULL)
    }

    # Reaction pair edges from the pinned reference (no KEGG queries).
    message("Building network from pinned reaction-pair reference (",
            nrow(reaction_pairs), " reference rows)")
    edges <- filter_reaction_pairs_to_features(reaction_pairs, kegg_mets$KEGG)
    message("  ", nrow(edges), " undirected reaction-pair edge(s) among DE metabolites")

    # Build node data
    nodes <- data.frame(
        name       = kegg_mets$feature_id,
        kegg_id    = kegg_mets$KEGG,
        logFC      = kegg_mets$logFC,
        pvalue     = kegg_mets$P.Value,
        neg_log10p = -log10(kegg_mets$P.Value),
        direction  = ifelse(kegg_mets$logFC > 0, "Up", "Down"),
        stringsAsFactors = FALSE
    )
    nodes$color <- ifelse(nodes$direction == "Up", "#E74C3C", "#3498DB")

    # Scale node size by -log10(p-value), with min/max bounds
    size_raw <- nodes$neg_log10p
    size_min <- 8
    size_max <- 30
    if (max(size_raw) > min(size_raw)) {
        nodes$size <- size_min + (size_max - size_min) *
            (size_raw - min(size_raw)) / (max(size_raw) - min(size_raw))
    } else {
        nodes$size <- (size_min + size_max) / 2
    }

    # Build igraph
    # Map edges from KEGG IDs to metabolite names
    kegg_to_name <- setNames(nodes$name, nodes$kegg_id)

    if (nrow(edges) > 0) {
        edges$from_name <- kegg_to_name[edges$from]
        edges$to_name   <- kegg_to_name[edges$to]
        # Drop edges where mapping failed (shouldn't happen but safety check)
        edges <- edges[!is.na(edges$from_name) & !is.na(edges$to_name), ]
    }

    if (nrow(edges) > 0) {
        # Columns 1-2 are the endpoints; reaction_ids / n_reactions ride along as
        # edge attributes so reaction provenance survives the edge collapse.
        g <- igraph::graph_from_data_frame(
            edges[, c("from_name", "to_name", "reaction_ids", "n_reactions")],
            directed = FALSE,
            vertices = nodes$name
        )
    } else {
        g <- igraph::make_empty_graph(n = nrow(nodes), directed = FALSE)
        igraph::V(g)$name <- nodes$name
    }

    # Set vertex attributes
    idx <- match(igraph::V(g)$name, nodes$name)
    igraph::V(g)$kegg_id    <- nodes$kegg_id[idx]
    igraph::V(g)$logFC      <- nodes$logFC[idx]
    igraph::V(g)$pvalue     <- nodes$pvalue[idx]
    igraph::V(g)$neg_log10p <- nodes$neg_log10p[idx]
    igraph::V(g)$direction  <- nodes$direction[idx]
    igraph::V(g)$color      <- nodes$color[idx]
    igraph::V(g)$size       <- nodes$size[idx]

    n_connected <- sum(igraph::degree(g) > 0)
    n_isolated  <- igraph::vcount(g) - n_connected

    # Remove isolated nodes (no edges) if requested
    if (remove_isolated && n_isolated > 0) {
        isolated_verts <- igraph::V(g)[igraph::degree(g) == 0]
        g <- igraph::delete_vertices(g, isolated_verts)
        nodes <- nodes[nodes$name %in% igraph::V(g)$name, ]
        message("  Removed ", n_isolated, " isolated nodes (no KEGG reaction pair edges)")
    }

    message("  Network: ", igraph::vcount(g), " nodes, ", igraph::ecount(g),
            " edges, ", n_connected, " connected, ", n_isolated, " isolated removed")

    list(
        graph            = g,
        nodes            = nodes,
        edges            = edges,
        n_total_de       = nrow(de_sig),
        n_with_kegg      = nrow(kegg_mets),
        n_connected      = n_connected,
        n_isolated       = n_isolated
    )
}


# ==== INTERACTIVE VISUALIZATION ===============================================

#' Create an interactive plotly network visualization
#'
#' @param network_result  Output from build_de_metabolite_network().
#' @param output_file     Path for output HTML file.
#' @param title           Plot title.
#' @return Invisible path to output file.
plot_metabolite_network_interactive <- function(network_result, output_file,
                                                title = "DE Metabolite Network") {
    if (!requireNamespace("plotly", quietly = TRUE))
        stop("plotly package required")
    if (!requireNamespace("htmlwidgets", quietly = TRUE))
        stop("htmlwidgets package required")

    g     <- network_result$graph
    nodes <- network_result$nodes

    # Layout
    set.seed(42)
    layout <- igraph::layout_with_fr(g)
    colnames(layout) <- c("x", "y")

    # Node coordinates
    node_x <- layout[, 1]
    node_y <- layout[, 2]

    # Edge traces
    edge_x <- c()
    edge_y <- c()
    el <- igraph::as_edgelist(g)
    if (nrow(el) > 0) {
        for (i in seq_len(nrow(el))) {
            src <- which(igraph::V(g)$name == el[i, 1])
            tgt <- which(igraph::V(g)$name == el[i, 2])
            edge_x <- c(edge_x, node_x[src], node_x[tgt], NA)
            edge_y <- c(edge_y, node_y[src], node_y[tgt], NA)
        }
    }

    # Build hover text
    idx <- match(igraph::V(g)$name, nodes$name)
    hover_text <- paste0(
        "<b>", nodes$name[idx], "</b><br>",
        "KEGG: ", nodes$kegg_id[idx], "<br>",
        "log2FC: ", round(nodes$logFC[idx], 2), "<br>",
        "p-value: ", signif(nodes$pvalue[idx], 3), "<br>",
        "Direction: ", nodes$direction[idx]
    )

    # Create plotly figure
    fig <- plotly::plot_ly()

    # Add edges
    if (length(edge_x) > 0) {
        fig <- plotly::add_trace(
            fig,
            x = edge_x, y = edge_y,
            type = "scatter", mode = "lines",
            line = list(color = "#CCCCCC", width = 1),
            hoverinfo = "none",
            showlegend = FALSE
        )
    }

    # Add nodes — separate traces for Up and Down for legend
    for (dir in c("Up", "Down")) {
        mask <- nodes$direction[idx] == dir
        if (!any(mask)) next
        fig <- plotly::add_trace(
            fig,
            x = node_x[mask], y = node_y[mask],
            type = "scatter", mode = "markers+text",
            marker = list(
                size = nodes$size[idx][mask],
                color = nodes$color[idx][mask],
                line = list(color = "#333333", width = 1)
            ),
            text = nodes$name[idx][mask],
            textposition = "top center",
            textfont = list(size = 8),
            hovertext = hover_text[mask],
            hoverinfo = "text",
            name = paste0(dir, "-regulated")
        )
    }

    # Layout
    n_info <- network_result
    subtitle <- sprintf(
        "%d DE metabolites | %d with KEGG IDs | %d connected | %d edges (KEGG reaction pairs)",
        n_info$n_total_de, n_info$n_with_kegg, n_info$n_connected,
        igraph::ecount(g)
    )

    fig <- plotly::layout(
        fig,
        title = list(
            text = paste0(title, "<br><sub>", subtitle, "</sub>"),
            font = list(size = 16)
        ),
        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE,
                     title = ""),
        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE,
                     title = ""),
        plot_bgcolor = "white",
        paper_bgcolor = "white",
        legend = list(x = 0.02, y = 0.98)
    )

    # Save
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    htmlwidgets::saveWidget(fig, output_file, selfcontained = TRUE)
    message("Network saved to: ", output_file)
    invisible(output_file)
}
