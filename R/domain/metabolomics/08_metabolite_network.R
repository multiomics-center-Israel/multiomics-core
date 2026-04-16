# R/domain/metabolomics/08_metabolite_network.R
#
# Metabolite network construction using KEGG reaction pairs.
# Builds interactive HTML networks of DE metabolites connected by known
# biochemical reactions (substrate-product pairs from KEGG).
#
# Dependencies: KEGGREST, igraph, plotly, htmlwidgets


# ==== KEGG REACTION PAIR FETCHING ============================================

#' Fetch KEGG reaction pairs for a set of compound IDs
#'
#' Queries KEGG for all reactions involving the given compounds and extracts
#' substrate-product pairs. Returns only edges where BOTH compounds are in
#' the input set.
#'
#' @param kegg_ids Character vector of KEGG compound IDs (e.g., "C00022").
#' @return data.frame with columns: from, to (KEGG compound IDs).
fetch_kegg_reaction_pairs <- function(kegg_ids) {
    if (!requireNamespace("KEGGREST", quietly = TRUE))
        stop("KEGGREST package required for KEGG reaction pair fetching")

    kegg_ids <- unique(kegg_ids)
    message("Fetching KEGG reactions for ", length(kegg_ids), " compounds...")

    # Step 1: Find all reactions involving our compounds
    # Use KEGG LINK to get compound -> reaction mappings
    all_edges <- data.frame(from = character(0), to = character(0),
                            stringsAsFactors = FALSE)

    # Batch queries: KEGGREST::keggLink returns compound-reaction links
    # Process in batches of 10 to respect API limits
    batch_size <- 10
    reaction_ids <- character(0)

    for (i in seq(1, length(kegg_ids), by = batch_size)) {
        batch <- kegg_ids[i:min(i + batch_size - 1, length(kegg_ids))]
        tryCatch({
            # keggLink("reaction", "cpd:CXXXXX") returns reaction IDs for compound
            links <- KEGGREST::keggLink("reaction", paste0("cpd:", batch))
            if (length(links) > 0) {
                reaction_ids <- c(reaction_ids, unique(unname(links)))
            }
            Sys.sleep(0.3)  # rate limit
        }, error = function(e) {
            message("  Warning: KEGG query failed for batch starting at ", batch[1],
                    ": ", e$message)
        })
    }

    reaction_ids <- unique(reaction_ids)
    # Strip "rn:" prefix if present
    reaction_ids <- sub("^rn:", "", reaction_ids)
    message("  Found ", length(reaction_ids), " reactions involving input compounds")

    if (length(reaction_ids) == 0) return(all_edges)

    # Step 2: For each reaction, get substrate and product compounds
    # Process in batches of 10
    for (i in seq(1, length(reaction_ids), by = batch_size)) {
        batch <- reaction_ids[i:min(i + batch_size - 1, length(reaction_ids))]
        tryCatch({
            rxn_data <- KEGGREST::keggGet(batch)
            for (rxn in rxn_data) {
                # Extract compound IDs from EQUATION field
                substrates <- character(0)
                products   <- character(0)

                if (!is.null(rxn$EQUATION)) {
                    # EQUATION format: "C00001 + C00002 <=> C00003 + C00004"
                    parts <- strsplit(rxn$EQUATION, "<=>")[[1]]
                    if (length(parts) == 2) {
                        substrates <- regmatches(parts[1],
                                                 gregexpr("C[0-9]{5}", parts[1]))[[1]]
                        products   <- regmatches(parts[2],
                                                 gregexpr("C[0-9]{5}", parts[2]))[[1]]
                    }
                }

                # Create edges: every substrate-product pair where BOTH are in our set
                sub_in <- substrates[substrates %in% kegg_ids]
                prod_in <- products[products %in% kegg_ids]

                if (length(sub_in) > 0 && length(prod_in) > 0) {
                    pairs <- expand.grid(from = sub_in, to = prod_in,
                                         stringsAsFactors = FALSE)
                    # Remove self-loops
                    pairs <- pairs[pairs$from != pairs$to, ]
                    if (nrow(pairs) > 0) {
                        all_edges <- rbind(all_edges, pairs)
                    }
                }
            }
            Sys.sleep(0.3)  # rate limit
        }, error = function(e) {
            message("  Warning: KEGG reaction fetch failed for batch: ", e$message)
        })
    }

    # Deduplicate edges (treat as undirected)
    if (nrow(all_edges) > 0) {
        all_edges <- unique(all_edges)
        # Normalize: always put smaller ID first for undirected dedup
        normalized <- t(apply(all_edges, 1, function(row) sort(row)))
        all_edges <- unique(data.frame(from = normalized[, 1],
                                       to = normalized[, 2],
                                       stringsAsFactors = FALSE))
    }

    message("  Found ", nrow(all_edges), " unique reaction pair edges")
    all_edges
}


# ==== NETWORK CONSTRUCTION ====================================================

#' Build a DE metabolite network with KEGG reaction pair edges
#'
#' @param de_res           data.frame with columns: feature_id, logFC, P.Value.
#' @param feature_annotations data.frame with columns: feature_id, KEGG, plus others.
#' @param p_cutoff         P-value cutoff for DE significance (default 0.05).
#' @return list with:
#'   - graph: igraph object
#'   - nodes: data.frame (name, kegg_id, logFC, pvalue, neg_log10p, direction, color, size)
#'   - edges: data.frame (from, to)
#'   - n_total_de: total DE metabolites
#'   - n_with_kegg: DE metabolites with valid KEGG IDs
#'   - n_connected: DE metabolites with at least one edge
build_de_metabolite_network <- function(de_res, feature_annotations, p_cutoff = 0.05,
                                        remove_isolated = TRUE) {
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

    # Fetch reaction pair edges
    edges <- fetch_kegg_reaction_pairs(kegg_mets$KEGG)

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
        g <- igraph::graph_from_data_frame(
            edges[, c("from_name", "to_name")],
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
