# R/domain/lipidomics/20_links_enzyme_lipid_viz.R
#
# LipidONE / LipidLynxX-style enzyme-lipid pathway figure.
#
# Builds a directed metabolic reaction graph where:
#   * nodes  = lipid classes, coloured by the mean lipid logFC (per contrast)
#              of species in that class with sign retained;
#   * edges  = enzymatic reactions, coloured + thickness-weighted by the
#              gene/protein logFC of the catalysing enzyme family
#              (max-magnitude across worm orthologs, retrieved from the
#              external RNAseq / proteomics DE tables);
#   * labels = enzyme family name + representative worm gene symbol.
#
# The reaction backbone is the human reaction database from
# `lipid_reaction_database("human")` in 09_lipid_pathway.R, re-mapped to
# C. elegans via a curated human-symbol -> worm-symbol family table.
#
# One figure is rendered per contrast.


# ---------------------------------------------------------------------------
# Curated human-enzyme -> C. elegans gene-family map
# ---------------------------------------------------------------------------
.enzyme_human_to_cel <- function() {
    list(
        PLA2G2A = c("pla-1", "pla-2", "ipla-3", "ipla-7"),
        PLD1    = c("pld-1"),
        PLCG1   = c("plc-1", "plc-2", "plc-3", "plc-4", "egl-8", "itr-1"),
        LPIN1   = c("lpin-1"),
        DGAT1   = c("dgat-2", "mboa-2", "mboa-3"),
        DGKA    = c("dgk-1", "dgk-2", "dgk-3", "dgk-4", "dgk-5"),
        MOGAT1  = c("mboa-1", "mboa-2", "mboa-3"),
        LPCAT1  = c("mboa-7", "lpcat-1", "lcat-1", "lcat-3"),
        LPEAT1  = c("mboa-6", "lpcat-1"),
        PEMT    = c("pmt-1", "pmt-2"),
        CHPT1   = c("cept-1", "cept-2", "pcyt-1"),
        SELENOI = c("ept-1"),
        CPT1A   = c("cpt-1", "cpt-2", "cpt-3", "cpt-4", "cpt-5", "cpt-6"),
        SGMS1   = c("sms-1", "sms-2", "sms-3", "sms-4", "sms-5"),
        SMPD1   = c("asm-1", "asm-2", "asm-3"),
        UGCG    = c("cgt-1", "cgt-2", "cgt-3"),
        B4GALT5 = c("bre-1", "bre-3", "bre-5"),
        ST3GAL5 = c("st3gal-1"),
        CERS2   = c("hyl-1", "hyl-2", "lagr-1"),
        ASAH1   = c("asah-1", "asah-2"),
        SPHK1   = c("sphk-1", "sphk-2"),
        DEGS1   = c("ttm-5"),
        PI4KA   = c("pifk-1"),
        PIP5K1A = c("ppk-1", "ppk-3"),
        PLA2G7  = c("dpf-1", "ipla-3"),
        LIPE    = c("lid-1", "hsl-1", "lipl-1", "lipl-2", "lipl-3",
                    "lipl-4", "lipl-5", "lipl-7"),
        PNPLA2  = c("atgl-1"),
        MGLL    = c("faah-1", "mgl-1")
    )
}


# ---------------------------------------------------------------------------
# Build the reaction-table overlay with DE values per node + per edge.
# Returns a list:
#   $nodes : data.frame (lipid_class, mean_logFC, n_features, has_data)
#   $edges : data.frame (from, to, enzyme, gene_symbol_canonical,
#            cel_gene_symbols, enzyme_logFC, enzyme_layer,
#            enzyme_padj, reaction_type)
# ---------------------------------------------------------------------------
build_enzyme_lipid_overlay <- function(lipid_long, ext_long_list, contrast,
                                         organism = "cel") {
    # Reaction database — canonical (human-named) source
    rxn <- lipid_reaction_database("human")

    # ---- Nodes: aggregate lipid logFC per class for this contrast ----
    lp <- lipid_long[lipid_long$contrast == contrast, , drop = FALSE]
    if (nrow(lp) == 0) {
        return(NULL)
    }
    if (!"ClassKey" %in% colnames(lp)) lp$ClassKey <- NA_character_

    cls_agg <- stats::aggregate(
        cbind(logFC = lp$logFC) ~ ClassKey,
        data = lp,
        FUN = function(x) mean(x, na.rm = TRUE)
    )
    cls_n <- stats::aggregate(
        feature_id ~ ClassKey,
        data = lp,
        FUN = length
    )
    cls_n_sig <- stats::aggregate(
        is_sig ~ ClassKey,
        data = lp,
        FUN = function(x) sum(x, na.rm = TRUE)
    )
    node_tbl <- merge(cls_agg, cls_n,   by = "ClassKey", all = TRUE)
    node_tbl <- merge(node_tbl, cls_n_sig, by = "ClassKey", all = TRUE)
    names(node_tbl)[names(node_tbl) == "ClassKey"]  <- "lipid_class"
    names(node_tbl)[names(node_tbl) == "feature_id"] <- "n_features"
    names(node_tbl)[names(node_tbl) == "is_sig"]    <- "n_sig"
    node_tbl$has_data <- TRUE

    # Build superset of all classes referenced by the reactions
    rxn_classes <- unique(c(rxn$reactant_class, rxn$product_class))
    extra <- setdiff(rxn_classes, node_tbl$lipid_class)
    if (length(extra) > 0) {
        node_tbl <- rbind(node_tbl, data.frame(
            lipid_class = extra,
            logFC       = NA_real_,
            n_features  = 0L,
            n_sig       = 0L,
            has_data    = FALSE,
            stringsAsFactors = FALSE
        ))
    }
    names(node_tbl)[names(node_tbl) == "logFC"] <- "mean_logFC"

    # ---- Edges: enzyme DE overlay ----
    cel_map <- if (organism == "cel") .enzyme_human_to_cel() else list()

    # Combine rnaseq + proteomics into one feature_name-keyed table
    rna <- ext_long_list$rnaseq
    prt <- ext_long_list$proteomics

    de_pool <- list()
    if (!is.null(rna)) {
        rna_c <- rna[rna$contrast == contrast, , drop = FALSE]
        if (nrow(rna_c) > 0) {
            de_pool$rnaseq <- data.frame(
                gene_symbol = tolower(rna_c$feature_name),
                logFC       = rna_c$logFC,
                p_adj       = rna_c$p_adj,
                layer       = "rnaseq",
                stringsAsFactors = FALSE
            )
        }
    }
    if (!is.null(prt)) {
        prt_c <- prt[prt$contrast == contrast, , drop = FALSE]
        if (nrow(prt_c) > 0) {
            de_pool$proteomics <- data.frame(
                gene_symbol = tolower(prt_c$feature_name),
                logFC       = prt_c$logFC,
                p_adj       = prt_c$p_adj,
                layer       = "proteomics",
                stringsAsFactors = FALSE
            )
        }
    }
    de_all <- if (length(de_pool) > 0) do.call(rbind, de_pool) else NULL

    pick_best <- function(syms) {
        if (is.null(de_all) || length(syms) == 0) {
            return(list(logFC = NA_real_, gene = NA_character_,
                         layer = NA_character_, p_adj = NA_real_))
        }
        keys <- tolower(syms)
        hits <- de_all[de_all$gene_symbol %in% keys, , drop = FALSE]
        if (nrow(hits) == 0) {
            return(list(logFC = NA_real_, gene = NA_character_,
                         layer = NA_character_, p_adj = NA_real_))
        }
        # Take max-magnitude logFC across the family
        best_idx <- which.max(abs(hits$logFC))
        list(
            logFC = hits$logFC[best_idx],
            gene  = hits$gene_symbol[best_idx],
            layer = hits$layer[best_idx],
            p_adj = hits$p_adj[best_idx]
        )
    }

    edges <- lapply(seq_len(nrow(rxn)), function(i) {
        hsym <- rxn$gene_symbol[i]
        syms <- cel_map[[hsym]] %||% character(0)
        best <- pick_best(syms)
        data.frame(
            from       = rxn$reactant_class[i],
            to         = rxn$product_class[i],
            enzyme     = rxn$enzyme[i],
            human_gene = hsym,
            cel_genes  = paste(syms, collapse = ","),
            best_gene  = best$gene,
            enzyme_logFC = best$logFC,
            enzyme_layer = best$layer,
            enzyme_padj  = best$p_adj,
            reaction_type = rxn$reaction_type[i],
            stringsAsFactors = FALSE
        )
    })
    edge_tbl <- do.call(rbind, edges)

    list(nodes = node_tbl, edges = edge_tbl, contrast = contrast)
}


# ---------------------------------------------------------------------------
# Render the enzyme-lipid network image for one contrast.
# Returns a ggplot object, or NULL if igraph unavailable or no usable rows.
# ---------------------------------------------------------------------------
plot_enzyme_lipid_network <- function(overlay,
                                        keep_edges = c("any", "with_data"),
                                        title = NULL) {
    if (is.null(overlay)) return(NULL)
    nodes <- overlay$nodes
    edges <- overlay$edges
    if (is.null(nodes) || nrow(nodes) == 0 ||
        is.null(edges) || nrow(edges) == 0) {
        return(NULL)
    }
    keep_edges <- match.arg(keep_edges)

    # Restrict node set to those referenced by any retained edge plus those
    # with lipidomics data, then trim isolated nodes that have neither data
    # nor an edge.
    if (keep_edges == "with_data") {
        edges <- edges[!is.na(edges$enzyme_logFC), , drop = FALSE]
    }
    used_nodes <- unique(c(edges$from, edges$to,
                            nodes$lipid_class[nodes$has_data]))
    nodes <- nodes[nodes$lipid_class %in% used_nodes, , drop = FALSE]
    edges <- edges[edges$from %in% nodes$lipid_class &
                    edges$to   %in% nodes$lipid_class, , drop = FALSE]
    if (nrow(edges) == 0 || nrow(nodes) == 0) return(NULL)

    if (!requireNamespace("igraph", quietly = TRUE)) {
        message("[links/viz] igraph not available — falling back to circle")
        n_nodes <- nrow(nodes)
        angles  <- seq(0, 2 * pi, length.out = n_nodes + 1)[seq_len(n_nodes)]
        nodes$x  <- cos(angles)
        nodes$y  <- sin(angles)
    } else {
        g <- igraph::graph_from_data_frame(
            edges[, c("from", "to"), drop = FALSE],
            directed = TRUE,
            vertices = data.frame(name = nodes$lipid_class)
        )
        lay <- igraph::layout_with_fr(g, niter = 500)
        nodes$x <- lay[, 1]
        nodes$y <- lay[, 2]
    }

    # Edge segment coordinates
    map_xy <- function(cls, key) {
        node_idx <- match(cls, nodes$lipid_class)
        nodes[[key]][node_idx]
    }
    edges$x    <- map_xy(edges$from, "x")
    edges$y    <- map_xy(edges$from, "y")
    edges$xend <- map_xy(edges$to,   "x")
    edges$yend <- map_xy(edges$to,   "y")

    # Edge midpoint for labels
    edges$mid_x <- (edges$x + edges$xend) / 2
    edges$mid_y <- (edges$y + edges$yend) / 2

    # Visual encodings
    edges$has_de <- !is.na(edges$enzyme_logFC)
    edges$thickness <- ifelse(edges$has_de,
                               pmin(2.5, abs(edges$enzyme_logFC) * 1.5 + 0.5),
                               0.4)
    edges$alpha <- ifelse(edges$has_de, 0.85, 0.25)
    edges$label <- ifelse(edges$has_de,
                           paste0(edges$enzyme, ": ", edges$best_gene),
                           edges$enzyme)
    nodes$size <- 5 + 2 * pmin(5, log1p(nodes$n_features))

    p <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = edges,
            ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                          color = enzyme_logFC, linewidth = thickness,
                          alpha = alpha),
            arrow = grid::arrow(length = grid::unit(0.12, "inches"),
                                 type = "closed", angle = 25)
        ) +
        ggplot2::geom_point(
            data = nodes,
            ggplot2::aes(x = x, y = y, fill = mean_logFC, size = size),
            shape = 21, color = "black", stroke = 0.7
        ) +
        ggrepel::geom_text_repel(
            data = nodes,
            ggplot2::aes(x = x, y = y, label = lipid_class),
            size = 3.6, fontface = "bold", max.overlaps = Inf,
            point.padding = 0.3, box.padding = 0.4, seed = 1
        ) +
        ggrepel::geom_text_repel(
            data = edges[edges$has_de, ],
            ggplot2::aes(x = mid_x, y = mid_y, label = label),
            size = 2.8, color = "grey20", max.overlaps = Inf,
            min.segment.length = 0, seed = 2,
            box.padding = 0.2, point.padding = 0.1
        ) +
        ggplot2::scale_fill_gradient2(
            low = "#2166AC", mid = "white", high = "#B2182B",
            midpoint = 0, na.value = "grey85",
            name = "Lipid class\nmean log2FC"
        ) +
        ggplot2::scale_color_gradient2(
            low = "#2166AC", mid = "grey60", high = "#B2182B",
            midpoint = 0, na.value = "grey85",
            name = "Enzyme\nlog2FC"
        ) +
        ggplot2::scale_size_identity() +
        ggplot2::scale_linewidth_identity() +
        ggplot2::scale_alpha_identity() +
        ggplot2::labs(
            title = title %||% paste0("Enzyme ↔ Lipid network: ",
                                       overlay$contrast),
            subtitle = "Nodes: lipid classes — colour = mean log2FC of class species  •  Edges: KEGG reactions — colour = best worm-ortholog enzyme log2FC"
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.title    = ggplot2::element_text(face = "bold", size = 13,
                                                    hjust = 0.5),
            plot.subtitle = ggplot2::element_text(size  = 9, hjust = 0.5),
            legend.position = "right"
        )
    p
}


# ---------------------------------------------------------------------------
# Top-level: produce one PNG per contrast under <out_dir>/figures/ and
# also write the underlying overlay tables.
# ---------------------------------------------------------------------------
write_enzyme_lipid_figures <- function(lipid_long, ext_long_list,
                                         out_dir, organism = "cel",
                                         width = 11, height = 8, dpi = 200) {
    if (is.null(lipid_long) || nrow(lipid_long) == 0) return(character(0))
    fig_dir   <- file.path(out_dir, "figures")
    table_dir <- file.path(out_dir, "tables")
    ensure_dir(fig_dir); ensure_dir(table_dir)

    contrasts <- unique(lipid_long$contrast)
    out_files <- character(0)
    for (ctr in contrasts) {
        overlay <- build_enzyme_lipid_overlay(lipid_long, ext_long_list,
                                                contrast = ctr,
                                                organism = organism)
        if (is.null(overlay)) next

        nodes_file <- file.path(table_dir,
                                 paste0("enzyme_lipid_nodes.", ctr, ".tsv"))
        edges_file <- file.path(table_dir,
                                 paste0("enzyme_lipid_edges.", ctr, ".tsv"))
        utils::write.table(overlay$nodes, nodes_file, sep = "\t",
                           quote = FALSE, row.names = FALSE)
        utils::write.table(overlay$edges, edges_file, sep = "\t",
                           quote = FALSE, row.names = FALSE)
        out_files <- c(out_files, nodes_file, edges_file)

        for (mode in c("any", "with_data")) {
            p <- plot_enzyme_lipid_network(overlay, keep_edges = mode)
            if (is.null(p)) next
            png_file <- file.path(
                fig_dir,
                paste0("enzyme_lipid_network.", ctr, ".", mode, ".png")
            )
            tryCatch({
                ggplot2::ggsave(png_file, p, width = width, height = height,
                                 dpi = dpi)
                out_files <- c(out_files, png_file)
            }, error = function(e) {
                message("[links/viz] ggsave failed (", mode, "/", ctr,
                        "): ", conditionMessage(e))
            })
        }
    }
    out_files
}
