#' Tri-Omics Feature Network
#'
#' Builds an integrated network connecting top features from transcriptomics,
#' proteomics, and metabolomics based on MOFA/DIABLO weights, co-loading,
#' and expression correlation.
#'
#' @name feature_network
NULL

#' Build Tri-Omics Feature Network
#'
#' Selects top features per omics from MOFA2 or DIABLO results, computes
#' cross-omics edges from cosine similarity, co-loading, and expression
#' correlation, then builds and plots an igraph network.
#'
#' @param integration_res Integration results (MOFA2/DIABLO/SNF).
#' @param harmonization_res Harmonization results including MAE and gene-protein mapping.
#' @param config Pipeline configuration list.
#' @param out_dir Output directory.
#' @return List with graph, nodes, and edges data frames.
#' @export
build_triomics_feature_network <- function(integration_res,
                                           harmonization_res,
                                           config,
                                           out_dir = NULL) {
    message("=== Building Tri-Omics Feature Network ===")

    if (!requireNamespace("igraph", quietly = TRUE)) {
        message("Package 'igraph' not installed. Skipping feature network.")
        return(NULL)
    }

    net_config <- config$modes$multiomics$feature_network %||% list()
    top_n <- net_config$top_n %||% 30
    mofa_cos_thresh <- net_config$mofa_cosine_threshold %||% 0.5
    diablo_load_thresh <- net_config$diablo_loading_threshold %||% 0.3
    corr_thresh <- net_config$correlation_threshold %||% 0.6

    # Create output directories
    if (!is.null(out_dir)) {
        dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    }

    # =========================================================================
    # Step A: Select top features per omics
    # =========================================================================
    top_features <- .select_top_features(integration_res, top_n)

    if (is.null(top_features) || nrow(top_features) == 0) {
        message("No top features found. Skipping feature network.")
        return(NULL)
    }

    message("  Selected ", nrow(top_features), " features across ",
            length(unique(top_features$view)), " omics views")

    # =========================================================================
    # Step B: Compute cross-omics edges
    # =========================================================================
    edges <- list()

    # B1: MOFA cosine similarity edges
    mofa_edges <- .compute_mofa_cosine_edges(
        integration_res$mofa_results, top_features, mofa_cos_thresh
    )
    if (!is.null(mofa_edges) && nrow(mofa_edges) > 0) {
        mofa_edges$type <- "mofa"
        edges[["mofa"]] <- mofa_edges
        message("  MOFA cosine edges: ", nrow(mofa_edges))
    }

    # B2: DIABLO co-loading edges
    diablo_edges <- .compute_diablo_coloading_edges(
        integration_res$diablo_results, top_features, diablo_load_thresh
    )
    if (!is.null(diablo_edges) && nrow(diablo_edges) > 0) {
        diablo_edges$type <- "diablo"
        edges[["diablo"]] <- diablo_edges
        message("  DIABLO co-loading edges: ", nrow(diablo_edges))
    }

    # B3: Expression correlation edges
    corr_edges <- .compute_correlation_edges(
        harmonization_res, top_features, corr_thresh
    )
    if (!is.null(corr_edges) && nrow(corr_edges) > 0) {
        corr_edges$type <- "correlation"
        edges[["correlation"]] <- corr_edges
        message("  Correlation edges: ", nrow(corr_edges))
    }

    if (length(edges) == 0) {
        message("No edges found above thresholds. Skipping feature network.")
        return(NULL)
    }

    all_edges <- do.call(rbind, edges)
    # Deduplicate: keep edge with highest weight for each pair
    all_edges$pair_key <- apply(
        all_edges[, c("from", "to")], 1,
        function(x) paste(sort(x), collapse = "___")
    )
    all_edges <- all_edges[order(all_edges$weight, decreasing = TRUE), ]
    all_edges <- all_edges[!duplicated(all_edges$pair_key), ]
    all_edges$pair_key <- NULL

    message("  Total unique edges: ", nrow(all_edges))

    # =========================================================================
    # Step C: Build igraph and plot
    # =========================================================================
    g <- igraph::graph_from_data_frame(
        all_edges[, c("from", "to")],
        directed = FALSE
    )

    # Node attributes from top_features
    node_info <- top_features[match(igraph::V(g)$name, top_features$feature_id), ]

    # Color by omics view
    view_colors <- c(
        transcriptomics = "#d73027",
        proteomics = "#4575b4",
        metabolomics = "#1b7837"
    )
    node_views <- node_info$view
    node_views[is.na(node_views)] <- "unknown"
    igraph::V(g)$color <- vapply(
        node_views,
        function(v) if (v %in% names(view_colors)) view_colors[[v]] else "gray60",
        character(1)
    )
    igraph::V(g)$view <- node_views

    # Node size by importance (rescale abs_weight to 4-15 range)
    importance <- node_info$importance
    importance[is.na(importance)] <- 0
    mn <- min(importance, na.rm = TRUE)
    mx <- max(importance, na.rm = TRUE)
    if (mx > mn) {
        scaled <- (importance - mn) / (mx - mn)
    } else {
        scaled <- rep(0.5, length(importance))
    }
    igraph::V(g)$size <- 4 + scaled * 11

    # Edge color by type
    edge_type_colors <- c(
        mofa = "gray70",
        diablo = "#e66101",
        correlation = "#5e3c99"
    )
    edge_types <- all_edges$type[match(
        apply(igraph::as_edgelist(g), 1, function(x) paste(sort(x), collapse = "___")),
        apply(all_edges[, c("from", "to")], 1, function(x) paste(sort(x), collapse = "___"))
    )]
    edge_types[is.na(edge_types)] <- "mofa"
    igraph::E(g)$color <- vapply(
        edge_types,
        function(et) if (et %in% names(edge_type_colors)) edge_type_colors[[et]] else "gray70",
        character(1)
    )

    # Resolve GENE_N to original IDs for labels
    gpm <- harmonization_res$gene_protein_mapping
    resolved_names <- .resolve_feature_ids(igraph::V(g)$name, igraph::V(g)$view, gpm)
    igraph::V(g)$original_id <- resolved_names

    # Labels: show top nodes only (top 15 by importance)
    sorted_idx <- order(importance, decreasing = TRUE)
    show_label <- rep(FALSE, length(igraph::V(g)))
    show_label[utils::head(sorted_idx, 15)] <- TRUE
    igraph::V(g)$label <- ifelse(show_label, resolved_names, NA)

    # Plot
    if (!is.null(out_dir)) {
        grDevices::png(
            file.path(out_dir, "plots", "triomics_feature_network.png"),
            width = 14, height = 14, units = "in", res = 300
        )

        tryCatch({
            layout <- igraph::layout_with_fr(g)

            igraph::plot.igraph(
                g,
                vertex.size = igraph::V(g)$size,
                vertex.color = igraph::V(g)$color,
                vertex.label = igraph::V(g)$label,
                vertex.label.cex = 0.65,
                vertex.label.color = "black",
                vertex.frame.color = NA,
                edge.color = igraph::E(g)$color,
                edge.width = 1.2,
                layout = layout,
                main = "Tri-Omics Feature Network"
            )

            # Legend
            legend("bottomleft",
                   legend = c("Transcriptomics", "Proteomics", "Metabolomics"),
                   col = c("#d73027", "#4575b4", "#1b7837"),
                   pch = 19, pt.cex = 2, bty = "n", title = "Omics Layer")
            legend("bottomright",
                   legend = c("MOFA cosine", "DIABLO co-loading", "Correlation"),
                   col = c("gray70", "#e66101", "#5e3c99"),
                   lwd = 2, bty = "n", title = "Edge Source")

        }, error = function(e) {
            message("  Network plot failed: ", e$message)
        })

        grDevices::dev.off()
        message("  Saved triomics_feature_network.png")
    }

    # =========================================================================
    # Step D: Save tables
    # =========================================================================
    node_degrees <- igraph::degree(g)
    nodes_df <- data.frame(
        feature_id = igraph::V(g)$name,
        original_id = igraph::V(g)$original_id,
        view = igraph::V(g)$view,
        importance = importance,
        degree = node_degrees[igraph::V(g)$name],
        stringsAsFactors = FALSE
    )

    # Resolve edge IDs too
    edge_from_views <- igraph::V(g)$view[match(all_edges$from, igraph::V(g)$name)]
    edge_to_views <- igraph::V(g)$view[match(all_edges$to, igraph::V(g)$name)]
    edges_df <- data.frame(
        from = all_edges$from,
        from_original = .resolve_feature_ids(all_edges$from, edge_from_views, gpm),
        to = all_edges$to,
        to_original = .resolve_feature_ids(all_edges$to, edge_to_views, gpm),
        weight = all_edges$weight,
        type = all_edges$type,
        stringsAsFactors = FALSE
    )

    if (!is.null(out_dir)) {
        write.csv(nodes_df, file.path(out_dir, "tables", "triomics_network_nodes.csv"),
                  row.names = FALSE)
        write.csv(edges_df, file.path(out_dir, "tables", "triomics_network_edges.csv"),
                  row.names = FALSE)
        message("  Saved node/edge tables")
    }

    list(
        graph = g,
        nodes = nodes_df,
        edges = edges_df
    )
}


# =============================================================================
# Internal helpers
# =============================================================================

#' Resolve GENE_N identifiers back to original gene/protein IDs
#'
#' MOFA/harmonization creates GENE_1, GENE_2, ... from gene_protein_mapping rows.
#' Feature IDs may also carry view suffixes like _transcriptomics, _proteomics.
#' This function maps them back to original identifiers.
#'
#' @param feature_ids Character vector of feature IDs (e.g., "GENE_3978_transcriptomics")
#' @param views Character vector of omics view names (same length as feature_ids)
#' @param gene_protein_mapping Data frame with gene_id, protein_id columns (row N = GENE_N)
#' @return Character vector of resolved IDs (original IDs where possible)
#' @noRd
.resolve_feature_ids <- function(feature_ids, views = NULL, gene_protein_mapping = NULL) {
    if (is.null(gene_protein_mapping) || !is.data.frame(gene_protein_mapping) ||
        nrow(gene_protein_mapping) == 0) {
        return(feature_ids)
    }

    resolved <- feature_ids
    view_suffixes <- c("_transcriptomics", "_proteomics", "_metabolomics",
                       "_rna", "_protein", "_metab")

    for (i in seq_along(feature_ids)) {
        fid <- feature_ids[i]
        if (is.na(fid)) next

        # Strip view suffix to get base ID
        base_id <- fid
        detected_view <- views[i]
        for (sfx in view_suffixes) {
            if (endsWith(base_id, sfx)) {
                base_id <- sub(paste0(sfx, "$"), "", base_id)
                if (is.null(detected_view) || is.na(detected_view)) {
                    detected_view <- gsub("^_", "", sfx)
                }
                break
            }
        }

        # Check if it's a GENE_N pattern
        if (!grepl("^GENE_\\d+$", base_id)) next

        row_n <- as.integer(sub("^GENE_", "", base_id))
        if (is.na(row_n) || row_n < 1 || row_n > nrow(gene_protein_mapping)) next

        # Resolve based on view
        vw <- tolower(detected_view %||% "")
        if (grepl("prot", vw)) {
            resolved[i] <- gene_protein_mapping$protein_id[row_n]
        } else {
            # Default to gene_id for transcriptomics/unknown
            # Use gene_symbol if available and non-empty
            sym <- gene_protein_mapping$gene_symbol[row_n]
            gid <- gene_protein_mapping$gene_id[row_n]
            if (!is.null(sym) && !is.na(sym) && nzchar(sym) && sym != gid) {
                resolved[i] <- sym
            } else {
                resolved[i] <- gid
            }
        }
    }

    resolved
}

#' Select top features per omics from MOFA or DIABLO
#' @noRd
.select_top_features <- function(integration_res, top_n) {

    # Try MOFA top_features first
    if (!is.null(integration_res$mofa_results$top_features)) {
        tf <- integration_res$mofa_results$top_features
        if (is.data.frame(tf) && nrow(tf) > 0 &&
            all(c("feature_id", "abs_weight", "view") %in% colnames(tf))) {

            # Aggregate max abs_weight per feature
            agg <- stats::aggregate(abs_weight ~ feature_id + view, data = tf, FUN = max)

            # Take top N per view
            result <- do.call(rbind, lapply(split(agg, agg$view), function(vdf) {
                vdf <- vdf[order(vdf$abs_weight, decreasing = TRUE), ]
                utils::head(vdf, top_n)
            }))

            if (nrow(result) > 0) {
                result$importance <- result$abs_weight
                rownames(result) <- NULL
                return(result[, c("feature_id", "view", "importance")])
            }
        }
    }

    # Fallback: DIABLO loadings
    if (!is.null(integration_res$diablo_results$loadings)) {
        loadings <- integration_res$diablo_results$loadings
        loadings <- loadings[names(loadings) != "Y"]

        if (length(loadings) > 0) {
            rows <- list()
            for (view_name in names(loadings)) {
                mat <- loadings[[view_name]]
                if (!is.matrix(mat) && !is.data.frame(mat)) next
                mat <- as.matrix(mat)
                # Max abs loading across components per feature
                max_load <- apply(abs(mat), 1, max, na.rm = TRUE)
                df <- data.frame(
                    feature_id = rownames(mat),
                    view = view_name,
                    importance = max_load,
                    stringsAsFactors = FALSE
                )
                df <- df[order(df$importance, decreasing = TRUE), ]
                rows[[view_name]] <- utils::head(df, top_n)
            }
            result <- do.call(rbind, rows)
            if (!is.null(result) && nrow(result) > 0) {
                rownames(result) <- NULL
                return(result)
            }
        }
    }

    NULL
}


#' Compute MOFA cosine similarity edges between features across views
#' @noRd
.compute_mofa_cosine_edges <- function(mofa_results, top_features, threshold) {
    if (is.null(mofa_results$weights)) return(NULL)

    weights <- mofa_results$weights
    if (!is.list(weights) || length(weights) < 2) return(NULL)

    # Build weight vectors for selected features
    feature_vectors <- list()
    for (view_name in names(weights)) {
        mat <- weights[[view_name]]
        if (!is.matrix(mat) && !is.data.frame(mat)) next
        mat <- as.matrix(mat)
        view_features <- top_features$feature_id[top_features$view == view_name]
        common <- intersect(view_features, rownames(mat))
        if (length(common) > 0) {
            for (feat in common) {
                feature_vectors[[feat]] <- list(
                    vec = mat[feat, ],
                    view = view_name
                )
            }
        }
    }

    if (length(feature_vectors) < 2) return(NULL)

    feat_names <- names(feature_vectors)
    edges <- list()

    for (i in seq_len(length(feat_names) - 1)) {
        fi <- feat_names[i]
        vi <- feature_vectors[[fi]]
        for (j in (i + 1):length(feat_names)) {
            fj <- feat_names[j]
            vj <- feature_vectors[[fj]]

            # Only cross-omics edges
            if (vi$view == vj$view) next

            # Cosine similarity
            dot <- sum(vi$vec * vj$vec, na.rm = TRUE)
            norm_i <- sqrt(sum(vi$vec^2, na.rm = TRUE))
            norm_j <- sqrt(sum(vj$vec^2, na.rm = TRUE))

            if (norm_i == 0 || norm_j == 0) next

            cos_sim <- dot / (norm_i * norm_j)

            if (abs(cos_sim) >= threshold) {
                edges[[length(edges) + 1]] <- data.frame(
                    from = fi, to = fj, weight = abs(cos_sim),
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    if (length(edges) == 0) return(NULL)
    do.call(rbind, edges)
}


#' Compute DIABLO co-loading edges
#' @noRd
.compute_diablo_coloading_edges <- function(diablo_results, top_features, threshold) {
    if (is.null(diablo_results$loadings)) return(NULL)

    loadings <- diablo_results$loadings
    loadings <- loadings[names(loadings) != "Y"]
    if (length(loadings) < 2) return(NULL)

    # For each component, find features from different omics with abs loading > threshold
    view_names <- names(loadings)
    edges <- list()

    # Determine number of components (min across views)
    n_comp <- min(vapply(loadings, ncol, integer(1)))

    for (comp in seq_len(n_comp)) {
        # Collect features with high loading on this component from each view
        view_feats <- list()
        for (vn in view_names) {
            mat <- as.matrix(loadings[[vn]])
            if (comp > ncol(mat)) next

            selected <- top_features$feature_id[top_features$view == vn]
            common <- intersect(selected, rownames(mat))
            if (length(common) == 0) next

            vals <- mat[common, comp]
            high <- names(vals[abs(vals) >= threshold])
            if (length(high) > 0) {
                view_feats[[vn]] <- data.frame(
                    feature_id = high,
                    loading = vals[high],
                    view = vn,
                    stringsAsFactors = FALSE
                )
            }
        }

        if (length(view_feats) < 2) next

        # Cross-view edges
        view_pairs <- utils::combn(names(view_feats), 2, simplify = FALSE)
        for (vp in view_pairs) {
            df1 <- view_feats[[vp[1]]]
            df2 <- view_feats[[vp[2]]]
            for (r1 in seq_len(nrow(df1))) {
                for (r2 in seq_len(nrow(df2))) {
                    w <- (abs(df1$loading[r1]) + abs(df2$loading[r2])) / 2
                    edges[[length(edges) + 1]] <- data.frame(
                        from = df1$feature_id[r1],
                        to = df2$feature_id[r2],
                        weight = w,
                        stringsAsFactors = FALSE
                    )
                }
            }
        }
    }

    if (length(edges) == 0) return(NULL)
    do.call(rbind, edges)
}


#' Compute expression correlation edges for RNA-protein pairs
#' @noRd
.compute_correlation_edges <- function(harmonization_res, top_features, threshold) {
    mapping <- harmonization_res$gene_protein_mapping
    mae <- harmonization_res$mae

    if (is.null(mapping) || is.null(mae)) return(NULL)
    if (!is.data.frame(mapping) || nrow(mapping) == 0) return(NULL)
    if (!all(c("gene_id", "protein_id") %in% colnames(mapping))) return(NULL)

    # Get experiment names
    exp_names <- names(mae@ExperimentList)
    rna_name <- intersect(c("transcriptomics", "rna"), exp_names)
    prot_name <- intersect(c("proteomics", "protein"), exp_names)

    if (length(rna_name) == 0 || length(prot_name) == 0) return(NULL)
    rna_name <- rna_name[1]
    prot_name <- prot_name[1]

    rna_mat <- as.matrix(SummarizedExperiment::assay(mae@ExperimentList[[rna_name]]))
    prot_mat <- as.matrix(SummarizedExperiment::assay(mae@ExperimentList[[prot_name]]))

    # Common samples
    common_samples <- intersect(colnames(rna_mat), colnames(prot_mat))
    if (length(common_samples) < 5) {
        message("  <5 common samples for correlation edges, skipping")
        return(NULL)
    }

    rna_mat <- rna_mat[, common_samples, drop = FALSE]
    prot_mat <- prot_mat[, common_samples, drop = FALSE]

    # Filter to top features
    rna_features <- top_features$feature_id[top_features$view %in% c("transcriptomics", "rna")]
    prot_features <- top_features$feature_id[top_features$view %in% c("proteomics", "protein")]

    # Filter mapping to selected features
    mapping_sub <- mapping[
        mapping$gene_id %in% intersect(rna_features, rownames(rna_mat)) &
        mapping$protein_id %in% intersect(prot_features, rownames(prot_mat)),
    ]

    if (nrow(mapping_sub) == 0) return(NULL)

    edges <- list()
    for (i in seq_len(nrow(mapping_sub))) {
        gid <- mapping_sub$gene_id[i]
        pid <- mapping_sub$protein_id[i]

        rna_vals <- rna_mat[gid, ]
        prot_vals <- prot_mat[pid, ]

        if (stats::sd(rna_vals, na.rm = TRUE) == 0 ||
            stats::sd(prot_vals, na.rm = TRUE) == 0) next

        r <- stats::cor(rna_vals, prot_vals, use = "pairwise.complete.obs",
                        method = "spearman")
        if (!is.na(r) && abs(r) >= threshold) {
            edges[[length(edges) + 1]] <- data.frame(
                from = gid, to = pid, weight = abs(r),
                stringsAsFactors = FALSE
            )
        }
    }

    if (length(edges) == 0) return(NULL)
    do.call(rbind, edges)
}


# =============================================================================
# Extended Network Functions
# =============================================================================

#' Build Pathway Crosstalk Network
#'
#' Connects KEGG pathways that are co-significant across omics layers.
#' Nodes = pathways, colored by shared/specific category. Edges = co-significance.
#'
#' @param enrichment_res Enrichment results from mod_multiomics_enrichment.
#' @param foundational_res Foundational results from mod_multiomics_foundational.
#' @param config Pipeline configuration.
#' @param out_dir Output directory.
#' @return List with graph, nodes, edges data frames, or NULL.
#' @export
build_pathway_crosstalk_network <- function(enrichment_res, foundational_res,
                                            config, out_dir) {
    message("=== Building Pathway Crosstalk Network ===")

    if (!requireNamespace("igraph", quietly = TRUE)) {
        message("  Package 'igraph' not installed. Skipping.")
        return(NULL)
    }

    net_cfg <- config$modes$multiomics$extended_networks$pathway_crosstalk %||% list()
    max_pathways <- net_cfg$max_pathways %||% 60
    padj_thresh  <- net_cfg$padj_threshold %||% 0.05

    # ---- Collect per-omics pathway significance ----
    per_omics <- enrichment_res$per_omics
    if (is.null(per_omics) || length(per_omics) == 0) {
        message("  No per-omics enrichment data. Skipping pathway crosstalk.")
        return(NULL)
    }

    # Build pathway x omics significance matrix
    pw_sig <- list()
    for (om in names(per_omics)) {
        df <- per_omics[[om]]
        if (!is.data.frame(df) || nrow(df) == 0) next
        if (!"padj" %in% colnames(df)) next
        sig <- df[!is.na(df$padj) & df$padj < padj_thresh, , drop = FALSE]
        if (nrow(sig) == 0) next
        pw_key <- if ("pathway" %in% colnames(sig)) sig$pathway else sig$ID
        pw_sig[[om]] <- unique(pw_key)
    }

    if (length(pw_sig) == 0) {
        message("  No significant pathways found. Skipping pathway crosstalk.")
        return(NULL)
    }

    all_pathways <- unique(unlist(pw_sig))
    if (length(all_pathways) == 0) return(NULL)

    # Classification: which omics is each pathway significant in?
    omics_names <- names(pw_sig)
    pw_omics <- data.frame(pathway = all_pathways, stringsAsFactors = FALSE)
    for (om in omics_names) {
        pw_omics[[om]] <- pw_omics$pathway %in% pw_sig[[om]]
    }
    pw_omics$n_omics <- rowSums(pw_omics[, omics_names, drop = FALSE])

    n_total_omics <- length(omics_names)
    pw_omics$category <- ifelse(
        pw_omics$n_omics == n_total_omics, "shared_all",
        ifelse(pw_omics$n_omics > 1, "shared_some", "omics_specific")
    )

    # Get combined p-values from meta-analysis if available
    meta <- enrichment_res$cross_omics$meta_analysis
    if (!is.null(meta) && is.data.frame(meta) && "pathway" %in% colnames(meta)) {
        pw_omics <- merge(pw_omics, meta[, c("pathway", "combined_pval"), drop = FALSE],
                          by = "pathway", all.x = TRUE)
    } else {
        pw_omics$combined_pval <- NA_real_
    }

    # Also try foundational pathway overlap classification
    fo_class <- tryCatch(
        foundational_res$foundational_results$pathway_overlap$shared_specific$classification,
        error = function(e) NULL
    )
    if (!is.null(fo_class) && is.data.frame(fo_class) && "pathway" %in% colnames(fo_class)) {
        fo_sub <- fo_class[, c("pathway", "category"), drop = FALSE]
        colnames(fo_sub)[2] <- "fo_category"
        pw_omics <- merge(pw_omics, fo_sub, by = "pathway", all.x = TRUE)
        has_fo <- !is.na(pw_omics$fo_category)
        pw_omics$category[has_fo] <- pw_omics$fo_category[has_fo]
        pw_omics$fo_category <- NULL
    }

    # Limit to top pathways
    pw_omics <- pw_omics[order(-pw_omics$n_omics, pw_omics$combined_pval, na.last = TRUE), ]
    pw_omics <- utils::head(pw_omics, max_pathways)

    if (nrow(pw_omics) < 3) {
        message("  Too few pathways (", nrow(pw_omics), "). Skipping pathway crosstalk.")
        return(NULL)
    }

    # ---- Build edges: co-significant pathways within same omics ----
    pw_edges <- list()
    pw_names <- pw_omics$pathway

    for (om in omics_names) {
        om_pws <- pw_names[pw_omics[[om]]]
        if (length(om_pws) < 2) next
        pairs <- utils::combn(om_pws, 2, simplify = FALSE)
        for (p in pairs) {
            key <- paste(sort(p), collapse = "___")
            if (is.null(pw_edges[[key]])) {
                pw_edges[[key]] <- data.frame(
                    from = p[1], to = p[2], weight = 1, stringsAsFactors = FALSE
                )
            } else {
                pw_edges[[key]]$weight <- pw_edges[[key]]$weight + 1
            }
        }
    }

    if (length(pw_edges) == 0) {
        message("  No pathway edges found. Skipping pathway crosstalk.")
        return(NULL)
    }

    edge_df <- do.call(rbind, pw_edges)
    rownames(edge_df) <- NULL

    # ---- Build igraph ----
    g <- igraph::graph_from_data_frame(edge_df[, c("from", "to")], directed = FALSE)
    igraph::E(g)$weight <- edge_df$weight

    node_info <- pw_omics[match(igraph::V(g)$name, pw_omics$pathway), ]

    # Color by category
    view_colors <- c(transcriptomics = "#d73027", proteomics = "#4575b4",
                     metabolomics = "#1b7837")
    node_colors <- character(nrow(node_info))
    for (i in seq_len(nrow(node_info))) {
        cat <- node_info$category[i]
        if (!is.na(cat) && cat == "shared_all") {
            node_colors[i] <- "#FFB300"
        } else if (!is.na(cat) && cat == "shared_some") {
            node_colors[i] <- "#FB8C00"
        } else {
            specific_om <- omics_names[unlist(node_info[i, omics_names, drop = TRUE])]
            if (length(specific_om) == 1 && specific_om %in% names(view_colors)) {
                node_colors[i] <- view_colors[[specific_om]]
            } else {
                node_colors[i] <- "gray70"
            }
        }
    }
    igraph::V(g)$color <- node_colors
    igraph::V(g)$category <- node_info$category

    # Node size by -log10(combined_pval), rescaled 5-18
    log_pval <- -log10(pmax(node_info$combined_pval, 1e-50, na.rm = TRUE))
    log_pval[is.na(log_pval) | is.infinite(log_pval)] <- 5
    mn <- min(log_pval, na.rm = TRUE)
    mx <- max(log_pval, na.rm = TRUE)
    scaled <- if (mx > mn) (log_pval - mn) / (mx - mn) else rep(0.5, length(log_pval))
    igraph::V(g)$size <- 5 + scaled * 13

    # Labels: show top 20
    show_label <- rep(FALSE, length(igraph::V(g)))
    top_idx <- order(node_info$combined_pval, na.last = TRUE)
    show_label[utils::head(top_idx, 20)] <- TRUE
    short_names <- substr(igraph::V(g)$name, 1, 35)
    igraph::V(g)$label <- ifelse(show_label, short_names, NA)

    igraph::E(g)$width <- 0.5 + igraph::E(g)$weight * 0.8
    igraph::E(g)$color <- "gray60"

    # ---- Plot ----
    if (!is.null(out_dir)) {
        grDevices::png(
            file.path(out_dir, "plots", "pathway_crosstalk_network.png"),
            width = 14, height = 14, units = "in", res = 300
        )
        tryCatch({
            layout <- igraph::layout_with_fr(g)
            igraph::plot.igraph(
                g, vertex.size = igraph::V(g)$size, vertex.color = igraph::V(g)$color,
                vertex.label = igraph::V(g)$label, vertex.label.cex = 0.55,
                vertex.label.color = "black", vertex.frame.color = NA,
                edge.color = igraph::E(g)$color, edge.width = igraph::E(g)$width,
                layout = layout, main = "Pathway Crosstalk Network"
            )
            legend("bottomleft",
                   legend = c("Shared (all omics)", "Shared (some omics)",
                              names(view_colors)),
                   col = c("#FFB300", "#FB8C00", unname(view_colors)),
                   pch = 19, pt.cex = 2, bty = "n", title = "Pathway Category")
        }, error = function(e) message("  Pathway crosstalk plot failed: ", e$message))
        grDevices::dev.off()
        message("  Saved pathway_crosstalk_network.png")
    }

    # ---- Tables ----
    node_degrees <- igraph::degree(g)
    nodes_df <- data.frame(
        pathway = igraph::V(g)$name, category = node_info$category,
        n_omics = node_info$n_omics, combined_pval = node_info$combined_pval,
        degree = node_degrees[igraph::V(g)$name], stringsAsFactors = FALSE
    )
    if (!is.null(out_dir)) {
        write.csv(nodes_df, file.path(out_dir, "tables", "pathway_crosstalk_nodes.csv"),
                  row.names = FALSE)
        write.csv(edge_df, file.path(out_dir, "tables", "pathway_crosstalk_edges.csv"),
                  row.names = FALSE)
        message("  Saved pathway crosstalk tables")
    }

    list(graph = g, nodes = nodes_df, edges = edge_df)
}


#' Build Per-Factor Feature Networks
#'
#' For each MOFA factor (or DIABLO component), builds a network of top features
#' across omics views, connected by their shared factor loading.
#'
#' @param integration_res Integration results (MOFA2/DIABLO).
#' @param config Pipeline configuration.
#' @param out_dir Output directory.
#' @return List of per-factor results, or NULL.
#' @export
build_per_factor_networks <- function(integration_res, config, out_dir,
                                      gene_protein_mapping = NULL) {
    message("=== Building Per-Factor Feature Networks ===")

    if (!requireNamespace("igraph", quietly = TRUE)) {
        message("  Package 'igraph' not installed. Skipping.")
        return(NULL)
    }

    net_cfg <- config$modes$multiomics$extended_networks$per_factor %||% list()
    top_n_per_view <- net_cfg$top_n_per_view %||% 25
    max_factors    <- net_cfg$max_factors %||% 5

    # ---- Try MOFA weights ----
    weights <- integration_res$mofa_results$weights
    top_features <- integration_res$mofa_results$top_features
    var_explained <- integration_res$mofa_results$variance_explained

    if (!is.null(weights) && is.list(weights) && length(weights) >= 2 &&
        !is.null(top_features) && is.data.frame(top_features)) {

        factor_col <- "factor"
        view_col <- "view"
        weight_col <- "abs_weight"

        factor_names <- unique(top_features[[factor_col]])
        factor_names <- utils::head(factor_names, max_factors)

        factor_labels <- list()
        if (!is.null(var_explained) && is.data.frame(var_explained)) {
            for (fn in factor_names) {
                if (fn %in% colnames(var_explained)) {
                    total_var <- sum(var_explained[[fn]], na.rm = TRUE)
                    factor_labels[[fn]] <- sprintf("%s (%.1f%% variance)", fn, total_var)
                } else {
                    factor_labels[[fn]] <- fn
                }
            }
        } else {
            for (fn in factor_names) factor_labels[[fn]] <- fn
        }

    } else if (!is.null(integration_res$diablo_results$loadings)) {
        loadings <- integration_res$diablo_results$loadings
        loadings <- loadings[names(loadings) != "Y"]
        if (length(loadings) < 2) {
            message("  No multi-view loadings found. Skipping per-factor networks.")
            return(NULL)
        }

        n_comp <- min(vapply(loadings, ncol, integer(1)))
        factor_names <- paste0("Component", seq_len(min(n_comp, max_factors)))
        rows <- list()
        for (vn in names(loadings)) {
            mat <- as.matrix(loadings[[vn]])
            for (comp in seq_len(min(n_comp, max_factors))) {
                vals <- mat[, comp]
                ord <- order(abs(vals), decreasing = TRUE)
                sel <- utils::head(ord, top_n_per_view)
                rows[[paste0(vn, "_", comp)]] <- data.frame(
                    feature_id = rownames(mat)[sel], weight = vals[sel],
                    abs_weight = abs(vals[sel]),
                    factor = paste0("Component", comp), view = vn,
                    stringsAsFactors = FALSE
                )
            }
        }
        top_features <- do.call(rbind, rows)
        rownames(top_features) <- NULL
        factor_col <- "factor"
        view_col <- "view"
        weight_col <- "abs_weight"
        factor_labels <- list()
        for (fn in factor_names) factor_labels[[fn]] <- fn
    } else {
        message("  No MOFA or DIABLO results found. Skipping per-factor networks.")
        return(NULL)
    }

    view_colors <- c(transcriptomics = "#d73027", proteomics = "#4575b4",
                     metabolomics = "#1b7837")

    all_nodes <- list()
    all_edges <- list()
    factor_results <- list()

    for (fn in factor_names) {
        fdf <- top_features[top_features[[factor_col]] == fn, , drop = FALSE]
        if (nrow(fdf) == 0) next

        fdf_sel <- do.call(rbind, lapply(split(fdf, fdf[[view_col]]), function(vdf) {
            vdf <- vdf[order(vdf[[weight_col]], decreasing = TRUE), ]
            utils::head(vdf, top_n_per_view)
        }))
        rownames(fdf_sel) <- NULL

        views_present <- unique(fdf_sel[[view_col]])
        if (length(views_present) < 2) next

        # Build cross-view edges
        factor_edges <- list()
        view_pairs <- utils::combn(views_present, 2, simplify = FALSE)
        for (vp in view_pairs) {
            df1 <- fdf_sel[fdf_sel[[view_col]] == vp[1], , drop = FALSE]
            df2 <- fdf_sel[fdf_sel[[view_col]] == vp[2], , drop = FALSE]
            for (r1 in seq_len(nrow(df1))) {
                for (r2 in seq_len(nrow(df2))) {
                    w1 <- df1[[weight_col]][r1]
                    w2 <- df2[[weight_col]][r2]
                    geo_mean <- sqrt(w1 * w2)
                    if (geo_mean > 0) {
                        factor_edges[[length(factor_edges) + 1]] <- data.frame(
                            from = df1$feature_id[r1], to = df2$feature_id[r2],
                            weight = geo_mean, stringsAsFactors = FALSE
                        )
                    }
                }
            }
        }

        if (length(factor_edges) == 0) next
        f_edge_df <- do.call(rbind, factor_edges)
        f_edge_df <- f_edge_df[order(f_edge_df$weight, decreasing = TRUE), ]
        f_edge_df <- utils::head(f_edge_df, 500)

        g <- igraph::graph_from_data_frame(f_edge_df[, c("from", "to")], directed = FALSE)

        node_match <- fdf_sel[match(igraph::V(g)$name, fdf_sel$feature_id), ]
        node_views <- node_match[[view_col]]
        node_views[is.na(node_views)] <- "unknown"
        igraph::V(g)$color <- vapply(
            node_views,
            function(v) if (v %in% names(view_colors)) view_colors[[v]] else "gray60",
            character(1)
        )
        igraph::V(g)$view <- node_views

        raw_weight <- node_match$weight
        raw_weight[is.na(raw_weight)] <- 0
        igraph::V(g)$shape <- ifelse(raw_weight >= 0, "circle", "square")

        abs_w <- node_match[[weight_col]]
        abs_w[is.na(abs_w)] <- 0
        mn <- min(abs_w, na.rm = TRUE)
        mx <- max(abs_w, na.rm = TRUE)
        scaled <- if (mx > mn) (abs_w - mn) / (mx - mn) else rep(0.5, length(abs_w))
        igraph::V(g)$size <- 4 + scaled * 11

        # Resolve GENE_N to original IDs
        resolved <- .resolve_feature_ids(igraph::V(g)$name, igraph::V(g)$view,
                                         gene_protein_mapping)
        igraph::V(g)$original_id <- resolved
        igraph::V(g)$label <- resolved
        igraph::V(g)$label.cex <- 0.5

        # Plot
        if (!is.null(out_dir)) {
            fn_clean <- gsub("[^a-zA-Z0-9]", "_", fn)
            grDevices::png(
                file.path(out_dir, "plots",
                          paste0("factor_", fn_clean, "_feature_network.png")),
                width = 10, height = 10, units = "in", res = 300
            )
            tryCatch({
                layout <- igraph::layout_with_fr(g)
                title <- factor_labels[[fn]] %||% fn
                igraph::plot.igraph(
                    g, vertex.size = igraph::V(g)$size,
                    vertex.color = igraph::V(g)$color,
                    vertex.shape = igraph::V(g)$shape,
                    vertex.label = igraph::V(g)$label,
                    vertex.label.cex = 0.5, vertex.label.color = "black",
                    vertex.frame.color = NA,
                    edge.color = "gray75", edge.width = 0.8,
                    layout = layout, main = title
                )
                avail_views <- intersect(views_present, names(view_colors))
                if (length(avail_views) > 0) {
                    legend("bottomleft", legend = avail_views,
                           col = view_colors[avail_views],
                           pch = 19, pt.cex = 2, bty = "n", title = "Omics Layer")
                }
                legend("bottomright",
                       legend = c("Positive weight", "Negative weight"),
                       pch = c(21, 22), pt.bg = "gray80",
                       pt.cex = 2, bty = "n", title = "Weight Sign")
            }, error = function(e) message("  Factor ", fn, " plot failed: ", e$message))
            grDevices::dev.off()
            message("  Saved factor_", fn_clean, "_feature_network.png")
        }

        node_deg <- igraph::degree(g)
        nodes_out <- data.frame(
            feature_id = igraph::V(g)$name,
            original_id = igraph::V(g)$original_id,
            view = igraph::V(g)$view,
            factor = fn, abs_weight = abs_w,
            degree = node_deg[igraph::V(g)$name], stringsAsFactors = FALSE
        )
        # Resolve edge IDs
        f_edge_df$from_original <- .resolve_feature_ids(
            f_edge_df$from, rep(NA, nrow(f_edge_df)), gene_protein_mapping)
        f_edge_df$to_original <- .resolve_feature_ids(
            f_edge_df$to, rep(NA, nrow(f_edge_df)), gene_protein_mapping)
        f_edge_df$factor <- fn

        all_nodes[[fn]] <- nodes_out
        all_edges[[fn]] <- f_edge_df
        factor_results[[fn]] <- list(graph = g, nodes = nodes_out, edges = f_edge_df)
    }

    if (length(factor_results) == 0) {
        message("  No multi-view factors found. Skipping per-factor networks.")
        return(NULL)
    }

    if (!is.null(out_dir)) {
        combined_nodes <- do.call(rbind, all_nodes)
        combined_edges <- do.call(rbind, all_edges)
        rownames(combined_nodes) <- NULL
        rownames(combined_edges) <- NULL
        write.csv(combined_nodes, file.path(out_dir, "tables", "factor_feature_nodes.csv"),
                  row.names = FALSE)
        write.csv(combined_edges, file.path(out_dir, "tables", "factor_feature_edges.csv"),
                  row.names = FALSE)
        message("  Saved factor feature tables")
    }

    factor_results
}


#' Build Cross-Omics Module Network
#'
#' Visualizes WGCNA cross-omics modules showing features from different omics
#' that co-vary together.
#'
#' @param foundational_res Foundational results from mod_multiomics_foundational.
#' @param config Pipeline configuration.
#' @param out_dir Output directory.
#' @return List with graph, nodes, edges data frames, or NULL.
#' @export
build_crossomics_module_network <- function(foundational_res, config, out_dir,
                                             gene_protein_mapping = NULL) {
    message("=== Building Cross-Omics Module Network ===")

    if (!requireNamespace("igraph", quietly = TRUE)) {
        message("  Package 'igraph' not installed. Skipping.")
        return(NULL)
    }

    net_cfg <- config$modes$multiomics$extended_networks$crossomics_modules %||% list()
    max_modules            <- net_cfg$max_modules %||% 5
    max_features_per_omics <- net_cfg$max_features_per_omics %||% 30
    eigengene_cor_thresh   <- net_cfg$eigengene_cor_threshold %||% 0.5

    modules <- tryCatch(
        foundational_res$foundational_results$crossomics_modules,
        error = function(e) NULL
    )
    if (is.null(modules)) {
        message("  No cross-omics module data found. Skipping.")
        return(NULL)
    }

    mod_assign <- modules$module_assignments
    mod_eigen  <- modules$module_eigengenes

    if (!is.data.frame(mod_assign) || nrow(mod_assign) == 0) {
        message("  No module assignments found. Skipping.")
        return(NULL)
    }

    required_cols <- c("feature_id", "module_number", "module_color", "omics")
    if (!all(required_cols %in% colnames(mod_assign))) {
        message("  Module assignments missing required columns. Skipping.")
        return(NULL)
    }

    # Identify cross-omics modules
    mod_summary <- stats::aggregate(
        omics ~ module_number, data = mod_assign,
        FUN = function(x) length(unique(x))
    )
    colnames(mod_summary)[2] <- "n_omics_types"
    crossomics_mods <- mod_summary$module_number[mod_summary$n_omics_types >= 2]

    if (length(crossomics_mods) == 0) {
        message("  No cross-omics modules found. Skipping.")
        return(NULL)
    }

    mod_sizes <- table(mod_assign$module_number[mod_assign$module_number %in% crossomics_mods])
    top_mods <- as.integer(names(sort(mod_sizes, decreasing = TRUE)))
    top_mods <- utils::head(top_mods, max_modules)

    # Subset features
    sub_assign <- mod_assign[mod_assign$module_number %in% top_mods, , drop = FALSE]
    sub_list <- list()
    for (mn in top_mods) {
        mod_df <- sub_assign[sub_assign$module_number == mn, , drop = FALSE]
        for (om in unique(mod_df$omics)) {
            om_df <- mod_df[mod_df$omics == om, , drop = FALSE]
            sub_list[[paste0(mn, "_", om)]] <- utils::head(om_df, max_features_per_omics)
        }
    }
    sub_assign <- do.call(rbind, sub_list)
    rownames(sub_assign) <- NULL

    if (nrow(sub_assign) < 4) {
        message("  Too few features for module network. Skipping.")
        return(NULL)
    }

    # Build edges
    mod_edges <- list()

    # Intra-module: connect cross-omics features within same module
    for (mn in top_mods) {
        mod_feats <- sub_assign[sub_assign$module_number == mn, , drop = FALSE]
        omics_in_mod <- unique(mod_feats$omics)
        if (length(omics_in_mod) < 2) next

        pairs <- utils::combn(omics_in_mod, 2, simplify = FALSE)
        for (vp in pairs) {
            f1 <- mod_feats$feature_id[mod_feats$omics == vp[1]]
            f2 <- mod_feats$feature_id[mod_feats$omics == vp[2]]
            f1_top <- utils::head(f1, 5)
            f2_top <- utils::head(f2, 5)
            for (a in f1_top) {
                for (b in f2_top) {
                    mod_edges[[length(mod_edges) + 1]] <- data.frame(
                        from = a, to = b, weight = 0.5,
                        module = mn, type = "intra", stringsAsFactors = FALSE
                    )
                }
            }
        }
    }

    # Inter-module edges via eigengene correlation
    if (!is.null(mod_eigen) && ncol(mod_eigen) >= 2) {
        mod_cols <- colnames(mod_eigen)
        for (i in seq_len(length(top_mods) - 1)) {
            for (j in (i + 1):length(top_mods)) {
                col_i <- paste0("module_", top_mods[i])
                col_j <- paste0("module_", top_mods[j])
                if (!(col_i %in% mod_cols && col_j %in% mod_cols)) next
                r <- stats::cor(mod_eigen[, col_i], mod_eigen[, col_j],
                               use = "pairwise.complete.obs")
                if (!is.na(r) && abs(r) > eigengene_cor_thresh) {
                    feats_i <- sub_assign$feature_id[sub_assign$module_number == top_mods[i]]
                    feats_j <- sub_assign$feature_id[sub_assign$module_number == top_mods[j]]
                    if (length(feats_i) > 0 && length(feats_j) > 0) {
                        mod_edges[[length(mod_edges) + 1]] <- data.frame(
                            from = feats_i[1], to = feats_j[1],
                            weight = abs(r), module = 0, type = "inter",
                            stringsAsFactors = FALSE
                        )
                    }
                }
            }
        }
    }

    if (length(mod_edges) == 0) {
        message("  No module edges built. Skipping.")
        return(NULL)
    }

    edge_df <- do.call(rbind, mod_edges)
    rownames(edge_df) <- NULL
    edge_df <- edge_df[edge_df$from != edge_df$to, , drop = FALSE]
    edge_df$pair_key <- apply(edge_df[, c("from", "to")], 1,
                              function(x) paste(sort(x), collapse = "___"))
    edge_df <- edge_df[!duplicated(edge_df$pair_key), ]
    edge_df$pair_key <- NULL

    if (nrow(edge_df) == 0) {
        message("  No valid edges after dedup. Skipping.")
        return(NULL)
    }

    # Build igraph
    g <- igraph::graph_from_data_frame(edge_df[, c("from", "to")], directed = FALSE)
    valid_nodes <- igraph::V(g)$name %in% sub_assign$feature_id
    if (!all(valid_nodes)) {
        g <- igraph::induced_subgraph(g, igraph::V(g)$name[valid_nodes])
    }

    node_info <- sub_assign[match(igraph::V(g)$name, sub_assign$feature_id), ]
    view_colors <- c(transcriptomics = "#d73027", proteomics = "#4575b4",
                     metabolomics = "#1b7837")
    node_omics <- node_info$omics
    node_omics[is.na(node_omics)] <- "unknown"
    igraph::V(g)$color <- vapply(
        node_omics,
        function(v) if (v %in% names(view_colors)) view_colors[[v]] else "gray60",
        character(1)
    )
    igraph::V(g)$view <- node_omics

    mod_colors <- node_info$module_color
    mod_colors[is.na(mod_colors)] <- "gray50"
    igraph::V(g)$frame.color <- mod_colors
    igraph::V(g)$size <- 6

    # Resolve GENE_N to original IDs
    resolved <- .resolve_feature_ids(igraph::V(g)$name, igraph::V(g)$view,
                                     gene_protein_mapping)
    igraph::V(g)$original_id <- resolved

    # Labels: top 10 highest-degree per module
    node_deg <- igraph::degree(g)
    show_label <- rep(FALSE, length(igraph::V(g)))
    for (mn in top_mods) {
        mod_idx <- which(node_info$module_number == mn)
        if (length(mod_idx) == 0) next
        deg_order <- order(node_deg[mod_idx], decreasing = TRUE)
        top_idx <- mod_idx[utils::head(deg_order, 10)]
        show_label[top_idx] <- TRUE
    }
    igraph::V(g)$label <- ifelse(show_label, resolved, NA)

    # Edge appearance
    edge_match <- edge_df[match(
        apply(igraph::as_edgelist(g), 1, function(x) paste(sort(x), collapse = "___")),
        apply(edge_df[, c("from", "to")], 1, function(x) paste(sort(x), collapse = "___"))
    ), ]
    edge_types <- edge_match$type
    edge_types[is.na(edge_types)] <- "intra"
    igraph::E(g)$color <- ifelse(edge_types == "inter", "gray30", "gray70")
    igraph::E(g)$lty <- ifelse(edge_types == "inter", 2, 1)

    # Plot
    if (!is.null(out_dir)) {
        grDevices::png(
            file.path(out_dir, "plots", "crossomics_module_network.png"),
            width = 14, height = 14, units = "in", res = 300
        )
        tryCatch({
            layout <- igraph::layout_with_fr(g)
            igraph::plot.igraph(
                g, vertex.size = igraph::V(g)$size,
                vertex.color = igraph::V(g)$color,
                vertex.frame.color = igraph::V(g)$frame.color,
                vertex.frame.width = 2,
                vertex.label = igraph::V(g)$label,
                vertex.label.cex = 0.5, vertex.label.color = "black",
                edge.color = igraph::E(g)$color,
                edge.lty = igraph::E(g)$lty, edge.width = 1,
                layout = layout, main = "Cross-Omics Module Network (WGCNA)"
            )
            avail_views <- intersect(unique(node_omics), names(view_colors))
            if (length(avail_views) > 0) {
                legend("bottomleft", legend = avail_views,
                       col = view_colors[avail_views],
                       pch = 19, pt.cex = 2, bty = "n", title = "Omics Layer")
            }
            mod_colors_unique <- unique(node_info$module_color[!is.na(node_info$module_color)])
            if (length(mod_colors_unique) > 0) {
                mod_nums <- vapply(mod_colors_unique, function(mc) {
                    paste0("Module ", node_info$module_number[node_info$module_color == mc][1])
                }, character(1))
                legend("bottomright", legend = mod_nums,
                       col = mod_colors_unique, pch = 0, pt.cex = 2,
                       lwd = 2, bty = "n", title = "Module (border)")
            }
        }, error = function(e) message("  Module network plot failed: ", e$message))
        grDevices::dev.off()
        message("  Saved crossomics_module_network.png")
    }

    # Tables
    nodes_df <- data.frame(
        feature_id = igraph::V(g)$name,
        original_id = igraph::V(g)$original_id,
        omics = igraph::V(g)$view,
        module_number = node_info$module_number, module_color = node_info$module_color,
        degree = node_deg[igraph::V(g)$name], stringsAsFactors = FALSE
    )
    if (!is.null(out_dir)) {
        write.csv(nodes_df, file.path(out_dir, "tables", "module_network_nodes.csv"),
                  row.names = FALSE)
        write.csv(edge_df, file.path(out_dir, "tables", "module_network_edges.csv"),
                  row.names = FALSE)
        message("  Saved module network tables")
    }

    list(graph = g, nodes = nodes_df, edges = edge_df)
}


#' Include Mechanistic Regulatory Network
#'
#' Copies pre-computed mechanistic regulatory network plot and hub regulators
#' table into the feature_network output directory for unified reporting.
#'
#' @param mechanistic_res Mechanistic analysis results.
#' @param out_dir Output directory (feature_network dir).
#' @return List with paths of copied files, or NULL.
#' @export
include_mechanistic_network <- function(mechanistic_res, out_dir) {
    message("=== Including Mechanistic Regulatory Network ===")

    if (is.null(out_dir)) return(NULL)

    mech_dir <- file.path(dirname(out_dir), "mechanistic", "mechanistic")
    src_plot <- file.path(mech_dir, "plots", "mech_regulatory_network.png")
    src_hubs <- file.path(mech_dir, "tables", "mech_hub_regulators.csv")

    copied <- list()

    if (file.exists(src_plot)) {
        dst <- file.path(out_dir, "plots", "mech_regulatory_network.png")
        file.copy(src_plot, dst, overwrite = TRUE)
        copied$network_plot <- dst
        message("  Copied mech_regulatory_network.png")
    } else {
        message("  mech_regulatory_network.png not found at: ", src_plot)
    }

    if (file.exists(src_hubs)) {
        dst <- file.path(out_dir, "tables", "mech_hub_regulators.csv")
        file.copy(src_hubs, dst, overwrite = TRUE)
        copied$hub_regulators <- dst
        message("  Copied mech_hub_regulators.csv")
    } else {
        message("  mech_hub_regulators.csv not found at: ", src_hubs)
    }

    if (length(copied) == 0) {
        message("  No mechanistic files found to copy.")
        return(NULL)
    }

    copied
}


#' Build Extended Networks (Wrapper)
#'
#' Calls all extended network functions (pathway crosstalk, per-factor,
#' cross-omics modules, mechanistic) in tryCatch blocks.
#'
#' @param enrichment_res Enrichment results.
#' @param foundational_res Foundational analysis results.
#' @param mechanistic_res Mechanistic analysis results.
#' @param integration_res Integration results (MOFA2/DIABLO).
#' @param harmonization_res Harmonization results.
#' @param config Pipeline configuration.
#' @param out_dir Output directory (feature_network dir).
#' @return List of results from each network function.
#' @export
build_extended_networks <- function(enrichment_res, foundational_res,
                                    mechanistic_res, integration_res,
                                    harmonization_res, config, out_dir) {
    message("=== Building Extended Multi-Omics Networks ===")

    if (!is.null(out_dir)) {
        dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    }

    results <- list()
    gpm <- harmonization_res$gene_protein_mapping

    results$pathway_crosstalk <- tryCatch(
        build_pathway_crosstalk_network(
            enrichment_res = enrichment_res, foundational_res = foundational_res,
            config = config, out_dir = out_dir
        ),
        error = function(e) { warning("Pathway crosstalk network: ", e$message); NULL }
    )

    results$per_factor <- tryCatch(
        build_per_factor_networks(
            integration_res = integration_res, config = config, out_dir = out_dir,
            gene_protein_mapping = gpm
        ),
        error = function(e) { warning("Per-factor networks: ", e$message); NULL }
    )

    results$crossomics_modules <- tryCatch(
        build_crossomics_module_network(
            foundational_res = foundational_res, config = config, out_dir = out_dir,
            gene_protein_mapping = gpm
        ),
        error = function(e) { warning("Cross-omics module network: ", e$message); NULL }
    )

    results$mechanistic <- tryCatch(
        include_mechanistic_network(
            mechanistic_res = mechanistic_res, out_dir = out_dir
        ),
        error = function(e) { warning("Mechanistic network inclusion: ", e$message); NULL }
    )

    n_success <- sum(!vapply(results, is.null, logical(1)))
    message("=== Extended Networks Complete: ", n_success, "/4 built ===")

    results
}
