# =============================================================================
# COSMOS: Causal Oriented Search of Multi-Omics Space
# =============================================================================
#
# Uses the cosmosR Bioconductor package to find directed mechanistic subnetworks
# connecting transcription factors, kinases, and metabolites through a prior
# knowledge network (OmniPath + STITCH + Recon3D, ~117K edges).
#
# COSMOS uses CARNIVAL's Integer Linear Programming (ILP) to find the smallest
# coherent subnetwork causally connecting deregulated molecular features.
#
# Requirements (not installed by default):
#   BiocManager::install(c("cosmosR", "CARNIVAL", "decoupleR"))
#   install.packages("lpSolve")  # ILP solver backend
#
# Reference:
#   Dugourd et al. (2021) Mol Syst Biol 17:e9730
#   https://saezlab.github.io/cosmosR/
# =============================================================================


#' Run COSMOS causal network analysis
#'
#' Estimates TF activities from transcriptomics DE results, maps metabolite
#' features to HMDB/KEGG IDs, and uses COSMOS to find a directed subnetwork
#' connecting deregulated TFs and metabolites through known signaling and
#' metabolic interactions.
#'
#' @param mae MultiAssayExperiment object
#' @param de_results Named list of DE results per omics
#' @param config Full config object
#' @param out_dir Output directory
#' @return List with network edges, nodes, and plot paths, or NULL on failure
run_cosmos_causal_network <- function(mae, de_results = NULL, config,
                                      out_dir = NULL) {

    message("Running COSMOS causal network analysis...")

    # Check dependencies
    if (!requireNamespace("cosmosR", quietly = TRUE)) {
        message("  cosmosR not installed. Install with: BiocManager::install('cosmosR')")
        return(NULL)
    }
    if (!requireNamespace("CARNIVAL", quietly = TRUE)) {
        message("  CARNIVAL not installed. Install with: BiocManager::install('CARNIVAL')")
        return(NULL)
    }
    if (!requireNamespace("decoupleR", quietly = TRUE)) {
        message("  decoupleR not installed. Install with: BiocManager::install('decoupleR')")
        return(NULL)
    }

    # Check for a solver backend
    has_solver <- requireNamespace("lpSolve", quietly = TRUE) ||
                  requireNamespace("cbc", quietly = TRUE)
    if (!has_solver) {
        message("  No ILP solver found. Install lpSolve: install.packages('lpSolve')")
        return(NULL)
    }

    # --- Step 1: Extract TF activities from transcriptomics DE ---
    tf_activities <- estimate_tf_activities_for_cosmos(mae, de_results)
    if (is.null(tf_activities) || length(tf_activities) == 0) {
        message("  Could not estimate TF activities. Skipping COSMOS.")
        return(NULL)
    }
    message("  Estimated activities for ", length(tf_activities), " TFs")

    # --- Step 2: Extract metabolite measurements ---
    metab_measurements <- extract_metabolite_measurements(mae, de_results)
    if (is.null(metab_measurements) || length(metab_measurements) == 0) {
        message("  Could not extract metabolite measurements. Skipping COSMOS.")
        return(NULL)
    }
    message("  Mapped ", length(metab_measurements), " metabolites to COSMOS IDs")

    # --- Step 3: Load prior knowledge network ---
    message("  Loading COSMOS prior knowledge network...")
    pkn <- tryCatch({
        cosmosR::load_tf_regulon_dorothea()
    }, error = function(e) {
        message("  Failed to load PKN: ", e$message)
        NULL
    })

    # Use the full meta-PKN from cosmosR
    meta_network <- tryCatch({
        cosmosR::meta_network
    }, error = function(e) {
        message("  Could not load meta_network from cosmosR: ", e$message)
        NULL
    })

    if (is.null(meta_network)) {
        message("  COSMOS meta-network not available. Skipping.")
        return(NULL)
    }
    message("  Meta-network: ", nrow(meta_network), " edges")

    # --- Step 4: Run COSMOS signaling to metabolism ---
    message("  Running COSMOS optimization (this may take a few minutes)...")

    cosmos_result <- tryCatch({
        # Preprocess inputs
        cosmos_input <- cosmosR::preprocess_COSMOS_signaling_to_metabolism(
            meta_network = meta_network,
            signaling_data = tf_activities,
            metabolic_data = metab_measurements,
            diff_expression_data = extract_de_t_values(mae, de_results),
            maximum_network_depth = 8,
            remove_unexpressed_nodes = TRUE,
            filter_tf_gene_interaction_by_tf = TRUE
        )

        # Run COSMOS
        cosmos_res <- cosmosR::run_COSMOS_signaling_to_metabolism(
            cosmos_data = cosmos_input,
            max_time = 300,  # 5 minute time limit
            solver_name = if (requireNamespace("lpSolve", quietly = TRUE)) "lpSolve" else "cbc"
        )

        cosmos_res
    }, error = function(e) {
        message("  COSMOS optimization failed: ", e$message)
        NULL
    })

    if (is.null(cosmos_result)) return(NULL)

    # --- Step 5: Extract and format results ---
    result <- format_cosmos_results(cosmos_result, out_dir)
    message("  COSMOS complete: ", result$n_nodes, " nodes, ", result$n_edges, " edges")

    result
}


#' Estimate TF activities for COSMOS input using decoupleR
#'
#' @param mae MultiAssayExperiment
#' @param de_results DE results list
#' @return Named numeric vector of TF activity scores (t-values)
estimate_tf_activities_for_cosmos <- function(mae, de_results) {

    # Get transcriptomics DE t-values or log2FC
    if (is.null(de_results$transcriptomics)) {
        message("  No transcriptomics DE results for TF estimation.")
        return(NULL)
    }

    # Extract a named vector: gene -> statistic (t-value or logFC)
    de <- de_results$transcriptomics
    de_vec <- NULL

    if (!is.null(de$tables) && is.list(de$tables)) {
        tbl <- de$tables[[1]]
        stat_col <- intersect(c("stat", "t", "log2FoldChange", "logFC"), names(tbl))[1]
        id_col <- intersect(c("FeatureID", "gene_id", "feature_id"), names(tbl))[1]
        if (!is.na(stat_col)) {
            ids <- if (is.na(id_col)) rownames(tbl) else tbl[[id_col]]
            de_vec <- setNames(tbl[[stat_col]], ids)
        }
    } else if (!is.null(de$summary_df)) {
        df <- de$summary_df
        lfc_cols <- grep("^(logFC\\.|log2FoldChange\\.|logFC_)", names(df), value = TRUE)
        id_col <- intersect(c("FeatureID", "feature_id"), names(df))[1]
        if (length(lfc_cols) > 0) {
            ids <- if (is.na(id_col)) rownames(df) else df[[id_col]]
            de_vec <- setNames(df[[lfc_cols[1]]], ids)
        }
    }

    if (is.null(de_vec) || length(de_vec) == 0) return(NULL)
    de_vec <- de_vec[!is.na(de_vec)]

    # Run decoupleR to get TF activities
    tryCatch({
        # Load dorothea regulons via decoupleR
        regulons <- decoupleR::get_dorothea(organism = "human", levels = c("A", "B", "C"))

        # Convert de_vec to a matrix (1 column)
        mat <- matrix(de_vec, ncol = 1, dimnames = list(names(de_vec), "contrast"))

        # Run weighted mean activity estimation
        tf_acts <- decoupleR::run_wmean(
            mat = mat,
            net = regulons,
            .source = "source",
            .target = "target",
            .mor = "mor",
            times = 100,
            minsize = 5
        )

        # Extract scores for the "wmean" statistic
        tf_scores <- tf_acts[tf_acts$statistic == "norm_wmean", ]
        result <- setNames(tf_scores$score, tf_scores$source)
        # Filter to significant TFs (|score| > 1.5)
        result <- result[abs(result) > 1.5]
        result
    }, error = function(e) {
        message("  decoupleR TF activity estimation failed: ", e$message)
        NULL
    })
}


#' Extract metabolite measurements for COSMOS
#'
#' Maps metabolite feature IDs to HMDB or KEGG compound IDs using rowData,
#' then extracts log2FC values as the measurement vector.
#'
#' @param mae MultiAssayExperiment
#' @param de_results DE results list
#' @return Named numeric vector: HMDB/KEGG ID -> log2FC
extract_metabolite_measurements <- function(mae, de_results) {

    if (!"metabolomics" %in% names(mae@ExperimentList)) return(NULL)

    # Get metabolite rowData for ID mapping
    rd <- as.data.frame(SummarizedExperiment::rowData(mae[["metabolomics"]]))

    # Look for HMDB or KEGG compound IDs in rowData
    hmdb_col <- intersect(c("HMDB", "HMDB_ID", "hmdb_id", "hmdb"), colnames(rd))[1]
    kegg_col <- intersect(c("KEGG", "KEGG_ID", "kegg_id", "Compound_ID"), colnames(rd))[1]

    id_col <- if (!is.na(hmdb_col)) hmdb_col else if (!is.na(kegg_col)) kegg_col else NA

    if (is.na(id_col)) {
        message("  No HMDB or KEGG compound IDs in metabolomics rowData. ",
                "Available columns: ", paste(colnames(rd), collapse = ", "))
        return(NULL)
    }

    # Get log2FC for metabolites
    metab_de <- de_results$metabolomics
    if (is.null(metab_de)) return(NULL)

    de_norm <- tryCatch(
        extract_de_contrast(metab_de, 1, "metabolomics"),
        error = function(e) NULL
    )
    if (is.null(de_norm)) return(NULL)

    # Map feature_id -> compound_id -> logFC
    feature_ids <- de_norm$feature_id
    matched_rows <- match(feature_ids, rownames(rd))
    compound_ids <- rd[[id_col]][matched_rows]

    # Build named vector: compound_id -> logFC
    valid <- !is.na(compound_ids) & compound_ids != "" & !is.na(de_norm$logFC)
    if (sum(valid) == 0) return(NULL)

    result <- setNames(de_norm$logFC[valid], compound_ids[valid])

    # For COSMOS, HMDB IDs should be prefixed with "Metab__HMDB:"
    # and KEGG compound IDs with "Metab__"
    if (!is.na(hmdb_col)) {
        names(result) <- paste0("Metab__HMDB:", names(result))
    } else {
        names(result) <- paste0("Metab__", names(result))
    }

    # Remove duplicates (keep first)
    result <- result[!duplicated(names(result))]
    result
}


#' Extract DE t-values for COSMOS diff_expression_data input
#'
#' @param mae MultiAssayExperiment
#' @param de_results DE results list
#' @return Named numeric vector of gene-level t-values or logFC
extract_de_t_values <- function(mae, de_results) {
    if (is.null(de_results$transcriptomics)) return(numeric(0))

    de <- de_results$transcriptomics
    if (!is.null(de$tables) && is.list(de$tables)) {
        tbl <- de$tables[[1]]
        stat_col <- intersect(c("stat", "t", "log2FoldChange", "logFC"), names(tbl))[1]
        id_col <- intersect(c("FeatureID", "gene_id", "feature_id"), names(tbl))[1]
        if (is.na(stat_col)) return(numeric(0))
        ids <- if (is.na(id_col)) rownames(tbl) else tbl[[id_col]]
        result <- setNames(tbl[[stat_col]], ids)
        return(result[!is.na(result)])
    }
    numeric(0)
}


#' Format COSMOS results into edges, nodes, and plots
#'
#' @param cosmos_result Output from run_COSMOS_signaling_to_metabolism()
#' @param out_dir Output directory
#' @return List with edges, nodes, n_nodes, n_edges, plot paths
format_cosmos_results <- function(cosmos_result, out_dir = NULL) {

    # Extract the solution network
    sif <- tryCatch({
        cosmosR::format_COSMOS_res(cosmos_result)
    }, error = function(e) {
        message("  Could not format COSMOS result: ", e$message)
        NULL
    })

    if (is.null(sif) || nrow(sif) == 0) {
        return(list(n_nodes = 0, n_edges = 0, edges = NULL, nodes = NULL))
    }

    edges <- sif
    nodes <- unique(c(sif$source, sif$target))

    # Classify node types
    node_df <- data.frame(
        node = nodes,
        type = ifelse(grepl("^Metab__", nodes), "metabolite",
               ifelse(grepl("^Enzyme", nodes), "enzyme", "gene/TF")),
        stringsAsFactors = FALSE
    )

    # Extract activity scores if available
    if (!is.null(cosmos_result$activity_scores)) {
        node_df$activity <- cosmos_result$activity_scores[node_df$node]
    }

    result <- list(
        edges = edges,
        nodes = node_df,
        n_nodes = length(nodes),
        n_edges = nrow(edges)
    )

    if (!is.null(out_dir)) {
        dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

        write.csv(edges,
                  file.path(out_dir, "tables", "cosmos_network_edges.csv"),
                  row.names = FALSE)
        write.csv(node_df,
                  file.path(out_dir, "tables", "cosmos_network_nodes.csv"),
                  row.names = FALSE)

        # Plot network
        result$plots <- plot_cosmos_network(edges, node_df, out_dir)
    }

    result
}


#' Plot COSMOS causal network
#'
#' @param edges Edge data frame from format_cosmos_results
#' @param node_df Node data frame
#' @param out_dir Output directory
#' @return List of plot file paths
plot_cosmos_network <- function(edges, node_df, out_dir) {

    if (!requireNamespace("igraph", quietly = TRUE)) {
        message("  igraph not available for COSMOS network plot")
        return(NULL)
    }

    plots <- list()

    # Build igraph
    g <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = node_df)

    # Color by node type
    type_colors <- c(
        "gene/TF" = "#E64B35",
        "metabolite" = "#4DAF4A",
        "enzyme" = "#377EB8"
    )
    node_colors <- type_colors[igraph::V(g)$type]
    node_colors[is.na(node_colors)] <- "grey70"

    # Simplify labels (strip prefixes)
    labels <- igraph::V(g)$name
    labels <- sub("^Metab__HMDB:", "", labels)
    labels <- sub("^Metab__", "", labels)
    labels <- sub("^Enzyme__", "", labels)

    # Truncate long labels
    labels <- ifelse(nchar(labels) > 15,
                     paste0(substr(labels, 1, 12), "..."), labels)

    # Edge sign coloring
    edge_colors <- ifelse(edges$interaction > 0, "#E64B35", "#4DBBD5")
    edge_colors[is.na(edge_colors)] <- "grey50"

    # Layout
    n_nodes <- igraph::vcount(g)
    layout <- if (n_nodes > 50) {
        igraph::layout_with_fr(g)
    } else {
        igraph::layout_with_kk(g)
    }

    # Save plot
    out_png <- file.path(out_dir, "plots", "cosmos_causal_network.png")
    png(out_png, width = 14, height = 10, units = "in", res = 150)
    par(mar = c(1, 1, 3, 1))
    igraph::plot.igraph(
        g,
        layout = layout,
        vertex.color = node_colors,
        vertex.size = pmin(8 + igraph::degree(g) * 0.5, 20),
        vertex.label = labels,
        vertex.label.cex = 0.6,
        vertex.label.color = "black",
        vertex.frame.color = NA,
        edge.arrow.size = 0.4,
        edge.color = edge_colors,
        edge.width = 0.8,
        main = "COSMOS Causal Network: Signaling to Metabolism"
    )
    legend("bottomleft",
           legend = names(type_colors),
           fill = type_colors,
           border = NA,
           bty = "n",
           cex = 0.9,
           title = "Node Type")
    legend("bottomright",
           legend = c("Activating", "Inhibiting"),
           col = c("#E64B35", "#4DBBD5"),
           lwd = 2,
           bty = "n",
           cex = 0.9,
           title = "Edge Sign")
    dev.off()

    plots$network <- out_png
    message("  COSMOS network plot saved: ", out_png)
    plots
}
