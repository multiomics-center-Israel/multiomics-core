# =============================================================================
# Protein-Protein Interaction Network Analysis (multiomics-core)
# =============================================================================
# Ported from multiomics/R/07b_ppi_networks.R with config adaptations:
#   log_message(...) → message(...)
#   config$ppi_analysis$... → cfg$ppi$... (cfg = config$modes$proteomics)
#   config$output$output_dir → out_dir parameter
#   DA data bridged via extract_de_table_for_pathway() + map_proteins_to_gene_symbols()

# -----------------------------------------------------------------------------
# Main Orchestrator
# -----------------------------------------------------------------------------

#' Run PPI Network Analysis
#'
#' @param da_results  Named list of per-contrast DE tables (FeatureID-based, gene-mapped)
#' @param metadata    Sample metadata data.frame
#' @param config      Full pipeline config list
#' @param out_dir     Output root directory
#' @return List with network, topology, communities, active_subnetwork, complexes,
#'         subcellular, enrichment, figures, summary
#' @export
run_ppi_network_analysis <- function(da_results, metadata, config, out_dir) {
    message("=== Running PPI Network Analysis ===")

    cfg <- config$modes$proteomics
    ppi_cfg <- cfg$ppi %||% list()

    # Resolve organism name for OrgDb dispatch (used by subcellular + enrichment)
    organism <- normalize_organism_name(cfg$annotation$organism %||% "Homo sapiens")

    # Set up output directories
    output_dir <- file.path(out_dir, "ppi_networks")
    plots_dir  <- file.path(output_dir, "plots")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (!dir.exists(plots_dir))  dir.create(plots_dir, recursive = TRUE)

    sig_threshold <- ppi_cfg$significance_threshold %||% 0.05
    lfc_threshold <- ppi_cfg$lfc_threshold %||% 0

    # Merge all contrast tables into one combined DA data frame
    # Use the first contrast's table as the primary
    da_df <- NULL
    if (is.data.frame(da_results)) {
        da_df <- da_results
    } else if (is.list(da_results) && length(da_results) > 0) {
        da_df <- da_results[[1]]
    }

    if (is.null(da_df) || !is.data.frame(da_df)) {
        message("WARNING: Could not extract DE data frame. Skipping PPI analysis.")
        return(NULL)
    }

    # Standardise column names
    if (!"padj" %in% colnames(da_df) && "adj.P.Val" %in% colnames(da_df))
        da_df$padj <- da_df$adj.P.Val
    if (!"log2FoldChange" %in% colnames(da_df) && "logFC" %in% colnames(da_df))
        da_df$log2FoldChange <- da_df$logFC
    if (!"protein_id" %in% colnames(da_df)) {
        if ("FeatureID" %in% colnames(da_df)) {
            da_df$protein_id <- da_df$FeatureID
        } else if (!is.null(rownames(da_df))) {
            da_df$protein_id <- rownames(da_df)
        }
    }
    if (!"pvalue" %in% colnames(da_df) && "P.Value" %in% colnames(da_df))
        da_df$pvalue <- da_df$P.Value

    # Filter significant proteins
    sig_proteins <- tryCatch({
        idx <- !is.na(da_df$padj) & da_df$padj < sig_threshold &
               abs(da_df$log2FoldChange) > lfc_threshold
        da_df$protein_id[idx]
    }, error = function(e) character(0))

    message("  Found ", length(sig_proteins), " significant proteins for network analysis")

    if (length(sig_proteins) < 5) {
        message("WARNING: Too few significant proteins for network analysis. Skipping.")
        return(NULL)
    }

    # Build PPI network from STRING
    ppi_network <- build_string_network(
        proteins        = sig_proteins,
        species         = ppi_cfg$species %||% 9606,
        score_threshold = ppi_cfg$string_score_threshold %||% 400,
        config          = config
    )

    if (is.null(ppi_network) || igraph::ecount(ppi_network$graph) == 0) {
        message("WARNING: No PPI edges found. Skipping network analysis.")
        return(NULL)
    }

    # Network topology analysis
    topology_results <- analyze_network_topology(
        graph = ppi_network$graph, da_results = da_df, output_dir = output_dir
    )

    # Community detection
    community_results <- detect_network_communities(
        graph = ppi_network$graph, da_results = da_df, output_dir = output_dir
    )

    # Active subnetwork identification
    active_subnet <- NULL
    if (isTRUE(ppi_cfg$active_subnetwork)) {
        active_subnet <- tryCatch(
            identify_active_subnetworks(ppi_network$graph, da_df, config, output_dir),
            error = function(e) { message("  Active subnetwork failed: ", e$message); NULL }
        )
    }

    # Protein complex analysis (CORUM)
    complex_results <- NULL
    if (isTRUE(ppi_cfg$complex_analysis)) {
        complex_results <- analyze_protein_complexes(
            sig_proteins, da_df, ppi_cfg$species %||% 9606, output_dir
        )
    }

    # Subcellular localization
    subcell_results <- analyze_subcellular_localization(sig_proteins, da_df, output_dir, organism = organism)

    # Network-based enrichment
    network_enrichment <- run_network_enrichment(
        ppi_network$graph, community_results, config, output_dir
    )

    # STRING functional enrichment (diseases, tissues, compartments, literature)
    string_enrichment <- get_string_functional_enrichment(
        sig_proteins, species_id, da_df, output_dir
    )

    # Visualizations
    figures <- generate_ppi_plots(
        ppi_network, topology_results, community_results,
        complex_results, subcell_results, da_df, config, plots_dir
    )

    results <- list(
        network               = ppi_network,
        topology              = topology_results,
        communities           = community_results,
        active_subnetwork     = active_subnet,
        complexes             = complex_results,
        subcellular           = subcell_results,
        enrichment            = network_enrichment,
        string_enrichment     = string_enrichment,
        figures               = figures,
        summary               = create_ppi_summary(ppi_network, topology_results,
                                                   community_results, complex_results)
    )

    save_ppi_outputs(results, output_dir)
    message("PPI network analysis complete!")
    results
}

# -----------------------------------------------------------------------------
# STRING Network Construction
# -----------------------------------------------------------------------------

build_string_network <- function(proteins, species = 9606,
                                  score_threshold = 400, config = NULL) {
    message("Building STRING PPI network...")
    message("  Species: ", species, " (", get_species_name(species), ")")
    message("  Score threshold: ", score_threshold)

    if (!requireNamespace("STRINGdb", quietly = TRUE)) {
        message("WARNING: STRINGdb not available. Attempting API fallback...")
        return(build_string_network_api(proteins, species, score_threshold))
    }

    tryCatch({
        string_db <- STRINGdb::STRINGdb$new(
            version = "12.0", species = species,
            score_threshold = score_threshold, input_directory = ""
        )

        proteins_df <- data.frame(protein = proteins, stringsAsFactors = FALSE)
        mapped <- string_db$map(proteins_df, "protein", removeUnmappedRows = TRUE)

        message("  Mapped ", nrow(mapped), " of ", length(proteins), " proteins to STRING")

        if (nrow(mapped) < 2) return(NULL)

        interactions <- string_db$get_interactions(mapped$STRING_id)

        message("  Found ", nrow(interactions), " interactions")
        if (nrow(interactions) == 0) return(NULL)

        graph <- igraph::graph_from_data_frame(
            interactions[, c("from", "to", "combined_score")], directed = FALSE
        )

        string_to_protein <- setNames(mapped$protein, mapped$STRING_id)
        igraph::V(graph)$protein_name <- string_to_protein[igraph::V(graph)$name]

        graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)

        list(graph = graph, edges = interactions, mapping = mapped,
             n_nodes = igraph::vcount(graph), n_edges = igraph::ecount(graph))

    }, error = function(e) {
        message("ERROR in STRING network: ", e$message, ". Attempting API fallback...")
        build_string_network_api(proteins, species, score_threshold)
    })
}

build_string_network_api <- function(proteins, species, score_threshold,
                                      chunk_size = 500) {
    message("  Using STRING API fallback...")

    # Clean protein IDs: remove any containing "|" (empty gene symbols) and whitespace
    proteins <- proteins[!grepl("\\|", proteins) & nzchar(trimws(proteins))]
    if (length(proteins) == 0) {
        message("WARNING: No valid protein IDs for STRING API query")
        return(NULL)
    }

    # Do NOT chunk: the /network endpoint returns edges only among submitted identifiers,
    # so independent per-chunk queries silently drop all inter-chunk edges.
    # Cap explicitly instead so truncation is visible and the user can install STRINGdb.
    if (length(proteins) > chunk_size) {
        message("WARNING: ", length(proteins), " proteins exceed the STRING API fallback ",
                "limit (", chunk_size, "). Capping to ", chunk_size, " proteins. ",
                "Install STRINGdb for complete cross-protein results.")
        proteins <- proteins[seq_len(chunk_size)]
    }

    base_url <- "https://string-db.org/api/tsv/network"
    url <- paste0(base_url,
                  "?identifiers=", paste(proteins, collapse = "%0a"),
                  "&species=", species,
                  "&required_score=", score_threshold,
                  "&caller_identity=multiomics_pipeline")
    combined <- tryCatch(
        utils::read.delim(url(url), stringsAsFactors = FALSE),
        error = function(e) {
            message("WARNING: STRING API request failed: ", e$message)
            NULL
        }
    )
    if (is.null(combined) || nrow(combined) == 0) return(NULL)

    required_cols <- c("preferredName_A", "preferredName_B", "score")
    missing_cols  <- setdiff(required_cols, colnames(combined))
    if (length(missing_cols) > 0) {
        message("WARNING: STRING API response is missing expected network columns (",
                paste(missing_cols, collapse = ", "), "). ",
                "Got: ", paste(colnames(combined), collapse = ", "), ". ",
                "Response may be a non-network payload (e.g. error text). Skipping.")
        return(NULL)
    }

    graph <- igraph::graph_from_data_frame(
        combined[, required_cols], directed = FALSE
    )
    graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)

    message("  API returned ", igraph::ecount(graph), " unique edges")

    list(graph = graph, edges = combined, mapping = NULL,
         n_nodes = igraph::vcount(graph), n_edges = igraph::ecount(graph))
}

get_species_name <- function(taxid) {
    species_map <- c(
        "9606" = "Homo sapiens", "10090" = "Mus musculus",
        "10116" = "Rattus norvegicus", "7955" = "Danio rerio",
        "6239" = "Caenorhabditis elegans", "7227" = "Drosophila melanogaster"
    )
    species_map[as.character(taxid)] %||% paste0("Species ", taxid)
}

#' Map organism name to NCBI taxonomy ID
#' @param organism Character, organism name (e.g. "Homo sapiens")
#' @return Integer taxonomy ID (e.g. 9606), or NULL if unknown
organism_to_taxid <- function(organism) {
    if (is.null(organism) || is.na(organism)) return(NULL)
    taxid_map <- c(
        "Homo sapiens"              = 9606L,
        "Mus musculus"              = 10090L,
        "Rattus norvegicus"         = 10116L,
        "Danio rerio"              = 7955L,
        "Caenorhabditis elegans"   = 6239L,
        "Drosophila melanogaster"  = 7227L,
        "Saccharomyces cerevisiae" = 4932L,
        "Sus scrofa"               = 9823L,
        "Bos taurus"               = 9913L,
        "Gallus gallus"            = 9031L,
        "Ovis aries"               = 9940L,
        "Canis familiaris"         = 9615L,
        "Felis catus"              = 9685L,
        "Pan troglodytes"          = 9598L,
        "Macaca mulatta"           = 9544L,
        "Equus caballus"           = 9796L,
        "Oryctolagus cuniculus"    = 9986L,
        "Arabidopsis thaliana"     = 3702L,
        "Giardia lamblia"          = 5741L
    )
    organism_norm <- normalize_organism_name(organism)
    tid <- taxid_map[organism_norm]
    if (!is.na(tid)) return(as.integer(tid))
    NULL
}

# -----------------------------------------------------------------------------
# Network Topology
# -----------------------------------------------------------------------------

analyze_network_topology <- function(graph, da_results, output_dir) {
    message("Analyzing network topology...")

    n_nodes <- igraph::vcount(graph)
    n_edges <- igraph::ecount(graph)
    density <- igraph::edge_density(graph)
    diameter <- igraph::diameter(graph)
    avg_path_length <- igraph::mean_distance(graph)
    clustering_coef <- igraph::transitivity(graph, type = "global")

    degree <- igraph::degree(graph)
    betweenness <- igraph::betweenness(graph)
    closeness <- igraph::closeness(graph)
    eigenvector <- tryCatch(igraph::eigen_centrality(graph)$vector,
                            error = function(e) rep(NA, n_nodes))

    node_names <- igraph::V(graph)$name
    protein_names <- igraph::V(graph)$protein_name
    if (is.null(protein_names)) protein_names <- node_names

    node_metrics <- data.frame(
        node_id = node_names, protein_name = protein_names,
        degree = degree, betweenness = betweenness,
        closeness = closeness, eigenvector = eigenvector,
        stringsAsFactors = FALSE
    )

    if (!is.null(da_results) && "protein_id" %in% colnames(da_results)) {
        node_metrics <- merge(
            node_metrics,
            da_results[, intersect(c("protein_id", "log2FoldChange", "padj"), colnames(da_results))],
            by.x = "protein_name", by.y = "protein_id", all.x = TRUE
        )
    }

    degree_threshold <- quantile(degree, 0.9)
    hub_proteins <- node_metrics[node_metrics$degree >= degree_threshold, ]
    hub_proteins <- hub_proteins[order(-hub_proteins$degree), ]

    betweenness_threshold <- quantile(betweenness, 0.9)
    bottleneck_proteins <- node_metrics[node_metrics$betweenness >= betweenness_threshold, ]
    bottleneck_proteins <- bottleneck_proteins[order(-bottleneck_proteins$betweenness), ]

    write.csv(node_metrics, file.path(output_dir, "network_node_metrics.csv"), row.names = FALSE)
    write.csv(hub_proteins, file.path(output_dir, "hub_proteins.csv"), row.names = FALSE)

    list(
        global_metrics = list(
            n_nodes = n_nodes, n_edges = n_edges, density = density,
            diameter = diameter, avg_path_length = avg_path_length,
            clustering_coefficient = clustering_coef
        ),
        node_metrics = node_metrics,
        hub_proteins = hub_proteins,
        bottleneck_proteins = bottleneck_proteins
    )
}

# -----------------------------------------------------------------------------
# Community Detection
# -----------------------------------------------------------------------------

detect_network_communities <- function(graph, da_results, output_dir) {
    message("Detecting network communities...")

    communities <- list()

    tryCatch({
        communities$louvain <- igraph::cluster_louvain(graph)
    }, error = function(e) message("  Louvain failed: ", e$message))

    tryCatch({
        communities$walktrap <- igraph::cluster_walktrap(graph)
    }, error = function(e) message("  Walktrap failed: ", e$message))

    tryCatch({
        communities$fast_greedy <- igraph::cluster_fast_greedy(graph)
    }, error = function(e) message("  Fast greedy failed: ", e$message))

    primary_comm <- communities$louvain %||% communities$walktrap
    if (is.null(primary_comm)) {
        message("WARNING: No community detection succeeded")
        return(NULL)
    }

    membership <- igraph::membership(primary_comm)
    modularity_val <- igraph::modularity(primary_comm)

    node_names <- igraph::V(graph)$name
    protein_names <- igraph::V(graph)$protein_name
    if (is.null(protein_names)) protein_names <- node_names

    community_df <- data.frame(
        node_id = node_names, protein_name = protein_names,
        community = as.integer(membership), stringsAsFactors = FALSE
    )

    if (!is.null(da_results) && "protein_id" %in% colnames(da_results)) {
        community_df <- merge(
            community_df,
            da_results[, intersect(c("protein_id", "log2FoldChange", "padj"), colnames(da_results))],
            by.x = "protein_name", by.y = "protein_id", all.x = TRUE
        )
    }

    write.csv(community_df, file.path(output_dir, "community_membership.csv"), row.names = FALSE)

    list(
        primary = primary_comm, all_methods = communities,
        membership = community_df, modularity = modularity_val,
        n_communities = max(membership)
    )
}

# -----------------------------------------------------------------------------
# Active Subnetwork Identification
# -----------------------------------------------------------------------------

identify_active_subnetworks <- function(graph, da_results, config, output_dir) {
    message("Identifying active subnetworks...")

    node_names <- igraph::V(graph)$protein_name
    if (is.null(node_names)) node_names <- igraph::V(graph)$name

    scores <- sapply(node_names, function(p) {
        idx <- which(da_results$protein_id == p)
        if (length(idx) == 0) return(0)
        pval <- da_results$pvalue[idx[1]]
        lfc  <- da_results$log2FoldChange[idx[1]]
        if (is.na(pval) || is.na(lfc)) return(0)
        -log10(pval + 1e-300) * sign(lfc)
    })

    igraph::V(graph)$score <- scores

    abs_scores <- abs(scores)
    abs_scores[is.na(abs_scores)] <- 0
    score_threshold <- quantile(abs_scores, 0.75, na.rm = TRUE)
    if (is.na(score_threshold) || score_threshold == 0) return(NULL)

    active_nodes <- which(abs_scores >= score_threshold)
    if (length(active_nodes) < 3) return(NULL)

    active_subgraph <- igraph::induced_subgraph(graph, active_nodes)
    components <- igraph::components(active_subgraph)
    comp_sizes <- table(components$membership)
    large_comps <- as.integer(names(comp_sizes)[comp_sizes >= 3])

    active_subnets <- lapply(large_comps, function(comp_id) {
        comp_nodes <- which(components$membership == comp_id)
        subnet <- igraph::induced_subgraph(active_subgraph, comp_nodes)
        node_names_sub <- igraph::V(subnet)$protein_name %||% igraph::V(subnet)$name
        scores_sub <- igraph::V(subnet)$score
        list(graph = subnet, proteins = node_names_sub, scores = scores_sub,
             mean_score = mean(scores_sub, na.rm = TRUE),
             n_nodes = igraph::vcount(subnet), n_edges = igraph::ecount(subnet))
    })

    active_subnets <- active_subnets[order(
        sapply(active_subnets, function(x) abs(x$mean_score)), decreasing = TRUE
    )]

    message("  Found ", length(active_subnets), " active subnetworks")

    subnet_summary <- data.frame(
        subnet_id  = seq_along(active_subnets),
        n_nodes    = sapply(active_subnets, `[[`, "n_nodes"),
        n_edges    = sapply(active_subnets, `[[`, "n_edges"),
        mean_score = sapply(active_subnets, `[[`, "mean_score"),
        proteins   = sapply(active_subnets, function(x) paste(head(x$proteins, 10), collapse = ", ")),
        stringsAsFactors = FALSE
    )
    write.csv(subnet_summary, file.path(output_dir, "active_subnetworks.csv"), row.names = FALSE)

    list(subnetworks = active_subnets, summary = subnet_summary)
}

# -----------------------------------------------------------------------------
# Protein Complex Analysis (CORUM)
# -----------------------------------------------------------------------------

analyze_protein_complexes <- function(proteins, da_results, species, output_dir) {
    message("Analyzing protein complexes (CORUM)...")

    corum_db <- get_corum_database(species)
    if (is.null(corum_db)) {
        message("  CORUM database not available. Skipping.")
        return(NULL)
    }

    proteins_upper <- toupper(proteins)

    complex_hits <- lapply(seq_len(nrow(corum_db)), function(i) {
        complex_proteins <- toupper(unlist(strsplit(corum_db$subunits_gene_name[i], ";")))
        overlap <- intersect(proteins_upper, complex_proteins)
        if (length(overlap) >= 2) {
            data.frame(
                complex_id = corum_db$complex_id[i],
                complex_name = corum_db$complex_name[i],
                n_total_subunits = length(complex_proteins),
                n_detected = length(overlap),
                detection_ratio = length(overlap) / length(complex_proteins),
                detected_proteins = paste(overlap, collapse = ";"),
                stringsAsFactors = FALSE
            )
        } else NULL
    })

    complex_df <- do.call(rbind, complex_hits[!sapply(complex_hits, is.null)])
    if (is.null(complex_df) || nrow(complex_df) == 0) {
        message("  No protein complexes found with >= 2 detected subunits")
        return(NULL)
    }

    complex_df <- complex_df[order(complex_df$detection_ratio, decreasing = TRUE), ]

    # Hypergeometric enrichment
    total_in_corum <- length(unique(unlist(strsplit(corum_db$subunits_gene_name, ";"))))
    n_query <- length(proteins)
    complex_df$pvalue <- sapply(seq_len(nrow(complex_df)), function(i) {
        phyper(complex_df$n_detected[i] - 1, complex_df$n_total_subunits[i],
               total_in_corum - complex_df$n_total_subunits[i], n_query, lower.tail = FALSE)
    })
    complex_df$padj <- p.adjust(complex_df$pvalue, method = "BH")

    write.csv(complex_df, file.path(output_dir, "protein_complex_analysis.csv"), row.names = FALSE)

    list(complexes = complex_df, n_total = nrow(complex_df),
         n_significant = sum(complex_df$padj < 0.05, na.rm = TRUE))
}

get_corum_database <- function(species) {
    corum_file <- system.file("extdata", "corum_core.csv", package = "multiomics")
    if (file.exists(corum_file)) return(read.csv(corum_file, stringsAsFactors = FALSE))

    tryCatch({
        message("  Downloading CORUM database...")
        url <- "https://mips.helmholtz-muenchen.de/corum/download/coreComplexes.txt"
        corum <- read.delim(url(url), stringsAsFactors = FALSE)
        if (species == 9606) corum <- corum[corum$Organism == "Human", ]
        else if (species == 10090) corum <- corum[corum$Organism == "Mouse", ]
        names(corum) <- gsub(" ", "_", tolower(names(corum)))
        corum$complex_id <- corum$complexid %||% corum$complex_id
        corum$complex_name <- corum$complexname %||% corum$complex_name
        corum$subunits_gene_name <- corum$`subunits(gene_name)` %||% corum$subunits_gene_name
        corum
    }, error = function(e) { message("WARNING: Could not load CORUM: ", e$message); NULL })
}

# -----------------------------------------------------------------------------
# Subcellular Localization
# -----------------------------------------------------------------------------

analyze_subcellular_localization <- function(proteins, da_results, output_dir,
                                             organism = "Homo sapiens") {
    message("Analyzing subcellular localization...")

    orgdb_pkg <- get_orgdb_package(organism)
    if (is.null(orgdb_pkg) || !requireNamespace(orgdb_pkg, quietly = TRUE)) {
        message("  OrgDb package for '", organism, "' not available. Skipping subcellular localization.")
        return(NULL)
    }
    org_db <- getExportedValue(orgdb_pkg, orgdb_pkg)

    tryCatch({
        go_cc <- AnnotationDbi::select(org_db,
            keys = proteins, keytype = "SYMBOL", columns = c("GOALL", "ONTOLOGYALL"))
        go_cc <- go_cc[go_cc$ONTOLOGYALL == "CC" & !is.na(go_cc$GOALL), ]
        if (nrow(go_cc) == 0) return(NULL)

        go_terms <- AnnotationDbi::select(GO.db::GO.db,
            keys = unique(go_cc$GOALL), keytype = "GOID", columns = "TERM")
        go_cc <- merge(go_cc, go_terms, by.x = "GOALL", by.y = "GOID")

        location_counts <- aggregate(SYMBOL ~ TERM, data = go_cc,
                                     FUN = function(x) length(unique(x)))
        colnames(location_counts) <- c("TERM", "n_proteins")
        location_counts <- location_counts[order(-location_counts$n_proteins), ]

        major_locs <- c("nucleus", "cytoplasm", "membrane", "mitochondrion",
                        "endoplasmic reticulum", "Golgi", "lysosome", "ribosome",
                        "extracellular", "cytoskeleton", "plasma membrane")
        location_summary <- location_counts[
            grepl(paste(major_locs, collapse = "|"), location_counts$TERM, ignore.case = TRUE), ]
        location_summary <- head(location_summary, 20)

        write.csv(location_summary, file.path(output_dir, "subcellular_localization.csv"),
                  row.names = FALSE)

        list(locations = location_summary, all_annotations = go_cc)

    }, error = function(e) { message("  Subcellular analysis failed: ", e$message); NULL })
}

# -----------------------------------------------------------------------------
# Network Enrichment
# -----------------------------------------------------------------------------

run_network_enrichment <- function(graph, community_results, config, output_dir) {
    message("Running network-based enrichment...")

    if (is.null(community_results)) return(NULL)
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        message("  clusterProfiler not available. Skipping.")
        return(NULL)
    }

    cfg_prot <- config$modes$proteomics
    organism <- normalize_organism_name(cfg_prot$annotation$organism %||% "Homo sapiens")
    orgdb_pkg <- get_orgdb_package(organism)
    if (is.null(orgdb_pkg) || !requireNamespace(orgdb_pkg, quietly = TRUE)) {
        message("  OrgDb package for '", organism, "' not available. Skipping network enrichment.")
        return(NULL)
    }
    org_db <- getExportedValue(orgdb_pkg, orgdb_pkg)

    enrichment_results <- list()
    community_df <- community_results$membership
    communities_ids <- unique(community_df$community)

    # Batch all symbol → Entrez lookups in a single OrgDb call before the loop
    all_comm_symbols <- unique(community_df$protein_name)
    entrez_batch <- tryCatch(
        AnnotationDbi::mapIds(org_db, keys = all_comm_symbols,
                              keytype = "SYMBOL", column = "ENTREZID",
                              multiVals = "first"),
        error = function(e) {
            message("  mapIds batch failed: ", e$message)
            setNames(rep(NA_character_, length(all_comm_symbols)), all_comm_symbols)
        }
    )

    for (comm in communities_ids) {
        comm_proteins <- community_df$protein_name[community_df$community == comm]
        if (length(comm_proteins) < 3) next

        tryCatch({
            entrez <- entrez_batch[comm_proteins]
            entrez <- entrez[!is.na(entrez)]
            if (length(entrez) < 3) next

            go_result <- clusterProfiler::enrichGO(
                gene = entrez, OrgDb = org_db,
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1
            )

            if (!is.null(go_result) && nrow(go_result@result) > 0) {
                enrichment_results[[paste0("community_", comm)]] <- list(
                    go = go_result@result, n_proteins = length(comm_proteins), proteins = comm_proteins
                )
            }
        }, error = function(e) NULL)
    }

    if (length(enrichment_results) == 0) return(NULL)

    all_enrichments <- lapply(names(enrichment_results), function(comm) {
        df <- enrichment_results[[comm]]$go
        df$community <- comm
        df
    })
    combined <- do.call(rbind, all_enrichments)
    write.csv(combined, file.path(output_dir, "community_enrichment.csv"), row.names = FALSE)

    enrichment_results
}

# -----------------------------------------------------------------------------
# STRING Functional Enrichment (Diseases, Tissues, Compartments, Literature)
# -----------------------------------------------------------------------------

#' Query STRING API for functional enrichment across multiple categories
#'
#' Retrieves enrichment results from STRING for: Process (GO), Component,
#' Function, KEGG, Reactome, DISEASES, COMPARTMENTS, TISSUES, and PMID.
#'
#' @param proteins Character vector of protein/gene symbols
#' @param species Integer NCBI taxonomy ID
#' @param da_results DA data frame (used to rank proteins by significance)
#' @param output_dir Output directory for CSV files
#' @return List with enrichment data frames per category
get_string_functional_enrichment <- function(proteins, species, da_results, output_dir) {
    message("=== STRING Functional Enrichment ===")

    species <- resolve_string_taxid(species)

    # Rank proteins by significance (stat = sign(lfc) * -log10(p)) for prioritized input
    if (!is.null(da_results) && "pvalue" %in% colnames(da_results) &&
        "log2FoldChange" %in% colnames(da_results)) {
        da_sub <- da_results[da_results$protein_id %in% proteins, ]
        da_sub$rank_stat <- sign(da_sub$log2FoldChange) * -log10(da_sub$pvalue + 1e-300)
        da_sub <- da_sub[order(abs(da_sub$rank_stat), decreasing = TRUE), ]
        proteins <- da_sub$protein_id
        message("  Proteins ranked by fGSEA-style stat for STRING query")
    }

    # Clean protein IDs
    proteins <- proteins[!grepl("\\|", proteins) & nzchar(trimws(proteins))]
    if (length(proteins) == 0) {
        message("  No valid protein IDs. Skipping STRING enrichment.")
        return(NULL)
    }
    if (length(proteins) > 500) {
        message("  Limiting to top 500 ranked proteins for STRING enrichment")
        proteins <- proteins[1:500]
    }

    # Query STRING enrichment API
    protein_list <- paste(proteins, collapse = "%0D")
    url <- paste0("https://string-db.org/api/tsv/enrichment",
                  "?identifiers=", protein_list,
                  "&species=", species,
                  "&caller_identity=multiomics_pipeline")

    enrichment_df <- tryCatch({
        resp <- utils::read.delim(url(url), stringsAsFactors = FALSE)
        if (nrow(resp) == 0) {
            message("  STRING enrichment returned no results.")
            return(NULL)
        }
        resp
    }, error = function(e) {
        message("  STRING enrichment API failed: ", e$message)

        # Fallback: try STRINGdb package
        if (requireNamespace("STRINGdb", quietly = TRUE)) {
            tryCatch({
                string_db <- STRINGdb::STRINGdb$new(
                    version = "12.0", species = species,
                    score_threshold = 0, input_directory = ""
                )
                proteins_df <- data.frame(protein = proteins, stringsAsFactors = FALSE)
                mapped <- string_db$map(proteins_df, "protein", removeUnmappedRows = TRUE)
                if (nrow(mapped) >= 2) {
                    enrich <- string_db$get_enrichment(mapped$STRING_id)
                    return(enrich)
                }
                NULL
            }, error = function(e2) {
                message("  STRINGdb fallback also failed: ", e2$message)
                NULL
            })
        } else NULL
    })

    if (is.null(enrichment_df) || nrow(enrichment_df) == 0) return(NULL)

    # Split by category
    categories <- c("Process", "Component", "Function", "KEGG", "Reactome",
                     "DISEASES", "COMPARTMENTS", "TISSUES", "PMID",
                     "Keyword", "InterPro", "Pfam", "NetworkNeighborAL")

    string_dir <- file.path(output_dir, "string_enrichment")
    dir.create(string_dir, recursive = TRUE, showWarnings = FALSE)

    results <- list()
    cat_col <- if ("category" %in% colnames(enrichment_df)) "category" else NULL

    if (!is.null(cat_col)) {
        for (cat in unique(enrichment_df[[cat_col]])) {
            cat_df <- enrichment_df[enrichment_df[[cat_col]] == cat, ]
            cat_df <- cat_df[order(cat_df$fdr), ]
            clean_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
            results[[cat]] <- cat_df

            out_file <- file.path(string_dir, paste0("string_", clean_cat, ".csv"))
            write.csv(cat_df, out_file, row.names = FALSE)
            n_sig <- sum(cat_df$fdr < 0.05, na.rm = TRUE)
            message("  ", cat, ": ", nrow(cat_df), " terms (", n_sig, " FDR < 0.05)")
        }
    } else {
        results$all <- enrichment_df
        write.csv(enrichment_df, file.path(string_dir, "string_all.csv"), row.names = FALSE)
    }

    # Summary
    write.csv(enrichment_df, file.path(string_dir, "string_enrichment_full.csv"), row.names = FALSE)
    message("  STRING enrichment saved to: ", string_dir)

    results
}

# -----------------------------------------------------------------------------
# Visualizations
# -----------------------------------------------------------------------------

generate_ppi_plots <- function(ppi_network, topology_results, community_results,
                                complex_results, subcell_results, da_results,
                                config, plots_dir) {
    message("Generating PPI network plots...")
    figures <- list()

    figures$network <- tryCatch(
        plot_ppi_network(ppi_network$graph, community_results, da_results, plots_dir),
        error = function(e) { message("  Network plot failed: ", e$message); NULL }
    )

    figures$topology <- tryCatch(
        plot_network_topology(topology_results, plots_dir),
        error = function(e) { message("  Topology plot failed: ", e$message); NULL }
    )

    if (!is.null(community_results)) {
        figures$communities <- tryCatch(
            plot_community_network(ppi_network$graph, community_results, plots_dir),
            error = function(e) { message("  Community plot failed: ", e$message); NULL }
        )
    }

    figures$hubs <- tryCatch(
        plot_hub_proteins(topology_results, da_results, plots_dir),
        error = function(e) { message("  Hub plot failed: ", e$message); NULL }
    )

    if (!is.null(subcell_results)) {
        figures$subcellular <- tryCatch(
            plot_subcellular_enrichment(subcell_results, plots_dir),
            error = function(e) { message("  Subcellular plot failed: ", e$message); NULL }
        )
    }

    figures
}

plot_ppi_network <- function(graph, community_results, da_results, plots_dir) {
    node_names <- igraph::V(graph)$protein_name %||% igraph::V(graph)$name

    log2fc <- sapply(node_names, function(p) {
        idx <- which(da_results$protein_id == p)
        if (length(idx) > 0) da_results$log2FoldChange[idx[1]] else 0
    })

    # Simple color scale (blue → white → red)
    col_vals <- pmax(pmin(log2fc, 2), -2)
    col_scaled <- (col_vals + 2) / 4  # 0 to 1
    node_colors <- grDevices::rgb(col_scaled, 1 - abs(col_scaled - 0.5) * 2, 1 - col_scaled)

    degrees <- igraph::degree(graph)
    node_sizes <- 3 + 5 * (degrees / max(degrees, 1))
    layout <- igraph::layout_with_fr(graph)

    fig_path <- file.path(plots_dir, "ppi_network.png")
    grDevices::png(fig_path, width = 12, height = 10, units = "in", res = 150)
    par(mar = c(0, 0, 2, 0))
    plot(graph, layout = layout, vertex.color = node_colors, vertex.size = node_sizes,
         vertex.label = if (igraph::vcount(graph) <= 50) node_names else NA,
         vertex.label.cex = 0.6, vertex.label.color = "black", vertex.frame.color = NA,
         edge.width = 0.5, edge.color = grDevices::adjustcolor("grey50", alpha.f = 0.5),
         main = "Protein-Protein Interaction Network")
    legend("bottomleft", legend = c("Down-regulated", "No change", "Up-regulated"),
           fill = c("blue", "white", "red"), border = "black", title = "log2FC", cex = 0.8)
    grDevices::dev.off()
    fig_path
}

plot_network_topology <- function(topology_results, plots_dir) {
    node_metrics <- topology_results$node_metrics

    p1 <- ggplot2::ggplot(node_metrics, ggplot2::aes(x = degree)) +
        ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "black") +
        ggplot2::labs(title = "Degree Distribution", x = "Degree", y = "Count") +
        ggplot2::theme_minimal()

    p2 <- ggplot2::ggplot(node_metrics, ggplot2::aes(x = degree, y = betweenness)) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::labs(title = "Degree vs Betweenness", x = "Degree", y = "Betweenness") +
        ggplot2::theme_minimal()

    combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
        patchwork::wrap_plots(p1, p2, ncol = 2)
    } else {
        p1
    }

    fig_path <- file.path(plots_dir, "network_topology.png")
    ggplot2::ggsave(fig_path, combined, width = 12, height = 5)
    fig_path
}

plot_community_network <- function(graph, community_results, plots_dir) {
    membership <- community_results$membership
    node_names <- igraph::V(graph)$protein_name %||% igraph::V(graph)$name
    community_map <- setNames(membership$community, membership$protein_name)
    node_communities <- community_map[node_names]
    n_comm <- max(node_communities, na.rm = TRUE)

    palette <- if (requireNamespace("RColorBrewer", quietly = TRUE)) {
        RColorBrewer::brewer.pal(min(12, max(3, n_comm)), "Set3")
    } else {
        grDevices::rainbow(max(3, n_comm))
    }
    node_colors <- palette[(node_communities %% length(palette)) + 1]

    fig_path <- file.path(plots_dir, "network_communities.png")
    grDevices::png(fig_path, width = 12, height = 10, units = "in", res = 150)
    par(mar = c(0, 0, 2, 0))
    plot(graph, layout = igraph::layout_with_fr(graph), vertex.color = node_colors,
         vertex.size = 5, vertex.label = NA, vertex.frame.color = NA,
         edge.width = 0.3, edge.color = grDevices::adjustcolor("grey50", alpha.f = 0.3),
         main = paste0("Network Communities (n=", n_comm,
                       ", modularity=", round(community_results$modularity, 3), ")"))
    grDevices::dev.off()
    fig_path
}

plot_hub_proteins <- function(topology_results, da_results, plots_dir) {
    hub_df <- topology_results$hub_proteins
    if (nrow(hub_df) == 0) return(NULL)

    hub_df <- head(hub_df, 20)
    hub_df$color <- ifelse(is.na(hub_df$log2FoldChange), "grey",
                           ifelse(hub_df$log2FoldChange > 0, "red", "blue"))

    p <- ggplot2::ggplot(hub_df,
             ggplot2::aes(x = reorder(protein_name, degree), y = degree, fill = color)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_identity() +
        ggplot2::coord_flip() +
        ggplot2::labs(title = "Top Hub Proteins by Degree", x = "Protein",
                      y = "Degree (Number of Interactions)") +
        ggplot2::theme_minimal()

    fig_path <- file.path(plots_dir, "hub_proteins.png")
    ggplot2::ggsave(fig_path, p, width = 8, height = 6)
    fig_path
}

plot_subcellular_enrichment <- function(subcell_results, plots_dir) {
    locations <- subcell_results$locations
    if (is.null(locations) || nrow(locations) == 0) return(NULL)
    locations <- head(locations, 15)

    p <- ggplot2::ggplot(locations,
             ggplot2::aes(x = reorder(TERM, n_proteins), y = n_proteins)) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        ggplot2::coord_flip() +
        ggplot2::labs(title = "Subcellular Localization of Significant Proteins",
                      x = "Cellular Compartment", y = "Number of Proteins") +
        ggplot2::theme_minimal()

    fig_path <- file.path(plots_dir, "subcellular_enrichment.png")
    ggplot2::ggsave(fig_path, p, width = 10, height = 8)
    fig_path
}

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

create_ppi_summary <- function(ppi_network, topology_results, community_results,
                                complex_results) {
    summary <- list(
        n_nodes = ppi_network$n_nodes, n_edges = ppi_network$n_edges,
        density = topology_results$global_metrics$density,
        n_hub_proteins = nrow(topology_results$hub_proteins),
        top_hubs = head(topology_results$hub_proteins$protein_name, 5)
    )
    if (!is.null(community_results)) {
        summary$n_communities <- community_results$n_communities
        summary$modularity <- community_results$modularity
    }
    if (!is.null(complex_results)) {
        summary$n_complexes_detected <- complex_results$n_total
        summary$n_significant_complexes <- complex_results$n_significant
    }
    summary
}

save_ppi_outputs <- function(results, output_dir) {
    message("Saving PPI outputs...")
    if (!is.null(results$network$graph)) {
        tryCatch(
            igraph::write_graph(results$network$graph,
                                file.path(output_dir, "ppi_network.graphml"), format = "graphml"),
            error = function(e) message("  Could not save graphml: ", e$message)
        )
    }
    if (!is.null(results$network$edges)) {
        write.csv(results$network$edges, file.path(output_dir, "ppi_network_edges.csv"),
                  row.names = FALSE)
    }
    saveRDS(results, file.path(output_dir, "ppi_analysis_results.rds"))
    message("  PPI outputs saved to: ", output_dir)
}
