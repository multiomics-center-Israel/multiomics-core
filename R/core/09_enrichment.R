#' Core Enrichment Utilities
#'
#' Shared functions for gene set / pathway analysis across omics types.
#' Provides GMT reading, gene set loading, and ORA (Fisher's exact test).

# ==============================================================================
# GMT READING
# ==============================================================================

#' Read GMT file for custom gene sets
#'
#' @param gmt_file Path to GMT file
#' @return Named list of gene sets (character vectors)
#' @export
read_gmt <- function(gmt_file) {
    if (!file.exists(gmt_file)) {
        stop("GMT file not found: ", gmt_file)
    }

    lines <- readLines(gmt_file)
    gene_sets <- list()
    descriptions <- character()

    for (line in lines) {
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) >= 3) {
            set_name <- parts[1]
            genes <- parts[3:length(parts)]
            genes <- genes[genes != ""]
            gene_sets[[set_name]] <- genes
            descriptions[set_name] <- parts[2]
        }
    }

    message("Loaded ", length(gene_sets), " gene sets from GMT file")
    attr(gene_sets, "descriptions") <- descriptions
    gene_sets
}

# ==============================================================================
# GENE SET LOADING
# ==============================================================================

#' Load gene sets from various sources
#'
#' @param organism Organism name (used for OrgDb/KEGG lookups)
#' @param pathway_database Character vector of databases to use (e.g. "GO", "KEGG")
#' @param gmt_file Optional custom GMT file path
#' @param annotation Gene annotation data frame (with gene_id and entrez_id columns)
#' @return Named list of gene set collections (each a named list of character vectors)
#' @export
load_gene_sets <- function(organism,
                           pathway_database = c("GO", "KEGG"),
                           gmt_file = NULL,
                           annotation = NULL,
                           target_id_type = "ensembl") {

    gene_sets <- list()
    org_info <- get_organism_info(organism)

    # Custom GMT takes priority
    if (!is.null(gmt_file) && file.exists(gmt_file)) {
        gene_sets$custom <- read_gmt(gmt_file)
        message("Loaded custom gene sets from: ", gmt_file)

        # Validate GMT coverage against annotation features if available
        if (!is.null(annotation) && "gene_id" %in% colnames(annotation)) {
            feature_ids <- annotation$gene_id
        } else if (!is.null(annotation)) {
            feature_ids <- annotation[[1]]
        } else {
            feature_ids <- NULL
        }

        if (!is.null(feature_ids) && length(feature_ids) > 0) {
            gmt_val <- tryCatch(
                validate_gmt(gene_sets$custom, feature_ids, verbose = TRUE),
                error = function(e) {
                    warning("GMT validation failed: ", e$message)
                    NULL
                }
            )
            if (!is.null(gmt_val) && length(gmt_val$filtered_pathways) > 0) {
                gene_sets$custom <- gmt_val$filtered_pathways
                message("GMT filtered to ", length(gmt_val$filtered_pathways),
                        " pathways with coverage in data")
            }
        }
    }

    # For non-model organisms with a KEGG code, auto-fetch KEGG pathways
    if (!org_info$supported && length(gene_sets) == 0 && !is.na(org_info$kegg)) {
        message("Non-model organism with KEGG code '", org_info$kegg,
                "' — attempting to fetch KEGG pathways via KEGGREST...")
        kegg_sets <- tryCatch(
            fetch_kegg_via_rest(org_info$kegg),
            error = function(e) {
                warning("Failed to fetch KEGG pathways via REST: ", e$message)
                list()
            }
        )
        if (length(kegg_sets) > 0) {
            gene_sets$KEGG <- kegg_sets
            message("Loaded ", length(kegg_sets), " KEGG pathways via KEGGREST")
        }
    }

    if (!org_info$supported && length(gene_sets) == 0) {
        # Fallback: try to generate GMT from GO or KEGG via biomaRt/KEGGREST
        message("Attempting to generate gene sets for non-model organism: ", organism)

        go_sets <- tryCatch({
            generate_gmt_from_go(organism, id_type = "symbol", verbose = TRUE)
        }, error = function(e) {
            message("  GO GMT generation failed: ", e$message)
            list()
        })

        if (length(go_sets) > 0) {
            gene_sets$GO <- go_sets
            message("Generated ", length(go_sets), " GO gene sets via biomaRt")
        } else if (!is.na(org_info$kegg)) {
            kegg_sets <- tryCatch({
                generate_gmt_from_kegg(org_info$kegg, verbose = TRUE)
            }, error = function(e) {
                message("  KEGG GMT generation failed: ", e$message)
                list()
            })
            if (length(kegg_sets) > 0) {
                gene_sets$KEGG <- kegg_sets
                message("Generated ", length(kegg_sets), " KEGG gene sets")
            }
        }

        if (length(gene_sets) == 0) {
            warning("Organism '", organism, "' not supported for standard pathway databases. ",
                    "Provide a GMT file for pathway analysis.")
            return(gene_sets)
        }
    }

    if (org_info$supported) {

        # GO terms via OrgDb
        if ("GO" %in% pathway_database) {
            tryCatch({
                if (!is.na(org_info$orgdb) && requireNamespace(org_info$orgdb, quietly = TRUE)) {
                    orgdb <- getExportedValue(org_info$orgdb, org_info$orgdb)

                    # Use SYMBOL keytype for proteomics (gene symbols),
                    # ENSEMBL keytype for RNA-seq (Ensembl gene IDs)
                    go_keytype <- if (target_id_type == "symbol") "SYMBOL" else "ENSEMBL"

                    # Verify the keytype is available in this OrgDb
                    available_keytypes <- AnnotationDbi::keytypes(orgdb)
                    if (!go_keytype %in% available_keytypes) {
                        message("Keytype '", go_keytype, "' not available in OrgDb. ",
                                "Falling back to ENSEMBL.")
                        go_keytype <- "ENSEMBL"
                    }

                    go_bp <- AnnotationDbi::select(
                        orgdb,
                        keys = AnnotationDbi::keys(orgdb, keytype = go_keytype),
                        columns = c(go_keytype, "GO"),
                        keytype = go_keytype
                    )

                    go_bp <- go_bp[!is.na(go_bp$GO), ]
                    go_sets <- split(go_bp[[go_keytype]], go_bp$GO)
                    go_sets <- go_sets[lengths(go_sets) >= 10 & lengths(go_sets) <= 500]

                    gene_sets$GO <- go_sets
                    message("Loaded ", length(go_sets), " GO gene sets (keytype: ",
                            go_keytype, ")")
                }
            }, error = function(e) {
                warning("Failed to load GO terms: ", e$message)
            })
        }

        # KEGG pathways
        if ("KEGG" %in% pathway_database && !is.na(org_info$kegg)) {
            tryCatch({
                # Cache clusterProfiler KEGG download (network call)
                kegg_cache_dir <- file.path(tempdir(), "kegg_cache")
                kegg_cp_cache <- file.path(kegg_cache_dir, paste0("kegg_cp_", org_info$kegg, ".rds"))
                if (file.exists(kegg_cp_cache) &&
                    difftime(Sys.time(), file.mtime(kegg_cp_cache), units = "days") < 7) {
                    kegg_gene <- readRDS(kegg_cp_cache)
                    message("Using cached clusterProfiler KEGG data for '", org_info$kegg, "'")
                } else {
                    kegg_gene <- clusterProfiler::download_KEGG(org_info$kegg, keggType = "KEGG")
                    dir.create(kegg_cache_dir, recursive = TRUE, showWarnings = FALSE)
                    tryCatch(saveRDS(kegg_gene, kegg_cp_cache), error = function(e) NULL)
                }

                if (!is.null(kegg_gene) && nrow(kegg_gene$KEGGPATHID2EXTID) > 0) {
                    kegg_df <- kegg_gene$KEGGPATHID2EXTID

                    if (!is.null(annotation) && "entrez_id" %in% colnames(annotation)) {
                        entrez_to_gene <- annotation[!is.na(annotation$entrez_id),
                                                     c("gene_id", "entrez_id")]
                        kegg_df <- merge(kegg_df, entrez_to_gene,
                                         by.x = "to", by.y = "entrez_id")
                        kegg_sets <- split(kegg_df$gene_id, kegg_df$from)
                    } else {
                        kegg_sets <- split(kegg_df$to, kegg_df$from)
                    }

                    kegg_sets <- kegg_sets[lengths(kegg_sets) >= 10 & lengths(kegg_sets) <= 500]
                    gene_sets$KEGG <- kegg_sets
                    message("Loaded ", length(kegg_sets), " KEGG pathways")
                }
            }, error = function(e) {
                warning("Failed to load KEGG pathways via clusterProfiler: ", e$message)
            })

            # Fallback: fetch KEGG directly via KEGGREST (works for non-model organisms
            # where clusterProfiler fails due to missing Entrez ID mappings)
            if (is.null(gene_sets$KEGG) || length(gene_sets$KEGG) == 0) {
                message("Trying KEGGREST fallback for organism code '", org_info$kegg, "'...")
                kegg_sets <- tryCatch(
                    fetch_kegg_via_rest(org_info$kegg),
                    error = function(e) {
                        warning("KEGGREST fallback also failed: ", e$message)
                        list()
                    }
                )
                if (length(kegg_sets) > 0) {
                    gene_sets$KEGG <- kegg_sets
                    message("Loaded ", length(kegg_sets), " KEGG pathways via KEGGREST")
                }
            }
        }
    }

    if (length(gene_sets) == 0) {
        warning("No gene sets loaded. Pathway analysis will be skipped.")
    }

    gene_sets
}

# ==============================================================================
# OVER-REPRESENTATION ANALYSIS
# ==============================================================================

#' Run over-representation analysis (Fisher's exact test)
#'
#' @param sig_genes Significant gene IDs
#' @param gene_sets Named list of gene sets
#' @param background Background gene IDs
#' @param min_size Minimum gene set size
#' @param max_size Maximum gene set size
#' @return Data frame of ORA results
#' @export
run_ora <- function(sig_genes, gene_sets, background, min_size = 10, max_size = 500) {

    gs_sizes <- lengths(gene_sets)
    gs_filtered <- gene_sets[gs_sizes >= min_size & gs_sizes <= max_size]

    if (length(gs_filtered) == 0) return(data.frame())

    results <- lapply(names(gs_filtered), function(gs_name) {
        gs_genes <- intersect(gs_filtered[[gs_name]], background)
        sig_in_gs <- length(intersect(sig_genes, gs_genes))
        sig_not_gs <- length(sig_genes) - sig_in_gs
        gs_not_sig <- length(gs_genes) - sig_in_gs
        neither <- length(background) - length(sig_genes) - gs_not_sig

        mat <- matrix(c(sig_in_gs, gs_not_sig, sig_not_gs, neither), nrow = 2)
        if (any(mat < 0)) return(NULL)

        test <- fisher.test(mat, alternative = "greater")

        data.frame(
            pathway = gs_name,
            size = length(gs_genes),
            overlap = sig_in_gs,
            pvalue = test$p.value,
            odds_ratio = test$estimate,
            stringsAsFactors = FALSE
        )
    })

    results_df <- do.call(rbind, results[!sapply(results, is.null)])

    if (nrow(results_df) > 0) {
        results_df$padj <- p.adjust(results_df$pvalue, method = "BH")
        results_df <- results_df[order(results_df$pvalue), ]
    }

    results_df
}

# ==============================================================================
# KEGG VIA KEGGREST (FOR NON-MODEL ORGANISMS)
# ==============================================================================

#' Fetch KEGG pathways directly via KEGGREST for any organism with a KEGG code
#'
#' @param kegg_code Three-letter KEGG organism code (e.g. "gla", "hsa")
#' @param min_size Minimum genes per pathway (default 3)
#' @return Named list of gene sets (character vectors of gene IDs)
fetch_kegg_via_rest <- function(kegg_code, min_size = 3, cache_dir = NULL,
                                cache_days = 7) {
    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop("KEGGREST package required. Install with: BiocManager::install('KEGGREST')")
    }

    # --- Local cache: avoid hitting KEGG REST API every run ---
    if (is.null(cache_dir)) {
        cache_dir <- file.path(tempdir(), "kegg_cache")
    }
    cache_file <- file.path(cache_dir, paste0("kegg_", kegg_code, ".rds"))

    if (file.exists(cache_file)) {
        cache_age <- difftime(Sys.time(), file.mtime(cache_file), units = "days")
        if (cache_age < cache_days) {
            cached <- readRDS(cache_file)
            message("Using cached KEGG data for '", kegg_code,
                    "' (", round(as.numeric(cache_age), 1), " days old)")
            return(cached)
        }
    }

    # Bulk fetch: all gene <-> pathway links for this organism
    links <- KEGGREST::keggLink("pathway", kegg_code)
    if (length(links) == 0) return(list())

    # links: names = "org:GENEID", values = "path:orgNNNNN"
    gene_ids <- sub(paste0("^", kegg_code, ":"), "", names(links))
    pathway_ids <- sub("^path:", "", unname(links))

    pw_genes <- split(gene_ids, pathway_ids)
    pw_genes <- pw_genes[lengths(pw_genes) >= min_size]

    # Fetch pathway names for readable labels
    pw_names <- tryCatch({
        pw_list <- KEGGREST::keggList("pathway", kegg_code)
        setNames(sub(" - .*$", "", pw_list), names(pw_list))
    }, error = function(e) NULL)

    # Rename from IDs to readable names where possible
    if (!is.null(pw_names)) {
        new_names <- pw_names[names(pw_genes)]
        has_name <- !is.na(new_names)
        names(pw_genes)[has_name] <- paste0(names(pw_genes)[has_name], " ", new_names[has_name])
    }

    # Save to cache
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    tryCatch(saveRDS(pw_genes, cache_file), error = function(e) NULL)

    pw_genes
}

# ==============================================================================
# GO TERM SIMPLIFICATION (via rrvgo semantic similarity)
# ==============================================================================

#' Simplify GO enrichment results by clustering semantically similar terms
#'
#' Uses rrvgo to compute pairwise semantic similarity between GO terms,
#' clusters them, and retains only the most significant term per cluster.
#' Works on any data.frame with GO IDs and p-values (fGSEA or ORA output).
#'
#' @param results_df Data frame with GO enrichment results. Must contain a column
#'   with GO IDs and a significance column (pvalue or padj).
#' @param go_col Column name containing GO IDs (default "pathway")
#' @param pval_col Column name for significance scoring. The term with the
#'   lowest value in each cluster is kept (default "pvalue")
#' @param orgdb Character: OrgDb package name (e.g. "org.Hs.eg.db")
#' @param ontology Character: GO ontology — "BP", "MF", or "CC"
#' @param threshold Numeric 0–1: similarity cutoff for clustering.
#'   Higher = more aggressive merging (default 0.7)
#' @param measure Semantic similarity method: "Wang", "Rel", "Resnik",
#'   "Lin", or "Jiang" (default "Wang")
#' @return A list with components:
#'   \item{reduced}{Data frame: filtered to representative terms only,
#'     with added columns cluster, parent_term, parent_term_name}
#'   \item{full}{Original data frame with cluster/parent annotations added}
#'   \item{sim_matrix}{The pairwise semantic similarity matrix}
simplify_go_results <- function(results_df,
                                go_col    = "pathway",
                                pval_col  = "pvalue",
                                sig_col   = NULL,
                                sig_cutoff = 0.10,
                                max_terms  = 500,
                                orgdb     = "org.Hs.eg.db",
                                ontology  = "BP",
                                threshold = 0.7,
                                measure   = "Wang") {

    if (is.null(results_df) || nrow(results_df) == 0) {
        return(list(reduced = results_df, full = results_df, sim_matrix = NULL))
    }

    # Validate GO column exists
    if (!go_col %in% colnames(results_df)) {
        warning("simplify_go_results: column '", go_col, "' not found. Returning unmodified.")
        return(list(reduced = results_df, full = results_df, sim_matrix = NULL))
    }

    # Filter to valid GO IDs only
    go_ids <- results_df[[go_col]]
    is_go <- grepl("^GO:[0-9]+$", go_ids)
    if (sum(is_go) < 2) {
        message("simplify_go_results: fewer than 2 GO terms — skipping simplification.")
        return(list(reduced = results_df, full = results_df, sim_matrix = NULL))
    }

    # Check for rrvgo
    if (!requireNamespace("rrvgo", quietly = TRUE)) {
        warning("rrvgo package not installed. Install with: BiocManager::install('rrvgo'). ",
                "Returning unsimplified results.")
        return(list(reduced = results_df, full = results_df, sim_matrix = NULL))
    }

    # Check for GOSemSim (required by rrvgo)
    if (!requireNamespace("GOSemSim", quietly = TRUE)) {
        warning("GOSemSim package not installed. Install with: BiocManager::install('GOSemSim'). ",
                "Returning unsimplified results.")
        return(list(reduced = results_df, full = results_df, sim_matrix = NULL))
    }

    # Load OrgDb
    if (!requireNamespace(orgdb, quietly = TRUE)) {
        warning("OrgDb package '", orgdb, "' not available. Returning unsimplified results.")
        return(list(reduced = results_df, full = results_df, sim_matrix = NULL))
    }

    orgdb_obj <- getExportedValue(orgdb, orgdb)

    # --- Performance optimisation: only simplify significant/top terms ---
    # Computing n×n Wang similarity is O(n²) and very slow for thousands of
    # terms. Pre-filter to significant terms, then cap at max_terms.
    go_subset <- results_df[is_go, , drop = FALSE]

    # Determine significance column (padj preferred, fall back to pval_col)
    eff_sig_col <- sig_col
    if (is.null(eff_sig_col)) {
        eff_sig_col <- intersect(c("padj", "p.adjust", pval_col), colnames(go_subset))[1]
    }
    if (!is.null(eff_sig_col) && eff_sig_col %in% colnames(go_subset)) {
        sig_vals <- as.numeric(go_subset[[eff_sig_col]])
        keep_sig <- !is.na(sig_vals) & sig_vals < sig_cutoff
        if (sum(keep_sig) >= 2) {
            go_subset <- go_subset[keep_sig, , drop = FALSE]
            message(sprintf("GO simplification: pre-filtered to %d significant terms (%s < %.2f)",
                            nrow(go_subset), eff_sig_col, sig_cutoff))
        }
    }

    # Cap at max_terms (take top by pvalue)
    if (nrow(go_subset) > max_terms) {
        go_subset <- go_subset[order(go_subset[[pval_col]]), ]
        go_subset <- go_subset[seq_len(max_terms), , drop = FALSE]
        message(sprintf("GO simplification: capped at %d terms for similarity computation", max_terms))
    }

    go_terms <- go_subset[[go_col]]

    if (length(go_terms) < 2) {
        message("simplify_go_results: fewer than 2 GO terms after filtering — skipping.")
        return(list(reduced = results_df, full = results_df, sim_matrix = NULL))
    }

    # Build score vector (lower pvalue = higher priority, so use -log10)
    scores <- setNames(-log10(go_subset[[pval_col]] + 1e-300), go_terms)

    # Compute semantic similarity matrix
    message(sprintf("Computing %s similarity for %d GO terms...", measure, length(go_terms)))
    sim_matrix <- tryCatch({
        rrvgo::calculateSimMatrix(
            go_terms,
            orgdb     = orgdb_obj,
            ont       = ontology,
            method    = measure
        )
    }, error = function(e) {
        warning("Semantic similarity computation failed: ", e$message)
        NULL
    })

    if (is.null(sim_matrix) || nrow(sim_matrix) < 2) {
        message("simplify_go_results: similarity matrix too small — skipping.")
        return(list(reduced = results_df, full = results_df, sim_matrix = NULL))
    }

    # Tag the matrix with orgdb name so downstream plotting can reuse it
    attr(sim_matrix, "orgdb") <- orgdb

    # Reduce: cluster terms and pick representatives
    reduced_terms <- tryCatch({
        rrvgo::reduceSimMatrix(
            sim_matrix,
            scores    = scores[rownames(sim_matrix)],
            threshold = threshold,
            orgdb     = orgdb_obj
        )
    }, error = function(e) {
        warning("reduceSimMatrix failed: ", e$message)
        NULL
    })

    if (is.null(reduced_terms)) {
        return(list(reduced = results_df, full = results_df, sim_matrix = sim_matrix))
    }

    # reduced_terms has columns: go, cluster, parent, parentTerm, score, ...
    # Merge cluster info back into original results
    cluster_map <- data.frame(
        go_id        = reduced_terms$go,
        cluster      = reduced_terms$cluster,
        parent_term  = reduced_terms$parent,
        parent_term_name = reduced_terms$parentTerm,
        stringsAsFactors = FALSE
    )

    results_annotated <- merge(
        results_df,
        cluster_map,
        by.x    = go_col,
        by.y    = "go_id",
        all.x   = TRUE,
        sort    = FALSE
    )

    # Restore original row order
    results_annotated <- results_annotated[order(results_annotated[[pval_col]]), ]

    # Keep only the best (lowest pval) term per cluster
    go_annotated <- results_annotated[!is.na(results_annotated$cluster), , drop = FALSE]
    non_go       <- results_annotated[is.na(results_annotated$cluster), , drop = FALSE]

    # Within each cluster, keep the row with the minimum pvalue
    representatives <- do.call(rbind, lapply(
        split(go_annotated, go_annotated$cluster),
        function(grp) grp[which.min(grp[[pval_col]]), , drop = FALSE]
    ))

    reduced_df <- rbind(representatives, non_go)
    reduced_df <- reduced_df[order(reduced_df[[pval_col]]), ]

    n_before <- sum(is_go)
    n_after  <- nrow(representatives)
    message(sprintf("GO simplification: %d terms -> %d representative terms (threshold=%.2f, %s)",
                    n_before, n_after, threshold, measure))

    list(
        reduced    = reduced_df,
        full       = results_annotated,
        sim_matrix = sim_matrix
    )
}

#' Generate treemap plot of simplified GO term clusters
#'
#' @param reduced_terms Output from rrvgo::reduceSimMatrix()
#' @param output_file Path for the output PNG file
#' @param title Plot title
#' @return The output file path (invisibly)
plot_go_treemap <- function(sim_matrix, scores, threshold, orgdb_obj,
                            output_file, title = "GO Term Clusters") {

    if (!requireNamespace("rrvgo", quietly = TRUE)) return(invisible(NULL))

    reduced_terms <- tryCatch(
        rrvgo::reduceSimMatrix(sim_matrix, scores = scores,
                               threshold = threshold, orgdb = orgdb_obj),
        error = function(e) NULL
    )
    if (is.null(reduced_terms)) return(invisible(NULL))

    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

    tryCatch({
        png(output_file, width = 12, height = 8, units = "in", res = 150)
        rrvgo::treemapPlot(reduced_terms)
        title(main = title, cex.main = 1.2)
        dev.off()
        message("  Saved GO treemap: ", output_file)
    }, error = function(e) {
        tryCatch(dev.off(), error = function(x) NULL)
        warning("Treemap plot failed: ", e$message)
    })

    invisible(output_file)
}
