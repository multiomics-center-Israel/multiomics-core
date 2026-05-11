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

    # Expand "all" shorthand to all supported databases
    if ("all" %in% pathway_database) {
        pathway_database <- c("GO", "KEGG", "Reactome")
    }

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

        # GO terms via OrgDb — split by ontology (BP, CC, MF)
        # Accept "GO" (loads all 3 ontologies separately) or specific "GO_BP", "GO_CC", "GO_MF"
        go_requested <- any(grepl("^GO", pathway_database))
        if (go_requested) {
            tryCatch({
                if (!is.na(org_info$orgdb) && requireNamespace(org_info$orgdb, quietly = TRUE)) {
                    orgdb <- getExportedValue(org_info$orgdb, org_info$orgdb)

                    go_keytype <- if (target_id_type == "symbol") "SYMBOL" else "ENSEMBL"

                    available_keytypes <- AnnotationDbi::keytypes(orgdb)
                    if (!go_keytype %in% available_keytypes) {
                        message("Keytype '", go_keytype, "' not available in OrgDb. ",
                                "Falling back to ENSEMBL.")
                        go_keytype <- "ENSEMBL"
                    }

                    go_all <- AnnotationDbi::select(
                        orgdb,
                        keys = AnnotationDbi::keys(orgdb, keytype = go_keytype),
                        columns = c(go_keytype, "GO", "ONTOLOGY"),
                        keytype = go_keytype
                    )
                    go_all <- go_all[!is.na(go_all$GO) & !is.na(go_all$ONTOLOGY), ]

                    # Determine which GO ontologies to load
                    if ("GO" %in% pathway_database) {
                        go_onts <- c("BP", "CC", "MF")
                    } else {
                        go_onts <- sub("^GO_", "", grep("^GO_", pathway_database, value = TRUE))
                    }

                    for (ont in go_onts) {
                        go_sub <- go_all[go_all$ONTOLOGY == ont, ]
                        if (nrow(go_sub) == 0) next
                        go_sets <- split(go_sub[[go_keytype]], go_sub$GO)
                        go_sets <- go_sets[lengths(go_sets) >= 10 & lengths(go_sets) <= 500]
                        gs_name <- paste0("GO_", ont)
                        gene_sets[[gs_name]] <- go_sets
                        message("Loaded ", length(go_sets), " ", gs_name,
                                " gene sets (keytype: ", go_keytype, ")")
                    }
                }
            }, error = function(e) {
                warning("Failed to load GO terms: ", e$message)
            })
        }

        # Reactome pathways via ReactomePA
        if ("Reactome" %in% pathway_database) {
            tryCatch({
                if (requireNamespace("ReactomePA", quietly = TRUE)) {
                    orgdb <- getExportedValue(org_info$orgdb, org_info$orgdb)

                    # ReactomePA needs Entrez IDs — build gene set lists keyed by target_id_type
                    all_keys <- AnnotationDbi::keys(orgdb, keytype = if (target_id_type == "symbol") "SYMBOL" else "ENSEMBL")
                    id_map <- tryCatch(
                        AnnotationDbi::select(orgdb,
                            keys = all_keys,
                            keytype = if (target_id_type == "symbol") "SYMBOL" else "ENSEMBL",
                            columns = c(if (target_id_type == "symbol") "SYMBOL" else "ENSEMBL", "ENTREZID")),
                        error = function(e) NULL
                    )

                    if (!is.null(id_map) && nrow(id_map) > 0) {
                        id_map <- id_map[!is.na(id_map$ENTREZID), ]
                        entrez_to_target <- setNames(
                            id_map[[if (target_id_type == "symbol") "SYMBOL" else "ENSEMBL"]],
                            id_map$ENTREZID
                        )

                        # Fetch Reactome pathway-gene mappings
                        reactome_organism <- gsub(" ", "_", organism)
                        reactome_db <- tryCatch(
                            ReactomePA:::get_Reactome_Env(),
                            error = function(e) NULL
                        )
                        if (is.null(reactome_db)) {
                            # Fallback: use enrichPathway on a dummy gene to trigger DB download,
                            # then extract the gene sets
                            reactome_gene2path <- tryCatch({
                                rpa <- reactome.db::reactomeEXTID2PATHID
                                rpa_df <- AnnotationDbi::toTable(rpa)
                                rpa_df
                            }, error = function(e) NULL)
                        } else {
                            reactome_gene2path <- tryCatch(
                                AnnotationDbi::toTable(reactome_db),
                                error = function(e) NULL
                            )
                        }

                        # Build Reactome gene sets from reactome.db directly
                        reactome_sets <- tryCatch({
                            if (requireNamespace("reactome.db", quietly = TRUE)) {
                                path2ext <- AnnotationDbi::toTable(reactome.db::reactomePATHID2EXTID)
                                path2name <- AnnotationDbi::toTable(reactome.db::reactomePATHID2NAME)
                                # Filter to organism
                                org_key <- organism
                                path2name_org <- path2name[grepl(org_key, path2name$path_name, ignore.case = TRUE), ]
                                if (nrow(path2name_org) == 0) path2name_org <- path2name
                                path2ext_org <- path2ext[path2ext$DB_ID %in% path2name_org$DB_ID, ]

                                # Map Entrez to target ID type
                                path2ext_org$target_id <- entrez_to_target[path2ext_org$gene_id]
                                path2ext_org <- path2ext_org[!is.na(path2ext_org$target_id), ]

                                # Create named gene sets
                                r_sets <- split(path2ext_org$target_id, path2ext_org$DB_ID)
                                r_sets <- r_sets[lengths(r_sets) >= 10 & lengths(r_sets) <= 500]

                                # Use pathway names instead of IDs
                                name_map <- setNames(path2name_org$path_name, path2name_org$DB_ID)
                                name_map <- sub(paste0("^", org_key, ": "), "", name_map)
                                new_names <- name_map[names(r_sets)]
                                has_name <- !is.na(new_names) & nzchar(new_names)
                                names(r_sets)[has_name] <- new_names[has_name]
                                r_sets
                            } else NULL
                        }, error = function(e) {
                            message("reactome.db extraction failed: ", e$message)
                            NULL
                        })

                        if (!is.null(reactome_sets) && length(reactome_sets) > 0) {
                            gene_sets$Reactome <- reactome_sets
                            message("Loaded ", length(reactome_sets), " Reactome pathways")
                        }
                    }
                } else {
                    message("ReactomePA not installed. Install with: BiocManager::install('ReactomePA')")
                }
            }, error = function(e) {
                warning("Failed to load Reactome pathways: ", e$message)
            })
        }

        # KEGG pathways
        if ("KEGG" %in% pathway_database && !is.na(org_info$kegg)) {
            tryCatch({
                kegg_gene <- clusterProfiler::download_KEGG(org_info$kegg, keggType = "KEGG")

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
fetch_kegg_via_rest <- function(kegg_code, min_size = 3) {
    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop("KEGGREST package required. Install with: BiocManager::install('KEGGREST')")
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

    pw_genes
}

# ==============================================================================
# PATHWAY ANALYSIS (fGSEA + ORA + Plots)
# Shared pathway helpers used by rnaseq and proteomics pipelines.
# Expects generic DE table interface (FeatureID, log2FoldChange, pvalue, padj, stat).
# ==============================================================================

#' Look up GO term names from GO IDs
#'
#' @param go_ids Character vector of GO IDs (e.g., "GO:0000001")
#' @return Named vector with GO IDs as names and term names as values
lookup_go_term_names <- function(go_ids) {
    if (length(go_ids) == 0) return(character(0))

    # Initialize with IDs as fallback
    term_names <- setNames(go_ids, go_ids)

    # Try GO.db package first
    if (requireNamespace("GO.db", quietly = TRUE)) {
        tryCatch({
            # Filter to valid GO IDs present in GO.db (some may be obsolete)
            valid_keys <- AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID")
            valid_ids <- go_ids[go_ids %in% valid_keys]
            if (length(valid_ids) > 0) {
                go_terms <- AnnotationDbi::select(
                    GO.db::GO.db,
                    keys = valid_ids,
                    columns = "TERM",
                    keytype = "GOID"
                )
                # Update names where we got results
                matched <- match(go_ids, go_terms$GOID)
                has_term <- !is.na(matched) & !is.na(go_terms$TERM[matched])
                term_names[has_term] <- go_terms$TERM[matched[has_term]]
            }
        }, error = function(e) {
            message("Could not look up GO terms: ", e$message)
        })
    }

    term_names
}

#' Add pathway names to fGSEA/ORA results
#'
#' @param pathway_df Data frame with pathway analysis results (has 'pathway' column)
#' @param database Character string indicating the database ("GO", "KEGG", etc.)
#' @return Data frame with added 'pathway_name' column
add_pathway_names <- function(pathway_df, database, gene_sets = NULL) {
    if (is.null(pathway_df) || nrow(pathway_df) == 0) return(pathway_df)
    if (!"pathway" %in% colnames(pathway_df)) return(pathway_df)

    pathway_ids <- pathway_df$pathway

    if (database == "GO" || grepl("^GO", database, ignore.case = TRUE)) {
        # Look up GO term names
        names_vec <- lookup_go_term_names(pathway_ids)
        pathway_df$pathway_name <- unname(names_vec[pathway_ids])
    } else if (database == "KEGG" || grepl("KEGG", database, ignore.case = TRUE)) {
        # For KEGG, the pathway names are usually already descriptive
        # but we'll keep the ID for now as KEGG lookup is more complex
        pathway_df$pathway_name <- pathway_ids
    } else {
        # For custom GMT, use descriptions attribute if available
        descriptions <- if (!is.null(gene_sets)) attr(gene_sets, "descriptions") else NULL
        if (!is.null(descriptions)) {
            pathway_df$pathway_name <- unname(descriptions[pathway_ids])
        } else {
            pathway_df$pathway_name <- pathway_ids
        }
        # For entries where description is still a GO ID, look up the term name
        is_go_id <- grepl("^GO:[0-9]+$", pathway_df$pathway_name)
        if (any(is_go_id)) {
            go_names <- lookup_go_term_names(pathway_df$pathway_name[is_go_id])
            pathway_df$pathway_name[is_go_id] <- unname(go_names)
        }
    }

    # Reorder columns to put pathway_name right after pathway
    cols <- colnames(pathway_df)
    pathway_idx <- which(cols == "pathway")
    new_order <- c(
        cols[1:pathway_idx],
        "pathway_name",
        setdiff(cols[(pathway_idx + 1):length(cols)], "pathway_name")
    )
    pathway_df <- pathway_df[, new_order, drop = FALSE]

    pathway_df
}

#' Run pathway analysis on DE results
#'
#' @param de_tables Named list of DE result data frames (from run_deseq2_de()$tables).
#'   Each table has columns: FeatureID, log2FoldChange, pvalue, padj, stat, etc.
#' @param gene_sets Named list of gene set collections from load_gene_sets()
#' @param annotation Gene annotation data frame (gene_id, symbol, entrez_id)
#' @param method "fgsea", "ora", or "both"
#' @param min_size Minimum gene set size for fGSEA
#' @param max_size Maximum gene set size for fGSEA
#' @return Named list (by contrast) of named lists (by db+method) of result data frames
#' @export
run_pathway_analysis <- function(de_tables,
                                  gene_sets,
                                  annotation = NULL,
                                  method = "fgsea",
                                  min_size = 10,
                                  max_size = 500) {

    if (length(gene_sets) == 0) {
        message("No gene sets available. Skipping pathway analysis.")
        return(list())
    }

    pathway_results <- list()

    for (contrast_name in names(de_tables)) {
        res <- de_tables[[contrast_name]]
        if (is.null(res) || nrow(res) == 0) next

        message("Running pathway analysis for: ", contrast_name)
        contrast_results <- list()

        for (db_name in names(gene_sets)) {
            gs <- gene_sets[[db_name]]
            if (length(gs) == 0) next

            message("  Database: ", db_name, " (", length(gs), " gene sets)")

            tryCatch({
                # ---- fGSEA ----
                if (method %in% c("fgsea", "both")) {

                    # Build ranked gene list from DE table
                    # Prefer stat column (Wald statistic); fallback to sign(lfc)*-log10(p)
                    if ("stat" %in% colnames(res)) {
                        ranks <- setNames(res$stat, res$FeatureID)
                    } else {
                        ranks <- setNames(
                            sign(res$log2FoldChange) * -log10(res$pvalue + 1e-300),
                            res$FeatureID
                        )
                    }

                    ranks <- ranks[!is.na(ranks)]
                    ranks <- sort(ranks, decreasing = TRUE)

                    fgsea_res <- fgsea::fgsea(
                        pathways = gs,
                        stats = ranks,
                        minSize = min_size,
                        maxSize = max_size,
                        nPermSimple = 10000
                    )

                    fgsea_df <- as.data.frame(fgsea_res)

                    if (nrow(fgsea_df) == 0) {
                        message("    fgsea: no gene set overlap — skipping ", db_name)
                    } else {
                        fgsea_df$contrast <- contrast_name
                        fgsea_df$database <- db_name
                        fgsea_df$method <- "fgsea"

                        # Collapse leading edge to comma-separated string
                        fgsea_df$leadingEdge <- sapply(fgsea_df$leadingEdge, paste, collapse = ",")

                        # Add pathway names (GO term names, etc.)
                        fgsea_df <- add_pathway_names(fgsea_df, db_name, gs)

                        contrast_results[[paste0(db_name, "_fgsea")]] <- fgsea_df

                        n_sig <- sum(fgsea_df$padj < 0.05, na.rm = TRUE)
                        message("    fgsea: ", n_sig, " significant pathways (padj < 0.05)")
                    }
                }

                # ---- ORA ----
                if (method %in% c("ora", "both")) {

                    # Identify significant up/down genes
                    de_cfg_padj <- 0.05
                    de_cfg_lfc  <- log2(1.5)
                    sig_up   <- res$FeatureID[!is.na(res$padj) &
                                              res$padj < de_cfg_padj &
                                              res$log2FoldChange > de_cfg_lfc]
                    sig_down <- res$FeatureID[!is.na(res$padj) &
                                              res$padj < de_cfg_padj &
                                              res$log2FoldChange < -de_cfg_lfc]
                    background <- res$FeatureID

                    if (length(sig_up) >= 5) {
                        ora_up <- run_ora(sig_up, gs, background, min_size, max_size)
                        if (nrow(ora_up) > 0) {
                            ora_up$contrast <- contrast_name
                            ora_up$database <- db_name
                            ora_up$method <- "ora"
                            ora_up$direction <- "up"
                            ora_up <- add_pathway_names(ora_up, db_name, gs)
                            contrast_results[[paste0(db_name, "_ora_up")]] <- ora_up
                        }
                    }

                    if (length(sig_down) >= 5) {
                        ora_down <- run_ora(sig_down, gs, background, min_size, max_size)
                        if (nrow(ora_down) > 0) {
                            ora_down$contrast <- contrast_name
                            ora_down$database <- db_name
                            ora_down$method <- "ora"
                            ora_down$direction <- "down"
                            ora_down <- add_pathway_names(ora_down, db_name, gs)
                            contrast_results[[paste0(db_name, "_ora_down")]] <- ora_down
                        }
                    }
                }
            }, error = function(e) {
                warning("Pathway analysis failed for ", db_name, ": ", e$message)
            })
        }

        pathway_results[[contrast_name]] <- contrast_results
    }

    pathway_results
}

#' Save pathway results to CSV files
#'
#' @param pathway_results Output of run_pathway_analysis()
#' @param output_dir Directory to write CSV files
#' @return Character vector of output file paths
#' @export
save_pathway_results <- function(pathway_results, output_dir) {

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_files <- c()

    for (contrast_name in names(pathway_results)) {
        contrast_res <- pathway_results[[contrast_name]]

        for (analysis_name in names(contrast_res)) {
            res_df <- contrast_res[[analysis_name]]
            if (is.null(res_df) || nrow(res_df) == 0) next

            clean_contrast <- gsub("[^a-zA-Z0-9_-]", "_", contrast_name)
            clean_analysis <- gsub("[^a-zA-Z0-9_-]", "_", analysis_name)

            output_file <- file.path(output_dir,
                                     paste0("pathway_", clean_contrast, "_", clean_analysis, ".csv"))
            write.csv(res_df, output_file, row.names = FALSE)
            output_files <- c(output_files, output_file)
        }
    }

    message("Saved ", length(output_files), " pathway result file(s) to ", output_dir)
    output_files
}

# ==============================================================================
# SEMANTIC CLUSTERING OF ENRICHMENT TERMS (rrvgo)
# ==============================================================================

#' Cluster enrichment terms to reduce redundancy
#'
#' For GO terms, uses rrvgo (semantic similarity via GOSemSim).
#' For KEGG/custom terms, uses Jaccard similarity on gene overlap.
#'
#' @param enrichment_df Data frame with enrichment results (must have 'pathway' and 'padj' columns)
#' @param database Character: "GO", "KEGG", or "custom"
#' @param gene_sets Named list of gene sets (needed for Jaccard clustering of non-GO terms)
#' @param organism Character: organism name for OrgDb lookup (needed for GO clustering)
#' @param threshold Numeric: similarity threshold for merging (0-1, default 0.7). Lower = more aggressive merging.
#' @param ont Character: GO ontology for rrvgo ("BP", "MF", "CC"). Default "BP".
#' @return Data frame with added columns: cluster, parentTerm, parentPadj
#' @export
cluster_enrichment_terms <- function(enrichment_df,
                                      database,
                                      gene_sets = NULL,
                                      organism = "Homo sapiens",
                                      threshold = 0.7,
                                      ont = "BP") {

    if (is.null(enrichment_df) || nrow(enrichment_df) == 0) return(enrichment_df)
    if (!"pathway" %in% colnames(enrichment_df)) return(enrichment_df)

    # Only cluster significant terms
    sig <- enrichment_df[!is.na(enrichment_df$padj) & enrichment_df$padj < 0.05, ]
    if (nrow(sig) < 2) return(NULL)

    is_go <- grepl("^GO", database, ignore.case = TRUE) ||
        all(grepl("^GO:[0-9]+", sig$pathway))

    if (is_go) {
        clustered <- .cluster_go_terms(sig, organism, threshold, ont)
    } else {
        clustered <- .cluster_by_jaccard(sig, gene_sets, threshold)
    }

    clustered
}

#' Cluster GO terms using rrvgo semantic similarity
#' @keywords internal
.cluster_go_terms <- function(sig_df, organism, threshold, ont) {
    if (!requireNamespace("rrvgo", quietly = TRUE)) {
        message("rrvgo package not installed. Install with: BiocManager::install('rrvgo')")
        message("Falling back to Jaccard-based clustering.")
        return(NULL)
    }

    org_info <- get_organism_info(organism)
    if (is.na(org_info$orgdb)) {
        message("No OrgDb available for '", organism, "'. Cannot compute GO semantic similarity.")
        return(NULL)
    }

    if (!requireNamespace(org_info$orgdb, quietly = TRUE)) {
        message("OrgDb package '", org_info$orgdb, "' not installed.")
        return(NULL)
    }

    go_ids <- sig_df$pathway
    # Ensure all are valid GO IDs
    valid_go <- grepl("^GO:[0-9]+", go_ids)
    if (sum(valid_go) < 2) {
        message("Fewer than 2 valid GO IDs found. Skipping GO clustering.")
        return(NULL)
    }

    sig_df <- sig_df[valid_go, ]
    go_ids <- sig_df$pathway
    scores <- setNames(-log10(sig_df$padj + 1e-300), go_ids)

    sim_matrix <- tryCatch({
        rrvgo::calculateSimMatrix(
            go_ids,
            orgdb = org_info$orgdb,
            ont = ont,
            method = "Rel"
        )
    }, error = function(e) {
        message("rrvgo similarity matrix computation failed: ", e$message)
        NULL
    })

    if (is.null(sim_matrix) || nrow(sim_matrix) < 2) return(NULL)

    # Only keep GO IDs that appear in the similarity matrix
    common_ids <- intersect(go_ids, rownames(sim_matrix))
    if (length(common_ids) < 2) return(NULL)
    scores <- scores[common_ids]

    reduced <- tryCatch({
        rrvgo::reduceSimMatrix(sim_matrix, scores, threshold = threshold,
                               orgdb = org_info$orgdb)
    }, error = function(e) {
        message("rrvgo reduceSimMatrix failed: ", e$message)
        NULL
    })

    if (is.null(reduced) || nrow(reduced) == 0) return(NULL)

    # Build clustered result: merge reduced info back into enrichment data
    # reduced has columns: go, cluster, parent, parentSimScore, parentTerm, score
    sig_df <- sig_df[sig_df$pathway %in% reduced$go, ]
    sig_df$cluster <- reduced$cluster[match(sig_df$pathway, reduced$go)]
    sig_df$parentTerm <- reduced$parentTerm[match(sig_df$pathway, reduced$go)]
    sig_df$parent <- reduced$parent[match(sig_df$pathway, reduced$go)]

    # Build cluster summary: one row per cluster (representative = parent term)
    cluster_summary <- .build_cluster_summary(sig_df)

    cluster_summary
}

#' Cluster non-GO terms by Jaccard similarity on gene overlap
#' @keywords internal
.cluster_by_jaccard <- function(sig_df, gene_sets, threshold) {
    if (is.null(gene_sets) || length(gene_sets) == 0) {
        message("No gene sets provided for Jaccard clustering.")
        return(NULL)
    }

    pathways <- sig_df$pathway
    # Only use pathways that exist in gene_sets
    pathways <- pathways[pathways %in% names(gene_sets)]
    if (length(pathways) < 2) return(NULL)

    sig_df <- sig_df[sig_df$pathway %in% pathways, ]

    # Compute pairwise Jaccard similarity
    n <- length(pathways)
    sim_mat <- matrix(0, nrow = n, ncol = n, dimnames = list(pathways, pathways))

    for (i in seq_len(n)) {
        for (j in i:n) {
            a <- gene_sets[[pathways[i]]]
            b <- gene_sets[[pathways[j]]]
            inter <- length(intersect(a, b))
            union_size <- length(union(a, b))
            sim <- if (union_size > 0) inter / union_size else 0
            sim_mat[i, j] <- sim
            sim_mat[j, i] <- sim
        }
    }

    # Hierarchical clustering with complete linkage
    dist_mat <- as.dist(1 - sim_mat)
    hc <- hclust(dist_mat, method = "complete")
    clusters <- cutree(hc, h = 1 - threshold)

    sig_df$cluster <- clusters[sig_df$pathway]

    # Within each cluster, pick the term with the lowest padj as the representative
    sig_df$parentTerm <- NA_character_
    sig_df$parent <- NA_character_
    for (cl in unique(sig_df$cluster)) {
        cl_rows <- which(sig_df$cluster == cl)
        best <- cl_rows[which.min(sig_df$padj[cl_rows])]
        label <- if ("pathway_name" %in% colnames(sig_df) && nzchar(sig_df$pathway_name[best] %||% "")) {
            sig_df$pathway_name[best]
        } else {
            sig_df$pathway[best]
        }
        sig_df$parentTerm[cl_rows] <- label
        sig_df$parent[cl_rows] <- sig_df$pathway[best]
    }

    .build_cluster_summary(sig_df)
}

#' Build cluster summary table from clustered enrichment data
#' @keywords internal
.build_cluster_summary <- function(clustered_df) {
    clusters <- split(clustered_df, clustered_df$cluster)

    summaries <- lapply(clusters, function(cl_df) {
        # Representative is the parent term (lowest padj in cluster)
        rep_idx <- which.min(cl_df$padj)

        rep_name <- if ("pathway_name" %in% colnames(cl_df) && nzchar(cl_df$pathway_name[rep_idx] %||% "")) {
            cl_df$pathway_name[rep_idx]
        } else {
            cl_df$parentTerm[rep_idx]
        }

        # Collect member term names
        member_names <- if ("pathway_name" %in% colnames(cl_df)) {
            cl_df$pathway_name[cl_df$pathway_name != rep_name]
        } else {
            cl_df$pathway[cl_df$pathway != cl_df$pathway[rep_idx]]
        }

        # Carry over method-specific columns from the representative row
        rep_row <- cl_df[rep_idx, , drop = FALSE]

        data.frame(
            cluster_id = cl_df$cluster[1],
            representative_term = rep_name,
            representative_id = cl_df$pathway[rep_idx],
            n_members = nrow(cl_df),
            member_terms = paste(head(member_names, 10), collapse = "; "),
            padj = rep_row$padj,
            NES = if ("NES" %in% colnames(rep_row)) rep_row$NES else NA_real_,
            size = if ("size" %in% colnames(rep_row)) rep_row$size else NA_integer_,
            pvalue = if ("pvalue" %in% colnames(rep_row)) rep_row$pvalue else
                if ("pval" %in% colnames(rep_row)) rep_row$pval else NA_real_,
            stringsAsFactors = FALSE
        )
    })

    result <- do.call(rbind, summaries)
    rownames(result) <- NULL
    result <- result[order(result$padj), ]
    result
}

#' Generate pathway visualization dotplots
#'
#' @param pathway_results Output of run_pathway_analysis()
#' @param output_dir Directory for PNG files
#' @return Named list of plot file paths
#' @export
generate_pathway_plots <- function(pathway_results, output_dir) {

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    plot_files <- list()

    for (contrast_name in names(pathway_results)) {
        contrast_res <- pathway_results[[contrast_name]]
        clean_contrast <- gsub("[^a-zA-Z0-9_-]", "_", contrast_name)

        for (analysis_name in names(contrast_res)) {
            res_df <- contrast_res[[analysis_name]]
            if (is.null(res_df) || nrow(res_df) == 0) next

            clean_analysis <- gsub("[^a-zA-Z0-9_-]", "_", analysis_name)

            if ("padj" %in% colnames(res_df)) {
                top_pathways <- head(res_df[order(res_df$padj), ], 20)
            } else if ("pvalue" %in% colnames(res_df)) {
                top_pathways <- head(res_df[order(res_df$pvalue), ], 20)
            } else {
                next
            }

            if (nrow(top_pathways) < 3) next

            # Use pathway_name if available, otherwise fall back to pathway ID
            label_col <- if ("pathway_name" %in% colnames(top_pathways)) "pathway_name" else "pathway"
            # Truncate long names for display (max 60 characters)
            top_pathways$pathway_label <- substr(top_pathways[[label_col]], 1, 60)

            if ("NES" %in% colnames(top_pathways)) {
                p <- ggplot2::ggplot(top_pathways,
                                     ggplot2::aes(x = NES, y = reorder(pathway_label, NES))) +
                    ggplot2::geom_point(ggplot2::aes(size = size, color = -log10(padj))) +
                    ggplot2::scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
                    ggplot2::labs(
                        title = paste("Top Pathways:", contrast_name),
                        subtitle = analysis_name,
                        x = "Normalized Enrichment Score",
                        y = "",
                        size = "Gene Set Size"
                    ) +
                    ggplot2::theme_minimal() +
                    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
            } else {
                p <- ggplot2::ggplot(top_pathways,
                                     ggplot2::aes(x = -log10(pvalue),
                                                  y = reorder(pathway_label, -log10(pvalue)))) +
                    ggplot2::geom_point(ggplot2::aes(size = overlap, color = -log10(padj))) +
                    ggplot2::scale_color_gradient(low = "blue", high = "red", name = "-log10(padj)") +
                    ggplot2::labs(
                        title = paste("Top Pathways:", contrast_name),
                        subtitle = analysis_name,
                        x = "-log10(p-value)",
                        y = "",
                        size = "Overlap"
                    ) +
                    ggplot2::theme_minimal() +
                    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
            }

            plot_file <- file.path(output_dir,
                                   paste0("pathway_", clean_contrast, "_", clean_analysis, ".png"))
            ggplot2::ggsave(plot_file, p, width = 10, height = 8)
            plot_files[[paste0(clean_contrast, "_", clean_analysis)]] <- plot_file
        }
    }

    message("Generated ", length(plot_files), " pathway plots in ", output_dir)
    plot_files
}

#' Generate dotplots from clustered enrichment CSVs
#'
#' Reads clustered result CSVs and produces one dotplot per file showing
#' cluster representatives instead of redundant terms.
#'
#' @param clustered_dir Directory containing clustered CSV files
#' @param output_dir Directory for PNG output
#' @return Named list of plot file paths
#' @export
generate_clustered_dotplots <- function(clustered_dir, output_dir) {

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    plot_files <- list()

    csv_files <- list.files(clustered_dir, pattern = "\\.csv$", full.names = TRUE)
    if (length(csv_files) == 0) return(plot_files)

    for (cf in csv_files) {
        cl_df <- read.csv(cf, stringsAsFactors = FALSE)
        if (is.null(cl_df) || nrow(cl_df) < 2) next

        label_raw <- gsub("^clustered_|\\.csv$", "", basename(cf))
        top <- head(cl_df[order(cl_df$padj), ], 25)

        # Build display label: "term (n members)" for clusters with >1 member
        top$display_label <- ifelse(
            top$n_members > 1,
            paste0(substr(top$representative_term, 1, 50), " (", top$n_members, ")"),
            substr(top$representative_term, 1, 55)
        )

        has_nes <- "NES" %in% colnames(top) && !all(is.na(top$NES))

        if (has_nes) {
            p <- ggplot2::ggplot(top,
                    ggplot2::aes(x = NES, y = reorder(display_label, NES))) +
                ggplot2::geom_point(ggplot2::aes(
                    size = n_members,
                    color = -log10(padj))) +
                ggplot2::scale_color_gradient(low = "blue", high = "red",
                                              name = "-log10(padj)") +
                ggplot2::labs(
                    title = paste("Clustered Pathways:", gsub("_", " ", label_raw)),
                    x = "NES (representative term)",
                    y = "",
                    size = "# Merged\nTerms"
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
        } else {
            pval_col <- if ("pvalue" %in% colnames(top)) top$pvalue else top$padj
            top$.neg_log10_p <- -log10(pval_col + 1e-300)

            p <- ggplot2::ggplot(top,
                    ggplot2::aes(x = .neg_log10_p,
                                 y = reorder(display_label, .neg_log10_p))) +
                ggplot2::geom_point(ggplot2::aes(
                    size = n_members,
                    color = -log10(padj))) +
                ggplot2::scale_color_gradient(low = "blue", high = "red",
                                              name = "-log10(padj)") +
                ggplot2::labs(
                    title = paste("Clustered Pathways:", gsub("_", " ", label_raw)),
                    x = "-log10(p-value)",
                    y = "",
                    size = "# Merged\nTerms"
                ) +
                ggplot2::theme_minimal() +
                ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
        }

        plot_file <- file.path(output_dir, paste0(label_raw, ".png"))
        ggplot2::ggsave(plot_file, p, width = 10, height = 8)
        plot_files[[label_raw]] <- plot_file
    }

    message("Generated ", length(plot_files), " clustered dotplots in ", output_dir)
    plot_files
}
