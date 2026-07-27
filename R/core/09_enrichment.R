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


# ==============================================================================
# PHASE 1 — LOCAL OFFLINE ENRICHMENT (enrichment migration v2)
# Local table-driven KEGG/GO loading, ranked gene lists, GSEA, cluster-based ORA.
# Re-applied from reference branch claude/enrichment-migration-plan-cuduc.
# No online resources. GO simplify is optional/gated (see run_cluster_ora).
# ==============================================================================

#' Read a two-column term table, tolerant of a missing header row
#'
#' The offline TERM2GENE/TERM2NAME tables are documented as header-carrying, and
#' read.delim(header = TRUE) is correct for those. A genuinely headerless file with
#' >= 2 rows would otherwise silently lose its first mapping (row 1 consumed as
#' column names, and a >0-row result never triggers a fallback). read.delim mangles
#' an ID-shaped header cell (e.g. "GO:0008150" -> "GO.0008150", "00010" -> "X00010")
#' while a real text header ("term", "GO_ID") stays a word, so an ID-shaped first
#' column name signals a headerless file: re-read with header = FALSE to keep all rows.
#'
#' @param path Path to a tab-separated table.
#' @return A data.frame with every data row preserved.
.read_term_table <- function(path) {
    x <- read.delim(path, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE, row.names = NULL)
    first_name <- if (ncol(x) >= 1) colnames(x)[1] else ""
    looks_like_id <- grepl("^(GO\\.[0-9]|X[0-9]{4,}|[A-Za-z]{2,4}[0-9]{4,})", first_name)
    if (isTRUE(looks_like_id) || nrow(x) == 0) {
        x <- read.delim(path, sep = "\t", header = FALSE,
                        stringsAsFactors = FALSE, row.names = NULL)
    }
    x
}

load_local_pathway_tables <- function(annotation_dir,
                                      databases = c("KEGG", "GO_BP", "GO_MF", "GO_CC"),
                                      feature_ids = NULL) {

    if (!dir.exists(annotation_dir)) {
        stop("Annotation directory not found: ", annotation_dir)
    }

    # Map database names to file pairs
    file_map <- list(
        KEGG  = list(gene = "KEGG_pathway2gene.tab", name = "KEGG_pathway2name.tab"),
        GO_BP = list(gene = "GO2gene_BP.tab",        name = "GO2name_BP.tab"),
        GO_MF = list(gene = "GO2gene_MF.tab",        name = "GO2name_MF.tab"),
        GO_CC = list(gene = "GO2gene_CC.tab",        name = "GO2name_CC.tab")
    )

    result <- list()

    for (db in databases) {
        if (!db %in% names(file_map)) {
            warning("Unknown database '", db, "' — skipping. ",
                    "Valid: ", paste(names(file_map), collapse = ", "))
            next
        }

        fmap <- file_map[[db]]
        gene_file <- file.path(annotation_dir, fmap$gene)
        name_file <- file.path(annotation_dir, fmap$name)

        if (!file.exists(gene_file)) {
            warning("TERM2GENE file not found for ", db, ": ", gene_file, " — skipping")
            next
        }
        if (!file.exists(name_file)) {
            warning("TERM2NAME file not found for ", db, ": ", name_file, " — skipping")
            next
        }

        # Read two-column tab files. Header-carrying is the documented format;
        # .read_term_table() also tolerates a genuinely headerless file (which
        # would otherwise silently drop its first mapping row).
        term2gene <- .read_term_table(gene_file)
        term2name <- .read_term_table(name_file)

        # Validate column count
        if (ncol(term2gene) < 2) {
            warning(db, " TERM2GENE file has fewer than 2 columns: ", gene_file, " — skipping")
            next
        }
        if (ncol(term2name) < 2) {
            warning(db, " TERM2NAME file has fewer than 2 columns: ", name_file, " — skipping")
            next
        }

        # Keep only first two columns, standardize names for clusterProfiler
        term2gene <- term2gene[, 1:2, drop = FALSE]
        term2name <- term2name[, 1:2, drop = FALSE]
        colnames(term2gene) <- c("term", "gene")
        colnames(term2name) <- c("term", "name")

        # Remove NAs and empty values
        term2gene <- term2gene[!is.na(term2gene$term) & !is.na(term2gene$gene) &
                               nzchar(term2gene$term) & nzchar(term2gene$gene), , drop = FALSE]
        term2name <- term2name[!is.na(term2name$term) & !is.na(term2name$name) &
                               nzchar(term2name$term) & nzchar(term2name$name), , drop = FALSE]

        if (nrow(term2gene) == 0) {
            warning(db, " TERM2GENE table is empty after cleaning: ", gene_file, " — skipping")
            next
        }

        # Overlap check against pipeline feature IDs
        if (!is.null(feature_ids) && length(feature_ids) > 0) {
            genes_in_db <- unique(term2gene$gene)
            overlap <- length(intersect(genes_in_db, feature_ids))
            pct <- round(100 * overlap / length(feature_ids), 1)
            message("  ", db, ": ", nrow(term2gene), " gene-term pairs, ",
                    length(unique(term2gene$term)), " terms, ",
                    overlap, "/", length(feature_ids), " features overlap (", pct, "%)")
            if (pct < 5) {
                warning(db, ": very low overlap (", pct, "%) between TERM2GENE genes and ",
                        "pipeline feature IDs. Check that gene ID types match.")
            }
        } else {
            message("  ", db, ": ", nrow(term2gene), " gene-term pairs, ",
                    length(unique(term2gene$term)), " terms")
        }

        result[[db]] <- list(TERM2GENE = term2gene, TERM2NAME = term2name)
    }

    if (length(result) == 0) {
        warning("No local pathway tables loaded from: ", annotation_dir)
    }

    result
}

# ==============================================================================
# RANKED GENE LIST BUILDERS (Phase 1 — enrichment migration)
# ==============================================================================
# These functions build named numeric vectors suitable for clusterProfiler::GSEA().
# Each implements a specific ranking strategy from the legacy enrichment workflow.
# Input: per-contrast DE tables from the current pipeline (FeatureID, log2FoldChange, pvalue).

#' Build ranked gene list: -log10(pvalue), no direction
#'
#' Ranking value is always positive. Genes with the most significant p-values rank highest.
#'
#' @param de_table Data.frame with FeatureID and pvalue columns
#' @return Named numeric vector (gene IDs as names, ranking values as elements), sorted descending
rank_by_pval_wo_direction <- function(de_table) {
    df <- data.frame(
        gene = de_table$FeatureID,
        pval = de_table$pvalue,
        stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$pval), , drop = FALSE]
    # Floor p at the smallest positive double so an underflowed p == 0 yields a
    # large FINITE score and stays at the TOP of the ranking. (The legacy
    # Inf -> 0 replacement pushed the strongest signal to the bottom — see PR #128.)
    df$rank_val <- -log10(pmax(df$pval, .Machine$double.xmin))
    df <- df[order(df$rank_val, decreasing = TRUE), , drop = FALSE]

    ranks <- df$rank_val
    names(ranks) <- df$gene
    ranks
}

#' Build ranked gene list: sign(FC) * -log10(pvalue)
#'
#' Signed ranking: positive values = upregulated and significant,
#' negative values = downregulated and significant.
#'
#' @param de_table Data.frame with FeatureID, log2FoldChange, and pvalue columns
#' @return Named numeric vector sorted descending
rank_by_pval_with_direction <- function(de_table) {
    df <- data.frame(
        gene = de_table$FeatureID,
        lfc  = de_table$log2FoldChange,
        pval = de_table$pvalue,
        stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$lfc) & !is.na(df$pval), , drop = FALSE]
    # Floor p at the smallest positive double so an underflowed p == 0 yields a
    # large FINITE magnitude, kept at the significant end after signing rather
    # than collapsing to 0 as the legacy Inf -> 0 replacement did — see PR #128.
    df$neg_log_p <- -log10(pmax(df$pval, .Machine$double.xmin))
    # Apply direction: positive LFC -> positive rank, negative LFC -> negative rank
    # Legacy logic: if fc is NA, rank = 0; if fc > 0, rank = pval; else rank = -pval
    # Since we already filtered NA lfc, just apply sign
    df$rank_val <- ifelse(df$lfc > 0, df$neg_log_p, -df$neg_log_p)
    # Note: lfc == 0 → treated as downregulated (-neg_log_p), matching legacy behavior
    df <- df[order(df$rank_val, decreasing = TRUE), , drop = FALSE]

    ranks <- df$rank_val
    names(ranks) <- df$gene
    ranks
}

#' Build ranked gene list: log2 of signed fold change
#'
#' Legacy behavior: converts linear FC via ifelse(fc > 0, fc, -1/fc), then log2,
#' then signif(digits=4). The current pipeline provides log2FoldChange, so we
#' recover linear FC first: linearFC = 2^log2FoldChange.
#'
#' @param de_table Data.frame with FeatureID, log2FoldChange, and pvalue columns
#' @return Named numeric vector sorted descending
rank_by_fc <- function(de_table) {
    df <- data.frame(
        gene = de_table$FeatureID,
        lfc  = de_table$log2FoldChange,
        pval = de_table$pvalue,
        stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$lfc) & !is.na(df$pval), , drop = FALSE]

    # Recover linear fold change from log2FC
    # Legacy input is linearFC where: up = positive value > 1, down = negative value (e.g. -2)
    # Convention: linearFC > 0 means upregulated, linearFC < 0 means downregulated
    # From log2FC: if lfc >= 0, linearFC = 2^lfc (positive, > 1)
    #              if lfc < 0,  linearFC = -(2^(-lfc)) = -(1/(2^lfc)) to get negative linear FC
    # This matches the legacy convention where downregulated genes have negative linearFC.
    linear_fc <- ifelse(df$lfc >= 0, 2^df$lfc, -(2^(-df$lfc)))

    # Legacy transform: ifelse(fc > 0, fc, -1/fc) then log2
    # Maps signed linear FC into a symmetric log2 scale:
    #   fc = +2 → log2(2) = 1;  fc = -2 → log2(-1/(-2)) = log2(0.5) = -1
    # Net effect is equivalent to log2FC, but we apply it to match legacy exactly.
    fc_transformed <- ifelse(linear_fc > 0, linear_fc, -1 / linear_fc)
    df$rank_val <- log2(fc_transformed)
    # Legacy: signif(digits = 4)
    df$rank_val <- signif(df$rank_val, digits = 4)

    df <- df[order(df$rank_val, decreasing = TRUE), , drop = FALSE]

    ranks <- df$rank_val
    names(ranks) <- df$gene
    ranks
}

#' Build ranked gene list: minimum p-value across all contrasts
#'
#' For each gene, takes the minimum raw p-value across all provided contrasts,
#' then ranks by -log10(min_pvalue). Always positive (no direction).
#'
#' @param de_tables Named list of DE tables (each with FeatureID and pvalue columns)
#' @return Named numeric vector sorted descending
rank_by_min_pval_any_contrast <- function(de_tables) {
    if (length(de_tables) == 0) return(numeric(0))

    # Collect pvalue columns from all contrasts, aligned by gene ID
    pval_list <- lapply(de_tables, function(dt) {
        df <- data.frame(
            gene = dt$FeatureID,
            pval = dt$pvalue,
            stringsAsFactors = FALSE
        )
        df[!duplicated(df$gene), , drop = FALSE]
    })

    # Merge all contrasts by gene, keeping all genes (full outer join)
    merged <- pval_list[[1]]
    colnames(merged)[2] <- "pval_1"
    if (length(pval_list) > 1) {
        for (i in 2:length(pval_list)) {
            p <- pval_list[[i]]
            colnames(p)[2] <- paste0("pval_", i)
            merged <- merge(merged, p, by = "gene", all = TRUE)
        }
    }

    # Compute row-wise minimum pvalue (matching legacy: min(x, na.rm = TRUE), NA if all NA)
    pval_cols <- grep("^pval_", colnames(merged), value = TRUE)
    if (length(pval_cols) == 1) {
        min_pval <- merged[[pval_cols]]
    } else {
        pval_mat <- as.matrix(merged[, pval_cols, drop = FALSE])
        min_pval <- apply(pval_mat, 1, function(x) {
            if (all(is.na(x))) NA else min(x, na.rm = TRUE)
        })
    }

    df <- data.frame(
        gene = merged$gene,
        min_pval = min_pval,
        stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$min_pval), , drop = FALSE]
    # Floor p at the smallest positive double (see PR #128): an underflowed
    # min p == 0 stays at the TOP instead of collapsing to 0 via Inf -> 0.
    df$rank_val <- -log10(pmax(df$min_pval, .Machine$double.xmin))
    df <- df[order(df$rank_val, decreasing = TRUE), , drop = FALSE]

    ranks <- df$rank_val
    names(ranks) <- df$gene
    ranks
}

#' Build all ranked gene lists for GSEA
#'
#' Convenience wrapper that builds all four ranking variants from per-contrast DE tables.
#'
#' @param de_tables Named list of DE tables (each with FeatureID, log2FoldChange, pvalue)
#' @return Nested list: ranking_method -> contrast_name -> named numeric vector.
#'   The "any_contrast" method has a single entry keyed "any_contrast".
build_ranked_gene_lists <- function(de_tables) {
    ranked <- list(
        pval_wo_direction   = list(),
        pval_with_direction = list(),
        fc                  = list()
    )

    for (contrast in names(de_tables)) {
        dt <- de_tables[[contrast]]
        ranked[["pval_wo_direction"]][[contrast]]   <- rank_by_pval_wo_direction(dt)
        ranked[["pval_with_direction"]][[contrast]] <- rank_by_pval_with_direction(dt)
        ranked[["fc"]][[contrast]]                  <- rank_by_fc(dt)
    }

    # Cross-contrast: minimum pvalue across all contrasts
    ranked[["pval_wo_direction"]][["any_contrast"]] <- rank_by_min_pval_any_contrast(de_tables)

    ranked
}

# ==============================================================================
# LOCAL GSEA (Phase 1 — enrichment migration)
# ==============================================================================

#' Run GSEA using local TERM2GENE/TERM2NAME tables
#'
#' Wraps clusterProfiler::GSEA() with the legacy enrichment parameters.
#' maxGSSize is set to the total number of unique genes in TERM2GENE (legacy behavior).
#'
#' @param ranked_genes Named numeric vector (gene IDs as names, ranking metric as values).
#'   Must be sorted descending.
#' @param term2gene Two-column data.frame: term ID, gene ID
#' @param term2name Two-column data.frame: term ID, term name
#' @param pvalueCutoff Adjusted p-value cutoff (default 0.05)
#' @param pAdjustMethod P-value adjustment method (default "fdr")
#' @return gseaResult object, or NULL if GSEA fails or produces no results
run_gsea_local <- function(ranked_genes,
                           term2gene,
                           term2name,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "fdr") {

    if (length(ranked_genes) == 0) {
        message("    Empty ranked gene list — skipping GSEA")
        return(NULL)
    }

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        warning("clusterProfiler package required for GSEA. ",
                "Install with: BiocManager::install('clusterProfiler')")
        return(NULL)
    }

    # Legacy behavior: maxGSSize = total unique genes in the pathway database
    nr_total_genes <- length(unique(term2gene[, 2]))

    res <- tryCatch({
        clusterProfiler::GSEA(
            geneList      = ranked_genes,
            TERM2GENE     = term2gene,
            TERM2NAME     = term2name,
            minGSSize     = 4,
            maxGSSize     = nr_total_genes,
            pAdjustMethod = pAdjustMethod,
            pvalueCutoff  = pvalueCutoff,
            # Run fgsea serially: our parallelism is the outer future.apply job
            # fan-out, so fgsea's default (Windows) SnowParam SOCK cluster would
            # be nested parallelism — fragile ("error writing to connection")
            # and redundant. SerialParam keeps each GSEA job single-process.
            BPPARAM       = BiocParallel::SerialParam()
        )
    }, error = function(e) {
        message("    GSEA failed: ", e$message)
        NULL
    })

    res
}

#' Run independent enrichment jobs, in parallel when workers > 1
#'
#' Generic, method-agnostic orchestration layer for the enrichment engine. A
#' "job" is any independent unit of enrichment work (e.g. one ORA over one gene
#' set and one database, or one GSEA over one ranked list and one database).
#' The same mechanism serves GO, KEGG, ORA, GSEA, and any future method
#' (Reactome, WikiPathways, MSigDB, ...) — callers just supply their own flat
#' job list and a pure worker function; the parallel logic lives here, once.
#'
#' Workers must be PURE compute: no file writing, no plotting, no messages
#' (those belong in the caller's serial assembly step). Results are returned in
#' the SAME order as `jobs`, so downstream assembly is deterministic regardless
#' of worker count or scheduling.
#'
#' @param jobs List of opaque job descriptors.
#' @param fun  Worker function mapping one job -> one result (pure compute).
#' @param workers Integer. Controls only the backend, never the results:
#'   `<= 1` uses a `future::sequential` plan (one job at a time, in-process);
#'   `> 1` uses `future::multisession` with that many workers (Windows-safe,
#'   separate processes). Both paths go through `future.apply::future_lapply()`.
#' @param seed Integer reproducibility seed (the project's `params$seed`). Passed
#'   as `future.seed = seed`, which derives each job's independent L'Ecuyer-CMRG
#'   RNG stream **from this fixed integer** — so results depend only on `seed` and
#'   job position, never on the ambient RNG state, the backend, or the worker
#'   count. This is the single source of truth for enrichment reproducibility:
#'   permutation-based methods (GSEA) return identical results for `workers = 1`,
#'   `4`, or any N, and across independent pipeline rebuilds. (Using
#'   `future.seed = TRUE` instead would key off the ambient RNG at call time, which
#'   drifts between builds — the bug this replaces.) Method-agnostic: any future
#'   RNG-using enrichment method inherits reproducibility for free. If
#'   `future`/`future.apply` are not installed, it degrades to plain `lapply()`.
#' @return List of results, one per job, in input order.
run_enrichment_jobs <- function(jobs, fun, workers = 1L, seed = 1L) {
    if (length(jobs) == 0) return(list())

    have_future <- requireNamespace("future", quietly = TRUE) &&
        requireNamespace("future.apply", quietly = TRUE)

    if (!have_future) {
        if (workers > 1) {
            message("  future/future.apply not available — running sequentially. ",
                    "Install with: renv::install(c('future', 'future.apply'))")
        }
        # Reproducibility in the no-future fallback: seed the sequential run from
        # the same project `seed` the future path uses. withr::with_seed sets the
        # RNG deterministically and RESTORES the caller's global RNG state on exit
        # (no leak). Streams are not byte-identical to future's L'Ecuyer-CMRG
        # per-job streams, but results are reproducible across independent runs
        # with the same seed — which is what matters for permutation-based GSEA.
        return(withr::with_seed(seed, lapply(jobs, fun)))
    }

    # NB: workers must capture ONLY the data they use. Callers build worker
    # functions with a minimal environment (see .make_ora_worker()), which keeps
    # exported globals tiny (~5 MiB here) — well under future's default 500 MiB
    # guard. That guard is intentionally left at its default: it is a useful
    # early warning if a future method ever starts broadcasting large objects.
    # Route EVERY worker count through future_lapply with an EXPLICIT integer
    # future.seed (see @param seed): RNG streams depend only on `seed` + job
    # position — not on ambient RNG, backend, or worker count — so results are
    # worker-count-invariant AND identical across independent rebuilds. Sequential
    # plan for workers <= 1 keeps one-job-at-a-time, in-process semantics.
    old_plan <- if (workers > 1) {
        future::plan(future::multisession, workers = workers)
    } else {
        future::plan(future::sequential)
    }
    on.exit(future::plan(old_plan), add = TRUE)

    future.apply::future_lapply(jobs, fun, future.seed = seed)
}

# ==============================================================================
# ENRICHMENT OUTPUT LAYOUT — centralized path builders
# ==============================================================================
# All enrichment output paths are constructed HERE (one place per subsystem) so
# the directory layout is defined once, not scattered across file.path() calls.
# The analysis context lives in the DIRECTORY path; filenames within a unit are
# short and fixed (results.csv, dotplot.pdf/png, ...) — no db/method/contrast/
# collection tokens repeated in the filename. GO_BP/GO_MF/GO_CC/KEGG all use the
# identical structure (db is just a directory level).

#' Directory for one GSEA result unit: <db>/ranking_by_<method>/<contrast>
#'
#' @param gsea_root The GSEA root (e.g. Enrichment/GSEA).
#' @param db_name Database (GO_BP/GO_MF/GO_CC/KEGG).
#' @param ranking_method Ranking method (fc / pval_with_direction / pval_wo_direction).
#' @param contrast Contrast (or the `any_contrast` pseudo-contrast).
#' @return The unit directory path (not created).
#' @noRd
gsea_unit_dir <- function(gsea_root, db_name, ranking_method, contrast) {
    file.path(gsea_root, db_name,
              paste0("ranking_by_", ranking_method), contrast)
}

#' Directory for one ORA result unit under Enrichment/ORA/<db>/
#'
#' Maps the gene-list collection + round to the nested layout:
#'   all_DE               -> <db>/all_DE/any_contrast
#'   contrasts            -> <db>/contrasts/with_direction/<round>
#'   contrasts_wo_direction -> <db>/contrasts/without_direction/<round>
#'   partition            -> <db>/clustering/partition
#'   binary_patterns      -> <db>/clustering/binary_patterns
#' Unknown collections fall back to <db>/<clust_method>/<clust_round> so nothing
#' is ever written outside ORA/<db>/.
#'
#' @param ora_root The ORA root (e.g. Enrichment/ORA).
#' @param db_name Database (GO_BP/GO_MF/GO_CC/KEGG).
#' @param clust_method Collection: all_DE / contrasts / contrasts_wo_direction /
#'   partition / binary_patterns.
#' @param clust_round Round within the collection (contrast name, "any_contrast",
#'   "k", or "best").
#' @return The unit directory path (not created).
#' @noRd
ora_unit_dir <- function(ora_root, db_name, clust_method, clust_round) {
    db_dir <- file.path(ora_root, db_name)
    switch(clust_method,
        all_DE                 = file.path(db_dir, "all_DE", clust_round),
        contrasts              = file.path(db_dir, "contrasts", "with_direction", clust_round),
        contrasts_wo_direction = file.path(db_dir, "contrasts", "without_direction", clust_round),
        partition              = file.path(db_dir, "clustering", "partition"),
        binary_patterns        = file.path(db_dir, "clustering", "binary_patterns"),
        # fallback: keep any unforeseen collection inside ORA/<db>/
        file.path(db_dir, clust_method, clust_round)
    )
}

#' Build a pure-compute GSEA worker with a minimal captured environment
#'
#' Returns a `function(job)` that runs GSEA for one job. Defining it here (not
#' nested inside run_gsea_all) bounds the closure's environment to just the
#' arguments below — so future.apply serializes only `local_tables` + scalars,
#' never the `run_gsea_all` frame (which holds the large `jobs` list of per-job
#' ranked vectors and `ranked_genes`). The per-job ranked vector arrives in
#' `job$ranked` (built in run_gsea_all); the worker does pure computation only
#' (no file I/O, no messages), with fgsea forced serial via SerialParam.
#' Analogous to .make_ora_worker().
#'
#' @param local_tables Output of load_local_pathway_tables().
#' @param pvalueCutoff GSEA adjusted-p cutoff.
#' @param pAdjustMethod P-value adjustment method.
#' @return A function(job) -> list(ranking_method, contrast, db_name, gsea_result).
#' @noRd
.make_gsea_worker <- function(local_tables, pvalueCutoff, pAdjustMethod) {
    force(local_tables); force(pvalueCutoff); force(pAdjustMethod)
    function(job) {
        ranked    <- job$ranked
        term2gene <- local_tables[[job$db_name]]$TERM2GENE
        term2name <- local_tables[[job$db_name]]$TERM2NAME

        res <- tryCatch({
            clusterProfiler::GSEA(
                geneList      = ranked,
                TERM2GENE     = term2gene,
                TERM2NAME     = term2name,
                minGSSize     = 4,
                maxGSSize     = length(unique(term2gene[, 2])),
                pAdjustMethod = pAdjustMethod,
                pvalueCutoff  = pvalueCutoff,
                # Serial fgsea inside each job. Outer parallelism is future.apply
                # over jobs (see run_enrichment_jobs); fgsea's default (Windows)
                # SnowParam spawns a 14-worker SOCK cluster PER job, which is
                # nested parallelism — the source of "error writing to
                # connection" failures under RStudio — and redundant here.
                BPPARAM       = BiocParallel::SerialParam()
            )
        }, error = function(e) {
            structure(list(message = e$message), class = "gsea_error")
        })
        list(
            ranking_method = job$ranking_method,
            contrast       = job$contrast,
            db_name        = job$db_name,
            gsea_result    = res
        )
    }
}

#' Run GSEA across all ranking methods, contrasts, and databases
#'
#' Orchestrator that runs run_gsea_local() for every combination of
#' ranking method x contrast x database loaded from local tables.
#' Jobs are independent and run in parallel when future.apply is available.
#'
#' @param ranked_genes Output of build_ranked_gene_lists()
#' @param local_tables Output of load_local_pathway_tables()
#' @param pvalueCutoff Adjusted p-value cutoff
#' @param pAdjustMethod P-value adjustment method
#' @param output_dir Directory for GSEA result CSVs
#' @param workers Number of parallel workers (default 1 = sequential).
#'   Parallelization uses future::plan(multisession) which is Windows-safe.
#' @return Nested list compatible with downstream consumers:
#'   contrast -> db_method_key -> data.frame with padj, NES, etc.
run_gsea_all <- function(ranked_genes,
                         local_tables,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "fdr",
                         output_dir = NULL,
                         workers = 1,
                         per_pathway_artifacts = FALSE,
                         max_terms_in_dotplot = 20,
                         dotplot = TRUE,
                         ridgeplot = FALSE,
                         ridgeplot_all_genes = FALSE,
                         pathway_heatmaps = FALSE,
                         gene_context = NULL,
                         expr_mat = NULL,
                         annotation_col = NULL,
                         seed = 1L) {

    if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # ------------------------------------------------------------------
    # 1. Build flat job list. Each job carries its OWN ranked vector
    #    (`ranked`), not a reference into `ranked_genes`. This keeps the whole
    #    `ranked_genes` structure OUT of the worker closure's environment, so it
    #    is never broadcast as a future global — the per-job vectors ride in the
    #    `jobs` iteration list (sent one-at-a-time to workers, held once on the
    #    master), keeping exported globals ~= local_tables regardless of the
    #    number of contrasts. (`local_tables` stays a captured global: small.)
    # ------------------------------------------------------------------
    jobs <- list()
    for (ranking_method in names(ranked_genes)) {
        for (contrast in names(ranked_genes[[ranking_method]])) {
            if (length(ranked_genes[[ranking_method]][[contrast]]) == 0) next

            for (db_name in names(local_tables)) {
                jobs[[length(jobs) + 1]] <- list(
                    ranking_method = ranking_method,
                    contrast       = contrast,
                    db_name        = db_name,
                    ranked         = ranked_genes[[ranking_method]][[contrast]]
                )
            }
        }
    }

    if (length(jobs) == 0) {
        return(list(results = list(), plot_files = list()))
    }

    message("  ", length(jobs), " GSEA jobs to run",
            if (workers > 1) paste0(" (", workers, " workers)") else " (sequential)")

    # ------------------------------------------------------------------
    # 2. Run GSEA computation (parallel or sequential)
    # ------------------------------------------------------------------
    # Build the worker via a factory (NOT a nested closure). A closure defined
    # inside run_gsea_all() would carry this whole execution frame — including
    # the large `jobs` list (which holds every per-job ranked vector) and
    # `ranked_genes` — as its environment, and future would serialize all of it
    # (the 500 MiB-globals failure). The factory bounds the worker's environment
    # to only `local_tables` + scalars; the ranked vectors ride in `jobs` (the
    # future_lapply iteration list, sent per-job, not a broadcast global).
    run_one_gsea_job <- .make_gsea_worker(local_tables, pvalueCutoff, pAdjustMethod)

    # Dispatch pure GSEA compute through the generic parallel orchestration
    # layer. Assembly + all file I/O happen serially below (deterministic).
    # `seed` (the project's params$seed) fixes the per-job RNG streams so GSEA
    # is reproducible across worker counts and independent rebuilds.
    job_results <- run_enrichment_jobs(jobs, run_one_gsea_job, workers, seed = seed)

    # ------------------------------------------------------------------
    # 3. Assemble results and write files (serial, deterministic)
    # ------------------------------------------------------------------
    results <- list()
    plot_files <- list()

    for (jr in job_results) {
        db_name        <- jr$db_name
        ranking_method <- jr$ranking_method
        contrast       <- jr$contrast
        res            <- jr$gsea_result
        result_key     <- paste0(db_name, "_gsea_", ranking_method)

        # Handle failed jobs
        if (inherits(res, "gsea_error")) {
            message("  GSEA failed: ", db_name, " | ", ranking_method, " | ",
                    contrast, " — ", res$message)
            next
        }

        if (is.null(res) || nrow(as.data.frame(res)) == 0) {
            message("  ", db_name, " | ", ranking_method, " | ", contrast,
                    ": no results returned")
            next
        }

        # Convert to data.frame for storage and downstream compatibility
        res_df <- as.data.frame(res)
        res_df$contrast <- contrast
        res_df$database <- db_name
        res_df$ranking_method <- ranking_method

        # Ensure downstream-required columns exist
        if ("p.adjust" %in% colnames(res_df) && !"padj" %in% colnames(res_df)) {
            res_df$padj <- res_df$p.adjust
        }
        if ("Description" %in% colnames(res_df) && !"pathway" %in% colnames(res_df)) {
            res_df$pathway <- res_df$Description
        }

        # Store in nested structure
        if (is.null(results[[contrast]])) results[[contrast]] <- list()
        results[[contrast]][[result_key]] <- res_df

        n_sig <- nrow(.gsea_significant_rows(res_df, pvalueCutoff))
        message("  ", db_name, " | ", ranking_method, " | ", contrast,
                ": ", n_sig, " significant (padj <= ", pvalueCutoff, ")")

        # Write CSV
        if (!is.null(output_dir)) {
            # Unit dir: GSEA/<db>/ranking_by_<method>/<contrast>/. Filenames are
            # short and fixed (context is already in the path).
            gsea_sub_dir <- gsea_unit_dir(output_dir, db_name, ranking_method, contrast)
            dir.create(gsea_sub_dir, recursive = TRUE, showWarnings = FALSE)
            csv_file <- file.path(gsea_sub_dir, "results.csv")
            write.csv(res_df, file = csv_file, row.names = FALSE)

            # Generate dotplot if significant results exist (config-gated).
            # Primary: enrichplot::dotplot() on the gseaResult object.
            # Fallback: basic ggplot2 scatterplot (visual approximation only).
            if (isTRUE(dotplot) && n_sig >= 3) {
                plot_file <- file.path(gsea_sub_dir, "dotplot.png")
                plot_key <- paste0(db_name, "_", ranking_method, "_", contrast)
                show_n <- min(max_terms_in_dotplot, n_sig)

                plotted <- FALSE

                # Primary: enrichplot::dotplot on gseaResult S4 object
                if (requireNamespace("enrichplot", quietly = TRUE)) {
                    tryCatch({
                        p <- enrichplot::dotplot(res, showCategory = show_n)
                        ggplot2::ggsave(plot_file, p, width = 10, height = 8)
                        plot_files[[plot_key]] <- plot_file
                        plotted <- TRUE
                    }, error = function(e) {
                        message("    enrichplot::dotplot() failed: ", e$message,
                                " — falling back to ggplot2")
                    })
                }

                # Fallback: basic ggplot2 scatterplot (not equivalent to dotplot)
                if (!plotted) {
                    tryCatch({
                        top <- head(res_df[order(res_df$padj), ], show_n)
                        top$pathway_label <- substr(top$pathway, 1, 60)
                        p <- ggplot2::ggplot(
                            top,
                            ggplot2::aes(x = NES, y = reorder(pathway_label, NES))
                        ) +
                            ggplot2::geom_point(
                                ggplot2::aes(size = setSize, color = -log10(padj))
                            ) +
                            ggplot2::scale_color_gradient(
                                low = "blue", high = "red", name = "-log10(padj)"
                            ) +
                            ggplot2::labs(
                                title = paste("GSEA:", db_name, "|", ranking_method),
                                subtitle = contrast,
                                x = "Normalized Enrichment Score",
                                y = "",
                                size = "Gene Set Size"
                            ) +
                            ggplot2::theme_minimal() +
                            ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
                        ggplot2::ggsave(plot_file, p, width = 10, height = 8)
                        plot_files[[plot_key]] <- plot_file
                    }, error = function(e) {
                        message("    Fallback plot also failed: ", e$message)
                    })
                }
            }

            # ----------------------------------------------------------
            # Per-pathway GSEA artifacts (gseaplot2 PNGs + rich core-gene CSVs +
            # optional expression heatmaps). The per-pathway loop runs when EITHER
            # gsea.per_pathway_artifacts (plots + core-gene tables) OR
            # plots.pathway_heatmaps is enabled; each artifact is gated inside.
            # `gene_context`/`expr_mat`/`annotation_col` are used only here in the
            # serial stage — never captured by the bounded GSEA worker.
            # ----------------------------------------------------------
            if (isTRUE(per_pathway_artifacts) || isTRUE(pathway_heatmaps)) {
                save_gsea_per_pathway_artifacts(
                    gsea_result    = res,
                    res_df         = res_df,
                    output_dir     = file.path(gsea_sub_dir, "per_pathway"),
                    gene_context   = gene_context,
                    expr_mat       = expr_mat,
                    annotation_col = annotation_col,
                    plots          = isTRUE(per_pathway_artifacts),
                    heatmaps       = isTRUE(pathway_heatmaps),
                    sig_cutoff     = pvalueCutoff
                )
            }

            # GSEA ridgeplots (legacy ridgeplot_edited / ridgeplot_edited1).
            # Additive, config-gated, fail-soft; serial stage only (never in a
            # worker). Leading-edge -> ridgeplot/plot.png + data.csv; all-genes ->
            # ridgeplot/plot_all_genes.png + data_all_genes.csv (same dir).
            if (isTRUE(ridgeplot) || isTRUE(ridgeplot_all_genes)) {
                ridge_key <- paste0(db_name, "_", ranking_method, "_", contrast)
                ridge_dir <- file.path(gsea_sub_dir, "ridgeplot")
                if (isTRUE(ridgeplot)) {
                    rp <- plot_gsea_ridgeplot(
                        gsea_result     = res,
                        out_dir         = ridge_dir,
                        show_category   = max_terms_in_dotplot,
                        x_axis_title    = paste0("ranking (", ranking_method, ")"),
                        core_enrichment = TRUE
                    )
                    if (!is.null(rp)) plot_files[[paste0("ridge_", ridge_key)]] <- rp
                }
                if (isTRUE(ridgeplot_all_genes)) {
                    rpa <- plot_gsea_ridgeplot(
                        gsea_result     = res,
                        out_dir         = ridge_dir,
                        show_category   = max_terms_in_dotplot,
                        x_axis_title    = paste0("ranking (", ranking_method, ")"),
                        core_enrichment = FALSE
                    )
                    if (!is.null(rpa)) plot_files[[paste0("ridge_allG_", ridge_key)]] <- rpa
                }
            }
        }
    }

    list(results = results, plot_files = plot_files)
}

# ==============================================================================
# PER-PATHWAY GSEA ARTIFACTS (Phase 3 — legacy outputs)
# ==============================================================================

#' Select the significant rows of a GSEA result at the effective cutoff
#'
#' clusterProfiler::GSEA() already filters its `@result` to rows with both
#' `pvalue <= pvalueCutoff` AND `p.adjust <= pvalueCutoff`, so a valid GSEA
#' data.frame's rows are all significant at the configured cutoff. This helper
#' re-selects on the ADJUSTED p-value (`padj`) using the SAME configured cutoff
#' (`gsea_pvalue_cutoff` if set, else `pvalue_cutoff`) so every downstream
#' artifact decision (dotplot gate, per-pathway plots/tables/heatmaps, counts)
#' tracks the config instead of a hard-coded 0.05. Uses `<=` to match GSEA's
#' filtering contract. NA `padj` rows are excluded.
#'
#' @param res_df GSEA result data.frame (must have a `padj` column).
#' @param sig_cutoff Effective adjusted-p significance cutoff.
#' @return The subset of `res_df` that is significant (0-row frame if none/NA-only).
#' @noRd
.gsea_significant_rows <- function(res_df, sig_cutoff = 0.05) {
    if (is.null(res_df) || !is.data.frame(res_df) ||
        !"padj" %in% names(res_df) || nrow(res_df) == 0) {
        return(res_df[0, , drop = FALSE])
    }
    res_df[!is.na(res_df$padj) & res_df$padj <= sig_cutoff, , drop = FALSE]
}

#' Save per-pathway GSEA artifacts: enrichment plots, core-gene tables, heatmaps
#'
#' For each significant pathway in the GSEA result, produces (each gated):
#'   - plots/{pathway_id}.png            enrichplot::gseaplot2  (needs `plots`)
#'   - core_genes/{pathway_id}.csv       leading-edge genes with stats
#'   - heatmaps_all_genes/{pathway_id}.png   expression heatmap, all pathway
#'                                           genes (needs `heatmaps` + `expr_mat`)
#'   - heatmaps_core_genes/{pathway_id}.png  expression heatmap, leading-edge only
#'
#' The core-gene CSV is enriched by left-joining `gene_context` (a per-gene table
#' of expression / DE stats / z-scores, built omics-specifically upstream and
#' keyed by its first column). When `gene_context` is NULL it falls back to the
#' minimal `gene, rank_value` table. All artifacts are fail-soft; this runs only
#' in the serial output stage, never inside a parallel worker.
#'
#' File writing vs computation: gseaplot2 PNGs AND the core-gene CSVs are written
#' only when `plots` (i.e. gsea.per_pathway_artifacts) is TRUE. Heatmaps are an
#' independent output (`heatmaps` / plots.pathway_heatmaps): they may run with
#' `plots = FALSE`, deriving their gene sets in memory without writing any
#' core-gene CSV.
#'
#' @param gsea_result gseaResult S4 object from clusterProfiler::GSEA()
#' @param res_df Data.frame version of the result (with padj, pathway, etc.)
#' @param output_dir The per_pathway/ directory for this GSEA unit.
#' @param gene_context Optional per-gene data.frame (first column = gene ID)
#'   with annotation / expression / DE / z-score columns to enrich core-gene CSVs.
#' @param expr_mat Optional feature x sample expression matrix for heatmaps.
#' @param annotation_col Optional pheatmap column-annotation data.frame.
#' @param plots Logical: emit gseaplot2 PNGs + core-gene CSVs (default TRUE).
#' @param heatmaps Logical: emit per-pathway expression heatmaps (default FALSE).
#' @param sig_cutoff Effective adjusted-p cutoff selecting which pathways get
#'   artifacts (the configured GSEA cutoff; default 0.05).
#' @param max_heatmap_genes Cap on genes drawn per heatmap (guards huge sets).
#' @noRd
save_gsea_per_pathway_artifacts <- function(gsea_result, res_df, output_dir,
                                            gene_context = NULL, expr_mat = NULL,
                                            annotation_col = NULL,
                                            plots = TRUE, heatmaps = FALSE,
                                            sig_cutoff = 0.05,
                                            max_heatmap_genes = 500) {

    if (is.null(gsea_result) || is.null(output_dir)) return(invisible(NULL))

    sig_rows <- .gsea_significant_rows(res_df, sig_cutoff)
    if (nrow(sig_rows) == 0) return(invisible(NULL))

    # Directories (create only those we will populate). The core_genes CSV
    # directory is created only when `plots` is on — heatmaps never write CSVs.
    plots_dir <- file.path(output_dir, "plots")
    excel_dir <- file.path(output_dir, "core_genes")
    hm_all_dir  <- file.path(output_dir, "heatmaps_all_genes")
    hm_core_dir <- file.path(output_dir, "heatmaps_core_genes")
    if (isTRUE(plots)) {
        dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(excel_dir, recursive = TRUE, showWarnings = FALSE)
    }
    do_heatmaps <- isTRUE(heatmaps) && !is.null(expr_mat)
    if (do_heatmaps) {
        expr_mat <- as.matrix(expr_mat)
        dir.create(hm_all_dir,  recursive = TRUE, showWarnings = FALSE)
        dir.create(hm_core_dir, recursive = TRUE, showWarnings = FALSE)
    }

    has_enrichplot <- requireNamespace("enrichplot", quietly = TRUE)
    ctx_key <- if (!is.null(gene_context) && ncol(gene_context) > 0) {
        names(gene_context)[1]
    } else NULL
    gene_list <- gsea_result@geneList  # full ranked vector (gene -> value)

    # Match the pathway's genes to expression rows on the `Gene:`-normalized key
    # (DE/pathway IDs may have the prefix stripped while expr rownames retain it —
    # see the ID-normalization note in mod_rnaseq_pathway). Select by POSITION so
    # the matrix keeps its original rownames for the plot. Ambiguity guard: keys
    # shared by >1 expr row (e.g. both "Gene:X" and "X") are EXCLUDED from heatmaps
    # so an ambiguous gene is never plotted from an arbitrary row.
    expr_keys <- sub("^Gene:", "", rownames(expr_mat))
    ambig_keys <- unique(expr_keys[duplicated(expr_keys)])
    if (length(ambig_keys) > 0) {
        warning("pathway heatmaps: ", length(ambig_keys), " ambiguous gene ID(s) ",
                "after Gene: normalization (e.g. ",
                paste(utils::head(ambig_keys, 5L), collapse = ", "),
                "); excluded from heatmaps.")
    }

    # Helper: write one expression heatmap for a gene set (fail-soft).
    write_pathway_heatmap <- function(genes, out_png, title) {
        want <- setdiff(sub("^Gene:", "", genes), ambig_keys)  # drop ambiguous
        sel  <- which(expr_keys %in% want)
        if (length(sel) < 2) return(invisible(NULL))  # need >=2 rows to cluster
        if (length(sel) > max_heatmap_genes) sel <- sel[seq_len(max_heatmap_genes)]
        tryCatch({
            hm <- plot_heatmap_core(
                expr_mat[sel, , drop = FALSE],
                annotation_col = annotation_col,
                title          = title,
                scale_rows     = TRUE,
                silent         = TRUE
            )
            save_heatmap_to_file(hm, out_png)
        }, error = function(e) {
            message("      pathway heatmap failed for ", title, ": ", e$message)
        })
    }

    for (i in seq_len(nrow(sig_rows))) {
        pathway_id <- sig_rows$ID[i]
        safe_id    <- gsub("[^a-zA-Z0-9_.-]", "_", pathway_id)  # safe file name
        pw_title   <- sig_rows$Description[i]

        # --- GSEA plot (enrichplot::gseaplot2) ---
        if (isTRUE(plots) && has_enrichplot) {
            tryCatch({
                plot_file <- file.path(plots_dir, paste0(safe_id, ".png"))
                png(plot_file, width = 800, height = 600, res = 120)
                print(enrichplot::gseaplot2(gsea_result, geneSetID = pathway_id,
                                            title = pw_title))
                dev.off()
            }, error = function(e) {
                tryCatch(dev.off(), error = function(e2) NULL)
                message("      gseaplot2 failed for ", pathway_id, ": ", e$message)
            })
        }

        # Leading-edge (core) genes for this pathway.
        core_genes_str <- sig_rows$core_enrichment[i]
        core_genes <- if (is.null(core_genes_str) || is.na(core_genes_str) ||
                          !nzchar(core_genes_str)) {
            character(0)
        } else {
            g <- trimws(strsplit(core_genes_str, "/", fixed = TRUE)[[1]])
            g[nzchar(g)]
        }

        # --- Core-gene CSV (rich when gene_context supplied) ---
        # Written only when per-pathway artifacts are enabled (`plots`); heatmaps
        # alone never write these tables.
        if (isTRUE(plots) && length(core_genes) > 0) {
            tryCatch({
                if (!is.null(ctx_key)) {
                    idx <- match(core_genes, gene_context[[ctx_key]])
                    core_df <- gene_context[idx, , drop = FALSE]
                    core_df$rank_value <- as.numeric(gene_list[core_genes])
                    # Keep genes with no context row but record their id/rank.
                    core_df[[ctx_key]] <- core_genes
                } else {
                    core_df <- data.frame(
                        gene       = core_genes,
                        rank_value = as.numeric(gene_list[core_genes]),
                        stringsAsFactors = FALSE
                    )
                }
                write.csv(core_df, file = file.path(excel_dir, paste0(safe_id, ".csv")),
                          row.names = FALSE)
            }, error = function(e) {
                message("      Core genes CSV failed for ", pathway_id, ": ", e$message)
            })
        }

        # --- Per-pathway expression heatmaps (all genes + core genes) ---
        if (do_heatmaps) {
            all_genes <- tryCatch(as.character(gsea_result@geneSets[[pathway_id]]),
                                  error = function(e) character(0))
            write_pathway_heatmap(all_genes, file.path(hm_all_dir, paste0(safe_id, ".png")),
                                  paste0(pw_title, " (all genes)"))
            if (length(core_genes) > 0) {
                write_pathway_heatmap(core_genes, file.path(hm_core_dir, paste0(safe_id, ".png")),
                                      paste0(pw_title, " (leading edge)"))
            }
        }
    }

    invisible(NULL)
}

#' Plot a GSEA ridgeplot (per-pathway distribution of gene ranking values)
#'
#' Faithful, robust adaptation of the legacy `ridgeplot_edited()`. For the top
#' `show_category` pathways of a GSEA result, draws ridgeline densities of the
#' per-gene ranking values (the leading-edge / core genes' values from the
#' gseaResult's geneList), filled by `-log10(p.adjust)`, ordered by NES. Writes
#' a PNG and a CSV of the underlying data next to it.
#'
#' Pure plotting — no enrichment computation, no RNG (deterministic given the
#' result), and safe to call only in the serial output stage (never in a
#' worker). Fail-soft: on any missing dependency, empty/invalid result, or
#' plotting error it emits a warning and returns NULL without stopping the
#' pipeline.
#'
#' @param gsea_result A clusterProfiler `gseaResult` S4 object.
#' @param out_dir Output directory (the unit's `ridgeplot/` folder). Writes
#'   `plot.png` and `data.csv` here.
#' @param show_category Max pathways to show (reuses `max_terms_in_dotplot`).
#' @param fill Fill column: "p.adjust" (default), "pvalue", or "qvalue".
#' @param x_axis_title X-axis label (e.g. the ranking method).
#' @param core_enrichment When TRUE (default, legacy `ridgeplot_edited`) the
#'   density uses only each pathway's leading-edge/core genes and writes
#'   `plot.png`/`data.csv`. When FALSE (legacy `ridgeplot_edited1(core=FALSE)`)
#'   it uses every ranked gene in the pathway and writes
#'   `plot_all_genes.png`/`data_all_genes.csv`.
#' @return Invisibly the PNG path if written, else NULL.
#' @noRd
plot_gsea_ridgeplot <- function(gsea_result, out_dir, show_category = 20,
                                fill = "p.adjust", x_axis_title = "ranking value",
                                core_enrichment = TRUE) {
    if (!requireNamespace("ggridges", quietly = TRUE)) {
        warning("ggridges not installed; skipping GSEA ridgeplot: ", out_dir)
        return(invisible(NULL))
    }

    valid <- tryCatch(
        !is.null(gsea_result) && methods::is(gsea_result, "gseaResult") &&
            nrow(gsea_result@result) > 0,
        error = function(e) FALSE
    )
    if (!isTRUE(valid)) {
        warning("GSEA ridgeplot skipped (empty/invalid result): ", out_dir)
        return(invisible(NULL))
    }
    # Two legacy variants (ridgeplot_edited1 with core_enrichment T/F): the
    # leading-edge ridgeplot uses only each pathway's core genes; the "all genes"
    # variant uses every ranked gene in the pathway. Distinct short filenames so
    # both live side-by-side in the same ridgeplot/ dir.
    if (isTRUE(core_enrichment)) {
        out_file  <- file.path(out_dir, "plot.png")
        data_file <- file.path(out_dir, "data.csv")
    } else {
        out_file  <- file.path(out_dir, "plot_all_genes.png")
        data_file <- file.path(out_dir, "data_all_genes.csv")
    }

    out <- tryCatch({
        x <- gsea_result
        if (fill == "qvalue") fill <- "qvalues"
        if (!fill %in% colnames(x@result)) fill <- "p.adjust"

        n   <- min(show_category, nrow(x@result))
        rdf <- x@result[seq_len(n), , drop = FALSE]

        # Per-pathway gene sets: leading-edge (core_enrichment column) or, for the
        # all-genes variant, the pathway's full ranked membership (@geneSets).
        if (isTRUE(core_enrichment)) {
            gs2id <- strsplit(as.character(rdf$core_enrichment), "/", fixed = TRUE)
        } else {
            gs2id <- x@geneSets[as.character(rdf$ID)]
        }
        names(gs2id) <- rdf$ID

        genes2values <- x@geneList  # named ranking vector (gene -> value)
        gs2val <- lapply(gs2id, function(g) {
            v <- genes2values[g]; v[!is.na(v)]
        })
        keep   <- vapply(gs2val, length, integer(1)) > 0
        gs2val <- gs2val[keep]
        rdf    <- rdf[keep, , drop = FALSE]
        if (nrow(rdf) < 1) stop("no ranking values for any pathway's core genes")

        len <- vapply(gs2val, length, integer(1))
        ord <- order(rdf$NES, decreasing = FALSE)

        df <- data.frame(
            category = rep(rdf$Description, times = len),
            fillval  = rep(-log10(rdf[[fill]]), times = len),
            value    = unlist(gs2val, use.names = FALSE),
            stringsAsFactors = FALSE
        )
        df$category <- factor(df$category, levels = rdf$Description[ord])

        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        # Underlying-data CSV (legacy parity).
        utils::write.csv(df, file = data_file, row.names = FALSE)

        fill_lab <- paste0("-log10(", fill, ")")
        p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = category, fill = fillval)) +
            ggridges::geom_density_ridges() +
            ggplot2::scale_fill_continuous(low = "blue", high = "red", name = fill_lab) +
            ggplot2::labs(x = x_axis_title, y = NULL) +
            ggplot2::theme_minimal()

        ggplot2::ggsave(out_file, p, width = 8, height = 8, dpi = 150)
        out_file
    }, error = function(e) {
        warning("GSEA ridgeplot failed for ", out_dir, ": ", conditionMessage(e))
        NULL
    })
    invisible(out)
}

# ==============================================================================
# CLUSTER-BASED ORA (Phase 2 — enrichment migration)
# ==============================================================================
# Reproduces the legacy Clusters_Enrichment_Test() behavior:
#   - per-cluster enricher() with minGSSize=0, maxGSSize=10000, qvalueCutoff=1
#   - merge_result() into compareClusterResult
#   - GO simplify (cutoff=0.7, by="p.adjust", select_fun=min)
#   - process_enrichment_table() with fold enrichment
#   - enrichplot::dotplot() on the merged result

#' Run cluster-based ORA for a single database
#'
#' Reproduces legacy Clusters_Enrichment_Test() behavior exactly.
#'
#' @param clusters Named vector: gene IDs as names, cluster labels as values.
#' @param TERM2GENE Two-column data.frame (term ID, gene ID).
#' @param TERM2NAME Two-column data.frame (term ID, term name).
#' @param type "KEGG" or "GO". Controls simplify and @fun slot.
#' @param pvalueCutoff Adjusted p-value cutoff for filtering (default 0.05).
#' @param pAdjustMethod P-value adjustment method (default "fdr").
#' @param outDir ORA unit directory to write into (results.csv / dotplot.pdf /
#'   simplify.csv). NULL skips file writing.
#' @param maxCategory Max categories to show in dotplot (default 1000).
#' @return List of 4 elements (matching legacy):
#'   [[1]] allRes (compareClusterResult or NULL),
#'   [[2]] allRes_simplify (compareClusterResult or NULL, GO only),
#'   [[3]] enrichment_table (data.frame or NULL),
#'   [[4]] enrichment_table_simplify (data.frame or NULL, GO only).
#'   Returns list() if no clusters have significant enrichment.
#' @export
run_cluster_ora <- function(clusters,
                            TERM2GENE,
                            TERM2NAME,
                            type = "KEGG",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "fdr",
                            outDir = NULL,
                            maxCategory = 1000,
                            orgdb = NULL,
                            ont = NULL) {

    # Thin wrapper: pure compute (parallel-safe) + serial file I/O. Kept for
    # backward compatibility and direct callers/tests. The module orchestration
    # calls run_cluster_ora_compute() and write_cluster_ora_outputs() separately
    # so compute can run inside parallel workers while I/O stays serial. GO simplify
    # is always applied for GO results (see run_cluster_ora_compute()).
    res <- run_cluster_ora_compute(
        clusters      = clusters,
        TERM2GENE     = TERM2GENE,
        TERM2NAME     = TERM2NAME,
        type          = type,
        pvalueCutoff  = pvalueCutoff,
        pAdjustMethod = pAdjustMethod,
        orgdb         = orgdb,
        ont           = ont
    )
    if (length(res) == 0) return(list())
    if (!is.null(outDir)) {
        write_cluster_ora_outputs(res, unit_dir = outDir,
                                  type = type, maxCategory = maxCategory)
    }
    res
}

# Per-process cache of GO semantic-similarity data, keyed by ontology. The Wang
# measure (legacy simplify default) derives term similarity from the GO DAG in
# GO.db and needs NO organism OrgDb, so one semData per ontology (BP/MF/CC) is
# reused across all ORA jobs in a process instead of being rebuilt each call
# (godata() is slow). In multisession workers each process keeps its own cache.
.enrich_semdata_cache <- new.env(parent = emptyenv())

#' Build (and cache) GO Wang semantic-similarity data for an ontology — offline
#'
#' Reproduces the legacy `simplify()` dependency: `GOSemSim::godata(ont=...)`
#' over `GO.db` only, with no OrgDb and no network. Returns NULL (never errors)
#' when GO.db / GOSemSim are unavailable so callers can skip-and-warn.
#'
#' @param ont GO ontology: "BP", "MF", or "CC" (NULL -> "BP").
#' @return A `GOSemSimDATA` object, or NULL if it could not be built.
#' @noRd
.go_semdata <- function(ont) {
    key <- as.character(ont %||% "BP")
    cached <- .enrich_semdata_cache[[key]]
    if (!is.null(cached)) return(cached)
    if (!requireNamespace("GOSemSim", quietly = TRUE) ||
        !requireNamespace("GO.db", quietly = TRUE)) {
        return(NULL)
    }
    sem <- tryCatch(
        GOSemSim::godata(ont = key),   # Wang: GO DAG from GO.db, no OrgDb
        error = function(e) {
            message("    GOSemSim::godata(ont=", key, ") failed: ", e$message)
            NULL
        }
    )
    if (!is.null(sem)) .enrich_semdata_cache[[key]] <- sem
    sem
}

#' Cluster ORA — pure computation (no file I/O), parallel-worker-safe
#'
#' The compute half of run_cluster_ora(): per-cluster `enricher()`, merge into a
#' compareClusterResult, fold-enrichment processing, and the optional (gated) GO
#' simplify. Writes no files, draws no plots — safe to run inside future workers.
#' `enricher()` is deterministic (hypergeometric, no RNG), so results are
#' identical regardless of worker count.
#'
#' @param clusters Named vector: gene IDs as names, cluster labels as values.
#' @param TERM2GENE Two-column data.frame (term ID, gene ID).
#' @param TERM2NAME Two-column data.frame (term ID, term name).
#' @param type "KEGG" or "GO". For "GO", redundant terms are ALWAYS simplified
#'   (mandatory; not configurable). For "KEGG" no simplification is applied.
#' @param pvalueCutoff Adjusted p-value cutoff for filtering (default 0.05).
#' @param pAdjustMethod P-value adjustment method (default "fdr").
#' @param orgdb Reserved for future IC-based similarity measures (which require an
#'   organism OrgDb). Not used by the default Wang simplify. (or NULL).
#' @param ont GO ontology ("BP"/"MF"/"CC") for simplify semData (or NULL).
#' @return List of 4 elements (matching legacy):
#'   [[1]] allRes (compareClusterResult or NULL),
#'   [[2]] allRes_simplify (compareClusterResult or NULL, GO only),
#'   [[3]] enrichment_table (data.frame or NULL),
#'   [[4]] enrichment_table_simplify (data.frame or NULL, GO only).
#'   Returns list() if no clusters have significant enrichment.
#' @noRd
run_cluster_ora_compute <- function(clusters,
                                    TERM2GENE,
                                    TERM2NAME,
                                    type = "KEGG",
                                    pvalueCutoff = 0.05,
                                    pAdjustMethod = "fdr",
                                    orgdb = NULL,
                                    ont = NULL) {

    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        warning("clusterProfiler required for cluster ORA. ",
                "Install with: BiocManager::install('clusterProfiler')")
        return(list())
    }

    # ------------------------------------------------------------------
    # Input validation: clusters must be a named atomic vector
    # ------------------------------------------------------------------
    if (!is.atomic(clusters) || is.null(names(clusters))) {
        warning("run_cluster_ora(): clusters must be a named atomic vector ",
                "(gene IDs as names, cluster labels as values). Got: ",
                class(clusters)[1])
        return(list())
    }
    # Remove entries with NA names or NA labels
    valid <- !is.na(names(clusters)) & nzchar(names(clusters)) & !is.na(clusters)
    if (sum(valid) == 0) {
        warning("run_cluster_ora(): no valid gene-cluster assignments after cleaning")
        return(list())
    }
    clusters <- clusters[valid]

    # Per-cluster enrichment
    # Legacy: sort(unique(clusters)). Preserve numeric sort for integer labels,
    # alphabetic sort for character labels.
    allRes0 <- list()
    genes_having_pathway <- unique(TERM2GENE[, 2])
    cluster_labels <- sort(unique(clusters))

    for (cluster_name in cluster_labels) {
        genes_in_cluster <- names(clusters[clusters == cluster_name])
        Genes <- intersect(genes_in_cluster, genes_having_pathway)

        if (length(Genes) == 0) next

        res <- tryCatch({
            clusterProfiler::enricher(
                Genes,
                TERM2GENE     = TERM2GENE,
                TERM2NAME     = TERM2NAME,
                minGSSize     = 0,
                maxGSSize     = 10000,
                pAdjustMethod = pAdjustMethod,
                pvalueCutoff  = pvalueCutoff,
                qvalueCutoff  = 1
            )
        }, error = function(e) {
            message("    enricher() failed for cluster ", cluster_name, ": ", e$message)
            NULL
        })

        if (!is.null(res) &&
            nrow(res@result) > 0 &&
            nrow(res@result[res@result$p.adjust < pvalueCutoff, , drop = FALSE]) > 0) {
            allRes0[[as.character(cluster_name)]] <- res
        }
    }

    if (length(allRes0) == 0) {
        return(list())
    }

    # Merge per-cluster results into a compareClusterResult
    allRes <- clusterProfiler::merge_result(enrichResultList = allRes0)

    # Process the enrichment table (fold enrichment, expanded ratios)
    enrichment_table <- process_enrichment_table(allRes@compareClusterResult)

    # Set @fun slot for enrichplot/simplify dispatch (legacy does this before CSV write)
    # enricher() sets @fun = "enricher" but simplify() and dotplot() dispatch
    # differently for enrichGO/enrichKEGG
    if (type == "GO") {
        allRes@fun <- "enrichGO"
    } else {
        allRes@fun <- "enrichKEGG"
    }

    # GO simplify — legacy-equivalent, fully offline.
    # GO simplify is MANDATORY for GO results (not configurable). The ORA above
    # (enricher over local TERM2GENE/TERM2NAME) always produces the UNSIMPLIFIED GO
    # result; simplify() then collapses redundant GO terms with the DEFAULT Wang
    # measure, whose semantic similarity is derived from the GO DAG in GO.db —
    # organism-agnostic, needing NO OrgDb (exactly what the legacy
    # Clusters_Enrichment_Test() did: a bare simplify(allRes, cutoff=0.7, ...)).
    # KEGG (type != "GO") is never simplified. The `orgdb` argument is retained
    # ONLY for future IC-based measures (Resnik/Lin/Jiang/Rel), which additionally
    # need an organism OrgDb; it is not used for the Wang default. Fail-soft: on
    # any missing package / failure, warn and keep the unsimplified result — the
    # pipeline never fails.
    allRes_simplify <- NULL
    enrichment_table_simplify <- NULL

    if (type == "GO") {
        sem_data <- .go_semdata(ont)
        if (is.null(sem_data)) {
            warning("GO simplify could not run: GO.db/GOSemSim are not installed. ",
                    "Retaining the unsimplified GO ORA result.")
        } else {
            allRes_simplify <- tryCatch({
                clusterProfiler::simplify(
                    allRes,
                    cutoff     = 0.7,
                    by         = "p.adjust",
                    select_fun = min,
                    measure    = "Wang",
                    semData    = sem_data
                )
            }, error = function(e) {
                message("    simplify() failed: ", e$message)
                NULL
            })
            if (!is.null(allRes_simplify)) {
                enrichment_table_simplify <- process_enrichment_table(
                    allRes_simplify@compareClusterResult
                )
            }
        }
    }

    list(allRes, allRes_simplify, enrichment_table, enrichment_table_simplify)
}

#' Write cluster ORA outputs to disk (serial, deterministic — never in a worker)
#'
#' The I/O half of run_cluster_ora(): writes the enrichment table, the optional
#' simplify table, and the dotplot into the ORA unit directory, plus (optionally)
#' the shared-gene heatmaps under `unit_dir/shared_genes/`. Split out so it can
#' run in the serial assembly step while compute is parallelized. Filenames are
#' short and fixed — the analysis context (db / collection / round) is in the
#' unit directory path (see ora_unit_dir()).
#'
#' @param ora_result The 4-element list returned by run_cluster_ora_compute().
#' @param unit_dir The ORA unit directory (from ora_unit_dir()).
#' @param type "KEGG" or "GO".
#' @param maxCategory Max categories to show in dotplot (default 1000).
#' @param dotplot Logical; when TRUE write the ORA dotplot PDF (the config
#'   `plots.dotplot` toggle — controls ORA and GSEA dotplots alike). Default TRUE.
#' @param shared_genes Logical; when TRUE also emit the legacy ORA shared-gene
#'   heatmaps (gene<->term / term<->term) via plot_ora_shared_genes(). Additive
#'   and fail-soft; default FALSE (callers pass the config toggle).
#' @return Invisibly NULL (called for its file-writing side effects).
#' @noRd
write_cluster_ora_outputs <- function(ora_result, unit_dir,
                                      type = "KEGG", maxCategory = 1000,
                                      dotplot = TRUE,
                                      shared_genes = FALSE) {
    if (length(ora_result) == 0 || is.null(unit_dir)) return(invisible(NULL))

    allRes                    <- ora_result[[1]]
    enrichment_table          <- ora_result[[3]]
    enrichment_table_simplify <- ora_result[[4]]

    dir.create(unit_dir, recursive = TRUE, showWarnings = FALSE)

    # Enrichment result table
    if (!is.null(enrichment_table)) {
        write.csv(x = enrichment_table, file = file.path(unit_dir, "results.csv"),
                  quote = TRUE, row.names = TRUE)
    }

    # Simplify table (only when GO simplify ran and produced a table)
    if (!is.null(enrichment_table_simplify)) {
        write.csv(x = enrichment_table_simplify, file = file.path(unit_dir, "simplify.csv"),
                  quote = FALSE, row.names = TRUE)
    }

    # Dotplot (config-gated: plots.dotplot). ORA and GSEA share this toggle.
    if (isTRUE(dotplot) && !is.null(allRes) && nrow(allRes@compareClusterResult) > 0) {
        # Legacy font size heuristic
        font.size <- if (nrow(allRes@compareClusterResult) > 50) 4 else 9

        dot_plot_file <- file.path(unit_dir, "dotplot.pdf")
        tryCatch({
            if (requireNamespace("enrichplot", quietly = TRUE)) {
                p <- enrichplot::dotplot(allRes, showCategory = maxCategory,
                                         font.size = font.size)
                ggplot2::ggsave(filename = dot_plot_file, plot = p,
                                dpi = 600, device = "pdf",
                                width = 20, height = 20)
            }
        }, error = function(e) {
            message("    Dotplot generation failed: ", e$message)
        })
    }

    # Legacy ORA shared-gene heatmaps (additive, config-gated, fail-soft) →
    # unit_dir/shared_genes/cluster_<label>_{genes_to_terms,terms_to_terms}.{csv,pdf}
    if (isTRUE(shared_genes) && !is.null(allRes)) {
        plot_ora_shared_genes(allRes@compareClusterResult,
                              file.path(unit_dir, "shared_genes"))
    }

    invisible(NULL)
}

#' Plot ORA shared-gene heatmaps (gene<->term and term<->term overlap)
#'
#' Faithful, robust port of the legacy `plot_shared_genes()` (file outputs only).
#' For each cluster in a compareClusterResult table, produces two views:
#'   - gene<->term: a binary term x gene membership matrix (which genes drive
#'     which enriched terms);
#'   - term<->term: a term x term Jaccard overlap matrix (percent shared genes,
#'     `100 * 2|A∩B| / (|A|+|B|)`).
#' Each is rendered with `pheatmap` (PDF) and written as a CSV, reordered by the
#' heatmap's row/column clustering. Large clusters (legacy gates: >20 terms, or
#' >200 genes for the gene<->term view) skip the PDF and write the CSV only.
#'
#' Pure plotting — no enrichment computation, no RNG (deterministic). Call only
#' in the serial output stage (never in a worker). Fail-soft: missing deps, too
#' few terms/genes, or any plotting error → warn and skip, never stop the
#' pipeline. The legacy `plotly`/`htmltools` interactive return is intentionally
#' NOT ported (that is deferred report/UI work).
#'
#' Note (corrected legacy quirk): the term<->term matrix holds overlap *percent*
#' (0–100); legacy hard-coded `breaks = seq(0, 1, 0.01)` (a 0–1 range) which does
#' not cover the data. Here the RdYlBu gradient is mapped across the actual 0–100
#' range so the heatmap renders correctly — a display-only correction.
#'
#' @param ora_df Data.frame with `Cluster`, `ID`, `Description`, `geneID`
#'   (slash-separated genes) — i.e. `allRes@compareClusterResult`.
#' @param out_dir Output directory (the unit's `shared_genes/` folder). Files are
#'   `cluster_<label>_genes_to_terms.{csv,pdf}` and `..._terms_to_terms.{csv,pdf}`
#'   — the analysis context is already in the parent path, so it is not repeated.
#' @return Invisibly a character vector of files written (may be empty).
#' @noRd
plot_ora_shared_genes <- function(ora_df, out_dir) {
    if (is.null(ora_df) || !is.data.frame(ora_df) || nrow(ora_df) == 0 ||
        is.null(out_dir)) {
        return(invisible(character(0)))
    }
    if (!all(c("Cluster", "ID", "Description", "geneID") %in% colnames(ora_df))) {
        warning("plot_ora_shared_genes(): missing required columns; skipping ",
                out_dir)
        return(invisible(character(0)))
    }
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
        warning("pheatmap not installed; skipping ORA shared-gene plots: ", out_dir)
        return(invisible(character(0)))
    }

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Render one matrix view: heatmap PDF (small case only) + a CSV reordered by
    # the heatmap's clustering. Always writes the CSV (falls back to the matrix's
    # original order if the heatmap can't be drawn, e.g. a constant/degenerate
    # matrix), so the data is preserved even when the plot is skipped. `big`
    # suppresses the PDF (legacy behavior for large clusters). Returns the paths
    # written. Fail-soft.
    render_view <- function(mat, view, cl, big, small_args) {
        base  <- paste0("cluster_", cl, "_", view)
        pdf_f <- file.path(out_dir, paste0(base, ".pdf"))
        csv_f <- file.path(out_dir, paste0(base, ".csv"))
        constant <- length(unique(as.vector(mat))) <= 1  # clustering meaningless

        hm <- NULL
        if (!constant) {
            hm <- tryCatch({
                if (big) {
                    pheatmap::pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
                                       silent = TRUE, legend = isTRUE(small_args$legend))
                } else {
                    do.call(pheatmap::pheatmap, c(list(
                        mat = mat, cluster_rows = TRUE, cluster_cols = TRUE,
                        treeheight_col = 0, treeheight_row = 0,
                        width = 5, height = 5, border_color = "white",
                        main = paste("Cluster", cl), filename = pdf_f, silent = TRUE
                    ), small_args))
                }
            }, error = function(e) {
                message("      ", view, " heatmap failed (cluster ", cl, "): ", e$message)
                NULL
            })
        }

        ord_r <- if (!is.null(hm)) hm$tree_row$order else seq_len(nrow(mat))
        ord_c <- if (!is.null(hm)) hm$tree_col$order else seq_len(ncol(mat))
        utils::write.csv(mat[ord_r, ord_c, drop = FALSE], quote = TRUE,
                         row.names = TRUE, file = csv_f)

        if (!big && !is.null(hm) && file.exists(pdf_f)) c(pdf_f, csv_f) else csv_f
    }

    written <- character(0)
    for (cl in unique(ora_df$Cluster)) {
        cl_res <- tryCatch({
            # Keep rows valid for the shared-gene matrices: drop only rows whose
            # REQUIRED fields (Description label + geneID membership) are missing.
            # Unrelated optional columns (e.g. qvalue) may be NA and must NOT cause
            # an otherwise-enriched term to be discarded (a whole-row na.omit did).
            cl_sub <- ora_df[ora_df$Cluster == cl, , drop = FALSE]
            keep   <- !is.na(cl_sub$Description) & !is.na(cl_sub$geneID) &
                      nzchar(as.character(cl_sub$geneID))
            TEMP   <- cl_sub[keep, , drop = FALSE]
            term_genes <- strsplit(as.character(TEMP$geneID), "/", fixed = TRUE)
            all_genes  <- unique(unlist(term_genes))
            n_terms <- nrow(TEMP)

            # Legacy guard: need >1 term and >1 gene to build meaningful matrices.
            # (if/else, not return() — a bare return() inside tryCatch would exit
            # the whole function, skipping the remaining clusters.)
            if (n_terms <= 1 || length(all_genes) <= 1) {
                character(0)
            } else {
                # ---- gene <-> term: binary membership (rows = terms, cols = genes) ----
                g2t <- matrix(0L, nrow = n_terms, ncol = length(all_genes),
                              dimnames = list(TEMP$Description, all_genes))
                for (i in seq_len(n_terms)) g2t[i, term_genes[[i]]] <- 1L
                g2t_color <- if (length(unique(as.vector(g2t))) == 2) {
                    grDevices::colorRampPalette(c("white", "pink"))(2)
                } else {
                    grDevices::colorRampPalette(c("pink"))(2)
                }
                f_g2t <- render_view(
                    g2t, "genes_to_terms", cl,
                    big = n_terms > 20 || length(all_genes) > 200,
                    small_args = list(color = g2t_color, legend = FALSE,
                                      cellheight = 10, fontsize_row = 5, fontsize_col = 1)
                )

                # ---- term <-> term: Jaccard overlap percent (0-100) ----
                t2t <- matrix(0, nrow = n_terms, ncol = n_terms,
                              dimnames = list(TEMP$Description, TEMP$Description))
                for (i in seq_len(n_terms)) for (j in seq_len(n_terms)) {
                    gi <- term_genes[[i]]; gj <- term_genes[[j]]
                    t2t[i, j] <- 100 * 2 * length(intersect(gi, gj)) / (length(gi) + length(gj))
                }
                t2t_color <- grDevices::colorRampPalette(
                    rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
                f_t2t <- render_view(
                    t2t, "terms_to_terms", cl,
                    big = n_terms > 20,
                    small_args = list(color = t2t_color,
                                      breaks = seq(0, 100, length.out = 101),  # corrected range (see @note)
                                      legend = TRUE, fontsize_row = 5, fontsize_col = 5)
                )

                c(f_g2t, f_t2t)
            }
        }, error = function(e) {
            warning("plot_ora_shared_genes(): cluster ", cl, " failed in ",
                    out_dir, ": ", conditionMessage(e))
            character(0)
        })
        written <- c(written, cl_res)
    }

    invisible(written)
}


#' Process clusterProfiler results table with fold enrichment
#'
#' Port of legacy process_clusterprofiler_results_table().
#' Expands GeneRatio and BgRatio into numeric components and computes
#' Fold_enrichment = (in_cluster_in_term / in_cluster) / (in_term / in_genome).
#'
#' @param clusterprofiler_results_table Data.frame from
#'   compareClusterResult@@compareClusterResult or enrichResult@@result
#' @return Data.frame with expanded ratio columns and Fold_enrichment
#' @export
process_enrichment_table <- function(clusterprofiler_results_table) {
    MAX_NR_GENES_TO_SHOW <- 1000
    text_to_show <- paste0("Too many to show (>", MAX_NR_GENES_TO_SHOW, ")")

    et <- clusterprofiler_results_table

    if (nrow(et) == 0) return(et)

    # Split GeneRatio "k/n" into two numeric columns
    if (!"GeneRatio" %in% colnames(et) || !"BgRatio" %in% colnames(et)) {
        warning("process_enrichment_table: missing GeneRatio or BgRatio columns")
        return(et)
    }

    gr_parts <- strsplit(as.character(et$GeneRatio), "/")
    et$in_cluster_in_term <- as.numeric(vapply(gr_parts, `[`, character(1), 1))
    et$in_cluster         <- as.numeric(vapply(gr_parts, `[`, character(1), 2))

    # Split BgRatio "M/N" into two numeric columns
    bg_parts <- strsplit(as.character(et$BgRatio), "/")
    et$in_term    <- as.numeric(vapply(bg_parts, `[`, character(1), 1))
    et$in_genome  <- as.numeric(vapply(bg_parts, `[`, character(1), 2))

    # Fold enrichment: (k/n) / (M/N)
    # Guard against division by zero (produces NaN/Inf)
    denom <- (et$in_term / et$in_genome)
    denom[denom == 0] <- NA
    et$Fold_enrichment <- signif(
        (et$in_cluster_in_term / et$in_cluster) / denom,
        digits = 2
    )

    # Truncate geneID for very large gene lists
    if ("geneID" %in% colnames(et) && "Count" %in% colnames(et)) {
        et$geneID <- ifelse(et$Count <= MAX_NR_GENES_TO_SHOW,
                            et$geneID, text_to_show)
    }

    et
}

# ==============================================================================
# GENE LIST BUILDER FOR ORA (legacy orchestration layer)
# ==============================================================================
# The legacy enrichment workflow iterates over gene_lists[[clust_method]][[clust_round]],
# where each entry is a named vector mapping gene IDs to cluster labels.
# This includes contrast-derived "clusters" (up/down per contrast, all DE per contrast)
# and actual clustering-derived assignments (partition, hierarchical, binary patterns).

#' Build gene_lists structure for cluster-based ORA
#'
#' Constructs a nested list gene_lists[[method]][[round]] where each leaf is
#' a named character/integer vector (gene IDs as names, cluster labels as values).
#' This matches the legacy enrichment orchestration structure.
#'
#' @param de_tables Named list of per-contrast DE tables
#'   (each with FeatureID, log2FoldChange, padj columns)
#' @param clustering_res Result from mod_rnaseq_clustering(), or NULL
#' @param p_cutoff Adjusted p-value cutoff for DE significance (default 0.05)
#' @param lfc_cutoff log2 fold change cutoff for DE significance (default log2(1.5))
#' @return Named list: method -> round -> named vector (gene ID -> cluster label).
#'   Returns empty list if no gene lists can be built.
#' @export
build_gene_lists <- function(de_tables,
                             clustering_res = NULL,
                             p_cutoff = 0.05,
                             lfc_cutoff = log2(1.5)) {

    gene_lists <- list()

    # Match the pipeline's canonical DE significance rule (see
    # R/domain/rnaseq/04_de_summary.R:167-180) exactly, so ORA operates on the same
    # gene set the summary reports as DE: padj <= p_cutoff AND
    # abs(signif(linearFC, 3)) >= linear cutoff, where
    # linearFC = ifelse(lfc >= 0, 2^lfc, -2^-lfc). The caller passes lfc_cutoff in
    # log2 units (log2(linear_fc_cutoff)), so recover the linear cutoff here.
    # Direction uses the sign of the rounded linear FC (also matching the summary).
    linear_cut <- 2 ^ lfc_cutoff
    .sig_rows <- function(dt) {
        lin_fc     <- ifelse(dt$log2FoldChange >= 0, 2 ^ dt$log2FoldChange,
                             -1 * (2 ^ -dt$log2FoldChange))
        rounded_fc <- signif(lin_fc, 3)
        keep <- !is.na(dt$padj) & dt$padj <= p_cutoff &
                !is.na(dt$log2FoldChange) & abs(rounded_fc) >= linear_cut
        out <- dt[keep, , drop = FALSE]
        out$.rounded_fc <- rounded_fc[keep]
        out
    }

    # ------------------------------------------------------------------
    # 1. Contrast-based gene lists (always available when DE tables exist)
    # ------------------------------------------------------------------
    if (length(de_tables) > 0) {

        # "contrasts": per contrast, genes assigned to "up" or "down"
        for (cn in names(de_tables)) {
            dt <- de_tables[[cn]]
            sig <- .sig_rows(dt)
            if (nrow(sig) == 0) next
            # Deduplicate by FeatureID (keep first occurrence)
            sig <- sig[!duplicated(sig$FeatureID), , drop = FALSE]

            labels <- ifelse(sig$.rounded_fc > 0, "up", "down")
            names(labels) <- sig$FeatureID
            gene_lists[["contrasts"]][[cn]] <- labels
        }

        # "contrasts_wo_direction": per contrast, all DE genes in one cluster "all"
        for (cn in names(de_tables)) {
            dt <- de_tables[[cn]]
            sig <- .sig_rows(dt)
            if (nrow(sig) == 0) next
            sig <- sig[!duplicated(sig$FeatureID), , drop = FALSE]

            labels <- rep("all", nrow(sig))
            names(labels) <- sig$FeatureID
            gene_lists[["contrasts_wo_direction"]][[cn]] <- labels
        }

        # "all_DE": union of DE genes across ALL contrasts, single cluster "all"
        # (legacy gene_lists[["all_DE"]][["any_contrast"]]). One ORA run over the
        # combined DE set, independent of contrast or direction.
        all_de_ids <- character(0)
        for (cn in names(de_tables)) {
            dt <- de_tables[[cn]]
            sig <- .sig_rows(dt)
            all_de_ids <- c(all_de_ids, sig$FeatureID)
        }
        all_de_ids <- unique(all_de_ids[!is.na(all_de_ids)])
        if (length(all_de_ids) > 0) {
            gene_lists[["all_DE"]][["any_contrast"]] <-
                setNames(rep("all", length(all_de_ids)), all_de_ids)
        }
    }

    # ------------------------------------------------------------------
    # 2. Clustering-derived gene lists (when clustering results available)
    #
    # ID-space alignment (M1): the DE-derived lists above use de_tables$FeatureID,
    # which the pathway module may have stripped of a leading "Gene:" prefix.
    # Clustering IDs come from clustering_res (unstripped), so we strip the same
    # prefix here to keep partition/binary genes in the SAME ID space as the DE
    # lists and the local TERM2GENE tables. sub() is a no-op when no prefix exists.
    # ------------------------------------------------------------------
    if (!is.null(clustering_res) && is.list(clustering_res$objects)) {
        objs <- clustering_res$objects
        eo   <- clustering_res$excel_order

        # Partition clusters — sourced UNAMBIGUOUSLY (M2).
        # objects$clusters is overwritten (hierarchical -> partition), so it is not
        # a reliable "partition" signal by itself. excel_order$partition_clusters is
        # set ONLY when partition ran (alongside hierarchical) -> prefer it. With no
        # excel_order at all, hierarchical did not run and objects$clusters then holds
        # partition clusters (partition-only run). If excel_order exists but
        # partition_clusters is NULL, partition did NOT run and objects$clusters holds
        # HIERARCHICAL cuts -> do NOT treat them as "partition" (legacy never ran ORA
        # on hierarchical cuts).
        partition_clusters <- NULL
        if (!is.null(eo) && !is.null(eo$partition_clusters)) {
            partition_clusters <- eo$partition_clusters
        } else if (is.null(eo) && !is.null(objs$clusters)) {
            partition_clusters <- objs$clusters
        }
        if (!is.null(partition_clusters) && length(partition_clusters) > 0) {
            names(partition_clusters) <- sub("^Gene:", "", names(partition_clusters))
            gene_lists[["partition"]][["k"]] <- partition_clusters
        }

        # Binary patterns: objects$patterns is a DATA.FRAME with feature_id +
        # best_pattern columns (run_binary_patterns()$best). (M3: guard on nrow(),
        # not length(), since length() of a data.frame counts columns, not rows.)
        if (!is.null(objs$patterns) && is.data.frame(objs$patterns) &&
            nrow(objs$patterns) > 0) {
            bin_pattern <- objs$patterns[!is.na(objs$patterns$best_pattern), , drop = FALSE]
            if (nrow(bin_pattern) > 0) {
                clusters <- setNames(
                    as.character(bin_pattern$best_pattern),
                    sub("^Gene:", "", bin_pattern$feature_id)
                )
                gene_lists[["binary_patterns"]][["best"]] <- clusters
            }
        }
    }

    gene_lists
}
