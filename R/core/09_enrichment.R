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
