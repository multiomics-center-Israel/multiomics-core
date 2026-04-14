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

# ==============================================================================
# LOCAL PATHWAY TABLE LOADING (Phase 1 — enrichment migration)
# ==============================================================================

#' Load local precomputed KEGG/GO pathway tables
#'
#' Reads TERM2GENE and TERM2NAME data.frames from tab-separated files
#' in a local annotation directory. No online resources are contacted.
#'
#' Expected files under annotation_dir:
#'   KEGG:  KEGG_pathway2gene.tab, KEGG_pathway2name.tab
#'   GO:    GO2gene_{BP,MF,CC}.tab, GO2name_{BP,MF,CC}.tab
#'
#' Each file is two-column tab-separated (term ID, gene ID or term name).
#'
#' @param annotation_dir Path to directory containing the pathway table files
#' @param databases Character vector of databases to load.
#'   Valid values: "KEGG", "GO_BP", "GO_MF", "GO_CC"
#' @param feature_ids Optional character vector of pipeline feature IDs.
#'   If provided, overlap with TERM2GENE gene columns is checked and warned.
#' @return Named list where each element is list(TERM2GENE = df, TERM2NAME = df)
#' @export
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

        # Read two-column tab files.
        # Legacy uses read.delim() with default header=TRUE, so files may have headers.
        # We read with header=TRUE (matching legacy), then fall back to header=FALSE
        # if the result has fewer than 2 rows (suggesting no header was present).
        term2gene <- read.delim(gene_file, sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE, row.names = NULL)
        term2name <- read.delim(name_file, sep = "\t", header = TRUE,
                                stringsAsFactors = FALSE, row.names = NULL)
        # If header=TRUE produced zero rows, retry without header
        if (nrow(term2gene) == 0) {
            term2gene <- read.delim(gene_file, sep = "\t", header = FALSE,
                                    stringsAsFactors = FALSE, row.names = NULL)
        }
        if (nrow(term2name) == 0) {
            term2name <- read.delim(name_file, sep = "\t", header = FALSE,
                                    stringsAsFactors = FALSE, row.names = NULL)
        }

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
    df$rank_val <- -log10(df$pval)
    # Legacy: replace any NaN/Inf from log10(0) with 0
    df$rank_val[!is.finite(df$rank_val)] <- 0
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
    df$neg_log_p <- -log10(df$pval)
    # Legacy: replace any NaN/Inf with 0
    df$neg_log_p[!is.finite(df$neg_log_p)] <- 0
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
    df$rank_val <- -log10(df$min_pval)
    df$rank_val[!is.finite(df$rank_val)] <- 0
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
            pvalueCutoff  = pvalueCutoff
        )
    }, error = function(e) {
        message("    GSEA failed: ", e$message)
        NULL
    })

    res
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
                         workers = 1) {

    if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # ------------------------------------------------------------------
    # 1. Build flat job list (lightweight identifiers only — data looked
    #    up by reference inside the worker to avoid copying large tables
    #    across serialized futures on Windows/multisession)
    # ------------------------------------------------------------------
    jobs <- list()
    for (ranking_method in names(ranked_genes)) {
        for (contrast in names(ranked_genes[[ranking_method]])) {
            if (length(ranked_genes[[ranking_method]][[contrast]]) == 0) next

            for (db_name in names(local_tables)) {
                jobs[[length(jobs) + 1]] <- list(
                    ranking_method = ranking_method,
                    contrast       = contrast,
                    db_name        = db_name
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
    run_one_gsea_job <- function(job) {
        # Pure computation — no file I/O, no message() (avoids interleaved output).
        # Looks up data from ranked_genes / local_tables by identifier.
        ranked    <- ranked_genes[[job$ranking_method]][[job$contrast]]
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
                pvalueCutoff  = pvalueCutoff
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

    use_parallel <- workers > 1 &&
        requireNamespace("future", quietly = TRUE) &&
        requireNamespace("future.apply", quietly = TRUE)

    if (use_parallel) {
        # Save current plan, set multisession, restore on exit
        old_plan <- future::plan()
        future::plan(future::multisession, workers = workers)
        on.exit(future::plan(old_plan), add = TRUE)

        job_results <- future.apply::future_lapply(
            jobs, run_one_gsea_job,
            future.seed = TRUE
        )
    } else {
        if (workers > 1) {
            message("  future/future.apply not available — running sequentially. ",
                    "Install with: install.packages(c('future', 'future.apply'))")
        }
        job_results <- lapply(jobs, function(job) {
            message("  GSEA: ", job$db_name, " | ", job$ranking_method, " | ", job$contrast)
            run_one_gsea_job(job)
        })
    }

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

        n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
        message("  ", db_name, " | ", ranking_method, " | ", contrast,
                ": ", n_sig, " significant (padj < 0.05)")

        # Write CSV
        if (!is.null(output_dir)) {
            gsea_sub_dir <- file.path(output_dir, db_name,
                                      paste0("ranking_by_", ranking_method),
                                      contrast)
            dir.create(gsea_sub_dir, recursive = TRUE, showWarnings = FALSE)
            csv_file <- file.path(gsea_sub_dir,
                                  paste0("GSEA_results_", contrast, ".csv"))
            write.csv(res_df, file = csv_file, row.names = FALSE)

            # Generate dotplot if significant results exist.
            # Primary: enrichplot::dotplot() on the gseaResult object.
            # Fallback: basic ggplot2 scatterplot (visual approximation only).
            if (n_sig >= 3) {
                plot_file <- file.path(gsea_sub_dir,
                                       paste0("GSEA_dotplot_", contrast, ".png"))
                plot_key <- paste0(db_name, "_", ranking_method, "_", contrast)
                show_n <- min(20, n_sig)

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
        }
    }

    list(results = results, plot_files = plot_files)
}
