#' Cross-omics pathway enrichment analysis
#'
#' Combines enrichment results from multiple omics layers to identify
#' pathways that are consistently dysregulated across modalities.


# =============================================================================
# Build per-omics enrichment from DE results
# =============================================================================

#' Build per-omics enrichment results
#'
#' Uses pre-computed enrichment if available. Otherwise runs KEGG enrichment
#' from DE tables for omics layers that have DE results.
#'
#' @param enrichment_results Pre-computed enrichment results (may be NULL/incomplete)
#' @param de_results Named list of DE results per omics
#' @param harmonization_res Harmonization result with MAE and pre-processing data
#' @param config Full config object
#' @param out_dir Output directory
#' @return Named list of enrichment data frames per omics (pathway, pvalue columns)
build_per_omics_enrichment <- function(enrichment_results, de_results,
                                        harmonization_res, config, out_dir) {

    organism <- config$global$organism
    omics_present <- config$global$omics_present
    kegg_org <- get_kegg_organism(organism)

    # Pre-compute KEGG ID conversion once (expensive KEGG REST API call)
    kegg_conv_cache <- NULL
    if (!is.null(kegg_org)) {
        message("  Building KEGG ID conversion table (one-time)...")
        kegg_conv_cache <- tryCatch(
            build_kegg_conversion_cache(de_results, harmonization_res, organism, kegg_org),
            error = function(e) {
                message("  KEGG conversion cache failed: ", e$message)
                NULL
            }
        )
    }

    per_omics <- list()

    for (om in omics_present) {
        # Check if pre-computed enrichment has usable results
        precomp <- enrichment_results[[om]]
        precomp_df <- extract_enrichment_df(precomp)

        if (!is.null(precomp_df) && nrow(precomp_df) > 0) {
            message("  ", om, ": using pre-computed enrichment (", nrow(precomp_df), " pathways)")
            per_omics[[om]] <- precomp_df
            next
        }

        # Run enrichment from DE results
        de_data <- de_results[[om]]
        if (is.null(de_data)) {
            message("  ", om, ": no DE results available, skipping enrichment")
            next
        }

        message("  ", om, ": running KEGG enrichment from DE results...")
        enrich_df <- tryCatch(
            run_kegg_enrichment_for_omics(
                de_data = de_data,
                omics_type = om,
                harmonization_res = harmonization_res,
                organism = organism,
                config = config,
                out_dir = file.path(out_dir, om),
                kegg_conv_cache = kegg_conv_cache
            ),
            error = function(e) {
                warning("  ", om, " enrichment failed: ", e$message)
                NULL
            }
        )

        if (!is.null(enrich_df) && nrow(enrich_df) > 0) {
            message("    Found ", nrow(enrich_df), " enriched pathways")
            per_omics[[om]] <- enrich_df
        } else {
            message("    No enriched pathways found")
        }
    }

    per_omics
}


#' Build KEGG ID conversion cache (ENTREZID -> KEGG gene ID)
#'
#' Downloads the conversion table once for all omics.
build_kegg_conversion_cache <- function(de_results, harmonization_res, organism, kegg_org) {
    org_db <- get_organism_db(organism)
    if (is.null(org_db)) return(NULL)

    # Collect all ENTREZID from all omics
    all_entrez <- c()

    for (om in names(de_results)) {
        de_data <- de_results[[om]]
        de_tables <- extract_de_tables(de_data, om, harmonization_res)
        if (length(de_tables) == 0) next

        id_map <- tryCatch(
            map_feature_ids_to_entrez(de_tables, om, harmonization_res, org_db),
            error = function(e) NULL
        )
        if (!is.null(id_map)) {
            all_entrez <- c(all_entrez, id_map$ENTREZID)
        }
    }

    all_entrez <- unique(all_entrez[!is.na(all_entrez)])
    if (length(all_entrez) == 0) return(NULL)

    message("    Converting ", length(all_entrez), " unique ENTREZID to KEGG IDs...")
    convert_entrez_to_kegg(all_entrez, kegg_org)
}


#' Extract a usable enrichment data frame from various result formats
extract_enrichment_df <- function(enrich_res) {
    if (is.null(enrich_res)) return(NULL)

    # Direct data frame
    if (is.data.frame(enrich_res)) return(enrich_res)

    # $enrichment_df slot
    if (!is.null(enrich_res$enrichment_df) && is.data.frame(enrich_res$enrichment_df)) {
        return(enrich_res$enrichment_df)
    }

    # Contrast-keyed collection. The per-contrast results live either under a
    # $pathway_results slot (RNA) or at the top level of the list keyed by
    # contrast name (proteomics: list(<contrast> = list(<db>_fgsea = df, ...))).
    # A top-level list only counts as contrast-keyed when every element is a list
    # of data frames. The metabolomics wrapper (qea, ora, ..., plots, files) fails
    # that on purpose: its method tables use other p-value columns, and taking
    # them here would skip the KEGG-from-DE fallback for that layer.
    is_contrast_keyed <- is.list(enrich_res) && length(enrich_res) > 0 &&
        all(vapply(enrich_res, function(x) {
            is.list(x) && !is.data.frame(x) &&
                all(vapply(x, function(y) is.null(y) || is.data.frame(y), logical(1)))
        }, logical(1)))

    contrast_list <- if (!is.null(enrich_res$pathway_results) &&
                         is.list(enrich_res$pathway_results)) {
        enrich_res$pathway_results
    } else if (is_contrast_keyed) {
        enrich_res
    } else {
        NULL
    }

    if (!is.null(contrast_list)) {
        dfs <- list()
        for (contrast_name in names(contrast_list)) {
            contrast_res <- contrast_list[[contrast_name]]
            if (is.data.frame(contrast_res) && nrow(contrast_res) > 0) {
                dfs[[contrast_name]] <- contrast_res
            } else if (is.list(contrast_res)) {
                # May have sub-results (e.g., custom_fgsea, KEGG, GO)
                for (sub_name in names(contrast_res)) {
                    sub_res <- contrast_res[[sub_name]]
                    if (is.data.frame(sub_res) && nrow(sub_res) > 0) {
                        dfs[[paste0(contrast_name, "_", sub_name)]] <- sub_res
                    }
                }
            }
        }

        # ORA and GSEA tables carry legitimately different columns (e.g. ORA has
        # Fold_enrichment/Count, GSEA has NES/core_enrichment). bind_rows() aligns
        # by name and NA-fills the missing method-specific columns, whereas rbind()
        # requires identical schemas and aborts on heterogeneous inputs.
        if (length(dfs) > 0) return(dplyr::bind_rows(dfs))
    }

    NULL
}


# =============================================================================
# Run KEGG enrichment from DE results
# =============================================================================

#' Run KEGG enrichment for a single omics layer
#'
#' Extracts gene-level stats from DE results, maps IDs to ENTREZ,
#' and runs ORA or GSEA via clusterProfiler.
run_kegg_enrichment_for_omics <- function(de_data, omics_type, harmonization_res,
                                           organism, config, out_dir,
                                           kegg_conv_cache = NULL) {

    kegg_org <- get_kegg_organism(organism)

    # For gene-based omics, need org_db and kegg_org
    if (omics_type != "metabolomics") {
        org_db <- get_organism_db(organism)
        if (is.null(org_db) || is.null(kegg_org)) {
            message("    Organism annotation not available for: ", organism)
            return(NULL)
        }
    }

    # Extract DE tables with gene-level stats
    de_tables <- extract_de_tables(de_data, omics_type, harmonization_res)
    if (is.null(de_tables) || length(de_tables) == 0) {
        message("    No DE tables found for ", omics_type)
        return(NULL)
    }

    # Get ID mapping: gene-based or compound-based
    if (omics_type == "metabolomics") {
        id_map <- map_metabolite_ids_to_kegg(de_tables, harmonization_res)
        if (is.null(id_map) || nrow(id_map) == 0) {
            message("    Could not map metabolite IDs to KEGG compound IDs")
            return(NULL)
        }
        # KEGG_ID = compound ID directly
        id_map$KEGG_ID <- id_map$KEGG_CPD
    } else {
        id_map <- map_feature_ids_to_entrez(
            de_tables = de_tables,
            omics_type = omics_type,
            harmonization_res = harmonization_res,
            org_db = org_db
        )
        if (is.null(id_map) || nrow(id_map) == 0) {
            message("    Could not map feature IDs to ENTREZID for ", omics_type)
            return(NULL)
        }

        # Convert ENTREZID to KEGG gene IDs (needed for some organisms like C. elegans)
        if (!is.null(kegg_conv_cache)) {
            id_map$KEGG_ID <- kegg_conv_cache[id_map$ENTREZID]
        } else {
            all_entrez <- unique(id_map$ENTREZID[!is.na(id_map$ENTREZID)])
            kegg_conv <- convert_entrez_to_kegg(all_entrez, kegg_org)
            if (!is.null(kegg_conv)) {
                id_map$KEGG_ID <- kegg_conv[id_map$ENTREZID]
            } else {
                id_map$KEGG_ID <- id_map$ENTREZID
            }
        }
    }

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Run enrichment for each contrast
    all_results <- list()
    enrich_cfg <- config$modes$multiomics$enrichment
    method <- enrich_cfg$methods %||% "gsea"
    pval_cutoff <- enrich_cfg$ora_pvalue %||% 0.1
    min_gs <- enrich_cfg$min_set_size %||% 10
    max_gs <- enrich_cfg$max_set_size %||% 500

    for (contrast_name in names(de_tables)) {
        de_df <- de_tables[[contrast_name]]

        # Merge with KEGG IDs
        de_mapped <- merge(de_df, id_map, by = "feature_id")
        de_mapped <- de_mapped[!is.na(de_mapped$KEGG_ID) & !duplicated(de_mapped$KEGG_ID), ]

        if (omics_type == "metabolomics" && nrow(de_mapped) < 3) {
            message("    ", contrast_name, ": too few mapped metabolites (", nrow(de_mapped), ")")
            next
        } else if (omics_type != "metabolomics" && nrow(de_mapped) < 10) {
            message("    ", contrast_name, ": too few mapped features (", nrow(de_mapped), ")")
            next
        }

        message("    ", contrast_name, ": ", nrow(de_mapped), " mapped features")

        # For small datasets (< 500 features), prefer ORA over GSEA
        use_method <- method
        use_min_gs <- min_gs
        if (nrow(de_mapped) < 500) {
            use_method <- "ora"
            use_min_gs <- max(3, min(min_gs, 5))
        }

        result <- tryCatch({
            if (omics_type == "metabolomics") {
                # Use ALL measured KEGG-mapped metabolites as universe (not just DE table)
                full_universe <- unique(id_map$KEGG_ID[!is.na(id_map$KEGG_ID)])
                run_compound_ora(de_mapped, out_dir, use_min_gs, max_gs, pval_cutoff,
                                 universe = full_universe,
                                 exclude_classes = .excluded_pathway_classes(config))
            } else if (use_method == "gsea") {
                run_gsea_kegg(de_mapped, kegg_org, use_min_gs, max_gs, pval_cutoff)
            } else {
                run_ora_kegg(de_mapped, kegg_org, use_min_gs, max_gs, pval_cutoff)
            }
        }, error = function(e) {
            message("    ", contrast_name, " enrichment error: ", e$message)
            NULL
        })

        if (!is.null(result) && nrow(result) > 0) {
            result$contrast <- contrast_name
            result$omics <- omics_type
            all_results[[contrast_name]] <- result
        }
    }

    if (length(all_results) == 0) return(NULL)

    combined <- do.call(rbind, all_results)
    rownames(combined) <- NULL

    # Save per-omics results
    write.csv(combined, file.path(out_dir, paste0(omics_type, "_kegg_enrichment.csv")),
              row.names = FALSE)

    combined
}


#' Extract standardized DE tables from various omics result formats
#'
#' For proteomics with multi-imputation summary_df (numeric IDs), runs limma
#' on the actual expression matrix instead.
extract_de_tables <- function(de_data, omics_type, harmonization_res = NULL) {
    tables <- list()

    if (omics_type == "transcriptomics") {
        # RNA-seq: de_data$tables is a named list of data frames
        if (!is.null(de_data$tables)) {
            for (nm in names(de_data$tables)) {
                df <- de_data$tables[[nm]]
                if (!is.data.frame(df)) next

                # Standardize column names
                std <- data.frame(
                    feature_id = if ("FeatureID" %in% names(df)) df$FeatureID else rownames(df),
                    log2fc = if ("log2FoldChange" %in% names(df)) df$log2FoldChange
                             else if ("logFC" %in% names(df)) df$logFC else NA,
                    pvalue = if ("pvalue" %in% names(df)) df$pvalue
                             else if ("P.Value" %in% names(df)) df$P.Value else NA,
                    padj = if ("padj" %in% names(df)) df$padj
                           else if ("adj.P.Val" %in% names(df)) df$adj.P.Val else NA,
                    stringsAsFactors = FALSE
                )
                std <- std[!is.na(std$pvalue), ]
                if (nrow(std) > 0) tables[[nm]] <- std
            }
        }

    } else if (omics_type == "proteomics") {
        # Use precomputed DE from summary_df, mapping numeric IDs to UniProt
        tables <- extract_proteomics_de_tables(de_data, harmonization_res)

    } else if (omics_type == "metabolomics") {
        # Metabolomics: de_tables named list
        if (!is.null(de_data$de_tables)) {
            for (nm in names(de_data$de_tables)) {
                df <- de_data$de_tables[[nm]]
                if (!is.data.frame(df)) next

                std <- data.frame(
                    feature_id = if ("feature_id" %in% names(df)) df$feature_id else rownames(df),
                    log2fc = if ("logFC" %in% names(df)) df$logFC else NA,
                    pvalue = if ("P.Value" %in% names(df)) df$P.Value else NA,
                    padj = if ("adj.P.Val" %in% names(df)) df$adj.P.Val else NA,
                    stringsAsFactors = FALSE
                )
                std <- std[!is.na(std$pvalue), ]
                if (nrow(std) > 0) tables[[nm]] <- std
            }
        }
    }

    tables
}


#' Extract precomputed DE tables from proteomics results
#'
#' Parses the multi-imputation summary_df, resolves numeric FeatureIDs to
#' UniProt protein IDs via row_data, and converts linearFC to log2FC.
#' Falls back to re-running limma if summary_df is unavailable.
extract_proteomics_de_tables <- function(de_data, harmonization_res) {
    sdf <- de_data$summary_df
    if (is.null(sdf) || nrow(sdf) == 0) {
        return(run_limma_for_proteomics(harmonization_res))
    }

    # Build numeric-ID -> UniProt mapping from row_data
    id_map <- NULL
    prot_pre <- harmonization_res$inputs$proteomics
    if (!is.null(prot_pre) && !is.null(prot_pre$row_data)) {
        rd <- prot_pre$row_data
        # row_data rows correspond 1:1 with expr_work rows (both 7251)
        id_map <- setNames(rownames(prot_pre$expr_work), seq_len(nrow(rd)))
    }

    # Find contrast columns by their padj pattern
    padj_cols <- grep("^padj\\.imputs\\.", colnames(sdf), value = TRUE)
    if (length(padj_cols) == 0) {
        return(run_limma_for_proteomics(harmonization_res))
    }

    tables <- list()
    for (padj_col in padj_cols) {
        # Derive contrast name: "padj.imputs.1.56ppmvs.0ppm" -> "1.56ppm vs. 0ppm"
        contrast_key <- sub("^padj\\.imputs\\.", "", padj_col)
        # Insert space before "vs" and restore the dot: "1.56ppmvs.0ppm" -> "1.56ppm vs. 0ppm"
        contrast_name <- sub("vs\\.", " vs. ", contrast_key)

        pval_col <- sub("^padj\\.", "pvalue.", padj_col)
        fc_col <- sub("^padj\\.", "linearFC.", padj_col)

        if (!pval_col %in% colnames(sdf) || !fc_col %in% colnames(sdf)) next

        # Map numeric FeatureIDs to UniProt IDs
        feat_ids <- as.character(sdf$FeatureID)
        if (!is.null(id_map)) {
            resolved <- id_map[feat_ids]
            feat_ids <- ifelse(is.na(resolved), feat_ids, resolved)
        }

        # linearFC is signed linear fold change; convert to log2
        linear_fc <- sdf[[fc_col]]
        log2fc <- ifelse(linear_fc >= 0, log2(abs(linear_fc)), -log2(abs(linear_fc)))
        log2fc[is.na(linear_fc)] <- NA

        std <- data.frame(
            feature_id = feat_ids,
            log2fc = log2fc,
            pvalue = sdf[[pval_col]],
            padj = sdf[[padj_col]],
            stringsAsFactors = FALSE
        )
        std <- std[!is.na(std$pvalue), ]
        if (nrow(std) > 0) tables[[contrast_name]] <- std
    }

    if (length(tables) == 0) {
        return(run_limma_for_proteomics(harmonization_res))
    }

    tables
}


#' Run limma DE on proteomics expression matrix
#'
#' Fallback when precomputed DE is unavailable. Uses the actual protein
#' expression data rather than the multi-imputation summary.
run_limma_for_proteomics <- function(harmonization_res) {
    if (is.null(harmonization_res) || is.null(harmonization_res$inputs$proteomics)) {
        return(list())
    }

    prot <- harmonization_res$inputs$proteomics
    expr_mat <- prot$expr_work
    meta <- prot$meta

    if (is.null(expr_mat) || is.null(meta)) return(list())

    # Find condition column
    cond_col <- intersect(c("Treatment", "condition", "Condition", "group"), colnames(meta))
    if (length(cond_col) == 0) return(list())
    cond_col <- cond_col[1]

    conditions <- meta[[cond_col]]
    if (length(unique(conditions)) < 2) return(list())

    # Simple limma analysis
    if (!requireNamespace("limma", quietly = TRUE)) return(list())

    cond_factor <- factor(conditions)
    safe_levels <- make.names(levels(cond_factor))
    design <- stats::model.matrix(~ 0 + cond_factor)
    colnames(design) <- safe_levels

    fit <- limma::lmFit(expr_mat, design)

    # Map original levels to safe names
    level_map <- setNames(safe_levels, levels(cond_factor))
    ref_safe <- level_map[levels(cond_factor)[1]]
    other_safe <- setdiff(safe_levels, ref_safe)

    tables <- list()
    for (i in seq_along(other_safe)) {
        lv_safe <- other_safe[i]
        # Original level name for the contrast label
        orig_level <- names(level_map)[level_map == lv_safe]
        orig_ref <- names(level_map)[level_map == ref_safe]

        contrast_str <- paste0(lv_safe, " - ", ref_safe)
        contrast_name <- paste0(orig_level, " vs ", orig_ref)

        contrast_mat <- limma::makeContrasts(contrasts = contrast_str, levels = design)
        fit2 <- limma::contrasts.fit(fit, contrast_mat)
        fit2 <- limma::eBayes(fit2)

        tt <- limma::topTable(fit2, number = Inf, sort.by = "none")

        tables[[contrast_name]] <- data.frame(
            feature_id = rownames(tt),
            log2fc = tt$logFC,
            pvalue = tt$P.Value,
            padj = tt$adj.P.Val,
            stringsAsFactors = FALSE
        )
    }

    tables
}


#' Map feature IDs to KEGG gene IDs
#'
#' For C. elegans (and some other organisms), KEGG uses organism-specific gene
#' IDs (e.g. CELE_xxx) rather than NCBI ENTREZID. This function maps via:
#' feature_id -> ENTREZID -> KEGG gene ID (via bitr_kegg).
map_feature_ids_to_entrez <- function(de_tables, omics_type, harmonization_res, org_db) {

    # Collect all unique feature IDs
    all_ids <- unique(unlist(lapply(de_tables, function(df) df$feature_id)))

    if (omics_type == "transcriptomics") {
        # WBGene IDs -> ENTREZID directly via org.db
        entrez_df <- tryCatch({
            res <- AnnotationDbi::mapIds(
                org_db,
                keys = all_ids,
                keytype = "ENSEMBL",
                column = "ENTREZID",
                multiVals = "first"
            )
            data.frame(
                feature_id = names(res),
                ENTREZID = as.character(res),
                stringsAsFactors = FALSE
            )
        }, error = function(e) {
            tryCatch({
                res <- AnnotationDbi::mapIds(
                    org_db,
                    keys = all_ids,
                    keytype = "WORMBASE",
                    column = "ENTREZID",
                    multiVals = "first"
                )
                data.frame(
                    feature_id = names(res),
                    ENTREZID = as.character(res),
                    stringsAsFactors = FALSE
                )
            }, error = function(e2) NULL)
        })

        return(entrez_df)

    } else if (omics_type == "proteomics") {
        # Try direct UniProt -> ENTREZID mapping first (works for most organisms)
        entrez_df <- tryCatch({
            res <- AnnotationDbi::mapIds(
                org_db,
                keys = all_ids,
                keytype = "UNIPROT",
                column = "ENTREZID",
                multiVals = "first"
            )
            df <- data.frame(
                feature_id = names(res),
                ENTREZID = as.character(res),
                stringsAsFactors = FALSE
            )
            df <- df[!is.na(df$ENTREZID), ]
            if (nrow(df) > 0) {
                message("    Mapped ", nrow(df), "/", length(all_ids),
                        " UniProt IDs to ENTREZID directly")
            }
            df
        }, error = function(e) NULL)

        if (!is.null(entrez_df) && nrow(entrez_df) > 0) return(entrez_df)

        # Fallback: try via row_data WormBase/gene_id columns (C. elegans etc.)
        prot_pre <- harmonization_res$inputs$proteomics
        if (is.null(prot_pre) || is.null(prot_pre$row_data)) {
            message("    No proteomics row_data for ID mapping")
            return(NULL)
        }

        row_data <- prot_pre$row_data
        wbgene_col <- intersect(c("Wormbase_id", "wormbase_id", "gene_id"), colnames(row_data))
        if (length(wbgene_col) == 0) {
            message("    No WormBase/gene_id column in proteomics row_data")
            return(NULL)
        }

        prot_ids <- rownames(prot_pre$expr_work)
        wb_ids <- row_data[[wbgene_col[1]]]
        prot_to_wb <- data.frame(
            feature_id = prot_ids,
            WBGene = wb_ids,
            stringsAsFactors = FALSE
        )
        prot_to_wb <- prot_to_wb[!is.na(prot_to_wb$WBGene) & nzchar(prot_to_wb$WBGene), ]

        if (nrow(prot_to_wb) == 0) return(NULL)

        mapped <- tryCatch({
            res <- AnnotationDbi::mapIds(
                org_db,
                keys = unique(prot_to_wb$WBGene),
                keytype = "ENSEMBL",
                column = "ENTREZID",
                multiVals = "first"
            )
            wb_to_entrez <- data.frame(
                WBGene = names(res),
                ENTREZID = as.character(res),
                stringsAsFactors = FALSE
            )
            prot_mapped <- merge(prot_to_wb, wb_to_entrez, by = "WBGene")
            data.frame(
                feature_id = prot_mapped$feature_id,
                ENTREZID = prot_mapped$ENTREZID,
                stringsAsFactors = FALSE
            )
        }, error = function(e) NULL)

        return(mapped)

    } else if (omics_type == "metabolomics") {
        # Metabolomics: handled separately via map_metabolite_ids_to_kegg()
        message("    Metabolomics uses compound-level enrichment (handled separately)")
        return(NULL)
    }

    NULL
}


#' Convert ENTREZID to KEGG gene IDs
#'
#' For some organisms (e.g. C. elegans), KEGG uses organism-specific IDs
#' (CELE_xxx) rather than ENTREZID. Uses KEGG REST API directly.
convert_entrez_to_kegg <- function(entrez_ids, kegg_org) {
    # Try clusterProfiler first (fast, cached)
    conv <- tryCatch({
        res <- clusterProfiler::bitr_kegg(
            entrez_ids,
            fromType = "ncbi-geneid",
            toType = "kegg",
            organism = kegg_org
        )
        if (nrow(res) > 0) {
            message("    Converted ", nrow(res), "/", length(entrez_ids),
                    " ENTREZID to KEGG IDs (via clusterProfiler)")
            return(setNames(res$kegg, res[["ncbi-geneid"]]))
        }
        NULL
    }, error = function(e) NULL)

    if (!is.null(conv)) return(conv)

    # Fallback: KEGG REST API
    message("    Using KEGG REST API for ID conversion...")
    url <- paste0("https://rest.kegg.jp/conv/", kegg_org, "/ncbi-geneid")
    lines <- tryCatch(readLines(url, warn = FALSE), error = function(e) {
        message("    KEGG REST API error: ", e$message)
        NULL
    })
    if (is.null(lines) || length(lines) == 0) return(NULL)

    parts <- strsplit(lines, "\t")
    conv_df <- data.frame(
        ncbi = gsub("^ncbi-geneid:", "", vapply(parts, `[`, character(1), 1)),
        kegg = vapply(parts, `[`, character(1), 2),
        stringsAsFactors = FALSE
    )

    # Filter to our IDs
    conv_df <- conv_df[conv_df$ncbi %in% entrez_ids, ]
    if (nrow(conv_df) == 0) return(NULL)

    message("    Converted ", nrow(conv_df), "/", length(entrez_ids),
            " ENTREZID to KEGG IDs (via REST API)")
    setNames(conv_df$kegg, conv_df$ncbi)
}


# =============================================================================
# Metabolomics compound enrichment
# =============================================================================

#' Load HMDB-to-KEGG compound ID mapping
#'
#' Reads the tab-separated mapping file (columns: HMDB, KEGG) and returns
#' a named vector for lookup (HMDB -> KEGG compound ID).
#'
#' @param mapping_file Path to HMDB-to-KEGG mapping file. If NULL, looks for
#'   the default file at data/HMDB2kegg_cpd.Jan2026.v2.txt
#' @return Named character vector: names = HMDB IDs, values = KEGG compound IDs
load_hmdb_to_kegg_map <- function(mapping_file = NULL) {
    if (is.null(mapping_file)) {
        # Try standard locations relative to project root
        candidates <- c(
            "data/HMDB2kegg_cpd.Jan2026.v2.txt",
            file.path(Sys.getenv("PIPELINE_ROOT", "."), "data/HMDB2kegg_cpd.Jan2026.v2.txt")
        )
        for (f in candidates) {
            if (file.exists(f)) {
                mapping_file <- f
                break
            }
        }
    }

    if (is.null(mapping_file) || !file.exists(mapping_file)) {
        message("    HMDB-to-KEGG mapping file not found")
        return(NULL)
    }

    df <- utils::read.delim(mapping_file, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE, na.strings = c("NA", ""))
    df <- df[!is.na(df$HMDB) & !is.na(df$KEGG) & nzchar(df$KEGG), ]

    if (nrow(df) == 0) return(NULL)

    message("    Loaded ", nrow(df), " HMDB-to-KEGG compound mappings")
    stats::setNames(df$KEGG, df$HMDB)
}


#' Map metabolite feature IDs to KEGG compound IDs
#'
#' Uses row_data HMDB column and the HMDB-to-KEGG mapping file to convert
#' metabolite feature IDs to KEGG compound IDs (C00xxx format).
#'
#' @param de_tables Named list of DE tables
#' @param harmonization_res Harmonization result with inputs$metabolomics
#' @return data.frame with columns: feature_id, KEGG_CPD
map_metabolite_ids_to_kegg <- function(de_tables, harmonization_res) {
    metab_pre <- harmonization_res$inputs$metabolomics
    if (is.null(metab_pre) || is.null(metab_pre$row_data)) {
        message("    No metabolomics row_data for ID mapping")
        return(NULL)
    }

    row_data <- metab_pre$row_data

    # Find HMDB column
    hmdb_col <- intersect(c("HMDB", "hmdb", "HMDB_ID"), colnames(row_data))
    if (length(hmdb_col) == 0) {
        message("    No HMDB column in metabolomics row_data")
        return(NULL)
    }

    # Load HMDB -> KEGG compound mapping
    hmdb_kegg <- load_hmdb_to_kegg_map()
    if (is.null(hmdb_kegg)) return(NULL)

    # Collect all unique feature_ids from DE tables
    all_de_ids <- unique(unlist(lapply(de_tables, function(df) df$feature_id)))

    # Check if DE feature_ids are already HMDB IDs (direct mapping, no row_data needed)
    hmdb_pattern <- sum(grepl("^HMDB\\d+$", all_de_ids, ignore.case = TRUE), na.rm = TRUE)
    uses_hmdb_direct <- hmdb_pattern > length(all_de_ids) * 0.5

    hmdb_ids <- as.character(row_data[[hmdb_col[1]]])

    if (uses_hmdb_direct) {
        # DE tables already use HMDB IDs as feature_ids — direct mapping
        message("    DE tables use HMDB IDs as feature_ids (direct mapping)")
        feat_ids <- hmdb_ids
        kegg_cpds <- hmdb_kegg[hmdb_ids]
    } else {
        # Determine whether DE feature_ids are metabolite names or synthetic IDs.
        # DE tables may use metabolite names (e.g., "Palmitoleic acid") while
        # row_data$feature_id uses synthetic "feature_N" IDs.
        name_col <- intersect(c("Metabolite", "Name", "Molecule"), colnames(row_data))
        synthetic_ids <- if ("feature_id" %in% colnames(row_data)) row_data$feature_id else NULL

        # Check which key the DE tables are using
        uses_names <- FALSE
        uses_bare_numeric <- FALSE
        if (length(name_col) > 0) {
            metab_names <- as.character(row_data[[name_col[1]]])
            overlap_names <- sum(all_de_ids %in% metab_names, na.rm = TRUE)
            overlap_synth <- if (!is.null(synthetic_ids)) sum(all_de_ids %in% synthetic_ids, na.rm = TRUE) else 0
            uses_names <- overlap_names > overlap_synth
        }

        # DE tables may use bare row indices ("1","2",...) from limma
        if (!uses_names && !is.null(synthetic_ids) &&
            sum(all_de_ids %in% synthetic_ids, na.rm = TRUE) == 0) {
            bare_indices <- as.character(seq_len(nrow(row_data)))
            if (sum(all_de_ids %in% bare_indices, na.rm = TRUE) > length(all_de_ids) * 0.5) {
                uses_bare_numeric <- TRUE
                message("    DE tables use bare numeric row indices as feature_ids")
            }
        }

        if (uses_names) {
            # DE tables use metabolite names; build name -> HMDB -> KEGG
            feat_ids <- metab_names
            message("    DE tables use metabolite names as feature_ids")
        } else if (uses_bare_numeric) {
            # DE tables use bare numeric indices; use same for join
            feat_ids <- as.character(seq_len(nrow(row_data)))
        } else if (!is.null(synthetic_ids)) {
            feat_ids <- synthetic_ids
        } else {
            feat_ids <- rownames(row_data)
            if (is.null(feat_ids)) feat_ids <- paste0("feature_", seq_len(nrow(row_data)))
        }

        # Map HMDB -> KEGG compound
        kegg_cpds <- hmdb_kegg[hmdb_ids]
    }

    mapped <- data.frame(
        feature_id = feat_ids,
        HMDB = hmdb_ids,
        KEGG_CPD = as.character(kegg_cpds),
        stringsAsFactors = FALSE
    )
    mapped <- mapped[!is.na(mapped$KEGG_CPD) & nzchar(mapped$KEGG_CPD), ]

    # Deduplicate by feature_id (keep first)
    mapped <- mapped[!duplicated(mapped$feature_id), ]

    if (nrow(mapped) == 0) {
        message("    No metabolites mapped to KEGG compound IDs")
        return(NULL)
    }

    n_total <- length(unique(feat_ids))
    message("    Mapped ", nrow(mapped), "/", n_total,
            " metabolites to KEGG compound IDs")

    mapped[, c("feature_id", "KEGG_CPD")]
}


#' Get KEGG compound-pathway associations
#'
#' Downloads the compound-pathway mapping from KEGG REST API and caches
#' the result as an RDS file for subsequent runs.
#'
#' @param cache_dir Directory to cache the downloaded data (NULL to skip caching)
#' @return data.frame with columns: pathway, compound, name (pathway description)
get_kegg_compound_pathways <- function(cache_dir = NULL) {
    # Check cache first
    if (!is.null(cache_dir)) {
        cache_file <- file.path(cache_dir, "kegg_compound_pathways.rds")
        if (file.exists(cache_file)) {
            message("    Using cached compound-pathway associations")
            return(readRDS(cache_file))
        }
    }

    message("    Downloading KEGG compound-pathway associations...")

    # Download compound -> pathway links
    link_lines <- tryCatch(
        readLines("https://rest.kegg.jp/link/pathway/compound", warn = FALSE),
        error = function(e) {
            message("    KEGG REST API error: ", e$message)
            NULL
        }
    )

    if (is.null(link_lines) || length(link_lines) == 0) return(NULL)

    link_parts <- strsplit(link_lines, "\t")
    cpd_pathway <- data.frame(
        compound = gsub("^cpd:", "", vapply(link_parts, `[`, character(1), 1)),
        pathway = gsub("^path:", "", vapply(link_parts, `[`, character(1), 2)),
        stringsAsFactors = FALSE
    )

    # Keep only global "map" pathways (organism-independent compound pathways)
    cpd_pathway <- cpd_pathway[grepl("^map", cpd_pathway$pathway), ]

    if (nrow(cpd_pathway) == 0) return(NULL)

    # Download pathway names
    name_lines <- tryCatch(
        readLines("https://rest.kegg.jp/list/pathway", warn = FALSE),
        error = function(e) NULL
    )

    if (!is.null(name_lines) && length(name_lines) > 0) {
        name_parts <- strsplit(name_lines, "\t")
        pathway_names <- data.frame(
            pathway = gsub("^path:", "", vapply(name_parts, `[`, character(1), 1)),
            name = vapply(name_parts, `[`, character(1), 2),
            stringsAsFactors = FALSE
        )
        cpd_pathway <- merge(cpd_pathway, pathway_names, by = "pathway", all.x = TRUE)
    } else {
        cpd_pathway$name <- cpd_pathway$pathway
    }

    message("    Downloaded ", nrow(cpd_pathway), " compound-pathway associations (",
            length(unique(cpd_pathway$pathway)), " pathways)")

    # Cache for future use
    if (!is.null(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        saveRDS(cpd_pathway, file.path(cache_dir, "kegg_compound_pathways.rds"))
    }

    cpd_pathway
}


#' Run ORA for metabolite compounds against KEGG pathways
#'
#' Uses clusterProfiler::enricher() with custom TERM2GENE built from
#' KEGG compound-pathway associations.
#'
#' @param de_mapped DE table with KEGG_ID (compound IDs) column
#' @param cache_dir Directory for caching KEGG data
#' @param min_gs Minimum gene set size
#' @param max_gs Maximum gene set size
#' @param pval_cutoff P-value cutoff
#' @return data.frame with enrichment results, or NULL
run_compound_ora <- function(de_mapped, cache_dir, min_gs, max_gs, pval_cutoff,
                              universe = NULL, exclude_classes = NULL) {

    # Get compound-pathway associations
    cpd_pathways <- get_kegg_compound_pathways(cache_dir)
    if (is.null(cpd_pathways) || nrow(cpd_pathways) == 0) {
        message("    Could not retrieve KEGG compound-pathway associations")
        return(NULL)
    }

    # Universe: all measured metabolites with KEGG IDs (not just those in DE table)
    if (!is.null(universe)) {
        all_cpds <- unique(universe[!is.na(universe)])
    } else {
        all_cpds <- unique(de_mapped$KEGG_ID[!is.na(de_mapped$KEGG_ID)])
    }

    # Significant compounds from DE table
    sig_cpds <- de_mapped$KEGG_ID[!is.na(de_mapped$padj) & de_mapped$padj < 0.05]
    sig_cpds <- unique(sig_cpds[!is.na(sig_cpds)])

    if (length(sig_cpds) < 3) {
        # Relax to nominal p-value
        sig_cpds <- de_mapped$KEGG_ID[!is.na(de_mapped$pvalue) & de_mapped$pvalue < 0.05]
        sig_cpds <- unique(sig_cpds[!is.na(sig_cpds)])
    }

    if (length(sig_cpds) < 2) {
        message("    Too few significant compounds for ORA (", length(sig_cpds), ")")
        return(NULL)
    }

    message("    Running compound ORA: ", length(sig_cpds), " significant / ",
            length(all_cpds), " total compounds")

    # Build pathway -> compound sets
    pathway_sets <- split(cpd_pathways$compound, cpd_pathways$pathway)

    # Build pathway name lookup
    pathway_names <- stats::setNames(cpd_pathways$name, cpd_pathways$pathway)
    pathway_names <- pathway_names[!duplicated(names(pathway_names))]

    # Filter pathways by size (intersection with measured compounds)
    use_min_gs <- max(2, min(min_gs, 3))
    N <- length(all_cpds)  # total measured compounds
    k <- length(sig_cpds)  # significant compounds

    results <- list()
    for (pw in names(pathway_sets)) {
        pw_cpds <- pathway_sets[[pw]]
        # Compounds in this pathway that are in our measured set
        pw_measured <- intersect(pw_cpds, all_cpds)
        m <- length(pw_measured)

        if (m < use_min_gs || m > max_gs) next

        # Overlap: significant compounds in this pathway
        overlap <- intersect(sig_cpds, pw_measured)
        q <- length(overlap)

        if (q == 0) next

        # Fisher's exact test (hypergeometric)
        # q-1 because phyper uses P(X > q-1) = P(X >= q)
        pval <- stats::phyper(q - 1, m, N - m, k, lower.tail = FALSE)

        results[[pw]] <- data.frame(
            pathway = if (!is.null(pathway_names[pw]) && !is.na(pathway_names[pw]))
                          pathway_names[pw] else pw,
            ID = pw,
            pvalue = pval,
            GeneRatio = paste0(q, "/", k),
            BgRatio = paste0(m, "/", N),
            setSize = q,
            compounds = paste(overlap, collapse = "/"),
            stringsAsFactors = FALSE
        )
    }

    if (length(results) == 0) {
        message("    No enriched compound pathways found")
        return(NULL)
    }

    df <- do.call(rbind, results)
    rownames(df) <- NULL

    # Multiple testing correction
    df$padj <- stats::p.adjust(df$pvalue, method = "BH")

    # Filter by p-value cutoff and sort
    df <- df[df$pvalue < pval_cutoff, ]
    if (nrow(df) == 0) {
        message("    No enriched compound pathways after p-value filtering")
        return(NULL)
    }

    df <- df[order(df$pvalue), ]

    # Class exclusion is applied here, to the finished table: the tested
    # universe above and the BH adjustment over it are untouched, so every
    # pathway a project keeps carries the p-value it would have had anyway.
    if (length(unlist(exclude_classes)) > 0) {
        df <- df[keep_kegg_pathways(df$ID, exclude = exclude_classes,
                                    cache_dir = cache_dir,
                                    label = "compound pathways"), , drop = FALSE]
        if (nrow(df) == 0) {
            message("    No compound pathways left after KEGG class exclusion")
            return(NULL)
        }
    }

    # Compound ORA is over-representation and nothing else, so the producer says
    # so rather than leaving a consumer to infer it from column shape. Metadata
    # only: the universe, the p-values, the BH adjustment, the filtering and the
    # row order above are all untouched.
    df$method <- "ora"

    message("    Found ", nrow(df), " enriched compound pathways")
    df
}


#' Get organism annotation database
get_organism_db <- function(organism) {
    db_map <- list(
        c_elegans = "org.Ce.eg.db",
        "Caenorhabditis elegans" = "org.Ce.eg.db",
        human = "org.Hs.eg.db",
        "Homo sapiens" = "org.Hs.eg.db",
        mouse = "org.Mm.eg.db",
        "Mus musculus" = "org.Mm.eg.db",
        rat = "org.Rn.eg.db",
        "Rattus norvegicus" = "org.Rn.eg.db",
        zebrafish = "org.Dr.eg.db",
        "Danio rerio" = "org.Dr.eg.db",
        drosophila = "org.Dm.eg.db",
        "Drosophila melanogaster" = "org.Dm.eg.db"
    )

    pkg <- db_map[[organism]]
    if (is.null(pkg)) {
        message("No annotation database for organism: ", organism)
        return(NULL)
    }

    if (!requireNamespace(pkg, quietly = TRUE)) {
        message("Package ", pkg, " not installed")
        return(NULL)
    }

    get(pkg, envir = asNamespace(pkg))
}


#' Get KEGG organism code
#'
#' @param organism Organism name from \code{config$global$organism}, e.g.
#'   "human" or "Homo sapiens".
#' @return KEGG organism code (e.g. "hsa"), or NULL when the organism is
#'   missing, blank or not in the lookup.
get_kegg_organism <- function(organism) {
    # A missing/blank organism (e.g. no global.organism set) must return NULL, not
    # error: list[[character(0)]] throws "attempt to select less than one element".
    if (length(organism) == 0L) return(NULL)
    # Several organisms is a config mistake rather than a missing value; returning
    # NULL would silently drop KEGG enrichment instead of pointing at the cause.
    if (length(organism) > 1L) {
        stop("global.organism must be a single organism, but it has ",
             length(organism), " values (", paste(unlist(organism), collapse = ", "),
             "). Set one name, e.g. \"human\" or \"Homo sapiens\".", call. = FALSE)
    }
    if (is.na(organism) || !nzchar(organism)) return(NULL)
    kegg_map <- list(
        c_elegans = "cel",
        "Caenorhabditis elegans" = "cel",
        human = "hsa",
        "Homo sapiens" = "hsa",
        mouse = "mmu",
        "Mus musculus" = "mmu",
        rat = "rno",
        "Rattus norvegicus" = "rno",
        zebrafish = "dre",
        "Danio rerio" = "dre",
        drosophila = "dme",
        "Drosophila melanogaster" = "dme"
    )
    kegg_map[[organism]]
}


#' Run GSEA on KEGG pathways
#'
#' Tries clusterProfiler::gseKEGG, falls back to ORA if unavailable.
run_gsea_kegg <- function(de_mapped, kegg_org, min_gs, max_gs, pval_cutoff) {
    # Try clusterProfiler GSEA first
    # Use lenient cutoff (1.0) to retrieve all results, then filter manually
    # so we can fall back from padj to pvalue when padj is too strict
    gsea_res <- tryCatch({
        de_mapped$rank_stat <- -log10(de_mapped$pvalue + 1e-300) * sign(de_mapped$log2fc)
        de_mapped <- de_mapped[order(-de_mapped$rank_stat), ]
        gene_list <- setNames(de_mapped$rank_stat, de_mapped$KEGG_ID)

        res <- clusterProfiler::gseKEGG(
            geneList = gene_list,
            organism = kegg_org,
            keyType = "kegg",
            minGSSize = min_gs,
            maxGSSize = max_gs,
            pvalueCutoff = 1.0,
            verbose = FALSE
        )
        if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
            df <- as.data.frame(res)
            out <- data.frame(
                pathway = df$Description,
                ID = df$ID,
                pvalue = df$pvalue,
                padj = df$p.adjust,
                NES = df$NES,
                setSize = df$setSize,
                stringsAsFactors = FALSE
            )
            # Filter: prefer padj, fall back to pvalue < 0.05
            padj_hits <- out[!is.na(out$padj) & out$padj < pval_cutoff, ]
            if (nrow(padj_hits) > 0) return(padj_hits)
            pval_hits <- out[!is.na(out$pvalue) & out$pvalue < 0.05, ]
            if (nrow(pval_hits) > 0) {
                message("    GSEA: padj cutoff too strict, using pvalue < 0.05 (",
                        nrow(pval_hits), " pathways)")
                return(pval_hits)
            }
        }
        NULL
    }, error = function(e) {
        message("    clusterProfiler GSEA unavailable: ", e$message)
        message("    Falling back to ORA with Fisher's exact test")
        NULL
    })

    if (!is.null(gsea_res)) return(gsea_res)

    # Fallback: run ORA using Fisher's exact test
    run_ora_kegg(de_mapped, kegg_org, min_gs, max_gs, pval_cutoff)
}


#' Run ORA on KEGG pathways
#'
#' Tries clusterProfiler::enrichKEGG, falls back to Fisher's exact test.
run_ora_kegg <- function(de_mapped, kegg_org, min_gs, max_gs, pval_cutoff) {
    # Significant genes: prefer padj < 0.05, fall back to pvalue < 0.05
    sig_genes <- de_mapped$KEGG_ID[!is.na(de_mapped$padj) & de_mapped$padj < 0.05]
    all_genes <- unique(de_mapped$KEGG_ID[!is.na(de_mapped$KEGG_ID)])

    if (length(sig_genes) < 5) {
        sig_genes <- de_mapped$KEGG_ID[!is.na(de_mapped$pvalue) & de_mapped$pvalue < 0.05]
    }
    sig_genes <- unique(sig_genes[!is.na(sig_genes)])

    if (length(sig_genes) < 5) return(NULL)

    # Try clusterProfiler first — use lenient cutoff, filter manually after
    ora_res <- tryCatch({
        res <- clusterProfiler::enrichKEGG(
            gene = sig_genes,
            universe = all_genes,
            organism = kegg_org,
            keyType = "kegg",
            minGSSize = min_gs,
            maxGSSize = max_gs,
            pvalueCutoff = 1.0
        )
        if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
            df <- as.data.frame(res)
            out <- data.frame(
                pathway = df$Description,
                ID = df$ID,
                pvalue = df$pvalue,
                padj = df$p.adjust,
                GeneRatio = df$GeneRatio,
                setSize = df$Count,
                stringsAsFactors = FALSE
            )
            # Filter: prefer padj, fall back to pvalue < 0.05
            padj_hits <- out[!is.na(out$padj) & out$padj < pval_cutoff, ]
            if (nrow(padj_hits) > 0) return(padj_hits)
            pval_hits <- out[!is.na(out$pvalue) & out$pvalue < 0.05, ]
            if (nrow(pval_hits) > 0) {
                message("    ORA: padj cutoff too strict, using pvalue < 0.05 (",
                        nrow(pval_hits), " pathways)")
                return(pval_hits)
            }
        }
        NULL
    }, error = function(e) {
        message("    clusterProfiler ORA unavailable, using Fisher's exact test")
        NULL
    })

    if (!is.null(ora_res)) {
        # Say what this is. Only the producer knows for certain -- run_gsea_kegg()
        # calls this function as its own fallback, so a table reaching a consumer
        # from that direction is ORA despite the name of the function that asked
        # for it. Metadata only: no row, p-value or adjustment is touched.
        ora_res$method <- "ora"
        return(ora_res)
    }

    # Fallback: Fisher's exact test with KEGG REST pathway-gene links; it stamps
    # its own result the same way.
    run_ora_kegg_fisher(sig_genes, all_genes, kegg_org, min_gs, max_gs, pval_cutoff)
}


#' Run ORA using Fisher's exact test with KEGG gene-pathway data
#'
#' Self-contained implementation that doesn't depend on clusterProfiler.
run_ora_kegg_fisher <- function(sig_genes, all_genes, kegg_org,
                                 min_gs, max_gs, pval_cutoff) {

    # Download gene-pathway links from KEGG
    url <- paste0("https://rest.kegg.jp/link/pathway/", kegg_org)
    lines <- tryCatch(readLines(url, warn = FALSE), error = function(e) {
        message("    KEGG pathway link download failed: ", e$message)
        NULL
    })
    if (is.null(lines) || length(lines) == 0) return(NULL)

    parts <- strsplit(lines, "\t")
    gene_pathway <- data.frame(
        gene = vapply(parts, `[`, character(1), 1),
        pathway = gsub("^path:", "", vapply(parts, `[`, character(1), 2)),
        stringsAsFactors = FALSE
    )

    # Get pathway names
    url2 <- paste0("https://rest.kegg.jp/list/pathway/", kegg_org)
    name_lines <- tryCatch(readLines(url2, warn = FALSE), error = function(e) NULL)
    pathway_names <- NULL
    if (!is.null(name_lines) && length(name_lines) > 0) {
        name_parts <- strsplit(name_lines, "\t")
        pathway_names <- stats::setNames(
            gsub(" - .*$", "", vapply(name_parts, `[`, character(1), 2)),
            gsub("^path:", "", vapply(name_parts, `[`, character(1), 1))
        )
    }

    # Build pathway sets
    pathway_sets <- split(gene_pathway$gene, gene_pathway$pathway)

    N <- length(all_genes)
    k <- length(sig_genes)

    results <- list()
    for (pw in names(pathway_sets)) {
        pw_genes <- pathway_sets[[pw]]
        pw_measured <- intersect(pw_genes, all_genes)
        m <- length(pw_measured)

        if (m < min_gs || m > max_gs) next

        overlap <- intersect(sig_genes, pw_measured)
        q <- length(overlap)
        if (q == 0) next

        pval <- stats::phyper(q - 1, m, N - m, k, lower.tail = FALSE)

        pw_name <- if (!is.null(pathway_names) && pw %in% names(pathway_names))
                       pathway_names[pw] else pw

        results[[pw]] <- data.frame(
            pathway = pw_name,
            ID = pw,
            pvalue = pval,
            GeneRatio = paste0(q, "/", k),
            setSize = q,
            stringsAsFactors = FALSE
        )
    }

    if (length(results) == 0) return(NULL)

    df <- do.call(rbind, results)
    rownames(df) <- NULL
    df$padj <- stats::p.adjust(df$pvalue, method = "BH")
    df <- df[df$pvalue < pval_cutoff, ]
    if (nrow(df) == 0) return(NULL)

    df <- df[order(df$pvalue), ]

    # This function tests one way and only one way, so it can say so. Metadata
    # only: nothing above it is re-run, re-filtered or re-ordered.
    df$method <- "ora"

    message("    Found ", nrow(df), " enriched gene pathways (Fisher's test)")
    df
}


# =============================================================================
# Cross-omics analysis: combine per-omics enrichment
# =============================================================================

#' Analyze cross-omics pathway enrichment
#'
#' @param enrichment_results Named list of enrichment data frames per omics
#' @param config Full config object
#' @param out_dir Output directory for plots
#' @return List with: combined_pathways, meta_analysis, plots
analyze_cross_omics_enrichment <- function(enrichment_results, config, out_dir = NULL) {

    if (length(enrichment_results) < 2) {
        message("Cross-omics enrichment requires >= 2 omics layers with enrichment results")
        return(NULL)
    }

    message("Analyzing cross-omics pathway enrichment...")

    omics <- names(enrichment_results)

    # Extract pathway-level results from each omics
    pathway_tables <- list()
    for (om in omics) {
        enrich_res <- enrichment_results[[om]]

        if (is.data.frame(enrich_res) && nrow(enrich_res) > 0) {
            pathway_tables[[om]] <- enrich_res
        } else if (!is.null(enrich_res$enrichment_df)) {
            pathway_tables[[om]] <- enrich_res$enrichment_df
        } else {
            warning("Cannot extract pathway table from ", om, " enrichment results")
            next
        }
    }

    if (length(pathway_tables) < 2) {
        warning("Insufficient pathway tables for cross-omics enrichment")
        return(NULL)
    }

    # Join the layers on a stable identity rather than on whichever column each
    # of them happens to use. Without this a gene layer keyed on hsa00010 and a
    # compound layer keyed on map00010 are two different pathways, and a custom
    # collection keyed on map00010 never meets a gene layer keyed on the readable
    # KEGG description at all.
    kegg_org <- resolve_kegg_org_code(config$global$organism)

    all_pathways <- lapply(pathway_tables, function(df) {
        keys <- pathway_join_key(df, kegg_org)
        keys[!is.na(keys)]
    })

    union_pathways <- unique(unlist(all_pathways))
    common_pathways <- Reduce(intersect, all_pathways)

    if (length(union_pathways) == 0) {
        message("No pathways found across omics layers")
        return(NULL)
    }

    message(sprintf("  Found %d total pathways (%d in common) across %d omics layers",
                    length(union_pathways), length(common_pathways), length(omics)))

    # The candidate universe is the union of what the layers enriched, always.
    #
    # This is not a display choice: use_pathways is what merge_pathway_pvalues()
    # assembles and what stouffer_combined_pvalues() then scores, so it decides
    # which pathways reach the cross-omics meta-analysis at all. It used to
    # collapse to the intersection whenever five or more pathways were shared,
    # which let the narrowest layer decide eligibility for every other: a layer
    # contributing a handful of rows could cut the candidate set to what it
    # happened to have in common with the rest.
    #
    # How much multi-omics evidence a pathway has is the meta-analysis's to
    # record: stouffer_combined_pvalues() counts the layers that supplied a
    # p-value for each pathway and carries that as n_omics, so a reader, a
    # ranking or a filter downstream can require two. Collapsing upstream
    # removed the rows before that count could be taken at all.
    #
    # common_pathways is still reported, and is still in the returned result;
    # it is a description of the overlap, not a gate on it.
    use_pathways <- union_pathways

    # KEGG's reference maps are pan-species, so an organism with no KEGG code of
    # its own can score well on maps of organs it does not have. Excluding those
    # classes is a reporting decision, applied to finished results: the p-values
    # and the adjustment behind them are untouched, and nothing is excluded
    # unless the config asks. The per-omics tables get the same treatment --
    # they drive the per-layer barplots and CSVs, and filtering only the merged
    # selection would leave the excluded classes visible one section away.
    excl <- .excluded_pathway_classes(config)
    if (length(excl) > 0) {
        use_pathways <- use_pathways[
            keep_kegg_pathways(use_pathways, exclude = excl, kegg_org = kegg_org,
                               label = "cross-omics pathways")]
        pathway_tables <- lapply(pathway_tables, function(df) {
            col <- if ("ID" %in% names(df)) "ID"
                   else if ("pathway" %in% names(df)) "pathway" else NULL
            if (is.null(col) || nrow(df) == 0) return(df)
            df[keep_kegg_pathways(df[[col]], exclude = excl, kegg_org = kegg_org,
                                  label = "per-omics pathways"), , drop = FALSE]
        })
    }

    # Merge pathway p-values for meta-analysis
    merged_pathways <- merge_pathway_pvalues(pathway_tables, use_pathways, omics,
                                              kegg_org = kegg_org)

    # Combine p-values using Stouffer's method. Every candidate pathway gets a
    # row: the ones a single layer enriched are kept and carry n_omics = 1,
    # rather than being dropped here. Requiring two layers is a question for
    # whoever reads or ranks the table, and n_omics is what lets them ask it.
    meta_results <- stouffer_combined_pvalues(merged_pathways)

    # The key joined the layers; the label is what a reader sees. Both are kept,
    # so the table can be traced back to the accession that produced a row.
    meta_results <- attach_pathway_display_names(meta_results, pathway_tables,
                                                  kegg_org = kegg_org)

    # Sort by combined p-value
    meta_results <- meta_results[order(meta_results$combined_pval), ]

    # Generate plots
    plots <- list()
    if (!is.null(out_dir) && nrow(meta_results) > 0) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        # 1. Cross-omics heatmap
        plots$pathway_heatmap <- file.path(out_dir, "cross_omics_pathway_heatmap.png")
        png(plots$pathway_heatmap, width = 1200, height = 900, res = 120)
        tryCatch({
            plot_cross_omics_pathway_heatmap(meta_results, omics)
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Heatmap failed:", e$message), cex = 1.2)
        })
        dev.off()

        # 2. Dot plot of top pathways per omics
        plots$dot_plot <- file.path(out_dir, "cross_omics_enrichment_dotplot.png")
        png(plots$dot_plot, width = 1200, height = 800, res = 120)
        tryCatch({
            plot_enrichment_dotplot(meta_results, omics)
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste("Dot plot failed:", e$message), cex = 1.2)
        })
        dev.off()

        # 3. ORA evidence per layer, on each layer's own adjusted p-value.
        # Built before the device is opened so that a run with no adjusted ORA
        # p-value at all -- every layer GSEA, or every layer missing the columns
        # this needs -- leaves no figure behind for the report to show.
        ora_padj <- build_ora_adjusted_p_matrix(pathway_tables, use_pathways,
                                                 omics, kegg_org = kegg_org)
        ora_png <- file.path(out_dir, "cross_omics_ora_heatmap.png")
        if (any(!is.na(ora_padj))) {
            plots$ora_heatmap <- ora_png
            png(ora_png, width = 1200, height = 900, res = 120)
            tryCatch({
                plot_cross_omics_ora_heatmap(ora_padj, pathway_tables,
                                             kegg_org = kegg_org)
            }, error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("ORA heatmap failed:", e$message), cex = 1.2)
            })
            dev.off()
        } else {
            # Delete rather than merely skip. Every other figure here is written
            # on every run and so overwrites itself; this is the only one a run
            # can decline to produce, and the report includes it on file.exists()
            # alone. Left behind, the previous run's ORA evidence would be read
            # as this one's -- a wrong figure being worse than no figure. The
            # per-contrast directories get the same treatment because this whole
            # function runs again for each contrast, with out_dir pointing there.
            if (file.exists(ora_png)) unlink(ora_png)
            message("  No adjusted ORA p-values across layers; ",
                    "skipping the cross-omics ORA heatmap")
        }

        # 4. Per-omics enrichment bar plots
        for (om in names(pathway_tables)) {
            pt <- pathway_tables[[om]]
            plot_path <- file.path(out_dir, paste0(om, "_top_pathways.png"))
            plots[[paste0(om, "_barplot")]] <- plot_path
            png(plot_path, width = 1000, height = 700, res = 120)
            tryCatch({
                plot_per_omics_barplot(pt, om)
            }, error = function(e) {
                plot.new()
                text(0.5, 0.5, paste(om, "barplot failed:", e$message), cex = 1.2)
            })
            dev.off()
        }

        message("  Cross-omics enrichment plots saved to: ", out_dir)
    }

    list(
        common_pathways = common_pathways,
        union_pathways = union_pathways,
        meta_analysis = meta_results,
        pathway_tables = pathway_tables,
        plots = plots
    )
}


# =============================================================================
# Helper functions
# =============================================================================

#' Canonical key for a contrast name
#'
#' Different omics spell the same contrast differently -- `"1.56ppm_vs_0ppm"`
#' from RNA, `"1.56ppm vs. 0ppm"` from proteomics, `"1.56ppm - 0ppm"` from
#' metabolomics, and the ORA exports drop the spaces again. All of those name
#' one biological comparison and must reduce to one key, or per-contrast work
#' silently splits a contrast in two, or pairs the wrong halves.
#'
#' This is identity, not display: keep the original string for headings and
#' filenames, and use the key only to decide what belongs with what.
#'
#' Style is what gets dropped, not content. A decimal point survives, because
#' stripping every non-alphanumeric character made `"1.56ppm vs 0ppm"` and
#' `"15.6ppm vs 0ppm"` the same key -- two different doses merged into one, or,
#' where the key picks a table, one contrast's fold changes rendered under the
#' other's name.
#'
#' @param x Character vector of contrast names.
#' @return Character vector of canonical keys, same length as \code{x}.
#' @examples
#' normalize_contrast_key(c("A vs. B", "A_vs_B", "a - b"))   # all "avsb"
#' normalize_contrast_key(c("1.56ppm vs 0ppm", "15.6ppm vs 0ppm"))  # stay apart
normalize_contrast_key <- function(x) {
    x <- tolower(trimws(x))
    # Treat " - " / "-" between groups as an alias for "vs"
    x <- gsub("\\s*-\\s*", "vs", x)
    x <- gsub("\\s*vs\\.?\\s*", "vs", x)  # "vs." / " vs " / "vs" -> "vs"
    x <- gsub("[^a-z0-9.]", "", x)        # strip separators, dots decided below
    # A dot is only content when it sits between digits; everywhere else it is
    # punctuation ("vs.", a trailing stop, make.names padding) and goes.
    x <- gsub("(?<![0-9])\\.|\\.(?![0-9])", "", x, perl = TRUE)
    x
}


#' Normalize KEGG pathway IDs to bare numeric form
#'
#' Strips organism prefixes (hsa, mmu, cel, map, etc.) to allow
#' joining gene-based and compound-based pathway results.
#' @param ids Character vector of KEGG pathway IDs
#' @return Character vector of numeric-only pathway IDs (e.g., "00010")
normalize_kegg_pathway_id <- function(ids) {
    sub("^[a-zA-Z]+", "", ids)
}


#' Is this string a KEGG pathway accession?
#'
#' Deliberately narrow. The enrichment tables mix KEGG accessions with GO terms,
#' PFAM and InterPro accessions and bare gene-set names from custom GMTs, and a
#' shape-only rule such as "three letters then five digits" would quietly rewrite
#' a custom term that happened to look like one. Only the forms this pipeline can
#' actually produce are recognised: a bare map number, the species-neutral `map`
#' and `ko` prefixes, and the organism code of the run itself.
#'
#' @param ids Character vector of pathway identifiers.
#' @param kegg_org Active KEGG organism code for the run (e.g. "hsa"), or NULL
#'   when the organism has no KEGG code.
#' @return Logical vector, one per element of \code{ids}; NA input is FALSE.
#' @examples
#' is_kegg_pathway_accession(c("00010", "map00010", "GO:0006915"), "hsa")
#' # TRUE TRUE FALSE
is_kegg_pathway_accession <- function(ids, kegg_org = NULL) {
    ids <- as.character(ids)
    pattern <- paste0("^", .kegg_accession_regex(kegg_org), "$")
    !is.na(ids) & grepl(pattern, ids)
}


#' KEGG organism code for the run, from whichever registry knows the organism
#'
#' There are two organism registries. \code{get_kegg_organism()} in this file
#' knows six species and matches the name exactly; \code{get_organism_info()} in
#' \code{R/core/11_annotation.R} knows those plus yeast, Arabidopsis, chicken,
#' pig, cow and Giardia, and matches case-insensitively after trimming.
#'
#' The join key needs the wider of the two. A run on an organism only the core
#' registry knows still gets \code{<code>#####} keys out of
#' \code{fetch_kegg_via_rest()} -- Giardia is in that registry precisely because
#' it has a KEGG code and no OrgDb -- and without its prefix those keys never
#' reduce to the bare map number, so they never meet the compound layer.
#'
#' \code{get_kegg_organism()} is deliberately left alone: it also gates which
#' organisms get full KEGG enrichment, and widening that is a different question
#' from what the join key can recognise.
#'
#' @param organism Organism name from \code{config$global$organism}.
#' @return Three-letter KEGG organism code, or NULL when neither registry has one.
resolve_kegg_org_code <- function(organism) {
    # Called first so that its error on several organisms still fires -- that is a
    # config mistake worth reporting rather than quietly resolving to nothing.
    code <- get_kegg_organism(organism)
    if (!is.null(code) && !is.na(code) && nzchar(code)) return(code)

    # get_kegg_organism() returns NULL for a missing organism as readily as for an
    # unknown one, and get_organism_info() cannot be handed a zero-length value:
    # its `%in%` test yields logical(0) and `if` errors on that. A config with no
    # organism set is a legitimate non-KEGG run, not a crash.
    if (length(organism) != 1L) return(NULL)
    organism <- as.character(organism)
    if (is.na(organism) || !nzchar(trimws(organism))) return(NULL)

    code <- get_organism_info(organism)$kegg
    if (is.null(code) || is.na(code) || !nzchar(code)) return(NULL)
    code
}


#' Regex body matching a KEGG pathway accession
#'
#' Unanchored on purpose: callers add \code{^...$} to match a bare accession, or
#' \code{^...[[:space:]]} to find one at the head of a longer key.
#'
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return Single regex string, with no anchors.
#' @keywords internal
.kegg_accession_regex <- function(kegg_org = NULL) {
    prefixes <- c("map", "ko", kegg_org)
    prefixes <- prefixes[!is.na(prefixes) & nzchar(prefixes)]
    sprintf("(%s)?[0-9]{5}", paste(prefixes, collapse = "|"))
}


#' Normalize a pathway identifier for joining, KEGG accessions only
#'
#' Wraps \code{normalize_kegg_pathway_id()} so that its prefix stripping only
#' ever reaches strings that are genuinely KEGG accessions. Applied blindly it
#' would turn "GO:0006915" into ":0006915" and strip the letters off every custom
#' gene-set name, silently merging unrelated terms.
#'
#' Two shapes are recognised: a bare accession, and an accession at the head of a
#' longer key followed by whitespace, which is how \code{fetch_kegg_via_rest()}
#' names its gene sets.
#'
#' Detection and preservation are kept apart. Whether a value is a KEGG
#' accession is decided on a trimmed copy, so " map00010 " is still recognised,
#' and that trimmed copy is what gets normalized. Everything else is returned
#' exactly as it arrived, whitespace included: an identifier this function does
#' not understand is not one to tidy up. Two layers spelling one custom set with
#' different padding will not join, which is the correct trade -- silently
#' rewriting identifiers is the failure mode worth avoiding here.
#'
#' @param ids Character vector of pathway identifiers.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return Character vector the same length as \code{ids}: KEGG accessions
#'   reduced to their bare map number, everything else byte-identical.
normalize_pathway_join_key <- function(ids, kegg_org = NULL) {
    ids <- as.character(ids)
    trimmed <- trimws(ids)
    body <- .kegg_accession_regex(kegg_org)

    exact <- is_kegg_pathway_accession(trimmed, kegg_org)
    ids[exact] <- normalize_kegg_pathway_id(trimmed[exact])

    # A gene set from fetch_kegg_via_rest() is named "<accession> <readable name>"
    # (R/core/09_enrichment.R), which is the shape the non-model KEGG fallback
    # produces -- exactly the path this join most needs to work. Take the leading
    # accession; the rest of the string is a label, and pathway_display_label()
    # still has the whole of it to show.
    labelled <- !is.na(trimmed) & !exact &
        grepl(paste0("^", body, "[[:space:]]"), trimmed)
    ids[labelled] <- normalize_kegg_pathway_id(
        sub(paste0("^(", body, ")[[:space:]].*$"), "\\1", trimmed[labelled])
    )

    ids
}


#' Drop KEGG pathways whose BRITE class a project has excluded
#'
#' A biological and reporting exclusion, not a statistical one: it is applied to
#' finished results, so the tested universe, the p-values and the BH adjustment
#' behind them are exactly what they were. Excluding a class before testing
#' would move only that layer's BH denominator while the gene-based layers had
#' already been corrected over their full family.
#'
#' Identity comes from the existing contract rather than a second KEGG rule of
#' its own: \code{normalize_pathway_join_key()} reduces every spelling this
#' pipeline produces to the bare map number, and
#' \code{is_kegg_pathway_accession()} decides which values are KEGG accessions
#' at all. Anything that is not one -- a GO term, a custom gene-set name, a
#' novel id -- has no class and is therefore kept. So is any accession the
#' hierarchy does not list.
#'
#' Kept in this file, beside those two helpers, rather than in the core
#' enrichment utilities: core is loaded before domain, so a core function
#' calling these would invert the layer order.
#'
#' @param pathway_ids Character vector of pathway identifiers.
#' @param exclude Character vector of BRITE categories or subcategories to drop.
#'   NULL or empty keeps everything, which is the default for every project.
#' @param kegg_org Active KEGG organism code, or NULL. Lets an organism-prefixed
#'   accession be recognised as this run's own.
#' @param cache_dir Passed to \code{kegg_pathway_categories()}.
#' @param label Short context word for the message naming what was dropped.
#' @param classification The resolved class table. Defaults to fetching it, and
#'   the default is lazy, so a call with nothing to exclude never reaches the
#'   network. Passing it explicitly -- including as NULL -- lets the fail-open
#'   path be exercised without depending on whether a machine has network.
#' @return Logical vector, TRUE for the pathways to keep, one per element of
#'   \code{pathway_ids}.
#' @examples
#' keep_kegg_pathways(c("map00010", "GO:0006915"), exclude = "Human Diseases")
keep_kegg_pathways <- function(pathway_ids, exclude = NULL, kegg_org = NULL,
                               cache_dir = NULL, label = "pathways",
                               classification = kegg_pathway_categories(
                                   cache_dir = cache_dir)) {
    keep <- rep(TRUE, length(pathway_ids))
    exclude <- unlist(exclude, use.names = FALSE)
    if (length(pathway_ids) == 0 || is.null(exclude) || length(exclude) == 0) {
        return(keep)
    }

    cls <- classification
    if (is.null(cls)) {
        # Fail open, and say so: a silent empty result reads as "nothing was
        # enriched" rather than "the classification could not be reached".
        message("  KEGG classification unavailable; keeping all ", label)
        return(keep)
    }

    ids <- as.character(pathway_ids)
    normalized <- normalize_pathway_join_key(ids, kegg_org)
    classifiable <- is_kegg_pathway_accession(normalized, kegg_org)

    idx <- rep(NA_integer_, length(ids))
    idx[classifiable] <- match(normalized[classifiable], cls$pathway_id)
    hit_cat <- cls$category[idx]
    hit_sub <- cls$subcategory[idx]

    drop <- (!is.na(hit_cat) & hit_cat %in% exclude) |
            (!is.na(hit_sub) & hit_sub %in% exclude)
    keep <- !drop

    if (any(drop)) {
        by_class <- ifelse(!is.na(hit_cat[drop]) & hit_cat[drop] %in% exclude,
                           hit_cat[drop], hit_sub[drop])
        counts <- table(by_class)
        message("  Excluded ", sum(drop), " of ", length(ids), " ", label,
                " by KEGG class (",
                paste(sprintf("%s: %d", names(counts), as.integer(counts)),
                      collapse = "; "), ")")
    }
    keep
}


#' Classes this run excludes, from the config
#'
#' One reader for the key, so that every section filters on the same list and a
#' class removed from one report section cannot reappear in the next.
#'
#' @param config Full config object.
#' @return Character vector of excluded classes, empty when the key is absent.
#' @keywords internal
.excluded_pathway_classes <- function(config) {
    excl <- config$modes$multiomics$enrichment$exclude_pathway_classes
    excl <- unlist(excl, use.names = FALSE)
    if (is.null(excl)) character(0) else as.character(excl)
}


#' Stable join key for each row of a pathway table
#'
#' The tables reaching cross-omics analysis do not agree on where identity
#' lives. clusterProfiler-derived layers put the readable description in
#' `pathway` and the accession in `ID`; compound ORA does the same; tables bound
#' from custom GMT collections carry the gene-set name in `pathway` and often no
#' `ID` at all. Choosing one column for a whole table therefore joins some layers
#' on accessions and others on prose, and a table bound from several collections
#' can have an `ID` column that only some of its rows populate.
#'
#' The key is resolved per row instead: the first of `ID`, `pathway`,
#' `Description` that actually carries a value, normalized only where that value
#' is a KEGG accession. This is identity, not display -- see
#' \code{pathway_display_label()} for the readable side.
#'
#' @param df Enrichment data frame for one omics layer.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return Character vector of join keys, one per row of \code{df}; NA for a row
#'   where none of the three columns carries a value.
pathway_join_key <- function(df, kegg_org = NULL) {
    key <- rep(NA_character_, nrow(df))

    for (col in c("ID", "pathway", "Description")) {
        if (!col %in% names(df)) next
        vals <- as.character(df[[col]])
        # Trimming decides whether the cell is empty; the value carried forward is
        # the original, so a custom identifier reaches the key unaltered.
        fill <- is.na(key) & !is.na(vals) & nzchar(trimws(vals))
        key[fill] <- vals[fill]
    }

    normalize_pathway_join_key(key, kegg_org)
}


#' Readable label for each row of a pathway table
#'
#' The display counterpart of \code{pathway_join_key()}: same row-wise idea,
#' opposite preference. `pathway_name` is the column \code{add_pathway_names()}
#' fills, and where it is absent the readable text is whatever sits in `pathway`
#' or `Description`.
#'
#' Unlike the join key, this one does trim: padding is worth removing from a
#' label a reader sees, and nothing is matched against it. Keep the two that way
#' round -- trimming an identifier is a silent edit, trimming a label is not.
#'
#' @param df Enrichment data frame for one omics layer.
#' @return Character vector of labels, one per row; NA where the row carries no
#'   readable text at all.
pathway_display_label <- function(df) {
    label <- rep(NA_character_, nrow(df))

    for (col in c("pathway_name", "pathway", "Description")) {
        if (!col %in% names(df)) next
        vals <- trimws(as.character(df[[col]]))
        fill <- is.na(label) & !is.na(vals) & nzchar(vals)
        label[fill] <- vals[fill]
    }

    label
}


#' Attach readable pathway names to a table keyed on join identity
#'
#' The meta-analysis table is keyed on \code{norm_id}, which is an accession and
#' not something to show a reader. The per-omics tables already carry the names
#' the layers agreed on, so the label is joined back from them rather than
#' re-derived from a gene-set collection.
#'
#' Where two layers name one key differently, a genuinely readable label wins
#' over one that is only an accession -- a layer whose gene sets were keyed on
#' `map00010` should not fix that as the display name when another layer calls it
#' "Glycolysis / Gluconeogenesis". Beyond that the first label seen wins, so the
#' result is stable across reruns; the layers are visited in the order given.
#'
#' @param meta_results Data frame with a \code{norm_id} column.
#' @param pathway_tables Named list of per-omics enrichment data frames.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return \code{meta_results} with a \code{pathway} column holding the readable
#'   label, falling back to \code{norm_id} where no layer names that key.
attach_pathway_display_names <- function(meta_results, pathway_tables,
                                          kegg_org = NULL) {
    if (is.null(meta_results) || nrow(meta_results) == 0) return(meta_results)

    keys <- character(0)
    labels <- character(0)
    for (df in pathway_tables) {
        if (!is.data.frame(df) || nrow(df) == 0) next
        keys <- c(keys, pathway_join_key(df, kegg_org))
        labels <- c(labels, pathway_display_label(df))
    }

    keep <- !is.na(keys) & !is.na(labels) & nzchar(labels)
    keys <- keys[keep]
    labels <- labels[keep]

    # Sorting readable labels ahead of bare accessions before dropping duplicates
    # leaves one label per key, and the readable one wherever a layer offered it.
    readable <- !is_kegg_pathway_accession(labels, kegg_org)
    # seq_along breaks ties explicitly rather than leaning on order() being stable.
    ord <- order(!readable, seq_along(readable))
    keys <- keys[ord]
    labels <- labels[ord]

    lookup <- stats::setNames(labels, keys)
    lookup <- lookup[!duplicated(names(lookup))]

    matched <- unname(lookup[as.character(meta_results$norm_id)])
    meta_results$pathway <- ifelse(is.na(matched),
                                   as.character(meta_results$norm_id), matched)

    # Identity stays available, but the readable column reads first.
    meta_results[, c("pathway", setdiff(names(meta_results), "pathway")),
                 drop = FALSE]
}


#' Shorten a display label to fit a plot axis
#'
#' Pulled out because the order matters and is easy to get wrong: truncate, then
#' \code{disambiguate_pathway_labels()}. Done the other way round the appended
#' key is cut straight back off and the collision it was breaking returns.
#'
#' @param labels Character vector of display labels.
#' @param max_chars Longest label to keep whole; longer ones are cut and given an
#'   ellipsis, ending up \code{max_chars} characters long.
#' @return Character vector the same length as \code{labels}.
truncate_pathway_label <- function(labels, max_chars = 50) {
    labels <- as.character(labels)
    ifelse(nchar(labels) > max_chars,
           paste0(substr(labels, 1, max_chars - 3), "..."),
           labels)
}


#' Make display labels unique without losing which pathway each one is
#'
#' Rows are keyed on \code{norm_id} but labelled with a readable name, and two
#' keys can legitimately share a name -- or share the first 50 characters of one,
#' once a caller has truncated it. A plot that positions rows by label alone then
#' stacks distinct pathways on one axis slot, and \code{pheatmap} errors outright
#' on duplicate rownames. Appending the key to just the colliding labels keeps
#' them apart and still says which pathway each row is.
#'
#' @param labels Character vector of display labels.
#' @param keys Character vector of join keys aligned to \code{labels}, or NULL
#'   when the caller has none to fall back on.
#' @return Character vector the same length as \code{labels}, unique, and
#'   unchanged wherever there was no collision.
disambiguate_pathway_labels <- function(labels, keys = NULL) {
    labels <- as.character(labels)
    if (anyDuplicated(labels) == 0) return(labels)

    if (!is.null(keys)) {
        dup <- labels %in% labels[duplicated(labels)]
        labels[dup] <- paste0(labels[dup], " (", as.character(keys)[dup], ")")
    }

    # Two rows can still collide when they share a key as well as a label; a
    # numeric suffix is ugly but beats dropping one of them off the figure.
    make.unique(labels, sep = "_")
}


#' Merge pathway p-values from multiple omics
#'
#' @param pathway_tables Named list of per-omics enrichment data frames.
#' @param target_pathways Character vector of join keys to report on, as produced
#'   by \code{pathway_join_key()}.
#' @param omics Character vector naming which layers of \code{pathway_tables} to
#'   merge, in order.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return Data frame with one row per element of \code{target_pathways}: a
#'   \code{norm_id} column and one \code{pval_<omics>} column per layer.
merge_pathway_pvalues <- function(pathway_tables, target_pathways, omics,
                                   kegg_org = NULL) {

    merged <- data.frame(norm_id = target_pathways, stringsAsFactors = FALSE)

    for (om in omics) {
        df <- pathway_tables[[om]]

        # Identity is resolved per row, so a table that mixes collections -- some
        # rows carrying an ID, some only a gene-set name -- keys each row on what
        # that row actually has.
        keys <- pathway_join_key(df, kegg_org)

        if (all(is.na(keys))) {
            warning("Cannot identify pathway column in ", om, " enrichment table")
            next
        }

        # Identify p-value column
        pval_col <- if ("pvalue" %in% names(df)) "pvalue"
                    else if ("pval" %in% names(df)) "pval"
                    else if ("p.adjust" %in% names(df)) "p.adjust"
                    else if ("padj" %in% names(df)) "padj"
                    else NULL

        if (is.null(pval_col)) {
            warning("Cannot identify p-value column in ", om, " enrichment table")
            next
        }

        # One row per key, as before: contrasts and now also KEGG prefix variants
        # of one pathway collapse to their best p-value. The result having unique
        # keys is what keeps the merge below one-to-one, so no layer can multiply
        # another's rows. Rows with no key or no p-value are dropped first, which
        # is what the formula form of aggregate() used to do via na.omit.
        # Converted only when it is not already numeric: as.numeric() on a factor
        # returns level codes, and round-tripping a double through as.character()
        # would cost precision the p-values cannot spare.
        pvals <- df[[pval_col]]
        if (!is.numeric(pvals)) {
            pvals <- suppressWarnings(as.numeric(as.character(pvals)))
        }
        usable <- !is.na(keys) & !is.na(pvals)
        if (!any(usable)) next

        df_agg <- aggregate(list(pval = pvals[usable]),
                            by = list(norm_id = keys[usable]),
                            FUN = min)
        colnames(df_agg) <- c("norm_id", paste0("pval_", om))

        # Subset to target pathways
        df_sub <- df_agg[df_agg$norm_id %in% target_pathways, , drop = FALSE]

        # Merge
        merged <- merge(merged, df_sub, by = "norm_id", all.x = TRUE)
    }

    merged
}


#' Adjusted ORA p-value per pathway and omics layer
#'
#' Assembles what the ORA figure shows, and nothing else. Deliberately separate
#' from \code{merge_pathway_pvalues()}: that one feeds
#' \code{stouffer_combined_pvalues()} and so must stay a raw-p merger with one
#' contract, while this reads each layer's own adjusted value for display. The
#' mechanics overlap; the meanings do not, and one function serving both would
#' put the meta-analysis one default away from combining adjusted p-values.
#'
#' What lands in a cell is the adjusted p-value the layer's own ORA already
#' reported. Nothing is re-adjusted here, and no cutoff is re-applied.
#'
#' Two ways a layer is refused outright, both fail-closed:
#'
#' \itemize{
#'   \item No `method` column. ORA membership cannot be confirmed, and a GSEA
#'     row in an ORA figure is exactly what this filter exists to prevent -- so
#'     the layer contributes nothing rather than being guessed at. This is the
#'     same rule \code{.kegg_hits_by_contrast()} applies for the same reason.
#'   \item No adjusted p-value column. The figure is captioned as adjusted
#'     p-values, so falling back to `pvalue` would silently show a different
#'     statistic under that caption. \code{run_ora_kegg()} does relax to the raw
#'     p-value, but it does so as the producer, with the cutoff in view; a
#'     display has no standing to repeat that decision.
#' }
#'
#' Both refusals warn, naming the layer and what was missing, because they mean
#' a column the pipeline should be producing is absent. A layer that clears both
#' and simply has no ORA row for these pathways is not an error and says
#' nothing: its column is all NA, quietly.
#'
#' Where a layer holds several contrast-level ORA results for one pathway, the
#' smallest adjusted p-value among them is taken. That is a display reduction --
#' the strongest evidence available in that layer -- and not a combined FDR:
#' nothing here controls error across contrasts.
#'
#' @param pathway_tables Named list of per-omics enrichment data frames.
#' @param target_pathways Character vector of join keys to report on, as
#'   produced by \code{pathway_join_key()}.
#' @param omics Character vector naming which layers to read, in column order.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return Numeric matrix, one row per element of \code{target_pathways} and one
#'   column per element of \code{omics}, holding adjusted ORA p-values. NA means
#'   no adjusted ORA p-value for that pathway and layer reached the table behind
#'   the figure -- not that the pathway was untested there.
#' @examples
#' tabs <- list(transcriptomics = data.frame(
#'     ID = c("map00010", "map00020"), pvalue = c(1e-4, 1e-3),
#'     padj = c(1e-3, 1e-2), method = "ora", stringsAsFactors = FALSE))
#' build_ora_adjusted_p_matrix(tabs, c("00010", "00020"), "transcriptomics")
build_ora_adjusted_p_matrix <- function(pathway_tables, target_pathways, omics,
                                         kegg_org = NULL) {

    target_pathways <- as.character(target_pathways)

    out <- matrix(NA_real_,
                  nrow = length(target_pathways), ncol = length(omics),
                  dimnames = list(target_pathways, omics))

    if (length(target_pathways) == 0 || length(omics) == 0) return(out)

    for (om in omics) {
        df <- pathway_tables[[om]]
        # A layer that produced no table at all is absent, not malformed; the
        # caller already warned about that where it mattered.
        if (!is.data.frame(df) || nrow(df) == 0) next

        missing_cols <- character(0)
        if (!"method" %in% names(df)) {
            missing_cols <- c(missing_cols, "method (ORA membership)")
        }
        if (!any(.ORA_ADJUSTED_P_COLUMNS %in% names(df))) {
            missing_cols <- c(missing_cols, "adjusted p-value (padj)")
        }
        if (length(missing_cols) > 0) {
            warning("Skipping ", om, " in the cross-omics ORA figure: no ",
                    paste(missing_cols, collapse = " and no "), " column")
            next
        }

        # NA-safe and case-insensitive, matching .kegg_hits_by_contrast(): a row
        # whose method is unknown is not an ORA row. bind_rows() NA-fills this
        # column when a layer stacks tables that do not all carry it, so the NA
        # test is load-bearing rather than defensive.
        is_ora <- !is.na(df$method) & tolower(as.character(df$method)) == "ora"
        if (!any(is_ora)) next
        df <- df[is_ora, , drop = FALSE]

        # Identity is resolved per row, so a table mixing collections keys each
        # row on what that row actually carries.
        keys <- pathway_join_key(df, kegg_org)

        padj <- .ora_adjusted_p_values(df)

        usable <- !is.na(keys) & !is.na(padj)
        if (!any(usable)) next

        best <- tapply(padj[usable], keys[usable], min)
        out[, om] <- unname(best[target_pathways])
    }

    out
}


#' Column names this pipeline uses for an adjusted ORA p-value, best first
#'
#' A short list on purpose. `padj` is what every ORA producer here writes --
#' \code{run_ora()}, \code{run_ora_kegg()}, \code{run_ora_kegg_fisher()} and
#' \code{run_compound_ora()} -- and `p.adjust` is clusterProfiler's own name for
#' the same quantity, which reaches this layer wherever an enrichResult was
#' handed over without being renamed. Nothing else belongs here: `qvalue` is a
#' different adjustment, and `pvalue` is not an adjusted one at all.
#' @keywords internal
.ORA_ADJUSTED_P_COLUMNS <- c("padj", "p.adjust")


#' The adjusted ORA p-value for each row, whichever column carries it
#'
#' Resolved per row, not per table, for the same reason
#' \code{pathway_join_key()} resolves identity per row: the tables arriving here
#' are bound from heterogeneous sub-results. \code{extract_enrichment_df()} uses
#' \code{dplyr::bind_rows()}, which NA-fills, so one layer can hold fgsea rows
#' carrying `padj` beside clusterProfiler ORA rows carrying only `p.adjust`.
#' Both columns then exist, and picking one for the whole table reads NA for
#' every row that used the other -- silently, because the column check passed.
#'
#' Both names mean the same quantity, so this is choosing where to read it, not
#' choosing between statistics.
#'
#' @param df Enrichment data frame for one omics layer.
#' @return Numeric vector, one per row of \code{df}; NA where no accepted column
#'   carries a value for that row.
#' @examples
#' mixed <- data.frame(padj = c(0.01, NA), p.adjust = c(NA, 0.02))
#' .ora_adjusted_p_values(mixed)   # 0.01 0.02
.ora_adjusted_p_values <- function(df) {
    vals <- rep(NA_real_, nrow(df))

    for (col in .ORA_ADJUSTED_P_COLUMNS) {
        if (!col %in% names(df)) next
        # Converted only when it is not already numeric: as.numeric() on a factor
        # returns level codes, and a round trip through as.character() would cost
        # precision these p-values cannot spare.
        v <- df[[col]]
        if (!is.numeric(v)) v <- suppressWarnings(as.numeric(as.character(v)))
        fill <- is.na(vals) & !is.na(v)
        vals[fill] <- v[fill]
    }

    vals
}


#' Combine p-values across omics using Stouffer's method
#'
#' Stouffer's Z, not Fisher's — the two give different combined p-values and
#' the report used to name the wrong one.
stouffer_combined_pvalues <- function(merged_pathways) {

    pval_cols <- grep("^pval_", names(merged_pathways), value = TRUE)

    if (length(pval_cols) < 2) {
        # Single omics: just use the one column
        if (length(pval_cols) == 1) {
            merged_pathways$combined_pval <- merged_pathways[[pval_cols[1]]]
            merged_pathways$combined_padj <- p.adjust(merged_pathways$combined_pval, method = "BH")
            merged_pathways$n_omics <- as.integer(!is.na(merged_pathways[[pval_cols[1]]]))
            return(merged_pathways)
        }
        stop("Need at least 1 p-value column")
    }

    pval_matrix <- as.matrix(merged_pathways[, pval_cols])

    # Stouffer's method (Loughin 2004, PMC3653960):
    #   z_i = Φ^{-1}(p_i)  — small p → large negative z
    #   Z_S = Σ z_i / sqrt(k)  ~ N(0,1) under H0
    #   combined p = Φ(Z_S)    — left tail
    combined_pvals <- apply(pval_matrix, 1, function(pvals) {
        pvals <- pvals[!is.na(pvals) & pvals > 0]
        if (length(pvals) == 0) return(NA)

        # Clamp to avoid Inf from qnorm(0) or qnorm(1)
        pvals <- pmax(pvals, .Machine$double.xmin)
        pvals <- pmin(pvals, 1 - .Machine$double.eps)

        z_scores <- qnorm(pvals)
        z_combined <- sum(z_scores) / sqrt(length(z_scores))
        pnorm(z_combined)
    })

    n_omics <- rowSums(!is.na(pval_matrix))

    merged_pathways$combined_pval <- combined_pvals
    merged_pathways$combined_padj <- p.adjust(combined_pvals, method = "BH")
    merged_pathways$n_omics <- n_omics

    merged_pathways
}


# =============================================================================
# Plotting functions
# =============================================================================

#' A colour range that will not collapse on a degenerate matrix
#'
#' Both heatmap paths cut a value range into one interval per colour, and both
#' break when that range has zero width: \code{pheatmap} derives its breaks from
#' the minimum and maximum, so an identical pair yields duplicate breaks and
#' \code{cut()} refuses them, and \code{image()} wants a \code{zlim} whose ends
#' differ. Zero width is not exotic here -- one cell with evidence, or several
#' cells that happen to agree, is an ordinary sparse result, especially in a
#' per-contrast view. The figures' tryCatch would then paint a "failed"
#' placeholder over a perfectly drawable result.
#'
#' Widening a degenerate range is a display accommodation: the value shown is
#' unchanged, it simply lands mid-scale instead of at an end, which is honest
#' for a matrix that carries no contrast to show.
#'
#' @param m Numeric matrix, possibly holding NA.
#' @return Numeric pair, low then high, with the low strictly below the high.
.nondegenerate_range <- function(m) {
    rng <- suppressWarnings(range(m, na.rm = TRUE))
    # All NA gives c(Inf, -Inf); nothing is drawn from it, but it must not reach
    # a breaks calculation either.
    if (!all(is.finite(rng))) return(c(0, 1))
    if (rng[1] == rng[2]) return(rng + c(-0.5, 0.5))
    rng
}


#' Draw a heatmap too small for stats::heatmap()
#'
#' \code{heatmap()} rejects fewer than two rows or two columns even with both
#' dendrograms disabled, and one surviving pathway is a normal outcome in a
#' per-contrast view. It draws through \code{image()} anyway, so the degenerate
#' case goes straight there rather than padding the matrix -- a fabricated row
#' would put a pathway on the figure that the data does not have.
#'
#' @param m Numeric matrix to draw, with dimnames for the axes.
#' @param main Plot title.
#' @param col Colour vector.
#' @param zlim Value range, from \code{.nondegenerate_range()}.
#' @return Invisibly NULL; called for the plot it draws.
.draw_small_heatmap <- function(m, main, col, zlim) {
    # image() takes z indexed [x, y] -- columns then rows -- so the matrix is
    # transposed and the row labels follow it.
    #
    # Cell edges, not cell centres. Given centres, image() infers the edges from
    # the spacing between them, which needs at least two: a length-one x or y is
    # rejected with "dimensions of z are not length(x)(-1) times length(y)(-1)".
    # That is precisely the shape this function exists for, so the edges are
    # given outright and every size from 1x1 upward draws.
    graphics::image(x = seq(0.5, ncol(m) + 0.5, by = 1),
                    y = seq(0.5, nrow(m) + 0.5, by = 1), z = t(m),
                    col = col, zlim = zlim, axes = FALSE,
                    xlab = "", ylab = "", main = main)
    graphics::axis(1, at = seq_len(ncol(m)), labels = colnames(m),
                   las = 2, tick = FALSE)
    graphics::axis(2, at = seq_len(nrow(m)), labels = rownames(m),
                   las = 1, tick = FALSE)
    invisible(NULL)
}


#' Choose the rows a cross-omics figure shows
#'
#' The figures exist to show where the layers agree, and ordering by combined
#' p-value alone does not do that: the layer with the largest gene-set
#' collection contributes many single-layer pathways with very small p-values,
#' and they fill every slot. The pathways several layers support -- the point of
#' the figure -- rank below them and never appear.
#'
#' Rows are therefore ordered by how many layers contributed an enrichment
#' p-value for the pathway first, and by the combined p-value within that. Note
#' what that count is and is not: `n_omics` counts the layers whose enrichment
#' table carried a p-value for this pathway. Some of those tables reach here
#' already filtered, so a missing p-value can mean the pathway was never
#' testable in that layer, or that it was tested and did not survive into the
#' layer's result table. Nothing downstream can tell those apart.
#'
#' This is selection for display only. The meta-analysis table keeps its own
#' combined_pval ordering, and no p-value, adjustment or membership is
#' touched.
#'
#' The two ranking columns are named rather than fixed, because the ordering is
#' the reusable part and the quantities are not. The defaults are the
#' meta-analysis pair and every existing caller keeps them. The ORA figure,
#' whose evidence is a different table with different missingness, passes its
#' own pair: it must not be ranked on raw-p meta-analysis columns.
#'
#' @param meta_results Meta-analysis table, as
#'   \code{stouffer_combined_pvalues()} returns it, or any table carrying
#'   \code{count_col} and \code{score_col}.
#' @param top_n Number of rows to keep.
#' @param count_col Column holding the number of contributing layers; more is
#'   better. Absent, every row counts as one and the order is left alone.
#' @param score_col Column holding the score that breaks ties within a count;
#'   smaller is better. Absent, nothing breaks them but the incoming order.
#' @return The selected rows of \code{meta_results}, in display order.
#' @examples
#' meta <- data.frame(norm_id = c("00010", "00020"),
#'                    n_omics = c(1L, 2L), combined_pval = c(1e-9, 1e-3))
#' select_multi_omics_pathways(meta, top_n = 2)$norm_id   # "00020" first
select_multi_omics_pathways <- function(meta_results, top_n = 30,
                                         count_col = "n_omics",
                                         score_col = "combined_pval") {
    if (is.null(meta_results) || nrow(meta_results) == 0) return(meta_results)

    counts <- if (count_col %in% names(meta_results)) {
        as.numeric(meta_results[[count_col]])
    } else {
        # Nothing to rank on: leave the caller's order alone rather than invent
        # a preference between rows that carry no contributing-layer count.
        rep(1, nrow(meta_results))
    }
    scores <- if (score_col %in% names(meta_results)) {
        as.numeric(meta_results[[score_col]])
    } else {
        rep(NA_real_, nrow(meta_results))
    }

    # Ties on both keys fall back to the incoming order, which is itself sorted
    # by combined p-value, so the selection is reproducible run to run.
    ord <- order(-counts, scores, seq_len(nrow(meta_results)), na.last = TRUE)
    meta_results[utils::head(ord, min(top_n, nrow(meta_results))), , drop = FALSE]
}


#' Plot cross-omics pathway heatmap
plot_cross_omics_pathway_heatmap <- function(meta_results, omics, top_n = 30) {

    # Rows that several layers support lead; see select_multi_omics_pathways().
    top_pathways <- select_multi_omics_pathways(meta_results, top_n)

    pval_cols <- grep("^pval_", names(top_pathways), value = TRUE)
    pval_matrix <- as.matrix(top_pathways[, pval_cols, drop = FALSE])

    # Truncate first, then disambiguate -- see truncate_pathway_label().
    pathway_labels <- truncate_pathway_label(top_pathways$pathway, 50)
    pathway_labels <- disambiguate_pathway_labels(pathway_labels,
                                                   top_pathways$norm_id)
    rownames(pval_matrix) <- pathway_labels

    # Transform to -log10(p)
    log_pval_matrix <- -log10(pval_matrix + 1e-300)
    # Cap at 10 for display. The NA test is not decoration: a logical subscript
    # carrying NA is an error in `[<-`, and a missing layer p-value is NA -- the
    # normal case now that the candidate universe is the union of the layers.
    capped <- !is.na(log_pval_matrix) & log_pval_matrix > 10
    log_pval_matrix[capped] <- 10
    colnames(log_pval_matrix) <- gsub("^pval_", "", colnames(log_pval_matrix))

    # NA stays NA. "No enrichment p-value for this pathway in this layer" is
    # not "a p-value close to 1", and flattening the two to 0 rendered them the
    # same white -- which also left na_col below as dead configuration.

    # Breaks are passed rather than left to pheatmap, which derives them from the
    # minimum and maximum and so produces duplicates when every value agrees --
    # see .nondegenerate_range(). One interval per colour, hence 50 + 1.
    zl <- .nondegenerate_range(log_pval_matrix)

    # Heatmap
    if (requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(log_pval_matrix,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           main = "Cross-Omics Pathway Enrichment (-log10 p-value)",
                           color = colorRampPalette(c("white", "gold", "orange", "red"))(50),
                           breaks = seq(zl[1], zl[2], length.out = 51),
                           fontsize_row = 7, fontsize_col = 10,
                           angle_col = 45,
                           na_col = "grey90",
                           border_color = "grey80")
    } else {
        # Rowv = NA as well as Colv = NA, to match the pheatmap path above,
        # which disables both dendrograms. Leaving row clustering on was also a
        # failure waiting for the right input: two rows whose missing layer
        # p-values do not overlap share no observed cell, so dist() returns NA
        # between them and hclust() stops on it. Missing values stay missing --
        # nothing here fills them to make clustering possible.
        # image(), which heatmap() draws through, simply does not paint an NA
        # cell -- so the device background shows through, and on white that is
        # the same white as the low end of the scale. The legend promises grey
        # for a missing p-value, so the background is grey while this draws.
        # The value stays NA; only what shows behind it changes.
        withr::with_par(list(bg = "grey90"), {
            base_cols <- colorRampPalette(c("white", "orange", "red"))(50)
            if (nrow(log_pval_matrix) >= 2 && ncol(log_pval_matrix) >= 2) {
                heatmap(log_pval_matrix, scale = "none", Rowv = NA, Colv = NA,
                        main = "Cross-Omics Pathway Enrichment",
                        col = base_cols, zlim = zl)
            } else {
                # heatmap() refuses this shape outright; see .draw_small_heatmap().
                .draw_small_heatmap(log_pval_matrix,
                                    main = "Cross-Omics Pathway Enrichment",
                                    col = base_cols, zlim = zl)
            }
        })
    }
}


#' Plot the ORA evidence each omics layer holds for a pathway
#'
#' ORA is the one test all the layers run, so this is the figure that can show
#' metabolomics alongside the gene layers. Each cell is -log10 of the adjusted
#' p-value that layer's own ORA reported, as
#' \code{build_ora_adjusted_p_matrix()} assembled it.
#'
#' Read it down a column. Across columns the layers differ in gene-set universe,
#' coverage and upstream filtering, so a darker cell in one is not stronger
#' biology than a lighter cell in another; the figure shows what evidence each
#' layer holds, not a comparison of magnitudes between them.
#'
#' Rows with no adjusted ORA p-value anywhere are dropped before selection --
#' they would otherwise take slots from rows that have something to show. What
#' remains is ordered by how many layers contribute a value and then by the
#' smallest of them, the same shape the meta-analysis figures use, but computed
#' from this matrix rather than from the raw-p meta-analysis columns.
#'
#' @param padj_matrix Adjusted ORA p-values, as
#'   \code{build_ora_adjusted_p_matrix()} returns them.
#' @param pathway_tables Named list of per-omics enrichment data frames, used
#'   only to label rows with the names the layers agreed on. NULL labels rows
#'   with their join keys.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @param top_n Number of pathways to show.
#' @return Invisibly, the -log10 matrix that was drawn, in display order and
#'   with display labels for row names -- or NULL when there was nothing to
#'   draw. Returned so that which rows the figure leads with can be checked
#'   without reading pixels, or guessing which of the two drawing branches a
#'   machine took.
plot_cross_omics_ora_heatmap <- function(padj_matrix, pathway_tables = NULL,
                                          kegg_org = NULL, top_n = 30) {

    informative <- rowSums(!is.na(padj_matrix)) > 0
    padj_matrix <- padj_matrix[informative, , drop = FALSE]

    if (nrow(padj_matrix) == 0) {
        plot.new()
        text(0.5, 0.5, "No adjusted ORA p-values available")
        return(invisible(NULL))
    }

    # Ranked on this figure's own evidence. n_ora_layers counts the layers that
    # contributed an adjusted ORA p-value for the pathway, which is not the same
    # as the layers that could have tested it: some tables arrive already
    # filtered, and nothing here can tell a pathway that was never testable in a
    # layer from one that was tested and did not survive into its table.
    # unname() both: a named vector would hand data.frame() its row names, and
    # two layers can legitimately label one key the same way.
    # Every remaining row has at least one value, so min(na.rm = TRUE) is never
    # asked for the minimum of nothing.
    ranking <- data.frame(
        row           = seq_len(nrow(padj_matrix)),
        n_ora_layers  = unname(rowSums(!is.na(padj_matrix))),
        best_ora_padj = unname(apply(padj_matrix, 1, min, na.rm = TRUE)),
        stringsAsFactors = FALSE
    )
    selected <- select_multi_omics_pathways(ranking, top_n,
                                            count_col = "n_ora_layers",
                                            score_col = "best_ora_padj")
    padj_matrix <- padj_matrix[selected$row, , drop = FALSE]

    # Truncate first, then disambiguate -- see truncate_pathway_label().
    keys <- rownames(padj_matrix)
    labels <- keys
    if (!is.null(pathway_tables)) {
        labelled <- attach_pathway_display_names(
            data.frame(norm_id = keys, stringsAsFactors = FALSE),
            pathway_tables, kegg_org = kegg_org)
        labels <- labelled$pathway
    }
    labels <- truncate_pathway_label(labels, 50)
    rownames(padj_matrix) <- disambiguate_pathway_labels(labels, keys)

    log_padj_matrix <- -log10(padj_matrix + 1e-300)
    # Cap at 10 for display. The NA test is not decoration: a logical subscript
    # carrying NA is an error in `[<-`, and a layer holding no adjusted ORA
    # p-value for a pathway is the normal case here.
    capped <- !is.na(log_padj_matrix) & log_padj_matrix > 10
    log_padj_matrix[capped] <- 10

    # NA stays NA, as in plot_cross_omics_pathway_heatmap(): "no adjusted ORA
    # p-value reached this table" is not "an adjusted p-value close to 1".

    # Pink through magenta, deliberately unlike the white-gold-orange-red of the
    # meta-analysis heatmap: the two figures sit close together in the report and
    # are scored on different quantities.
    ora_palette <- c("#fff7f3", "#fcc5c0", "#f768a1", "#ae017e", "#7a0177")

    # A sparse result is the common case for this figure -- one layer, one
    # pathway, one cell -- and a zero-width range breaks both drawing paths.
    # See .nondegenerate_range(); one interval per colour, hence 50 + 1.
    zl <- .nondegenerate_range(log_padj_matrix)

    if (requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(log_padj_matrix,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           main = "Cross-Omics ORA Evidence (-log10 adjusted p-value)",
                           color = colorRampPalette(ora_palette)(50),
                           breaks = seq(zl[1], zl[2], length.out = 51),
                           fontsize_row = 7, fontsize_col = 10,
                           angle_col = 45,
                           na_col = "grey90",
                           border_color = "grey80")
    } else {
        # Both dendrograms off, matching the pheatmap path: two rows whose
        # missing layers do not overlap share no observed cell, so dist()
        # returns NA between them and hclust() stops on it. Nothing is filled in
        # to make clustering possible.
        # image(), which heatmap() draws through, does not paint an NA cell at
        # all, so the device background shows through it; on white that is the
        # same white as the low end of the scale, and the legend promises grey.
        withr::with_par(list(bg = "grey90"), {
            base_cols <- colorRampPalette(ora_palette)(50)
            if (nrow(log_padj_matrix) >= 2 && ncol(log_padj_matrix) >= 2) {
                heatmap(log_padj_matrix, scale = "none", Rowv = NA, Colv = NA,
                        main = "Cross-Omics ORA Evidence",
                        col = base_cols, zlim = zl)
            } else {
                # heatmap() refuses this shape outright; see .draw_small_heatmap().
                .draw_small_heatmap(log_padj_matrix,
                                    main = "Cross-Omics ORA Evidence",
                                    col = base_cols, zlim = zl)
            }
        })
    }

    invisible(log_padj_matrix)
}


#' Plot enrichment dot plot
plot_enrichment_dotplot <- function(meta_results, omics, top_n = 20) {

    top <- select_multi_omics_pathways(meta_results, top_n)

    pval_cols <- grep("^pval_", names(top), value = TRUE)

    # The axis is built from the label, so two keys sharing a name would land on
    # one position and hide each other. Truncate before disambiguating, not after:
    # the key suffix sits at the end of the string, so truncating second cuts it
    # straight back off and rebuilds the collision -- and the duplicate then
    # reaches factor(levels = ) below, which rejects duplicated levels outright.
    top$pathway <- truncate_pathway_label(top$pathway, 45)
    top$pathway <- disambiguate_pathway_labels(top$pathway, top$norm_id)

    # Build long-format data
    plot_data <- list()
    for (pc in pval_cols) {
        om_name <- sub("^pval_", "", pc)
        df <- data.frame(
            pathway = top$pathway,
            omics = om_name,
            neg_log10_p = -log10(top[[pc]] + 1e-300),
            stringsAsFactors = FALSE
        )
        df$neg_log10_p[is.na(top[[pc]])] <- NA
        plot_data[[om_name]] <- df
    }
    plot_df <- do.call(rbind, plot_data)
    plot_df <- plot_df[!is.na(plot_df$neg_log10_p), ]

    if (nrow(plot_df) == 0) {
        plot.new()
        text(0.5, 0.5, "No data for dot plot")
        return(invisible(NULL))
    }

    # Cap for display
    plot_df$neg_log10_p <- pmin(plot_df$neg_log10_p, 10)

    # Labels were truncated and made unique before plot_df was built, so the
    # levels are already the strings on the axis and are guaranteed distinct.
    pathway_order <- rev(unique(top$pathway))
    plot_df$pathway <- factor(plot_df$pathway, levels = pathway_order)

    omics_colors <- c(
        transcriptomics = "#E41A1C",
        proteomics = "#377EB8",
        metabolomics = "#4DAF4A"
    )
    available_colors <- omics_colors[intersect(names(omics_colors), unique(plot_df$omics))]

    # Plot
    par(mar = c(5, 15, 3, 2))
    plot(NULL, xlim = c(0, max(plot_df$neg_log10_p, na.rm = TRUE) * 1.1),
         ylim = c(0.5, length(pathway_order) + 0.5),
         xlab = "-log10(p-value)", ylab = "",
         yaxt = "n", main = "Cross-Omics Pathway Enrichment")
    axis(2, at = seq_along(pathway_order), labels = pathway_order, las = 1, cex.axis = 0.7)
    abline(h = seq_along(pathway_order), col = "grey90", lty = 2)

    omics_list <- unique(plot_df$omics)
    offsets <- seq(-0.15, 0.15, length.out = length(omics_list))

    for (i in seq_along(omics_list)) {
        om <- omics_list[i]
        sub_df <- plot_df[plot_df$omics == om, ]
        y_pos <- as.numeric(sub_df$pathway) + offsets[i]
        col <- if (om %in% names(available_colors)) available_colors[om] else i + 1
        points(sub_df$neg_log10_p, y_pos, pch = 19,
               cex = 1.2, col = col)
    }

    legend("bottomright", legend = omics_list, col = sapply(omics_list, function(om) {
        if (om %in% names(available_colors)) available_colors[om] else which(omics_list == om) + 1
    }), pch = 19, cex = 0.8, bg = "white")
}


#' Plot per-omics enrichment barplot
plot_per_omics_barplot <- function(pathway_table, omics_name, top_n = 15) {

    pval_col <- if ("pvalue" %in% names(pathway_table)) "pvalue"
                else if ("pval" %in% names(pathway_table)) "pval"
                else if ("padj" %in% names(pathway_table)) "padj"
                else NULL

    pathway_col <- if ("pathway" %in% names(pathway_table)) "pathway"
                   else if ("Description" %in% names(pathway_table)) "Description"
                   else if ("ID" %in% names(pathway_table)) "ID"
                   else NULL

    if (is.null(pval_col) || is.null(pathway_col)) {
        plot.new()
        text(0.5, 0.5, paste("Cannot identify columns for", omics_name))
        return(invisible(NULL))
    }

    # Aggregate across contrasts: take min p-value per pathway
    agg <- aggregate(
        stats::as.formula(paste(pval_col, "~", pathway_col)),
        data = pathway_table,
        FUN = min
    )
    colnames(agg) <- c("pathway", "pvalue")
    agg <- agg[order(agg$pvalue), ]
    agg <- agg[seq_len(min(top_n, nrow(agg))), ]

    # Truncate names
    agg$label <- ifelse(nchar(agg$pathway) > 45,
                        paste0(substr(agg$pathway, 1, 42), "..."),
                        agg$pathway)

    neg_log_p <- -log10(agg$pvalue + 1e-300)
    neg_log_p <- pmin(neg_log_p, 15)

    omics_colors <- c(
        transcriptomics = "#E41A1C",
        proteomics = "#377EB8",
        metabolomics = "#4DAF4A"
    )
    bar_col <- if (omics_name %in% names(omics_colors)) omics_colors[omics_name] else "steelblue"

    par(mar = c(5, 15, 3, 2))
    barplot(rev(neg_log_p), horiz = TRUE, names.arg = rev(agg$label),
            las = 1, cex.names = 0.7, col = bar_col,
            xlab = "-log10(p-value)",
            main = paste("Top Pathways -", omics_name))
    abline(v = -log10(0.05), col = "red", lty = 2)
}


#' Write cross-omics enrichment results
write_cross_omics_enrichment <- function(enrichment_res, out_dir) {

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    if (!is.null(enrichment_res$meta_analysis)) {
        write.csv(enrichment_res$meta_analysis,
                  file.path(out_dir, "cross_omics_pathways_meta_analysis.csv"),
                  row.names = FALSE)
    }

    # Write per-omics tables
    if (!is.null(enrichment_res$pathway_tables)) {
        for (om in names(enrichment_res$pathway_tables)) {
            write.csv(enrichment_res$pathway_tables[[om]],
                      file.path(out_dir, paste0(om, "_enriched_pathways.csv")),
                      row.names = FALSE)
        }
    }

    message("Cross-omics enrichment results written to: ", out_dir)
    invisible(NULL)
}


# =============================================================================
# Loadings-based enrichment (DIABLO / MOFA2 top features)
# =============================================================================

#' Run geneset enrichment on integration loadings
#'
#' Takes top features from DIABLO loadings or MOFA2 weights and runs
#' ORA enrichment (KEGG) for each component/factor per omics view.
#'
#' @param integration_res Output from mod_multiomics_integration()
#' @param harmonization_res Output from mod_multiomics_harmonization()
#' @param config Full config object
#' @param out_dir Output directory for results
#' @param top_n Number of top features per component/factor to use (default 50)
#' @return List with diablo and mofa enrichment results
run_loadings_enrichment <- function(integration_res, harmonization_res,
                                     config, out_dir, top_n = 50) {

    message("\n=== Loadings-based Geneset Enrichment ===\n")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    organism <- config$global$organism
    kegg_org <- get_kegg_organism(organism)
    org_db <- get_organism_db(organism)

    results <- list()

    # --- DIABLO loadings enrichment ---
    if (!is.null(integration_res$diablo_results)) {
        message("Running enrichment on DIABLO loadings...")
        diablo_dir <- file.path(out_dir, "diablo_loadings")
        dir.create(diablo_dir, showWarnings = FALSE)

        results$diablo <- tryCatch(
            run_diablo_loadings_enrichment(
                exclude_classes = .excluded_pathway_classes(config),
                diablo_results = integration_res$diablo_results,
                harmonization_res = harmonization_res,
                organism = organism,
                kegg_org = kegg_org,
                org_db = org_db,
                out_dir = diablo_dir,
                top_n = top_n
            ),
            error = function(e) {
                message("  DIABLO loadings enrichment failed: ", conditionMessage(e))
                NULL
            }
        )
    }

    # --- MOFA2 weights enrichment ---
    if (!is.null(integration_res$mofa_results)) {
        message("Running enrichment on MOFA2 weights...")
        mofa_dir <- file.path(out_dir, "mofa_loadings")
        dir.create(mofa_dir, showWarnings = FALSE)

        results$mofa <- tryCatch(
            run_mofa_weights_enrichment(
                exclude_classes = .excluded_pathway_classes(config),
                mofa_results = integration_res$mofa_results,
                harmonization_res = harmonization_res,
                organism = organism,
                kegg_org = kegg_org,
                org_db = org_db,
                out_dir = mofa_dir,
                top_n = top_n
            ),
            error = function(e) {
                message("  MOFA2 weights enrichment failed: ", conditionMessage(e))
                NULL
            }
        )
    }

    message("Loadings enrichment complete")
    results
}


#' Run enrichment on DIABLO top loadings per component
run_diablo_loadings_enrichment <- function(diablo_results, harmonization_res,
                                            organism, kegg_org, org_db,
                                            out_dir, top_n = 50, exclude_classes = NULL) {

    top_features <- diablo_results$top_features
    if (is.null(top_features) || length(top_features) == 0) return(NULL)

    all_results <- list()

    for (om in names(top_features)) {
        if (om == "Y") next  # Skip outcome
        feat_df <- top_features[[om]]
        if (is.null(feat_df) || nrow(feat_df) == 0) next

        # Metabolomics: compound-based ORA instead of gene-based
        if (om == "metabolomics") {
            message("  Running metabolomics loadings enrichment (compound ORA)")
            components <- unique(feat_df$component)
            for (comp in components) {
                comp_feats <- feat_df[feat_df$component == comp, ]
                comp_feats <- comp_feats[order(-comp_feats$abs_loading), ]
                top_feat_ids <- head(comp_feats$feature, top_n)

                label <- paste0("DIABLO_metabolomics_", comp)
                message("  ", label, ": ", length(top_feat_ids), " features")

                metab_enrich <- run_metabolite_loadings_ora(
                    top_feat_ids, harmonization_res, out_dir, label,
                    exclude_classes = exclude_classes
                )
                if (!is.null(metab_enrich) && nrow(metab_enrich) > 0) {
                    metab_enrich$method <- "DIABLO"
                    metab_enrich$omics <- "metabolomics"
                    metab_enrich$component <- comp
                    all_results[[label]] <- metab_enrich
                }
            }
            next
        }

        components <- unique(feat_df$component)
        for (comp in components) {
            comp_feats <- feat_df[feat_df$component == comp, ]
            comp_feats <- comp_feats[order(-comp_feats$abs_loading), ]
            top_feat_ids <- head(comp_feats$feature, top_n)

            label <- paste0("DIABLO_", om, "_", comp)
            message("  ", label, ": ", length(top_feat_ids), " features")

            enrich_df <- enrich_feature_list(
                feature_ids = top_feat_ids,
                omics_type = om,
                harmonization_res = harmonization_res,
                organism = organism,
                kegg_org = kegg_org,
                org_db = org_db
            )

            if (!is.null(enrich_df) && nrow(enrich_df) > 0) {
                enrich_df$method <- "DIABLO"
                enrich_df$omics <- om
                enrich_df$component <- comp
                all_results[[label]] <- enrich_df

                write.csv(enrich_df,
                          file.path(out_dir, paste0(label, "_enrichment.csv")),
                          row.names = FALSE)

                # Barplot
                plot_loadings_enrichment_barplot(
                    enrich_df, label,
                    file.path(out_dir, paste0(label, "_enrichment.png"))
                )
            }
        }
    }

    if (length(all_results) == 0) return(NULL)
    combined <- .rbind_fill(all_results)
    rownames(combined) <- NULL
    write.csv(combined, file.path(out_dir, "diablo_loadings_enrichment_all.csv"),
              row.names = FALSE)
    combined
}


#' Row-bind data frames with differing column sets (union-of-columns).
#' Gene-based ORA and compound ORA return different columns, so a plain
#' do.call(rbind, ...) fails. This pads missing columns with NA.
.rbind_fill <- function(dfs) {
    dfs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, dfs)
    if (length(dfs) == 0) return(NULL)
    if (length(dfs) == 1) return(dfs[[1]])
    all_cols <- unique(unlist(lapply(dfs, colnames)))
    aligned <- lapply(dfs, function(df) {
        missing <- setdiff(all_cols, colnames(df))
        for (m in missing) df[[m]] <- NA
        df[, all_cols, drop = FALSE]
    })
    do.call(rbind, aligned)
}


#' Run enrichment on MOFA2 top weights per factor
run_mofa_weights_enrichment <- function(mofa_results, harmonization_res,
                                         organism, kegg_org, org_db,
                                         out_dir, top_n = 50, exclude_classes = NULL) {

    weights <- mofa_results$weights
    if (is.null(weights) || length(weights) == 0) return(NULL)

    all_results <- list()

    for (view in names(weights)) {
        w <- weights[[view]]
        if (view == "metabolomics") {
            message("  Running metabolomics weights enrichment (compound ORA)")
            n_factors <- min(ncol(w), 3)
            for (k in seq_len(n_factors)) {
                factor_name <- colnames(w)[k]
                loadings <- w[, k]
                ord <- order(abs(loadings), decreasing = TRUE)
                top_feat_ids <- rownames(w)[head(ord, top_n)]

                label <- paste0("MOFA_metabolomics_", factor_name)
                message("  ", label, ": ", length(top_feat_ids), " features")

                metab_enrich <- run_metabolite_loadings_ora(
                    top_feat_ids, harmonization_res, out_dir, label,
                    exclude_classes = exclude_classes
                )
                if (!is.null(metab_enrich) && nrow(metab_enrich) > 0) {
                    metab_enrich$method <- "MOFA2"
                    metab_enrich$view <- "metabolomics"
                    metab_enrich$factor <- factor_name
                    all_results[[label]] <- metab_enrich
                }
            }
            next
        }

        n_factors <- min(ncol(w), 3)  # Top 3 factors

        for (k in seq_len(n_factors)) {
            factor_name <- colnames(w)[k]
            loadings <- w[, k]
            ord <- order(abs(loadings), decreasing = TRUE)
            top_feat_ids <- rownames(w)[head(ord, top_n)]

            label <- paste0("MOFA_", view, "_", factor_name)
            message("  ", label, ": ", length(top_feat_ids), " features")

            enrich_df <- enrich_feature_list(
                feature_ids = top_feat_ids,
                omics_type = view,
                harmonization_res = harmonization_res,
                organism = organism,
                kegg_org = kegg_org,
                org_db = org_db
            )

            if (!is.null(enrich_df) && nrow(enrich_df) > 0) {
                enrich_df$method <- "MOFA2"
                enrich_df$view <- view
                enrich_df$factor <- factor_name
                all_results[[label]] <- enrich_df

                write.csv(enrich_df,
                          file.path(out_dir, paste0(label, "_enrichment.csv")),
                          row.names = FALSE)

                plot_loadings_enrichment_barplot(
                    enrich_df, label,
                    file.path(out_dir, paste0(label, "_enrichment.png"))
                )
            }
        }
    }

    if (length(all_results) == 0) return(NULL)
    combined <- .rbind_fill(all_results)
    rownames(combined) <- NULL
    write.csv(combined, file.path(out_dir, "mofa_weights_enrichment_all.csv"),
              row.names = FALSE)
    combined
}


#' Run compound ORA for metabolomics loadings/weights
#'
#' Maps top metabolomics features to KEGG compound IDs and runs compound ORA.
#' Used by loadings enrichment for DIABLO and MOFA2. Reuses
#' map_metabolite_ids_to_kegg() — the same mapper used by the per-omics
#' enrichment path — so loadings enrichment succeeds whenever metabolite-name
#' / synthetic-ID / bare-numeric feature IDs can be resolved to HMDB.
#'
#' @param feature_ids Character vector of top metabolomics feature IDs
#' @param harmonization_res Harmonization result (for HMDB mapping)
#' @param out_dir Output directory for CSVs and plots
#' @param label Label prefix for output files
#' @return data.frame of enriched pathways, or NULL
run_metabolite_loadings_ora <- function(feature_ids, harmonization_res, out_dir, label,
                                        exclude_classes = NULL) {
    # Build a pseudo-DE table: treat all top features as significant
    de_tbl <- data.frame(
        feature_id = feature_ids,
        pvalue = 0.01,
        padj = 0.01,
        stringsAsFactors = FALSE
    )

    # Use the same robust mapper as the per-omics enrichment path.
    id_map <- tryCatch(
        map_metabolite_ids_to_kegg(
            de_tables = list(loadings = de_tbl),
            harmonization_res = harmonization_res
        ),
        error = function(e) {
            message("    Metabolite ID mapping failed: ", e$message)
            NULL
        }
    )

    if (is.null(id_map) || nrow(id_map) == 0) {
        message("    Could not map metabolomics features to KEGG compound IDs")
        return(NULL)
    }

    de_mapped <- merge(de_tbl, id_map, by = "feature_id")
    de_mapped$KEGG_ID <- de_mapped$KEGG_CPD
    de_mapped <- de_mapped[!is.na(de_mapped$KEGG_ID) & !duplicated(de_mapped$KEGG_ID), ]

    if (nrow(de_mapped) < 2) {
        message("    Too few mapped metabolites for compound ORA (", nrow(de_mapped), ")")
        return(NULL)
    }

    # Use all measured KEGG-mapped metabolites as the background universe
    full_universe <- unique(id_map$KEGG_CPD[!is.na(id_map$KEGG_CPD)])

    enrich_df <- tryCatch(
        run_compound_ora(de_mapped, out_dir, 2, 500, 0.1, universe = full_universe,
                         exclude_classes = exclude_classes),
        error = function(e) {
            message("    Compound ORA failed: ", e$message)
            NULL
        }
    )

    if (!is.null(enrich_df) && nrow(enrich_df) > 0) {
        write.csv(enrich_df,
                  file.path(out_dir, paste0(label, "_enrichment.csv")),
                  row.names = FALSE)
        plot_loadings_enrichment_barplot(
            enrich_df, label,
            file.path(out_dir, paste0(label, "_enrichment.png"))
        )
        message("    Saved: ", file.path(out_dir, paste0(label, "_enrichment.png")))
    } else {
        message("    No enriched metabolomics pathways found for ", label)
    }

    enrich_df
}


#' Run ORA enrichment on a list of feature IDs
#'
#' Maps feature IDs to ENTREZ IDs and runs KEGG ORA via clusterProfiler.
#' Handles GENE_N synthetic IDs from the harmonized MAE by translating
#' them to WBGene IDs via the gene_protein_mapping table.
enrich_feature_list <- function(feature_ids, omics_type, harmonization_res,
                                 organism, kegg_org, org_db) {

    if (is.null(kegg_org) || is.null(org_db)) return(NULL)

    # Resolve IDs using the actual omics type
    resolved_ids <- resolve_gene_n_ids(feature_ids, harmonization_res, omics_type)

    # Map ALL features to ENTREZ (needed for both query and universe)
    # Then filter to only the query features for the enrichment test
    full_id_map <- tryCatch(
        map_feature_ids_to_entrez(
            de_tables = list(dummy = data.frame(feature_id = resolved_ids,
                                                 stringsAsFactors = FALSE)),
            omics_type = omics_type,
            harmonization_res = harmonization_res,
            org_db = org_db
        ),
        error = function(e) {
            message("    ID mapping failed: ", e$message)
            NULL
        }
    )

    if (is.null(full_id_map) || nrow(full_id_map) == 0) return(NULL)

    # Filter to only the query features (top-N loadings)
    query_map <- full_id_map[full_id_map$feature_id %in% resolved_ids, ]
    entrez_ids <- unique(query_map$ENTREZID[!is.na(query_map$ENTREZID)])
    if (length(entrez_ids) < 3) return(NULL)

    # Build universe from all measured features of the same omics type
    universe_entrez <- NULL
    om_key <- switch(omics_type,
                     "transcriptomics" = "transcriptomics",
                     "proteomics" = "proteomics",
                     omics_type)
    pre_data <- harmonization_res$inputs[[om_key]]
    if (!is.null(pre_data) && !is.null(pre_data$expr_work)) {
        all_features <- rownames(pre_data$expr_work)
        all_id_map <- tryCatch(
            map_feature_ids_to_entrez(
                de_tables = list(all = data.frame(feature_id = all_features,
                                                   stringsAsFactors = FALSE)),
                omics_type = omics_type,
                harmonization_res = harmonization_res,
                org_db = org_db
            ),
            error = function(e) NULL
        )
        if (!is.null(all_id_map)) {
            universe_entrez <- unique(all_id_map$ENTREZID[!is.na(all_id_map$ENTREZID)])
        }
    }

    message("    ORA (clusterProfiler): ", length(entrez_ids), " query ENTREZ IDs",
            if (!is.null(universe_entrez)) paste0(" / ", length(universe_entrez), " universe"))

    # Run KEGG ORA via clusterProfiler (handles KEGG ID conversion internally)
    enrich_res <- tryCatch({
        clusterProfiler::enrichKEGG(
            gene = entrez_ids,
            organism = kegg_org,
            keyType = "ncbi-geneid",
            universe = universe_entrez,
            pvalueCutoff = 0.1,
            minGSSize = 5,
            maxGSSize = 500
        )
    }, error = function(e) {
        message("    clusterProfiler::enrichKEGG failed: ", e$message)
        NULL
    })

    if (is.null(enrich_res)) return(NULL)

    df <- as.data.frame(enrich_res)
    if (nrow(df) == 0) return(NULL)

    # Standardize column names to match expected format
    result <- data.frame(
        pathway = df$Description,
        ID = df$ID,
        pvalue = df$pvalue,
        padj = df$p.adjust,
        GeneRatio = df$GeneRatio,
        setSize = df$Count,
        stringsAsFactors = FALSE
    )
    result <- result[order(result$pvalue), ]
    message("    Found ", nrow(result), " enriched KEGG pathways")
    result
}


#' Resolve GENE_N synthetic IDs to original feature IDs
#'
#' The harmonized MAE uses GENE_N IDs (where N = row in gene_protein_mapping).
#' This function translates them back to WBGene (for transcriptomics) or
#' protein IDs (for proteomics).
resolve_gene_n_ids <- function(feature_ids, harmonization_res, omics_type) {
    gpm <- harmonization_res$gene_protein_mapping
    if (is.null(gpm)) return(feature_ids)

    # Check if IDs look like GENE_N
    is_gene_n <- grepl("^GENE_\\d+$", feature_ids)
    if (!any(is_gene_n)) return(feature_ids)

    # Build lookup: GENE_N -> original ID
    if (omics_type == "transcriptomics") {
        id_col <- "gene_id"
    } else if (omics_type == "proteomics") {
        id_col <- "protein_id"
    } else {
        return(feature_ids)  # metabolomics uses feature_N, not GENE_N
    }

    lookup <- setNames(gpm[[id_col]], paste0("GENE_", seq_len(nrow(gpm))))

    resolved <- feature_ids
    resolved[is_gene_n] <- lookup[feature_ids[is_gene_n]]
    resolved <- resolved[!is.na(resolved)]

    n_mapped <- sum(is_gene_n) - sum(is.na(lookup[feature_ids[is_gene_n]]))
    message("    Resolved ", n_mapped, "/", sum(is_gene_n),
            " GENE_N IDs to ", id_col, " (", omics_type, ")")

    resolved
}


#' Plot barplot for loadings enrichment
plot_loadings_enrichment_barplot <- function(enrich_df, title, out_path, top_n = 15) {
    if (nrow(enrich_df) == 0) return(invisible(NULL))

    df <- enrich_df[order(enrich_df$pvalue), ]
    df <- df[seq_len(min(top_n, nrow(df))), ]

    df$label <- ifelse(nchar(df$pathway) > 45,
                        paste0(substr(df$pathway, 1, 42), "..."),
                        df$pathway)
    neg_log_p <- -log10(df$pvalue + 1e-300)
    neg_log_p <- pmin(neg_log_p, 15)

    png(out_path, width = 900, height = 600, res = 120)
    par(mar = c(5, 15, 3, 2))
    barplot(rev(neg_log_p), horiz = TRUE, names.arg = rev(df$label),
            las = 1, cex.names = 0.65, col = "steelblue",
            xlab = "-log10(p-value)",
            main = paste("Loadings Enrichment:", title))
    abline(v = -log10(0.05), col = "red", lty = 2)
    dev.off()
    message("    Saved: ", out_path)
}
