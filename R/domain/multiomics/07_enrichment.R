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

    # $pathway_results slot (nested by contrast)
    if (!is.null(enrich_res$pathway_results) && is.list(enrich_res$pathway_results)) {
        df <- collect_contrast_enrichment_dfs(enrich_res$pathway_results)
        if (!is.null(df)) return(df)
    }

    # Bare contrast-keyed list. The proteomics pathway module returns
    # list(<contrast> = list(<db>_<method> = data.frame)) with no
    # $pathway_results wrapper, so without this branch the whole proteomics
    # layer was silently dropped from the cross-omics combination.
    if (is.list(enrich_res)) {
        df <- collect_contrast_enrichment_dfs(enrich_res)
        if (!is.null(df)) return(df)
    }

    NULL
}


#' Collect enrichment tables out of a contrast-keyed list
#'
#' Walks a \code{list(<contrast> = data.frame)} or
#' \code{list(<contrast> = list(<sub_result> = data.frame))} and row-binds every
#' non-empty table it finds.
#'
#' @param contrast_list Named list keyed by contrast.
#' @return A single data frame, or NULL when no table was found.
collect_contrast_enrichment_dfs <- function(contrast_list) {
    dfs <- list()
    for (contrast_name in names(contrast_list)) {
        contrast_res <- contrast_list[[contrast_name]]
        if (is.data.frame(contrast_res) && nrow(contrast_res) > 0) {
            dfs[[contrast_name]] <- contrast_res
        } else if (is.list(contrast_res)) {
            # May have sub-results (e.g., KEGG, GO)
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
    if (length(dfs) == 0) return(NULL)
    dplyr::bind_rows(dfs)
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
                                 universe = full_universe)
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
                              universe = NULL) {

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

        # `pathway` must hold the KEGG map id, not the human-readable name:
        # merge_pathway_pvalues() joins the omics layers on this column and the
        # gene layers emit map##### ids, so keying metabolites on the name made
        # the metabolite layer unable to ever intersect them. The name is kept
        # alongside in `pathway_name`, matching the gene-layer tables.
        results[[pw]] <- data.frame(
            pathway = pw,
            pathway_name = if (!is.null(pathway_names[pw]) && !is.na(pathway_names[pw]))
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
get_kegg_organism <- function(organism) {
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

    if (!is.null(ora_res)) return(ora_res)

    # Fallback: Fisher's exact test with KEGG REST pathway-gene links
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

    # Find common pathways across omics
    all_pathways <- lapply(pathway_tables, function(df) {
        if ("pathway" %in% names(df)) df$pathway
        else if ("ID" %in% names(df)) df$ID
        else if ("Description" %in% names(df)) df$Description
        else rownames(df)
    })

    # Use union of all pathways (not just intersection) for broader view
    union_pathways <- unique(unlist(all_pathways))
    common_pathways <- Reduce(intersect, all_pathways)

    if (length(union_pathways) == 0) {
        message("No pathways found across omics layers")
        return(NULL)
    }

    message(sprintf("  Found %d total pathways (%d in common) across %d omics layers",
                    length(union_pathways), length(common_pathways), length(omics)))

    # Use union for the heatmap (show all enriched), common for meta-analysis
    use_pathways <- if (length(common_pathways) >= 5) common_pathways else union_pathways

    # Merge pathway p-values for meta-analysis
    merged_pathways <- merge_pathway_pvalues(pathway_tables, use_pathways, omics)

    # Combine p-values using Stouffer's method (only for pathways with >= 2 p-values)
    meta_results <- stouffer_combined_pvalues(merged_pathways)

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

        # 3. Per-omics enrichment bar plots
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

#' Normalize KEGG pathway IDs to bare numeric form
#'
#' Strips organism prefixes (hsa, mmu, cel, map, etc.) to allow
#' joining gene-based and compound-based pathway results.
#' @param ids Character vector of KEGG pathway IDs
#' @return Character vector of numeric-only pathway IDs (e.g., "00010")
normalize_kegg_pathway_id <- function(ids) {
    sub("^[a-zA-Z]+", "", ids)
}


#' Merge pathway p-values from multiple omics
merge_pathway_pvalues <- function(pathway_tables, target_pathways, omics) {

    merged <- data.frame(pathway = target_pathways, stringsAsFactors = FALSE)

    for (om in omics) {
        df <- pathway_tables[[om]]

        # Identify pathway column
        pathway_col <- if ("pathway" %in% names(df)) "pathway"
                       else if ("ID" %in% names(df)) "ID"
                       else if ("Description" %in% names(df)) "Description"
                       else NULL

        if (is.null(pathway_col)) {
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

        # For multiple contrasts, take the minimum p-value per pathway
        df_agg <- aggregate(
            stats::as.formula(paste(pval_col, "~", pathway_col)),
            data = df,
            FUN = min
        )
        colnames(df_agg) <- c("pathway", paste0("pval_", om))

        # Subset to target pathways
        df_sub <- df_agg[df_agg$pathway %in% target_pathways, , drop = FALSE]

        # Merge
        merged <- merge(merged, df_sub, by = "pathway", all.x = TRUE)
    }

    merged
}


#' Combine p-values using Fisher's method
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

#' Plot cross-omics pathway heatmap
plot_cross_omics_pathway_heatmap <- function(meta_results, omics, top_n = 30) {

    # Select top N pathways by combined p-value
    top_pathways <- meta_results[seq_len(min(top_n, nrow(meta_results))), ]

    pval_cols <- grep("^pval_", names(top_pathways), value = TRUE)
    pval_matrix <- as.matrix(top_pathways[, pval_cols, drop = FALSE])

    # Truncate long pathway names
    pathway_labels <- top_pathways$pathway
    pathway_labels <- ifelse(nchar(pathway_labels) > 50,
                             paste0(substr(pathway_labels, 1, 47), "..."),
                             pathway_labels)
    rownames(pval_matrix) <- pathway_labels

    # Transform to -log10(p)
    log_pval_matrix <- -log10(pval_matrix + 1e-300)
    # Cap at 10 for display
    log_pval_matrix[log_pval_matrix > 10] <- 10
    colnames(log_pval_matrix) <- gsub("^pval_", "", colnames(log_pval_matrix))

    # Replace NA with 0
    log_pval_matrix[is.na(log_pval_matrix)] <- 0

    # Heatmap
    if (requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(log_pval_matrix,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           main = "Cross-Omics Pathway Enrichment (-log10 p-value)",
                           color = colorRampPalette(c("white", "gold", "orange", "red"))(50),
                           fontsize_row = 7, fontsize_col = 10,
                           angle_col = 45,
                           na_col = "grey90",
                           border_color = "grey80")
    } else {
        heatmap(log_pval_matrix, scale = "none", Colv = NA,
                main = "Cross-Omics Pathway Enrichment",
                col = colorRampPalette(c("white", "orange", "red"))(50))
    }
}


#' Plot enrichment dot plot
plot_enrichment_dotplot <- function(meta_results, omics, top_n = 20) {

    top <- meta_results[seq_len(min(top_n, nrow(meta_results))), ]

    pval_cols <- grep("^pval_", names(top), value = TRUE)

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

    # Truncate long names
    plot_df$pathway <- ifelse(nchar(plot_df$pathway) > 45,
                              paste0(substr(plot_df$pathway, 1, 42), "..."),
                              plot_df$pathway)

    # Reverse pathway order for bottom-to-top display
    pathway_order <- rev(unique(top$pathway))
    pathway_order <- ifelse(nchar(pathway_order) > 45,
                            paste0(substr(pathway_order, 1, 42), "..."),
                            pathway_order)
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
                                            out_dir, top_n = 50) {

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
                    top_feat_ids, harmonization_res, out_dir, label
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
                                         out_dir, top_n = 50) {

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
                    top_feat_ids, harmonization_res, out_dir, label
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
run_metabolite_loadings_ora <- function(feature_ids, harmonization_res, out_dir, label) {
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
        run_compound_ora(de_mapped, out_dir, 2, 500, 0.1, universe = full_universe),
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
