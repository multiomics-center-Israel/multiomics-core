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
    excl_classes <- config$modes$multiomics$enrichment$exclude_pathway_classes

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
                ora_res <- run_compound_ora(de_mapped, out_dir, use_min_gs, max_gs,
                                            pval_cutoff,
                                            exclude_classes = excl_classes,
                                            universe = full_universe,
                                            label = contrast_name)
                # filter_ora_rows() switches to a literal `method` match the
                # moment any row in the table carries one, so the ORA half has
                # to name itself now that GSEA rows sit beside it.
                if (is.data.frame(ora_res) && nrow(ora_res) > 0 &&
                    !"method" %in% colnames(ora_res)) {
                    ora_res$method <- "ora"
                }
                # Metabolites get both: ORA for the crowding question, GSEA for
                # an NES the cross-omics heatmap can colour a cell with.
                gsea_res <- run_compound_gsea(de_mapped, out_dir, use_min_gs, max_gs,
                                              seed = config$params$seed %||% 1L,
                                              exclude_classes = excl_classes,
                                              label = contrast_name)
                .rbind_fill(list(ora = ora_res, gsea = gsea_res))
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

    # .rbind_fill(), not rbind(): a contrast can contribute ORA rows only (no
    # compound pathway large enough to score) while another contributes both,
    # and rbind() aborts on the differing column sets.
    combined <- .rbind_fill(all_results)
    if (is.null(combined) || nrow(combined) == 0) return(NULL)
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
                    # Kept for compound GSEA: the moderated t is the ranking
                    # statistic, and it is unrecoverable from log2fc + p alone
                    # once the table has been standardised.
                    statistic = if ("statistic" %in% names(df)) df$statistic
                                else if ("t" %in% names(df)) df$t
                                else NA_real_,
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
#' Two routes are combined. The HMDB route translates the row_data HMDB column
#' through the HMDB-to-KEGG mapping file; the direct route reads a `KEGG` column
#' when the annotation table already carries one. Many vendor tables annotate
#' compounds KEGG-side that have no HMDB accession, so consulting only HMDB
#' silently drops them from every downstream compound test.
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
    # Find a pre-annotated KEGG compound column
    kegg_col <- intersect(c("KEGG", "kegg", "KEGG_ID", "KEGG_CPD", "KEGG_COMPOUND"),
                          colnames(row_data))

    if (length(hmdb_col) == 0 && length(kegg_col) == 0) {
        message("    No HMDB or KEGG column in metabolomics row_data")
        return(NULL)
    }

    # Load HMDB -> KEGG compound mapping. Missing is only fatal when it is the
    # only route available.
    hmdb_kegg <- if (length(hmdb_col) > 0) load_hmdb_to_kegg_map() else NULL
    if (is.null(hmdb_kegg) && length(kegg_col) == 0) return(NULL)

    # Collect all unique feature_ids from DE tables
    all_de_ids <- unique(unlist(lapply(de_tables, function(df) df$feature_id)))

    # Check if DE feature_ids are already HMDB IDs (direct mapping, no row_data needed)
    hmdb_pattern <- sum(grepl("^HMDB\\d+$", all_de_ids, ignore.case = TRUE), na.rm = TRUE)
    uses_hmdb_direct <- !is.null(hmdb_kegg) &&
        hmdb_pattern > length(all_de_ids) * 0.5

    hmdb_ids <- if (length(hmdb_col) > 0) {
        as.character(row_data[[hmdb_col[1]]])
    } else {
        rep(NA_character_, nrow(row_data))
    }

    if (uses_hmdb_direct) {
        # DE tables already use HMDB IDs as feature_ids — direct mapping
        message("    DE tables use HMDB IDs as feature_ids (direct mapping)")
        feat_ids <- hmdb_ids
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
    }

    # Route 1: HMDB accession -> KEGG compound via the mapping file.
    kegg_from_hmdb <- if (!is.null(hmdb_kegg)) {
        as.character(hmdb_kegg[hmdb_ids])
    } else {
        rep(NA_character_, length(hmdb_ids))
    }

    # Route 2: a KEGG compound column the annotation table already carries.
    # Take the first id when a row lists several.
    kegg_direct <- if (length(kegg_col) > 0) {
        raw <- trimws(as.character(row_data[[kegg_col[1]]]))
        raw[!nzchar(raw) | raw %in% c("NA", "-")] <- NA_character_
        sub("[;,|/[:space:]].*$", "", raw)
    } else {
        rep(NA_character_, nrow(row_data))
    }

    # HMDB wins where both routes fire, so nothing that maps today changes id;
    # the direct column only fills the gaps.
    kegg_cpds <- kegg_from_hmdb
    fill <- is.na(kegg_cpds) | !nzchar(kegg_cpds)
    kegg_cpds[fill] <- kegg_direct[fill]

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
    n_hmdb_route <- sum(!is.na(kegg_from_hmdb) & nzchar(kegg_from_hmdb))
    n_direct_only <- sum(fill & !is.na(kegg_direct) & nzchar(kegg_direct))
    message("    Mapped ", nrow(mapped), "/", n_total,
            " metabolites to KEGG compound IDs (", n_hmdb_route, " via HMDB, ",
            n_direct_only, " added from the KEGG column)")

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
#' Hypergeometric ORA computed directly with `stats::phyper()` against
#' pathway-compound sets built from the KEGG compound-pathway associations —
#' `clusterProfiler` is not involved. Compound sets are restricted to the
#' measured background before sizing, so the universe is the measured
#' metabolome, not all of KEGG.
#'
#' The returned table is filtered to `pvalue < pval_cutoff`, which is a
#' candidate list rather than a set of FDR-significant pathways. The full
#' tested table (every pathway that passed the size filter, unfiltered by
#' p-value, with its BH-adjusted p) is written to `out_dir` and attached to the
#' return value as the `all_tested` attribute so the correction stays auditable.
#'
#' @param de_mapped DE table with KEGG_ID (compound IDs) column
#' @param cache_dir Directory for caching KEGG data
#' @param min_gs Minimum compound set size, honoured as given (floored at 2,
#'   below which the hypergeometric test is degenerate)
#' @param max_gs Maximum compound set size
#' @param pval_cutoff Nominal p-value cutoff applied to the returned table
#' @param universe Character vector of all measured KEGG compound IDs. When
#'   NULL, the compounds present in `de_mapped` are used.
#' @param out_dir Directory to write the full tested table to; NULL to skip the
#'   write. Defaults to `cache_dir`, which is what every caller passes.
#' @param label Optional prefix for the full-table filename, so several
#'   contrasts writing to one directory do not overwrite each other.
#' @return data.frame of candidate pathways (`pvalue < pval_cutoff`) with an
#'   `all_tested` attribute holding every tested pathway, or NULL.
run_compound_ora <- function(de_mapped, cache_dir, min_gs, max_gs, pval_cutoff,
                              universe = NULL, out_dir = cache_dir,
                              label = NULL, exclude_classes = NULL) {

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
    query_rule <- "padj<0.05"
    sig_cpds <- de_mapped$KEGG_ID[!is.na(de_mapped$padj) & de_mapped$padj < 0.05]
    sig_cpds <- unique(sig_cpds[!is.na(sig_cpds)])

    if (length(sig_cpds) < 3) {
        # Relax to nominal p-value. Recorded on every row: the relaxed query is
        # a different, weaker question than the FDR one and the reader has no
        # other way to tell which was asked.
        query_rule <- "pvalue<0.05 (relaxed: <3 compounds at padj<0.05)"
        sig_cpds <- de_mapped$KEGG_ID[!is.na(de_mapped$pvalue) & de_mapped$pvalue < 0.05]
        sig_cpds <- unique(sig_cpds[!is.na(sig_cpds)])
    }

    if (length(sig_cpds) < 2) {
        message("    Too few significant compounds for ORA (", length(sig_cpds), ")")
        return(NULL)
    }

    message("    Running compound ORA: ", length(sig_cpds), " significant / ",
            length(all_cpds), " total compounds [query rule: ", query_rule, "]")

    # Build pathway -> compound sets
    pathway_sets <- split(cpd_pathways$compound, cpd_pathways$pathway)

    # Build pathway name lookup
    pathway_names <- stats::setNames(cpd_pathways$name, cpd_pathways$pathway)
    pathway_names <- pathway_names[!duplicated(names(pathway_names))]

    # Filter pathways by size (intersection with measured compounds).
    # The configured minimum is honoured: silently clamping it to 3 let
    # 3-compound pathways top the table under a config asking for >= 10.
    use_min_gs <- max(2, min_gs)
    N <- length(all_cpds)  # total measured compounds
    k <- length(sig_cpds)  # significant compounds

    # Drop KEGG classes this organism cannot have before testing them, so they
    # never reach the reported table or the multiple-testing family.
    if (!is.null(exclude_classes) && length(exclude_classes) > 0) {
        keep_pw <- keep_kegg_pathways(names(pathway_sets),
                                      exclude = unlist(exclude_classes),
                                      cache_dir = cache_dir,
                                      label = "compound pathways")
        pathway_sets <- pathway_sets[keep_pw]
    }

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
            query_rule = query_rule,
            stringsAsFactors = FALSE
        )
    }

    if (length(results) == 0) {
        message("    No compound pathways passed the size filter (min ",
                use_min_gs, ", max ", max_gs, ")")
        return(NULL)
    }

    df <- do.call(rbind, results)
    rownames(df) <- NULL

    # Multiple testing correction over every pathway actually tested
    df$padj <- stats::p.adjust(df$pvalue, method = "BH")
    df <- df[order(df$pvalue), ]

    all_tested <- df
    message("    Tested ", nrow(all_tested), " compound pathways (min set size ",
            use_min_gs, "); ", sum(all_tested$padj < 0.05, na.rm = TRUE),
            " at BH padj < 0.05")

    if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        fname <- if (is.null(label)) {
            "compound_ora_all_tested.csv"
        } else {
            paste0(label, "_compound_ora_all_tested.csv")
        }
        utils::write.csv(all_tested, file.path(out_dir, fname), row.names = FALSE)
    }

    # The returned table is the nominal-p candidate list the report consumes.
    # It is deliberately NOT an FDR-significant set — `padj` travels with every
    # row and `all_tested` holds what it was corrected against.
    df <- df[df$pvalue < pval_cutoff, ]
    if (nrow(df) == 0) {
        message("    No compound pathways below nominal p < ", pval_cutoff)
        return(NULL)
    }

    message("    Found ", nrow(df), " candidate compound pathways at nominal p < ",
            pval_cutoff)
    attr(df, "all_tested") <- all_tested
    df
}


#' Rank mapped metabolites for compound GSEA
#'
#' Ranks on the moderated t when the DE table carries usable values, falling
#' back to sign(log2FC) * -log10(p). The test is for finite VALUES rather than
#' for the column: the standardised metabolomics DE table always carries a
#' `statistic` column and fills it with NA when the source export has none, and
#' an all-NA rank vector makes fgsea report "no gene set overlap" — a message
#' that sends the reader looking for an id-mapping bug that is not there.
#'
#' @param de_mapped DE table already merged with KEGG compound ids: columns
#'   `KEGG_ID`, `pvalue`, `log2fc`, and optionally `statistic`.
#' @return Named numeric vector of ranks (names are KEGG compound ids) sorted
#'   decreasing, one entry per compound; `numeric(0)` when nothing is rankable.
rank_compounds_for_gsea <- function(de_mapped) {
    if (!is.data.frame(de_mapped) || nrow(de_mapped) == 0) return(numeric(0))
    if (!"KEGG_ID" %in% colnames(de_mapped)) return(numeric(0))

    stat <- if ("statistic" %in% colnames(de_mapped)) {
        suppressWarnings(as.numeric(de_mapped$statistic))
    } else {
        NULL
    }

    if (!is.null(stat) && any(is.finite(stat))) {
        ranks <- stat
    } else {
        lfc <- suppressWarnings(as.numeric(de_mapped$log2fc))
        pval <- suppressWarnings(as.numeric(de_mapped$pvalue))
        ranks <- sign(lfc) * -log10(pval + 1e-300)
    }

    names(ranks) <- as.character(de_mapped$KEGG_ID)
    keep <- is.finite(ranks) & !is.na(names(ranks)) & nzchar(names(ranks))
    ranks <- ranks[keep]
    if (length(ranks) == 0) return(numeric(0))

    # Several metabolites can annotate to the same compound; fgsea would score
    # that compound twice. The strongest |rank| wins, ties broken on the compound
    # id so the vector does not depend on the DE table's row order.
    ranks <- ranks[order(-abs(ranks), names(ranks))]
    ranks <- ranks[!duplicated(names(ranks))]
    sort(ranks, decreasing = TRUE)
}


#' Run GSEA for metabolite compounds against KEGG compound pathways
#'
#' The rank-based companion to `run_compound_ora()`, not a replacement for it.
#' ORA asks whether the significant compounds crowd into a pathway; this asks
#' whether a pathway's compounds sit at one end of the whole measured ranking.
#' It needs no significance cutoff and it yields an NES — the signed value the
#' cross-omics heatmap draws for the gene layers and, with ORA alone, could
#' never draw for metabolites.
#'
#' Compound sets are intersected with the measured, KEGG-mapped metabolites
#' before the size filter, so a reported `size` is the number of compounds
#' actually scored rather than KEGG's full membership. With a few hundred
#' mapped metabolites most compound pathways fall below any sensible minimum;
#' the count that survives is reported rather than worked around.
#'
#' @param de_mapped DE table merged with KEGG compound ids (see
#'   `rank_compounds_for_gsea()` for the columns used).
#' @param cache_dir Directory holding the cached KEGG compound-pathway and
#'   pathway-class tables.
#' @param min_gs Minimum measured compounds per pathway, honoured as given and
#'   floored at 2, below which a set is not a set.
#' @param max_gs Maximum measured compounds per pathway.
#' @param seed Integer seed for fgsea's stochastic multilevel step.
#' @param exclude_classes KEGG pathway classes to drop before testing, from
#'   `enrichment.exclude_pathway_classes`; NULL tests every class.
#' @param out_dir Directory to write the full scored table to; NULL to skip the
#'   write. Defaults to `cache_dir`, matching `run_compound_ora()`.
#' @param label Optional filename prefix so several contrasts writing into one
#'   directory do not overwrite each other.
#' @return data.frame of every scored pathway in the fgsea column shape the gene
#'   layers emit (pathway, pathway_name, pval, padj, NES, size, database,
#'   method), or NULL when no pathway was testable.
run_compound_gsea <- function(de_mapped, cache_dir, min_gs, max_gs,
                               seed = 1L, exclude_classes = NULL,
                               out_dir = cache_dir, label = NULL) {

    cpd_pathways <- get_kegg_compound_pathways(cache_dir)
    if (is.null(cpd_pathways) || nrow(cpd_pathways) == 0) {
        message("    Could not retrieve KEGG compound-pathway associations")
        return(NULL)
    }

    ranks <- rank_compounds_for_gsea(de_mapped)
    if (length(ranks) < 3) {
        message("    Too few rankable compounds for GSEA (", length(ranks), ")")
        return(NULL)
    }

    pathway_sets <- split(cpd_pathways$compound, cpd_pathways$pathway)
    pathway_names <- stats::setNames(cpd_pathways$name, cpd_pathways$pathway)
    pathway_names <- pathway_names[!duplicated(names(pathway_names))]

    # The same class exclusion the other layers get, applied before testing so
    # the dropped maps never enter the multiple-testing family.
    if (!is.null(exclude_classes) && length(exclude_classes) > 0) {
        keep_pw <- keep_kegg_pathways(names(pathway_sets),
                                      exclude = unlist(exclude_classes),
                                      cache_dir = cache_dir,
                                      label = "compound pathways")
        pathway_sets <- pathway_sets[keep_pw]
    }

    pathway_sets <- lapply(pathway_sets,
                           function(cpds) unique(intersect(cpds, names(ranks))))
    use_min_gs <- max(2, min_gs)
    set_sizes <- lengths(pathway_sets)
    testable <- set_sizes >= use_min_gs & set_sizes <= max_gs

    message("    Compound GSEA: ", length(ranks), " ranked compounds; ",
            sum(testable), "/", length(set_sizes), " pathways carry ",
            use_min_gs, "-", max_gs, " measured compounds")

    if (!any(testable)) {
        message("    No compound pathway is large enough to score")
        return(NULL)
    }
    pathway_sets <- pathway_sets[testable]

    # fgseaMultilevel is stochastic: unseeded, borderline pathways cross
    # padj = 0.05 between otherwise identical runs.
    fgsea_res <- withr::with_seed(seed, fgsea::fgsea(
        pathways = pathway_sets,
        stats = ranks,
        minSize = use_min_gs,
        maxSize = max_gs,
        nPermSimple = 10000
    ))

    df <- as.data.frame(fgsea_res)
    if (nrow(df) == 0) {
        message("    Compound GSEA scored no pathways")
        return(NULL)
    }

    pw_ids <- as.character(df$pathway)
    pw_names <- unname(pathway_names[pw_ids])
    pw_names[is.na(pw_names) | !nzchar(pw_names)] <- pw_ids[is.na(pw_names) | !nzchar(pw_names)]

    # `pathway` holds the KEGG map id, matching run_compound_ora() and the gene
    # layers: merge_pathway_pvalues() joins the layers on this column, and a
    # human-readable name here would leave metabolites unable to intersect them.
    out <- data.frame(
        pathway = pw_ids,
        pathway_name = pw_names,
        ID = pw_ids,
        pval = df$pval,
        padj = df$padj,
        ES = df$ES,
        NES = df$NES,
        size = df$size,
        leadingEdge = vapply(df$leadingEdge, paste, character(1), collapse = ","),
        database = "KEGG",
        method = "fgsea",
        stringsAsFactors = FALSE
    )
    out <- out[order(out$pval), ]
    rownames(out) <- NULL

    finite_nes <- out$NES[is.finite(out$NES)]
    message("    Compound GSEA: ", sum(out$padj < 0.05, na.rm = TRUE),
            " of ", nrow(out), " pathways at BH padj < 0.05",
            if (length(finite_nes) > 0) {
                paste0(" (NES ", round(min(finite_nes), 2), " to ",
                       round(max(finite_nes), 2), ")")
            } else "")

    # Every scored pathway is returned, matching what the gene layers hand the
    # cross-omics step; the same table goes to disk so the NES a figure shows
    # can be traced back without re-running fgsea.
    if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        fname <- if (is.null(label)) {
            "compound_gsea_all_tested.csv"
        } else {
            paste0(label, "_compound_gsea_all_tested.csv")
        }
        utils::write.csv(out, file.path(out_dir, fname), row.names = FALSE)
    }

    out
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
#' @param organism Organism identifier from `config$global$organism`, e.g.
#'   "human" or "Caenorhabditis elegans".
#' @return Three-letter KEGG organism code, or NULL when the organism is
#'   unknown, absent, or not a single usable string. Callers treat NULL as
#'   "no KEGG for this organism" and skip the KEGG path.
get_kegg_organism <- function(organism) {
    # A missing/blank organism must not reach kegg_map[[...]]: a zero-length or
    # NA subscript there is an error, not a lookup miss, and non-model configs
    # legitimately leave this unset.
    if (is.null(organism) || length(organism) != 1L) return(NULL)
    if (is.na(organism) || !nzchar(trimws(as.character(organism)))) return(NULL)

    organism <- as.character(organism)

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

    # Always work from the union. The intersection was used whenever it reached 5
    # pathways, which made the whole analysis depend on the narrowest layer: when
    # metabolomics carried 6 compound-ORA rows the union was taken (4,471
    # pathways), and adding compound GSEA to that layer silently collapsed the
    # table to the 41 pathways all three happen to share. Row *selection* for the
    # figures already prefers pathways with several layers, so the ranking shows
    # where the omics meet without the table having to hide everything else.
    use_pathways <- union_pathways

    # KEGG reference maps are pan-species. For an organism that has no KEGG code
    # the vertebrate organ and human-disease maps light up purely because the
    # underlying orthologs are generic, so they are dropped rather than reported
    # as findings. Configurable; nothing is excluded unless the config asks.
    excl <- config$modes$multiomics$enrichment$exclude_pathway_classes
    if (!is.null(excl) && length(excl) > 0) {
        keep <- keep_kegg_pathways(use_pathways, exclude = unlist(excl),
                                   cache_dir = out_dir,
                                   label = "cross-omics pathways")
        use_pathways <- use_pathways[keep]

        # The per-omics tables are written to disk and drive the per-layer
        # barplots, so they need the same treatment; filtering only the merged
        # selection would leave the excluded classes visible one section away.
        pathway_tables <- lapply(pathway_tables, function(df) {
            col <- if ("pathway" %in% names(df)) "pathway"
                   else if ("ID" %in% names(df)) "ID" else NULL
            if (is.null(col)) return(df)
            df[keep_kegg_pathways(df[[col]], exclude = unlist(excl),
                                  cache_dir = out_dir,
                                  label = "per-omics pathways"), , drop = FALSE]
        })
    }

    # One significance statistic for the whole section — figures and meta table
    # alike — so a reader never has to guess whether a bar is raw or adjusted.
    pvalue_type <- resolve_pvalue_type(config)

    # Merge pathway significance values for meta-analysis
    merged_pathways <- merge_pathway_pvalues(pathway_tables, use_pathways, omics,
                                             pvalue_type = pvalue_type)
    types_used <- attr(merged_pathways, "pvalue_type_used")
    value_label <- describe_pvalue_axis(pvalue_type, types_used)

    # Combine across layers using Stouffer's method (only for pathways with >= 2
    # values). With pvalue_type = "padj" the inputs are already BH-adjusted, so
    # combined_padj is conservative — treat it as a ranking key, not a calibrated
    # test.
    meta_results <- stouffer_combined_pvalues(merged_pathways)

    # Sort by combined p-value. This ordering is the contract of the meta CSV —
    # the plots re-order their own rows (see select_multi_omics_pathways()).
    meta_results <- meta_results[order(meta_results$combined_pval), ]

    # Carry the readable term names over from the per-omics tables, otherwise
    # every label downstream is a bare accession (GO:0034329, map01523).
    meta_results <- attach_pathway_names_from_tables(meta_results, pathway_tables)

    # The pval_<omics> column names are a fixed contract, so the table has to say
    # in-band which statistic they actually hold.
    fell_back <- names(types_used)[types_used != pvalue_type]
    meta_results$pvalue_type <- if (length(fell_back) == 0) pvalue_type else
        paste0(pvalue_type, " (", if (identical(pvalue_type, "padj")) "raw p" else "padj",
               " for ", paste(fell_back, collapse = ", "), ")")

    # Generate plots
    plots <- list()
    if (!is.null(out_dir) && nrow(meta_results) > 0) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        # 1. One heatmap per gene-set collection, not one combined figure.
        #    The collections do not cover the same layers or the same ID space:
        #    GO is annotated for the gene layers only, while the KEGG map#####
        #    space is the one all three share. Mixing them into a single figure
        #    also let the largest collection win every row -- GO contributes
        #    4,107 pathways here against KEGG's 195, so a flat top-30 was
        #    entirely GO and the first KEGG map sat at rank 39.
        #
        #    Each figure keeps only the layers that actually scored something in
        #    that collection, so an empty column means "this layer has no
        #    annotation here", not "this layer found nothing".
        coll <- if (exists("classify_pathway_collection")) {
            classify_pathway_collection(as.character(meta_results$pathway))
        } else {
            rep("all", nrow(meta_results))
        }

        for (cl in unique(coll[!is.na(coll)])) {
            sub <- meta_results[coll == cl, , drop = FALSE]
            if (nrow(sub) == 0) next

            # Layers with at least one value in this collection.
            keep_omics <- omics[vapply(omics, function(om) {
                cn <- paste0("pval_", om)
                cn %in% names(sub) && any(!is.na(sub[[cn]]))
            }, logical(1))]
            if (length(keep_omics) == 0) next

            slug <- gsub("[^A-Za-z0-9]+", "_", cl)
            key  <- paste0("pathway_heatmap_", slug)
            plots[[key]] <- file.path(
                out_dir, paste0("cross_omics_pathway_heatmap_", slug, ".png"))
            png(plots[[key]], width = 1200, height = 900, res = 120)
            tryCatch({
                plot_cross_omics_pathway_heatmap(sub, keep_omics,
                                                 pathway_tables = pathway_tables,
                                                 value_label = value_label,
                                                 pvalue_type = pvalue_type,
                                                 collection = cl)
            }, error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("Heatmap failed:", e$message), cex = 1.2)
            })
            dev.off()

            # The same collection seen through over-representation only. The NES
            # figure can only colour a layer that ran GSEA, so a layer tested end
            # to end by over-representation has no column there; here it does.
            key_ora <- paste0("pathway_heatmap_ora_", slug)
            plots[[key_ora]] <- file.path(
                out_dir, paste0("cross_omics_pathway_heatmap_ora_", slug, ".png"))
            png(plots[[key_ora]], width = 1200, height = 900, res = 120)
            tryCatch({
                plot_cross_omics_ora_heatmap(sub, keep_omics,
                                             pathway_tables = pathway_tables,
                                             collection = cl)
            }, error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("ORA heatmap failed:", e$message), cex = 1.2)
            })
            dev.off()
        }

        # 2. Dot plot of top pathways per omics
        plots$dot_plot <- file.path(out_dir, "cross_omics_enrichment_dotplot.png")
        png(plots$dot_plot, width = 1200, height = 800, res = 120)
        tryCatch({
            plot_enrichment_dotplot(meta_results, omics, value_label = value_label)
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
                plot_per_omics_barplot(pt, om, pvalue_type = pvalue_type)
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
        pvalue_type = pvalue_type,
        pvalue_type_used = types_used,
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


#' Attach readable pathway names to a cross-omics meta-analysis table
#'
#' The per-omics enrichment tables already carry a `pathway_name` column;
#' `merge_pathway_pvalues()` keeps only the accession, so the names have to be
#' joined back on. This deliberately does not reuse `add_pathway_names()` from
#' `R/core/09_enrichment.R`: that helper re-derives names from a gene-set
#' collection or GO.db, which is neither available nor needed here — the names
#' the omics layers already agreed on are the ones we want.
#'
#' @param meta_results Meta-analysis data frame with a `pathway` column.
#' @param pathway_tables Named list of per-omics enrichment data frames, each
#'   optionally carrying `pathway` and `pathway_name` columns.
#' @return `meta_results` with a `pathway_name` column added (NA where no source
#'   table names that accession). Returned unchanged if no table carries names.
attach_pathway_names_from_tables <- function(meta_results, pathway_tables) {
    if (is.null(meta_results) || nrow(meta_results) == 0) return(meta_results)

    name_lookup <- do.call(rbind, lapply(pathway_tables, function(df) {
        if (is.data.frame(df) && all(c("pathway", "pathway_name") %in% names(df))) {
            unique(data.frame(pathway = as.character(df$pathway),
                              pathway_name = as.character(df$pathway_name),
                              stringsAsFactors = FALSE))
        } else {
            NULL
        }
    }))

    if (is.null(name_lookup) || nrow(name_lookup) == 0) return(meta_results)

    # An accession can be named by more than one layer; first name wins so the
    # label is stable across reruns.
    name_lookup <- name_lookup[!duplicated(name_lookup$pathway), , drop = FALSE]
    meta_results$pathway_name <- name_lookup$pathway_name[
        match(as.character(meta_results$pathway), name_lookup$pathway)]

    meta_results
}


#' Select the cross-omics pathways to show in a plot
#'
#' Ranking by combined p-value alone lets whichever layer contributed the
#' largest collection monopolise the plot — a genome-wide GO GMT on one omics
#' produces thousands of tiny p-values and buries every pathway the layers
#' actually share, which is the opposite of what a cross-omics figure is for.
#' Ranking by the number of contributing layers first surfaces the overlap.
#'
#' @param meta_results Meta-analysis data frame from `stouffer_combined_pvalues()`.
#' @param top_n Maximum number of rows to return.
#' @return Data frame with at most `top_n` rows, ordered by `n_omics` descending
#'   then `combined_pval` ascending.
select_multi_omics_pathways <- function(meta_results, top_n = 30) {
    if (is.null(meta_results) || nrow(meta_results) == 0) return(meta_results)

    # Rank on how many layers are ENRICHED, not how many were merely tested.
    # n_omics counts non-missing p-values, so ordering by it alone floated
    # pathways that every layer measured and none found anything in above
    # pathways genuinely significant in two layers.
    pval_cols <- grep("^pval_", names(meta_results), value = TRUE)
    n_sig <- if (length(pval_cols) > 0) {
        pm <- suppressWarnings(vapply(meta_results[pval_cols], as.numeric,
                                      numeric(nrow(meta_results))))
        pm <- matrix(pm, nrow = nrow(meta_results))
        rowSums(!is.na(pm) & pm < 0.05)
    } else {
        rep(0L, nrow(meta_results))
    }
    n_omics <- if ("n_omics" %in% names(meta_results)) {
        as.numeric(meta_results$n_omics)
    } else {
        rep(1, nrow(meta_results))
    }
    combined <- if ("combined_pval" %in% names(meta_results)) {
        as.numeric(meta_results$combined_pval)
    } else {
        rep(NA_real_, nrow(meta_results))
    }

    ord <- order(-n_sig, -n_omics, combined, na.last = TRUE)
    keep <- min(top_n, nrow(meta_results))

    # Draw in turn from each gene-set collection rather than taking a flat top-N.
    # The collections differ by an order of magnitude in size -- 4,107 GO terms
    # against 195 KEGG maps in this run -- so a flat cut is won outright by the
    # largest one: every row was GO and the first KEGG map sat at rank 39, which
    # made the figure look like a GO-only analysis. Ranking within a collection
    # is unchanged; only the interleaving is new.
    coll <- if (exists("classify_pathway_collection")) {
        classify_pathway_collection(as.character(meta_results$pathway))[ord]
    } else {
        rep("all", length(ord))
    }
    by_coll <- split(ord, coll)
    # Best-ranked collection first, so the strongest evidence still leads.
    by_coll <- by_coll[order(vapply(by_coll, function(i) match(i[1], ord), numeric(1)))]

    picked <- integer(0)
    rank_i <- 1L
    while (length(picked) < keep && any(lengths(by_coll) >= rank_i)) {
        for (idx in by_coll) {
            if (length(picked) >= keep) break
            if (length(idx) >= rank_i) picked <- c(picked, idx[rank_i])
        }
        rank_i <- rank_i + 1L
    }
    picked <- picked[order(match(picked, ord))]
    meta_results[picked, , drop = FALSE]
}


#' Build readable pathway labels from accessions and names
#'
#' @param ids Character vector of pathway accessions (GO:xxxxxxx, map#####, ...).
#' @param names_vec Character vector of pathway names aligned to `ids`, or NULL
#'   when the source table carries no name column.
#' @param max_chars Maximum label length before the name is truncated.
#' @param append_id Append " (<accession>)" after the name. The heatmap needs
#'   this: it both keeps the accession visible and guarantees the unique
#'   rownames `pheatmap` requires.
#' @return Character vector of labels, one per element of `ids`, falling back to
#'   the bare accession wherever the name is missing or blank.
format_pathway_labels <- function(ids, names_vec = NULL, max_chars = 50,
                                  append_id = FALSE) {
    ids <- as.character(ids)
    labels <- if (is.null(names_vec)) {
        rep(NA_character_, length(ids))
    } else {
        as.character(names_vec)
    }

    unnamed <- is.na(labels) | !nzchar(trimws(labels))
    labels[unnamed] <- ids[unnamed]

    labels <- ifelse(nchar(labels) > max_chars,
                     paste0(substr(labels, 1, max_chars - 3), "..."),
                     labels)

    if (append_id) labels <- paste0(labels, " (", ids, ")")

    labels
}


#' Pull the pathway-name column out of an enrichment table, if it has one
#'
#' @param df Enrichment data frame.
#' @return Character vector of names aligned to `df` rows, or NULL.
pathway_name_column <- function(df) {
    if (is.data.frame(df) && "pathway_name" %in% names(df)) {
        as.character(df$pathway_name)
    } else {
        NULL
    }
}


#' Read the significance statistic the cross-omics figures should use
#'
#' Historically every cross-omics figure plotted the raw p-value even where the
#' per-omic table carried a BH-adjusted one, which over-states how many pathways
#' clear the 0.05 line. The statistic is now a config choice so a run cannot
#' drift between figures.
#'
#' @param config Full pipeline config; the key read is
#'   `modes.multiomics.enrichment.pvalue_type`.
#' @return "padj" (the default) or "pvalue".
resolve_pvalue_type <- function(config) {
    requested <- config$modes$multiomics$enrichment$pvalue_type
    if (is.null(requested) || !nzchar(as.character(requested)[1])) return("padj")
    requested <- tolower(trimws(as.character(requested)[1]))
    if (requested %in% c("padj", "pvalue")) return(requested)
    warning("Unknown modes.multiomics.enrichment.pvalue_type '", requested,
            "'; expected \"padj\" or \"pvalue\". Using \"padj\".")
    "padj"
}


#' Pick the significance column of an enrichment table
#'
#' Enrichment tables reach the cross-omics layer from four different producers
#' (fgsea, clusterProfiler ORA, the Fisher fallback, compound ORA) and spell the
#' same two quantities differently. This resolves the requested statistic and
#' reports what it actually found, so a fallback can be shown on the figure
#' rather than silently swapped in.
#'
#' @param df Enrichment data frame.
#' @param pvalue_type "padj" or "pvalue" — the statistic asked for.
#' @return `NULL` when the table has neither, otherwise a list with `column`
#'   (the column name) and `type` ("padj" or "pvalue" — what that column holds).
pick_pvalue_column <- function(df, pvalue_type = "padj") {
    if (!is.data.frame(df)) return(NULL)
    adj_cols <- c("padj", "p.adjust", "FDR", "qvalue")
    raw_cols <- c("pvalue", "pval", "PValue", "p_value")

    first_present <- function(candidates) {
        hit <- candidates[candidates %in% names(df)]
        if (length(hit) == 0) return(NULL)
        # A column of all-NA is present but unusable, and the merge downstream
        # cannot tell the difference between that and "not tested".
        for (h in hit) if (any(!is.na(df[[h]]))) return(h)
        NULL
    }

    order_pref <- if (identical(pvalue_type, "pvalue")) {
        list(c(raw_cols, "pvalue"), adj_cols)
    } else {
        list(adj_cols, raw_cols)
    }
    wanted <- if (identical(pvalue_type, "pvalue")) "pvalue" else "padj"
    other  <- if (identical(pvalue_type, "pvalue")) "padj" else "pvalue"

    col <- first_present(order_pref[[1]])
    if (!is.null(col)) return(list(column = col, type = wanted))
    col <- first_present(order_pref[[2]])
    if (!is.null(col)) return(list(column = col, type = other))
    NULL
}


#' Describe the statistic on an axis, naming any layer that fell back
#'
#' @param pvalue_type The statistic requested ("padj" or "pvalue").
#' @param types_used Named character vector: layer -> statistic actually used.
#' @return A one-line label such as "-log10(adjusted p)" or
#'   "-log10(adjusted p); raw p for: transcriptomics".
describe_pvalue_axis <- function(pvalue_type = "padj", types_used = NULL) {
    base <- if (identical(pvalue_type, "pvalue")) "-log10(p-value)" else "-log10(adjusted p)"
    if (is.null(types_used) || length(types_used) == 0) return(base)
    fallback <- names(types_used)[types_used != pvalue_type]
    if (length(fallback) == 0) return(base)
    other <- if (identical(pvalue_type, "pvalue")) "adjusted p" else "raw p"
    paste0(base, "; ", other, " for: ", paste(fallback, collapse = ", "))
}


#' Merge pathway significance values from multiple omics
#'
#' Column names stay `pval_<omics>` whichever statistic is selected: they are the
#' contract `select_multi_omics_pathways()`, the plots and the report all key on.
#' The statistic actually taken per layer is returned as the `pvalue_type_used`
#' attribute so the figures can say so.
#'
#' @param pathway_tables Named list of per-omics enrichment data frames.
#' @param target_pathways Pathway accessions to keep.
#' @param omics Character vector of omics layer names, in column order.
#' @param pvalue_type "padj" (default) or "pvalue".
#' @return Data frame with one `pval_<omics>` column per layer, carrying a
#'   `pvalue_type_used` attribute (named character vector, layer -> statistic).
merge_pathway_pvalues <- function(pathway_tables, target_pathways, omics,
                                  pvalue_type = "padj") {

    merged <- data.frame(pathway = target_pathways, stringsAsFactors = FALSE)
    types_used <- character(0)

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

        picked <- pick_pvalue_column(df, pvalue_type)

        if (is.null(picked)) {
            warning("Cannot identify a significance column in ", om,
                    " enrichment table")
            next
        }
        if (picked$type != pvalue_type) {
            message("  ", om, ": no usable ", pvalue_type, " column, using ",
                    picked$column, " instead")
        }
        types_used[om] <- picked$type

        # For multiple contrasts, take the most significant value per pathway
        df_agg <- aggregate(
            stats::as.formula(paste(picked$column, "~", pathway_col)),
            data = df,
            FUN = min
        )
        colnames(df_agg) <- c("pathway", paste0("pval_", om))

        # Subset to target pathways
        df_sub <- df_agg[df_agg$pathway %in% target_pathways, , drop = FALSE]

        # Merge
        merged <- merge(merged, df_sub, by = "pathway", all.x = TRUE)
    }

    attr(merged, "pvalue_type_used") <- types_used
    merged
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

#' Build a per-layer NES matrix for a set of pathways
#'
#' Pulls the signed GSEA normalised enrichment score out of each layer's
#' enrichment table. A layer may hold several rows per pathway (one per contrast,
#' database or method), so the NES is taken from the most significant row that
#' has one — the same "best row wins" rule `merge_pathway_pvalues()` applies to
#' the p-values, so a cell's colour and its p-value describe the same test.
#'
#' @param pathway_tables Named list of per-omics enrichment data frames.
#' @param pathways Pathway accessions, in the row order wanted.
#' @param omic_names Layer names, in the column order wanted.
#' @param pvalue_type "padj" or "pvalue" — which statistic ranks the rows.
#' @return Numeric matrix (pathways x layers) of NES values with NA where the
#'   pathway was not scored, or NULL when no layer reports NES at all.
build_layer_nes_matrix <- function(pathway_tables, pathways, omic_names,
                                   pvalue_type = "padj") {
    if (is.null(pathway_tables) || length(pathways) == 0) return(NULL)

    nes_mat <- matrix(NA_real_, nrow = length(pathways), ncol = length(omic_names),
                      dimnames = list(pathways, omic_names))
    any_nes <- FALSE

    for (om in omic_names) {
        df <- pathway_tables[[om]]
        if (!is.data.frame(df) || !("NES" %in% names(df))) next

        pathway_col <- if ("pathway" %in% names(df)) "pathway"
                       else if ("ID" %in% names(df)) "ID"
                       else if ("Description" %in% names(df)) "Description"
                       else NULL
        if (is.null(pathway_col)) next

        keep <- !is.na(df$NES)
        if (!any(keep)) next
        sub <- df[keep, , drop = FALSE]

        picked <- pick_pvalue_column(sub, pvalue_type)
        rank_by <- if (is.null(picked)) rep(0, nrow(sub)) else
            suppressWarnings(as.numeric(sub[[picked$column]]))
        rank_by[is.na(rank_by)] <- Inf

        # order() then !duplicated() keeps exactly one row per pathway, the one
        # with the smallest p; ties resolve on the table's own row order so the
        # figure is reproducible across runs.
        sub <- sub[order(rank_by), , drop = FALSE]
        sub <- sub[!duplicated(as.character(sub[[pathway_col]])), , drop = FALSE]

        hit <- match(pathways, as.character(sub[[pathway_col]]))
        nes_mat[, om] <- suppressWarnings(as.numeric(sub$NES[hit]))
        any_nes <- TRUE
    }

    if (!any_nes || all(is.na(nes_mat))) return(NULL)
    nes_mat
}


#' Plot cross-omics pathway heatmap
#'
#' Rows are the pathways seen by the most omics layers first, then the most
#' significant — not simply the most significant, which would fill the figure
#' with single-layer rows and leave the other columns blank.
#'
#' Cells show the signed GSEA NES: red positive, blue negative, white near zero,
#' on a scale symmetric about zero. Over-representation layers (compound ORA,
#' Fisher ORA) have no NES at all; such a column would otherwise be a block of
#' grey reading as "not tested", so its header is marked "(no NES)" and the
#' caption says what that grey means. With no NES anywhere the figure falls back
#' to the unsigned -log10 scale it used before.
#'
#' @param meta_results Meta-analysis data frame from `stouffer_combined_pvalues()`,
#'   optionally carrying `pathway_name`.
#' @param omics Character vector of omics layer names (kept for signature
#'   compatibility with the other cross-omics plotters).
#' @param top_n Number of pathways to show.
#' @param pathway_tables Named list of per-omics enrichment data frames, the
#'   source of the NES values. NULL keeps the -log10 scale.
#' @param value_label Axis/title wording for the significance statistic, from
#'   `describe_pvalue_axis()`.
#' @param pvalue_type "padj" or "pvalue" — which statistic picks a pathway's NES
#'   when a layer scored it more than once.
#' @return Invisibly NULL; called for the plot it draws on the active device.
plot_cross_omics_pathway_heatmap <- function(meta_results, omics, top_n = 30,
                                             pathway_tables = NULL,
                                             value_label = "-log10(p-value)",
                                             pvalue_type = "padj",
                                             collection = NULL) {

    top_pathways <- select_multi_omics_pathways(meta_results, top_n)

    pval_cols <- grep("^pval_", names(top_pathways), value = TRUE)
    omic_names <- sub("^pval_", "", pval_cols)
    pval_matrix <- as.matrix(top_pathways[, pval_cols, drop = FALSE])

    pathway_labels <- format_pathway_labels(
        top_pathways$pathway,
        pathway_name_column(top_pathways),
        append_id = TRUE
    )
    rownames(pval_matrix) <- pathway_labels
    colnames(pval_matrix) <- omic_names

    nes_matrix <- build_layer_nes_matrix(pathway_tables,
                                         as.character(top_pathways$pathway),
                                         omic_names, pvalue_type)
    no_nes <- if (is.null(nes_matrix)) omic_names else
        omic_names[colSums(!is.na(nes_matrix)) == 0]

    # pheatmap draws `main` at full size and clips it against the device rather
    # than wrapping, so every title line below is kept short deliberately.
    stat_word <- if (identical(pvalue_type, "pvalue")) "p" else "adjusted p"
    rank_note <- paste0("rows ranked on layers with ", stat_word,
                        " < 0.05, then combined p")
    col_labels <- omic_names

    if (!is.null(nes_matrix)) {
        plot_matrix <- nes_matrix
        rownames(plot_matrix) <- pathway_labels
        # A layer that never reports NES is grey top to bottom, which on its own
        # reads as "not tested". Marking its header is what keeps the two kinds
        # of grey apart; a bare asterisk rather than a word because pheatmap
        # clips a long rotated label of the leftmost column against the device.
        marked <- omic_names %in% no_nes
        col_labels[marked] <- paste0(omic_names[marked], " *")
        # Symmetric limits so zero sits on the white midpoint of the ramp and a
        # +1.8 in one layer is the same colour distance from white as a -1.8.
        lim <- max(abs(plot_matrix), na.rm = TRUE)
        if (!is.finite(lim) || lim <= 0) lim <- 1
        # Mid-tone ends rather than the full RdBu extremes: the per-cell NES has
        # to stay readable in black on top of the darkest colour in the ramp.
        heat_colors <- colorRampPalette(c("#4393C3", "#F7F7F7", "#D6604D"))(51)
        heat_breaks <- seq(-lim, lim, length.out = 52)
        cell_labels <- matrix(formatC(plot_matrix, format = "f", digits = 2),
                              nrow = nrow(plot_matrix),
                              dimnames = dimnames(plot_matrix))
        cell_labels[is.na(plot_matrix)] <- ""
        grey_note <- if (length(no_nes) == 0) {
            "grey = pathway not scored in that layer"
        } else {
            "grey = not scored; * = over-representation layer, no NES"
        }
        heatmap_title <- paste0(
        if (!is.null(collection)) paste0(collection, " gene sets\n"),
            "Cross-Omics Pathway Enrichment (GSEA NES)\n",
            "red = positive, blue = negative\n",
            rank_note, "\n",
            grey_note
        )
    } else {
        plot_matrix <- -log10(pval_matrix + 1e-300)
        plot_matrix[] <- pmin(plot_matrix, 10)  # pmin keeps NA as NA
        heat_colors <- colorRampPalette(c("white", "gold", "orange", "red"))(50)
        heat_breaks <- NULL
        cell_labels <- FALSE
        heatmap_title <- paste0(
        if (!is.null(collection)) paste0(collection, " gene sets\n"),
            "Cross-Omics Pathway Enrichment\n",
            value_label, "; no GSEA NES in any layer\n",
            rank_note, "\n",
            "grey = not tested in that layer"
        )
    }

    # Heatmap. NAs stay NA so na_col separates "not tested in this layer" from
    # "tested, no effect" — coercing them to 0 painted both like a real result.
    if (requireNamespace("pheatmap", quietly = TRUE)) {
        args <- list(plot_matrix,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     main = heatmap_title,
                     color = heat_colors,
                     labels_col = col_labels,
                     fontsize = 9,
                     fontsize_row = 7, fontsize_col = 10,
                     fontsize_number = 6, number_color = "black",
                     angle_col = 45,
                     na_col = "grey90",
                     border_color = "grey80",
                     display_numbers = cell_labels)
        if (!is.null(heat_breaks)) args$breaks <- heat_breaks
        do.call(pheatmap::pheatmap, args)
    } else {
        # stats::heatmap() has no na_col, so the fallback keeps the old coercion
        plot_matrix[is.na(plot_matrix)] <- 0
        heatmap(plot_matrix, scale = "none", Colv = NA,
                main = heatmap_title, col = heat_colors)
    }

    invisible(NULL)
}


#' Pull the pathway-accession column out of an enrichment table, if it has one
#'
#' The four enrichment producers that feed this section (fgsea, clusterProfiler
#' ORA, the Fisher fallback, compound ORA) each spell the accession column
#' differently, so the lookup is centralised here rather than repeated at every
#' call site.
#'
#' @param df Enrichment data frame.
#' @return Name of the accession column, or NULL when the table has none.
pathway_id_column <- function(df) {
    if (!is.data.frame(df)) return(NULL)
    for (candidate in c("pathway", "ID", "Description")) {
        if (candidate %in% names(df)) return(candidate)
    }
    NULL
}


#' Keep only the over-representation rows of an enrichment table
#'
#' Most layers write a `method` column and the filter is literal. The compound
#' ORA producer does not: its table carries neither `method` nor `NES`, and
#' dropping it on a missing label would hide the one layer this figure exists
#' for. A table with no NES column cannot hold GSEA output, so it is taken as
#' over-representation throughout.
#'
#' @param df Enrichment data frame, or NULL.
#' @return `df` restricted to its ORA rows (possibly zero rows); NULL stays NULL.
filter_ora_rows <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(df)

    if ("method" %in% names(df)) {
        return(df[tolower(trimws(as.character(df$method))) %in% "ora", , drop = FALSE])
    }
    if ("NES" %in% names(df)) return(df[0, , drop = FALSE])
    df
}


#' Build a per-layer significance matrix from over-representation rows only
#'
#' A layer can report the same pathway several times (one row per contrast,
#' database, or ORA direction); the smallest value wins, the same "best row
#' wins" rule `build_layer_nes_matrix()` and `merge_pathway_pvalues()` apply, so
#' a cell means the same thing in all three places.
#'
#' Every requested layer gets a column even when it contributed nothing, so the
#' caller can tell "layer ran no ORA" from "pathway not tested in this layer" —
#' the two look identical once the matrix has been collapsed.
#'
#' @param ora_tables Named list of per-omics enrichment data frames already
#'   restricted to ORA rows by `filter_ora_rows()`.
#' @param pathways Pathway accessions, in the row order wanted.
#' @param omic_names Layer names, in the column order wanted.
#' @param pvalue_type "padj" (default) or "pvalue" — the statistic requested.
#' @return Numeric matrix (pathways x layers) of significance values, NA where
#'   the pathway was not tested in that layer, carrying a `pvalue_type_used`
#'   attribute (named character vector, layer -> statistic actually taken, NA
#'   for a layer with no usable ORA rows).
build_layer_ora_matrix <- function(ora_tables, pathways, omic_names,
                                   pvalue_type = "padj") {
    pmat <- matrix(NA_real_, nrow = length(pathways), ncol = length(omic_names),
                   dimnames = list(pathways, omic_names))
    types_used <- stats::setNames(rep(NA_character_, length(omic_names)), omic_names)

    for (om in omic_names) {
        df <- ora_tables[[om]]
        if (!is.data.frame(df) || nrow(df) == 0) next

        pathway_col <- pathway_id_column(df)
        if (is.null(pathway_col)) next

        picked <- pick_pvalue_column(df, pvalue_type)
        if (is.null(picked)) next
        types_used[om] <- picked$type

        vals <- suppressWarnings(as.numeric(df[[picked$column]]))
        keep <- !is.na(vals)
        if (!any(keep)) next

        best <- tapply(vals[keep], as.character(df[[pathway_col]])[keep], min)
        pmat[, om] <- as.numeric(best[match(pathways, names(best))])
    }

    attr(pmat, "pvalue_type_used") <- types_used
    pmat
}


#' Plot the cross-omics over-representation heatmap
#'
#' Companion to `plot_cross_omics_pathway_heatmap()`, built from the ORA rows
#' alone. The NES figure can only colour a layer that ran GSEA, so a layer whose
#' enrichment is over-representation end to end (compound ORA on metabolites,
#' the Fisher fallback on a non-model organism) is structurally absent from it.
#' Here every layer is on the same footing because -log10(FDR) exists for all of
#' them.
#'
#' Rows are ranked by `select_multi_omics_pathways()`, as in the NES figure, but
#' over an ORA-only meta-analysis: ranking on the shared table would fill the
#' figure with terms that only GSEA ever scored, leaving little but grey.
#'
#' @param meta_results Meta-analysis data frame from `stouffer_combined_pvalues()`
#'   (used only to name layers when `pathway_tables` is NULL).
#' @param omics Character vector of omics layer names, in column order.
#' @param top_n Number of pathways to show.
#' @param pathway_tables Named list of per-omics enrichment data frames — the
#'   source of the ORA rows. Without it there is nothing to draw.
#' @param pvalue_type "padj" (default) or "pvalue". The figure is defined on the
#'   FDR; the argument exists so a layer that carries no usable adjusted p can
#'   still be plotted on its raw p, which the title then says.
#' @return Invisibly NULL; called for the plot it draws on the active device.
plot_cross_omics_ora_heatmap <- function(meta_results, omics, top_n = 30,
                                         pathway_tables = NULL,
                                         pvalue_type = "padj",
                                         collection = NULL) {

    omic_names <- if (!is.null(pathway_tables) && length(pathway_tables) > 0) {
        names(pathway_tables)
    } else if (!is.null(meta_results)) {
        sub("^pval_", "", grep("^pval_", names(meta_results), value = TRUE))
    } else {
        omics
    }

    draw_note <- function(msg) {
        plot.new()
        text(0.5, 0.5, msg, cex = 1.1)
        invisible(NULL)
    }

    if (length(omic_names) == 0 || is.null(pathway_tables)) {
        return(draw_note("No enrichment tables available for the ORA heatmap"))
    }

    ora_tables <- lapply(pathway_tables[omic_names], filter_ora_rows)
    names(ora_tables) <- omic_names

    # This figure rebuilds its own ORA-only ranking from pathway_tables rather
    # than reading meta_results, so a caller that subsets meta_results to one
    # gene-set collection has no effect unless the same subset is applied here.
    # Without this every per-collection file came out byte-identical, each one
    # showing GO, KEGG and Pfam rows together.
    if (!is.null(collection) && exists("classify_pathway_collection")) {
        ora_tables <- lapply(ora_tables, function(d) {
            if (!is.data.frame(d) || nrow(d) == 0) return(d)
            idc <- pathway_id_column(d)
            if (is.na(idc)) return(d)
            d[classify_pathway_collection(as.character(d[[idc]])) %in% collection, ,
              drop = FALSE]
        })
    }
    n_ora <- vapply(ora_tables,
                    function(d) if (is.data.frame(d)) nrow(d) else 0L, integer(1))

    if (all(n_ora == 0)) {
        return(draw_note("No over-representation results in any omics layer"))
    }

    populated <- ora_tables[n_ora > 0]
    ora_pathways <- unique(unlist(lapply(populated, function(d) {
        col <- pathway_id_column(d)
        if (is.null(col)) NULL else as.character(d[[col]])
    })))
    if (length(ora_pathways) == 0) {
        return(draw_note("No over-representation results in any omics layer"))
    }

    ora_meta <- merge_pathway_pvalues(populated, ora_pathways, names(populated),
                                      pvalue_type = pvalue_type)
    if (length(grep("^pval_", names(ora_meta))) == 0) {
        return(draw_note("No usable significance column in the ORA results"))
    }
    ora_meta <- stouffer_combined_pvalues(ora_meta)
    ora_meta <- attach_pathway_names_from_tables(ora_meta, populated)

    top_pathways <- select_multi_omics_pathways(ora_meta, top_n)
    if (is.null(top_pathways) || nrow(top_pathways) == 0) {
        return(draw_note("No over-representation results to show"))
    }

    pmat <- build_layer_ora_matrix(ora_tables, as.character(top_pathways$pathway),
                                   omic_names, pvalue_type)
    types_used <- attr(pmat, "pvalue_type_used")

    # 1e-300 guards -log10(0) the same way the -log10 fallback of the NES figure
    # does; it is far below any FDR a real test returns.
    score <- -log10(pmax(pmat, 1e-300))
    # -log10(1) is IEEE negative zero and formatC prints that as "-0.00", which
    # reads like a sign on a magnitude scale. Re-assign it as a plain zero.
    score[!is.na(score) & score == 0] <- 0
    rownames(score) <- format_pathway_labels(
        top_pathways$pathway,
        pathway_name_column(top_pathways),
        append_id = TRUE
    )

    cell_labels <- matrix(formatC(score, format = "f", digits = 2),
                          nrow = nrow(score), dimnames = dimnames(score))
    cell_labels[is.na(score)] <- ""

    # Two ways a column ends up blank, and they mean different things: the layer
    # ran no ORA at all, or it did but none of its pathways reached these rows.
    no_ora <- omic_names[n_ora == 0]
    not_shown <- omic_names[n_ora > 0 & colSums(!is.na(score)) == 0]
    col_labels <- omic_names
    col_labels[omic_names %in% no_ora] <- paste0(no_ora, " *")
    col_labels[omic_names %in% not_shown] <- paste0(not_shown, " +")

    # pheatmap draws `main` at full size and clips it against the device rather
    # than wrapping, so every line below is kept to roughly the width of the
    # figure's own title -- the markers are terse for that reason, not for style.
    marker_note <- c(
        if (length(no_ora) > 0) "* = no ORA",
        if (length(not_shown) > 0) "+ = ORA run, none shown"
    )

    stat_word <- if (identical(pvalue_type, "pvalue")) "p-value" else "FDR"
    fell_back <- names(types_used)[!is.na(types_used) & types_used != pvalue_type]
    fallback_line <- if (length(fell_back) == 0) NULL else {
        other <- if (identical(pvalue_type, "pvalue")) "adjusted p" else "raw p"
        paste0(other, " for: ", paste(fell_back, collapse = ", "))
    }

    heatmap_title <- paste(c(
        "Cross-Omics Pathway Over-Representation (ORA only)",
        paste0("colour = -log10(", stat_word, "); 1.3 = ", stat_word, " 0.05"),
        fallback_line,
        paste0("rows ranked on layers with ", stat_word,
               " < 0.05, then combined p"),
        paste(c("grey = not tested", marker_note), collapse = "; ")
    ), collapse = "\n")

    # Sequential ramp: this is a magnitude, not a direction, so it runs pale pink
    # to dark red with no meaningful midpoint. The dark end stops short of black
    # so the per-cell number stays readable on top of it.
    heat_colors <- colorRampPalette(c("#FFF0F2", "#FBC4C9", "#F08D8D",
                                      "#DD4B3E", "#A50F15"))(51)
    lim <- suppressWarnings(max(score, na.rm = TRUE))
    if (!is.finite(lim) || lim <= 0) lim <- 1
    heat_breaks <- seq(0, lim, length.out = 52)

    # NAs stay NA so na_col separates "not tested in this layer" from "tested,
    # FDR near 1" — coercing them to 0 paints both the palest pink.
    if (requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(score,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           main = heatmap_title,
                           color = heat_colors,
                           breaks = heat_breaks,
                           labels_col = col_labels,
                           fontsize = 9,
                           fontsize_row = 7, fontsize_col = 10,
                           fontsize_number = 6, number_color = "black",
                           angle_col = 45,
                           na_col = "grey90",
                           border_color = "grey80",
                           display_numbers = cell_labels)
    } else {
        # stats::heatmap() has no na_col, so the fallback keeps the old coercion
        score[is.na(score)] <- 0
        heatmap(score, scale = "none", Colv = NA,
                main = heatmap_title, col = heat_colors)
    }

    invisible(NULL)
}


#' Plot enrichment dot plot
#'
#' Shares the row selection of `plot_cross_omics_pathway_heatmap()`: layers
#' first, significance second, so the figure shows where the omics meet.
#'
#' @param meta_results Meta-analysis data frame from `stouffer_combined_pvalues()`,
#'   optionally carrying `pathway_name`.
#' @param omics Character vector of omics layer names (kept for signature
#'   compatibility with the other cross-omics plotters).
#' @param top_n Number of pathways to show.
#' @param value_label Axis wording for the significance statistic, from
#'   `describe_pvalue_axis()`.
#' @return Invisibly NULL; called for the plot it draws on the active device.
plot_enrichment_dotplot <- function(meta_results, omics, top_n = 20,
                                    value_label = "-log10(p-value)") {

    top <- select_multi_omics_pathways(meta_results, top_n)

    pval_cols <- grep("^pval_", names(top), value = TRUE)

    # One label per row, computed once so the long-format frame and the axis
    # ordering below cannot drift apart.
    # make.unique(): two accessions can share a term name, and a duplicated
    # label would collapse two rows onto one axis position.
    top$plot_label <- make.unique(format_pathway_labels(
        top$pathway, pathway_name_column(top), max_chars = 45
    ))

    # Build long-format data
    plot_data <- list()
    for (pc in pval_cols) {
        om_name <- sub("^pval_", "", pc)
        df <- data.frame(
            pathway = top$plot_label,
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

    # Reverse pathway order for bottom-to-top display
    pathway_order <- rev(unique(top$plot_label))
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
         xlab = value_label, ylab = "",
         yaxt = "n", main = "Cross-Omics Pathway Enrichment")
    mtext("ordered by number of omics layers, then combined p-value",
          side = 3, line = 0.2, cex = 0.7)
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
#'
#' @param pathway_table Per-omics enrichment data frame, optionally carrying
#'   `pathway_name`.
#' @param omics_name Name of the omics layer, used for the title and bar colour.
#' @param top_n Number of pathways to show.
#' @param pvalue_type "padj" (default) or "pvalue". When the table has no usable
#'   column of the requested kind the other one is plotted and the axis label
#'   says so, so the bars are never ambiguous.
#' @return Invisibly NULL; called for the plot it draws on the active device.
plot_per_omics_barplot <- function(pathway_table, omics_name, top_n = 15,
                                   pvalue_type = "padj") {

    picked <- pick_pvalue_column(pathway_table, pvalue_type)
    pval_col <- picked$column

    pathway_col <- if ("pathway" %in% names(pathway_table)) "pathway"
                   else if ("Description" %in% names(pathway_table)) "Description"
                   else if ("ID" %in% names(pathway_table)) "ID"
                   else NULL

    if (is.null(pval_col) || is.null(pathway_col)) {
        plot.new()
        text(0.5, 0.5, paste("Cannot identify columns for", omics_name))
        return(invisible(NULL))
    }

    x_label <- describe_pvalue_axis(
        pvalue_type,
        stats::setNames(picked$type, omics_name)
    )

    # Aggregate across contrasts: take the most significant value per pathway
    agg <- aggregate(
        stats::as.formula(paste(pval_col, "~", pathway_col)),
        data = pathway_table,
        FUN = min
    )
    colnames(agg) <- c("pathway", "pvalue")
    agg <- agg[order(agg$pvalue), ]
    agg <- agg[seq_len(min(top_n, nrow(agg))), ]

    # aggregate() drops every column but the two it keyed on, so re-join the
    # names rather than labelling from the accession.
    names_vec <- NULL
    all_names <- pathway_name_column(pathway_table)
    if (!is.null(all_names)) {
        lookup <- stats::setNames(all_names, as.character(pathway_table[[pathway_col]]))
        lookup <- lookup[!duplicated(names(lookup))]
        names_vec <- unname(lookup[as.character(agg$pathway)])
    }
    agg$label <- make.unique(format_pathway_labels(agg$pathway, names_vec, max_chars = 45))

    neg_log_p <- -log10(agg$pvalue + 1e-300)
    neg_log_p <- pmin(neg_log_p, 15)

    omics_colors <- c(
        transcriptomics = "#E41A1C",
        proteomics = "#377EB8",
        metabolomics = "#4DAF4A"
    )
    bar_col <- if (omics_name %in% names(omics_colors)) omics_colors[omics_name] else "steelblue"

    # Anchor the axis so the 0.05 reference line is always on scale. Under
    # adjusted p a layer can have nothing significant at all, and barplot()'s
    # automatic range would then draw a bare -1..1 axis with the reference line
    # off-screen — a reader could not tell "nothing survives correction" from a
    # broken figure.
    x_max <- max(c(neg_log_p, -log10(0.05) * 1.15), na.rm = TRUE)

    par(mar = c(5, 15, 3, 2))
    barplot(rev(neg_log_p), horiz = TRUE, names.arg = rev(agg$label),
            las = 1, cex.names = 0.7, col = bar_col,
            xlab = x_label, xlim = c(0, x_max),
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
    excl_classes <- config$modes$multiomics$enrichment$exclude_pathway_classes

    results <- list()

    # --- DIABLO loadings enrichment ---
    if (!is.null(integration_res$diablo_results)) {
        message("Running enrichment on DIABLO loadings...")
        diablo_dir <- file.path(out_dir, "diablo_loadings")
        dir.create(diablo_dir, showWarnings = FALSE)

        results$diablo <- tryCatch(
            run_diablo_loadings_enrichment(
                exclude_classes = excl_classes,
                diablo_results = integration_res$diablo_results,
                harmonization_res = harmonization_res,
                organism = organism,
                kegg_org = kegg_org,
                org_db = org_db,
                config = config,
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
                exclude_classes = excl_classes,
                mofa_results = integration_res$mofa_results,
                harmonization_res = harmonization_res,
                organism = organism,
                kegg_org = kegg_org,
                org_db = org_db,
                config = config,
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
#'
#' @param diablo_results DIABLO integration result carrying `$top_features`.
#' @param harmonization_res Harmonization result (ID resolution + universe).
#' @param organism Organism name from the config.
#' @param kegg_org KEGG organism code, or NULL when the organism has none.
#' @param org_db OrgDb package name, or NULL when the organism has none.
#' @param config Full pipeline config; needed for the per-omic custom GMTs that
#'   gene views fall back to when there is no KEGG organism / OrgDb.
#' @param out_dir Output directory for CSVs and barplots.
#' @param top_n Number of top-loading features per component to enrich.
#' @param exclude_classes KEGG pathway classes to drop (compound ORA only).
#' @return Row-bound data frame of enriched pathways over all components, or NULL.
run_diablo_loadings_enrichment <- function(diablo_results, harmonization_res,
                                            organism, kegg_org, org_db,
                                            out_dir, top_n = 50,
                                            exclude_classes = NULL,
                                            config = NULL) {

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
                org_db = org_db,
                config = config
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
#'
#' @param mofa_results MOFA2 integration result carrying `$weights` per view.
#' @param harmonization_res Harmonization result (ID resolution + universe).
#' @param organism Organism name from the config.
#' @param kegg_org KEGG organism code, or NULL when the organism has none.
#' @param org_db OrgDb package name, or NULL when the organism has none.
#' @param config Full pipeline config; needed for the per-omic custom GMTs that
#'   gene views fall back to when there is no KEGG organism / OrgDb.
#' @param out_dir Output directory for CSVs and barplots.
#' @param top_n Number of top-weight features per factor to enrich.
#' @param exclude_classes KEGG pathway classes to drop (compound ORA only).
#' @return Row-bound data frame of enriched pathways over all factors, or NULL.
run_mofa_weights_enrichment <- function(mofa_results, harmonization_res,
                                         organism, kegg_org, org_db,
                                         out_dir, top_n = 50,
                                         exclude_classes = NULL,
                                         config = NULL) {

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
                org_db = org_db,
                config = config
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
                         exclude_classes = exclude_classes,
                         label = label),
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
#' Maps feature IDs to ENTREZ IDs and runs KEGG ORA via clusterProfiler. Handles
#' GENE_N synthetic IDs from the harmonized MAE by translating them back to the
#' native gene / protein namespace via the gene_protein_mapping table.
#'
#' @param feature_ids Top-loading feature IDs of one view (GENE_N or native).
#' @param omics_type One of "transcriptomics" / "proteomics".
#' @param harmonization_res Harmonization result (ID resolution + universe).
#' @param organism Organism name from the config.
#' @param kegg_org KEGG organism code, or NULL when the organism has none.
#' @param org_db OrgDb package name, or NULL when the organism has none.
#' @param config Full pipeline config, used only by the GMT fallback below.
#' @return data.frame(pathway, ID, pvalue, padj, GeneRatio, setSize), or NULL.
enrich_feature_list <- function(feature_ids, omics_type, harmonization_res,
                                 organism, kegg_org, org_db, config = NULL) {

    # Resolve IDs using the actual omics type. Needed by both branches below,
    # so it happens before the KEGG/OrgDb decision rather than after it.
    resolved_ids <- resolve_gene_n_ids(feature_ids, harmonization_res, omics_type)

    # Non-model organism: clusterProfiler has neither a KEGG code nor an OrgDb
    # to map against, so the whole gene view used to return NULL and the
    # loadings-enrichment folder came out metabolomics-only. The per-omic custom
    # GMTs are the same gene sets the DE-driven enrichment already uses, so fall
    # back to those instead of skipping the view.
    if (is.null(kegg_org) || is.null(org_db)) {
        return(enrich_feature_list_gmt(resolved_ids, omics_type,
                                       harmonization_res, config))
    }

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


#' GMT-based ORA over a loadings feature list (non-model fallback)
#'
#' Over-representation of one view's top-loading features against that omic's
#' custom gene sets (\code{modes.<omic>.pathway.gmt_file}), for organisms with
#' no KEGG code / OrgDb. Mirrors \code{run_multi_ora_gmt()}: same path
#' resolution, same GMT reader, same \code{clusterProfiler::enricher} wrapper,
#' so a pathway reported here means the same thing as one reported by the
#' DE-driven Multi-ORA.
#'
#' @param resolved_ids Native-namespace feature IDs of the top loadings
#'   (already through \code{resolve_gene_n_ids()}).
#' @param omics_type One of "transcriptomics" / "proteomics"; other views have
#'   no gene-set collection here and return NULL.
#' @param harmonization_res Harmonization result; its per-omic \code{expr_work}
#'   supplies the tested-feature background.
#' @param config Full pipeline config, for the per-omic \code{gmt_file}.
#' @return data.frame(pathway, ID, pvalue, padj, GeneRatio, setSize), or NULL.
enrich_feature_list_gmt <- function(resolved_ids, omics_type,
                                    harmonization_res, config) {
    if (is.null(config)) return(NULL)

    cfg_key <- c(transcriptomics = "rna", proteomics = "proteomics")[omics_type]
    if (is.na(cfg_key)) return(NULL)

    # gmt_file may be a single path or a YAML list (GO + KEGG + PFAM);
    # resolve_input_path() vectorises and leaves absolute paths alone, and
    # read_gmt() (via gmt_to_term2gene) merges several files into one collection.
    gmt_path <- unlist(config$modes[[cfg_key]]$pathway$gmt_file, use.names = FALSE)
    if (length(gmt_path) == 0 || !any(nzchar(gmt_path))) return(NULL)
    gmt_abs <- resolve_input_path(config, gmt_path)
    if (any(!file.exists(gmt_abs))) {
        message("    Loadings ORA (GMT): ", omics_type, " gmt_file not found: ",
                paste(gmt_abs[!file.exists(gmt_abs)], collapse = ", "))
        return(NULL)
    }

    gs <- gmt_to_term2gene(gmt_abs)
    if (is.null(gs) || nrow(gs$t2g) == 0) return(NULL)

    # Background = every feature of this view that survived preprocessing, in the
    # same namespace as the query. Without it enricher() would silently use the
    # gene sets themselves as the universe and inflate every p-value.
    universe <- NULL
    pre_data <- harmonization_res$inputs[[omics_type]]
    if (!is.null(pre_data) && !is.null(pre_data$expr_work)) {
        universe <- unique(resolve_gene_n_ids(rownames(pre_data$expr_work),
                                              harmonization_res, omics_type))
    }

    sig_genes <- unique(resolved_ids[!is.na(resolved_ids)])
    message("    Loadings ORA (GMT) ", omics_type, ": ", length(sig_genes),
            " query / ", length(universe), " background features, ",
            length(unique(gs$t2g$term)), " gene sets")

    ora <- run_multi_ora_enricher(
        sig_genes = sig_genes,
        universe  = universe,
        term2gene = gs$t2g,
        term2name = gs$t2n,
        label     = paste0(omics_type, " loadings")
    )
    if (is.null(ora) || nrow(ora) == 0) return(NULL)

    # Same column shape as the KEGG branch above, so .rbind_fill() can stack
    # gene-view and compound-view results into one loadings table.
    data.frame(
        pathway   = ora$pathway,
        ID        = ora$ID,
        pvalue    = ora$pvalue,
        padj      = ora$padj,
        GeneRatio = ora$GeneRatio,
        setSize   = ora$Count,
        stringsAsFactors = FALSE
    )
}


#' Resolve GENE_N synthetic IDs to original feature IDs
#'
#' The harmonized MAE uses GENE_N IDs (where N = row in gene_protein_mapping).
#' This function translates them back to the native gene ID (transcriptomics) or
#' protein ID (proteomics). MOFA2 requires feature names to be unique across
#' views and appends the view name to the ones it finds in more than one view,
#' so "GENE_12" and "GENE_12_transcriptomics" both have to resolve — matching
#' only the bare form left 41% of the MOFA weights unresolved.
#'
#' @param feature_ids Character vector of feature IDs, GENE_N or native.
#' @param harmonization_res Harmonization result carrying `gene_protein_mapping`.
#' @param omics_type One of "transcriptomics" / "proteomics"; anything else is
#'   returned untouched (metabolomics does not use GENE_N).
#' @return Character vector of resolved IDs; entries with no mapping are dropped.
resolve_gene_n_ids <- function(feature_ids, harmonization_res, omics_type) {
    gpm <- harmonization_res$gene_protein_mapping
    if (is.null(gpm)) return(feature_ids)

    gene_n_pattern <- "^GENE_\\d+(_.+)?$"
    is_gene_n <- grepl(gene_n_pattern, feature_ids)
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
    keys <- sub("^(GENE_\\d+)(_.+)?$", "\\1", feature_ids[is_gene_n])

    resolved <- feature_ids
    resolved[is_gene_n] <- lookup[keys]
    resolved <- resolved[!is.na(resolved)]

    n_mapped <- sum(!is.na(lookup[keys]))
    message("    Resolved ", n_mapped, "/", sum(is_gene_n),
            " GENE_N IDs to ", id_col, " (", omics_type, ")")

    resolved
}


#' Plot barplot for loadings enrichment
#'
#' @param enrich_df Enrichment data frame with `pathway` and `pvalue` columns,
#'   optionally carrying `pathway_name`.
#' @param title Label for the plot title (component or factor name).
#' @param out_path PNG path to write.
#' @param top_n Number of pathways to show.
#' @return Invisibly NULL; called for the PNG it writes.
plot_loadings_enrichment_barplot <- function(enrich_df, title, out_path, top_n = 15) {
    if (nrow(enrich_df) == 0) return(invisible(NULL))

    df <- enrich_df[order(enrich_df$pvalue), ]
    df <- df[seq_len(min(top_n, nrow(df))), ]

    # The compound-ORA tables key `pathway` on the KEGG map id, so labelling
    # from it alone renders these barplots as a column of bare accessions.
    df$label <- make.unique(format_pathway_labels(
        df$pathway, pathway_name_column(df), max_chars = 45
    ))
    neg_log_p <- -log10(df$pvalue + 1e-300)
    neg_log_p <- pmin(neg_log_p, 15)

    png(out_path, width = 900, height = 600, res = 120)
    par(mar = c(5, 15, 3, 2))
    # The component/factor labels run long ("MOFA_transcriptomics_Factor1"), so
    # the title needs shrinking to survive the 900px device.
    barplot(rev(neg_log_p), horiz = TRUE, names.arg = rev(df$label),
            las = 1, cex.names = 0.65, col = "steelblue",
            xlab = "-log10(p-value)",
            main = paste("Loadings Enrichment:", title), cex.main = 0.9)
    abline(v = -log10(0.05), col = "red", lty = 2)
    dev.off()
    message("    Saved: ", out_path)
}
