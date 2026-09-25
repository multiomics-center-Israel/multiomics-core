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

    # .rbind_fill() rather than rbind(): use_method is chosen per contrast on the
    # mapped-feature count above, so one layer's contrasts can straddle the
    # threshold and return different schemas -- ORA carries Fold_enrichment and
    # Count, GSEA carries NES and core_enrichment. rbind() aborts on that, and
    # since the caller catches the error and warns, the whole layer would come
    # back NULL and read as "no enriched pathways" rather than as a failure.
    combined <- .rbind_fill(all_results)
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
        # Which test produced `statistic` decides whether it can be carried at
        # all. limma and the two t-tests store a signed, zero-centred t; the
        # supported `wilcoxon` method stores wilcox.test()'s W, which is
        # non-negative and centred at n1*n2/2. Handing W to a ranker that treats
        # its input as a signed score gives fgsea a list with no low end, so
        # decreased metabolites get small positive weights instead of strong
        # negative ranks -- and the duplicate collapse then prefers the largest
        # W rather than the strongest change. W is not converted into something
        # signed here: it is simply not carried, and the ranker falls back to
        # sign(log2fc) * -log10(p), which is correct for any method.
        #
        # `precomputed`, a missing method and anything unrecognised are treated
        # the same way, because none of them tells us what the column holds.
        signed_statistic <- tolower(as.character(de_data$method %||% "")) %in%
            c("limma", "t_test", "t_test_equal")

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
                    # Carried for compound GSEA, which ranks on the moderated
                    # statistic where there is one. Both metabolomics DE paths
                    # name it `statistic` -- the limma path renames topTable()'s
                    # `t` to that on the way out -- and `t` is accepted only for
                    # a table handed in from elsewhere. It has to survive here or
                    # not at all: sign(log2fc) and a p-value cannot reconstruct a
                    # moderated t, so a ranker downstream would have no way to
                    # prefer it and would silently always take the fallback.
                    # Gated on the method, per signed_statistic above.
                    statistic = if (!signed_statistic) NA_real_
                                else if ("statistic" %in% names(df)) df$statistic
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


#' Split protein-group IDs into their member accessions
#'
#' A protein group ("P1;P2;P3") names every protein its peptides could belong
#' to, leading protein first. None of those strings is itself a UniProt key, so
#' anything that looks accessions up has to take the group apart first.
#'
#' @param ids Character vector of protein-group IDs.
#' @return Data frame with one row per member: \code{feature_id} (the group ID
#'   as given), \code{accession} (the trimmed member) and \code{position}
#'   (1 for the leading protein). Missing and empty IDs are dropped.
#' @examples
#' protein_group_members(c("P1;P2", "P3"))
protein_group_members <- function(ids) {
    ids <- unique(as.character(ids))
    ids <- ids[!is.na(ids) & nzchar(trimws(ids))]
    parts <- strsplit(ids, ";", fixed = TRUE)
    n <- lengths(parts)
    out <- data.frame(
        feature_id = rep(ids, n),
        accession  = trimws(as.character(unlist(parts, use.names = FALSE))),
        position   = as.integer(unlist(lapply(n, seq_len), use.names = FALSE)),
        stringsAsFactors = FALSE
    )
    out[nzchar(out$accession), , drop = FALSE]
}


#' Canonical UniProt accession of an isoform accession
#'
#' UniProt writes isoforms as the canonical accession plus "-<n>" ("P12345-2").
#' Annotation packages key genes on the canonical accession, so an isoform with
#' no entry of its own can still be placed through it.
#'
#' @param accessions Character vector of accessions.
#' @return The accessions with a trailing "-<digits>" removed; others unchanged.
#' @examples
#' canonical_uniprot_accession(c("P12345-2", "P12345", "Q9ABC1"))
canonical_uniprot_accession <- function(accessions) {
    sub("-[0-9]+$", "", as.character(accessions))
}


#' First annotated gene symbol of each proteomics feature
#'
#' Reads the gene annotation that came in with the protein table (DIA-NN's
#' \code{Genes}, or the other columns \code{extract_protein_symbols()} knows)
#' and keeps its first symbol -- the convention harmonization already uses when
#' it matches proteins to transcripts. The annotation is joined to the features
#' by value, through the \code{row_data} column that holds the feature IDs, so a
#' reordered \code{row_data} still lines up. That column must be a 1:1 key --
#' as many rows as features, no duplicates, the same set of IDs -- or the
#' annotation cannot be trusted to line up and none is returned.
#'
#' @param prot_pre The proteomics entry of \code{harmonization_res$inputs}:
#'   a list carrying \code{row_data} and \code{expr_work}.
#' @return Character vector of first gene symbols (NA where a feature has none),
#'   named by feature ID; NULL when there is no usable, aligned annotation.
#' @examples
#' pre <- list(expr_work = matrix(0, 2, 1, dimnames = list(c("P1;P2", "P3"), "S1")),
#'             row_data = data.frame(Protein.Group = c("P3", "P1;P2"),
#'                                   Genes = c("GENE3", "GENEA;GENEB")))
#' protein_group_gene_symbols(pre)   # "P3" = "GENE3", "P1;P2" = "GENEA"
protein_group_gene_symbols <- function(prot_pre) {
    row_data <- prot_pre$row_data
    feature_ids <- rownames(prot_pre$expr_work)
    if (is.null(row_data) || nrow(row_data) == 0 || is.null(feature_ids)) return(NULL)

    # A 1:1 key: one row per feature and nothing else, in any order.
    holds_ids <- vapply(colnames(row_data), function(col) {
        v <- as.character(row_data[[col]])
        length(v) == length(feature_ids) && !anyDuplicated(v) && setequal(v, feature_ids)
    }, logical(1))
    if (!any(holds_ids)) {
        message("    No row_data column is a 1:1 key for the proteomics feature IDs; ",
                "protein groups are mapped through their accessions only")
        return(NULL)
    }

    symbols <- extract_protein_symbols(row_data, NULL)
    if (is.null(symbols)) return(NULL)
    stats::setNames(as.character(symbols),
                    as.character(row_data[[names(which(holds_ids))[1]]]))
}


#' Map proteomics features to Entrez IDs, one gene per feature
#'
#' A single-accession feature is looked up as a UniProt accession, and an
#' isoform accession with no entry of its own is tried again as its canonical
#' accession.
#'
#' A protein group ("P1;P2") is one measurement, so it gets one representative
#' gene, never one per member: expanding it would turn a single statistical
#' observation into several gene-level ones in ORA and GSEA. The representative
#' is the group's first annotated gene symbol when that symbol maps. When it is
#' missing or does not map, the group falls back to its accessions in order and
#' takes the first one that resolves in the annotation database -- an
#' annotation fallback, not a claim about which member the group "is".
#'
#' @param ids Character vector of feature IDs: UniProt accessions or protein
#'   groups of them.
#' @param lookup Function taking a character vector of unique accessions and
#'   returning Entrez IDs named by accession, NA where unmapped -- the shape
#'   \code{AnnotationDbi::mapIds()} returns.
#' @param symbols Optional character vector of first gene symbols named by
#'   feature ID, as \code{protein_group_gene_symbols()} returns. Read for
#'   protein groups only.
#' @param symbol_lookup Function like \code{lookup}, keyed by gene symbol.
#'   Required for \code{symbols} to be used.
#' @return Data frame with one row per mapped feature, ordered by
#'   \code{feature_id}: \code{feature_id}, \code{ENTREZID}, \code{source}
#'   (\code{"accession"}, \code{"canonical_accession"} or \code{"gene_symbol"})
#'   and \code{matched_key} (the accession or symbol that resolved). A feature
#'   nothing resolves for is left out.
#' @examples
#' acc <- function(keys) setNames(c(P1 = "101", P2 = "102", P3 = NA)[keys], keys)
#' map_protein_groups_to_entrez(c("P1;P2", "P3"), acc)   # P1;P2 -> 101 only
map_protein_groups_to_entrez <- function(ids, lookup, symbols = NULL,
                                         symbol_lookup = NULL) {
    empty <- data.frame(feature_id = character(0), ENTREZID = character(0),
                        source = character(0), matched_key = character(0),
                        stringsAsFactors = FALSE)
    members <- protein_group_members(ids)
    if (nrow(members) == 0) return(empty)
    n_members <- stats::ave(members$position, members$feature_id, FUN = length)
    group_ids <- unique(members$feature_id[n_members > 1])

    by_symbol <- empty
    if (length(group_ids) > 0 && !is.null(symbols) && is.function(symbol_lookup)) {
        sym <- trimws(unname(as.character(symbols[group_ids])))
        has_sym <- !is.na(sym) & nzchar(sym)
        if (any(has_sym)) {
            hits <- symbol_lookup(unique(sym[has_sym]))
            entrez <- unname(as.character(hits[sym[has_sym]]))
            ok <- !is.na(entrez) & nzchar(entrez)
            by_symbol <- data.frame(feature_id = group_ids[has_sym][ok],
                                    ENTREZID = entrez[ok],
                                    source = rep("gene_symbol", sum(ok)),
                                    matched_key = sym[has_sym][ok],
                                    stringsAsFactors = FALSE)
        }
    }

    rest <- members[!members$feature_id %in% by_symbol$feature_id, , drop = FALSE]
    by_accession <- empty
    if (nrow(rest) > 0) {
        rest$canonical <- canonical_uniprot_accession(rest$accession)
        hits <- lookup(unique(c(rest$accession, rest$canonical)))
        hit_of <- function(keys) unname(as.character(hits[keys]))
        direct <- hit_of(rest$accession)
        via_canonical <- hit_of(rest$canonical)
        direct_ok <- !is.na(direct) & nzchar(direct)
        canonical_ok <- !is.na(via_canonical) & nzchar(via_canonical)
        rest$ENTREZID <- ifelse(direct_ok, direct, via_canonical)
        rest$source <- ifelse(direct_ok, "accession", "canonical_accession")
        rest$matched_key <- ifelse(direct_ok, rest$accession, rest$canonical)
        rest <- rest[direct_ok | canonical_ok, , drop = FALSE]
        # First member that resolves, in the group's own order.
        rest <- rest[order(rest$feature_id, rest$position), , drop = FALSE]
        rest <- rest[!duplicated(rest$feature_id), , drop = FALSE]
        by_accession <- rest[, names(empty), drop = FALSE]
    }

    out <- rbind(by_symbol, by_accession)
    out <- out[order(out$feature_id), , drop = FALSE]
    rownames(out) <- NULL
    out
}


#' Map feature IDs to Entrez gene IDs
#'
#' Transcriptomics IDs are looked up as ENSEMBL, then WORMBASE, keys.
#' Proteomics features go through \code{map_protein_groups_to_entrez()}; when
#' no accession resolves at all, the whole layer falls back to the WormBase /
#' gene_id column of \code{row_data} (C. elegans and similar). When that
#' fallback is unavailable or maps nothing, the groups already resolved by
#' gene symbol are returned rather than nothing.
#'
#' @param de_tables Named list of DE data frames, each with a \code{feature_id}
#'   column; the IDs of every table are mapped together.
#' @param omics_type One of \code{"transcriptomics"}, \code{"proteomics"},
#'   \code{"metabolomics"}.
#' @param harmonization_res Harmonization result; its
#'   \code{inputs$proteomics} supplies the gene annotation and the fallback
#'   gene IDs for proteomics.
#' @param org_db OrgDb annotation object.
#' @return Data frame with \code{feature_id} and \code{ENTREZID}, \strong{at
#'   most one row per \code{feature_id}}; NULL when nothing can be mapped (and
#'   always for metabolomics). Proteomics: one representative gene per
#'   feature, as \code{map_protein_groups_to_entrez()} describes; unmapped
#'   features are omitted. Transcriptomics: one row per input ID, with NA
#'   \code{ENTREZID} where unmapped. Several features may share an
#'   \code{ENTREZID}; collapsing by gene is the caller's job.
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
        # Try direct UniProt -> ENTREZID mapping first (works for most organisms).
        # A group string such as "P1;P2" is never a UniProt key itself, so groups
        # are resolved to one representative gene each. mapIds() errors when
        # none of the keys is valid; that reads as "nothing mapped" here, as it
        # did when the whole lookup sat in one tryCatch.
        entrez_via <- function(keytype) function(keys) tryCatch(
            AnnotationDbi::mapIds(org_db, keys = keys, keytype = keytype,
                                  column = "ENTREZID", multiVals = "first"),
            error = function(e) stats::setNames(rep(NA_character_, length(keys)), keys))
        entrez_df <- tryCatch({
            df <- map_protein_groups_to_entrez(
                all_ids,
                lookup = entrez_via("UNIPROT"),
                symbols = protein_group_gene_symbols(harmonization_res$inputs$proteomics),
                symbol_lookup = entrez_via("SYMBOL")
            )
            n_src <- table(factor(df$source, levels = c("accession", "canonical_accession",
                                                        "gene_symbol")))
            message("    Mapped ", nrow(df), "/", length(all_ids),
                    " proteomics features to ENTREZID (", n_src[["accession"]],
                    " by accession, ", n_src[["canonical_accession"]],
                    " through the canonical accession of an isoform, ",
                    n_src[["gene_symbol"]], " protein groups by their first gene symbol)")
            df
        }, error = function(e) NULL)

        # The WormBase fallback below stays a whole-layer decision, taken when no
        # accession resolved -- the same condition as before protein groups were
        # handled. Groups resolved by gene symbol alone do not count, or a layer
        # that needs the fallback would skip it on the strength of a few groups.
        if (!is.null(entrez_df) &&
            any(entrez_df$source %in% c("accession", "canonical_accession"))) {
            return(entrez_df[, c("feature_id", "ENTREZID"), drop = FALSE])
        }
        # Groups already resolved by gene symbol are what this layer returns if
        # the fallback below is unavailable or maps nothing -- not NULL.
        symbol_only <- if (!is.null(entrez_df) && nrow(entrez_df) > 0) {
            entrez_df[, c("feature_id", "ENTREZID"), drop = FALSE]
        }

        # Fallback: try via row_data WormBase/gene_id columns (C. elegans etc.)
        prot_pre <- harmonization_res$inputs$proteomics
        if (is.null(prot_pre) || is.null(prot_pre$row_data)) {
            message("    No proteomics row_data for ID mapping")
            return(symbol_only)
        }

        row_data <- prot_pre$row_data
        wbgene_col <- intersect(c("Wormbase_id", "wormbase_id", "gene_id"), colnames(row_data))
        if (length(wbgene_col) == 0) {
            message("    No WormBase/gene_id column in proteomics row_data")
            return(symbol_only)
        }

        prot_ids <- rownames(prot_pre$expr_work)
        wb_ids <- row_data[[wbgene_col[1]]]
        prot_to_wb <- data.frame(
            feature_id = prot_ids,
            WBGene = wb_ids,
            stringsAsFactors = FALSE
        )
        prot_to_wb <- prot_to_wb[!is.na(prot_to_wb$WBGene) & nzchar(prot_to_wb$WBGene), ]

        if (nrow(prot_to_wb) == 0) return(symbol_only)

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

        if (is.null(mapped) || !any(!is.na(mapped$ENTREZID))) {
            return(symbol_only %||% mapped)
        }
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


# =============================================================================
# Compound GSEA (metabolomics)
# =============================================================================
#
# Deliberately outside the ORA path. Compound GSEA is never added to
# `pathway_tables`, because merge_pathway_pvalues() aggregates a layer's rows
# with FUN = min and applies no method filter: a GSEA row sitting beside an ORA
# row would make the metabolomics p-value entering Stouffer the smaller of the
# two, silently changing the cross-omics meta-analysis. It therefore has its own
# orchestrator, its own CSV and its own report section, and run_compound_ora()
# and run_kegg_enrichment_for_omics() are untouched.

#' The cached KEGG compound-pathway table, and only the cache
#'
#' \code{get_kegg_compound_pathways()} falls back to the KEGG REST API when its
#' cache is cold. Compound GSEA runs after compound ORA has already populated
#' that cache, and must never be the call that goes out to the network -- so it
#' reads the cache file directly and gives up when it is not there. Structural
#' rather than conditional: there is no branch here that could fetch.
#'
#' @param cache_dir Directory \code{get_kegg_compound_pathways()} caches into,
#'   or NULL.
#' @return The cached data frame of compound-pathway associations, or NULL when
#'   no usable cache exists.
.cached_compound_pathways <- function(cache_dir) {
    if (is.null(cache_dir) || !nzchar(cache_dir)) return(NULL)
    # Same filename get_kegg_compound_pathways() writes.
    f <- file.path(cache_dir, "kegg_compound_pathways.rds")
    if (!file.exists(f)) return(NULL)

    cpd <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.data.frame(cpd) || nrow(cpd) == 0) return(NULL)
    if (!all(c("compound", "pathway") %in% names(cpd))) return(NULL)
    cpd
}


#' The cached KEGG BRITE classification, and only the cache
#'
#' The companion to \code{.cached_compound_pathways()}, for the other KEGG
#' resource this path can reach. \code{keep_kegg_pathways()} defaults its
#' \code{classification} argument to \code{kegg_pathway_categories()}, which
#' downloads br08901 whenever its own cache is cold or a week old -- so leaving
#' that default in place would have made the compound GSEA path fetch after all,
#' indirectly, the moment a project configured \code{exclude_pathway_classes}.
#'
#' Reading the cache and stopping there keeps the guarantee literal.
#' \code{keep_kegg_pathways()} already fails open on a NULL classification, and
#' says so, which is the right outcome: a class exclusion that cannot be
#' resolved keeps every pathway rather than silently dropping some.
#'
#' @param cache_dir Directory \code{kegg_pathway_categories()} caches into;
#'   NULL falls back to the same default that function uses.
#' @return The cached classification table, or NULL when there is no usable one.
.cached_pathway_categories <- function(cache_dir = NULL) {
    if (is.null(cache_dir) || !nzchar(cache_dir)) {
        cache_dir <- file.path(tempdir(), "kegg_cache")
    }
    # Same filename kegg_pathway_categories() writes.
    f <- file.path(cache_dir, "kegg_pathway_categories.rds")
    if (!file.exists(f)) return(NULL)

    cls <- tryCatch(readRDS(f), error = function(e) NULL)
    # Validated with the producer's own check, not a second opinion about shape.
    if (!.is_kegg_category_table(cls)) return(NULL)
    cls
}


#' Rank mapped metabolites for compound GSEA
#'
#' Two statistics, in order, and no third: the moderated statistic where the DE
#' table carries usable values, and \code{sign(log2fc) * -log10(pvalue)}
#' otherwise.
#'
#' The choice is made on the VALUES, not on the column. The standardised
#' metabolomics DE table always carries a `statistic` column and fills it with
#' NA where the source had none, so testing for the column would pick an all-NA
#' vector, and fgsea would then report no gene-set overlap -- sending the reader
#' after an ID-mapping bug that is not there.
#'
#' The p-value gets the same 1e-300 floor used elsewhere in this file, so a
#' p-value that underflowed to zero ranks very high rather than infinite. A zero
#' log2 fold change ranks at zero and is kept: no change is evidence of no
#' change, not missing evidence.
#'
#' Duplicates are collapsed here rather than upstream. Several metabolites can
#' annotate to one KEGG compound, and fgsea would otherwise score that compound
#' more than once. The strongest absolute rank wins. Two rows of one compound
#' can still tie on magnitude while differing in sign -- +2 against -2 -- and
#' the compound id cannot separate them, so the more positive rank is taken.
#' That last rule is arbitrary and deliberately so: it is there to make the
#' outcome independent of the DE table's row order, not to express a preference
#' about direction. Without it the first row to arrive won.
#'
#' This is GSEA's own rule -- the ORA path keeps its own upstream
#' de-duplication, which this does not touch.
#'
#' @param de_mapped DE table merged with KEGG compound ids: `KEGG_ID`,
#'   `log2fc`, `pvalue`, and optionally `statistic`.
#' @return Named numeric vector of ranks, names being KEGG compound ids, sorted
#'   decreasing and unique; \code{numeric(0)} when nothing is rankable.
#' @examples
#' de <- data.frame(KEGG_ID = c("C00031", "C00022"), statistic = c(3.1, -2.4),
#'                  log2fc = c(1, -1), pvalue = c(0.01, 0.02))
#' rank_compounds_for_gsea(de)
rank_compounds_for_gsea <- function(de_mapped) {
    if (!is.data.frame(de_mapped) || nrow(de_mapped) == 0) return(numeric(0))
    if (!"KEGG_ID" %in% names(de_mapped)) return(numeric(0))

    stat <- if ("statistic" %in% names(de_mapped)) {
        suppressWarnings(as.numeric(de_mapped$statistic))
    } else {
        NULL
    }

    if (!is.null(stat) && any(is.finite(stat))) {
        ranks <- stat
    } else {
        if (!all(c("log2fc", "pvalue") %in% names(de_mapped))) return(numeric(0))
        lfc  <- suppressWarnings(as.numeric(de_mapped$log2fc))
        pval <- suppressWarnings(as.numeric(de_mapped$pvalue))
        ranks <- sign(lfc) * -log10(pval + 1e-300)
    }

    names(ranks) <- as.character(de_mapped$KEGG_ID)
    keep <- is.finite(ranks) & !is.na(names(ranks)) & nzchar(names(ranks))
    ranks <- ranks[keep]
    if (length(ranks) == 0) return(numeric(0))

    # Three keys, and the third is load-bearing: magnitude, then the compound
    # id, then the signed rank. Without the last one, two rows of one compound
    # at +x and -x tie on both earlier keys and order() falls back to the
    # arrival index, which is precisely the row-order dependence this collapse
    # exists to remove.
    ranks <- ranks[order(-abs(ranks), names(ranks), -ranks)]
    ranks <- ranks[!duplicated(names(ranks))]

    # Ordered explicitly rather than through sort(), so that two compounds with
    # the same rank land in one fixed order instead of relying on sort stability.
    ranks[order(-ranks, names(ranks))]
}


#' Run GSEA over KEGG compound pathways
#'
#' The companion to \code{run_compound_ora()}, not a replacement: ORA asks which
#' pathways are crowded with significantly changed metabolites, this asks
#' whether a pathway's metabolites sit systematically high or low in the ranked
#' list, and reports a direction with it.
#'
#' \code{fgsea::fgseaMultilevel()} is called by name rather than through
#' \code{fgsea::fgsea()}, which dispatches to it only while no argument says
#' otherwise. Which statistical procedure runs should not depend on dispatch
#' behaviour that a future version could change underneath the pipeline. It is
#' stochastic either way, so the call is seeded -- the same
#' \code{withr::with_seed()} treatment \code{run_pathway_analysis()} gives its
#' own fgsea call, for the same reason: unseeded, pathways near the cutoff cross
#' it between otherwise identical runs.
#'
#' Nothing is filtered on significance. fgsea's `padj` is a BH family over the
#' pathways it scored, and the complete scored table is what reaches the export;
#' choosing what to show is the report's business. Class exclusion is applied
#' after scoring, so a project's reporting choice cannot move that family --
#' the rule \code{run_compound_ora()} follows.
#'
#' @param de_mapped DE table merged with KEGG compound ids; see
#'   \code{rank_compounds_for_gsea()} for the columns used.
#' @param cache_dir Directory holding the cached compound-pathway table.
#' @param min_gs,max_gs Pathway size bounds, counted on measured compounds.
#' @param seed Seed for the stochastic scoring.
#' @param exclude_classes KEGG BRITE classes to drop after scoring, or NULL.
#' @param cpd_pathways Compound-pathway associations. Defaults to the cache and
#'   never to a download; passing it explicitly lets a caller score many
#'   contrasts from one read, and lets a test supply its own.
#' @return Data frame of every pathway fgsea scored, or NULL.
run_compound_gsea <- function(de_mapped, cache_dir, min_gs, max_gs, seed = 1L,
                               exclude_classes = NULL,
                               cpd_pathways = .cached_compound_pathways(cache_dir)) {

    if (is.null(cpd_pathways)) {
        message("    No cached KEGG compound-pathway associations; ",
                "skipping compound GSEA")
        return(NULL)
    }

    ranks <- rank_compounds_for_gsea(de_mapped)
    if (length(ranks) < 3) {
        message("    Too few rankable compounds for GSEA (", length(ranks), ")")
        return(NULL)
    }

    pathway_names <- stats::setNames(
        if ("name" %in% names(cpd_pathways)) cpd_pathways$name else cpd_pathways$pathway,
        cpd_pathways$pathway)
    pathway_names <- pathway_names[!duplicated(names(pathway_names))]

    # Sets are restricted to what was measured, so a size is a count of
    # compounds that could actually contribute -- the same basis
    # run_compound_ora() sizes on.
    pathway_sets <- split(cpd_pathways$compound, cpd_pathways$pathway)
    pathway_sets <- lapply(pathway_sets,
                           function(cpds) unique(intersect(cpds, names(ranks))))

    # The same compound-specific floor run_compound_ora() applies, and for the
    # same reason: min_set_size is configured on a gene-set scale, where 10 is
    # modest, but a KEGG compound pathway rarely has ten MEASURED members in one
    # experiment. Clamping up to the configured value rather than down would
    # leave compound GSEA testing almost nothing under the shipped default --
    # working, reporting no error, and finding nothing. The configured maximum
    # is honoured as given.
    use_min_gs <- max(2, min(min_gs, 3))
    set_sizes <- lengths(pathway_sets)
    testable <- set_sizes >= use_min_gs & set_sizes <= max_gs

    if (!any(testable)) {
        message("    No compound pathway carries between ", use_min_gs, " and ",
                max_gs, " measured compounds; skipping compound GSEA")
        return(NULL)
    }
    pathway_sets <- pathway_sets[testable]

    # Checked here rather than on entry, so that a run without the package still
    # reports which of the inputs was the problem instead of blaming fgsea for
    # an empty ranking or a mapping that matched nothing.
    if (!requireNamespace("fgsea", quietly = TRUE)) {
        message("    fgsea not available; skipping compound GSEA")
        return(NULL)
    }

    fgsea_res <- tryCatch(
        withr::with_seed(seed, fgsea::fgseaMultilevel(
            pathways = pathway_sets,
            stats = ranks,
            minSize = use_min_gs,
            maxSize = max_gs,
            nPermSimple = 10000
        )),
        error = function(e) {
            warning("Compound GSEA failed: ", e$message)
            NULL
        }
    )

    if (is.null(fgsea_res)) return(NULL)
    df <- as.data.frame(fgsea_res)
    if (nrow(df) == 0) {
        message("    Compound GSEA scored no pathways")
        return(NULL)
    }

    pw_ids <- as.character(df$pathway)
    labels <- unname(pathway_names[pw_ids])
    unnamed <- is.na(labels) | !nzchar(labels)
    labels[unnamed] <- pw_ids[unnamed]

    # `ID` carries the accession and `pathway` the readable name, which is what
    # pathway_join_key() and pathway_display_label() each read first.
    out <- data.frame(
        pathway     = labels,
        ID          = pw_ids,
        pvalue      = df$pval,
        padj        = df$padj,
        NES         = df$NES,
        ES          = df$ES,
        setSize     = df$size,
        leadingEdge = if ("leadingEdge" %in% names(df)) {
            vapply(df$leadingEdge, paste, character(1), collapse = ",")
        } else {
            NA_character_
        },
        database    = "KEGG",
        # Said outright, so no consumer has to infer the method from which
        # columns happen to be present.
        method      = "fgsea",
        stringsAsFactors = FALSE
    )

    if (length(unlist(exclude_classes)) > 0) {
        # classification passed explicitly: the default would resolve through
        # kegg_pathway_categories(), which downloads br08901 on a cold cache.
        # See .cached_pathway_categories(). A NULL here fails open, keeping
        # every pathway and saying so, which is the safe direction for an
        # exclusion that cannot be resolved.
        out <- out[keep_kegg_pathways(out$ID, exclude = exclude_classes,
                                      cache_dir = cache_dir,
                                      label = "compound GSEA pathways",
                                      classification =
                                          .cached_pathway_categories(cache_dir)), ,
                   drop = FALSE]
        if (nrow(out) == 0) {
            message("    No compound GSEA pathways left after KEGG class exclusion")
            return(NULL)
        }
    }

    out <- out[order(out$pvalue, out$ID), , drop = FALSE]
    rownames(out) <- NULL
    message("    Compound GSEA scored ", nrow(out), " pathways")
    out
}


#' Run compound GSEA for every metabolomics contrast
#'
#' The one caller of \code{run_compound_gsea()} in the pipeline, kept separate
#' from \code{run_kegg_enrichment_for_omics()} on purpose: that function's
#' return becomes `pathway_tables`, and everything in `pathway_tables` reaches
#' the cross-omics meta-analysis. Scoring here instead is what keeps the ORA and
#' Stouffer paths unable to see a GSEA row at all.
#'
#' The compound-pathway table is read once, from the cache the ORA path filled,
#' and shared by every contrast. A cold cache means no GSEA rather than a
#' download.
#'
#' @param de_results Named list of per-omics DE results.
#' @param harmonization_res Harmonization result, for the ID mapping.
#' @param config Full config object.
#' @param out_dir The cross-enrichment output directory, whose `metabolomics`
#'   subdirectory holds the cache \code{run_compound_ora()} writes.
#' @return One data frame for all contrasts, carrying `contrast` and `omics`
#'   columns, or NULL when nothing could be scored.
run_compound_gsea_for_contrasts <- function(de_results, harmonization_res,
                                             config, out_dir) {

    de_data <- de_results[["metabolomics"]]
    if (is.null(de_data)) return(NULL)

    de_tables <- extract_de_tables(de_data, "metabolomics", harmonization_res)
    if (length(de_tables) == 0) {
        message("  No metabolomics DE tables for compound GSEA")
        return(NULL)
    }

    cache_dir <- file.path(out_dir, "metabolomics")
    cpd_pathways <- .cached_compound_pathways(cache_dir)
    if (is.null(cpd_pathways)) {
        # Named precisely, because the cause is not obvious from the symptom:
        # this table is written by the compound ORA step, and a run where that
        # step never executed -- metabolomics enrichment supplied pre-computed,
        # for instance -- leaves no cache for GSEA to read. Deliberately not
        # fetched here; a display-side analysis is not the right place to start
        # a KEGG download.
        message("  Compound GSEA skipped: no cached KEGG compound-pathway ",
                "table at ", cache_dir, ". It is written by the compound ORA ",
                "step, which has not run for this output directory.")
        return(NULL)
    }

    id_map <- tryCatch(
        map_metabolite_ids_to_kegg(de_tables, harmonization_res),
        error = function(e) {
            message("  Metabolite ID mapping failed: ", e$message)
            NULL
        }
    )
    if (is.null(id_map) || nrow(id_map) == 0) {
        message("  Could not map metabolite IDs to KEGG compounds for GSEA")
        return(NULL)
    }
    id_map$KEGG_ID <- id_map$KEGG_CPD

    enrich_cfg <- config$modes$multiomics$enrichment
    min_gs <- enrich_cfg$min_set_size %||% 10
    max_gs <- enrich_cfg$max_set_size %||% 500
    seed   <- config$params$seed %||% 1L
    excl   <- .excluded_pathway_classes(config)

    results <- list()
    for (cname in names(de_tables)) {
        de_mapped <- merge(de_tables[[cname]], id_map, by = "feature_id")
        de_mapped <- de_mapped[!is.na(de_mapped$KEGG_ID), , drop = FALSE]

        # Not de-duplicated here: rank_compounds_for_gsea() collapses compounds
        # on the strongest rank, which is a better rule than first-row-wins and
        # is the one GSEA should use.
        if (nrow(de_mapped) < 3) {
            message("    ", cname, ": too few mapped metabolites for GSEA (",
                    nrow(de_mapped), ")")
            next
        }

        res <- run_compound_gsea(de_mapped, cache_dir = cache_dir,
                                 min_gs = min_gs, max_gs = max_gs, seed = seed,
                                 exclude_classes = excl,
                                 cpd_pathways = cpd_pathways)
        if (is.null(res) || nrow(res) == 0) next

        res$contrast <- cname
        res$omics <- "metabolomics"
        results[[cname]] <- res
    }

    if (length(results) == 0) return(NULL)

    combined <- .rbind_fill(results)
    if (is.null(combined) || nrow(combined) == 0) return(NULL)
    rownames(combined) <- NULL
    combined
}


#' Write the compound GSEA export, and only for the run that produced it
#'
#' The report includes this file on \code{file.exists()} alone, and a run can
#' legitimately produce no GSEA at all -- no metabolomics DE, no compound
#' mapping, a cold cache, nothing scorable. Writing conditionally but never
#' clearing would leave the previous run's result in place for the report to
#' present as this run's.
#'
#' So the file is removed before the decision, not instead of it: after this
#' returns, the export exists if and only if this invocation produced rows. The
#' same invariant the cross-omics ORA figure was given.
#'
#' @param compound_gsea Result of \code{run_compound_gsea_for_contrasts()}, or
#'   NULL.
#' @param out_dir Cross-enrichment output directory.
#' @return Invisibly, the path when one was written, otherwise NULL.
write_compound_gsea_export <- function(compound_gsea, out_dir) {
    path <- file.path(out_dir, "metabolomics_compound_gsea.csv")

    if (file.exists(path)) unlink(path)

    if (is.null(compound_gsea) || !is.data.frame(compound_gsea) ||
        nrow(compound_gsea) == 0) {
        return(invisible(NULL))
    }

    write.csv(compound_gsea, path, row.names = FALSE)
    invisible(path)
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
                # Said outright, as the other producers do: the cross-omics
                # merge chooses rank-based rows by this column, and an unlabelled
                # GSEA row would be treated as a method it cannot identify.
                method = "gsea",
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
                # Stamped here, where the frame is built, and not after the
                # tryCatch() below: both filtered subsets inherit it, and both
                # leave this function from inside the tryCatch expression. See
                # the note above that call for why nothing after it is reached.
                method = "ora",
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

    # Nothing here sees a successful clusterProfiler result. `return()` inside a
    # tryCatch() expression leaves the ENCLOSING function, because the expression
    # is a promise evaluated in this frame -- so both hit branches above exit
    # run_ora_kegg() outright and ora_res is only ever NULL. That is why the
    # method stamp lives in the frame construction and not on ora_res: put here,
    # it would be dead code that reads as though it were doing the job.
    # The guard stays, so that a future edit which stops returning early still
    # behaves, and run_gsea_kegg() has the same shape for the same reason.
    if (!is.null(ora_res)) return(ora_res)

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
#' @param rank_tables Optional named list of rank-based result tables per omics
#'   that are not part of \code{enrichment_results} -- compound GSEA for
#'   metabolomics, which is kept out of the per-omics tables on purpose. Their
#'   rows join that layer's evidence for the meta-analysis only, where
#'   \code{merge_pathway_pvalues()} picks one method per layer; the per-omics
#'   bar plots, CSVs and the ORA figure are drawn from \code{enrichment_results}
#'   alone. A layer present only here still counts towards the two needed.
#' @return List with: common_pathways, union_pathways, meta_analysis,
#'   pathway_tables, layer_methods (the test each layer contributed to the
#'   meta-analysis, named by layer), lookup_files (see
#'   \code{write_cross_omics_lookups()}) and plots
analyze_cross_omics_enrichment <- function(enrichment_results, config, out_dir = NULL,
                                           rank_tables = NULL) {

    # Cleared before anything can return, not beside the code that writes the
    # figures. These two are now written one per gene-set collection, and which
    # collections exist depends on the run -- so a previous run's files, under
    # either the old combined names or a collection this run does not produce,
    # would otherwise survive every path that declines to draw: too few layers,
    # no pathway keys, nothing left after filtering. The report finds these by
    # glob, so a survivor is shown as this run's evidence.
    #
    # Scoped to these two figure families and to this out_dir, which covers the
    # per-contrast directories too because this function runs again for each.
    if (!is.null(out_dir) && dir.exists(out_dir)) {
        .clear_collection_heatmaps(out_dir)
        .clear_cross_lookup_outputs(out_dir)
    }

    rank_tables <- .usable_rank_tables(rank_tables)

    if (length(union(names(enrichment_results), names(rank_tables))) < 2) {
        message("Cross-omics enrichment requires >= 2 omics layers with enrichment results")
        return(NULL)
    }

    message("Analyzing cross-omics pathway enrichment...")

    omics <- union(names(enrichment_results), names(rank_tables))

    # Extract pathway-level results from each omics
    pathway_tables <- list()
    for (om in names(enrichment_results)) {
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

    # What the meta-analysis reads: each layer's own table plus any rank-based
    # rows supplied beside it. pathway_tables stays exactly what the layers
    # handed over, because it drives the per-layer figures and CSVs.
    merge_tables <- pathway_tables
    for (om in names(rank_tables)) {
        merge_tables[[om]] <- .rbind_fill(list(merge_tables[[om]], rank_tables[[om]]))
    }
    omics <- intersect(omics, names(merge_tables))

    if (length(merge_tables) < 2) {
        warning("Insufficient pathway tables for cross-omics enrichment")
        return(NULL)
    }

    # Join the layers on a stable identity rather than on whichever column each
    # of them happens to use. Without this a gene layer keyed on hsa00010 and a
    # compound layer keyed on map00010 are two different pathways, and a custom
    # collection keyed on map00010 never meets a gene layer keyed on the readable
    # KEGG description at all.
    kegg_org <- resolve_kegg_org_code(config$global$organism)

    # KEGG's reference maps are pan-species, so an organism with no KEGG code of
    # its own can score well on maps of organs it does not have. Excluding those
    # classes is a reporting decision, applied to finished results: the p-values
    # and the adjustment behind them are untouched, and nothing is excluded
    # unless the config asks. The per-omics tables get the same treatment --
    # they drive the per-layer barplots and CSVs, and filtering only the merged
    # selection would leave the excluded classes visible one section away.
    #
    # Applied before the candidates are found, not after: which method a layer
    # contributes is chosen from the rows it still has, so filtering later could
    # change that choice under candidates already taken from the other method.
    excl <- .excluded_pathway_classes(config)
    if (length(excl) > 0) {
        drop_excluded <- function(df) {
            col <- if ("ID" %in% names(df)) "ID"
                   else if ("pathway" %in% names(df)) "pathway" else NULL
            if (is.null(col) || nrow(df) == 0) return(df)
            df[keep_kegg_pathways(df[[col]], exclude = excl, kegg_org = kegg_org,
                                  label = "per-omics pathways"), , drop = FALSE]
        }
        pathway_tables <- lapply(pathway_tables, drop_excluded)
        merge_tables <- lapply(merge_tables, drop_excluded)
    }

    # Candidates come from the rows that can contribute -- the same rows
    # merge_pathway_pvalues() aggregates, through the same helper -- so a
    # pathway that only an unselected method's rows support is not a candidate,
    # and cannot reach the meta-analysis with no layer behind it.
    contrib <- lapply(merge_tables[omics], .layer_contribution, kegg_org = kegg_org)
    all_pathways <- lapply(contrib, function(cb) {
        if (!is.null(cb$problem)) return(character(0))
        unique(cb$keys[cb$keep])
    })

    # A layer contributes one method; say when that leaves a contrast or a
    # collection of it out of the run-level meta-analysis entirely.
    for (om in omics) {
        lost <- .method_selection_losses(merge_tables[[om]], contrib[[om]])
        if (lost$n_rows == 0) next
        msg <- sprintf("  %s: the meta-analysis uses its %s rows only; %d row(s) of other methods are left out",
                       om, contrib[[om]]$method, lost$n_rows)
        if (length(lost$contrasts) > 0) {
            msg <- paste0(msg, ", including every row for contrast(s) ",
                          paste(lost$contrasts, collapse = ", "))
        }
        if (length(lost$collections) > 0) {
            msg <- paste0(msg, if (length(lost$contrasts) > 0) " and" else ", including",
                          " every row for collection(s) ",
                          paste(lost$collections, collapse = ", "))
        }
        message(msg)
    }

    union_pathways <- unique(unlist(all_pathways))
    common_pathways <- Reduce(intersect, all_pathways)

    if (length(union_pathways) == 0) {
        message("No pathways found across omics layers")
        return(NULL)
    }

    message(sprintf("  Found %d total pathways (%d in common) across %d omics layers",
                    length(union_pathways), length(common_pathways), length(omics)))

    # The candidate universe is the union of what the layers contribute, always.
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
    if (length(excl) > 0) {
        use_pathways <- use_pathways[
            keep_kegg_pathways(use_pathways, exclude = excl, kegg_org = kegg_org,
                               label = "cross-omics pathways")]
    }

    # Merge pathway p-values for meta-analysis, one method per layer -- see
    # merge_pathway_pvalues(). The method each layer contributed is recorded
    # beside its p-values, so the table says what was combined.
    merged_pathways <- merge_pathway_pvalues(merge_tables, use_pathways, omics,
                                              kegg_org = kegg_org)
    layer_methods <- vapply(omics, function(om) {
        col <- paste0("method_", om)
        vals <- if (col %in% names(merged_pathways)) merged_pathways[[col]] else NA
        vals <- unique(vals[!is.na(vals)])
        if (length(vals) == 0) NA_character_ else vals[1]
    }, character(1))
    message("  Meta-analysis p-values per layer: ",
            paste(sprintf("%s = %s", omics, ifelse(is.na(layer_methods), "none",
                                                   layer_methods)),
                  collapse = "; "))

    # Combine p-values using Stouffer's method. Every candidate pathway gets a
    # row: the ones a single layer enriched are kept and carry n_omics = 1,
    # rather than being dropped here. Requiring two layers is a question for
    # whoever reads or ranks the table, and n_omics is what lets them ask it.
    meta_results <- stouffer_combined_pvalues(merged_pathways)

    # The key joined the layers; the label is what a reader sees. Both are kept,
    # so the table can be traced back to the accession that produced a row.
    meta_results <- attach_pathway_display_names(meta_results, merge_tables,
                                                  kegg_org = kegg_org)

    # Sort by combined p-value
    meta_results <- meta_results[order(meta_results$combined_pval), ]

    # Generate plots
    plots <- list()
    if (!is.null(out_dir) && nrow(meta_results) > 0) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

        # 1. Cross-omics heatmaps, one per gene-set collection.
        #
        # The collections cover neither the same layers nor the same identifier
        # space -- GO is annotated for the gene layers only, while the KEGG map
        # space is the one every layer can share -- and they differ by an order
        # of magnitude in size. A single combined figure therefore gave its rows
        # to whichever collection was largest, and read as an analysis of that
        # collection alone.
        #
        # Each figure keeps only the layers that scored something in its
        # collection, so an absent column means "this layer carries no
        # annotation here", not "this layer found nothing".
        # Keyed on norm_id, passed explicitly. The classifier's default would
        # take pathway_join_key(), which prefers `pathway` here -- and by this
        # point attach_pathway_display_names() has made that a readable label,
        # so every row would classify as "Other". norm_id is the accession this
        # table is actually keyed on.
        collections <- classify_pathway_collection(
            meta_results, kegg_org, keys = as.character(meta_results$norm_id))
        for (cl in sort(unique(collections[!is.na(collections)]))) {
            sub <- meta_results[collections == cl, , drop = FALSE]
            if (nrow(sub) == 0) next

            keep_omics <- .layers_with_values(sub, omics)
            if (length(keep_omics) == 0) next

            slug <- .collection_slug(cl)
            key <- paste0("pathway_heatmap_", slug)
            plots[[key]] <- file.path(
                out_dir, paste0("cross_omics_pathway_heatmap_", slug, ".png"))
            png(plots[[key]], width = 1200, height = 900, res = 120)
            tryCatch({
                plot_cross_omics_pathway_heatmap(sub, keep_omics, collection = cl)
            }, error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("Heatmap failed:", e$message), cex = 1.2)
            })
            dev.off()
        }

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

        # 3. ORA evidence per layer, on each layer's own adjusted p-value, and
        # likewise one figure per collection.
        #
        # The matrix is rebuilt from the collection's own pathways rather than
        # subset afterwards: this figure ranks on its own ORA evidence, so a
        # matrix built over every pathway would rank against rows the figure
        # does not show. Built before each device is opened so that a
        # collection with no adjusted ORA p-value at all leaves no figure
        # behind for the report to find.
        #
        # Its candidates come from the per-layer tables it draws from, not from
        # the meta-analysis: the meta-analysis keeps one method per layer, so a
        # layer contributing its rank-based rows there still has ORA evidence --
        # GO terms only ORA scored, say -- that belongs in this figure.
        ora_pathways <- .ora_figure_candidates(pathway_tables, kegg_org, exclude = excl)
        ora_collections <- classify_pathway_collection(NULL, kegg_org, keys = ora_pathways)
        any_ora <- FALSE
        for (cl in sort(unique(ora_collections))) {
            cl_pathways <- ora_pathways[ora_collections == cl]
            if (length(cl_pathways) == 0) next

            ora_padj <- build_ora_adjusted_p_matrix(pathway_tables, cl_pathways,
                                                    omics, kegg_org = kegg_org)
            # The builder returns a column per requested layer, all-NA for a
            # layer that contributed nothing to this collection. Dropping those
            # keeps the figure consistent with its legend: a layer is absent,
            # not shown empty.
            informative_cols <- colSums(!is.na(ora_padj)) > 0
            ora_padj <- ora_padj[, informative_cols, drop = FALSE]
            if (ncol(ora_padj) == 0 || !any(!is.na(ora_padj))) next

            any_ora <- TRUE
            slug <- .collection_slug(cl)
            key <- paste0("ora_heatmap_", slug)
            plots[[key]] <- file.path(
                out_dir, paste0("cross_omics_ora_heatmap_", slug, ".png"))
            png(plots[[key]], width = 1200, height = 900, res = 120)
            tryCatch({
                plot_cross_omics_ora_heatmap(ora_padj, pathway_tables,
                                             kegg_org = kegg_org, collection = cl)
            }, error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("ORA heatmap failed:", e$message), cex = 1.2)
            })
            dev.off()
        }
        if (!any_ora) {
            message("  No adjusted ORA p-values across layers; ",
                    "skipping the cross-omics ORA heatmaps")
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

    # 5. Each layer's top pathways, looked up in every other layer. Read from
    # merge_tables, so the rows each layer contributes are the ones the
    # meta-analysis used -- rank-based where the layer has them.
    lookup_files <- list()
    lookup_cfg <- .cross_lookup_config(config)
    if (!is.null(out_dir) && lookup_cfg$enabled) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        lookup_files <- tryCatch(
            write_cross_omics_lookups(merge_tables, omics, out_dir,
                                      top_n = lookup_cfg$top_n, kegg_org = kegg_org),
            error = function(e) {
                message("  Cross-omics lookup failed: ", e$message)
                list()
            })
        for (nm in names(lookup_files)) {
            pngs <- lookup_files[[nm]][names(lookup_files[[nm]]) != "tsv"]
            for (rk in names(pngs)) plots[[paste0("lookup_", nm, "_", rk)]] <- pngs[[rk]]
        }
    }

    list(
        common_pathways = common_pathways,
        union_pathways = union_pathways,
        meta_analysis = meta_results,
        pathway_tables = pathway_tables,
        layer_methods = layer_methods,
        lookup_files = lookup_files,
        plots = plots
    )
}


#' Keep the rank-based rows of the supplementary tables, and only those
#'
#' \code{analyze_cross_omics_enrichment()} accepts these tables precisely
#' because they are rank-based, so that is checked rather than assumed. A table
#' with no `method` column cannot be confirmed and is refused with a warning --
#' the same fail-closed rule the ORA figure applies to ORA membership. Rows of
#' any other method are dropped quietly; a table left with none is dropped.
#'
#' @param rank_tables Named list of data frames, or NULL.
#' @return Named list holding only the usable tables, each reduced to its
#'   rank-based rows; an empty list when nothing is usable.
#' @keywords internal
.usable_rank_tables <- function(rank_tables) {
    if (is.null(rank_tables) || length(rank_tables) == 0) return(list())

    out <- list()
    for (om in names(rank_tables)) {
        df <- rank_tables[[om]]
        if (!is.data.frame(df) || nrow(df) == 0) next
        if (!"method" %in% names(df)) {
            warning("Ignoring the rank-based table supplied for ", om,
                    ": it has no method column, so its rows cannot be ",
                    "confirmed as rank-based")
            next
        }
        is_rank <- !is.na(df$method) &
            tolower(as.character(df$method)) %in% .RANK_BASED_METHODS
        if (any(is_rank)) out[[om]] <- df[is_rank, , drop = FALSE]
    }
    out
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
#' Each layer contributes one kind of test, chosen by
#' \code{select_layer_method_rows()}: its rank-based rows (fgsea / GSEA) where it
#' has any, otherwise its ORA rows. Choosing one column for the whole table used
#' to decide this by accident -- `pvalue` was taken before `pval`, so a layer
#' holding ORA rows in `pvalue` beside fgsea rows in `pval` contributed its ORA
#' rows only, while a layer scored by GSEA alone contributed GSEA, and Stouffer
#' then combined unlike tests.
#'
#' @param pathway_tables Named list of per-omics enrichment data frames.
#' @param target_pathways Character vector of join keys to report on, as produced
#'   by \code{pathway_join_key()}.
#' @param omics Character vector naming which layers of \code{pathway_tables} to
#'   merge, in order.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return Data frame with one row per element of \code{target_pathways}: a
#'   \code{norm_id} column, and per layer a \code{pval_<omics>} column holding the
#'   raw p-value and a \code{method_<omics>} column naming the test it came from
#'   (NA where that layer has no p-value for the pathway).
merge_pathway_pvalues <- function(pathway_tables, target_pathways, omics,
                                   kegg_org = NULL) {

    merged <- data.frame(norm_id = target_pathways, stringsAsFactors = FALSE)

    for (om in omics) {
        df <- pathway_tables[[om]]

        contrib <- .layer_contribution(df, kegg_org)
        if (!is.null(contrib$problem)) {
            warning("Cannot identify ", contrib$problem, " column in ", om,
                    " enrichment table")
            next
        }

        # One row per key, as before: contrasts and now also KEGG prefix variants
        # of one pathway collapse to their best p-value -- within the one method
        # chosen, never across methods. The result having unique keys is what
        # keeps the merge below one-to-one, so no layer can multiply another's
        # rows.
        usable <- contrib$keep
        if (!any(usable)) next

        df_agg <- aggregate(list(pval = contrib$pvals[usable]),
                            by = list(norm_id = contrib$keys[usable]),
                            FUN = min)
        colnames(df_agg) <- c("norm_id", paste0("pval_", om))
        df_agg[[paste0("method_", om)]] <- contrib$method

        # Subset to target pathways
        df_sub <- df_agg[df_agg$norm_id %in% target_pathways, , drop = FALSE]

        # Merge
        merged <- merge(merged, df_sub, by = "norm_id", all.x = TRUE)
    }

    merged
}


#' Methods the cross-omics merge treats as rank-based
#'
#' `fgsea` is what \code{run_pathway_analysis()} and \code{run_compound_gsea()}
#' write; `gsea` is what \code{run_gsea_kegg()} writes. Both score a ranked list
#' of every measured feature, with no significance cutoff upstream, which is why
#' they are preferred: a layer can carry evidence here without any feature
#' passing its DE threshold.
#' @keywords internal
.RANK_BASED_METHODS <- c("fgsea", "gsea")


#' Raw p-value for each row, whichever column carries it
#'
#' Resolved per row, like \code{.ora_adjusted_p_values()} and for the same
#' reason: \code{extract_enrichment_df()} binds heterogeneous sub-results with
#' \code{dplyr::bind_rows()}, so one layer can hold fgsea rows carrying `pval`
#' beside ORA rows carrying `pvalue`. Picking one column for the table reads NA
#' for every row that used the other.
#'
#' Only when a table carries neither raw column does this fall back to an
#' adjusted one, which is what the merge did before; that fallback is kept for
#' tables from producers that export nothing else, not as a preference.
#'
#' @param df Enrichment data frame for one omics layer.
#' @return Numeric vector, one per row of \code{df} (NA where no accepted column
#'   carries a value), or NULL when the table has no p-value column at all.
#' @examples
#' .raw_p_values(data.frame(pvalue = c(0.01, NA), pval = c(NA, 0.02)))  # 0.01 0.02
#' @keywords internal
.raw_p_values <- function(df) {
    cols <- intersect(c("pvalue", "pval"), names(df))
    if (length(cols) == 0) cols <- intersect(c("p.adjust", "padj"), names(df))[1]
    if (length(cols) == 0 || all(is.na(cols))) return(NULL)

    vals <- rep(NA_real_, nrow(df))
    for (col in cols) {
        # Converted only when it is not already numeric: as.numeric() on a factor
        # returns level codes, and round-tripping a double through as.character()
        # would cost precision the p-values cannot spare.
        v <- df[[col]]
        if (!is.numeric(v)) v <- suppressWarnings(as.numeric(as.character(v)))
        fill <- is.na(vals) & !is.na(v)
        vals[fill] <- v[fill]
    }
    vals
}


#' The rows of one layer that can contribute to the meta-analysis
#'
#' The single place the cross-omics analysis decides which of a layer's rows
#' count: \code{merge_pathway_pvalues()} aggregates them and
#' \code{analyze_cross_omics_enrichment()} builds its candidate pathways from
#' them, so a pathway only rows of an unselected method support can never
#' become a candidate with no layer behind it.
#'
#' Identity is resolved per row (\code{pathway_join_key()}), so a table mixing
#' collections keys each row on what that row carries; the raw p-value is read
#' per row (\code{.raw_p_values()}); one method is chosen for the layer
#' (\code{select_layer_method_rows()}). A row contributes when it belongs to
#' that method and has both a key and a p-value.
#'
#' @param df Enrichment data frame for one omics layer.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return List with \code{problem} -- NULL, or \code{"pathway"} /
#'   \code{"p-value"} when the table has no usable column of that kind, in
#'   which case nothing else is returned -- and otherwise \code{keys},
#'   \code{pvals} (one per row), \code{keep} (logical, the contributing rows)
#'   and \code{method} (the test they come from).
#' @examples
#' df <- data.frame(ID = c("map00010", "map00020"), method = c("fgsea", "ora"),
#'                  pval = c(0.2, NA), pvalue = c(NA, 0.001))
#' .layer_contribution(df)$keep   # TRUE FALSE: map00020 has ORA rows only
#' @keywords internal
.layer_contribution <- function(df, kegg_org = NULL) {
    keys <- pathway_join_key(df, kegg_org)
    if (all(is.na(keys))) return(list(problem = "pathway"))

    pvals <- .raw_p_values(df)
    if (is.null(pvals)) return(list(problem = "p-value"))

    chosen <- select_layer_method_rows(df, pvals)
    list(problem = NULL, keys = keys, pvals = pvals,
         keep = chosen$keep & !is.na(keys) & !is.na(pvals),
         method = chosen$method)
}


#' What choosing one method left out of a layer
#'
#' Rows of the methods not chosen are left out on purpose -- a layer
#' contributes one kind of test -- but when that removes every row of a
#' contrast or of a gene-set collection, the layer's evidence for it is gone
#' from the run-level meta-analysis, and that is said rather than left for the
#' reader to infer from a missing cell.
#'
#' @param df Enrichment data frame for one omics layer.
#' @param contrib Its \code{.layer_contribution()} result.
#' @return List with \code{n_rows} (rows with a key and p-value left out),
#'   \code{contrasts} and \code{collections} (values of the \code{contrast} /
#'   \code{database} columns that only left-out rows carry; empty when the
#'   column is absent).
#' @keywords internal
.method_selection_losses <- function(df, contrib) {
    none <- list(n_rows = 0L, contrasts = character(0), collections = character(0))
    if (!is.null(contrib$problem)) return(none)

    left_out <- !contrib$keep & !is.na(contrib$keys) & !is.na(contrib$pvals)
    only_left_out <- function(col) {
        if (!col %in% names(df)) return(character(0))
        v <- as.character(df[[col]])
        lost <- unique(v[left_out & !is.na(v)])
        sort(setdiff(lost, unique(v[contrib$keep])))
    }
    list(n_rows = sum(left_out),
         contrasts = only_left_out("contrast"),
         collections = only_left_out("database"))
}


#' Choose one kind of test for a layer's rows
#'
#' Rank-based rows (see \code{.RANK_BASED_METHODS}) win when the layer has any
#' with a p-value; otherwise ORA rows; otherwise every row, as before. Rows whose
#' `method` is missing are left out whenever a known method was chosen, because
#' they cannot be shown to belong to it -- \code{dplyr::bind_rows()} NA-fills the
#' column when a layer stacks tables that do not all carry it.
#'
#' A table with no `method` column at all is used whole and labelled
#' "unspecified". Refusing it, as the ORA figure does, would drop layers that
#' custom producers hand over with p-values and nothing else, and the merge has
#' accepted those since it was written.
#'
#' The choice is made per layer, not per contrast: the run-level merge collapses
#' contrasts to their best p-value, and that minimum must not range over two
#' different tests either. A layer whose contrasts were scored by different
#' methods therefore contributes its rank-based contrasts only at run level; the
#' per-contrast calls still see each contrast's own rows.
#'
#' @param df Enrichment data frame for one omics layer.
#' @param pvals Raw p-values for its rows, from \code{.raw_p_values()}.
#' @return List with `keep` (logical, one per row of \code{df}) and `method`
#'   (single string naming the test the kept rows come from; several unknown
#'   methods are joined with "+").
#' @examples
#' df <- data.frame(ID = c("map00010", "map00010"), method = c("ora", "fgsea"),
#'                  pvalue = c(0.001, NA), pval = c(NA, 0.2))
#' select_layer_method_rows(df, .raw_p_values(df))$method   # "fgsea"
select_layer_method_rows <- function(df, pvals) {
    n <- nrow(df)
    if (!"method" %in% names(df)) {
        return(list(keep = rep(TRUE, n), method = "unspecified"))
    }

    m <- tolower(trimws(as.character(df$method)))
    m[!is.na(m) & !nzchar(m)] <- NA_character_
    has_p <- !is.na(pvals)

    rank_rows <- !is.na(m) & m %in% .RANK_BASED_METHODS & has_p
    if (any(rank_rows)) {
        return(list(keep = rank_rows,
                    method = paste(sort(unique(m[rank_rows])), collapse = "+")))
    }

    ora_rows <- !is.na(m) & m == "ora" & has_p
    if (any(ora_rows)) {
        return(list(keep = ora_rows, method = "ora"))
    }

    known <- sort(unique(m[!is.na(m)]))
    list(keep = rep(TRUE, n),
         method = if (length(known) == 0) "unspecified" else paste(known, collapse = "+"))
}


#' Candidate pathways for the cross-omics ORA figure
#'
#' Every pathway the per-layer tables carry, whatever the meta-analysis chose
#' to combine: this figure shows each layer's own adjusted ORA p-value, and
#' \code{build_ora_adjusted_p_matrix()} decides which rows of a table are ORA.
#' Rows with no ORA value for any layer are dropped by the figure itself, so
#' the candidates need not pre-filter on method. The KEGG class exclusion is
#' applied to the keys as it is to the meta-analysis candidates.
#'
#' @param pathway_tables Named list of per-omics enrichment data frames, with
#'   the class exclusion already applied to their rows.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @param exclude Excluded KEGG classes (\code{.excluded_pathway_classes()}).
#' @return Character vector of unique join keys, as \code{pathway_join_key()}
#'   produces them.
#' @examples
#' .ora_figure_candidates(list(proteomics = data.frame(
#'     ID = c("map00010", "GO:0006096"), method = c("fgsea", "ora"))))
#' @keywords internal
.ora_figure_candidates <- function(pathway_tables, kegg_org = NULL,
                                   exclude = character(0)) {
    keys <- unique(unlist(lapply(pathway_tables, function(df) {
        if (!is.data.frame(df) || nrow(df) == 0) return(NULL)
        k <- pathway_join_key(df, kegg_org)
        k[!is.na(k)]
    }), use.names = FALSE))
    if (length(keys) == 0) return(character(0))
    if (length(exclude) > 0) {
        keys <- keys[keep_kegg_pathways(keys, exclude = exclude, kegg_org = kegg_org,
                                        label = "cross-omics ORA pathways")]
    }
    keys
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


#' Count the layers that nominally support each pathway
#'
#' Nominal support is how many layers gave the pathway a raw p-value below
#' \code{alpha}. It is a display-ranking quantity and nothing else: it defines
#' no significance threshold, belongs to no tested family, and never reaches a
#' p-value, an adjustment or the exported table.
#'
#' The rule is exactly \code{p < alpha}, so a p-value sitting on the threshold
#' is not support. A layer with no p-value for the pathway contributes zero:
#' absent evidence is not evidence, and the table cannot distinguish a pathway
#' that was untestable in a layer from one tested and dropped before it got
#' here.
#'
#' Counted from the same `pval_*` columns the figures then display, so the
#' ranking and the picture are derived from one thing.
#'
#' @param meta_results Meta-analysis table carrying `pval_<layer>` columns.
#' @param alpha Raw-p threshold for nominal support. Local to display ranking
#'   and deliberately not configurable.
#' @return Integer vector, one count per row; all zero when the table carries
#'   no `pval_*` columns.
#' @keywords internal
.nominal_support <- function(meta_results, alpha = 0.05) {
    pval_cols <- grep("^pval_", names(meta_results), value = TRUE)
    support <- integer(nrow(meta_results))

    for (col in pval_cols) {
        p <- suppressWarnings(as.numeric(meta_results[[col]]))
        # Written out rather than leaning on na.rm: "NA contributes zero" is
        # the contract, and it should be visible in the line that implements it.
        support <- support + as.integer(!is.na(p) & p < alpha)
    }

    support
}


#' Choose the rows a cross-omics figure shows
#'
#' The figures exist to show where the layers agree, and ordering by combined
#' p-value alone does not do that: the layer with the largest gene-set
#' collection contributes many single-layer pathways with very small p-values,
#' and they fill every slot. The pathways several layers support -- the point of
#' the figure -- rank below them and never appear.
#'
#' Counting the layers that merely *held* a p-value does not fix it either. A
#' pathway every layer measured and none found anything in outranks one
#' genuinely significant in two, because the count rewards coverage rather than
#' agreement. The meta figures therefore rank on **nominal support** -- the
#' layers whose raw p-value for the pathway is below 0.05, per
#' \code{.nominal_support()} -- and break ties on the combined p-value, then on
#' a deterministic identity.
#'
#' Nominal support is a display-ranking definition, not a significance
#' threshold: it is counted on the raw per-layer p-values the figure itself
#' shows, so what puts a row at the top is what the reader can see in it.
#'
#' This is selection for display only. The meta-analysis table keeps its own
#' combined_pval ordering, and no p-value, adjustment or membership is
#' touched -- nor is the caller's table, which is copied on assignment rather
#' than edited in place.
#'
#' The ranking columns are named rather than fixed, because the ordering is
#' the reusable part and the quantities are not. The ORA figure, whose evidence
#' is a different table with different missingness, passes its own pair: it
#' must not be ranked on raw-p meta-analysis columns, and it does not supply an
#' identity, so it keeps the incoming-order tie-break it has always had.
#'
#' @param meta_results Meta-analysis table, as
#'   \code{stouffer_combined_pvalues()} returns it, or any table carrying
#'   \code{count_col} and \code{score_col}.
#' @param top_n Number of rows to keep.
#' @param count_col Column holding the support count; more is better. Absent,
#'   every row counts as one and the order is left alone.
#' @param score_col Column holding the score that breaks ties within a count;
#'   smaller is better. Absent, nothing breaks them but the incoming order.
#' @param id_col Column holding a stable identity, breaking ties the first two
#'   keys leave. NULL, or naming a column the table does not carry, falls back
#'   to the incoming order -- which is what every caller predating this
#'   argument relies on.
#' @return The selected rows of \code{meta_results}, in display order.
#' @examples
#' # Two layers nominally support 00020; 00010 has a far smaller combined p but
#' # only one layer under 0.05, so it does not lead the figure.
#' meta <- data.frame(norm_id = c("00010", "00020"),
#'                    pval_transcriptomics = c(1e-9, 0.01),
#'                    pval_proteomics = c(0.9, 0.02),
#'                    combined_pval = c(1e-9, 1e-3))
#' meta$n_nominal_support <- .nominal_support(meta)
#' select_multi_omics_pathways(meta, top_n = 2,
#'                             count_col = "n_nominal_support",
#'                             id_col = "norm_id")$norm_id   # "00020" first
select_multi_omics_pathways <- function(meta_results, top_n = 30,
                                         count_col = "n_omics",
                                         score_col = "combined_pval",
                                         id_col = NULL) {
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

    # An identity breaks what the first two keys leave, so the figure does not
    # depend on the order rows happened to be bound in. Without one -- the ORA
    # caller supplies none -- the incoming order decides, as it always has;
    # that order is itself sorted by combined p-value, so it is reproducible
    # run to run.
    last_key <- if (!is.null(id_col) && id_col %in% names(meta_results)) {
        as.character(meta_results[[id_col]])
    } else {
        seq_len(nrow(meta_results))
    }

    ord <- order(-counts, scores, last_key, na.last = TRUE)
    meta_results[utils::head(ord, min(top_n, nrow(meta_results))), , drop = FALSE]
}


#' Filename-safe slug for a gene-set collection
#'
#' Must be injective over the collection names \code{classify_pathway_collection()}
#' returns, or two collections would write one file and the second would be read
#' as the first. Those names are drawn from a fixed vocabulary of alphanumeric
#' words, so collapsing runs of non-alphanumerics is one-to-one over them.
#'
#' @param collection Collection name.
#' @return Filename-safe string.
#' @keywords internal
.collection_slug <- function(collection) {
    gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", collection))
}


#' Which layers scored anything in this slice of the meta table
#'
#' A per-collection figure shows only the layers that carry a p-value in that
#' collection. Drawing an all-empty column instead would invite the reader to
#' read it as "tested here and found nothing", which is the one inference these
#' figures must not support.
#'
#' Note what this can and cannot distinguish: a layer drops out both when it
#' has no annotation in the collection and when it was scored there but no
#' result was retained into these tables. The producers keep only their own
#' hits, so absence here cannot establish absence of annotation -- the figure
#' legends say so rather than claiming the stronger reading.
#'
#' @param meta_slice Meta-analysis rows for one collection.
#' @param omics Character vector of layer names, in display order.
#' @return The subset of \code{omics} with at least one non-NA p-value here.
#' @keywords internal
.layers_with_values <- function(meta_slice, omics) {
    omics[vapply(omics, function(om) {
        cn <- paste0("pval_", om)
        cn %in% names(meta_slice) && any(!is.na(meta_slice[[cn]]))
    }, logical(1))]
}


#' Remove cross-omics heatmaps a previous run left in this directory
#'
#' These two figures used to be written once each under a fixed name, so every
#' run overwrote its predecessor. They are now written once per gene-set
#' collection, and which collections a run produces depends on its data -- so a
#' collection that has since disappeared, or the pre-existing combined figure,
#' would survive and still be discovered by the report's glob and shown as this
#' run's evidence.
#'
#' Deliberately narrow: only these two families, only in the directory given.
#' The wider question of output lifecycle across this module is tracked
#' separately and is not settled here.
#'
#' @param out_dir Directory the current run is about to write into.
#' @return Invisibly, the paths removed.
#' @keywords internal
.clear_collection_heatmaps <- function(out_dir) {
    stale <- c(
        file.path(out_dir, c("cross_omics_pathway_heatmap.png",
                             "cross_omics_ora_heatmap.png")),
        list.files(out_dir,
                   pattern = "^cross_omics_(pathway|ora)_heatmap_.*\\.png$",
                   full.names = TRUE)
    )
    stale <- unique(stale[file.exists(stale)])
    if (length(stale) > 0) unlink(stale)
    invisible(stale)
}


#' Plot cross-omics pathway heatmap
#'
#' @param meta_results Meta-analysis rows to draw, already restricted to one
#'   gene-set collection where the caller drew per collection.
#' @param omics Character vector of layer names, in display order. Only these
#'   layers become columns, so a caller passing the layers that scored
#'   something keeps empty columns out of the figure.
#' @param top_n Number of pathways to show.
#' @param collection Collection name to name in the title, or NULL for a figure
#'   drawn across collections.
#' @return Whatever the drawing call returned; called for the figure it puts on
#'   the active device, not for its value.
plot_cross_omics_pathway_heatmap <- function(meta_results, omics, top_n = 30,
                                             collection = NULL) {

    # Named in the title where the caller drew one collection, so a reader
    # holding two of these figures can tell which gene sets each one scored.
    main_title <- function(base) {
        if (is.null(collection)) base else paste0(collection, " gene sets\n", base)
    }

    # Rows that several layers nominally support lead. The count comes off the
    # full table's pval_* columns, before selection -- the matrix below is built
    # from the rows already chosen, so it cannot be what chooses them.
    meta_results$n_nominal_support <- .nominal_support(meta_results)
    top_pathways <- select_multi_omics_pathways(
        meta_results, top_n,
        count_col = "n_nominal_support",
        id_col = "norm_id")

    # Restricted to the layers the caller named. A per-collection figure passes
    # only the layers that scored something in that collection, and taking every
    # pval_* column instead would put back the all-empty column whose absence is
    # the whole point -- an empty column reads as "tested here and found
    # nothing", which is what these figures must not imply.
    pval_cols <- intersect(paste0("pval_", omics), names(top_pathways))
    if (length(pval_cols) == 0) {
        pval_cols <- grep("^pval_", names(top_pathways), value = TRUE)
    }
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
                           main = main_title("Cross-Omics Pathway Enrichment (-log10 p-value)"),
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
                        main = main_title("Cross-Omics Pathway Enrichment"),
                        col = base_cols, zlim = zl)
            } else {
                # heatmap() refuses this shape outright; see .draw_small_heatmap().
                .draw_small_heatmap(log_pval_matrix,
                                    main = main_title("Cross-Omics Pathway Enrichment"),
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
#' @param collection Collection name to name in the title, or NULL for a figure
#'   drawn across collections.
#' @return Invisibly, the -log10 matrix that was drawn, in display order and
#'   with display labels for row names -- or NULL when there was nothing to
#'   draw. Returned so that which rows the figure leads with can be checked
#'   without reading pixels, or guessing which of the two drawing branches a
#'   machine took.
plot_cross_omics_ora_heatmap <- function(padj_matrix, pathway_tables = NULL,
                                          kegg_org = NULL, top_n = 30,
                                          collection = NULL) {

    main_title <- function(base) {
        if (is.null(collection)) base else paste0(collection, " gene sets\n", base)
    }

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
                           main = main_title("Cross-Omics ORA Evidence (-log10 adjusted p-value)"),
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
                        main = main_title("Cross-Omics ORA Evidence"),
                        col = base_cols, zlim = zl)
            } else {
                # heatmap() refuses this shape outright; see .draw_small_heatmap().
                .draw_small_heatmap(log_padj_matrix,
                                    main = main_title("Cross-Omics ORA Evidence"),
                                    col = base_cols, zlim = zl)
            }
        })
    }

    invisible(log_padj_matrix)
}


#' Plot enrichment dot plot
plot_enrichment_dotplot <- function(meta_results, omics, top_n = 20) {

    # Same ranking as the heatmap: nominal support counted on the full table's
    # pval_* columns before any row is dropped. See .nominal_support().
    meta_results$n_nominal_support <- .nominal_support(meta_results)
    top <- select_multi_omics_pathways(
        meta_results, top_n,
        count_col = "n_nominal_support",
        id_col = "norm_id")

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
#' By default, GSEA over every feature ranked by its signed DIABLO loading or
#' MOFA2 weight (\code{run_loadings_gsea()}, in 07e_loadings_gsea.R). With
#' `modes.multiomics.enrichment.loadings.method: "ora"` it instead takes the
#' top features by absolute loading or weight and runs ORA (KEGG) for each
#' component/factor per omics view, as it always did.
#'
#' @param integration_res Output from mod_multiomics_integration()
#' @param harmonization_res Output from mod_multiomics_harmonization()
#' @param config Full config object
#' @param out_dir Output directory for results
#' @param top_n Number of top features per component/factor to use (default 50);
#'   ORA only.
#' @return List with diablo and mofa enrichment results
run_loadings_enrichment <- function(integration_res, harmonization_res,
                                     config, out_dir, top_n = 50) {

    message("\n=== Loadings-based Geneset Enrichment ===\n")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # GSEA unless the config asks for ORA. The record tells the report which
    # test's files in these directories belong to this run.
    method <- loadings_enrichment_method(config)
    record_loadings_method(out_dir, method)
    if (identical(method, "gsea")) {
        return(run_loadings_gsea(integration_res, harmonization_res, config, out_dir))
    }

    # The record now says "ora", so a previous run's ORA files -- for an
    # integration that is absent this time, or components that no longer
    # produce rows -- would be reported as this run's.
    .clear_loadings_ora_all(out_dir)

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
                config = config,
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
                config = config,
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
                                            out_dir, top_n = 50, exclude_classes = NULL,
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
run_mofa_weights_enrichment <- function(mofa_results, harmonization_res,
                                         organism, kegg_org, org_db,
                                         out_dir, top_n = 50, exclude_classes = NULL,
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
                                 organism, kegg_org, org_db, config = NULL) {

    # Resolve IDs using the actual omics type
    resolved_ids <- resolve_gene_n_ids(feature_ids, harmonization_res, omics_type)

    # No KEGG code or OrgDb -- true of every non-model organism. The per-omic
    # GMTs are the same gene sets the DE-driven enrichment already uses, so fall
    # back to those rather than skipping the view entirely. With no config there
    # is nothing to fall back to and this returns NULL, as it always did.
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


#' Loadings enrichment against the per-omic custom GMTs
#'
#' The fallback for \code{enrich_feature_list()} when the organism has no KEGG
#' code or OrgDb, which is the case for every non-model organism. Those runs
#' produced loadings enrichment for metabolomics only: the gene views returned
#' NULL and the section simply had nothing in it, with no indication that the
#' gene sets to test against were configured and sitting unused.
#'
#' Runs the same over-representation test the DE-driven enrichment runs, on the
#' same gene sets, through \code{gmt_to_term2gene()} and
#' \code{run_multi_ora_enricher()} -- no enrichment logic of its own.
#'
#' @param resolved_ids Top-loading feature IDs of one view, already through
#'   \code{resolve_gene_n_ids()}.
#' @param omics_type One of "transcriptomics" / "proteomics"; anything else
#'   returns NULL, since only those two carry a configured GMT.
#' @param harmonization_res Harmonization result, for the background feature set.
#' @param config Full pipeline config, or NULL to decline.
#' @return Data frame in the same column shape as the KEGG branch above, so
#'   \code{.rbind_fill()} can stack gene-view and compound-view results into one
#'   loadings table; NULL when there is nothing to test against.
enrich_feature_list_gmt <- function(resolved_ids, omics_type,
                                    harmonization_res, config) {
    if (is.null(config)) return(NULL)

    cfg_key <- c(transcriptomics = "rna", proteomics = "proteomics")[omics_type]
    if (is.na(cfg_key)) return(NULL)

    # gmt_file may be a single path or a YAML list (GO + KEGG + Pfam);
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

    # Background = every feature of this view that survived preprocessing, put
    # through the same resolution as the query so the two share a namespace.
    universe <- NULL
    pre_data <- harmonization_res$inputs[[omics_type]]
    if (!is.null(pre_data) && !is.null(pre_data$expr_work)) {
        universe <- unique(resolve_gene_n_ids(rownames(pre_data$expr_work),
                                              harmonization_res, omics_type))
    }

    # Fail closed. Handing run_multi_ora_enricher() a NULL universe is not a
    # degraded run, it is a different test: enricher() then takes its background
    # from TERM2GENE, so the universe becomes the union of the gene sets rather
    # than what this view measured, and every p-value comes out inflated. The
    # result would look like an ordinary enrichment table. One guard covers all
    # the ways the background can come up empty -- no preprocessed data for the
    # view, no expr_work, unnamed rows, or nothing left after resolution.
    if (length(universe) == 0) {
        message("    Loadings ORA (GMT): ", omics_type,
                " has no usable background feature set ",
                "(harmonization_res$inputs[[\"", omics_type,
                "\"]]$expr_work), skipping rather than testing against the ",
                "gene sets themselves")
        return(NULL)
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
#' views and appends the view name to any it finds in more than one, so
#' "GENE_12" and "GENE_12_transcriptomics" both have to resolve; matching only
#' the bare form left a large share of the MOFA weights unresolved and the
#' enrichment that follows working from a partial feature list.
#'
#' @param feature_ids Character vector of feature IDs, GENE_N or native.
#' @param harmonization_res Harmonization result carrying
#'   \code{gene_protein_mapping}.
#' @param omics_type One of "transcriptomics" / "proteomics"; anything else is
#'   returned untouched (metabolomics uses feature_N, not GENE_N).
#' @return Character vector of resolved IDs; entries with no mapping are dropped.
resolve_gene_n_ids <- function(feature_ids, harmonization_res, omics_type) {
    gpm <- harmonization_res$gene_protein_mapping
    if (is.null(gpm)) return(feature_ids)

    # Check if IDs look like GENE_N, with or without a MOFA view suffix
    is_gene_n <- grepl("^GENE_\\d+(_.+)?$", feature_ids)
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
    # Strip any view suffix before the lookup; the mapping is keyed on the bare
    # GENE_N, which is what the harmonized MAE assigned.
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
