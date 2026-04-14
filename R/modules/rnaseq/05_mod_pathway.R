#' RNA-seq Pathway Analysis Module
#'
#' Orchestrates pathway/enrichment analysis for RNA-seq data.
#'
#' Two modes:
#'   1. Local enrichment (if config$modes$rna$enrichment$annotation_dir is set):
#'      Loads precomputed KEGG/GO tables from local files.
#'      Runs cluster-based ORA (legacy behavior) when clustering results are available.
#'      Runs GSEA with multiple ranking methods as supplementary analysis.
#'   2. Online fallback (if annotation_dir is not set):
#'      Original behavior — organism detection, gene annotation, online gene
#'      set loading (OrgDb/KEGG/biomaRt), fGSEA + ORA per contrast.
#'
#' @param de_res  Result from mod_rnaseq_de() — list(dds, tables)
#' @param pre     Preprocessed data list (with expr_filt, meta)
#' @param config  Full pipeline config
#' @param out_dir Output directory for the RNA mode (e.g. .../rna)
#' @param clustering_res Result from mod_rnaseq_clustering(), or NULL
#' @return List with annotation, pathway_results, and plot_files
#' @export
mod_rnaseq_pathway <- function(de_res, pre, config, out_dir, clustering_res = NULL) {

    rna_cfg <- config$modes$rna
    pw_cfg  <- rna_cfg$pathway    %||% list()
    enr_cfg <- rna_cfg$enrichment %||% list()

    # Skip entirely if pathway analysis is explicitly disabled
    if (isFALSE(pw_cfg$enabled)) {
        message("Pathway analysis disabled in config (pathway.enabled: false)")
        return(list(annotation = NULL, pathway_results = list(), plot_files = list()))
    }

    de_tables <- de_res$tables
    if (is.null(de_tables) || length(de_tables) == 0) {
        warning("No DE tables available for pathway analysis")
        return(list(annotation = NULL, pathway_results = list(), plot_files = list()))
    }

    # ------------------------------------------------------------------
    # Collect all gene IDs from DE tables (shared by both paths)
    # ------------------------------------------------------------------
    all_gene_ids <- unique(unlist(lapply(de_tables, function(x) x$FeatureID)))

    # Strip "Gene:" prefix if present (e.g. "Gene:WBGene00000001" -> "WBGene00000001")
    has_prefix <- grepl("^Gene:", all_gene_ids)
    if (mean(has_prefix) > 0.5) {
        message("Stripping 'Gene:' prefix from feature IDs for annotation/pathway")
        for (cn in names(de_tables)) {
            de_tables[[cn]]$FeatureID <- sub("^Gene:", "", de_tables[[cn]]$FeatureID)
        }
        all_gene_ids <- sub("^Gene:", "", all_gene_ids)
    }

    # ------------------------------------------------------------------
    # Route: local enrichment vs. online fallback
    # ------------------------------------------------------------------
    annotation_dir <- enr_cfg$annotation_dir
    use_local <- !is.null(annotation_dir) && nzchar(annotation_dir)

    if (use_local) {
        return(.run_local_enrichment(
            de_tables      = de_tables,
            feature_ids    = all_gene_ids,
            enr_cfg        = enr_cfg,
            config         = config,
            out_dir        = out_dir,
            clustering_res = clustering_res
        ))
    }

    # ==================================================================
    # ONLINE FALLBACK — existing behavior, unchanged
    # ==================================================================
    ann_cfg <- rna_cfg$annotation %||% list()

    # Auto-detect organism
    organism <- ann_cfg$organism %||% "auto"
    if (tolower(organism) == "auto") {
        org_detect <- detect_organism_from_ids(all_gene_ids)
        if (org_detect$confidence %in% c("high", "medium")) {
            organism <- org_detect$organism
            message("Auto-detected organism: ", organism,
                    " (", org_detect$confidence, " confidence)")
        } else {
            warning("Could not auto-detect organism. Pathway analysis may be limited.")
            organism <- "unknown"
        }
    } else {
        organism <- normalize_organism_name(organism)
    }

    # Annotate genes (custom -> biomaRt -> OrgDb fallback)
    annotation_result <- NULL
    if (!isTRUE(ann_cfg$skip_annotation)) {
        anno_config <- list(
            organism = organism,
            annotation = list(
                skip_annotation = FALSE,
                custom_mapping_file = ann_cfg$custom_mapping_file,
                id_type = ann_cfg$id_type %||% "auto",
                fallback_chain = c("custom", "biomart", "orgdb", "keggrest")
            )
        )
        annotation_result <- annotate_genes_v2(all_gene_ids, anno_config, verbose = TRUE)
    }

    annotation_df <- if (!is.null(annotation_result)) annotation_result$annotation else NULL

    # Load gene sets (GO, KEGG, and/or custom GMT)
    databases <- pw_cfg$databases %||% c("GO", "KEGG")
    gmt_file  <- pw_cfg$gmt_file

    gene_sets <- load_gene_sets(
        organism         = organism,
        pathway_database = databases,
        gmt_file         = gmt_file,
        annotation       = annotation_df
    )

    if (length(gene_sets) == 0) {
        message("No gene sets loaded. Skipping pathway analysis.")
        return(list(annotation = annotation_result,
                    pathway_results = list(),
                    plot_files = list()))
    }

    # Run fGSEA / ORA per contrast
    pw_method <- pw_cfg$method   %||% "fgsea"
    pw_min    <- pw_cfg$min_size %||% 10
    pw_max    <- pw_cfg$max_size %||% 500

    pathway_results <- run_pathway_analysis(
        de_tables  = de_tables,
        gene_sets  = gene_sets,
        annotation = annotation_df,
        method     = pw_method,
        min_size   = pw_min,
        max_size   = pw_max
    )

    # Save CSV results + dotplot PNGs
    enrich_dir <- file.path(out_dir, "Enrichment")
    plot_dir   <- file.path(enrich_dir, "plots")

    saved_files <- save_pathway_results(pathway_results, enrich_dir)
    plot_files  <- generate_pathway_plots(pathway_results, plot_dir)

    if (!is.null(annotation_df)) {
        anno_file <- file.path(enrich_dir, "gene_annotation.csv")
        write.csv(annotation_df, anno_file, row.names = FALSE)
        message("Saved gene annotation to: ", anno_file)
    }

    list(
        annotation      = annotation_result,
        pathway_results = pathway_results,
        plot_files      = plot_files
    )
}


# ==============================================================================
# LOCAL ENRICHMENT PATH
# ==============================================================================

#' Internal: run enrichment using local precomputed tables
#'
#' @param de_tables Named list of DE tables
#' @param feature_ids All unique feature IDs
#' @param enr_cfg enrichment config section
#' @param config Full config
#' @param out_dir Output directory
#' @param clustering_res Clustering result from mod_rnaseq_clustering(), or NULL
#' @return List with annotation, pathway_results, plot_files
#' @noRd
.run_local_enrichment <- function(de_tables, feature_ids, enr_cfg, config,
                                  out_dir, clustering_res = NULL) {

    message("\n=== Local Enrichment (offline, table-driven) ===\n")

    annotation_dir <- enr_cfg$annotation_dir

    # Resolve relative path against project dir if needed
    if (!dir.exists(annotation_dir)) {
        project_dir <- config$project$dir
        if (!is.null(project_dir) && nzchar(project_dir)) {
            candidate <- file.path(project_dir, annotation_dir)
            if (dir.exists(candidate)) {
                annotation_dir <- candidate
            }
        }
    }

    # ------------------------------------------------------------------
    # 1. Load local pathway tables
    # ------------------------------------------------------------------
    databases <- enr_cfg$databases %||% c("KEGG", "GO_BP", "GO_MF", "GO_CC")

    message("Loading local pathway tables from: ", annotation_dir)
    local_tables <- load_local_pathway_tables(
        annotation_dir = annotation_dir,
        databases      = databases,
        feature_ids    = feature_ids
    )

    if (length(local_tables) == 0) {
        warning("No local pathway tables loaded. Skipping enrichment.")
        return(list(annotation = NULL, pathway_results = list(), plot_files = list()))
    }

    enrich_dir <- file.path(out_dir, "Enrichment")
    pval_cutoff <- enr_cfg$pvalue_cutoff %||% 0.05
    padj_method <- enr_cfg$padj_method   %||% "fdr"

    pathway_results <- list()
    plot_files      <- list()

    # ------------------------------------------------------------------
    # 2. Cluster-based ORA (primary legacy behavior)
    # ------------------------------------------------------------------
    clusters <- .extract_clusters(clustering_res)

    if (!is.null(clusters) && length(clusters) > 0) {
        message("\n--- Cluster-based ORA ---")
        message("  ", length(unique(clusters)), " clusters, ",
                length(clusters), " genes")

        ora_results <- list()

        for (db_name in names(local_tables)) {
            tbl <- local_tables[[db_name]]

            # Determine type for simplify dispatch: GO or KEGG
            db_type <- if (grepl("^GO", db_name)) "GO" else "KEGG"

            ora_dir <- file.path(enrich_dir, "ORA", db_name)
            ora_file_name <- paste0(db_name, "_cluster_ora")

            message("  ORA: ", db_name, " (type=", db_type, ")")

            allRes <- run_cluster_ora(
                clusters      = clusters,
                TERM2GENE     = tbl$TERM2GENE,
                TERM2NAME     = tbl$TERM2NAME,
                type          = db_type,
                pvalueCutoff  = pval_cutoff,
                pAdjustMethod = padj_method,
                outDir        = ora_dir,
                file_name     = ora_file_name,
                maxCategory   = enr_cfg$max_terms_in_dotplot %||% 1000
            )

            if (length(allRes) == 0) {
                message("    No significant enrichment")
                next
            }

            # allRes is a 4-element list:
            #   [[1]] compareClusterResult, [[2]] simplified (or NULL),
            #   [[3]] enrichment_table, [[4]] enrichment_table_simplify (or NULL)

            # Store processed tables in pathway_results for downstream consumers.
            # Use a structure where the enrichment table (data.frame with padj)
            # is accessible to collect_pipeline_stats() and extract_enrichment_df().
            if (!is.null(allRes[[3]]) && nrow(allRes[[3]]) > 0) {
                ora_df <- allRes[[3]]
                # Ensure downstream-compatible column names
                if ("p.adjust" %in% colnames(ora_df) && !"padj" %in% colnames(ora_df)) {
                    ora_df$padj <- ora_df$p.adjust
                }
                if ("Description" %in% colnames(ora_df) && !"pathway" %in% colnames(ora_df)) {
                    ora_df$pathway <- ora_df$Description
                }
                ora_results[[paste0(db_name, "_ora")]] <- ora_df
            }

            if (!is.null(allRes[[4]]) && nrow(allRes[[4]]) > 0) {
                simp_df <- allRes[[4]]
                if ("p.adjust" %in% colnames(simp_df) && !"padj" %in% colnames(simp_df)) {
                    simp_df$padj <- simp_df$p.adjust
                }
                if ("Description" %in% colnames(simp_df) && !"pathway" %in% colnames(simp_df)) {
                    simp_df$pathway <- simp_df$Description
                }
                ora_results[[paste0(db_name, "_ora_simplify")]] <- simp_df
            }

            n_terms <- if (!is.null(allRes[[3]])) nrow(allRes[[3]]) else 0
            message("    ", n_terms, " enriched terms")
        }

        # Store ORA results under a dedicated key in pathway_results
        if (length(ora_results) > 0) {
            pathway_results[["cluster_ora"]] <- ora_results
        }
    } else {
        message("NOTE: Cluster-based ORA requires clustering results. ",
                "Enable clustering in config to run ORA. GSEA will proceed.")
    }

    # ------------------------------------------------------------------
    # 3. GSEA (supplementary — multiple ranking methods)
    # ------------------------------------------------------------------
    message("\n--- GSEA ---")
    message("Building ranked gene lists...")
    ranked_genes <- build_ranked_gene_lists(de_tables)

    n_contrasts <- length(de_tables)
    n_methods <- length(ranked_genes)
    n_total <- sum(vapply(ranked_genes, length, integer(1)))
    message("  Built ", n_total, " ranked lists (",
            n_methods, " methods x ", n_contrasts, " contrasts + any_contrast)")

    gsea_dir    <- file.path(enrich_dir, "GSEA")
    gsea_pval   <- enr_cfg$gsea_pvalue_cutoff %||% pval_cutoff
    gsea_padj   <- enr_cfg$gsea_padj_method   %||% padj_method
    workers     <- enr_cfg$workers %||% 1

    message("Running GSEA (pvalueCutoff=", gsea_pval,
            ", pAdjustMethod=", gsea_padj, ")...")

    gsea_out <- run_gsea_all(
        ranked_genes  = ranked_genes,
        local_tables  = local_tables,
        pvalueCutoff  = gsea_pval,
        pAdjustMethod = gsea_padj,
        output_dir    = gsea_dir,
        workers       = workers
    )

    # Merge GSEA results into pathway_results (contrast-keyed)
    for (contrast in names(gsea_out$results)) {
        if (is.null(pathway_results[[contrast]])) {
            pathway_results[[contrast]] <- gsea_out$results[[contrast]]
        } else {
            pathway_results[[contrast]] <- c(
                pathway_results[[contrast]], gsea_out$results[[contrast]]
            )
        }
    }
    plot_files <- c(plot_files, gsea_out$plot_files)

    # ------------------------------------------------------------------
    # 4. Summary
    # ------------------------------------------------------------------
    n_result_dfs <- sum(vapply(pathway_results, function(pr) {
        if (is.data.frame(pr)) 1L else if (is.list(pr)) length(pr) else 0L
    }, integer(1)))
    message("\n=== Local enrichment complete ===")
    message("  Result sets: ", n_result_dfs)
    message("  Plot files: ", length(plot_files))
    if (is.null(clusters) || length(clusters) == 0) {
        message("  ORA: skipped (no clustering results)")
    }

    list(
        annotation      = NULL,
        pathway_results = pathway_results,
        plot_files      = plot_files
    )
}


#' Extract cluster assignments from clustering result
#'
#' @param clustering_res Result from mod_rnaseq_clustering(), or NULL
#' @return Named integer vector (gene IDs -> cluster labels), or NULL
#' @noRd
.extract_clusters <- function(clustering_res) {
    if (is.null(clustering_res)) return(NULL)

    # mod_rnaseq_clustering returns list(plots, files, excel_order, objects)
    # objects$clusters is a named integer vector
    clusters <- clustering_res$objects$clusters
    if (is.null(clusters) || length(clusters) == 0) return(NULL)

    clusters
}
