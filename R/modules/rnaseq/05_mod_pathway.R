#' RNA-seq Pathway Analysis Module
#'
#' Orchestrates pathway/enrichment analysis for RNA-seq data.
#'
#' Two modes:
#'   1. Local enrichment (if config$modes$rna$enrichment$annotation_dir is set):
#'      Loads precomputed KEGG/GO tables from local files, runs GSEA with
#'      multiple ranking methods. Cluster-based ORA is deferred to Phase 2.
#'   2. Online fallback (if annotation_dir is not set):
#'      Original behavior — organism detection, gene annotation, online gene
#'      set loading (OrgDb/KEGG/biomaRt), fGSEA + ORA per contrast.
#'
#' @param de_res  Result from mod_rnaseq_de() — list(dds, tables)
#' @param pre     Preprocessed data list (with expr_filt, meta)
#' @param config  Full pipeline config
#' @param out_dir Output directory for the RNA mode (e.g. .../rna)
#' @return List with annotation, pathway_results, and plot_files
#' @export
mod_rnaseq_pathway <- function(de_res, pre, config, out_dir) {

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
            de_tables    = de_tables,
            feature_ids  = all_gene_ids,
            enr_cfg      = enr_cfg,
            config       = config,
            out_dir      = out_dir
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
# LOCAL ENRICHMENT PATH (Phase 1)
# ==============================================================================

#' Internal: run enrichment using local precomputed tables
#'
#' @param de_tables Named list of DE tables
#' @param feature_ids All unique feature IDs
#' @param enr_cfg enrichment config section
#' @param config Full config (for future Phase 2 clustering access)
#' @param out_dir Output directory
#' @return List with annotation, pathway_results, plot_files
#' @noRd
.run_local_enrichment <- function(de_tables, feature_ids, enr_cfg, config, out_dir) {

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

    # ------------------------------------------------------------------
    # 2. ORA: cluster-based (Phase 2 — deferred)
    # ------------------------------------------------------------------
    # Cluster-based ORA requires clustering results which are not yet wired
    # into this module (pipeline DAG change needed). Emit the agreed warning.
    message("NOTE: Cluster-based ORA requires clustering results. ",
            "Enable clustering and wire clustering_res into the pathway module ",
            "to run ORA (Phase 2). GSEA will proceed.")

    # ------------------------------------------------------------------
    # 3. Build ranked gene lists
    # ------------------------------------------------------------------
    message("Building ranked gene lists...")
    ranked_genes <- build_ranked_gene_lists(de_tables)

    n_contrasts <- length(de_tables)
    n_methods <- length(ranked_genes)
    n_total <- sum(vapply(ranked_genes, length, integer(1)))
    message("  Built ", n_total, " ranked lists (",
            n_methods, " methods x ", n_contrasts, " contrasts + any_contrast)")

    # ------------------------------------------------------------------
    # 4. Run GSEA across all combinations
    # ------------------------------------------------------------------
    enrich_dir <- file.path(out_dir, "Enrichment")
    gsea_dir   <- file.path(enrich_dir, "GSEA")

    pval_cutoff <- enr_cfg$gsea_pvalue_cutoff %||% enr_cfg$pvalue_cutoff %||% 0.05
    padj_method <- enr_cfg$gsea_padj_method   %||% enr_cfg$padj_method   %||% "fdr"

    message("Running GSEA (pvalueCutoff=", pval_cutoff,
            ", pAdjustMethod=", padj_method, ")...")

    gsea_out <- run_gsea_all(
        ranked_genes  = ranked_genes,
        local_tables  = local_tables,
        pvalueCutoff  = pval_cutoff,
        pAdjustMethod = padj_method,
        output_dir    = gsea_dir
    )

    pathway_results <- gsea_out$results
    plot_files      <- gsea_out$plot_files

    # ------------------------------------------------------------------
    # 5. Summary
    # ------------------------------------------------------------------
    n_result_dfs <- sum(vapply(pathway_results, function(pr) {
        if (is.list(pr)) length(pr) else 0L
    }, integer(1)))
    message("\n=== Local enrichment complete ===")
    message("  GSEA result sets: ", n_result_dfs)
    message("  Plot files: ", length(plot_files))
    message("  ORA: skipped (Phase 2 — requires clustering wiring)")

    list(
        annotation      = NULL,
        pathway_results = pathway_results,
        plot_files      = plot_files
    )
}
