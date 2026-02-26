#' RNA-seq Pathway Analysis Module
#'
#' Orchestrates organism detection, gene annotation, gene set loading,
#' and pathway/enrichment analysis (fGSEA + ORA).
#'
#' @param de_res  Result from mod_rnaseq_de() — list(dds, tables)
#' @param pre     Preprocessed data list (with expr_filt, meta)
#' @param config  Full pipeline config
#' @param out_dir Output directory for the RNA mode (e.g. .../rna)
#' @return List with annotation, pathway_results, and plot_files
#' @export
mod_rnaseq_pathway <- function(de_res, pre, config, out_dir) {

    rna_cfg <- config$modes$rna
    ann_cfg <- rna_cfg$annotation %||% list()
    pw_cfg  <- rna_cfg$pathway   %||% list()

    # Skip entirely if pathway analysis is disabled
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
    # 1. Collect all gene IDs from DE tables
    # ------------------------------------------------------------------
    all_gene_ids <- unique(unlist(lapply(de_tables, function(x) x$FeatureID)))

    # Strip "Gene:" prefix if present (e.g. "Gene:WBGene00000001" -> "WBGene00000001")
    has_prefix <- grepl("^Gene:", all_gene_ids)
    if (mean(has_prefix) > 0.5) {
        message("Stripping 'Gene:' prefix from feature IDs for annotation/pathway")
        # Also strip in each DE table for consistency
        for (cn in names(de_tables)) {
            de_tables[[cn]]$FeatureID <- sub("^Gene:", "", de_tables[[cn]]$FeatureID)
        }
        all_gene_ids <- sub("^Gene:", "", all_gene_ids)
    }

    # ------------------------------------------------------------------
    # 2. Auto-detect organism
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # 3. Annotate genes (custom -> biomaRt -> OrgDb fallback)
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # 4. Load gene sets (GO, KEGG, and/or custom GMT)
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # 5. Run fGSEA / ORA per contrast
    # ------------------------------------------------------------------
    pw_method  <- pw_cfg$method   %||% "fgsea"
    pw_min     <- pw_cfg$min_size %||% 10
    pw_max     <- pw_cfg$max_size %||% 500

    pathway_results <- run_pathway_analysis(
        de_tables  = de_tables,
        gene_sets  = gene_sets,
        annotation = annotation_df,
        method     = pw_method,
        min_size   = pw_min,
        max_size   = pw_max
    )

    # ------------------------------------------------------------------
    # 6. Save CSV results + dotplot PNGs
    # ------------------------------------------------------------------
    enrich_dir <- file.path(out_dir, "Enrichment")
    plot_dir   <- file.path(enrich_dir, "plots")

    saved_files <- save_pathway_results(pathway_results, enrich_dir)
    plot_files  <- generate_pathway_plots(pathway_results, plot_dir)

    # Save annotation table
    if (!is.null(annotation_df)) {
        anno_file <- file.path(enrich_dir, "gene_annotation.csv")
        write.csv(annotation_df, anno_file, row.names = FALSE)
        message("Saved gene annotation to: ", anno_file)
    }

    # Build pathway-colored volcano data
    pw_volcano_enabled <- isTRUE(pw_cfg$pathway_volcano) || is.null(pw_cfg$pathway_volcano)
    if (pw_volcano_enabled && length(pathway_results) > 0) {
        message("Building pathway-colored volcano data for RNA-seq...")
        for (cn in names(de_tables)) {
            volcano_data <- build_pathway_volcano_data(de_tables[[cn]], pathway_results)
            if (!is.null(volcano_data)) {
                volcano_file <- file.path(enrich_dir, sprintf("pathway_volcano_data_%s.csv", cn))
                write.csv(volcano_data, volcano_file, row.names = FALSE)
                message("  Saved pathway volcano data: ", volcano_file)
            }
        }
    }

    list(
        annotation      = annotation_result,
        pathway_results = pathway_results,
        plot_files      = plot_files
    )
}
