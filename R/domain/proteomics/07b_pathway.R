#' Proteomics Pathway Enrichment
#'
#' Bridges proteomics DE results to the shared pathway analysis
#' (run_pathway_analysis from R/core/09_enrichment.R).

# ==============================================================================
# DE TABLE EXTRACTION
# ==============================================================================

#' Extract a per-contrast DE table from the wide summary_df
#'
#' Converts from proteomics summary_df column naming
#' (padj.imputs.<cn>, pvalue.imputs.<cn>, linearFC.imputs.<cn>)
#' to the standard pathway input format:
#' FeatureID, log2FoldChange, pvalue, padj, stat
#'
#' @param summary_df Wide DE summary data.frame
#' @param contrast_name Character: contrast identifier
#' @param config Full config list
#' @return data.frame with FeatureID, log2FoldChange, pvalue, padj, stat
extract_de_table_for_pathway <- function(summary_df, contrast_name, config) {
    cfg <- config$modes$proteomics
    src_id_col <- cfg$de_table$id_col %||% "FeatureID"

    if (!src_id_col %in% colnames(summary_df))
        stop(sprintf("ID column '%s' not found in DE summary table. Available: %s",
                     src_id_col, paste(colnames(summary_df), collapse = ", ")))

    cn <- normalize_contrast_name(contrast_name)
    padj_col <- paste0("padj.imputs.", cn)
    pval_col <- paste0("pvalue.imputs.", cn)
    fc_col   <- paste0("linearFC.imputs.", cn)

    if (!padj_col %in% colnames(summary_df))
        stop(sprintf("Adjusted p-value column '%s' not found for contrast '%s'. Available: %s",
                     padj_col, contrast_name, paste(colnames(summary_df), collapse = ", ")))
    if (!fc_col %in% colnames(summary_df))
        stop(sprintf("Fold-change column '%s' not found for contrast '%s'. Available: %s",
                     fc_col, contrast_name, paste(colnames(summary_df), collapse = ", ")))

    padj_vals <- as.numeric(summary_df[[padj_col]])
    pval_vals <- as.numeric(summary_df[[pval_col]])
    lfc_vals  <- signed_fc_to_log2(as.numeric(summary_df[[fc_col]]))

    # stat = sign(log2FC) * -log10(pvalue)
    stat_vals <- sign(lfc_vals) * -log10(pval_vals + 1e-300)

    de_tbl <- data.frame(
        FeatureID       = summary_df[[src_id_col]],
        log2FoldChange  = lfc_vals,
        pvalue          = pval_vals,
        padj            = padj_vals,
        stat            = stat_vals,
        stringsAsFactors = FALSE
    )

    # Carry annotation columns if present
    for (a in intersect(c("Protein.Names", "Genes", "First.Protein.Description"),
                        colnames(summary_df))) {
        de_tbl[[a]] <- summary_df[[a]]
    }

    de_tbl
}

# ==============================================================================
# PROTEIN-TO-GENE MAPPING
# ==============================================================================

#' Map protein IDs to gene symbols for pathway analysis
#'
#' Uses the Genes column from summary_df (semicolon-separated,
#' takes first symbol per protein).
#'
#' @param de_table  DE table with FeatureID column
#' @param summary_df Original summary data.frame with Genes column
#' @param config    Full config list
#' @return Modified de_table with FeatureID remapped to gene symbols,
#'         original ID stored in ProteinID
map_proteins_to_gene_symbols <- function(de_table, summary_df, config) {
    cfg <- config$modes$proteomics
    src_id_col <- cfg$de_table$id_col %||% "FeatureID"

    if (!"Genes" %in% colnames(summary_df)) {
        message("No 'Genes' column in summary_df — cannot map to gene symbols.")
        return(de_table)
    }

    # Build mapping: protein ID → first gene symbol
    genes_raw <- summary_df[["Genes"]]
    protein_ids <- summary_df[[src_id_col]]

    gene_map <- vapply(genes_raw, function(g) {
        if (is.na(g) || !nzchar(g)) return(NA_character_)
        first <- trimws(strsplit(g, ";")[[1]][1])
        # Strip description after " | " (e.g. "GL50803_X | bZIP family protein" -> "GL50803_X")
        first <- trimws(strsplit(first, "\\|")[[1]][1])
        if (nzchar(first)) first else NA_character_
    }, character(1))
    names(gene_map) <- protein_ids

    # Remap
    de_table$ProteinID <- de_table$FeatureID
    mapped <- gene_map[de_table$FeatureID]
    de_table$FeatureID <- ifelse(is.na(mapped), de_table$FeatureID, mapped)

    # Remove duplicates (keep lowest padj per gene symbol)
    de_table <- de_table[order(de_table$padj, na.last = TRUE), ]
    de_table <- de_table[!duplicated(de_table$FeatureID), ]

    de_table
}

# ==============================================================================
# ORCHESTRATOR
# ==============================================================================

#' Run proteomics pathway enrichment (ORA + fGSEA)
#'
#' @param de_res   DE results list (with summary_df)
#' @param pre      Preprocessed proteomics object
#' @param config   Full config list
#' @param out_dir  Output root directory
#' @return Pathway results list (from run_pathway_analysis)
run_proteomics_pathway <- function(de_res, pre, config, out_dir) {
    cfg <- config$modes$proteomics
    pw_cfg <- cfg$pathway %||% list()

    message("=== Proteomics Pathway Enrichment ===")

    summary_df <- de_res$summary_df
    if (is.null(summary_df)) {
        message("No summary_df available. Skipping pathway analysis.")
        return(NULL)
    }

    # Identify contrasts
    contrasts <- sub("^padj\\.imputs\\.", "",
                     grep("^padj\\.imputs\\.", colnames(summary_df), value = TRUE))

    if (length(contrasts) == 0) {
        message("No contrasts found. Skipping pathway analysis.")
        return(NULL)
    }

    # Build per-contrast DE tables with gene symbols
    de_tables <- lapply(setNames(contrasts, contrasts), function(cn) {
        tbl <- extract_de_table_for_pathway(summary_df, cn, config)
        map_proteins_to_gene_symbols(tbl, summary_df, config)
    })

    # ------------------------------------------------------------------
    # Resolve organism from config (no auto-detect for proteomics)
    # ------------------------------------------------------------------
    ann_cfg  <- cfg$annotation %||% list()
    organism <- ann_cfg$organism %||% "Homo sapiens"
    organism <- normalize_organism_name(organism)
    message("Organism: ", organism)

    # ------------------------------------------------------------------
    # Annotate genes (custom -> biomaRt -> OrgDb fallback)
    # ------------------------------------------------------------------
    annotation_df <- NULL
    if (!isTRUE(ann_cfg$skip_annotation)) {
        all_gene_ids <- unique(unlist(lapply(de_tables, function(x) x$FeatureID)))
        anno_config <- list(
            organism = organism,
            annotation = list(
                skip_annotation = FALSE,
                custom_mapping_file = ann_cfg$custom_mapping_file,
                id_type = ann_cfg$id_type %||% "auto",
                fallback_chain = c("custom", "biomart", "orgdb", "keggrest")
            )
        )
        annotation_result <- tryCatch(
            annotate_genes_v2(all_gene_ids, anno_config, verbose = TRUE),
            error = function(e) {
                message("Gene annotation failed: ", e$message)
                NULL
            }
        )
        if (!is.null(annotation_result)) {
            annotation_df <- annotation_result$annotation
        }
    }

    # ------------------------------------------------------------------
    # Load gene sets (using SYMBOL keytype — proteomics IDs are gene symbols)
    # ------------------------------------------------------------------
    databases <- pw_cfg$databases %||% c("GO", "KEGG")
    gmt_file  <- pw_cfg$gmt_file

    gene_sets <- tryCatch(
        load_gene_sets(organism        = organism,
                       pathway_database = databases,
                       gmt_file         = gmt_file,
                       annotation       = annotation_df,
                       target_id_type   = "symbol"),
        error = function(e) {
            message("Failed to load gene sets: ", e$message)
            list()
        }
    )

    if (length(gene_sets) == 0) {
        message("No gene sets loaded. Skipping pathway analysis.")
        return(NULL)
    }

    # Run shared pathway analysis
    method   <- pw_cfg$method   %||% "both"
    min_size <- pw_cfg$min_size %||% 10
    max_size <- pw_cfg$max_size %||% 500

    pathway_results <- run_pathway_analysis(
        de_tables  = de_tables,
        gene_sets  = gene_sets,
        annotation = annotation_df,
        method     = method,
        min_size   = min_size,
        max_size   = max_size
    )

    # Save results and plots
    enrich_dir <- file.path(out_dir, "Enrichment")
    save_pathway_results(pathway_results, enrich_dir)

    plots_dir <- file.path(enrich_dir, "plots")
    generate_pathway_plots(pathway_results, plots_dir)

    # Save annotation table if generated
    if (!is.null(annotation_df)) {
        anno_file <- file.path(enrich_dir, "gene_annotation.csv")
        write.csv(annotation_df, anno_file, row.names = FALSE)
        message("Saved gene annotation to: ", anno_file)
    }

    message("Proteomics pathway analysis complete.")
    pathway_results
}
