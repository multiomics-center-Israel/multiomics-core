#' Proteomics PPI network analysis module
#'
#' Bridges DE results to run_ppi_network_analysis() by extracting
#' per-contrast tables and mapping to gene symbols.
#' Species ID for STRING is derived from annotation.organism in config.
#'
#' @param de_res  DE results list (with summary_df)
#' @param pre     Preprocessed proteomics object
#' @param config  Full config list
#' @param out_dir Output root directory
#' @return PPI analysis results list or NULL
mod_proteomics_ppi <- function(de_res, pre, config, out_dir) {
    cfg <- config$modes$proteomics
    if (!isTRUE(cfg$ppi$enabled)) {
        message("PPI analysis disabled.")
        return(NULL)
    }

    summary_df <- de_res$summary_df
    if (is.null(summary_df)) {
        message("No summary_df available. Skipping PPI analysis.")
        return(NULL)
    }

    # Identify contrasts from summary columns
    contrasts <- sub("^padj\\.imputs\\.", "",
                     grep("^padj\\.imputs\\.", colnames(summary_df), value = TRUE))

    if (length(contrasts) == 0) {
        message("No contrasts found. Skipping PPI analysis.")
        return(NULL)
    }

    # Build per-contrast DE tables with gene symbols
    de_tables <- lapply(setNames(contrasts, contrasts), function(cn) {
        tbl <- extract_de_table_for_pathway(summary_df, cn, config)
        map_proteins_to_gene_symbols(tbl, summary_df, config)
    })

    # Derive species ID from annotation.organism (single source of truth)
    ppi_cfg <- cfg$ppi %||% list()
    if (is.null(ppi_cfg$species)) {
        organism <- cfg$annotation$organism %||% "Homo sapiens"
        organism <- normalize_organism_name(organism)
        org_info <- get_organism_info(organism)
        if (!is.na(org_info$taxid)) {
            config$modes$proteomics$ppi$species <- org_info$taxid
            message("PPI: Species ", org_info$taxid, " (", organism, ")")
        } else {
            message("PPI: Unknown organism '", organism, "'. ",
                    "Set ppi.species in config to the NCBI taxonomy ID.")
            return(NULL)
        }
    }

    run_ppi_network_analysis(de_tables, pre$meta, config, out_dir)
}
