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

    # If a custom annotation mapping exists, remap IDs that the standard
    # gene-symbol mapping could not resolve (e.g. KAE* -> GL50803_* for
    # Giardia).  This keeps the pipeline annotation untouched — the
    # conversion is only applied to the IDs sent to STRING.
    custom_file <- cfg$annotation$custom_mapping_file
    if (!is.null(custom_file) && file.exists(custom_file)) {
        custom_map <- tryCatch(read.csv(custom_file, stringsAsFactors = FALSE),
                               error = function(e) NULL)
        if (!is.null(custom_map) &&
            "protein_id" %in% colnames(custom_map) &&
            "gene_name"  %in% colnames(custom_map)) {
            id_lookup <- setNames(custom_map$gene_name, custom_map$protein_id)
            # Also map versioned IDs (e.g. KAE8301209.1 -> KAE8301209)
            de_tables <- lapply(de_tables, function(tbl) {
                fids <- tbl$FeatureID
                # For protein groups (A;B;C), use first member for lookup
                fids_first <- sub(";.*$", "", fids)
                # Try exact match first, then without version suffix
                mapped <- id_lookup[fids_first]
                unversioned <- sub("\\.[0-9]+$", "", fids_first)
                still_na <- is.na(mapped)
                mapped[still_na] <- id_lookup[unversioned[still_na]]
                tbl$FeatureID <- ifelse(is.na(mapped), fids, mapped)
                # Deduplicate after remapping (keep lowest padj)
                tbl <- tbl[order(tbl$padj, na.last = TRUE), ]
                tbl[!duplicated(tbl$FeatureID), ]
            })
            n_remapped <- sum(sapply(de_tables, function(t)
                sum(grepl("^GL50803_|^ENS|^[A-Z][a-z]", t$FeatureID))))
            message("PPI: Remapped ", n_remapped,
                    " protein IDs via custom annotation for STRING query")
        }
    }

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
