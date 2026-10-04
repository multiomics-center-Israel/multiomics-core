#' Build gene-protein ID mapping for multi-omics harmonization
#'
#' Creates a mapping table between gene IDs (from RNA-seq) and protein IDs
#' (from proteomics) using multiple strategies:
#' 1. Gene symbols (if both have symbols)
#' 2. ENSEMBL/UniProt cross-reference via biomaRt
#' 3. Custom mapping file (if provided)
#'
#' @param rna_data Preprocessed RNA-seq data with row_data
#' @param prot_data Preprocessed proteomics data with row_data
#' @param config Full config object
#' @return Data frame with columns: gene_id, protein_id, mapping_source
build_gene_protein_mapping <- function(rna_data, prot_data, config) {

    # Extract gene and protein IDs
    gene_ids <- rownames(rna_data$expr_work)
    protein_ids <- rownames(prot_data$expr_work)

    # Check for custom mapping file (resolved from config)
    custom_map_file <- resolve_gene_protein_mapping_file(config)
    if (!is.null(custom_map_file) && file.exists(custom_map_file)) {
        message("Using gene-protein mapping file: ", basename(custom_map_file))
        return(load_custom_gene_protein_mapping(custom_map_file, gene_ids, protein_ids))
    }

    # Fallback: symbol-based mapping
    message("Building gene-protein mapping via symbols...")

    # Extract symbols from row_data if available
    gene_symbols <- extract_gene_symbols(rna_data$row_data, config$modes$rna)
    protein_symbols <- extract_protein_symbols(prot_data$row_data, config$modes$proteomics)

    # Map by matching symbols
    mapping_df <- map_by_symbols(gene_ids, protein_ids, gene_symbols, protein_symbols)

    message(sprintf(
        "Gene-protein mapping complete: %d pairs (%d genes, %d proteins)",
        nrow(mapping_df),
        length(unique(mapping_df$gene_id)),
        length(unique(mapping_df$protein_id))
    ))

    mapping_df
}


#' Resolve the custom gene-protein mapping file path from config
#'
#' Looks up the mapping file first under the multiomics mode, then the global
#' setting (resolved relative to the raw data directory). Returns NULL when no
#' path is configured.
#'
#' @param config Full config object.
#' @return Character path to the mapping file, or NULL if none is configured.
resolve_gene_protein_mapping_file <- function(config) {
    custom_map_file <- config$modes$multiomics$gene_protein_mapping_file
    if (!is.null(custom_map_file) && file.exists(custom_map_file)) {
        return(custom_map_file)
    }
    gpm <- config$global$gene_protein_mapping
    if (!is.null(gpm) && nzchar(gpm)) {
        return(file.path(config$project$dir, config$paths$raw, gpm))
    }
    NULL
}


#' Build a gene-protein mapping from scratch for a given set of IDs
#'
#' Reads the configured custom mapping file directly and filters it to the
#' supplied gene and protein IDs. Deliberately independent of the harmonized
#' MultiAssayExperiment, so downstream steps (e.g. concordance) can rebuild the
#' mapping from the original ID space rather than trusting harmonized
#' (\code{GENE_*}) feature IDs.
#'
#' When \code{modes$multiomics$require_one_to_one_mapping} is TRUE, the same
#' 1:1 filter that harmonization applies is applied here, at the same scope.
#' The order is: narrow to the experiment's measured IDs, judge ambiguity, then
#' narrow to the supplied (DE) IDs. Both ends of that order matter:
#' \itemize{
#'   \item Judging ambiguity on the whole file would drop a valid pair whenever
#'     a reusable mapping file lists an isoform that this experiment never
#'     measured — there is no ambiguity to resolve if only one side was observed.
#'   \item Judging it after narrowing to the DE IDs would call a gene
#'     unambiguous just because one of its proteins missed the DE cutoffs,
#'     making "1:1" drift with the cutoffs.
#' }
#'
#' @param gene_ids Character vector of RNA-seq gene IDs, in their original space.
#' @param protein_ids Character vector of proteomics protein IDs, original space.
#' @param config Full config object.
#' @param scope_gene_ids Gene IDs measured in this experiment (the RNA
#'   expression matrix rownames), defining the scope at which ambiguity is
#'   judged. NULL falls back to the whole mapping file.
#' @param scope_protein_ids Protein IDs measured in this experiment, likewise.
#' @return Data frame (gene_id, protein_id, mapping_source[, gene_symbol]), or
#'   NULL if no custom mapping file is configured or it cannot be read.
build_gene_protein_mapping_from_ids <- function(gene_ids, protein_ids, config,
                                                scope_gene_ids = NULL,
                                                scope_protein_ids = NULL) {
    custom_map_file <- resolve_gene_protein_mapping_file(config)
    if (is.null(custom_map_file) || !file.exists(custom_map_file)) {
        return(NULL)
    }
    mapping <- tryCatch(
        load_custom_gene_protein_mapping(custom_map_file),
        error = function(e) {
            warning("Could not read gene-protein mapping file: ", e$message)
            NULL
        }
    )
    if (is.null(mapping)) return(NULL)

    # Same scope harmonization uses (build_gene_protein_mapping() narrows to the
    # expression IDs, then 01_mod_harmonization.R applies the filter).
    mapping <- narrow_mapping_to_ids(
        mapping,
        gene_ids = scope_gene_ids,
        protein_ids = scope_protein_ids
    )

    # This mapping replaces the harmonized one for concordance, so it must honour
    # the same setting: otherwise a gene mapped to several proteins enters the
    # concordance table once per protein and gets extra weight in the reported r.
    if (isTRUE(config$modes$multiomics$require_one_to_one_mapping)) {
        mapping <- filter_to_one_to_one_mapping(mapping)
    }

    narrow_mapping_to_ids(mapping, gene_ids = gene_ids, protein_ids = protein_ids)
}


#' Narrow a gene-protein mapping to a given ID space
#'
#' @param mapping_df Mapping data frame with gene_id / protein_id columns.
#' @param gene_ids Gene IDs to keep, or NULL to keep every gene.
#' @param protein_ids Protein IDs to keep, or NULL to keep every protein.
#' @return The mapping restricted to rows whose gene AND protein are in scope.
narrow_mapping_to_ids <- function(mapping_df, gene_ids = NULL, protein_ids = NULL) {
    if (is.null(mapping_df) || nrow(mapping_df) == 0) return(mapping_df)

    keep <- rep(TRUE, nrow(mapping_df))
    if (!is.null(gene_ids))    keep <- keep & mapping_df$gene_id %in% gene_ids
    if (!is.null(protein_ids)) keep <- keep & mapping_df$protein_id %in% protein_ids

    mapping_df[keep, , drop = FALSE]
}


#' Load custom gene-protein mapping from file
#'
#' Accepts files with gene_id + protein_id columns, or gene_id + uniprot_id.
#' Also stores gene_symbol if present for downstream correlation analysis.
#'
#' @param file_path Path to the mapping CSV.
#' @param gene_ids Gene IDs to narrow to, or NULL to return the whole file.
#' @param protein_ids Protein IDs to narrow to, or NULL to return the whole file.
#' @return Data frame (gene_id, protein_id, mapping_source[, gene_symbol]).
load_custom_gene_protein_mapping <- function(file_path, gene_ids = NULL, protein_ids = NULL) {
    df <- read.csv(file_path, stringsAsFactors = FALSE)

    # Normalize: accept uniprot_id as protein_id alias
    if (!"protein_id" %in% colnames(df) && "uniprot_id" %in% colnames(df)) {
        df$protein_id <- df$uniprot_id
    }

    required_cols <- c("gene_id", "protein_id")
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
        stop(
            "Custom mapping file must contain columns: ",
            paste(required_cols, collapse = ", "),
            " (or uniprot_id instead of protein_id)",
            "\nMissing: ", paste(missing_cols, collapse = ", ")
        )
    }

    # Filter to present IDs
    df <- narrow_mapping_to_ids(df, gene_ids = gene_ids, protein_ids = protein_ids)
    df$mapping_source <- "custom_file"

    keep_cols <- c("gene_id", "protein_id", "mapping_source")
    if ("gene_symbol" %in% colnames(df)) keep_cols <- c(keep_cols, "gene_symbol")

    message(sprintf("  Loaded %d gene-protein pairs from mapping file", nrow(df)))
    df[, keep_cols, drop = FALSE]
}


#' Extract gene symbols from RNA-seq row_data
extract_gene_symbols <- function(row_data, rna_cfg) {
    if (is.null(row_data)) return(NULL)

    # Try common symbol column names
    symbol_cols <- c("symbol", "gene_name", "gene_symbol", "SYMBOL", "Gene_Symbol")
    found_col <- intersect(symbol_cols, colnames(row_data))

    if (length(found_col) > 0) {
        return(as.character(row_data[[found_col[1]]]))
    }

    # Fallback: use gene_id column
    gene_id_col <- rna_cfg$id_columns$gene_id %||% "gene_id"
    if (gene_id_col %in% colnames(row_data)) {
        return(as.character(row_data[[gene_id_col]]))
    }

    NULL
}


#' Extract protein symbols/genes from proteomics row_data
extract_protein_symbols <- function(row_data, prot_cfg) {
    if (is.null(row_data)) return(NULL)

    # Try Genes column (common in DIA-NN/MaxQuant)
    if ("Genes" %in% colnames(row_data)) {
        genes <- as.character(row_data$Genes)
        # Handle multi-gene proteins: take first gene
        genes <- sapply(strsplit(genes, ";"), function(x) trimws(x[1]))
        return(genes)
    }

    # Try other common columns
    symbol_cols <- c("Gene.Names", "gene_name", "Gene", "SYMBOL")
    found_col <- intersect(symbol_cols, colnames(row_data))

    if (length(found_col) > 0) {
        genes <- as.character(row_data[[found_col[1]]])
        genes <- sapply(strsplit(genes, ";"), function(x) trimws(x[1]))
        return(genes)
    }

    NULL
}


#' Map genes to proteins by matching symbols
map_by_symbols <- function(gene_ids, protein_ids, gene_symbols, protein_symbols) {

    if (is.null(gene_symbols) || is.null(protein_symbols)) {
        warning("Cannot map by symbols: missing symbol annotations")
        return(data.frame(
            gene_id = character(0),
            protein_id = character(0),
            mapping_source = character(0),
            stringsAsFactors = FALSE
        ))
    }

    # Normalize symbols (uppercase, trim whitespace)
    gene_symbols_norm <- toupper(trimws(gene_symbols))
    protein_symbols_norm <- toupper(trimws(protein_symbols))

    # Build mapping
    mapping_list <- list()
    for (i in seq_along(gene_ids)) {
        gene_sym <- gene_symbols_norm[i]
        if (is.na(gene_sym) || gene_sym == "") next

        # Find matching proteins
        matching_idx <- which(protein_symbols_norm == gene_sym)
        if (length(matching_idx) > 0) {
            for (j in matching_idx) {
                mapping_list[[length(mapping_list) + 1]] <- list(
                    gene_id = gene_ids[i],
                    protein_id = protein_ids[j],
                    mapping_source = "symbol"
                )
            }
        }
    }

    if (length(mapping_list) == 0) {
        warning("No gene-protein pairs found via symbol matching")
        return(data.frame(
            gene_id = character(0),
            protein_id = character(0),
            mapping_source = character(0),
            stringsAsFactors = FALSE
        ))
    }

    do.call(rbind, lapply(mapping_list, as.data.frame, stringsAsFactors = FALSE))
}


#' Filter mapping to retain only 1:1 matches (optional, for strict harmonization)
#'
#' Removes genes that map to multiple proteins and vice versa.
#' @param mapping_df Mapping data frame from build_gene_protein_mapping()
#' @return Filtered mapping with 1:1 relationships only
filter_to_one_to_one_mapping <- function(mapping_df) {
    # Count occurrences
    gene_counts <- table(mapping_df$gene_id)
    protein_counts <- table(mapping_df$protein_id)

    # Keep only 1:1
    one_to_one <- mapping_df[
        mapping_df$gene_id %in% names(gene_counts[gene_counts == 1]) &
        mapping_df$protein_id %in% names(protein_counts[protein_counts == 1]),
    ]

    removed <- nrow(mapping_df) - nrow(one_to_one)
    if (removed > 0) {
        message(sprintf(
            "Filtered to 1:1 mapping: removed %d ambiguous pairs (%d → %d)",
            removed, nrow(mapping_df), nrow(one_to_one)
        ))
    }

    one_to_one
}
