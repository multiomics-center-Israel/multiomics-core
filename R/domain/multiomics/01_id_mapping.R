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

    # Check for custom mapping file: modes > multiomics, then global
    custom_map_file <- config$modes$multiomics$gene_protein_mapping_file
    if (is.null(custom_map_file) || !file.exists(custom_map_file)) {
        # Fallback: global.gene_protein_mapping (relative to data dir)
        gpm <- config$global$gene_protein_mapping
        if (!is.null(gpm) && nzchar(gpm)) {
            custom_map_file <- file.path(config$project$dir, config$paths$raw, gpm)
        }
    }

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


#' Load custom gene-protein mapping from file
#'
#' Accepts files with gene_id + protein_id columns, or gene_id + uniprot_id.
#' Also stores gene_symbol if present for downstream correlation analysis.
load_custom_gene_protein_mapping <- function(file_path, gene_ids, protein_ids) {
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
    df <- df[df$gene_id %in% gene_ids & df$protein_id %in% protein_ids, ]
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

    # Fallback: use feature_id column or rownames (payload mode)
    if ("feature_id" %in% colnames(row_data)) {
        return(as.character(row_data$feature_id))
    }

    # Last resort: use rownames as gene symbols
    rn <- rownames(row_data)
    if (!is.null(rn) && length(rn) > 0) return(rn)

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
        # Handle "GL50803_xxx | gene_name" format (pipe-separated)
        genes <- sapply(strsplit(genes, "\\|"), function(x) trimws(x[1]))
        # Remove empty strings
        genes[genes == "" | genes == "NA"] <- NA
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
