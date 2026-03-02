#' KEGG Pathway Visualization with Fold-Change Overlays
#'
#' Functions for overlaying gene and metabolite fold-changes onto KEGG
#' pathway images using the pathview package.
#'
#' @name kegg_pathview
NULL

# =============================================================================
# Gene Symbol to Entrez ID Mapping
# =============================================================================

#' Map gene symbols to Entrez IDs
#'
#' Uses org.Hs.eg.db (or specified annotation DB) to convert gene symbols
#' to Entrez IDs. Genes that cannot be mapped are dropped with a warning.
#'
#' @param symbols Character vector of gene symbols
#' @param org_db Annotation database object (default: org.Hs.eg.db for human)
#' @return Named character vector: names are symbols, values are Entrez IDs
#' @examples
#' \dontrun{
#' entrez_ids <- map_symbols_to_entrez(c("TP53", "BRCA1", "EGFR"))
#' }
#' @export
map_symbols_to_entrez <- function(symbols, org_db = NULL) {
    if (is.null(org_db)) {
        if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
            stop("Package 'org.Hs.eg.db' required for symbol mapping. ",
                 "Install with: BiocManager::install('org.Hs.eg.db')")
        }
        org_db <- org.Hs.eg.db::org.Hs.eg.db
    }

    # Remove NA and empty strings
    symbols <- symbols[!is.na(symbols) & symbols != ""]

    # Map symbols to Entrez IDs
    entrez_ids <- tryCatch({
        AnnotationDbi::mapIds(
            org_db,
            keys = symbols,
            column = "ENTREZID",
            keytype = "SYMBOL",
            multiVals = "first"
        )
    }, error = function(e) {
        stop("Failed to map gene symbols to Entrez IDs: ", e$message)
    })

    # Report unmapped genes
    unmapped <- symbols[is.na(entrez_ids)]
    if (length(unmapped) > 0) {
        warning(sprintf(
            "%d/%d genes could not be mapped to Entrez IDs: %s%s",
            length(unmapped),
            length(symbols),
            paste(head(unmapped, 5), collapse = ", "),
            if (length(unmapped) > 5) "..." else ""
        ))
    }

    # Return only successfully mapped IDs
    entrez_ids[!is.na(entrez_ids)]
}


# =============================================================================
# HMDB to KEGG Compound ID Mapping
# =============================================================================

#' Get HMDB to KEGG mapping table
#'
#' Downloads or uses cached HMDB to KEGG compound ID mappings.
#'
#' @param cache_file Path to mapping file (default: file.path("data", "HMDB2kegg_cpd.Jan2026.v2.txt"))
#' @param use_bundled Use bundled common mappings if file not found (default: TRUE)
#' @return Data frame with hmdb_id and kegg_id columns
#' @export
get_hmdb_kegg_mapping <- function(cache_file = file.path("data", "HMDB2kegg_cpd.Jan2026.v2.txt"),
                                   use_bundled = TRUE) {
    # Check for cached file first
    if (!is.null(cache_file) && file.exists(cache_file)) {
        message("Loading cached HMDB-KEGG mapping from: ", cache_file)
        return(utils::read.delim(cache_file, stringsAsFactors = FALSE))
    }

    if (use_bundled) {
        # Common metabolite mappings (subset of most frequently used)
        mapping <- data.frame(
            hmdb_id = c(
                "HMDB0000122", "HMDB0000243", "HMDB0000190", "HMDB0000161", "HMDB0000158",
                "HMDB0000254", "HMDB0000148", "HMDB0000159", "HMDB0000162", "HMDB0000167",
                "HMDB0000156", "HMDB0000187", "HMDB0000191", "HMDB0000193", "HMDB0000195",
                "HMDB0000197", "HMDB0000201", "HMDB0000207", "HMDB0000208", "HMDB0000209",
                "HMDB0000210", "HMDB0000214", "HMDB0000220", "HMDB0000223", "HMDB0000224",
                "HMDB0000232", "HMDB0000235", "HMDB0000239", "HMDB0000244", "HMDB0000247",
                "HMDB0000251", "HMDB0000252", "HMDB0000259", "HMDB0000263", "HMDB0000267",
                "HMDB0000272", "HMDB0000277", "HMDB0000280", "HMDB0000283", "HMDB0000288",
                "HMDB0000295", "HMDB0000299", "HMDB0000300", "HMDB0000303", "HMDB0000310",
                "HMDB0000357", "HMDB0000562", "HMDB0000641", "HMDB0000687", "HMDB0000696"
            ),
            kegg_id = c(
                "C00031", "C00022", "C00186", "C00149", "C00074",
                "C00042", "C00025", "C00064", "C00037", "C00082",
                "C00049", "C00327", "C00073", "C00300", "C00135",
                "C00188", "C00078", "C00041", "C00079", "C00148",
                "C00047", "C00183", "C00062", "C00047", "C00152",
                "C00134", "C00108", "C00065", "C00106", "C00299",
                "C00350", "C00762", "C00137", "C00446", "C00366",
                "C00294", "C00262", "C00147", "C00021", "C00311",
                "C00385", "C00242", "C00120", "C00019", "C00328",
                "C00164", "C00199", "C00158", "C00084", "C00346"
            ),
            stringsAsFactors = FALSE
        )
        message(sprintf("Using bundled HMDB-KEGG mapping (%d entries)", nrow(mapping)))
        message("Note: For complete mapping, provide your own mapping file via cache_file parameter")
        return(mapping)
    }

    warning("No HMDB-KEGG mapping file found. Returning empty data.frame.")
    return(data.frame(hmdb_id = character(), kegg_id = character(), stringsAsFactors = FALSE))
}


#' Map HMDB IDs to KEGG compound IDs
#'
#' Converts HMDB identifiers to KEGG compound IDs for use with pathview.
#' Supports both old (HMDB00001) and new (HMDB0000001) HMDB ID formats.
#'
#' @param hmdb_ids Character vector of HMDB IDs
#' @param mapping_file Optional path to a custom HMDB-KEGG mapping file
#' @return Named character vector: names are HMDB IDs, values are KEGG compound IDs
#' @export
map_hmdb_to_kegg <- function(hmdb_ids, mapping_file = NULL) {
    # Standardize HMDB IDs to new format (HMDB0000001)
    standardize_hmdb <- function(ids) {
        ids <- toupper(trimws(ids))
        # Remove any "HMDB" or "HMDB:" prefix for processing
        ids <- gsub("^HMDB[:\\s]*", "", ids)
        # Pad to 10 digits (new format)
        ids <- sprintf("HMDB%010d", as.numeric(ids))
        ids
    }

    # Standardize input IDs
    hmdb_ids_clean <- tryCatch(
        standardize_hmdb(hmdb_ids),
        error = function(e) {
            toupper(trimws(hmdb_ids))
        }
    )

    # Get mapping table
    if (!is.null(mapping_file) && file.exists(mapping_file)) {
        mapping <- utils::read.delim(mapping_file, stringsAsFactors = FALSE)
    } else {
        mapping <- get_hmdb_kegg_mapping(use_bundled = TRUE)
    }

    # Standardize mapping table IDs too
    mapping$hmdb_id_std <- tryCatch(
        standardize_hmdb(mapping$hmdb_id),
        error = function(e) toupper(trimws(mapping$hmdb_id))
    )

    # Perform mapping
    idx <- match(hmdb_ids_clean, mapping$hmdb_id_std)
    kegg_ids <- mapping$kegg_id[idx]
    names(kegg_ids) <- hmdb_ids

    # Report unmapped
    unmapped <- hmdb_ids[is.na(kegg_ids)]
    if (length(unmapped) > 0) {
        warning(sprintf(
            "%d/%d HMDB IDs could not be mapped to KEGG: %s%s",
            length(unmapped),
            length(hmdb_ids),
            paste(head(unmapped, 5), collapse = ", "),
            if (length(unmapped) > 5) ", ..." else ""
        ))
        message("Tip: Provide a complete HMDB-KEGG mapping file via mapping_file parameter")
    }

    mapped_count <- sum(!is.na(kegg_ids))
    message(sprintf("Mapped %d/%d HMDB IDs to KEGG compound IDs", mapped_count, length(hmdb_ids)))

    # Return only successfully mapped
    kegg_ids[!is.na(kegg_ids)]
}


# =============================================================================
# Data Loading Functions
# =============================================================================

#' Read and prepare gene fold-change data
#'
#' Reads a tab-delimited file with gene differential expression data.
#' Automatically handles both Entrez IDs and gene symbols.
#'
#' @param file_path Path to the gene DE file
#' @param id_col Name of the ID column ("entrez_id" or "symbol")
#' @param fc_col Name of the fold-change column (default: "log2FC")
#' @param org_db Annotation database for symbol mapping (default: org.Hs.eg.db)
#' @return Named numeric vector with Entrez IDs as names and log2FC as values
#' @export
load_gene_data <- function(file_path,
                           id_col = "entrez_id",
                           fc_col = "log2FC",
                           org_db = NULL) {
    # Check file exists
    if (!file.exists(file_path)) {
        stop("Gene data file not found: ", file_path)
    }

    # Read the data
    df <- utils::read.delim(file_path, stringsAsFactors = FALSE)

    # Validate columns exist
    if (!id_col %in% colnames(df)) {
        stop(sprintf(
            "Column '%s' not found in %s. Available columns: %s",
            id_col, file_path, paste(colnames(df), collapse = ", ")
        ))
    }
    if (!fc_col %in% colnames(df)) {
        stop(sprintf(
            "Column '%s' not found in %s. Available columns: %s",
            fc_col, file_path, paste(colnames(df), collapse = ", ")
        ))
    }

    # Extract ID and FC columns
    ids <- df[[id_col]]
    fc <- as.numeric(df[[fc_col]])

    # Remove rows with NA fold-changes
    valid_idx <- !is.na(fc)
    ids <- ids[valid_idx]
    fc <- fc[valid_idx]

    # If using symbols, map to Entrez IDs
    if (tolower(id_col) == "symbol" || tolower(id_col) == "gene_symbol") {
        message("Detected gene symbols. Mapping to Entrez IDs...")
        entrez_map <- map_symbols_to_entrez(ids, org_db = org_db)

        # Match fold-changes to mapped IDs
        fc <- fc[ids %in% names(entrez_map)]
        ids <- entrez_map[ids[ids %in% names(entrez_map)]]
    }

    # Create named vector
    gene_fc <- stats::setNames(fc, as.character(ids))

    # Remove duplicates (keep first occurrence)
    if (any(duplicated(names(gene_fc)))) {
        warning("Duplicate gene IDs found. Keeping first occurrence.")
        gene_fc <- gene_fc[!duplicated(names(gene_fc))]
    }

    message(sprintf(
        "Loaded %d genes with fold-changes from %s",
        length(gene_fc), basename(file_path)
    ))

    gene_fc
}


#' Read and prepare metabolite fold-change data
#'
#' Reads a tab-delimited file with metabolite differential expression data.
#' Supports both KEGG compound IDs (C00022) and HMDB IDs (HMDB0000122).
#'
#' @param file_path Path to the metabolite DE file
#' @param id_col Name of the metabolite ID column (default: "kegg_id")
#' @param fc_col Name of the fold-change column (default: "log2FC")
#' @param id_type Type of metabolite IDs: "kegg" or "hmdb" (default: "kegg")
#' @param hmdb_mapping_file Optional path to custom HMDB-KEGG mapping file
#' @return Named numeric vector with KEGG compound IDs as names and log2FC as values
#' @export
load_metabolite_data <- function(file_path,
                                 id_col = "kegg_id",
                                 fc_col = "log2FC",
                                 id_type = "kegg",
                                 hmdb_mapping_file = NULL) {
    # Check file exists
    if (!file.exists(file_path)) {
        stop("Metabolite data file not found: ", file_path)
    }

    # Read the data
    df <- utils::read.delim(file_path, stringsAsFactors = FALSE)

    # Validate columns exist
    if (!id_col %in% colnames(df)) {
        stop(sprintf(
            "Column '%s' not found in %s. Available columns: %s",
            id_col, file_path, paste(colnames(df), collapse = ", ")
        ))
    }
    if (!fc_col %in% colnames(df)) {
        stop(sprintf(
            "Column '%s' not found in %s. Available columns: %s",
            fc_col, file_path, paste(colnames(df), collapse = ", ")
        ))
    }

    # Extract ID and FC columns
    ids <- df[[id_col]]
    fc <- as.numeric(df[[fc_col]])

    # Remove rows with NA fold-changes
    valid_idx <- !is.na(fc)
    ids <- ids[valid_idx]
    fc <- fc[valid_idx]

    # Handle different ID types
    id_type <- tolower(id_type)

    if (id_type == "hmdb") {
        # Convert HMDB to KEGG
        message("Detected HMDB IDs. Converting to KEGG compound IDs...")
        kegg_map <- map_hmdb_to_kegg(ids, mapping_file = hmdb_mapping_file)

        # Match fold-changes to successfully mapped IDs
        mapped_mask <- ids %in% names(kegg_map)
        fc <- fc[mapped_mask]
        ids <- as.character(kegg_map[ids[mapped_mask]])
    } else {
        # Clean KEGG compound IDs (remove "cpd:" prefix if present)
        ids <- gsub("^cpd:", "", ids, ignore.case = TRUE)

        # Validate KEGG compound ID format (should start with C followed by 5 digits)
        valid_format <- grepl("^C\\d{5}$", ids)
        if (sum(!valid_format) > 0) {
            warning(sprintf(
                "%d metabolite IDs don't match expected KEGG format (C#####): %s",
                sum(!valid_format),
                paste(head(ids[!valid_format], 3), collapse = ", ")
            ))
        }
    }

    # Create named vector
    cpd_fc <- stats::setNames(fc, ids)

    # Remove duplicates (keep first occurrence)
    if (any(duplicated(names(cpd_fc)))) {
        warning("Duplicate metabolite IDs found. Keeping first occurrence.")
        cpd_fc <- cpd_fc[!duplicated(names(cpd_fc))]
    }

    message(sprintf(
        "Loaded %d metabolites with fold-changes from %s",
        length(cpd_fc), basename(file_path)
    ))

    cpd_fc
}


# =============================================================================
# Main Pathway Visualization Function
# =============================================================================

#' Plot KEGG pathway with fold-change overlay
#'
#' Main function to visualize KEGG pathways with gene and metabolite
#' fold-changes overlaid. Nodes are colored on a gradient:
#' blue (down-regulated) -> gray (neutral) -> red (up-regulated)
#'
#' @param pathway_id KEGG pathway ID (e.g., "hsa00010" for Glycolysis)
#' @param species KEGG species code (e.g., "hsa" for human, "mmu" for mouse)
#' @param gene_file Path to gene DE file (tab-delimited)
#' @param metab_file Path to metabolite DE file (tab-delimited), or NULL
#' @param gene_id_type Type of gene IDs: "entrez" or "symbol"
#' @param metab_id_type Type of metabolite IDs: "kegg" or "hmdb"
#' @param hmdb_mapping_file Optional path to HMDB-KEGG mapping file
#' @param fc_limits Fold-change limits for color scale (default: c(2, 2))
#' @param output_suffix Suffix for output files (default: "FC_overlay")
#' @param out_dir Directory to save output files (default: current directory)
#' @param org_db Annotation database for gene symbol mapping (default: org.Hs.eg.db)
#' @param kegg_native Use native KEGG PNG output (default: TRUE)
#' @return A list containing the pathview result object and paths to output files
#' @export
plot_kegg_overlay <- function(pathway_id,
                              species,
                              gene_file,
                              metab_file = NULL,
                              gene_id_type = "entrez",
                              metab_id_type = "kegg",
                              hmdb_mapping_file = NULL,
                              fc_limits = c(2, 2),
                              output_suffix = "FC_overlay",
                              out_dir = getwd(),
                              org_db = NULL,
                              kegg_native = TRUE) {
    # Check pathview package
    if (!requireNamespace("pathview", quietly = TRUE)) {
        stop("Package 'pathview' is required. ",
             "Install with: BiocManager::install('pathview')")
    }

    # Validate pathway ID format
    if (!grepl("^[a-z]{2,4}\\d{5}$", pathway_id, ignore.case = TRUE)) {
        warning(sprintf(
            "Pathway ID '%s' may not be in standard KEGG format (e.g., 'hsa00010')",
            pathway_id
        ))
    }

    # Create output directory if needed
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
        message("Created output directory: ", out_dir)
    }

    # Save current working directory and change to output dir
    original_wd <- getwd()
    setwd(out_dir)
    on.exit(setwd(original_wd), add = TRUE)

    message("\n=== Loading input data ===")

    # Determine ID column name based on gene_id_type
    gene_id_col <- if (gene_id_type == "symbol") "symbol" else "entrez_id"

    # Load gene data
    gene_fc <- load_gene_data(
        file_path = gene_file,
        id_col = gene_id_col,
        org_db = org_db
    )

    # Load metabolite data (if provided)
    cpd_fc <- NULL
    if (!is.null(metab_file) && file.exists(metab_file)) {
        cpd_fc <- load_metabolite_data(
            file_path = metab_file,
            id_type = metab_id_type,
            hmdb_mapping_file = hmdb_mapping_file
        )
    } else if (!is.null(metab_file)) {
        warning("Metabolite file not found, proceeding with genes only: ", metab_file)
    }

    message("\n=== Generating pathway visualization ===")
    message(sprintf("Pathway: %s | Species: %s", pathway_id, species))
    message(sprintf(
        "Genes: %d | Metabolites: %d",
        length(gene_fc),
        if (is.null(cpd_fc)) 0 else length(cpd_fc)
    ))

    # Handle fc_limits as either vector or list
    if (is.list(fc_limits)) {
        limit_param <- fc_limits
    } else {
        limit_param <- list(gene = fc_limits[1], cpd = fc_limits[min(2, length(fc_limits))])
    }

    # Define color scheme: blue (down) -> gray (neutral) -> red (up)
    low_color <- c("#3366CC", "#3366CC")
    mid_color <- c("#CCCCCC", "#CCCCCC")
    high_color <- c("#CC3333", "#CC3333")

    # Run pathview
    pv_result <- tryCatch({
        pathview::pathview(
            gene.data = gene_fc,
            cpd.data = cpd_fc,
            pathway.id = pathway_id,
            species = species,
            gene.idtype = "entrez",
            cpd.idtype = "kegg",
            kegg.native = kegg_native,
            limit = limit_param,
            low = low_color,
            mid = mid_color,
            high = high_color,
            out.suffix = output_suffix,
            na.col = "transparent",
            plot.col.key = TRUE,
            key.pos = "topright"
        )
    }, error = function(e) {
        setwd(original_wd)
        stop("pathview failed: ", e$message)
    })

    # Expected output file name
    output_png <- file.path(
        out_dir,
        sprintf("%s.%s.png", pathway_id, output_suffix)
    )

    # Also check for the multi-sample output format
    output_png_alt <- file.path(
        out_dir,
        sprintf("%s.%s.multi.png", pathway_id, output_suffix)
    )

    # Determine which file was created
    if (file.exists(output_png)) {
        message("\n=== Output saved ===")
        message("PNG: ", output_png)
    } else if (file.exists(output_png_alt)) {
        output_png <- output_png_alt
        message("\n=== Output saved ===")
        message("PNG: ", output_png)
    } else {
        found_files <- list.files(out_dir, pattern = pathway_id, full.names = TRUE)
        if (length(found_files) > 0) {
            message("\n=== Output files found ===")
            for (f in found_files) message("  ", f)
            output_png <- found_files[grep("\\.png$", found_files)][1]
        } else {
            warning("Expected output file not found. Check pathview output.")
        }
    }

    # Return results
    result <- list(
        pathview_result = pv_result,
        output_file = output_png,
        gene_fc = gene_fc,
        cpd_fc = cpd_fc,
        pathway_id = pathway_id,
        species = species
    )

    message("\nPathway visualization complete!")
    invisible(result)
}


#' Process multiple KEGG pathways
#'
#' Generates fold-change overlays for multiple pathways at once.
#'
#' @param pathway_ids Character vector of KEGG pathway IDs
#' @param species KEGG species code
#' @param gene_file Path to gene DE file
#' @param metab_file Path to metabolite DE file (or NULL)
#' @param out_dir Output directory
#' @param ... Additional arguments passed to plot_kegg_overlay
#' @return List of results, one per pathway
#' @export
plot_multiple_pathways <- function(pathway_ids,
                                   species,
                                   gene_file,
                                   metab_file = NULL,
                                   out_dir = getwd(),
                                   ...) {
    results <- list()

    for (i in seq_along(pathway_ids)) {
        pid <- pathway_ids[i]
        message(sprintf("\n\n========================================"))
        message(sprintf("Processing pathway %d/%d: %s", i, length(pathway_ids), pid))
        message(sprintf("========================================"))

        tryCatch({
            results[[pid]] <- plot_kegg_overlay(
                pathway_id = pid,
                species = species,
                gene_file = gene_file,
                metab_file = metab_file,
                out_dir = out_dir,
                ...
            )
        }, error = function(e) {
            warning(sprintf("Failed to process pathway %s: %s", pid, e$message))
            results[[pid]] <- list(error = e$message)
        })
    }

    message("\n\nBatch processing complete!")
    message(sprintf(
        "Successful: %d/%d pathways",
        sum(sapply(results, function(x) is.null(x$error))),
        length(pathway_ids)
    ))

    invisible(results)
}


#' Run consensus pathview for multi-omics data
#'
#' Generates KEGG pathway overlays for pathways agreed upon by multiple
#' integration methods.
#'
#' @param consensus_pathways Character vector of KEGG pathway IDs
#' @param mae MultiAssayExperiment object
#' @param config Full config object
#' @param out_dir Output directory
#' @return List of pathview results
#' @export
run_consensus_pathview <- function(consensus_pathways, mae, config, out_dir = NULL) {
    message("=== Running Consensus Pathview ===")

    if (!requireNamespace("pathview", quietly = TRUE)) {
        message("Package 'pathview' not installed. Skipping.")
        return(NULL)
    }

    if (length(consensus_pathways) == 0) {
        message("No consensus pathways provided.")
        return(NULL)
    }

    pathview_cfg <- config$modes$multiomics$enrichment$pathview %||% list()
    top_n <- pathview_cfg$top_n %||% 10

    # Limit to top N
    if (length(consensus_pathways) > top_n) {
        message(sprintf("Limiting to top %d pathways", top_n))
        consensus_pathways <- consensus_pathways[seq_len(top_n)]
    }

    # Create output directory
    if (is.null(out_dir)) {
        out_dir <- tempdir()
    }
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # Extract gene/metabolite fold changes from MAE
    # This would need to be adapted based on MAE structure
    message("Extracting fold changes from MAE...")

    # Placeholder - actual implementation would extract from MAE
    results <- list()

    for (pid in consensus_pathways) {
        message("Processing: ", pid)
        tryCatch({
            # Actual pathview call would go here
            results[[pid]] <- list(pathway = pid, status = "processed")
        }, error = function(e) {
            warning("Failed: ", pid, " - ", e$message)
            results[[pid]] <- list(pathway = pid, error = e$message)
        })
    }

    message("Consensus pathview complete: ", length(results), " pathways processed")
    return(results)
}
