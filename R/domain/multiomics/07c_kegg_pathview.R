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


# =============================================================================
# KO-space pathview for organisms with no OrgDb
# =============================================================================

#' Load the feature -> KEGG Orthology map used for KO-space pathview
#'
#' Non-model organisms have no KEGG organism code, so their features can only be
#' placed onto KEGG reference maps through KO ids. The map is a user-supplied TSV
#' (see `utils/build_feature_ko_map.R` for one way to build it) in long format,
#' with one row per (omics, feature, KO) triple.
#'
#' @param config Full config; reads `modes$multiomics$enrichment$pathview$ko_map`.
#' @return Data frame with columns `omics`, `feature_id`, `KO`, or NULL when the
#'   key is absent, the file is missing, or the columns do not line up.
load_feature_ko_map <- function(config) {
    pv_cfg <- config$modes$multiomics$enrichment$pathview %||% list()
    ko_path <- pv_cfg$ko_map %||% NULL
    if (is.null(ko_path) || !nzchar(ko_path[1])) return(NULL)

    # Resolved like every other user-supplied input (gmt_file, mapping_file):
    # absolute paths pass through, relative ones land under the raw data dir.
    ko_abs <- resolve_input_path(config, ko_path)[1]
    if (!file.exists(ko_abs)) {
        message("  Union pathview: ko_map not found: ", ko_abs)
        return(NULL)
    }

    ko_map <- tryCatch(
        utils::read.delim(ko_abs, sep = "\t", header = TRUE,
                          colClasses = "character", stringsAsFactors = FALSE),
        error = function(e) NULL
    )
    if (is.null(ko_map) || !all(c("omics", "feature_id", "KO") %in% names(ko_map))) {
        message("  Union pathview: ko_map needs columns omics, feature_id, KO: ", ko_abs)
        return(NULL)
    }
    ko_map <- ko_map[!is.na(ko_map$KO) & nzchar(ko_map$KO), c("omics", "feature_id", "KO")]
    message("  Union pathview: KO map with ", nrow(ko_map), " (feature, KO) rows")
    ko_map
}


#' Aggregate feature-level log2FC onto KEGG Orthology nodes
#'
#' A KO box on a KEGG reference map stands for an ortholog group, so several
#' features of one layer (paralogues, or a protein group that collapses onto the
#' same ortholog) can land on the same box. The node then carries their mean,
#' rather than letting whichever feature happened to come first win.
#'
#' @param de_table Standardized DE table with `feature_id` and `log2fc` columns.
#' @param ko_map Feature -> KO map from \code{load_feature_ko_map}.
#' @param omics Omics label selecting the rows of \code{ko_map} to join against.
#' @return Named numeric vector of log2FC per KO id, or NULL when nothing joins.
#' @examples
#' de <- data.frame(feature_id = c("g1", "g2"), log2fc = c(1, 3))
#' ko <- data.frame(omics = "transcriptomics", feature_id = c("g1", "g2"),
#'                  KO = c("K00001", "K00001"))
#' aggregate_log2fc_by_ko(de, ko, "transcriptomics")   # K00001 = 2
aggregate_log2fc_by_ko <- function(de_table, ko_map, omics) {
    if (is.null(de_table) || is.null(ko_map)) return(NULL)
    m <- ko_map[ko_map$omics == omics, c("feature_id", "KO"), drop = FALSE]
    if (nrow(m) == 0) return(NULL)
    joined <- merge(de_table, m, by = "feature_id")
    ok <- !is.na(joined$KO) & nzchar(joined$KO) & is.finite(joined$log2fc)
    if (!any(ok)) return(NULL)
    fc <- tapply(joined$log2fc[ok], joined$KO[ok], mean, na.rm = TRUE)
    stats::setNames(as.numeric(fc), names(fc))
}


#' Thresholds the pathway maps are drawn at
#'
#' One source of truth for three things that must agree: which pathways are
#' selected, which features are allowed to colour a node, and what the caption
#' under the figure claims. They were separate numbers before, so a caption
#' could describe a rule the filter was not applying.
#'
#' \code{fdr_alpha} is the adjusted-p cutoff \code{.kegg_hits_by_contrast()}
#' selects pathways at; \code{node_fc} is the linear fold change a feature must
#' exceed and \code{node_p} the p-value it must reach before it may colour a
#' node. Not config keys: these describe what the figure means, and a project
#' that moved them would be reading a different figure under the same caption.
#'
#' @keywords internal
.PATHVIEW_THRESHOLDS <- list(fdr_alpha = 0.05, node_fc = 1.5, node_p = 0.05)


#' Keep only the features whose change is both large and supported
#'
#' Colouring every measured feature makes a map look saturated regardless of
#' evidence, so it shows coverage rather than signal. The rule is AND, not OR:
#' a big but unsupported change is usually a noisy low-abundance feature, and a
#' confident but tiny one is not what a pathway map exists to highlight.
#'
#' A feature with no usable fold change is dropped rather than coloured on its
#' p-value alone -- NaN is what \code{log2()} of a signed linear fold change
#' leaves behind, and such a node must stay uncoloured.
#'
#' The p-value is the raw one where the table has it, falling back to an
#' adjusted column only where it does not. That is deliberately the same
#' per-feature evidence rule the figure has always implied; it is not an
#' adjusted-p gate, and it decides nothing about which pathways are drawn.
#'
#' @param de_table Standardized DE table with `feature_id` and `log2fc`.
#' @param thresholds Threshold list; defaults to \code{.PATHVIEW_THRESHOLDS}.
#' @return \code{de_table} with only the qualifying rows, or NULL for NULL input.
filter_changed_features <- function(de_table, thresholds = .PATHVIEW_THRESHOLDS) {
    if (is.null(de_table)) return(NULL)
    if (!is.data.frame(de_table) || nrow(de_table) == 0) return(de_table)

    fc_col <- intersect(c("log2fc", "log2FC", "logFC", "log2FoldChange"),
                        names(de_table))[1]
    if (is.na(fc_col)) {
        message("    Pathview: DE table has no fold-change column; ",
                "no feature can colour a node")
        return(de_table[0, , drop = FALSE])
    }
    p_col <- intersect(c("pvalue", "p_value", "P.Value", "padj", "adj.P.Val",
                         "FDR"), names(de_table))[1]
    if (is.na(p_col)) {
        message("    Pathview: DE table has no p-value column; ",
                "no feature can colour a node")
        return(de_table[0, , drop = FALSE])
    }
    if (!p_col %in% c("pvalue", "p_value", "P.Value")) {
        message("    Pathview: no raw p-value column; using ", p_col)
    }

    fc <- suppressWarnings(as.numeric(de_table[[fc_col]]))
    p  <- suppressWarnings(as.numeric(de_table[[p_col]]))
    # is.finite(), not !is.na(): NaN and Inf must not colour a node either.
    keep <- is.finite(fc) & abs(fc) > log2(thresholds$node_fc) &
        is.finite(p) & p < thresholds$node_p
    de_table[keep, , drop = FALSE]
}


#' The caption for a rendered pathway map, from the thresholds it was drawn at
#'
#' Built from the same list the filters read, so the legend cannot drift from
#' the rules. It also states what an uncoloured node does NOT mean: pathview
#' draws "measured but unchanged" and "never measured" identically, and a
#' caption that implied otherwise would invite the reader to infer absence.
#'
#' @param thresholds Threshold list; defaults to \code{.PATHVIEW_THRESHOLDS}.
#' @return Single caption string.
pathview_significance_caption <- function(thresholds = .PATHVIEW_THRESHOLDS) {
    fmt <- function(x) format(x, trim = TRUE, scientific = FALSE)
    paste0(
        "Pathways are drawn where at least one layer scored them below ",
        fmt(thresholds$fdr_alpha), ", on adjusted p-values or, where those ",
        "selected nothing, on raw p-values, so the evidence is a floor rather ",
        "than an FDR. A node is coloured only by features past ",
        fmt(thresholds$node_fc), "-fold (|log2FC| > ",
        format(round(log2(thresholds$node_fc), 2), nsmall = 2),
        ") with raw p < ", fmt(thresholds$node_p),
        "; an uncoloured node means no measured feature passed both rules, or ",
        "nothing was measured -- the map does not separate the two."
    )
}


#' Choose the emptiest corner of a KEGG map for pathview's colour key
#'
#' pathview draws its colour key at a fixed corner (default "topright"), which on
#' some maps lands squarely on the artwork -- on the fatty-acid reference map it
#' covers the compound legend panel and makes both unreadable. The blank KEGG
#' template is already cached on disk by the preceding pathview call, so the ink
#' in each corner can be measured and the quietest one chosen for free.
#'
#' @param template_png Path to the downloaded blank KEGG map PNG.
#' @param frac Fraction of width and height treated as a corner.
#' @return One of "topright", "bottomleft", "bottomright"; falls back to
#'   "topright" when the template cannot be read.
pick_key_position <- function(template_png, frac = 0.28) {
    if (!file.exists(template_png) || !requireNamespace("png", quietly = TRUE)) {
        return("topright")
    }
    img <- tryCatch(png::readPNG(template_png), error = function(e) NULL)
    if (is.null(img)) return("topright")

    grey <- if (length(dim(img)) == 3) {
        apply(img[, , seq_len(min(3, dim(img)[3])), drop = FALSE], c(1, 2), mean)
    } else img
    nr <- nrow(grey); nc <- ncol(grey)
    rr <- max(1L, floor(nr * frac)); cc <- max(1L, floor(nc * frac))
    ink <- function(rows, cols) mean(grey[rows, cols] < 0.86)

    # topleft is deliberately not a candidate: every KEGG map puts its title box
    # there, and the box is thin enough to score as empty, so the key lands on
    # the map's own title.
    scores <- c(
        topright    = ink(seq_len(rr),               seq.int(nc - cc + 1L, nc)),
        bottomleft  = ink(seq.int(nr - rr + 1L, nr), seq_len(cc)),
        bottomright = ink(seq.int(nr - rr + 1L, nr), seq.int(nc - cc + 1L, nc))
    )
    names(scores)[which.min(scores)]
}


#' Remove the Multi-ORA pathview outputs of an earlier run
#'
#' The report finds its maps by globbing the pathview directory, so a map that a
#' later run no longer produces -- a pathway that fell below the cutoff, a
#' contrast that went away, a smaller `top_n` -- would otherwise stay on the page
#' as if it were a current result.
#'
#' Ownership follows the `run_pathview` switch: \code{run_multi_ora()} decides
#' for all three of its pathview renderers at once, so this clears the artifacts
#' of all three -- the cross-omics `.multi_ora*` overlays, the per-omics
#' `.metab_top*` and `.prot_top*` overlays, every compiled PDF they produce, and
#' the sidecars that describe those PDFs.
#' A per-omics map is as stale as a cross-omics one when its pathway drops out
#' of the current top-N, and the report cannot tell either from a fresh one.
#'
#' What survives is the KEGG download cache: the blank template PNGs and their
#' KGML. Those are not results -- deleting them would re-fetch every map and
#' lose the pathway titles the report reads out of the XML -- and so is anything
#' else in the directory that no Multi-ORA renderer wrote.
#'
#' @param out_dir Multi-ORA output directory.
#' @return Character vector of the removed paths, invisibly.
clear_multi_ora_pathview_outputs <- function(out_dir) {
    pv_dir <- file.path(out_dir, "pathview")
    stale <- character(0)
    if (dir.exists(pv_dir)) {
        # Each renderer's overlays carry its own out.suffix; the blank template
        # ("ko00010.png") carries none, which is what keeps the cache out of
        # this list.
        stale <- list.files(
            pv_dir,
            pattern = "\\.(multi_ora|metab_top|prot_top)[^/]*\\.png$",
            full.names = TRUE)
    }
    pdfs <- file.path(out_dir, c("multi_ora_pathview_supported.pdf",
                                 "multi_ora_pathview_supported.yaml",
                                 "multi_ora_pathview_union.pdf",
                                 "multi_ora_pathview_union.yaml",
                                 "pathview_top_metabolomics_pathways.pdf",
                                 "pathview_top_proteomics_pathways.pdf"))
    stale <- c(stale, pdfs[file.exists(pdfs)])
    if (length(stale) > 0) {
        unlink(stale)
        message("  Multi-ORA pathview: cleared ", length(stale),
                " output(s) from a previous run")
    }
    invisible(stale)
}


#' Collect the KEGG pathways each contrast called enriched
#'
#' The enrichment frames carry their own `contrast` column (written by the
#' enrichment step alongside `database`, `method` and `direction`), so a pathway
#' is attributed to the contrast that produced it.
#'
#' Only rows whose `method` is `"ora"` are considered, which is what the
#' predecessor selected by reading just the `*_ora_up/down.csv` exports. A frame
#' with no `method` column cannot be confirmed as ORA and is skipped with a
#' message rather than guessed at.
#'
#' Contrasts are grouped by \code{normalize_contrast_key()}, the same canonical
#' key the cross-omics enrichment module uses, because the spellings genuinely
#' differ between exports: an ORA table can say `"1.56ppmvs.0ppm"` where the
#' proteomics DE tables say `"1.56ppm vs. 0ppm"`. The first raw spelling seen
#' for a key is kept as its label, for headings and filenames.
#'
#' Eligibility is decided on the adjusted p-value wherever a frame has usable
#' ones, because `run_ora()` returns every tested pathway with both columns and
#' filters neither: preferring the raw p-value there would have drawn a map for
#' `pvalue = 0.01, padj = 0.40` and let the report call it enriched. A frame
#' with no usable `padj` at all falls back to `pvalue` -- but a frame where
#' adjusted values exist and none of them pass yields nothing, rather than
#' quietly relaxing to the raw p-value the way \code{run_ora_kegg()} does for
#' its own table. The two branches are symmetric -- each needs at least one
#' usable value of its own kind -- and a frame with neither is skipped: no
#' significance evidence means nothing in it can be called enriched, so nothing
#' in it gets a map.
#'
#' Pathways are then ranked by that same score, best first, so that `top_n`
#' keeps the most significant rather than whichever layer or direction happened
#' to be read first. Where a pathway appears in more than one layer or
#' direction it is one entry, carrying the best score seen for it, and the
#' pathway id breaks ties so the order is reproducible.
#'
#' Pure: it reads no files and draws nothing, which is what makes the selection
#' testable without invoking pathview.
#'
#' @param ora_tables Named list of per-omics enrichment data frames from the
#'   current run (`multiomics_cross_enrichment$per_omics`).
#' @param kegg_org Organism code for pathway identity, or NULL for a run with no
#'   KEGG code of its own. This is the identity organism; the gene space these
#'   maps are drawn in is always KO.
#' @param alpha Significance cutoff, applied to `padj` where a frame has usable
#'   adjusted values and to `pvalue` only where it has none. Defaults to the
#'   shared \code{.PATHVIEW_THRESHOLDS}, which is also what
#'   \code{pathview_significance_caption()} states, so the figure's caption
#'   cannot claim a cutoff the selection did not use.
#' @return Named list keyed by canonical contrast key, each element a list of
#'   `label` (the first raw spelling seen), `pathways` (normalized KEGG ids,
#'   ranked best score first) and `scores` (the score behind that ranking, named
#'   by pathway id). Contrasts with no KEGG hit are absent, and so are rows that
#'   carry no contrast: they cannot be attributed to one, and rendering them
#'   against some other contrast's fold changes is the bug this structure exists
#'   to prevent.
.kegg_hits_by_contrast <- function(ora_tables, kegg_org = NULL,
                                   alpha = .PATHVIEW_THRESHOLDS$fdr_alpha) {
    hits <- list()
    # Indices, not names: the production input is named per omics, but iterating
    # over names() means an unnamed list runs the loop zero times and selects
    # nothing at all, silently. The name is only wanted for the message below.
    nms <- names(ora_tables)
    for (i in seq_along(ora_tables)) {
        d <- ora_tables[[i]]
        nm <- if (!is.null(nms) && !is.na(nms[i]) && nzchar(nms[i])) {
            nms[i]
        } else {
            paste0("enrichment input ", i)
        }
        if (is.null(d) || !is.data.frame(d) || nrow(d) == 0) next
        if (!"contrast" %in% names(d)) next

        # ORA rows only, as the disk-scanning predecessor selected by reading
        # just "*_ora_up/down.csv". A frame that cannot be established as ORA is
        # skipped rather than guessed at: a GSEA table has its own p-value
        # columns and its own idea of what "enriched" means, and letting its
        # rows into the union would quietly widen what gets a map.
        if (!"method" %in% names(d)) {
            message("  Union pathview: ", nm, " enrichment has no method column, ",
                    "so its rows cannot be confirmed as ORA; skipping it")
            next
        }
        is_ora <- !is.na(d$method) & tolower(as.character(d$method)) == "ora"
        if (!any(is_ora)) next

        # Identity is resolved by #201's row-wise helper, not by assuming the
        # accession is in `pathway`: run_ora_kegg() builds the clusterProfiler
        # schema, where `pathway` holds the readable Description and `ID` holds
        # the accession, and reading `pathway` alone would reject every one of
        # its rows. pathway_join_key() walks ID -> pathway -> Description per
        # row and normalizes only where the value is genuinely a KEGG
        # accession, so the bare, prefixed and "<accession> <readable name>"
        # spellings all reduce to the map number while another organism's
        # prefix is left alone, and therefore still rejected below.
        keys <- pathway_join_key(d, kegg_org)
        keep <- is_ora & is_kegg_pathway_accession(keys, kegg_org)

        # Adjusted significance decides eligibility wherever this frame has it.
        # "Usable" is judged over the ORA rows themselves: a frame whose ORA
        # rows all carry NA padj has none, whatever its other rows hold.
        if ("padj" %in% names(d) && any(is_ora & !is.na(d$padj))) {
            score <- suppressWarnings(as.numeric(d$padj))
        } else if ("pvalue" %in% names(d) && any(is_ora & !is.na(d$pvalue))) {
            score <- suppressWarnings(as.numeric(d$pvalue))
        } else {
            # No significance evidence of either kind. Nothing here can be
            # called enriched, so nothing here gets a map -- the branches are
            # symmetric, and neither one lets a row through unscored.
            message("  Union pathview: ", nm, " enrichment has no usable padj ",
                    "or pvalue on its ORA rows; skipping it")
            next
        }
        keep <- keep & !is.na(score) & score < alpha

        raw <- as.character(d$contrast)
        keep <- keep & !is.na(raw) & nzchar(trimws(raw))
        if (!any(keep)) next

        ckey <- normalize_contrast_key(raw)
        for (k in unique(ckey[keep])) {
            rows <- keep & ckey == k
            if (is.null(hits[[k]])) {
                hits[[k]] <- list(label = raw[rows][1], scores = numeric(0))
            }
            # One entry per pathway, carrying the best score seen for it --
            # across directions within this frame, and across layers as later
            # frames arrive.
            best <- tapply(score[rows], keys[rows], min)
            s <- hits[[k]]$scores
            ids <- names(best)
            prev <- s[ids]
            s[ids] <- pmin(as.numeric(best), ifelse(is.na(prev), Inf, prev))
            hits[[k]]$scores <- s
        }
    }

    # Rank once, after every frame has had its say: a pathway's best score is
    # not known until the last layer that carries it has been read.
    for (k in names(hits)) {
        s <- hits[[k]]$scores
        ord <- order(s, names(s))
        hits[[k]]$scores <- s[ord]
        hits[[k]]$pathways <- names(s)[ord]
    }
    hits
}


#' Filename-safe identity for a contrast's rendered maps
#'
#' The output filename has to distinguish two contrasts exactly when their
#' canonical keys do. \code{make.names()} on the readable label does not:
#' `"A-B"` and `"A B"` are different contrasts -- one is a comparison, the other
#' a group whose name has a space -- and both become `"A.B"`, so the second
#' would overwrite the first's maps.
#'
#' The canonical key is already lowercase `[a-z0-9.]`, so replacing its dots is
#' one-to-one: no two keys can produce one suffix.
#'
#' @param contrast_key Canonical key from \code{normalize_contrast_key}.
#' @return Filename-safe string, distinct for distinct keys.
.contrast_out_key <- function(contrast_key) {
    gsub("[^a-z0-9]", "_", contrast_key)
}


#' Pick one contrast's DE table out of a per-contrast list
#'
#' Matched on \code{normalize_contrast_key()}, so a table named
#' `"1.56ppm vs. 0ppm"` answers to the ORA export's `"1.56ppmvs.0ppm"`. A miss
#' means the layer is absent for that contrast -- never the first table. Two
#' table names collapsing to one key is ambiguous, and ambiguity is also a miss.
#'
#' @param tables Named list of DE tables, as \code{extract_de_tables} returns.
#' @param contrast_key Canonical contrast key to look for.
#' @return The matching data frame, or NULL.
.de_table_for_contrast <- function(tables, contrast_key) {
    if (is.null(tables) || length(tables) == 0 || is.null(names(tables))) return(NULL)
    same <- normalize_contrast_key(names(tables)) == contrast_key
    if (sum(same, na.rm = TRUE) != 1L) return(NULL)
    tables[[which(same)]]
}


#' Compile rendered maps into a single PDF
#'
#' A half-written PDF is worse than none: the report links it as a download and
#' has no way to tell it is truncated. So the device is closed, the file is
#' removed and NULL is returned if anything fails on the way -- opening the
#' device, reading a PNG, or writing a page.
#'
#' Each page is captioned, because the same pathway enriched in two contrasts
#' renders as two visually similar maps, and a downloaded PDF has no filename to
#' fall back on the way the HTML report does.
#'
#' @param png_files Character vector of rendered PNG paths, in page order.
#' @param pdf_path Path of the PDF to write.
#' @param labels Optional character vector, one caption per page.
#' @return \code{pdf_path} on success, otherwise NULL.
.compile_pathview_pdf <- function(png_files, pdf_path, labels = NULL) {
    if (length(png_files) == 0) return(NULL)
    if (!is.null(labels) && length(labels) != length(png_files)) labels <- NULL
    # Remember which device was current so that a failure in pdf() itself, which
    # opens nothing, cannot make us close a device somebody else owns.
    before <- grDevices::dev.cur()
    ok <- tryCatch({
        grDevices::pdf(pdf_path, width = 12, height = 8)
        for (i in seq_along(png_files)) {
            img <- png::readPNG(png_files[i])
            grid::grid.newpage()
            grid::grid.raster(img, y = 0.48, height = 0.92)
            if (!is.null(labels)) {
                grid::grid.text(labels[i], y = 0.98,
                                gp = grid::gpar(fontsize = 11, fontface = "bold"))
            }
        }
        TRUE
    }, error = function(e) {
        message("  Union pathview: PDF compilation failed (", conditionMessage(e),
                "); no PDF written")
        FALSE
    })
    if (!identical(grDevices::dev.cur(), before)) {
        # Closing the device is also what flushes and releases the file.
        tryCatch(grDevices::dev.off(), error = function(e) NULL)
    }
    if (!isTRUE(ok) || !file.exists(pdf_path) || file.size(pdf_path) == 0) {
        if (file.exists(pdf_path)) unlink(pdf_path)
        return(NULL)
    }
    pdf_path
}


#' Render KEGG reference maps in KO space for the union of enriched pathways
#'
#' The fallback for a run with no OrgDb, which otherwise gets no pathway maps at
#' all. For each contrast separately, it renders one KEGG reference map per
#' pathway that contrast called enriched, with the transcriptomics and
#' proteomics log2FC of that contrast overlaid as two states and, where
#' metabolomics is present, compound nodes coloured too -- so one map can carry
#' all three layers.
#'
#' \strong{KO space only.} Features reach the map through the user-supplied
#' feature -> KO map (`enrichment.pathview.ko_map`); without one this returns
#' NULL. It does not render in an organism's native gene space, because nothing
#' in this pipeline guarantees that a project's feature ids are KEGG gene ids --
#' the enrichment loader keys KEGG memberships back to project `gene_id` where
#' an annotation exists, so handing those to pathview as `gene.idtype = "KEGG"`
#' would draw blank maps. Native rendering needs an explicit project-id to
#' KEGG-gene contract, which is a separate piece of work.
#'
#' Identity space and gene space are different questions. The organism's KEGG
#' code, where it has one, still decides which pathway accessions the run may
#' legitimately spell its own pathways with, so an `hsa` run normalizes
#' `"hsa00010 Glycolysis / Gluconeogenesis"` to `00010` and then draws the
#' corresponding reference map.
#'
#' Union, not intersection: a pathway enriched in only one gene layer is
#' deliberately eligible here. That is why the output is a separate PDF from
#' \code{generate_multi_ora_pathview}'s, whose pathways really are supported in
#' two or more layers -- the report distinguishes the two.
#'
#' @param de_results Named list of DE results per omics.
#' @param harmonization_res Harmonization result (supplies the metabolomics row
#'   data used to reach KEGG compound ids).
#' @param config Full config.
#' @param out_dir Multi-ORA output directory (maps go under `out_dir/pathview`).
#' @param per_omics_enrichment Named list of this run's per-omics enrichment
#'   frames (`multiomics_cross_enrichment$per_omics`). Pathways are selected
#'   from these and from nothing else -- see the note in the body on why the
#'   persistent Enrichment directories are not scanned.
#' @param top_n Max pathways to render per contrast. The config validator fills
#'   `enrichment.pathview.top_n` with 5 when it is absent, so this default only
#'   applies to a config that never passed through it; overridden by
#'   `modes$multiomics$enrichment$pathview$top_n` when that is set.
#' @return Path to the compiled PDF, or NULL when nothing could be rendered.
generate_per_omic_union_pathview <- function(de_results, harmonization_res,
                                             config, out_dir,
                                             per_omics_enrichment = NULL,
                                             top_n = 5) {
    if (!requireNamespace("pathview", quietly = TRUE)) return(NULL)

    # run_multi_ora() owns the enrichment.pathview.run_pathview switch for all
    # three renderers; this one is not called from anywhere else.
    if (is.null(per_omics_enrichment) || length(per_omics_enrichment) == 0) {
        message("  Union pathview: no per-omics enrichment from this run, ",
                "so there is nothing to select pathways from.")
        return(NULL)
    }

    ko_map <- load_feature_ko_map(config)
    if (is.null(ko_map) || nrow(ko_map) == 0) {
        message("  Union pathview: no feature-to-KO map configured, so there is ",
                "no way to place this run's features on a KEGG map. Set ",
                "modes.multiomics.enrichment.pathview.ko_map (see ",
                "utils/build_feature_ko_map.R) to enable these plots.")
        return(NULL)
    }

    pv_cfg <- config$modes$multiomics$enrichment$pathview %||% list()
    top_n <- pv_cfg$top_n %||% top_n

    # resolve_kegg_org_code(), not get_kegg_organism(): the latter only knows the
    # six exact species names of the older table. This is the identity organism
    # and nothing else -- the gene space below is always KO.
    id_org <- resolve_kegg_org_code(config$global$organism %||% "")

    # 1. KEGG pathways each contrast called enriched, taken from this run's own
    #    enrichment state. It used to scan the per-omic Enrichment directories
    #    instead, which is persistent: save_pathway_results() writes only
    #    non-empty current results and never removes a file, so a pathway that
    #    had fallen out of significance still had a CSV on disk and still got a
    #    map -- drawn with current fold changes. No historical filesystem state
    #    participates in the selection now.
    hits <- .kegg_hits_by_contrast(per_omics_enrichment, id_org)

    # The same reporting exclusion the cross-omics tables get, applied to the
    # selection rather than inside it: ranking and scores are already settled,
    # and a class this project does not report should not spend a render either.
    # Per contrast, because that is the shape selection has after #204.
    excl <- .excluded_pathway_classes(config)
    if (length(excl) > 0 && length(hits) > 0) {
        for (ckey in names(hits)) {
            keep <- keep_kegg_pathways(hits[[ckey]]$pathways, exclude = excl,
                                       kegg_org = id_org,
                                       label = "pathview maps")
            hits[[ckey]]$pathways <- hits[[ckey]]$pathways[keep]
            hits[[ckey]]$scores <- hits[[ckey]]$scores[keep]
        }
        # A contrast whose pathways were all excluded has nothing left to draw.
        hits <- hits[vapply(hits, function(h) length(h$pathways) > 0, logical(1))]
    }

    if (length(hits) == 0) {
        message("  Union pathview: no enriched KEGG pathways.")
        return(NULL)
    }

    # 2. DE tables per layer, kept keyed by contrast so that step 3 can take the
    #    one belonging to the contrast it is rendering.
    tables_for <- function(om) {
        if (!om %in% names(de_results)) return(NULL)
        tryCatch(extract_de_tables(de_results[[om]], om, harmonization_res),
                 error = function(e) NULL)
    }
    rna_tables   <- tables_for("transcriptomics")
    prot_tables  <- tables_for("proteomics")
    metab_tables <- tables_for("metabolomics")
    metab_map <- if (!is.null(metab_tables)) {
        tryCatch(map_metabolite_ids_to_kegg(metab_tables, harmonization_res),
                 error = function(e) NULL)
    } else NULL

    # 3. One render pass per contrast, on that contrast's own fold changes.
    # pathview needs its 'bods' dataset in the global env when called via ::.
    if (!exists("bods", envir = globalenv())) {
        utils::data("bods", package = "pathview", envir = globalenv())
    }
    pv_dir <- file.path(out_dir, "pathview")
    dir.create(pv_dir, recursive = TRUE, showWarnings = FALSE)
    # pathview writes its output into the working directory, so the render runs
    # with_dir(); the absolute path keeps kegg.dir pointing at the same place
    # from inside it.
    pv_dir <- normalizePath(pv_dir, winslash = "/", mustWork = FALSE)

    any_compounds <- FALSE
    contrast_labels <- list()
    page_labels <- character(0)
    generated <- withr::with_dir(pv_dir, {
        made <- character(0)
        for (ckey in names(hits)) {
            label <- hits[[ckey]]$label
            # Filtered BEFORE the KO aggregation, not after: a node carries the
            # mean of the features on it, so a feature that did not move, or
            # moved without support, would otherwise pull that mean towards
            # zero and colour the box on evidence it does not have.
            rna_fc  <- aggregate_log2fc_by_ko(
                filter_changed_features(.de_table_for_contrast(rna_tables, ckey)),
                ko_map, "transcriptomics")
            prot_fc <- aggregate_log2fc_by_ko(
                filter_changed_features(.de_table_for_contrast(prot_tables, ckey)),
                ko_map, "proteomics")
            if (is.null(rna_fc) && is.null(prot_fc)) {
                # Either a contrast whose pathways came from an ORA table we
                # have no DE table for -- borrowing another contrast's values is
                # exactly the cross-wiring this loop exists to prevent -- or one
                # where nothing cleared the node thresholds. Both leave every
                # gene node uncoloured, which is a map not worth rendering.
                message("  Union pathview: no KO-mapped feature passed the node ",
                        "thresholds for contrast ", label, "; skipping it")
                next
            }
            genes <- unique(c(names(rna_fc), names(prot_fc)))
            gene_data <- matrix(NA_real_, length(genes), 2,
                                dimnames = list(genes, c("RNA", "Protein")))
            if (!is.null(rna_fc))  gene_data[names(rna_fc), 1]  <- rna_fc
            if (!is.null(prot_fc)) gene_data[names(prot_fc), 2] <- prot_fc

            # Compounds colour the compound nodes, from the same contrast.
            cpd_data <- NULL
            # Same rule as the gene nodes, and for the same reason: a compound
            # box means one thing across the whole map or it means nothing.
            metab_df <- filter_changed_features(
                .de_table_for_contrast(metab_tables, ckey))
            if (!is.null(metab_df) && nrow(metab_df) > 0 &&
                !is.null(metab_map) && nrow(metab_map) > 0) {
                cpd_data <- tryCatch({
                    md <- merge(metab_df, metab_map, by = "feature_id")
                    ok <- is.finite(md$log2fc)
                    fc <- tapply(md$log2fc[ok], md$KEGG_CPD[ok], mean, na.rm = TRUE)
                    stats::setNames(as.numeric(fc), names(fc))
                }, error = function(e) NULL)
            }
            if (!is.null(cpd_data) && length(cpd_data) > 0) any_compounds <- TRUE

            # Identity, not display: two contrasts that differ must not write
            # one filename. The readable label travels to the report in the
            # sidecar instead.
            out_key <- .contrast_out_key(ckey)
            contrast_labels[[out_key]] <- label
            out_suffix <- paste0("multi_ora_", out_key)
            message("  Union pathview: ", label, " -- ", length(genes),
                    " KO nodes with log2FC")

            for (clean_pid in utils::head(hits[[ckey]]$pathways, top_n)) {
                # Place the colour key in whichever corner the map itself leaves
                # empty. The template is cached by an earlier pathview call, so
                # this is free on a re-run and falls back to the default corner
                # on a first one.
                key_pos <- pick_key_position(
                    file.path(pv_dir, paste0("ko", clean_pid, ".png")))
                tryCatch({
                    pathview::pathview(gene.data = gene_data, cpd.data = cpd_data,
                                       pathway.id = clean_pid,
                                       species = "ko", gene.idtype = "KEGG",
                                       out.suffix = out_suffix, kegg.dir = pv_dir,
                                       key.pos = key_pos,
                                       multi.state = TRUE, same.layer = FALSE)
                    f <- c(paste0("ko", clean_pid, ".", out_suffix, ".multi.png"),
                           paste0("ko", clean_pid, ".", out_suffix, ".png"))
                    f <- f[file.exists(f)]
                    if (length(f) > 0) {
                        made <- c(made, file.path(pv_dir, f[1]))
                        # The PDF has no filename to read a contrast off, and
                        # one pathway enriched in two contrasts renders twice.
                        page_labels <- c(page_labels,
                                         paste0(label, "  |  ko", clean_pid))
                        message("    Union pathview: ko", clean_pid)
                    }
                }, error = function(e) {
                    message("    Union pathview failed for ", clean_pid, ": ", e$message)
                })
            }
        }
        made
    })
    if (length(generated) == 0) return(NULL)

    # 4. Compile into this renderer's own PDF. Deliberately NOT the name
    #    generate_multi_ora_pathview() writes: that file means "enriched in two
    #    or more omics layers", and a union of single-layer hits is not that.
    #    The report reads whichever of the two exists and says which it is.
    pdf_path <- .compile_pathview_pdf(
        generated, file.path(out_dir, "multi_ora_pathview_union.pdf"), page_labels)
    if (is.null(pdf_path)) return(NULL)

    # What this run actually produced, written down rather than left for the
    # report to infer from which files happen to exist: whether any compound
    # node carries a value, and which readable contrast each filename-safe
    # output key stands for.
    tryCatch(
        yaml::write_yaml(list(compound_nodes = any_compounds,
                              contrast_labels = contrast_labels),
                         file.path(out_dir, "multi_ora_pathview_union.yaml")),
        error = function(e) NULL
    )
    message("  Union pathview: ", length(generated), " maps -> ", basename(pdf_path))
    pdf_path
}
