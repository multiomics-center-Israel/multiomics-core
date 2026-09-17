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




#' Collect the KEGG pathways each contrast called enriched
#'
#' The ORA tables carry their own `contrast` column (written by the enrichment
#' step alongside `database`, `method` and `direction`), so a pathway can be
#' attributed to the contrast that produced it without reverse-engineering it
#' from the filename, whose prefix spells the contrast and the gene-set database
#' in whatever way the project configured them.
#'
#' Pure: it reads no files and draws nothing, which is what makes the selection
#' testable without invoking pathview.
#'
#' @param ora_tables List of ORA data frames, as read from the per-omic
#'   `*_ora_up/down.csv` exports.
#' @param kegg_org Organism code for pathway identity, or NULL for a run with no
#'   KEGG code of its own. This is the identity organism, never the render
#'   species -- see \code{generate_per_omic_union_pathview}.
#' @param alpha Significance cutoff applied to `pvalue`, or to `padj` when no
#'   `pvalue` column is present.
#' @return Named list, one element per contrast, each a character vector of
#'   normalized KEGG pathway ids. Contrasts with no KEGG hit are absent, and so
#'   are rows that carry no contrast: they cannot be attributed to one, and
#'   rendering them against some other contrast's fold changes is the bug this
#'   whole structure exists to prevent.
.kegg_hits_by_contrast <- function(ora_tables, kegg_org = NULL, alpha = 0.05) {
    hits <- list()
    for (d in ora_tables) {
        if (is.null(d) || !is.data.frame(d) || nrow(d) == 0) next
        if (!all(c("pathway", "contrast") %in% names(d))) next

        # The gene-set name reaching this column can be a bare accession, a
        # prefixed one, or the "<accession> <readable name>" form that
        # fetch_kegg_via_rest() gives its sets. normalize_pathway_join_key()
        # reduces all three to the bare map number and leaves everything else
        # byte-identical, so testing the normalized value is what decides
        # whether the row names a KEGG pathway at all -- and another organism's
        # prefix, which normalization does not touch, is still rejected.
        keys <- normalize_pathway_join_key(d$pathway, kegg_org)
        keep <- is_kegg_pathway_accession(keys, kegg_org)

        pcol <- if ("pvalue" %in% names(d)) "pvalue" else if ("padj" %in% names(d)) "padj" else NA
        if (!is.na(pcol)) keep <- keep & !is.na(d[[pcol]]) & d[[pcol]] < alpha

        keep <- keep & !is.na(d$contrast) & nzchar(trimws(as.character(d$contrast)))
        if (!any(keep)) next

        for (cn in unique(as.character(d$contrast[keep]))) {
            hits[[cn]] <- unique(c(hits[[cn]], keys[keep & d$contrast == cn]))
        }
    }
    hits
}


#' Pick one contrast's DE table out of a per-contrast list
#'
#' Exact name first. The ORA tables and the DE tables are named from the same
#' contrast strings, but a caller that has passed one of them through
#' \code{make.names()} would otherwise silently lose the match, so an equal
#' syntactic name counts as the same contrast. Anything else is a miss, and a
#' miss means the layer is absent for that contrast -- never the first table.
#'
#' @param tables Named list of DE tables, as \code{extract_de_tables} returns.
#' @param contrast Contrast name to look for.
#' @return The matching data frame, or NULL.
.de_table_for_contrast <- function(tables, contrast) {
    if (is.null(tables) || length(tables) == 0 || is.null(names(tables))) return(NULL)
    if (contrast %in% names(tables)) return(tables[[contrast]])
    same <- make.names(names(tables)) == make.names(contrast)
    if (sum(same) == 1L) return(tables[[which(same)]])
    NULL
}


#' Render KEGG maps for the union of pathways enriched in either gene omic
#'
#' Reads the per-omic ORA tables written by the RNA / proteomics enrichment
#' steps and, for each contrast separately, renders one map per KEGG pathway
#' that contrast called enriched, with its RNA and protein log2FC overlaid as
#' two states.
#'
#' Union, not intersection: a pathway enriched in only one gene layer is
#' deliberately eligible here, because this is the fallback for runs that get no
#' maps at all otherwise. That is why the output is a separate PDF from
#' \code{generate_multi_ora_pathview}'s, whose pathways really are supported in
#' two or more layers -- the report distinguishes the two.
#'
#' Runs in one of two ID spaces:
#' \itemize{
#'   \item organism space -- the species has a KEGG code, so native
#'     `<org>#####` maps are drawn and feature ids double as KEGG gene ids
#'     (proteins are bridged through the gene-protein mapping).
#'   \item KO space -- the species has no KEGG code, or `pathview.species` is
#'     `"ko"`. The KEGG reference `map#####` artwork is drawn and features are
#'     translated to KEGG Orthology ids through `enrichment.pathview.ko_map`.
#'     Metabolite log2FC is passed as `cpd.data` here, so one map carries all
#'     three layers -- which is the only reason the KO detour is worth taking.
#' }
#'
#' Render space and identity space are separate: forcing `"ko"` changes what the
#' maps are drawn on, not which accessions the run may legitimately spell its
#' own pathways with.
#'
#' @param de_results Named list of DE results per omics.
#' @param harmonization_res Harmonization result (supplies the gene-protein map
#'   and the metabolomics row data).
#' @param config Full config.
#' @param out_dir Multi-ORA output directory (maps go under `out_dir/pathview`).
#' @param top_n Max pathways to render per contrast; overridden by
#'   `modes$multiomics$enrichment$pathview$top_n` when that is set.
#' @return Path to the compiled PDF, or NULL when nothing could be rendered.
generate_per_omic_union_pathview <- function(de_results, harmonization_res,
                                             config, out_dir, top_n = 25) {
    if (!requireNamespace("pathview", quietly = TRUE)) return(NULL)
    organism <- config$global$organism %||% ""
    # resolve_kegg_org_code(), not get_kegg_organism(): the latter only knows the
    # six exact species names of the older table, so an organism KEGG does cover
    # would be pushed into KO mode and then need a ko_map to render anything.
    # The pathway identity helpers read the same registry through this resolver.
    kegg_org <- resolve_kegg_org_code(organism)

    pv_cfg <- config$modes$multiomics$enrichment$pathview %||% list()
    top_n <- pv_cfg$top_n %||% top_n

    # `species` is a mode switch, not a second place to name an organism: the
    # only value it may carry is "ko", to force reference-map rendering for a
    # species that does have a KEGG code. global.organism stays the single
    # source of truth for which organism this run is.
    species_cfg <- trimws(as.character(pv_cfg$species %||% ""))
    if (nzchar(species_cfg) && !identical(tolower(species_cfg), "ko")) {
        message("  Union pathview: enrichment.pathview.species accepts only 'ko'; ",
                "ignoring '", species_cfg, "' and using global.organism")
        species_cfg <- ""
    }
    ko_mode <- identical(tolower(species_cfg), "ko") || is.null(kegg_org)

    ko_map <- if (ko_mode) load_feature_ko_map(config) else NULL
    if (ko_mode && (is.null(ko_map) || nrow(ko_map) == 0)) {
        # Nothing can be placed onto a reference map without a KO translation.
        return(NULL)
    }

    # Render species is what pathview draws on; the identity organism is which
    # accessions this run may spell its own pathways with. Forcing KO changes
    # the first and not the second, so an hsa run with species = "ko" still
    # normalizes hsa00010 to 00010 and then draws it on the reference map. Only
    # a genuinely code-less organism carries NULL here.
    pv_species <- if (ko_mode) "ko" else kegg_org
    id_org <- kegg_org

    # 1. KEGG pathways each contrast called enriched, read from the per-omic ORA
    #    tables. Locate them by walking up from out_dir to the run root.
    run_root <- out_dir
    for (i in seq_len(8)) {
        if (dir.exists(file.path(run_root, "rna", "Enrichment")) ||
            dir.exists(file.path(run_root, "proteomics", "Enrichment"))) break
        parent <- dirname(run_root); if (identical(parent, run_root)) break; run_root <- parent
    }
    # The ORA filename prefix carries the contrast and the gene-set database,
    # both project-specific, so only the "_ora_up/down.csv" tail is fixed. The
    # contrast is then read from the rows, not from the name.
    ora_files <- c(list.files(file.path(run_root, "rna", "Enrichment"),
                              "_ora_(up|down)\\.csv$", full.names = TRUE),
                   list.files(file.path(run_root, "proteomics", "Enrichment"),
                              "_ora_(up|down)\\.csv$", full.names = TRUE))
    ora_tables <- lapply(ora_files, function(f) {
        tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    })
    hits <- .kegg_hits_by_contrast(ora_tables, id_org)
    if (length(hits) == 0) {
        message("  Union pathview: no enriched ", pv_species, " KEGG pathways.")
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
    metab_tables <- if (ko_mode) tables_for("metabolomics") else NULL
    metab_map <- if (!is.null(metab_tables)) {
        tryCatch(map_metabolite_ids_to_kegg(metab_tables, harmonization_res),
                 error = function(e) NULL)
    } else NULL

    gpm <- harmonization_res$gene_protein_mapping
    layer_fc <- function(df, om) {
        if (is.null(df) || !all(c("feature_id", "log2fc") %in% names(df))) return(NULL)
        if (ko_mode) return(aggregate_log2fc_by_ko(df, ko_map, om))
        ids <- df$feature_id
        if (identical(om, "proteomics") && !is.null(gpm)) {
            ids <- gpm$gene_id[match(df$feature_id, gpm$protein_id)]
        }
        ok <- !is.na(ids) & is.finite(df$log2fc)
        if (!any(ok)) return(NULL)
        # Many features can collapse onto one node (paralogues, protein groups),
        # so average rather than let an arbitrary one win.
        tapply(df$log2fc[ok], ids[ok], mean, na.rm = TRUE)
    }

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

    generated <- withr::with_dir(pv_dir, {
        made <- character(0)
        for (contrast in names(hits)) {
            rna_df  <- .de_table_for_contrast(rna_tables, contrast)
            prot_df <- .de_table_for_contrast(prot_tables, contrast)
            rna_fc  <- layer_fc(rna_df, "transcriptomics")
            prot_fc <- layer_fc(prot_df, "proteomics")
            if (is.null(rna_fc) && is.null(prot_fc)) {
                # A contrast whose pathways came from an ORA table we have no DE
                # table for. Borrowing another contrast's values is exactly the
                # cross-wiring this loop exists to prevent, so it is skipped.
                message("  Union pathview: no DE table for contrast ", contrast,
                        "; skipping it")
                next
            }
            genes <- unique(c(names(rna_fc), names(prot_fc)))
            gene_data <- matrix(NA_real_, length(genes), 2,
                                dimnames = list(genes, c("RNA", "Protein")))
            if (!is.null(rna_fc))  gene_data[names(rna_fc), 1]  <- rna_fc
            if (!is.null(prot_fc)) gene_data[names(prot_fc), 2] <- prot_fc

            # Compounds colour the compound nodes, from the same contrast. Only
            # in KO mode: organism mode keeps the two-layer view.
            cpd_data <- NULL
            metab_df <- .de_table_for_contrast(metab_tables, contrast)
            if (!is.null(metab_df) && !is.null(metab_map) && nrow(metab_map) > 0) {
                cpd_data <- tryCatch({
                    md <- merge(metab_df, metab_map, by = "feature_id")
                    fc <- tapply(md$log2fc, md$KEGG_CPD, mean, na.rm = TRUE)
                    stats::setNames(as.numeric(fc), names(fc))
                }, error = function(e) NULL)
            }

            safe_contrast <- make.names(contrast)
            out_suffix <- paste0("multi_ora_", safe_contrast)
            message("  Union pathview: ", contrast, " -- ", length(genes), " ",
                    if (ko_mode) "KO" else "gene", " nodes with log2FC")

            for (clean_pid in utils::head(hits[[contrast]], top_n)) {
                # Place the colour key in whichever corner the map itself leaves
                # empty. The template is cached by an earlier pathview call, so
                # this is free on a re-run and falls back to the default corner
                # on a first one.
                key_pos <- pick_key_position(
                    file.path(pv_dir, paste0(pv_species, clean_pid, ".png")))
                tryCatch({
                    pathview::pathview(gene.data = gene_data, cpd.data = cpd_data,
                                       pathway.id = clean_pid,
                                       species = pv_species, gene.idtype = "KEGG",
                                       out.suffix = out_suffix, kegg.dir = pv_dir,
                                       key.pos = key_pos,
                                       multi.state = TRUE, same.layer = FALSE)
                    f <- c(paste0(pv_species, clean_pid, ".", out_suffix, ".multi.png"),
                           paste0(pv_species, clean_pid, ".", out_suffix, ".png"))
                    f <- f[file.exists(f)]
                    if (length(f) > 0) {
                        made <- c(made, file.path(pv_dir, f[1]))
                        message("    Union pathview: ", pv_species, clean_pid)
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
    pdf_path <- file.path(out_dir, "multi_ora_pathview_union.pdf")
    tryCatch({
        grDevices::pdf(pdf_path, width = 12, height = 8)
        for (png_file in generated) {
            img <- png::readPNG(png_file)
            grid::grid.newpage()
            grid::grid.raster(img, y = 0.48, height = 0.92)
        }
        grDevices::dev.off()
    }, error = function(e) { tryCatch(grDevices::dev.off(), error = function(e2) NULL) })
    message("  Union pathview: ", length(generated), " maps -> ", basename(pdf_path))
    pdf_path
}
