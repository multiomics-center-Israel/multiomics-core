#!/usr/bin/env Rscript
#' Generate GMT Files for Pathway Analysis
#'
#' This script creates GMT (Gene Matrix Transposed) files for any organism
#' that can be used for pathway enrichment analysis. Especially useful for
#' non-model organisms that lack pre-built pathway databases.
#'
#' Usage:
#'   Rscript scripts/generate_gmt.R --organism "Danio rerio" --output pathways.gmt
#'   Rscript scripts/generate_gmt.R --organism "Bos taurus" --database GO --ontology BP
#'   Rscript scripts/generate_gmt.R --kegg-org bta --output kegg_pathways.gmt
#'
#' @examples
#' # From R:
#' source("scripts/generate_gmt.R")
#' generate_gmt_for_organism("Sus scrofa", output_file = "pig_go_bp.gmt")
#' generate_kegg_gmt("ssc", output_file = "pig_kegg.gmt")

# ==============================================================================
# SETUP
# ==============================================================================

# Null coalescing (must be defined early, used throughout)
`%||%` <- function(x, y) if (is.null(x)) y else x

# Check for required packages
required_packages <- c("optparse")
missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
}

# Source shared utilities
script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) {
        # Fallback for Rscript invocation
        args <- commandArgs(trailingOnly = FALSE)
        file_arg <- grep("^--file=", args, value = TRUE)
        if (length(file_arg) > 0) {
            dirname(normalizePath(sub("^--file=", "", file_arg)))
        } else {
            "."
        }
    }
)
shared_dir <- file.path(dirname(script_dir), "shared", "R")

if (file.exists(file.path(shared_dir, "gmt_utils.R"))) {
    source(file.path(shared_dir, "gmt_utils.R"))
}
if (file.exists(file.path(shared_dir, "organism_detection.R"))) {
    source(file.path(shared_dir, "organism_detection.R"))
}

# ==============================================================================
# MAIN FUNCTIONS
# ==============================================================================

#' Generate GMT file for an organism using biomaRt
#'
#' @param organism Scientific name (e.g., "Danio rerio", "Sus scrofa")
#' @param database "GO" or "all"
#' @param ontology For GO: "BP", "MF", "CC", or "all"
#' @param id_type Gene ID type: "symbol", "ensembl", "entrez"
#' @param output_file Output GMT file path
#' @param min_genes Minimum genes per pathway
#' @param max_genes Maximum genes per pathway
#'
#' @export
generate_gmt_for_organism <- function(organism,
                                       database = "GO",
                                       ontology = "BP",
                                       id_type = "symbol",
                                       output_file = NULL,
                                       min_genes = 10,
                                       max_genes = 500) {

    message("\n=== Generating GMT for ", organism, " ===\n")

    # Check for biomaRt
    if (!requireNamespace("biomaRt", quietly = TRUE)) {
        stop("biomaRt package required. Install with: BiocManager::install('biomaRt')")
    }

    # Get Ensembl dataset name
    dataset <- get_ensembl_dataset_name(organism)
    if (is.null(dataset)) {
        stop("Could not find Ensembl dataset for: ", organism,
             "\nTry running: biomaRt::listDatasets(biomaRt::useMart('ensembl'))")
    }

    message("Using Ensembl dataset: ", dataset)

    # Connect to Ensembl
    mart <- tryCatch({
        biomaRt::useMart("ensembl", dataset = dataset)
    }, error = function(e) {
        stop("Failed to connect to Ensembl: ", e$message)
    })

    # Determine gene attribute
    gene_attr <- switch(id_type,
                        "symbol" = "external_gene_name",
                        "ensembl" = "ensembl_gene_id",
                        "entrez" = "entrezgene_id",
                        "external_gene_name")

    message("Fetching GO annotations (this may take a few minutes)...")

    # Fetch GO annotations
    go_data <- biomaRt::getBM(
        attributes = c(gene_attr, "go_id", "name_1006", "namespace_1003"),
        mart = mart
    )

    message("Retrieved ", nrow(go_data), " gene-GO associations")

    # Filter by ontology
    if (ontology != "all") {
        namespace_map <- c(
            "BP" = "biological_process",
            "MF" = "molecular_function",
            "CC" = "cellular_component"
        )
        if (ontology %in% names(namespace_map)) {
            go_data <- go_data[go_data$namespace_1003 == namespace_map[ontology], ]
            message("Filtered to ", ontology, ": ", nrow(go_data), " associations")
        }
    }

    # Remove empty entries
    go_data <- go_data[!is.na(go_data$go_id) & go_data$go_id != "", ]
    go_data <- go_data[!is.na(go_data[[gene_attr]]) & go_data[[gene_attr]] != "", ]

    # Build pathway list
    message("Building pathway list...")
    pathways <- split(go_data[[gene_attr]], go_data$go_id)
    pathways <- lapply(pathways, unique)

    # Create descriptions
    desc_df <- unique(go_data[, c("go_id", "name_1006")])
    descriptions <- setNames(desc_df$name_1006, desc_df$go_id)

    # Filter by size
    sizes <- sapply(pathways, length)
    keep <- sizes >= min_genes & sizes <= max_genes
    pathways <- pathways[keep]

    message("After size filtering (", min_genes, "-", max_genes, " genes): ",
            length(pathways), " pathways")

    # Add descriptions
    attr(pathways, "descriptions") <- descriptions[names(pathways)]

    # Generate default output filename
    if (is.null(output_file)) {
        org_short <- gsub(" ", "_", tolower(organism))
        output_file <- sprintf("%s_%s_%s.gmt", org_short, database, ontology)
    }

    # Write GMT
    write_gmt_file(pathways, output_file)

    message("\n✓ Successfully created: ", output_file)
    message("  Pathways: ", length(pathways))
    message("  Unique genes: ", length(unique(unlist(pathways))))

    invisible(pathways)
}

#' Generate GMT from KEGG pathways
#'
#' @param kegg_org KEGG organism code (e.g., "hsa", "mmu", "bta")
#' @param id_type "entrez" or "symbol"
#' @param output_file Output GMT file path
#'
#' @export
generate_kegg_gmt <- function(kegg_org,
                               id_type = "entrez",
                               output_file = NULL) {

    message("\n=== Generating KEGG GMT for organism: ", kegg_org, " ===\n")

    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop("KEGGREST package required. Install with: BiocManager::install('KEGGREST')")
    }

    # Get pathway list
    message("Fetching KEGG pathway list...")
    pathway_list <- tryCatch({
        KEGGREST::keggList("pathway", kegg_org)
    }, error = function(e) {
        stop("Failed to fetch KEGG pathways. Is '", kegg_org, "' a valid KEGG organism code?\n",
             "Check valid codes at: https://www.genome.jp/kegg/catalog/org_list.html")
    })

    message("Found ", length(pathway_list), " pathways")

    pathways <- list()
    descriptions <- character()

    message("Fetching genes for each pathway (this may take a while)...")

    pb <- txtProgressBar(min = 0, max = length(pathway_list), style = 3)

    for (i in seq_along(pathway_list)) {
        pathway_id <- names(pathway_list)[i]
        pathway_name <- pathway_list[i]
        pathway_name <- gsub(paste0(" - ", kegg_org, "$"), "", pathway_name)

        tryCatch({
            pathway_info <- KEGGREST::keggGet(pathway_id)[[1]]

            if (!is.null(pathway_info$GENE)) {
                if (id_type == "entrez") {
                    genes <- names(pathway_info$GENE)
                } else {
                    # Extract symbols
                    genes <- sapply(pathway_info$GENE, function(x) {
                        strsplit(x, ";")[[1]][1]
                    })
                    genes <- unname(genes)
                }

                if (length(genes) > 0) {
                    clean_id <- gsub("path:", "", pathway_id)
                    pathways[[clean_id]] <- genes
                    descriptions[clean_id] <- pathway_name
                }
            }
        }, error = function(e) {
            # Skip failed pathways
        })

        setTxtProgressBar(pb, i)

        # Rate limiting
        if (i %% 10 == 0) Sys.sleep(0.3)
    }

    close(pb)

    attr(pathways, "descriptions") <- descriptions

    # Generate default output filename
    if (is.null(output_file)) {
        output_file <- sprintf("kegg_%s.gmt", kegg_org)
    }

    write_gmt_file(pathways, output_file)

    message("\n✓ Successfully created: ", output_file)
    message("  Pathways: ", length(pathways))
    message("  Unique genes: ", length(unique(unlist(pathways))))

    invisible(pathways)
}

#' List available organisms in Ensembl
#'
#' @param pattern Optional pattern to filter organisms
#' @export
list_available_organisms <- function(pattern = NULL) {

    if (!requireNamespace("biomaRt", quietly = TRUE)) {
        stop("biomaRt package required")
    }

    message("Connecting to Ensembl...")
    mart <- biomaRt::useMart("ensembl")
    datasets <- biomaRt::listDatasets(mart)

    if (!is.null(pattern)) {
        matches <- grep(pattern, datasets$description, ignore.case = TRUE)
        datasets <- datasets[matches, ]
    }

    # Format nicely
    datasets$organism <- gsub(" genes.*", "", datasets$description)

    cat("\nAvailable organisms:\n")
    cat(sprintf("%-30s %-40s\n", "Dataset", "Organism"))
    cat(paste(rep("-", 72), collapse = ""), "\n")

    for (i in seq_len(min(50, nrow(datasets)))) {
        cat(sprintf("%-30s %-40s\n",
                    datasets$dataset[i],
                    substr(datasets$organism[i], 1, 40)))
    }

    if (nrow(datasets) > 50) {
        cat(sprintf("\n... and %d more. Use pattern to filter.\n", nrow(datasets) - 50))
    }

    invisible(datasets)
}

#' List available KEGG organisms
#'
#' @param pattern Optional pattern to filter
#' @export
list_kegg_organisms <- function(pattern = NULL) {

    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop("KEGGREST package required")
    }

    message("Fetching KEGG organism list...")
    orgs <- KEGGREST::keggList("organism")

    # Parse into data frame
    org_df <- data.frame(
        code = names(orgs),
        name = unname(orgs),
        stringsAsFactors = FALSE
    )

    if (!is.null(pattern)) {
        matches <- grep(pattern, org_df$name, ignore.case = TRUE)
        org_df <- org_df[matches, ]
    }

    cat("\nAvailable KEGG organisms:\n")
    cat(sprintf("%-6s %-60s\n", "Code", "Organism"))
    cat(paste(rep("-", 68), collapse = ""), "\n")

    for (i in seq_len(min(50, nrow(org_df)))) {
        cat(sprintf("%-6s %-60s\n",
                    org_df$code[i],
                    substr(org_df$name[i], 1, 60)))
    }

    if (nrow(org_df) > 50) {
        cat(sprintf("\n... and %d more. Use pattern to filter.\n", nrow(org_df) - 50))
    }

    invisible(org_df)
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Get Ensembl dataset name for an organism
get_ensembl_dataset_name <- function(organism) {
    organism_lower <- tolower(gsub("\\s+", "_", organism))

    # Common mappings
    known <- list(
        "homo_sapiens" = "hsapiens_gene_ensembl",
        "human" = "hsapiens_gene_ensembl",
        "mus_musculus" = "mmusculus_gene_ensembl",
        "mouse" = "mmusculus_gene_ensembl",
        "rattus_norvegicus" = "rnorvegicus_gene_ensembl",
        "danio_rerio" = "drerio_gene_ensembl",
        "zebrafish" = "drerio_gene_ensembl",
        "gallus_gallus" = "ggallus_gene_ensembl",
        "chicken" = "ggallus_gene_ensembl",
        "sus_scrofa" = "sscrofa_gene_ensembl",
        "pig" = "sscrofa_gene_ensembl",
        "bos_taurus" = "btaurus_gene_ensembl",
        "cow" = "btaurus_gene_ensembl",
        "ovis_aries" = "oaries_gene_ensembl",
        "sheep" = "oaries_gene_ensembl",
        "canis_familiaris" = "cfamiliaris_gene_ensembl",
        "canis_lupus_familiaris" = "cfamiliaris_gene_ensembl",
        "dog" = "cfamiliaris_gene_ensembl",
        "felis_catus" = "fcatus_gene_ensembl",
        "cat" = "fcatus_gene_ensembl",
        "equus_caballus" = "ecaballus_gene_ensembl",
        "horse" = "ecaballus_gene_ensembl",
        "oryctolagus_cuniculus" = "ocuniculus_gene_ensembl",
        "rabbit" = "ocuniculus_gene_ensembl",
        "pan_troglodytes" = "ptroglodytes_gene_ensembl",
        "macaca_mulatta" = "mmulatta_gene_ensembl",
        "drosophila_melanogaster" = "dmelanogaster_gene_ensembl",
        "caenorhabditis_elegans" = "celegans_gene_ensembl",
        "saccharomyces_cerevisiae" = "scerevisiae_gene_ensembl",
        "arabidopsis_thaliana" = "athaliana_gene_ensembl"
    )

    dataset <- known[[organism_lower]]

    # Try to construct if not found
    if (is.null(dataset)) {
        parts <- strsplit(organism_lower, "_")[[1]]
        if (length(parts) >= 2) {
            # Standard Ensembl naming: first letter of genus + species
            dataset <- paste0(substr(parts[1], 1, 1), parts[2], "_gene_ensembl")
        }
    }

    dataset
}

#' Write GMT file
write_gmt_file <- function(pathways, output_file) {
    descriptions <- attr(pathways, "descriptions") %||% rep("", length(pathways))

    lines <- character(length(pathways))
    for (i in seq_along(pathways)) {
        name <- names(pathways)[i]
        desc <- descriptions[name] %||% ""
        genes <- pathways[[i]]
        lines[i] <- paste(c(name, desc, genes), collapse = "\t")
    }

    writeLines(lines, output_file)
}

# Null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# ==============================================================================
# COMMAND LINE INTERFACE
# ==============================================================================

if (!interactive()) {

    library(optparse)

    option_list <- list(
        make_option(c("-o", "--organism"), type = "character", default = NULL,
                    help = "Organism scientific name (e.g., 'Danio rerio')"),
        make_option(c("-k", "--kegg-org"), type = "character", default = NULL,
                    help = "KEGG organism code (e.g., 'dre', 'bta')"),
        make_option(c("-d", "--database"), type = "character", default = "GO",
                    help = "Database: GO, KEGG [default: GO]"),
        make_option(c("--ontology"), type = "character", default = "BP",
                    help = "GO ontology: BP, MF, CC, all [default: BP]"),
        make_option(c("-i", "--id-type"), type = "character", default = "symbol",
                    help = "Gene ID type: symbol, ensembl, entrez [default: symbol]"),
        make_option(c("--output"), type = "character", default = NULL,
                    help = "Output GMT file path"),
        make_option(c("--min-genes"), type = "integer", default = 10,
                    help = "Minimum genes per pathway [default: 10]"),
        make_option(c("--max-genes"), type = "integer", default = 500,
                    help = "Maximum genes per pathway [default: 500]"),
        make_option(c("--list-organisms"), action = "store_true", default = FALSE,
                    help = "List available organisms in Ensembl"),
        make_option(c("--list-kegg"), action = "store_true", default = FALSE,
                    help = "List available KEGG organisms"),
        make_option(c("--search"), type = "character", default = NULL,
                    help = "Search pattern for organism lists")
    )

    opt <- parse_args(OptionParser(option_list = option_list))

    # Handle list commands
    if (opt$`list-organisms`) {
        list_available_organisms(opt$search)
        quit(status = 0)
    }

    if (opt$`list-kegg`) {
        list_kegg_organisms(opt$search)
        quit(status = 0)
    }

    # Generate GMT
    if (!is.null(opt$`kegg-org`)) {
        generate_kegg_gmt(
            kegg_org = opt$`kegg-org`,
            id_type = opt$`id-type`,
            output_file = opt$output
        )
    } else if (!is.null(opt$organism)) {
        generate_gmt_for_organism(
            organism = opt$organism,
            database = opt$database,
            ontology = opt$ontology,
            id_type = opt$`id-type`,
            output_file = opt$output,
            min_genes = opt$`min-genes`,
            max_genes = opt$`max-genes`
        )
    } else {
        cat("Usage: Rscript generate_gmt.R --organism 'Organism name' [options]\n")
        cat("       Rscript generate_gmt.R --kegg-org CODE [options]\n")
        cat("       Rscript generate_gmt.R --list-organisms [--search pattern]\n")
        cat("       Rscript generate_gmt.R --list-kegg [--search pattern]\n")
        cat("\nRun with --help for all options\n")
    }
}
