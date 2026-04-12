#!/usr/bin/env Rscript
#' Extract KEGG Compounds for an Organism
#'
#' This script extracts all compounds associated with a specific organism
#' from KEGG, organized by pathway. Outputs a TSV with compound IDs, names,
#' formulas, and pathway associations.
#'
#' Usage:
#'   Rscript extract_kegg_compounds.R --kegg-org hsa --output hsa_compounds.tsv
#'   Rscript extract_kegg_compounds.R --kegg-org cel --format gmt --output cel_compounds.gmt
#'   Rscript extract_kegg_compounds.R --list-kegg --search "elegans"
#'
#' @examples
#' # From R:
#' source("extract_kegg_compounds.R")
#' compounds <- extract_kegg_compounds("hsa")
#' compounds <- extract_kegg_compounds("cel", output_file = "cel_compounds.tsv")
#' compound_gmt <- extract_kegg_compound_gmt("rno", output_file = "rno_compound_pathways.gmt")

# ==============================================================================
# SETUP
# ==============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

required_packages <- c("optparse")
missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
}

# ==============================================================================
# MAIN FUNCTIONS
# ==============================================================================

#' Extract all KEGG compounds for an organism
#'
#' Retrieves compounds from all pathways of a given organism, then fetches
#' compound details (name, formula, exact mass, molecular weight).
#'
#' @param kegg_org KEGG organism code (e.g., "hsa", "mmu", "cel")
#' @param output_file Output TSV file path (NULL for auto-naming)
#' @param fetch_details If TRUE, fetch name/formula/mass for each compound
#' @return Data frame of compounds with pathway associations
#'
#' @export
extract_kegg_compounds <- function(kegg_org,
                                    output_file = NULL,
                                    fetch_details = TRUE) {

    message("\n=== Extracting KEGG Compounds for organism: ", kegg_org, " ===\n")

    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop("KEGGREST package required. Install with: BiocManager::install('KEGGREST')")
    }

    # Get pathway list for the organism
    message("Fetching KEGG pathway list...")
    pathway_list <- tryCatch({
        KEGGREST::keggList("pathway", kegg_org)
    }, error = function(e) {
        stop("Failed to fetch KEGG pathways. Is '", kegg_org, "' a valid KEGG organism code?\n",
             "Check valid codes at: https://www.genome.jp/kegg/catalog/org_list.html")
    })

    message("Found ", length(pathway_list), " pathways")

    # Collect compounds from each pathway
    compound_pathway <- data.frame(
        compound_id = character(),
        pathway_id = character(),
        pathway_name = character(),
        stringsAsFactors = FALSE
    )

    message("Fetching compounds for each pathway...")
    pb <- txtProgressBar(min = 0, max = length(pathway_list), style = 3)

    for (i in seq_along(pathway_list)) {
        pathway_id <- names(pathway_list)[i]
        pathway_name <- pathway_list[i]
        pathway_name <- gsub(paste0(" - .*$"), "", pathway_name)

        tryCatch({
            pathway_info <- KEGGREST::keggGet(pathway_id)[[1]]

            if (!is.null(pathway_info$COMPOUND)) {
                cpd_ids <- names(pathway_info$COMPOUND)
                cpd_names <- unname(pathway_info$COMPOUND)
                clean_pid <- gsub("path:", "", pathway_id)

                new_rows <- data.frame(
                    compound_id = cpd_ids,
                    compound_name_short = cpd_names,
                    pathway_id = clean_pid,
                    pathway_name = pathway_name,
                    stringsAsFactors = FALSE
                )
                compound_pathway <- rbind(compound_pathway, new_rows)
            }
        }, error = function(e) {
            # Skip failed pathways
        })

        setTxtProgressBar(pb, i)

        # Rate limiting for KEGG API
        if (i %% 10 == 0) Sys.sleep(0.3)
    }

    close(pb)

    if (nrow(compound_pathway) == 0) {
        message("No compounds found for organism: ", kegg_org)
        return(invisible(data.frame()))
    }

    # Get unique compounds
    unique_cpds <- unique(compound_pathway$compound_id)
    message("\nFound ", length(unique_cpds), " unique compounds across ",
            length(unique(compound_pathway$pathway_id)), " pathways")

    # Fetch compound details
    if (fetch_details) {
        message("Fetching compound details (name, formula, mass)...")
        cpd_details <- fetch_compound_details(unique_cpds)
    } else {
        cpd_details <- data.frame(
            compound_id = unique_cpds,
            name = NA_character_,
            formula = NA_character_,
            exact_mass = NA_character_,
            mol_weight = NA_character_,
            stringsAsFactors = FALSE
        )
    }

    # Merge details with pathway associations
    result <- merge(compound_pathway, cpd_details, by = "compound_id", all.x = TRUE)

    # Sort by pathway then compound
    result <- result[order(result$pathway_id, result$compound_id), ]
    rownames(result) <- NULL

    # Generate default output filename
    if (is.null(output_file)) {
        output_file <- sprintf("kegg_compounds_%s.tsv", kegg_org)
    }

    write.table(result, file = output_file, sep = "\t",
                row.names = FALSE, quote = FALSE)

    message("\nSuccessfully created: ", output_file)
    message("  Total compound-pathway associations: ", nrow(result))
    message("  Unique compounds: ", length(unique_cpds))
    message("  Pathways with compounds: ", length(unique(result$pathway_id)))

    # Also write a summary of unique compounds
    summary_file <- sub("\\.[^.]+$", "_unique.tsv", output_file)
    unique_result <- cpd_details[order(cpd_details$compound_id), ]
    write.table(unique_result, file = summary_file, sep = "\t",
                row.names = FALSE, quote = FALSE)
    message("  Unique compound list: ", summary_file)

    invisible(result)
}

#' Extract KEGG compounds as a GMT file (pathway -> compounds)
#'
#' Creates a GMT file where each line is a pathway and its associated compounds.
#' Useful for compound set enrichment analysis in metabolomics.
#'
#' @param kegg_org KEGG organism code
#' @param output_file Output GMT file path
#' @return List of pathways with compound members
#'
#' @export
extract_kegg_compound_gmt <- function(kegg_org, output_file = NULL) {

    message("\n=== Generating KEGG Compound GMT for organism: ", kegg_org, " ===\n")

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

    message("Fetching compounds for each pathway...")
    pb <- txtProgressBar(min = 0, max = length(pathway_list), style = 3)

    for (i in seq_along(pathway_list)) {
        pathway_id <- names(pathway_list)[i]
        pathway_name <- pathway_list[i]
        pathway_name <- gsub(paste0(" - .*$"), "", pathway_name)

        tryCatch({
            pathway_info <- KEGGREST::keggGet(pathway_id)[[1]]

            if (!is.null(pathway_info$COMPOUND)) {
                cpd_ids <- names(pathway_info$COMPOUND)
                clean_id <- gsub("path:", "", pathway_id)
                pathways[[clean_id]] <- cpd_ids
                descriptions[clean_id] <- pathway_name
            }
        }, error = function(e) {
            # Skip failed pathways
        })

        setTxtProgressBar(pb, i)
        if (i %% 10 == 0) Sys.sleep(0.3)
    }

    close(pb)

    if (length(pathways) == 0) {
        message("No compounds found for organism: ", kegg_org)
        return(invisible(list()))
    }

    attr(pathways, "descriptions") <- descriptions

    # Generate default output filename
    if (is.null(output_file)) {
        output_file <- sprintf("kegg_compounds_%s.gmt", kegg_org)
    }

    # Write GMT
    write_gmt_file(pathways, output_file)

    message("\nSuccessfully created: ", output_file)
    message("  Pathways with compounds: ", length(pathways))
    message("  Unique compounds: ", length(unique(unlist(pathways))))

    invisible(pathways)
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Fetch compound details from KEGG in batches
#'
#' @param compound_ids Vector of KEGG compound IDs (e.g., "C00001")
#' @return Data frame with compound details
fetch_compound_details <- function(compound_ids) {

    details <- data.frame(
        compound_id = character(),
        name = character(),
        formula = character(),
        exact_mass = character(),
        mol_weight = character(),
        stringsAsFactors = FALSE
    )

    # KEGG allows up to 10 IDs per request
    batch_size <- 10
    batches <- split(compound_ids, ceiling(seq_along(compound_ids) / batch_size))

    pb <- txtProgressBar(min = 0, max = length(batches), style = 3)

    for (b in seq_along(batches)) {
        batch <- batches[[b]]

        tryCatch({
            cpd_info <- KEGGREST::keggGet(batch)

            for (cpd in cpd_info) {
                name <- if (!is.null(cpd$NAME)) {
                    # Take first name, remove trailing semicolon
                    gsub(";\\s*$", "", cpd$NAME[1])
                } else {
                    NA_character_
                }

                new_row <- data.frame(
                    compound_id = cpd$ENTRY %||% NA_character_,
                    name = name,
                    formula = cpd$FORMULA %||% NA_character_,
                    exact_mass = as.character(cpd$EXACT_MASS %||% NA_character_),
                    mol_weight = as.character(cpd$MOL_WEIGHT %||% NA_character_),
                    stringsAsFactors = FALSE
                )
                details <- rbind(details, new_row)
            }
        }, error = function(e) {
            # Add NAs for failed batch
            for (cid in batch) {
                details <- rbind(details, data.frame(
                    compound_id = cid,
                    name = NA, formula = NA,
                    exact_mass = NA, mol_weight = NA,
                    stringsAsFactors = FALSE
                ))
            }
        })

        setTxtProgressBar(pb, b)
        Sys.sleep(0.3)
    }

    close(pb)
    details
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

#' Write GMT file
write_gmt_file <- function(pathways, output_file) {
    descriptions <- attr(pathways, "descriptions") %||% rep("", length(pathways))

    lines <- character(length(pathways))
    for (i in seq_along(pathways)) {
        name <- names(pathways)[i]
        desc <- descriptions[name] %||% ""
        members <- pathways[[i]]
        lines[i] <- paste(c(name, desc, members), collapse = "\t")
    }

    writeLines(lines, output_file)
}

# ==============================================================================
# COMMAND LINE INTERFACE
# ==============================================================================

if (!interactive()) {

    library(optparse)

    option_list <- list(
        make_option(c("-k", "--kegg-org"), type = "character", default = NULL,
                    help = "KEGG organism code (e.g., 'hsa', 'mmu', 'cel')"),
        make_option(c("-f", "--format"), type = "character", default = "tsv",
                    help = "Output format: tsv or gmt [default: tsv]"),
        make_option(c("--output"), type = "character", default = NULL,
                    help = "Output file path"),
        make_option(c("--no-details"), action = "store_true", default = FALSE,
                    help = "Skip fetching compound details (faster, TSV only)"),
        make_option(c("--list-kegg"), action = "store_true", default = FALSE,
                    help = "List available KEGG organisms"),
        make_option(c("--search"), type = "character", default = NULL,
                    help = "Search pattern for organism list")
    )

    opt <- parse_args(OptionParser(
        option_list = option_list,
        description = "Extract all KEGG compounds for a specific organism"
    ))

    # Handle list command
    if (opt$`list-kegg`) {
        list_kegg_organisms(opt$search)
        quit(status = 0)
    }

    # Extract compounds
    if (!is.null(opt$`kegg-org`)) {
        if (opt$format == "gmt") {
            extract_kegg_compound_gmt(
                kegg_org = opt$`kegg-org`,
                output_file = opt$output
            )
        } else {
            extract_kegg_compounds(
                kegg_org = opt$`kegg-org`,
                output_file = opt$output,
                fetch_details = !opt$`no-details`
            )
        }
    } else {
        cat("Usage: Rscript extract_kegg_compounds.R --kegg-org CODE [options]\n")
        cat("       Rscript extract_kegg_compounds.R --list-kegg [--search pattern]\n")
        cat("\nExamples:\n")
        cat("  Rscript extract_kegg_compounds.R --kegg-org hsa --output human_compounds.tsv\n")
        cat("  Rscript extract_kegg_compounds.R --kegg-org cel --format gmt\n")
        cat("  Rscript extract_kegg_compounds.R --list-kegg --search 'elegans'\n")
        cat("\nRun with --help for all options\n")
    }
}
