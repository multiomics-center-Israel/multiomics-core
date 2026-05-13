#' Unified Annotation Utilities for Multi-Omics Pipelines
#'
#' Organism-agnostic annotation functions with fallback chain:
#'   custom file -> biomaRt -> OrgDb -> none (original IDs)
#'
#' Also provides get_organism_info() lookup for database/package mappings.

# Cache for orgdb keys to avoid repeated expensive queries
.annotation_keys_cache <- new.env(parent = emptyenv())

#' Cached version of AnnotationDbi::keys()
#' @keywords internal
.cached_orgdb_keys <- function(orgdb, keytype) {
  org_name <- tryCatch(
    {
      md <- AnnotationDbi::metadata(orgdb)
      md$value[md$name == "ORGANISM"][1]
    },
    error = function(e) "unknown"
  )
  cache_id <- paste0(org_name, ":", keytype)
  if (!exists(cache_id, envir = .annotation_keys_cache)) {
    keys <- tryCatch(
      AnnotationDbi::keys(orgdb, keytype = keytype),
      error = function(e) character(0)
    )
    assign(cache_id, keys, envir = .annotation_keys_cache)
  }
  get(cache_id, envir = .annotation_keys_cache)
}
# ==============================================================================
# ORGANISM INFO LOOKUP
# ==============================================================================

#' Get organism database name for biomaRt/OrgDb/KEGG
#' @param organism User-provided organism string
#' @return List with dataset, orgdb, kegg, taxid, supported
#' @export
get_organism_info <- function(organism) {

    organism_map <- list(
        "human" = list(dataset = "hsapiens_gene_ensembl", orgdb = "org.Hs.eg.db",
                       kegg = "hsa", taxid = 9606),
        "homo sapiens" = list(dataset = "hsapiens_gene_ensembl", orgdb = "org.Hs.eg.db",
                              kegg = "hsa", taxid = 9606),
        "mouse" = list(dataset = "mmusculus_gene_ensembl", orgdb = "org.Mm.eg.db",
                       kegg = "mmu", taxid = 10090),
        "mus musculus" = list(dataset = "mmusculus_gene_ensembl", orgdb = "org.Mm.eg.db",
                              kegg = "mmu", taxid = 10090),
        "rat" = list(dataset = "rnorvegicus_gene_ensembl", orgdb = "org.Rn.eg.db",
                     kegg = "rno", taxid = 10116),
        "rattus norvegicus" = list(dataset = "rnorvegicus_gene_ensembl", orgdb = "org.Rn.eg.db",
                                   kegg = "rno", taxid = 10116),
        "zebrafish" = list(dataset = "drerio_gene_ensembl", orgdb = "org.Dr.eg.db",
                           kegg = "dre", taxid = 7955),
        "danio rerio" = list(dataset = "drerio_gene_ensembl", orgdb = "org.Dr.eg.db",
                              kegg = "dre", taxid = 7955),
        "fly" = list(dataset = "dmelanogaster_gene_ensembl", orgdb = "org.Dm.eg.db",
                     kegg = "dme", taxid = 7227),
        "drosophila melanogaster" = list(dataset = "dmelanogaster_gene_ensembl", orgdb = "org.Dm.eg.db",
                                         kegg = "dme", taxid = 7227),
        "worm" = list(dataset = "celegans_gene_ensembl", orgdb = "org.Ce.eg.db",
                      kegg = "cel", taxid = 6239),
        "caenorhabditis elegans" = list(dataset = "celegans_gene_ensembl", orgdb = "org.Ce.eg.db",
                                        kegg = "cel", taxid = 6239),
        "yeast" = list(dataset = "scerevisiae_gene_ensembl", orgdb = "org.Sc.sgd.db",
                       kegg = "sce", taxid = 559292),
        "saccharomyces cerevisiae" = list(dataset = "scerevisiae_gene_ensembl", orgdb = "org.Sc.sgd.db",
                                          kegg = "sce", taxid = 559292),
        "arabidopsis" = list(dataset = "athaliana_eg_gene", orgdb = "org.At.tair.db",
                             kegg = "ath", taxid = 3702),
        "arabidopsis thaliana" = list(dataset = "athaliana_eg_gene", orgdb = "org.At.tair.db",
                                      kegg = "ath", taxid = 3702),
        "chicken" = list(dataset = "ggallus_gene_ensembl", orgdb = "org.Gg.eg.db",
                         kegg = "gga", taxid = 9031),
        "gallus gallus" = list(dataset = "ggallus_gene_ensembl", orgdb = "org.Gg.eg.db",
                                kegg = "gga", taxid = 9031),
        "pig" = list(dataset = "sscrofa_gene_ensembl", orgdb = "org.Ss.eg.db",
                     kegg = "ssc", taxid = 9823),
        "sus scrofa" = list(dataset = "sscrofa_gene_ensembl", orgdb = "org.Ss.eg.db",
                             kegg = "ssc", taxid = 9823),
        "cow" = list(dataset = "btaurus_gene_ensembl", orgdb = "org.Bt.eg.db",
                     kegg = "bta", taxid = 9913),
        "bos taurus" = list(dataset = "btaurus_gene_ensembl", orgdb = "org.Bt.eg.db",
                             kegg = "bta", taxid = 9913),
        "giardia" = list(dataset = NA, orgdb = NA, kegg = "gla", taxid = 5741),
        "giardia lamblia" = list(dataset = NA, orgdb = NA, kegg = "gla", taxid = 5741)
    )

    org_lower <- tolower(trimws(organism))

    if (org_lower %in% names(organism_map)) {
        info <- organism_map[[org_lower]]
        info$supported <- TRUE
        return(info)
    }

    list(dataset = NA, orgdb = NA, kegg = NA, taxid = NA, supported = FALSE)
}

# ==============================================================================
# MAIN ANNOTATION FUNCTION
# ==============================================================================

#' Annotate features using a flexible fallback chain
#'
#' @param feature_ids Character vector of feature IDs to annotate
#' @param config List containing annotation configuration
#' @param verbose Logical, print progress messages
#' @return Data frame with columns: feature_id, symbol, description, entrez_id, annotation_source
#' @export
annotate_features <- function(feature_ids, config, verbose = TRUE) {

    result <- data.frame(
        feature_id = feature_ids,
        symbol = feature_ids,
        description = NA_character_,
        entrez_id = NA_character_,
        annotation_source = "none",
        stringsAsFactors = FALSE
    )

    if (isTRUE(config$annotation$skip_annotation)) {
        if (verbose) message("Annotation skipped as per config (skip_annotation: true)")
        return(create_annotation_result(result, verbose))
    }

    fallback_chain <- config$annotation$fallback_chain %||%
        c("custom", "biomart", "orgdb", "keggrest")

    annotated_mask <- rep(FALSE, length(feature_ids))

    for (method in fallback_chain) {
        if (all(annotated_mask)) break

        unannotated_ids <- feature_ids[!annotated_mask]

        if (verbose) {
            message(sprintf("Attempting annotation via '%s' for %d features...",
                            method, length(unannotated_ids)))
        }

        method_result <- tryCatch({
            switch(method,
                   "custom" = annotate_from_custom_file(unannotated_ids, config, verbose),
                   "biomart" = annotate_from_biomart(unannotated_ids, config, verbose),
                   "orgdb" = annotate_from_orgdb(unannotated_ids, config, verbose),
                   "keggrest" = annotate_from_keggrest(unannotated_ids, config, verbose),
                   NULL
            )
        }, error = function(e) {
            if (verbose) message(sprintf("  %s annotation failed: %s", method, e$message))
            NULL
        })

        if (!is.null(method_result) && nrow(method_result) > 0) {
            matched <- method_result$feature_id %in% unannotated_ids
            method_result <- method_result[matched, , drop = FALSE]

            if (nrow(method_result) > 0) {
                for (i in seq_len(nrow(method_result))) {
                    idx <- which(result$feature_id == method_result$feature_id[i])
                    if (length(idx) == 1) {
                        if (!is.na(method_result$symbol[i]) && method_result$symbol[i] != "") {
                            result$symbol[idx] <- method_result$symbol[i]
                            result$description[idx] <- method_result$description[i]
                            result$entrez_id[idx] <- method_result$entrez_id[i]
                            result$annotation_source[idx] <- method
                            annotated_mask[idx] <- TRUE
                        }
                    }
                }

                if (verbose) {
                    n_new <- sum(result$annotation_source == method)
                    message(sprintf("  Successfully annotated %d features via %s", n_new, method))
                }
            }
        }
    }

    create_annotation_result(result, verbose)
}

# ==============================================================================
# CUSTOM FILE ANNOTATION
# ==============================================================================

#' Annotate features from a custom CSV file
#' @param feature_ids Character vector of feature IDs
#' @param config Configuration list
#' @param verbose Print messages
#' @return Data frame with annotation columns
annotate_from_custom_file <- function(feature_ids, config, verbose = TRUE) {

    custom_file <- config$annotation$custom_mapping_file %||%
        config$annotation$custom_annotation_file %||%
        config$custom_gene_mapping

    if (is.null(custom_file) || !file.exists(custom_file)) {
        if (verbose) message("  No custom annotation file provided or file not found")
        return(NULL)
    }

    if (verbose) message(sprintf("  Reading custom annotation from: %s", custom_file))

    custom_data <- tryCatch({
        if (grepl("\\.csv$", custom_file, ignore.case = TRUE)) {
            read.csv(custom_file, stringsAsFactors = FALSE)
        } else if (grepl("\\.tsv$|\\.txt$", custom_file, ignore.case = TRUE)) {
            read.delim(custom_file, stringsAsFactors = FALSE)
        } else {
            read.csv(custom_file, stringsAsFactors = FALSE)
        }
    }, error = function(e) {
        if (verbose) message(sprintf("  Failed to read custom file: %s", e$message))
        return(NULL)
    })

    if (is.null(custom_data) || nrow(custom_data) == 0) return(NULL)

    id_col <- detect_id_column(custom_data, config)
    if (is.null(id_col)) {
        if (verbose) message("  Could not detect ID column in custom file")
        return(NULL)
    }

    symbol_col <- detect_column(custom_data, c("symbol", "gene_symbol", "gene_name",
                                                "name", "target_id", "gene"))
    desc_col <- detect_column(custom_data, c("description", "desc", "gene_description",
                                             "full_name", "definition"))
    entrez_col <- detect_column(custom_data, c("entrez_id", "entrezid", "entrez",
                                               "ncbi_gene_id", "gene_id"))

    result <- data.frame(
        feature_id = custom_data[[id_col]],
        symbol = if (!is.null(symbol_col)) custom_data[[symbol_col]] else NA_character_,
        description = if (!is.null(desc_col)) custom_data[[desc_col]] else NA_character_,
        entrez_id = if (!is.null(entrez_col)) as.character(custom_data[[entrez_col]]) else NA_character_,
        stringsAsFactors = FALSE
    )

    result <- result[result$feature_id %in% feature_ids, , drop = FALSE]

    if (verbose && nrow(result) > 0) {
        message(sprintf("  Found %d/%d features in custom file (%.1f%%)",
                        nrow(result), length(feature_ids),
                        100 * nrow(result) / length(feature_ids)))
    }

    result
}

# ==============================================================================
# BIOMART ANNOTATION
# ==============================================================================

#' Annotate features using biomaRt
#' @param feature_ids Character vector of feature IDs
#' @param config Configuration list
#' @param verbose Print messages
#' @return Data frame with annotation columns
annotate_from_biomart <- function(feature_ids, config, verbose = TRUE) {

    if (!requireNamespace("biomaRt", quietly = TRUE)) {
        if (verbose) message("  biomaRt package not installed")
        return(NULL)
    }

    organism <- config$organism %||% config$annotation$organism
    if (is.null(organism)) {
        if (verbose) message("  No organism specified for biomaRt query")
        return(NULL)
    }

    dataset <- get_ensembl_dataset(organism)
    if (is.null(dataset)) {
        if (verbose) message(sprintf("  Could not determine Ensembl dataset for: %s", organism))
        return(NULL)
    }

    id_type <- config$annotation$id_type %||% config$gene_id_type %||% "auto"
    if (id_type == "auto") {
        id_type <- detect_id_type(feature_ids)
    }

    biomart_filter <- get_biomart_filter(id_type)
    if (is.null(biomart_filter)) {
        if (verbose) message(sprintf("  Unsupported ID type for biomaRt: %s", id_type))
        return(NULL)
    }

    if (verbose) {
        message(sprintf("  Connecting to Ensembl (dataset: %s, filter: %s)...",
                        dataset, biomart_filter))
    }

    mart <- tryCatch({
        biomaRt::useMart("ensembl", dataset = dataset)
    }, error = function(e) {
        if (verbose) message(sprintf("  Failed to connect to Ensembl: %s", e$message))
        NULL
    })

    if (is.null(mart)) return(NULL)

    batch_size <- 500
    all_results <- list()

    for (i in seq(1, length(feature_ids), by = batch_size)) {
        batch_ids <- feature_ids[i:min(i + batch_size - 1, length(feature_ids))]

        batch_result <- tryCatch({
            biomaRt::getBM(
                attributes = c(biomart_filter, "external_gene_name",
                               "description", "entrezgene_id"),
                filters = biomart_filter,
                values = batch_ids,
                mart = mart
            )
        }, error = function(e) {
            if (verbose) message(sprintf("  Batch query failed: %s", e$message))
            NULL
        })

        if (!is.null(batch_result) && nrow(batch_result) > 0) {
            all_results[[length(all_results) + 1]] <- batch_result
        }
    }

    if (length(all_results) == 0) return(NULL)

    combined <- do.call(rbind, all_results)

    result <- data.frame(
        feature_id = combined[[biomart_filter]],
        symbol = combined$external_gene_name,
        description = gsub(" \\[Source:.*", "", combined$description),
        entrez_id = as.character(combined$entrezgene_id),
        stringsAsFactors = FALSE
    )

    result <- result[!duplicated(result$feature_id), , drop = FALSE]

    if (verbose) {
        message(sprintf("  Retrieved %d annotations from biomaRt (%.1f%%)",
                        nrow(result), 100 * nrow(result) / length(feature_ids)))
    }

    result
}

# ==============================================================================
# ORGDB ANNOTATION
# ==============================================================================

#' Annotate features using OrgDb packages
#' @param feature_ids Character vector of feature IDs
#' @param config Configuration list
#' @param verbose Print messages
#' @return Data frame with annotation columns
annotate_from_orgdb <- function(feature_ids, config, verbose = TRUE) {

    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
        if (verbose) message("  AnnotationDbi package not installed")
        return(NULL)
    }

    organism <- config$organism %||% config$annotation$organism
    if (is.null(organism)) {
        if (verbose) message("  No organism specified for OrgDb query")
        return(NULL)
    }

    orgdb_pkg <- get_orgdb_package(organism)
    if (is.null(orgdb_pkg)) {
        if (verbose) message(sprintf("  No OrgDb package available for: %s", organism))
        return(NULL)
    }

    if (!requireNamespace(orgdb_pkg, quietly = TRUE)) {
        if (verbose) message(sprintf("  OrgDb package '%s' not installed", orgdb_pkg))
        return(NULL)
    }

    orgdb <- tryCatch({
        getExportedValue(orgdb_pkg, orgdb_pkg)
    }, error = function(e) {
        if (verbose) message(sprintf("  Failed to load OrgDb: %s", e$message))
        NULL
    })

    if (is.null(orgdb)) return(NULL)

    id_type <- config$annotation$id_type %||% config$gene_id_type %||% "auto"
    if (id_type == "auto") {
        id_type <- detect_id_type(feature_ids)
    }

    keytype <- get_orgdb_keytype(id_type)
    if (is.null(keytype) || !(keytype %in% AnnotationDbi::keytypes(orgdb))) {
        if (verbose) message(sprintf("  Keytype '%s' not available in OrgDb", keytype))
        return(NULL)
    }

    if (verbose) message(sprintf("  Querying OrgDb (keytype: %s)...", keytype))

    # Pre-filter to keys that actually exist in the OrgDb to avoid
    # errors from invalid keys (e.g. unmapped protein IDs mixed with symbols)
    valid_keys <- .cached_orgdb_keys(orgdb, keytype)
    # Direct exact match
    exact_ids <- feature_ids[feature_ids %in% valid_keys]

    # Case-insensitive match for remaining IDs (handles rat Actb vs input ACTB)
    remaining <- feature_ids[!feature_ids %in% valid_keys]
    case_matched_ids <- character(0)
    if (length(remaining) > 0 && length(valid_keys) > 0) {
        valid_upper <- toupper(valid_keys)
        # Build lookup: uppercase -> first matching valid key
        dup_mask <- !duplicated(valid_upper)
        case_lookup <- setNames(valid_keys[dup_mask], valid_upper[dup_mask])
        rem_upper <- toupper(remaining)
        found <- rem_upper %in% names(case_lookup)
        if (any(found)) {
            case_matched_ids <- unname(case_lookup[rem_upper[found]])
        }
    }

    usable_ids <- unique(c(exact_ids, case_matched_ids))
    n_exact <- length(exact_ids)
    n_case  <- length(case_matched_ids)

    if (length(usable_ids) == 0) {
        if (verbose) message(sprintf("  0/%d feature IDs found in OrgDb %s keys",
                                     length(feature_ids), keytype))
        return(NULL)
    }
    if (verbose) {
        message(sprintf("  %d/%d feature IDs matched OrgDb %s keys (%d exact, %d case-insensitive)",
                        length(usable_ids), length(feature_ids), keytype, n_exact, n_case))
    }

    result <- tryCatch({
        AnnotationDbi::select(
            orgdb,
            keys = feature_ids,
            columns = c("SYMBOL", "GENENAME", "ENTREZID"),
            keytype = keytype
        )
    }, error = function(e) {
        if (verbose) message(sprintf("  OrgDb query failed: %s", e$message))
        NULL
    })

    if (is.null(result) || nrow(result) == 0) return(NULL)

    result <- data.frame(
        feature_id = result[[keytype]],
        symbol = result$SYMBOL,
        description = result$GENENAME,
        entrez_id = as.character(result$ENTREZID),
        stringsAsFactors = FALSE
    )

    result <- result[!duplicated(result$feature_id), , drop = FALSE]

    if (verbose) {
        message(sprintf("  Retrieved %d annotations from OrgDb (%.1f%%)",
                        nrow(result), 100 * nrow(result) / length(feature_ids)))
    }

    result
}

# ==============================================================================
# KEGGREST ANNOTATION (FOR NON-MODEL ORGANISMS)
# ==============================================================================

#' Annotate features using KEGG Orthology via KEGGREST
#'
#' For organisms without biomaRt/OrgDb support but with a KEGG code,
#' this fetches gene -> KO mappings and uses KO descriptions as gene names.
#'
#' @param feature_ids Character vector of feature IDs
#' @param config Configuration list
#' @param verbose Print messages
#' @return Data frame with annotation columns
annotate_from_keggrest <- function(feature_ids, config, verbose = TRUE) {

    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        if (verbose) message("  KEGGREST package not installed")
        return(NULL)
    }

    organism <- config$organism %||% config$annotation$organism
    if (is.null(organism)) {
        if (verbose) message("  No organism specified for KEGGREST query")
        return(NULL)
    }

    org_info <- get_organism_info(organism)
    if (is.na(org_info$kegg)) {
        if (verbose) message("  No KEGG code for organism: ", organism)
        return(NULL)
    }

    kegg_code <- org_info$kegg

    if (verbose) message(sprintf("  Fetching KO annotations from KEGG for '%s'...", kegg_code))

    # Get gene -> KO links
    ko_links <- tryCatch(
        KEGGREST::keggLink("ko", kegg_code),
        error = function(e) { if (verbose) message("  keggLink failed: ", e$message); NULL }
    )
    if (is.null(ko_links) || length(ko_links) == 0) return(NULL)

    gene_ids_kegg <- sub(paste0("^", kegg_code, ":"), "", names(ko_links))
    ko_ids <- sub("^ko:", "", unname(ko_links))

    # Get KO descriptions in bulk
    ko_descs <- tryCatch(
        KEGGREST::keggList("ko"),
        error = function(e) { if (verbose) message("  keggList(ko) failed: ", e$message); NULL }
    )
    if (is.null(ko_descs)) return(NULL)

    # Build gene -> description mapping
    ko_desc_map <- ko_descs[ko_ids]

    # Parse KO descriptions: "SYMBOL; description [EC:x.x.x.x]" or "SYMBOL, ALT; description"
    symbols <- vapply(ko_desc_map, function(d) {
        if (is.na(d)) return(NA_character_)
        # Take text before first ";"  as symbol, rest as description
        parts <- strsplit(d, ";\\s*", perl = TRUE)[[1]]
        if (length(parts) >= 1) {
            # First part may have "SYM1, SYM2" — take first symbol
            sym <- trimws(strsplit(parts[1], ",")[[1]][1])
            return(sym)
        }
        NA_character_
    }, character(1), USE.NAMES = FALSE)

    descriptions <- vapply(ko_desc_map, function(d) {
        if (is.na(d)) return(NA_character_)
        parts <- strsplit(d, ";\\s*", perl = TRUE)[[1]]
        if (length(parts) >= 2) {
            desc <- trimws(parts[2])
            # Remove trailing [EC:...] for cleaner description
            desc <- sub("\\s*\\[EC:[^]]+\\]$", "", desc)
            return(desc)
        }
        NA_character_
    }, character(1), USE.NAMES = FALSE)

    result <- data.frame(
        feature_id = gene_ids_kegg,
        symbol = symbols,
        description = descriptions,
        entrez_id = NA_character_,
        stringsAsFactors = FALSE
    )

    # Keep only features in our input set
    result <- result[result$feature_id %in% feature_ids, , drop = FALSE]
    # Remove rows with NA symbol
    result <- result[!is.na(result$symbol), , drop = FALSE]
    # Deduplicate (a gene might map to multiple KOs — keep first)
    result <- result[!duplicated(result$feature_id), , drop = FALSE]

    if (verbose && nrow(result) > 0) {
        message(sprintf("  Found %d/%d features via KEGGREST KO (%.1f%%)",
                        nrow(result), length(feature_ids),
                        100 * nrow(result) / length(feature_ids)))
    }

    result
}

# ==============================================================================
# V2 ANNOTATION WRAPPER (PREFERRED ENTRY POINT)
# ==============================================================================

#' Annotate genes using the fallback chain (preferred method)
#'
#' @param gene_ids Vector of gene IDs
#' @param config Pipeline configuration list
#' @param verbose Print progress messages
#' @return List with annotation data frame and statistics
#' @export
annotate_genes_v2 <- function(gene_ids, config, verbose = TRUE) {

    annotation_config <- list(
        organism = config$organism,
        annotation = list(
            skip_annotation = isTRUE(config$skip_annotation) ||
                              isTRUE(config$annotation$skip_annotation),
            custom_mapping_file = config$custom_gene_mapping %||%
                                  config$annotation$custom_mapping_file,
            custom_annotation_file = config$annotation$custom_annotation_file,
            fallback_chain = config$annotation$fallback_chain %||%
                             c("custom", "biomart", "orgdb", "keggrest"),
            id_type = config$gene_id_type %||%
                      config$annotation$id_type %||%
                      "auto"
        )
    )

    # Auto-detect ID type if needed
    if (annotation_config$annotation$id_type == "auto") {
        detected_type <- detect_id_type(gene_ids)
        if (verbose) message("Auto-detected ID type: ", detected_type)
        annotation_config$annotation$id_type <- detected_type
    }

    # Auto-detect organism if not specified or set to "auto"
    if (is.null(annotation_config$organism) ||
        tolower(annotation_config$organism) == "auto") {
        org_result <- detect_organism_from_ids(gene_ids)
        if (org_result$confidence %in% c("high", "medium")) {
            if (verbose) message("Auto-detected organism: ", org_result$organism,
                                 " (", org_result$confidence, " confidence)")
            annotation_config$organism <- org_result$organism
        }
    }

    annotation_result <- annotate_features(gene_ids, annotation_config, verbose = verbose)

    annotation <- data.frame(
        gene_id = annotation_result$feature_id,
        symbol = annotation_result$symbol,
        description = annotation_result$description,
        entrez_id = annotation_result$entrez_id,
        stringsAsFactors = FALSE
    )

    stats <- attr(annotation_result, "annotation_stats")
    if (is.null(stats)) {
        stats <- list(
            source = "shared_utils",
            total_genes = length(gene_ids),
            annotated_genes = sum(annotation_result$annotation_source != "none"),
            success_rate = mean(annotation_result$annotation_source != "none")
        )
    } else {
        stats <- list(
            source = paste(names(stats$sources), collapse = "+"),
            total_genes = stats$total_features,
            annotated_genes = stats$annotated_features,
            success_rate = stats$coverage
        )
    }

    list(annotation = annotation, stats = stats)
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Detect ID column in a data frame
detect_id_column <- function(df, config = NULL) {
    if (!is.null(config$annotation$id_column)) {
        if (config$annotation$id_column %in% names(df)) {
            return(config$annotation$id_column)
        }
    }

    id_patterns <- c("feature_id", "source_id", "gene_id", "ensembl_gene_id",
                     "ensembl_id", "uniprot", "protein_id", "id", "gene",
                     "ensembl", "entrez_id", "entrezid")

    for (pattern in id_patterns) {
        matches <- grep(paste0("^", pattern, "$"), names(df), ignore.case = TRUE)
        if (length(matches) > 0) return(names(df)[matches[1]])
    }

    if (ncol(df) > 0) return(names(df)[1])
    NULL
}

#' Detect a column by name patterns
detect_column <- function(df, patterns) {
    for (pattern in patterns) {
        matches <- grep(paste0("^", pattern, "$"), names(df), ignore.case = TRUE)
        if (length(matches) > 0) return(names(df)[matches[1]])
    }
    NULL
}

#' Get Ensembl dataset name for an organism
get_ensembl_dataset <- function(organism) {
    organism <- tolower(gsub("\\s+", "_", organism))

    datasets <- list(
        "homo_sapiens" = "hsapiens_gene_ensembl",
        "human" = "hsapiens_gene_ensembl",
        "mus_musculus" = "mmusculus_gene_ensembl",
        "mouse" = "mmusculus_gene_ensembl",
        "rattus_norvegicus" = "rnorvegicus_gene_ensembl",
        "rat" = "rnorvegicus_gene_ensembl",
        "danio_rerio" = "drerio_gene_ensembl",
        "zebrafish" = "drerio_gene_ensembl",
        "drosophila_melanogaster" = "dmelanogaster_gene_ensembl",
        "fly" = "dmelanogaster_gene_ensembl",
        "caenorhabditis_elegans" = "celegans_gene_ensembl",
        "worm" = "celegans_gene_ensembl",
        "saccharomyces_cerevisiae" = "scerevisiae_gene_ensembl",
        "yeast" = "scerevisiae_gene_ensembl",
        "arabidopsis_thaliana" = "athaliana_gene_ensembl",
        "gallus_gallus" = "ggallus_gene_ensembl",
        "sus_scrofa" = "sscrofa_gene_ensembl",
        "bos_taurus" = "btaurus_gene_ensembl",
        "canis_familiaris" = "cfamiliaris_gene_ensembl"
    )

    dataset <- datasets[[organism]]

    if (is.null(dataset)) {
        parts <- strsplit(organism, "_")[[1]]
        if (length(parts) >= 2) {
            dataset <- paste0(substr(parts[1], 1, 1), parts[2], "_gene_ensembl")
        }
    }

    dataset
}

#' Get OrgDb package name for an organism
get_orgdb_package <- function(organism) {
    organism <- tolower(gsub("\\s+", "_", organism))

    packages <- list(
        "homo_sapiens" = "org.Hs.eg.db",
        "human" = "org.Hs.eg.db",
        "mus_musculus" = "org.Mm.eg.db",
        "mouse" = "org.Mm.eg.db",
        "rattus_norvegicus" = "org.Rn.eg.db",
        "rat" = "org.Rn.eg.db",
        "danio_rerio" = "org.Dr.eg.db",
        "zebrafish" = "org.Dr.eg.db",
        "drosophila_melanogaster" = "org.Dm.eg.db",
        "fly" = "org.Dm.eg.db",
        "caenorhabditis_elegans" = "org.Ce.eg.db",
        "worm" = "org.Ce.eg.db",
        "saccharomyces_cerevisiae" = "org.Sc.sgd.db",
        "yeast" = "org.Sc.sgd.db",
        "arabidopsis_thaliana" = "org.At.tair.db",
        "gallus_gallus" = "org.Gg.eg.db",
        "sus_scrofa" = "org.Ss.eg.db",
        "bos_taurus" = "org.Bt.eg.db",
        "canis_familiaris" = "org.Cf.eg.db"
    )

    packages[[organism]]
}

#' Get biomaRt filter name for an ID type
get_biomart_filter <- function(id_type) {
    filters <- list(
        "ensembl_gene_id" = "ensembl_gene_id",
        "ensembl" = "ensembl_gene_id",
        "entrez" = "entrezgene_id",
        "entrez_id" = "entrezgene_id",
        "entrezid" = "entrezgene_id",
        "symbol" = "external_gene_name",
        "gene_symbol" = "external_gene_name",
        "uniprot" = "uniprotswissprot",
        "uniprot_id" = "uniprotswissprot",
        "refseq" = "refseq_mrna",
        "refseq_mrna" = "refseq_mrna",
        "wormbase" = "wormbase_gene"
    )

    filters[[tolower(id_type)]]
}

#' Get OrgDb keytype for an ID type
get_orgdb_keytype <- function(id_type) {
    keytypes <- list(
        "ensembl_gene_id" = "ENSEMBL",
        "ensembl" = "ENSEMBL",
        "entrez" = "ENTREZID",
        "entrez_id" = "ENTREZID",
        "entrezid" = "ENTREZID",
        "symbol" = "SYMBOL",
        "gene_symbol" = "SYMBOL",
        "uniprot" = "UNIPROT",
        "uniprot_id" = "UNIPROT",
        "refseq" = "REFSEQ",
        "wormbase" = "WORMBASE"
    )

    keytypes[[tolower(id_type)]]
}

#' Create final annotation result with statistics
create_annotation_result <- function(result, verbose = TRUE) {

    n_total <- nrow(result)
    n_annotated <- sum(result$annotation_source != "none")
    coverage <- n_annotated / n_total

    attr(result, "annotation_stats") <- list(
        total_features = n_total,
        annotated_features = n_annotated,
        coverage = coverage,
        sources = table(result$annotation_source)
    )

    if (verbose) {
        message(sprintf("\n=== Annotation Summary ==="))
        message(sprintf("Total features: %d", n_total))
        message(sprintf("Annotated: %d (%.1f%%)", n_annotated, 100 * coverage))
        message(sprintf("Annotation sources:"))
        for (src in names(attr(result, "annotation_stats")$sources)) {
            n <- attr(result, "annotation_stats")$sources[[src]]
            message(sprintf("  %s: %d", src, n))
        }
    }

    result
}

#' Get annotation coverage from result
#' @param annotation_result Result from annotate_features()
#' @return Numeric, proportion of features with annotation (0-1)
#' @export
get_annotation_coverage <- function(annotation_result) {
    stats <- attr(annotation_result, "annotation_stats")
    if (is.null(stats)) {
        return(sum(annotation_result$annotation_source != "none") / nrow(annotation_result))
    }
    stats$coverage
}

#' Check if annotation coverage is sufficient for enrichment
#' @param annotation_result Result from annotate_features() or annotate_genes_v2()
#' @param min_coverage Minimum required coverage (default 0.3)
#' @param verbose Print message
#' @return Logical
#' @export
check_enrichment_ready <- function(annotation_result, min_coverage = 0.3, verbose = TRUE) {
    if (is.list(annotation_result) && "stats" %in% names(annotation_result)) {
        coverage <- annotation_result$stats$success_rate
    } else {
        coverage <- get_annotation_coverage(annotation_result)
    }

    if (coverage < min_coverage) {
        if (verbose) {
            message(sprintf(
                "WARNING: Annotation coverage (%.1f%%) is below threshold (%.1f%%) for enrichment analysis.",
                100 * coverage, 100 * min_coverage
            ))
            message("Consider providing a custom GMT file or annotation mapping.")
        }
        return(FALSE)
    }

    TRUE
}
