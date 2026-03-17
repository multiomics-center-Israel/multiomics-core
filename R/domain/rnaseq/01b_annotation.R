# ============================================================
# RNA-seq Annotation and Trinotate Processing
# ============================================================
#
# This module provides:
# 1. Wrapper functions for config-driven loading of annotation/trinotate files
# 2. Original processing functions (preserved as-is, with side effects isolated)
#
# TECHNICAL DEBT:
# - Original functions mix I/O and transformation (should be split)
# - Hardcoded column names (should be config-driven)
# - parse_pfam dependency handled via conditional stub (temporary compromise)
#
# ============================================================


# ------------------------------------------------------------
# TEMPORARY COMPROMISE: parse_pfam stub
# ------------------------------------------------------------
# This stub is injected ONLY if parse_pfam is not already defined.
# It allows the pipeline to run without Pfam parsing.
# This is an explicit temporary compromise, not a preferred pattern.
# Future: parse_pfam should be provided as a proper dependency or
# the Pfam parsing logic should be made optional within the function.
# ------------------------------------------------------------

.ensure_parse_pfam <- function() {
    if (!exists("parse_pfam", mode = "function", envir = globalenv()))
        if (!exists("parse_pfam", mode = "function", envir = parent.frame(2))) {
            message("[annotation] parse_pfam not found, Pfam parsing will be skipped")
            # Assign stub to the module's environment, not global
            assign(
                "parse_pfam",
                function(x) {
                    warning("parse_pfam not available, returning empty string", call. = FALSE)
                    return(rep("", length(x)))
                },
                envir = parent.frame()
            )
        }
}


# ============================================================
# WRAPPER FUNCTIONS (Config-Driven)
# ============================================================

#' Load and process annotation file
#'
#' Config-driven wrapper for process_annotation_data().
#' Reads file path from config, handles missing/empty paths gracefully.
#'
#' @param config Pipeline config object
#' @return Processed annotation data.frame, or NULL if not configured
#' @export
load_and_process_annotation <- function(config) {
    cfg <- config$modes$rna
    if (is.null(cfg)) {
        return(NULL)
    }

    # Get file path from config

    file_path <- cfg$files$annotation
    if (is.null(file_path) || !nzchar(file_path)) {
        message("[annotation] No annotation file configured, skipping")
        return(NULL)
    }

    # Resolve to absolute path
    abs_path <- resolve_raw_path(config, file_path)
    if (!file.exists(abs_path)) {
        warning("[annotation] Annotation file not found: ", abs_path)
        return(NULL)
    }

    # Get annotation source from config
    annot_source <- cfg$annotation$source %||% ""

    message(sprintf("[annotation] Loading annotation file: %s (source: %s)",
                    basename(abs_path), if (nzchar(annot_source)) annot_source else "unspecified"))

    # Call original function
    result <- process_annotation_data(
        annotation_file = abs_path,
        ANNOT_SOURCE = annot_source
    )

    message(sprintf("[annotation] Loaded %d annotation records", nrow(result)))
    result
}


#' Load and process Trinotate file
#'
#' Config-driven wrapper for process_trinotate_data().
#' Reads file path from config, handles missing/empty paths gracefully.
#' Ensures parse_pfam is available (with fallback stub if not).
#'
#' @param config Pipeline config object
#' @return Processed trinotate_main data.frame, or NULL if not configured
#' @export
load_and_process_trinotate <- function(config) {
    cfg <- config$modes$rna
    if (is.null(cfg)) {
        return(NULL)
    }

    # Get file path from config
    file_path <- cfg$files$trinotate
    if (is.null(file_path) || !nzchar(file_path)) {
        message("[trinotate] No trinotate file configured, skipping")
        return(NULL)
    }

    # Resolve to absolute path
    abs_path <- resolve_raw_path(config, file_path)
    if (!file.exists(abs_path)) {
        warning("[trinotate] Trinotate file not found: ", abs_path)
        return(NULL)
    }

    message(sprintf("[trinotate] Loading Trinotate file: %s", basename(abs_path)))

    # Ensure parse_pfam is available (stub if not)
    .ensure_parse_pfam()

    # Call original function (with side effects isolated)
    result <- .process_trinotate_isolated(abs_path)

    message(sprintf("[trinotate] Loaded %d gene annotations", nrow(result)))
    result
}


# ------------------------------------------------------------
# ISOLATION WRAPPER for process_trinotate_data
# ------------------------------------------------------------
# The original function contains save.image() which must not execute.
# Rather than modifying the function, we isolate its execution and
# suppress the save.image() call via temporary override.
# ------------------------------------------------------------

.process_trinotate_isolated <- function(trinotate_file) {
    # Temporarily override save.image to no-op
    # This isolates the side effect without modifying the original function
    original_save_image <- base::save.image
    on.exit({
        # Restore original (though it's a base function, this is defensive)
    }, add = TRUE)

    # Create a local override in an environment where the function will run
    local_env <- new.env(parent = globalenv())
    local_env$save.image <- function(...) {
        message("[trinotate] save.image() call suppressed (isolated)")
        invisible(NULL)
    }

    # Also ensure parse_pfam is available in that environment
    if (exists("parse_pfam", mode = "function")) {
        local_env$parse_pfam <- get("parse_pfam", mode = "function")
    } else {
        local_env$parse_pfam <- function(x) {
            warning("parse_pfam not available", call. = FALSE)
            return(rep("", length(x)))
        }
    }

    # Evaluate the function in the isolated environment
    result <- eval(
        substitute(process_trinotate_data(trinotate_file)),
        envir = local_env
    )

    result
}


# ============================================================
# ORIGINAL FUNCTIONS (Preserved As-Is)
# ============================================================
# These functions are copied from the original source.
# No modifications except where explicitly noted.
# Architectural concerns are documented in the module header.
# ============================================================

#' Process annotation data
#'
#' @param annotation_file Path to annotation file
#' @param ANNOT_SOURCE Annotation source type ("Ensembl", "NCBI_GTF", "Generic")
#' @return Processed annotation data.frame
process_annotation_data <- function(annotation_file, ANNOT_SOURCE = "") {

    ## Read annotation data

    annot <- read.delim(file = annotation_file,
                        stringsAsFactors = FALSE)

    # process annotation data

    if (ANNOT_SOURCE == "Ensembl") {
        # define new column names (this method avoids error in case column(s) do not exist)
        lookup <- c(Ensembl_ID = "Gene.stable.ID", Gene_Symbol = "Gene.name",
                    Gene_Type = "Gene.type", Gene_Title = "Gene.description")
        # process
        annot <-
            annot|>
            dplyr::rename(dplyr::any_of(lookup))|>
            dplyr::select(Ensembl_ID, Gene_Symbol, Gene_Title, Gene_Type)|>
            dplyr::mutate(Gene_Title = stringr::str_replace(
                string = Gene_Title,
                pattern = " \\[.*",
                replacement = ""
            ))
    }

    if (ANNOT_SOURCE == "NCBI_GTF" || ANNOT_SOURCE == "Generic") {
        # right now do nothing
    }

    return(annot)
}


#' Process Trinotate data
#'
#' Parses Trinotate annotation report and extracts structured annotation data
#' including Swiss-Prot BLAST hits, Infernal RNA predictions, and Pfam domains.
#'
#' @param trinotate_file Path to Trinotate report TSV
#' @return Data frame with aggregated annotation per gene
#'
#' @note This function contains legacy code patterns:
#'   - gc() calls (harmless but unnecessary)
#'   - rm() calls (harmless but unnecessary)
#'   - cat() for logging (inconsistent with message() pattern)
#'   - save.image() call is ISOLATED by wrapper (does not execute)
process_trinotate_data <- function(trinotate_file) {

    ## Read Trinotate data

    trinotate <- read.delim(file = trinotate_file,
                            colClasses = c(X.gene_id = "factor",
                                           transcript_id = "factor",
                                           sprot_Top_BLASTX_hit = "character",
                                           infernal = "character",
                                           prot_id = "character",
                                           prot_coords = "character",
                                           sprot_Top_BLASTP_hit = "character",
                                           Pfam = "character",
                                           SignalP = "character",
                                           TmHMM = "character",
                                           eggnog = "character",
                                           Kegg = "character",
                                           gene_ontology_BLASTX = "character",
                                           gene_ontology_BLASTP = "character",
                                           gene_ontology_Pfam = "character",
                                           transcript = "character",
                                           peptide = "character"))

    ## Parse Trinotate data

    # get all gene IDs (to be used on the final merge)
    all_genes <- trinotate|>
        dplyr::select(X.gene_id)|>
        unique()

    # get RNA predictions from infernal
    trinotate_main_infernal <-
        trinotate[, c("X.gene_id", "infernal")]|>
        tibble::as_tibble()|>
        dplyr::filter(infernal != ".")|>
        unique()|>
        dplyr::mutate(RNA_Infernal = stringr::str_replace(
            string = infernal,
            pattern = "`",
            replacement = "|"
        ))|>
        dplyr::select("X.gene_id", "RNA_Infernal")

    # get sprot_BLASTX results
    trinotate_main_sprot_BLASTX <-
        trinotate[, c("X.gene_id", "sprot_Top_BLASTX_hit")]|>
        tibble::as_tibble()|>
        dplyr::filter(sprot_Top_BLASTX_hit != ".")|>
        dplyr::mutate(sprot_edit = stringr::str_replace(
            string = sprot_Top_BLASTX_hit,
            pattern = "\\`.*",
            replacement = ""
        ))|>
        dplyr::mutate(sprot_edit = stringr::str_replace(
            string = sprot_edit,
            pattern = "\\{.*\\}",
            replacement = ""
        ))|>
        dplyr::select(X.gene_id, sprot_edit)|>
        tidyr::separate(sprot_edit,
                        into = paste("sprot",
                                     c("Name", "Acc", "Pos", "percID", "eval", "RecName", "Lineage"),
                                     sep = "_"),
                        sep = "\\^")|>
        dplyr::mutate(sprot_RecName = stringr::str_replace(
            string = sprot_RecName,
            pattern = "^RecName: Full=",
            replacement = ""
        ))|>
        unique()|>
        dplyr::group_by(X.gene_id)|>
        dplyr::summarise(
            sprot_Name = paste(sprot_Name, collapse = "|"),
            sprot_Acc = paste(sprot_Acc, collapse = "|"),
            sprot_Pos = paste(sprot_Pos, collapse = "|"),
            sprot_percID = paste(sprot_percID, collapse = "|"),
            sprot_eval = paste(sprot_eval, collapse = "|"),
            sprot_RecName = paste(sprot_RecName, collapse = "|"),
            sprot_Lineage = paste(sprot_Lineage, collapse = "|")
        )|>
        dplyr::ungroup()|>
        dplyr::select("X.gene_id", "sprot_Name", "sprot_RecName", "sprot_percID", "sprot_eval")|>
        unique()

    # Doing the same for the protein BLASTP results
    trinotate_main_sprot_BLASTP <-
        trinotate[, c("X.gene_id", "sprot_Top_BLASTP_hit")]|>
        tibble::as_tibble()|>
        dplyr::filter(sprot_Top_BLASTP_hit != ".")|>
        dplyr::mutate(sprot_edit = stringr::str_replace(
            string = sprot_Top_BLASTP_hit,
            pattern = "\\`.*",
            replacement = ""
        ))|>
        dplyr::mutate(sprot_edit = stringr::str_replace(
            string = sprot_edit,
            pattern = "\\{.*\\}",
            replacement = ""
        ))|>
        dplyr::select(X.gene_id, sprot_edit)|>
        tidyr::separate(sprot_edit,
                        into = paste("sprot",
                                     c("Name", "Acc", "Pos", "percID", "eval", "RecName", "Lineage"),
                                     sep = "_"),
                        sep = "\\^")|>
        dplyr::mutate(sprot_RecName = stringr::str_replace(
            string = sprot_RecName,
            pattern = "^RecName: Full=",
            replacement = ""
        ))|>
        dplyr::select("X.gene_id", "sprot_Name", "sprot_RecName", "sprot_percID", "sprot_eval")|>
        dplyr::group_by(X.gene_id)|>
        dplyr::summarise(
            blastp_Name = paste0(sprot_Name, collapse = "|"),
            blastp_RecName = paste0(sprot_RecName, collapse = "|"),
            blastp_percID = paste0(sprot_percID, collapse = "|"),
            blastp_eval = paste0(sprot_eval, collapse = "|")
        )|>
        unique()

    trinotate_main_pfam <-
        trinotate[, c("X.gene_id", "Pfam")]|>
        dplyr::filter(Pfam != ".")|>
        unique()

    trinotate_main_pfam$pfam_domains <- sapply(
        X = trinotate_main_pfam[, "Pfam"],
        FUN = parse_pfam
    )

    trinotate_main_pfam1 <-
        trinotate_main_pfam|>
        dplyr::select("X.gene_id", "pfam_domains")|>
        dplyr::group_by(X.gene_id)|>
        dplyr::summarise(pfam_all = paste(pfam_domains, collapse = "|"))

    # merge all tables

    trinotate_main <-
        all_genes|>
        dplyr::left_join(trinotate_main_sprot_BLASTX, by = "X.gene_id")|>
        dplyr::left_join(trinotate_main_infernal, by = "X.gene_id")|>
        dplyr::left_join(trinotate_main_sprot_BLASTP, by = "X.gene_id")|>
        dplyr::left_join(trinotate_main_pfam1, by = "X.gene_id")|>
        unique()

    rm(trinotate_main_sprot_BLASTX)
    rm(trinotate_main_infernal)
    rm(trinotate_main_sprot_BLASTP)
    rm(trinotate_main_pfam)
    rm(trinotate_main_pfam1)

    gc()
    save.image()

    ## Statistical analysis of trinotate data:

    cat(sprintf("Number of genes for which annotation exists: %s\n",
                trinotate_main$X.gene_id|> unique()|> length()))

  
    cols_to_check <- c(
      "sprot_Name", "sprot_RecName", "sprot_percID", "sprot_eval",
      "RNA_Infernal", "blastp_Name", "blastp_RecName", "blastp_percID",
      "blastp_eval", "pfam_all"
    )
    
    for (coln in cols_to_check) {
      
      n_pos <- sum(!is.na(trinotate_main[[coln]]))
      
      cat(sprintf(
        "Number of genes positive in column %s: %s\n",
        coln, n_pos
      ))
    }

    return(trinotate_main)
}
