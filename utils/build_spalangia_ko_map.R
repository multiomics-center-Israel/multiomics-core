#!/usr/bin/env Rscript
#' Build a feature -> KEGG Orthology (KO) map for Spalangia cameroni
#'
#' Spalangia cameroni has no KEGG organism code and no OrgDb, so KEGG maps can
#' only be rendered in KO space (`pathview(species = "ko")`). This script turns
#' the eggNOG-mapper annotation of the funannotate proteome into the lookup that
#' the multiomics pathview step needs: one row per (omics, feature_id, KO).
#'
#' The three ID spaces in this project do not line up, so each omics layer needs
#' its own transform onto the eggNOG `#query` key:
#'   - transcriptomics `evm.TU.*`  -> `evm.model.*` (gene -> model naming)
#'   - transcriptomics `BRK_g*`    -> `BRK_g*.t1`   (gene -> first transcript)
#'   - proteomics `<id>|<species>` -> `<id>`, taking the first `;`-separated
#'     member of a protein group that resolves
#'
#' Usage:
#'   Rscript utils/build_spalangia_ko_map.R
#'   Rscript utils/build_spalangia_ko_map.R --data-dir data/Elad_Chiel_June2026 \
#'       --output data/Elad_Chiel_June2026/multiomics/feature_to_ko.tsv
#'
#' Optional: pass --de-ids <rna_ids.txt>,<prot_ids.txt> to also report coverage
#' restricted to the features that were actually DE-tested.
#'
#' @examples
#' # From R:
#' source("utils/build_spalangia_ko_map.R")
#' ko_by_query <- read_eggnog_ko_map("eggNOG_annotations.funannotate.txt")

# ==============================================================================
# SETUP
# ==============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

# ==============================================================================
# MAIN FUNCTIONS
# ==============================================================================

#' Read the KEGG_ko column of an eggNOG-mapper annotation file
#'
#' The emapper output carries `##` banner lines and a header line that starts
#' with `#query`; both are handled here rather than by `read.delim`, which would
#' otherwise treat every `#` line as a comment and lose the column names.
#'
#' @param path Path to the eggNOG-mapper annotation TSV.
#' @return Named list, one element per query, each a character vector of KO ids
#'   (`K#####`, no `ko:` prefix). Queries with no KO get a zero-length element.
read_eggnog_ko_map <- function(path) {
    lines <- readLines(path, warn = FALSE)
    hdr_i <- grep("^#query", lines)[1]
    if (is.na(hdr_i)) {
        stop("No '#query' header line found in ", path,
             " - is this an eggNOG-mapper annotation file?")
    }

    header <- strsplit(sub("^#", "", lines[hdr_i]), "\t", fixed = TRUE)[[1]]
    body <- lines[seq.int(hdr_i + 1L, length(lines))]
    body <- body[nzchar(body) & !grepl("^##", body)]
    if (length(body) == 0) stop("No annotation rows after the header in ", path)

    fields <- strsplit(body, "\t", fixed = TRUE)
    query <- vapply(fields, function(x) x[1], character(1))
    ko_col <- match("KEGG_ko", header)
    if (is.na(ko_col)) stop("No 'KEGG_ko' column in ", path)
    ko_raw <- vapply(fields, function(x) if (length(x) >= ko_col) x[ko_col] else "-",
                     character(1))

    ko_by_query <- lapply(strsplit(ko_raw, ",", fixed = TRUE), function(x) {
        x <- sub("^ko:", "", trimws(x))
        unique(x[nzchar(x) & x != "-"])
    })
    names(ko_by_query) <- query
    ko_by_query
}


#' Map transcriptomics feature ids onto eggNOG query keys
#'
#' @param feature_ids Character vector of RNA feature ids (counts-matrix Geneid).
#' @return Character vector of candidate eggNOG query keys, same length as input.
rna_feature_to_eggnog_key <- function(feature_ids) {
    keys <- feature_ids
    is_evm <- grepl("^evm\\.TU\\.", keys)
    keys[is_evm] <- sub("^evm\\.TU\\.", "evm.model.", keys[is_evm])
    is_brk <- grepl("^BRK_g[0-9]+$", keys)
    keys[is_brk] <- paste0(keys[is_brk], ".t1")
    keys
}


#' Map a proteomics protein group onto eggNOG query keys
#'
#' A DIA-NN protein group can list several proteins separated by `;`, each
#' carrying a `|<species>` suffix. Group members are ranked by the group's own
#' order, so the first member that carries a KO is the representative.
#'
#' @param protein_group Single protein-group string.
#' @return Character vector of candidate eggNOG query keys, in group order.
protein_group_to_eggnog_keys <- function(protein_group) {
    parts <- strsplit(protein_group, ";", fixed = TRUE)[[1]]
    sub("\\|.*$", "", parts)
}


#' Resolve feature ids to KO ids through the eggNOG lookup
#'
#' @param feature_ids Character vector of feature ids.
#' @param key_fun Function mapping one feature id to candidate eggNOG keys.
#' @param ko_by_query Named list from \code{read_eggnog_ko_map}.
#' @param omics Omics label written into the `omics` column.
#' @return Data frame with columns `omics`, `feature_id`, `KO`; long format, one
#'   row per (feature, KO) pair. Features with no KO contribute no rows.
resolve_features_to_ko <- function(feature_ids, key_fun, ko_by_query, omics) {
    feature_ids <- unique(feature_ids[!is.na(feature_ids) & nzchar(feature_ids)])
    hits <- lapply(feature_ids, function(fid) {
        for (key in key_fun(fid)) {
            ko <- ko_by_query[[key]]
            if (!is.null(ko) && length(ko) > 0) return(ko)
        }
        character(0)
    })
    n <- lengths(hits)
    data.frame(
        omics = rep(omics, sum(n)),
        feature_id = rep(feature_ids, n),
        # as.character() keeps the column a character vector when nothing maps;
        # unlist() on an all-empty list would return NULL and drop the column.
        KO = as.character(unlist(hits, use.names = FALSE)),
        stringsAsFactors = FALSE
    )
}


#' Build the full feature -> KO map for this project
#'
#' @param data_dir Project data directory (contains `multiomics/` and `proteomics/`).
#' @param output Path of the TSV to write.
#' @return The written data frame (invisibly), with columns `omics`, `feature_id`, `KO`.
build_spalangia_ko_map <- function(data_dir, output) {
    egg_path <- file.path(data_dir, "proteomics", "eggNOG_annotations.funannotate.txt")
    rna_path <- file.path(data_dir, "multiomics", "rna_counts_8samples.tsv")
    pg_path <- file.path(data_dir, "proteomics", "report.pg_matrix.tsv")
    for (p in c(egg_path, rna_path, pg_path)) {
        if (!file.exists(p)) stop("Required input not found: ", p)
    }

    ko_by_query <- read_eggnog_ko_map(egg_path)
    n_annotated <- sum(lengths(ko_by_query) > 0)
    message("eggNOG queries with >= 1 KO: ", n_annotated, " / ", length(ko_by_query))

    # colClasses = "character" keeps ids verbatim; the count columns are unused here.
    rna_ids <- read.delim(rna_path, sep = "\t", header = TRUE, check.names = FALSE,
                          colClasses = "character")$Geneid
    pg_ids <- read.delim(pg_path, sep = "\t", header = TRUE, check.names = FALSE,
                         colClasses = "character")$Protein.Group

    rna_map <- resolve_features_to_ko(rna_ids, rna_feature_to_eggnog_key,
                                      ko_by_query, "transcriptomics")
    prot_map <- resolve_features_to_ko(pg_ids, protein_group_to_eggnog_keys,
                                       ko_by_query, "proteomics")

    ko_map <- rbind(rna_map, prot_map)
    ko_map <- ko_map[order(ko_map$omics, ko_map$feature_id, ko_map$KO), ]

    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    write.table(ko_map, output, sep = "\t", quote = FALSE, row.names = FALSE)

    summary_df <- data.frame(
        layer = c("eggNOG queries", "transcriptomics", "proteomics"),
        with_ko = c(n_annotated,
                    length(unique(rna_map$feature_id)),
                    length(unique(prot_map$feature_id))),
        total = c(length(ko_by_query), length(unique(rna_ids)), length(unique(pg_ids))),
        expected = c(8625L, 6886L, 3210L),
        stringsAsFactors = FALSE
    )
    summary_df$match <- summary_df$with_ko == summary_df$expected
    print(summary_df, row.names = FALSE)
    message("Wrote ", nrow(ko_map), " (feature, KO) rows -> ", output)
    if (!all(summary_df$match)) {
        warning("Coverage does not match the expected counts - check the ID transforms.")
    }
    invisible(ko_map)
}


#' Report KO coverage restricted to the features actually DE-tested
#'
#' @param ko_map Data frame from \code{build_spalangia_ko_map}.
#' @param rna_de_ids Character vector of DE-tested RNA feature ids.
#' @param prot_de_ids Character vector of DE-tested proteomics feature ids.
#' @return Data frame of the per-layer coverage (invisibly).
report_de_ko_coverage <- function(ko_map, rna_de_ids, prot_de_ids) {
    rna_ko <- ko_map[ko_map$omics == "transcriptomics" & ko_map$feature_id %in% rna_de_ids, ]
    prot_ko <- ko_map[ko_map$omics == "proteomics" & ko_map$feature_id %in% prot_de_ids, ]
    cov <- data.frame(
        layer = c("transcriptomics (DE-tested)", "proteomics (DE-tested)"),
        with_ko = c(length(unique(rna_ko$feature_id)), length(unique(prot_ko$feature_id))),
        total = c(length(unique(rna_de_ids)), length(unique(prot_de_ids))),
        stringsAsFactors = FALSE
    )
    print(cov, row.names = FALSE)
    message("Shared unique KOs (RNA n proteomics, DE-tested): ",
            length(intersect(unique(rna_ko$KO), unique(prot_ko$KO))))
    invisible(cov)
}

# ==============================================================================
# CLI
# ==============================================================================

# sys.nframe() == 0 only at the top level of an Rscript call, so `source()`ing
# this file (from a test, or from an R session) gets the functions without
# triggering the build.
if (!interactive() && sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)
    get_arg <- function(flag, default) {
        i <- match(flag, args)
        if (is.na(i) || i == length(args)) default else args[i + 1L]
    }
    data_dir <- get_arg("--data-dir", "data/Elad_Chiel_June2026")
    output <- get_arg("--output", file.path(data_dir, "multiomics", "feature_to_ko.tsv"))
    de_ids <- get_arg("--de-ids", NULL)

    ko_map <- build_spalangia_ko_map(data_dir, output)

    if (!is.null(de_ids)) {
        paths <- strsplit(de_ids, ",", fixed = TRUE)[[1]]
        if (length(paths) == 2 && all(file.exists(paths))) {
            report_de_ko_coverage(ko_map, readLines(paths[1]), readLines(paths[2]))
        } else {
            message("--de-ids expects '<rna_ids.txt>,<prot_ids.txt>' with both files present")
        }
    }
}
