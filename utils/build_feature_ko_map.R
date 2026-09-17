#!/usr/bin/env Rscript
#' Build a feature -> KEGG Orthology (KO) map from an eggNOG-mapper annotation
#'
#' An organism with no KEGG organism code and no OrgDb can still be placed onto
#' KEGG's reference maps, because those maps label their gene boxes with KO ids
#' rather than with any one species' genes. What the multiomics pathview step
#' needs in order to do that is a lookup from this project's feature ids to KO
#' ids, and this script builds one from an eggNOG-mapper annotation of the
#' proteome: one row per (omics, feature_id, KO).
#'
#' The output is what `modes.multiomics.enrichment.pathview.ko_map` points at.
#'
#' Feature ids rarely match the annotation's query keys verbatim, so each layer
#' gets its own transform onto that key space (see the two `*_to_eggnog_key*`
#' functions). Both transforms are conventions of the tools that produced the
#' ids, not of any one project; check them against your own annotation before
#' trusting the coverage.
#'
#' Usage:
#'   Rscript utils/build_feature_ko_map.R \
#'     --annotation <emapper.annotations.tsv> \
#'     --rna <rna_counts.tsv> \
#'     --proteomics <protein_groups.tsv> \
#'     --output <feature_to_ko.tsv>
#'
#' At least one of --rna / --proteomics is required. Use --rna-column and
#' --proteomics-column when the feature ids do not sit in the default column
#' (`Geneid`, as featureCounts writes it, and `Protein.Group`, as DIA-NN does).
#'
#' Optional: --de-ids <rna_ids.txt>,<prot_ids.txt> also reports coverage
#' restricted to the features that were actually DE-tested, which is the number
#' that matters for a pathway map.
#'
#' @examples
#' # From R:
#' source("utils/build_feature_ko_map.R")
#' ko_by_query <- read_eggnog_ko_map("emapper.annotations.tsv")

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


#' Map transcriptomics feature ids onto candidate eggNOG query keys
#'
#' Two rewrites, both matching how the annotation tools name things rather than
#' anything project-specific:
#' \itemize{
#'   \item EVM/funannotate gene ids (`evm.TU.*`) become the model ids
#'     (`evm.model.*`) that the proteome, and therefore the annotation, is keyed
#'     on. This one is a naming rule and holds generally.
#'   \item A bare BRAKER gene id (`BRK_g123`) gains a `.t1` suffix to reach its
#'     first transcript. This is a matching heuristic, not a guarantee: BRAKER
#'     numbers transcripts per gene and `.t1` is only conventionally the one the
#'     proteome carries. A gene whose annotated protein comes from another
#'     isoform will simply not match, and shows up as missing coverage rather
#'     than as a wrong KO.
#' }
#'
#' @param feature_ids Character vector of RNA feature ids.
#' @return Character vector of candidate eggNOG query keys, same length as input.
rna_feature_to_eggnog_key <- function(feature_ids) {
    keys <- feature_ids
    is_evm <- grepl("^evm\\.TU\\.", keys)
    keys[is_evm] <- sub("^evm\\.TU\\.", "evm.model.", keys[is_evm])
    is_brk <- grepl("^BRK_g[0-9]+$", keys)
    keys[is_brk] <- paste0(keys[is_brk], ".t1")
    keys
}


#' Map a proteomics protein group onto candidate eggNOG query keys
#'
#' A protein group can list several proteins separated by `;`, each carrying a
#' `|<suffix>` tag. Group members are ranked by the group's own order, so the
#' first member that carries a KO is taken as the representative.
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


#' Read one column of feature ids from a delimited table
#'
#' @param path Path to a TSV whose rows are features.
#' @param column Name of the column holding the feature ids.
#' @return Character vector of ids.
read_feature_ids <- function(path, column) {
    # colClasses = "character" keeps ids verbatim -- without it a numeric-looking
    # id loses its leading zeros. The measurement columns are unused here.
    tbl <- utils::read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
                             colClasses = "character")
    if (!column %in% names(tbl)) {
        stop("Column '", column, "' not found in ", path,
             ". Available: ", paste(utils::head(names(tbl), 10), collapse = ", "),
             if (ncol(tbl) > 10) ", ..." else "")
    }
    tbl[[column]]
}


#' Build the feature -> KO map
#'
#' @param annotation Path to the eggNOG-mapper annotation TSV.
#' @param output Path of the TSV to write.
#' @param rna Path to a transcriptomics table, or NULL to skip that layer.
#' @param proteomics Path to a proteomics table, or NULL to skip that layer.
#' @param rna_column Column of \code{rna} holding the feature ids.
#' @param proteomics_column Column of \code{proteomics} holding the protein groups.
#' @return The written data frame (invisibly), with columns `omics`,
#'   `feature_id`, `KO`.
build_feature_ko_map <- function(annotation, output, rna = NULL, proteomics = NULL,
                                 rna_column = "Geneid",
                                 proteomics_column = "Protein.Group") {
    if (is.null(rna) && is.null(proteomics)) {
        stop("Nothing to map: pass at least one of --rna / --proteomics.")
    }
    for (p in c(annotation, rna, proteomics)) {
        if (!file.exists(p)) stop("Input not found: ", p)
    }

    ko_by_query <- read_eggnog_ko_map(annotation)
    message("eggNOG queries with >= 1 KO: ", sum(lengths(ko_by_query) > 0),
            " / ", length(ko_by_query))

    parts <- list()
    if (!is.null(rna)) {
        parts$transcriptomics <- resolve_features_to_ko(
            read_feature_ids(rna, rna_column),
            rna_feature_to_eggnog_key, ko_by_query, "transcriptomics")
    }
    if (!is.null(proteomics)) {
        parts$proteomics <- resolve_features_to_ko(
            read_feature_ids(proteomics, proteomics_column),
            protein_group_to_eggnog_keys, ko_by_query, "proteomics")
    }

    ko_map <- do.call(rbind, unname(parts))
    ko_map <- ko_map[order(ko_map$omics, ko_map$feature_id, ko_map$KO), ]

    if (nrow(ko_map) == 0) {
        warning("No feature resolved to a KO. The id transforms in this script ",
                "may not match how your annotation names its queries - compare a ",
                "few feature ids against the '#query' column before rerunning.")
    }

    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(ko_map, output, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Wrote ", nrow(ko_map), " (feature, KO) rows -> ", output)
    invisible(ko_map)
}


#' Report KO coverage restricted to the features actually DE-tested
#'
#' Coverage over the whole annotation flatters itself: what decides how much of
#' a pathway map gets coloured is coverage over the features that reached a DE
#' test.
#'
#' @param ko_map Data frame from \code{build_feature_ko_map}.
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
    message("Shared unique KOs (RNA and proteomics, DE-tested): ",
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
    get_arg <- function(flag, default = NULL) {
        i <- match(flag, args)
        if (is.na(i) || i == length(args)) default else args[i + 1L]
    }
    required <- function(flag) {
        v <- get_arg(flag)
        if (is.null(v)) {
            stop(flag, " is required. See the usage block at the top of this file.",
                 call. = FALSE)
        }
        v
    }

    ko_map <- build_feature_ko_map(
        annotation        = required("--annotation"),
        output            = required("--output"),
        rna               = get_arg("--rna"),
        proteomics        = get_arg("--proteomics"),
        rna_column        = get_arg("--rna-column", "Geneid"),
        proteomics_column = get_arg("--proteomics-column", "Protein.Group")
    )

    de_ids <- get_arg("--de-ids")
    if (!is.null(de_ids)) {
        paths <- strsplit(de_ids, ",", fixed = TRUE)[[1]]
        if (length(paths) == 2 && all(file.exists(paths))) {
            report_de_ko_coverage(ko_map, readLines(paths[1]), readLines(paths[2]))
        } else {
            message("--de-ids expects '<rna_ids.txt>,<prot_ids.txt>' with both files present")
        }
    }
}
