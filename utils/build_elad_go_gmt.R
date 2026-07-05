#!/usr/bin/env Rscript
#' Build a GO gene-set (GMT) for the Elad Chiel proteomics run
#'
#' Turns the eggNOG-mapper annotation of the *Spalangia cameroni* proteome (and,
#' when available, a matching *Sodalis praecaptivus* annotation) into two files
#' the proteomics pipeline can use for GO enrichment on this non-model,
#' combined-genome DIA-NN dataset:
#'
#'   1. protein_to_GO.tsv        long mapping, one (protein_id, GO_id) pair per row
#'   2. GO_spalangia_sodalis.gmt tab-separated GMT: GO_id <tab> "name [aspect]" <tab> members...
#'
#' Member IDs carry the species suffix (e.g. "BRK_g10146.t1|Spalangia_cameroni")
#' so they match the pipeline's DE FeatureID, which equals the raw DIA-NN
#' Protein.Group (the empty Genes column means no symbol remap happens).
#'
#' Usage (from the repo root, inside the pipeline env):
#'   conda activate mofa2_env
#'   Rscript utils/build_elad_go_gmt.R
#'
#' Re-run after dropping a Sodalis eggNOG annotation at
#' data/Elad_Chiel_June2026/proteomics/Sodalis_eggnog_annotations.tsv to fold the
#' bacterial proteins into the same GMT.

suppressWarnings(suppressMessages({
    library(GO.db)
    library(AnnotationDbi)
}))

`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------------------------
# Locate the repo and source the shared GMT helpers (write_gmt/validate_gmt/
# filter_gmt_by_size in R/core/13_gmt_utils.R; read_gmt in R/core/09_enrichment.R)
# ------------------------------------------------------------------------------
script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) {
        args <- commandArgs(trailingOnly = FALSE)
        file_arg <- grep("^--file=", args, value = TRUE)
        if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg))) else "."
    }
)
repo_root <- dirname(script_dir)
source(file.path(repo_root, "R", "core", "13_gmt_utils.R"))
source(file.path(repo_root, "R", "core", "09_enrichment.R"))

# ------------------------------------------------------------------------------
# Inputs / outputs (relative to repo root for portability)
# ------------------------------------------------------------------------------
data_dir     <- file.path(repo_root, "data", "Elad_Chiel_June2026", "proteomics")
wasp_eggnog  <- file.path(data_dir, "eggNOG_annotations.funannotate.txt")
sodalis_gaf  <- file.path(data_dir, "Sodalis_GO.gaf")                   # preferred, NCBI RefSeq GAF
sodalis_file <- file.path(data_dir, "Sodalis_eggnog_annotations.tsv")  # fallback (eggNOG / 2-col)
pg_matrix    <- file.path(data_dir, "report.pg_matrix.tsv")
out_map      <- file.path(data_dir, "protein_to_GO.tsv")
out_gmt      <- file.path(data_dir, "GO_spalangia_sodalis.gmt")

# Runtime enrichment size window, mirrored here only for the validation report.
MIN_SIZE <- 10
MAX_SIZE <- 500

# ------------------------------------------------------------------------------
# Parsers -> long (protein_id, GO_id) data frame
# ------------------------------------------------------------------------------

#' Extract (protein_id, GO_id) pairs from an eggNOG-mapper annotation table.
#'
#' @param path Path to an eggNOG-mapper output (## comment lines tolerated; the
#'   header row starts with "#query"; GO terms live in the comma-separated "GOs"
#'   column and "-" means none).
#' @param suffix Species suffix appended to each query ID so members match the
#'   DIA-NN Protein.Group (e.g. "|Spalangia_cameroni").
#' @return A data.frame with columns protein_id, GO_id (long form), or NULL if
#'   the file has no GO annotations.
parse_eggnog_go <- function(path, suffix) {
    raw <- readLines(path, warn = FALSE)
    raw <- raw[!startsWith(raw, "##")]               # drop eggNOG stat/comment lines
    hdr_i <- grep("^#?query\t", raw)[1]
    if (is.na(hdr_i)) stop(sprintf("No '#query' header found in %s", path))
    header <- strsplit(sub("^#", "", raw[hdr_i]), "\t", fixed = TRUE)[[1]]
    body <- raw[(hdr_i + 1):length(raw)]
    body <- body[nzchar(body)]
    df <- read.delim(text = body, header = FALSE, quote = "", comment.char = "",
                     stringsAsFactors = FALSE, col.names = header)

    q <- df[["query"]]
    gos <- df[["GOs"]]
    keep <- !is.na(gos) & gos != "-" & nzchar(gos)
    if (!any(keep)) return(NULL)

    terms_list <- strsplit(gos[keep], ",", fixed = TRUE)
    data.frame(
        protein_id = rep(paste0(q[keep], suffix), lengths(terms_list)),
        GO_id      = unlist(terms_list, use.names = FALSE),
        stringsAsFactors = FALSE
    )
}

#' Extract (protein_id, GO_id) pairs from a GAF 2.x file (NCBI RefSeq style).
#'
#' Uses the protein accession (column 2) and GO ID (column 5); skips "!" comment
#' lines and NOT-qualified (negated) annotations. The species suffix is appended
#' so IDs match the DIA-NN Protein.Group (WP_ accessions here have none).
#'
#' @param path Path to the GAF file.
#' @param suffix Species suffix, e.g. "|Sodalis_praecaptivus".
#' @return A long data.frame (protein_id, GO_id), or NULL if empty.
parse_gaf_go <- function(path, suffix) {
    raw <- readLines(path, warn = FALSE)
    raw <- raw[!startsWith(raw, "!") & nzchar(raw)]
    if (!length(raw)) return(NULL)
    parts <- strsplit(raw, "\t", fixed = TRUE)
    acc  <- vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_, "")
    qual <- vapply(parts, function(p) if (length(p) >= 4) p[4] else "",           "")
    go   <- vapply(parts, function(p) if (length(p) >= 5) p[5] else NA_character_, "")
    keep <- !grepl("NOT", qual, fixed = TRUE) & !is.na(acc) & !is.na(go) & nzchar(go)
    if (!any(keep)) return(NULL)
    data.frame(
        protein_id = paste0(acc[keep], suffix),
        GO_id      = go[keep],
        stringsAsFactors = FALSE
    )
}

#' Extract (protein_id, GO_id) pairs from a Sodalis annotation file.
#'
#' Accepts either an eggNOG-mapper table (same shape as the wasp file) or a plain
#' two-column mapping (protein_id, GO_id) e.g. reshaped from UniProt-GOA. The
#' species suffix is appended only when the IDs do not already carry a "|".
#'
#' @param path Path to the Sodalis annotation file.
#' @return A long data.frame (protein_id, GO_id), or NULL.
parse_sodalis_go <- function(path) {
    suffix <- "|Sodalis_praecaptivus"
    head_line <- readLines(path, n = 5, warn = FALSE)
    if (any(grepl("^#?query\t", head_line))) {
        return(parse_eggnog_go(path, suffix))
    }
    # Plain two-column table: protein_id <tab> GO_id (one pair per row, or GO
    # terms comma-separated in the second column).
    df <- read.delim(path, header = TRUE, quote = "", comment.char = "",
                     stringsAsFactors = FALSE)
    if (ncol(df) < 2) stop(sprintf("Sodalis file %s needs >=2 columns", path))
    pid <- df[[1]]; gos <- df[[2]]
    keep <- !is.na(gos) & gos != "-" & nzchar(gos)
    if (!any(keep)) return(NULL)
    terms_list <- strsplit(gos[keep], ",", fixed = TRUE)
    pid <- ifelse(grepl("|", pid[keep], fixed = TRUE), pid[keep], paste0(pid[keep], suffix))
    data.frame(
        protein_id = rep(pid, lengths(terms_list)),
        GO_id      = trimws(unlist(terms_list, use.names = FALSE)),
        stringsAsFactors = FALSE
    )
}

# ------------------------------------------------------------------------------
# 1. Build the long protein -> GO mapping
# ------------------------------------------------------------------------------
message("== Building protein -> GO mapping ==")
long <- parse_eggnog_go(wasp_eggnog, "|Spalangia_cameroni")
message(sprintf("  Spalangia: %d proteins with GO, %d pairs",
                length(unique(long$protein_id)), nrow(long)))

# Prefer the NCBI RefSeq GAF (WP_ accessions match the DIA-NN Sodalis IDs
# directly); fall back to an eggNOG / 2-column table if only that is present.
sod <- NULL; sod_src <- NULL
if (file.exists(sodalis_gaf)) {
    sod <- parse_gaf_go(sodalis_gaf, "|Sodalis_praecaptivus"); sod_src <- basename(sodalis_gaf)
} else if (file.exists(sodalis_file)) {
    sod <- parse_sodalis_go(sodalis_file); sod_src <- basename(sodalis_file)
}
if (!is.null(sod)) {
    message(sprintf("  Sodalis:   %d proteins with GO, %d pairs (from %s)",
                    length(unique(sod$protein_id)), nrow(sod), sod_src))
    long <- rbind(long, sod)
} else {
    message("  Sodalis:   no annotation found -> building Spalangia-only.")
    message("             (Add Sodalis_GO.gaf or Sodalis_eggnog_annotations.tsv to the ",
            "data folder and re-run to include the bacterial proteins.)")
}

long <- unique(long)

# Restrict to proteins actually measured in this run. Enrichment only ever
# considers genes in the universe (the DE FeatureID / Protein.Group set), so
# annotations for unmeasured proteins are inert. Restricting here keeps the
# mapping and GMT specific to this experiment and their size sane -- the eggNOG
# GOs are propagated to ancestors, so the unrestricted pair count is ~1.4M.
feature_ids <- read.delim(pg_matrix, header = TRUE, quote = "", comment.char = "",
                          check.names = FALSE)[[1]]
n_all <- nrow(long)
long <- long[long$protein_id %in% feature_ids, ]
message(sprintf("  Restricted to measured features: %d -> %d pairs (%d of %d proteins matched)",
                n_all, nrow(long), length(unique(long$protein_id)), length(feature_ids)))

long <- long[order(long$protein_id, long$GO_id), ]
write.table(long, out_map, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("  Wrote %s (%d unique pairs)", basename(out_map), nrow(long)))

# ------------------------------------------------------------------------------
# 2. Invert to gene sets and name the terms via GO.db
# ------------------------------------------------------------------------------
message("== Building GO gene sets ==")
gene_sets <- lapply(split(long$protein_id, long$GO_id), unique)
go_ids <- names(gene_sets)

go_keys <- keys(GO.db, "GOID")
known <- go_ids %in% go_keys
ann <- suppressMessages(AnnotationDbi::select(
    GO.db, keys = go_ids[known], columns = c("TERM", "ONTOLOGY"), keytype = "GOID"))
term_map <- setNames(ann$TERM, ann$GOID)
ont_map  <- setNames(ann$ONTOLOGY, ann$GOID)

descriptions <- ifelse(
    go_ids %in% names(term_map) & !is.na(term_map[go_ids]),
    sprintf("%s [%s]", term_map[go_ids], ont_map[go_ids]),
    go_ids  # fall back to the bare GO id for obsolete / secondary terms
)
names(descriptions) <- go_ids
attr(gene_sets, "descriptions") <- descriptions

n_unnamed <- sum(!(go_ids %in% names(term_map)) | is.na(term_map[go_ids]))
message(sprintf("  %d GO terms (%d without a current GO.db name -> kept as GO id)",
                length(gene_sets), n_unnamed))
sizes <- lengths(gene_sets)
message(sprintf("  Set sizes: min %d, median %d, max %d; %d within the %d-%d window",
                min(sizes), as.integer(median(sizes)), max(sizes),
                sum(sizes >= MIN_SIZE & sizes <= MAX_SIZE), MIN_SIZE, MAX_SIZE))

# ------------------------------------------------------------------------------
# 3. Trim to the enrichment size window and write the GMT
# ------------------------------------------------------------------------------
# Members are already restricted to the measured universe, so these set sizes are
# final: sets outside [MIN_SIZE, MAX_SIZE] are exactly the ones the pipeline would
# drop at runtime (fgsea min/maxSize; ORA gene-set size filter). Trimming here is
# result-neutral for this config and keeps the file small. NOTE: if the config's
# pathway.min_size/max_size change, re-run this builder with matching bounds.
kept <- filter_gmt_by_size(gene_sets, min_size = MIN_SIZE, max_size = MAX_SIZE)
message(sprintf("  Trimmed to size %d-%d: %d -> %d sets written (dropped %d out-of-window)",
                MIN_SIZE, MAX_SIZE, length(gene_sets), length(kept),
                length(gene_sets) - length(kept)))
write_gmt(kept, out_gmt, verbose = TRUE)

# ------------------------------------------------------------------------------
# 4. Coverage report against the actual DIA-NN features + round-trip check
# ------------------------------------------------------------------------------
message("== Validating against pg_matrix features ==")
invisible(validate_gmt(gene_sets, feature_ids,
                       min_coverage = 0.05,
                       min_pathway_size = MIN_SIZE,
                       max_pathway_size = MAX_SIZE,
                       verbose = TRUE))

rt <- read_gmt(out_gmt)
message(sprintf("Round-trip: read_gmt() loaded %d sets; first members carry the ",
                length(rt)), appendLF = FALSE)
message(sprintf("species suffix: %s",
                all(grepl("\\|", head(rt[[1]], 3)))))

message("Done.")
