# R/domain/lipidomics/15_links_kegg_fetch.R
#
# Fetch gene + compound members for each KEGG pathway via KEGGREST and cache
# results on disk so reruns are offline. The cache lives under
#   <out_dir>/lipid_links/cache/kegg/<pathway_id>.tsv
# and is reused if `force_refetch = FALSE`.
#
# Output: long table with columns
#   kegg_pathway_id, member_type ("gene" | "compound"),
#   member_id, member_name


fetch_kegg_pathway_members <- function(pathway_ids, cache_dir,
                                        force_refetch = FALSE,
                                        organism = "cel") {
    pathway_ids <- unique(stats::na.omit(pathway_ids))
    if (length(pathway_ids) == 0) return(NULL)
    ensure_dir(cache_dir)

    out_list <- list()
    for (pid in pathway_ids) {
        cache_file <- file.path(cache_dir, paste0(pid, ".tsv"))
        if (!force_refetch && file.exists(cache_file)) {
            out_list[[pid]] <- utils::read.table(cache_file, header = TRUE,
                                                  sep = "\t",
                                                  stringsAsFactors = FALSE,
                                                  quote = "")
            next
        }
        message("[links/kegg] fetching ", pid)
        members <- tryCatch(.kegg_one_pathway(pid),
                            error = function(e) {
                                message("[links/kegg] ", pid,
                                        " failed: ", conditionMessage(e))
                                NULL
                            })
        if (is.null(members) || nrow(members) == 0) {
            members <- data.frame(kegg_pathway_id = character(0),
                                   member_type    = character(0),
                                   member_id      = character(0),
                                   member_name    = character(0),
                                   stringsAsFactors = FALSE)
        }
        utils::write.table(members, cache_file, sep = "\t", quote = FALSE,
                           row.names = FALSE)
        out_list[[pid]] <- members
    }
    members_all <- do.call(rbind, out_list)
    if (is.null(members_all) || nrow(members_all) == 0) return(members_all)

    # Annotate gene members with NCBI Entrez IDs (cel uses CELE_* clone ids
    # natively, so a separate keggConv step is required).
    gene_idx <- members_all$member_type == "gene"
    if (any(gene_idx)) {
        e_map <- fetch_kegg_gene_id_map(organism = organism,
                                         cache_dir = cache_dir)
        # member_id is "CELE_*" — KEGG uses prefix "cel:" in keggConv keys
        kegg_keys <- paste0(organism, ":", members_all$member_id[gene_idx])
        members_all$entrez_id <- NA_character_
        members_all$entrez_id[gene_idx] <- unname(e_map[kegg_keys])
    } else {
        members_all$entrez_id <- NA_character_
    }
    members_all
}


# One-time fetch + cache of KEGG cel-id -> NCBI Entrez gene-id map.
# Returns a named character vector keyed by "cel:CELE_*".
fetch_kegg_gene_id_map <- function(organism = "cel", cache_dir,
                                     force_refetch = FALSE) {
    ensure_dir(cache_dir)
    map_file <- file.path(cache_dir, paste0(organism, "_to_entrez.tsv"))
    if (!force_refetch && file.exists(map_file)) {
        tbl <- utils::read.table(map_file, header = TRUE, sep = "\t",
                                  stringsAsFactors = FALSE, quote = "")
        return(stats::setNames(tbl$entrez_id, tbl$kegg_id))
    }
    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop("KEGGREST not available")
    }
    message("[links/kegg] fetching ", organism, " -> Entrez gene-id map")
    raw <- KEGGREST::keggConv("ncbi-geneid", organism)
    # raw is a named char vec: names = "cel:CELE_*", values = "ncbi-geneid:1234"
    entrez_clean <- sub("^ncbi-geneid:", "", unname(raw))
    tbl <- data.frame(kegg_id = names(raw),
                       entrez_id = entrez_clean,
                       stringsAsFactors = FALSE)
    utils::write.table(tbl, map_file, sep = "\t", quote = FALSE,
                       row.names = FALSE)
    stats::setNames(tbl$entrez_id, tbl$kegg_id)
}


# Internal: fetch a single pathway via KEGGREST::keggGet.
# Genes: pid is "cel00564" -> entry returns $GENE = chr vec of alternating
#   "12345" then "gene-symbol;long-desc"
# Compounds: $COMPOUND is a named char vec (names = C-ids, values = names).
.kegg_one_pathway <- function(pid) {
    if (!requireNamespace("KEGGREST", quietly = TRUE)) {
        stop("KEGGREST not available")
    }
    res <- KEGGREST::keggGet(pid)[[1]]
    rows <- list()

    g <- res$GENE
    if (!is.null(g) && length(g) > 0) {
        # Pairs: odd indices are Entrez IDs, even are description strings.
        idx_id   <- seq(1, length(g), by = 2)
        idx_desc <- seq(2, length(g), by = 2)
        # Truncate if uneven (some KEGG entries have stray rows)
        n_pairs <- min(length(idx_id), length(idx_desc))
        if (n_pairs > 0) {
            gid <- g[idx_id[seq_len(n_pairs)]]
            gdesc <- g[idx_desc[seq_len(n_pairs)]]
            # First token before ';' is the gene symbol when present.
            # Some entries start directly with "[KO:...]" — fall back to
            # the clone ID in that case.
            gsym <- vapply(seq_along(gdesc), function(i) {
                first <- trimws(strsplit(gdesc[i], ";", fixed = TRUE)[[1]][1])
                if (is.na(first) || !nzchar(first) ||
                    startsWith(first, "[")) {
                    sub("^CELE_", "", gid[i])
                } else {
                    first
                }
            }, character(1))
            rows$genes <- data.frame(
                kegg_pathway_id = pid,
                member_type    = "gene",
                member_id      = as.character(gid),
                member_name    = gsym,
                stringsAsFactors = FALSE
            )
        }
    }

    cp <- res$COMPOUND
    if (!is.null(cp) && length(cp) > 0) {
        rows$compounds <- data.frame(
            kegg_pathway_id = pid,
            member_type    = "compound",
            member_id      = names(cp),
            member_name    = unname(cp),
            stringsAsFactors = FALSE
        )
    }
    if (length(rows) == 0) return(NULL)
    do.call(rbind, rows)
}
