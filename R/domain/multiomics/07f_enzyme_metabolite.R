# =============================================================================
# Enzyme-metabolite pairs: KEGG access and parsing
# =============================================================================
#
# What an enzyme in the proteomics layer could act on among the metabolites the
# metabolomics layer measured. This file reaches KEGG and parses what comes
# back; 07g_enzyme_metabolite_pairs.R assembles the table from it.
#
# Every fetch here is a bulk endpoint cached as an RDS, so a rerun makes no
# network call at all, and a run with no network returns NULL rather than
# failing: the pairs are an annotation lookup, not a result the pipeline needs
# in order to finish.


#' Metabolites too common to make a pair informative
#'
#' A currency metabolite takes part in thousands of reactions, so pairing it
#' with an enzyme says nothing about that enzyme: ATP is a substrate of most
#' kinases in the table and of most of the rest too. Dropping them keeps the
#' pairs to compounds the enzyme is specific about. Named for readability in
#' the code that filters on them; only the KEGG ids are used.
#'
#' @format Named character vector: names are readable, values are KEGG compound
#'   ids.
ENZYME_METABOLITE_CURRENCY_COMPOUNDS <- c(
    water = "C00001", ATP = "C00002", ADP = "C00008", AMP = "C00020",
    NAD_plus = "C00003", NADH = "C00004", NADP_plus = "C00006", NADPH = "C00005",
    oxygen = "C00007", CO2 = "C00011", orthophosphate = "C00009",
    diphosphate = "C00013", proton = "C00080", CoA = "C00010",
    FAD = "C00016", FADH2 = "C01352", ammonia = "C00014"
)


#' Settings for the enzyme-metabolite pair table
#'
#' Defaults match \code{validate_multiomics_config()}; they are repeated here so
#' the builder can be called with a bare config in a test without the validator
#' having run first. The table is opt-in: it reaches KEGG, so it runs only when
#' `enabled` is TRUE. A flag that is not literally TRUE reads as FALSE, so a
#' string such as "no" can never switch anything on; the validator rejects
#' such values outright.
#'
#' @param config Full config object.
#' @return List with \code{enabled}, \code{enzyme_hits_only},
#'   \code{require_shared_pathway}, \code{drop_currency_metabolites} and
#'   \code{list_unpaired_enzymes}.
#' @examples
#' enzyme_metabolite_config(list())$enabled   # FALSE
enzyme_metabolite_config <- function(config) {
    cfg <- ((config$modes$multiomics$enrichment %||% list())$enzyme_metabolite) %||% list()
    flag <- function(x, default) if (is.null(x)) default else isTRUE(x)
    list(
        enabled = flag(cfg$enabled, FALSE),
        enzyme_hits_only = flag(cfg$enzyme_hits_only, TRUE),
        require_shared_pathway = flag(cfg$require_shared_pathway, TRUE),
        drop_currency_metabolites = flag(cfg$drop_currency_metabolites, TRUE),
        list_unpaired_enzymes = flag(cfg$list_unpaired_enzymes, TRUE)
    )
}


#' One KEGG link table, fetched in bulk and cached
#'
#' Wraps the \code{/link/<target>/<source>} endpoint, which returns every
#' association in one response -- one call for a whole organism's gene-to-enzyme
#' map, rather than one call per gene. The response is cached as an RDS under
#' \code{cache_dir}, so the only run that reaches the network is the first.
#'
#' @param target KEGG database to link to, e.g. "enzyme", "reaction",
#'   "compound", "pathway".
#' @param source KEGG database or organism code to link from, e.g. "enzyme",
#'   "reaction", or an organism code such as "hsa".
#' @param cache_dir Directory to cache into; NULL skips caching.
#' @return Data frame with columns \code{from} and \code{to}, prefixes stripped
#'   ("cpd:C00001" becomes "C00001"), or NULL when KEGG cannot be reached.
#' @examples
#' # kegg_link_table("enzyme", "hsa", cache_dir = tempdir())
kegg_link_table <- function(target, source, cache_dir = NULL) {
    stem <- paste0("kegg_link_", gsub("[^A-Za-z0-9]+", "_", target), "_",
                   gsub("[^A-Za-z0-9]+", "_", source), ".rds")
    cache_file <- if (!is.null(cache_dir)) file.path(cache_dir, stem) else NULL
    if (!is.null(cache_file) && file.exists(cache_file)) {
        cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
        if (is.data.frame(cached)) return(cached)
    }

    lines <- .kegg_rest_lines(paste0("https://rest.kegg.jp/link/", target, "/", source),
                              what = paste("link", source, "->", target))
    if (is.null(lines) || length(lines) == 0) return(NULL)

    parts <- strsplit(lines, "\t")
    keep <- lengths(parts) >= 2
    if (!any(keep)) return(NULL)
    out <- data.frame(
        from = .strip_kegg_prefix(vapply(parts[keep], `[`, character(1), 1)),
        to   = .strip_kegg_prefix(vapply(parts[keep], `[`, character(1), 2)),
        stringsAsFactors = FALSE
    )

    if (!is.null(cache_file)) {
        dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
        tryCatch(saveRDS(out, cache_file), error = function(e) NULL)
    }
    out
}


#' Read one KEGG REST response, retrying a failed request
#'
#' The link endpoints return large responses -- reaction to compound is hundreds
#' of thousands of lines -- and a single timeout used to be the end of it: the
#' table that needed it was skipped for the whole run, silently, while the other
#' link tables sat cached beside it. Each request is therefore tried a few
#' times, with a short pause between attempts, and says so when it gives up.
#'
#' @param url Full KEGG REST URL.
#' @param what Short description for the messages.
#' @param attempts How many times to try.
#' @param pause Seconds to wait after a failed attempt; doubles each time.
#' @return Character vector of lines, or NULL when every attempt failed.
#' @keywords internal
.kegg_rest_lines <- function(url, what, attempts = 3L, pause = 2) {
    for (i in seq_len(attempts)) {
        lines <- tryCatch(readLines(url, warn = FALSE), error = function(e) {
            message("    KEGG ", what, " attempt ", i, "/", attempts,
                    " failed: ", conditionMessage(e))
            NULL
        })
        if (!is.null(lines) && length(lines) > 0) return(lines)
        if (i < attempts) {
            Sys.sleep(pause)
            pause <- pause * 2
        }
    }
    message("    KEGG ", what, " unavailable after ", attempts, " attempts")
    NULL
}


#' KEGG's own name for each compound
#'
#' The metabolomics layer names a feature however its table did -- often an
#' m/z and retention time, which says nothing to a reader. KEGG carries a name
#' for every compound it knows, so a mapped metabolite can be shown by that
#' name as well as by the id it was measured under. One bulk call, cached like
#' the link tables.
#'
#' KEGG lists synonyms separated by "; "; the first is taken, which is the one
#' KEGG itself leads with.
#'
#' @param cache_dir Directory to cache into; NULL skips caching.
#' @return Named character vector, compound id to name, or NULL when KEGG
#'   cannot be reached.
#' @examples
#' # kegg_compound_names(tempdir())["C00026"]   # "2-Oxoglutarate"
kegg_compound_names <- function(cache_dir = NULL) {
    cache_file <- if (!is.null(cache_dir)) {
        file.path(cache_dir, "kegg_compound_names.rds")
    } else NULL
    if (!is.null(cache_file) && file.exists(cache_file)) {
        cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
        if (is.character(cached)) return(cached)
    }

    lines <- .kegg_rest_lines("https://rest.kegg.jp/list/compound",
                              what = "compound names")
    if (is.null(lines) || length(lines) == 0) return(NULL)

    parts <- strsplit(lines, "\t")
    keep <- lengths(parts) >= 2
    if (!any(keep)) return(NULL)
    ids <- .strip_kegg_prefix(vapply(parts[keep], `[`, character(1), 1))
    names_first <- trimws(vapply(strsplit(
        vapply(parts[keep], `[`, character(1), 2), ";", fixed = TRUE),
        `[`, character(1), 1))
    out <- stats::setNames(names_first, ids)

    if (!is.null(cache_file)) {
        dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
        tryCatch(saveRDS(out, cache_file), error = function(e) NULL)
    }
    out
}


#' Drop the database prefix KEGG puts on every identifier
#'
#' "cpd:C00001", "ec:1.1.1.1" and "path:map00010" all carry a prefix that no
#' other table in this pipeline uses. An organism-prefixed gene id ("hsa:1234")
#' keeps its prefix, because that prefix is part of the gene id.
#'
#' @param ids Character vector of KEGG identifiers.
#' @return The identifiers without their leading database prefix.
#' @keywords internal
.strip_kegg_prefix <- function(ids) {
    sub("^(cpd|ec|rn|path|gl|dr|br):", "", as.character(ids))
}


#' Compare KEGG gene ids that differ only by their organism prefix
#'
#' The two sides of the gene-to-EC join spell one gene differently.
#' \code{convert_entrez_to_kegg()} returns what \code{bitr_kegg()} gives it,
#' which for an organism whose KEGG gene ids are NCBI gene ids is the bare
#' number ("29740"), while \code{/link/enzyme/<org>} returns it prefixed
#' ("rno:29740"). Merging the two as they come matches nothing at all, and the
#' table then reports that no protein carries an EC number.
#'
#' @param ids Character vector of KEGG gene ids, prefixed or not.
#' @param kegg_org Organism code, e.g. "rno".
#' @return The ids without a leading \code{<org>:}.
#' @examples
#' kegg_gene_key(c("rno:29740", "29740"), "rno")   # both "29740"
kegg_gene_key <- function(ids, kegg_org) {
    ids <- as.character(ids)
    if (is.null(kegg_org) || !nzchar(kegg_org)) return(ids)
    sub(paste0("^", kegg_org, ":"), "", ids)
}


#' Substrates and products of one KEGG reaction equation
#'
#' A KEGG equation reads "C00001 + 2 C00002 <=> C00003". Everything left of the
#' arrow is a substrate and everything right of it a product, as KEGG writes
#' the reaction; most reactions are reversible in vivo, which is why the table
#' calls this a role rather than a direction.
#'
#' Stoichiometric coefficients, "(n)" polymer notation and glycan ids are
#' ignored: only compound ids are extracted.
#'
#' @param equation One equation string.
#' @return Data frame with \code{compound} and \code{role} ("substrate" or
#'   "product"); zero rows when the string holds no compound id. A compound on
#'   both sides appears twice.
#' @examples
#' parse_kegg_equation("C00001 + 2 C00002 <=> C00003")
parse_kegg_equation <- function(equation) {
    empty <- data.frame(compound = character(0), role = character(0),
                        stringsAsFactors = FALSE)
    if (length(equation) != 1 || is.na(equation) || !nzchar(equation)) return(empty)

    sides <- strsplit(equation, "<=>", fixed = TRUE)[[1]]
    if (length(sides) != 2) return(empty)

    ids <- function(side) unique(regmatches(side, gregexpr("C[0-9]{5}", side))[[1]])
    subs <- ids(sides[1])
    prods <- ids(sides[2])
    if (length(subs) == 0 && length(prods) == 0) return(empty)

    data.frame(
        compound = c(subs, prods),
        role = c(rep("substrate", length(subs)), rep("product", length(prods))),
        stringsAsFactors = FALSE
    )
}


#' Substrate and product roles for a set of KEGG reactions
#'
#' The only endpoint here that is not a single bulk call: \code{/get/} returns
#' at most ten entries at a time, so this is one call per ten reactions. It is
#' therefore called last, on the reactions that survived every filter, rather
#' than on every reaction an enzyme set touches.
#'
#' Results are cached per run directory and merged with what is already there,
#' so a rerun fetches only reactions it has not seen.
#'
#' @param reaction_ids Character vector of KEGG reaction ids ("R00623").
#' @param cache_dir Directory to cache into; NULL skips caching.
#' @return Data frame with \code{reaction}, \code{compound} and \code{role}, or
#'   NULL when nothing could be fetched.
fetch_reaction_roles <- function(reaction_ids, cache_dir = NULL) {
    reaction_ids <- unique(reaction_ids[!is.na(reaction_ids) & nzchar(reaction_ids)])
    if (length(reaction_ids) == 0) return(NULL)

    cache_file <- if (!is.null(cache_dir)) {
        file.path(cache_dir, "kegg_reaction_roles.rds")
    } else NULL
    cached <- if (!is.null(cache_file) && file.exists(cache_file)) {
        tryCatch(readRDS(cache_file), error = function(e) NULL)
    } else NULL
    if (!is.data.frame(cached)) {
        cached <- data.frame(reaction = character(0), compound = character(0),
                             role = character(0), stringsAsFactors = FALSE)
    }

    missing <- setdiff(reaction_ids, cached$reaction)
    if (length(missing) > 0) {
        fetched <- .fetch_reaction_equations(missing)
        if (!is.null(fetched) && nrow(fetched) > 0) {
            cached <- rbind(cached, fetched)
            if (!is.null(cache_file)) {
                dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
                tryCatch(saveRDS(cached, cache_file), error = function(e) NULL)
            }
        }
    }

    out <- cached[cached$reaction %in% reaction_ids, , drop = FALSE]
    if (nrow(out) == 0) return(NULL)
    out[order(out$reaction, out$role, out$compound), , drop = FALSE]
}


#' Read EQUATION records for up to ten KEGG reactions at a time
#'
#' @param reaction_ids Reaction ids to fetch.
#' @return Data frame with \code{reaction}, \code{compound}, \code{role}, or
#'   NULL when no request returned an equation.
#' @keywords internal
.fetch_reaction_equations <- function(reaction_ids) {
    batches <- split(reaction_ids, ceiling(seq_along(reaction_ids) / 10))
    message("    Reading ", length(reaction_ids), " KEGG reaction equation(s) in ",
            length(batches), " request(s)")

    rows <- list()
    for (batch in batches) {
        url <- paste0("https://rest.kegg.jp/get/",
                      paste0("rn:", batch, collapse = "+"))
        lines <- tryCatch(readLines(url, warn = FALSE), error = function(e) {
            message("    KEGG reaction fetch failed: ", e$message)
            NULL
        })
        if (is.null(lines)) next
        rows <- c(rows, list(.parse_reaction_records(lines)))
    }
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) return(NULL)
    out <- do.call(rbind, rows)
    if (is.null(out) || nrow(out) == 0) NULL else out
}


#' Parse a KEGG flat-file response into reaction roles
#'
#' Entries are separated by "///", the id sits on the ENTRY line and the
#' equation on the EQUATION line; KEGG wraps a long equation onto indented
#' continuation lines, which are joined back on.
#'
#' @param lines Lines of a \code{/get/} response.
#' @return Data frame with \code{reaction}, \code{compound} and \code{role};
#'   zero rows when the response holds no equation.
#' @keywords internal
.parse_reaction_records <- function(lines) {
    empty <- data.frame(reaction = character(0), compound = character(0),
                        role = character(0), stringsAsFactors = FALSE)
    if (length(lines) == 0) return(empty)

    breaks <- which(trimws(lines) == "///")
    starts <- c(1L, utils::head(breaks + 1L, -1L))
    ends <- if (length(breaks) > 0) breaks else length(lines)
    if (length(starts) == 0) return(empty)

    out <- list()
    for (i in seq_along(ends)) {
        block <- lines[starts[i]:ends[i]]
        entry <- grep("^ENTRY\\s+R[0-9]{5}", block, value = TRUE)
        if (length(entry) == 0) next
        rid <- regmatches(entry[1], regexpr("R[0-9]{5}", entry[1]))

        eq_start <- grep("^EQUATION", block)
        if (length(eq_start) == 0) next
        eq_lines <- block[eq_start[1]]
        j <- eq_start[1] + 1L
        while (j <= length(block) && grepl("^\\s{4,}\\S", block[j])) {
            eq_lines <- paste(eq_lines, trimws(block[j]))
            j <- j + 1L
        }
        roles <- parse_kegg_equation(sub("^EQUATION\\s+", "", eq_lines))
        if (nrow(roles) == 0) next
        out[[length(out) + 1]] <- data.frame(
            reaction = rid, compound = roles$compound, role = roles$role,
            stringsAsFactors = FALSE)
    }
    if (length(out) == 0) return(empty)
    do.call(rbind, out)
}
