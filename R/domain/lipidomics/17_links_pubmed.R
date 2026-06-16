# R/domain/lipidomics/17_links_pubmed.R
#
# Literature linker: for each (lipid_class, feature_name) pair from the
# pathway-link candidates, query NCBI PubMed via E-utilities, fetch
# abstracts, filter by journal-IF (>= min_journal_if) and year, then keep
# only papers where the two terms co-occur in the same sentence.
#
# Result schema (one row per kept paper):
#   layer, contrast, lipid_class, feature_name, pmid, year, journal,
#   journal_if, title, sentence


load_journal_if_table <- function(file = NULL, min_if = 15) {
    if (is.null(file)) {
        candidates <- c(
            file.path(getwd(), "data", "lipidomics", "journal_if_2024.tsv"),
            file.path("/home/ozsol/multiomics-core", "data", "lipidomics",
                      "journal_if_2024.tsv")
        )
        file <- candidates[file.exists(candidates)][1]
    }
    if (is.na(file) || !nzchar(file) || !file.exists(file)) {
        stop("journal_if_2024.tsv not found")
    }
    tbl <- utils::read.table(file, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE, check.names = FALSE,
                              quote = "")
    tbl <- tbl[!is.na(tbl$impact_factor) & tbl$impact_factor >= min_if, ,
               drop = FALSE]
    tbl$journal_full_lc <- tolower(tbl$journal_full)
    tbl$journal_iso_lc  <- tolower(tbl$journal_iso)
    tbl
}


# Run the literature linker against a candidate pair table.
# candidate_pairs: data.frame with columns layer, contrast, lipid_class,
#                  feature_name, [score] used only for ranking input
# Returns a list with $papers (long table of kept papers) and $summary
# (per-pair count of supporting papers).
run_pubmed_literature_links <- function(candidate_pairs, journal_if_tbl,
                                          cache_dir,
                                          year_min = 2021,
                                          api_key = NULL,
                                          sleep_seconds = 0.34,
                                          context_mode = "same_sentence",
                                          max_papers_per_pair = 50) {
    if (is.null(candidate_pairs) || nrow(candidate_pairs) == 0) return(NULL)
    ensure_dir(cache_dir)

    # De-duplicate pair queries (same class+feature can come from multiple
    # contrasts / pathway hits)
    keys <- unique(candidate_pairs[, c("lipid_class", "feature_name"),
                                    drop = FALSE])
    keys <- keys[!is.na(keys$lipid_class) & nzchar(keys$lipid_class) &
                  !is.na(keys$feature_name) & nzchar(keys$feature_name), ,
                  drop = FALSE]
    if (nrow(keys) == 0) return(NULL)
    message("[links/pubmed] ", nrow(keys),
            " unique (class, feature) pairs to query")

    paper_rows <- list()
    for (i in seq_len(nrow(keys))) {
        cls  <- keys$lipid_class[i]
        feat <- keys$feature_name[i]
        cache_file <- file.path(
            cache_dir,
            paste0(.safe_token(cls), "__", .safe_token(feat), ".json")
        )
        if (file.exists(cache_file)) {
            res <- tryCatch(jsonlite::fromJSON(cache_file,
                                                simplifyDataFrame = FALSE),
                            error = function(e) NULL)
        } else {
            res <- .pkg_one_pubmed_pair(cls, feat,
                                         year_min = year_min,
                                         api_key = api_key,
                                         max_papers = max_papers_per_pair)
            jsonlite::write_json(res %||% list(), cache_file,
                                  auto_unbox = TRUE, pretty = TRUE)
            Sys.sleep(sleep_seconds)
        }
        if (is.null(res) || length(res) == 0) next
        df <- .pkg_pubmed_to_df(res, cls, feat)
        paper_rows[[i]] <- df
    }

    if (length(paper_rows) == 0) return(list(papers = NULL, summary = NULL))
    all_papers <- do.call(rbind, paper_rows)
    if (is.null(all_papers) || nrow(all_papers) == 0) {
        return(list(papers = NULL, summary = NULL))
    }

    # IF + year filter
    jif_lc <- tolower(all_papers$journal)
    in_if <- jif_lc %in% c(journal_if_tbl$journal_full_lc,
                            journal_if_tbl$journal_iso_lc)
    all_papers <- all_papers[in_if, , drop = FALSE]
    if (nrow(all_papers) == 0) return(list(papers = NULL, summary = NULL))
    all_papers$journal_if <- .pkg_lookup_if(all_papers$journal,
                                              journal_if_tbl)

    # Context filter: same sentence (or same abstract as fallback)
    keep_rows <- list()
    for (i in seq_len(nrow(all_papers))) {
        cls  <- all_papers$lipid_class[i]
        feat <- all_papers$feature_name[i]
        text <- paste(all_papers$title[i], all_papers$abstract[i],
                      sep = ". ")
        ctx <- .pkg_context_match(text, cls, feat, mode = context_mode)
        if (length(ctx) > 0) {
            keep_rows[[i]] <- data.frame(
                all_papers[i, c("lipid_class", "feature_name", "pmid",
                                 "year", "journal", "journal_if", "title")],
                sentence = paste(head(ctx, 2), collapse = " | "),
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(keep_rows) == 0) return(list(papers = NULL, summary = NULL))
    kept <- do.call(rbind, keep_rows)

    # Re-attach to per-contrast candidate rows so downstream combine has
    # the same row granularity.
    kept_join <- merge(
        candidate_pairs[, c("layer", "contrast", "lipid_class",
                             "feature_name"), drop = FALSE],
        kept,
        by = c("lipid_class", "feature_name"), all = FALSE
    )

    summary_tbl <- aggregate(
        pmid ~ layer + contrast + lipid_class + feature_name,
        data = kept_join,
        FUN = function(x) length(unique(x))
    )
    names(summary_tbl)[names(summary_tbl) == "pmid"] <- "n_papers"

    list(papers = kept_join, summary = summary_tbl)
}


# -- helpers ---------------------------------------------------------------

.safe_token <- function(x) {
    gsub("[^A-Za-z0-9._-]+", "_", x)
}


# Single pair: esearch -> efetch -> parsed list. Returns list of paper dicts.
.pkg_one_pubmed_pair <- function(lipid_class, feature_name,
                                   year_min, api_key, max_papers) {
    if (!requireNamespace("xml2", quietly = TRUE)) {
        stop("xml2 not available")
    }
    base_search <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    base_fetch  <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    # Build query: class + feature in title/abstract, date range
    cls_q <- paste0("\"", lipid_class, "\"[Title/Abstract]")
    feat_q <- paste0("\"", feature_name, "\"[Title/Abstract]")
    date_q <- paste0(year_min, ":3000[dp]")
    term <- paste(cls_q, feat_q, date_q, sep = " AND ")

    params <- c(db = "pubmed", term = utils::URLencode(term),
                retmax = as.character(max_papers), retmode = "xml")
    if (!is.null(api_key) && nzchar(api_key)) {
        params <- c(params, api_key = api_key)
    }
    url <- paste0(base_search, "?",
                  paste(names(params), unname(params), sep = "=",
                        collapse = "&"))
    esearch_xml <- tryCatch(
        xml2::read_xml(url),
        error = function(e) {
            message("[pubmed] esearch failed: ", conditionMessage(e))
            NULL
        }
    )
    if (is.null(esearch_xml)) return(NULL)
    pmids <- xml2::xml_text(xml2::xml_find_all(esearch_xml, ".//Id"))
    if (length(pmids) == 0) return(list())

    # efetch
    fetch_params <- c(db = "pubmed", id = paste(pmids, collapse = ","),
                      rettype = "abstract", retmode = "xml")
    if (!is.null(api_key) && nzchar(api_key)) {
        fetch_params <- c(fetch_params, api_key = api_key)
    }
    fetch_url <- paste0(base_fetch, "?",
                        paste(names(fetch_params), unname(fetch_params),
                              sep = "=", collapse = "&"))
    efetch_xml <- tryCatch(
        xml2::read_xml(fetch_url),
        error = function(e) {
            message("[pubmed] efetch failed: ", conditionMessage(e))
            NULL
        }
    )
    if (is.null(efetch_xml)) return(NULL)

    articles <- xml2::xml_find_all(efetch_xml, ".//PubmedArticle")
    out <- lapply(articles, function(art) {
        pmid <- xml2::xml_text(xml2::xml_find_first(
            art, ".//MedlineCitation/PMID"
        ))
        title <- xml2::xml_text(xml2::xml_find_first(
            art, ".//Article/ArticleTitle"
        ))
        abs_nodes <- xml2::xml_find_all(art, ".//Abstract/AbstractText")
        abstract <- if (length(abs_nodes) > 0) {
            paste(xml2::xml_text(abs_nodes), collapse = " ")
        } else ""
        journal <- xml2::xml_text(xml2::xml_find_first(
            art, ".//Article/Journal/ISOAbbreviation"
        ))
        if (is.na(journal) || !nzchar(journal)) {
            journal <- xml2::xml_text(xml2::xml_find_first(
                art, ".//Article/Journal/Title"
            ))
        }
        year_node <- xml2::xml_find_first(art,
            ".//Article/Journal/JournalIssue/PubDate/Year")
        year <- xml2::xml_text(year_node)
        if (is.na(year) || !nzchar(year)) {
            md <- xml2::xml_text(xml2::xml_find_first(art,
                ".//Article/Journal/JournalIssue/PubDate/MedlineDate"))
            year <- regmatches(md, regexpr("\\d{4}", md))
            if (length(year) == 0) year <- NA_character_
        }
        list(pmid = pmid, title = title %||% "",
             abstract = abstract, journal = journal %||% "",
             year = suppressWarnings(as.integer(year)))
    })
    out
}


# Turn list-of-paper-dicts into a flat data.frame and tag with class/feature.
.pkg_pubmed_to_df <- function(res, lipid_class, feature_name) {
    if (length(res) == 0) return(NULL)
    do.call(rbind, lapply(res, function(p) {
        if (is.null(p) || is.null(p$pmid) || !nzchar(p$pmid)) return(NULL)
        data.frame(
            lipid_class  = lipid_class,
            feature_name = feature_name,
            pmid         = p$pmid,
            title        = p$title %||% "",
            abstract     = p$abstract %||% "",
            journal      = p$journal %||% "",
            year         = if (is.null(p$year)) NA_integer_ else p$year,
            stringsAsFactors = FALSE
        )
    }))
}


.pkg_lookup_if <- function(journals, journal_if_tbl) {
    lc <- tolower(journals)
    out <- rep(NA_real_, length(lc))
    # Match against ISO abbrev first, then full title
    iso_idx <- match(lc, journal_if_tbl$journal_iso_lc)
    out[!is.na(iso_idx)] <- journal_if_tbl$impact_factor[iso_idx[!is.na(iso_idx)]]
    miss <- is.na(out)
    full_idx <- match(lc[miss], journal_if_tbl$journal_full_lc)
    out[miss][!is.na(full_idx)] <-
        journal_if_tbl$impact_factor[full_idx[!is.na(full_idx)]]
    out
}


# Find sentences containing both terms (case-insensitive, word boundary).
# Returns vector of matched sentence strings.
.pkg_context_match <- function(text, term_a, term_b, mode = "same_sentence") {
    if (!nzchar(text)) return(character(0))
    # Crude sentence split: keep handles for "et al." and "e.g." minimally.
    s <- strsplit(text, "(?<=[.!?])\\s+(?=[A-Z0-9])", perl = TRUE)[[1]]
    if (length(s) == 0) s <- text
    pat_a <- paste0("(?i)\\b", regex_escape(term_a), "\\b")
    pat_b <- paste0("(?i)\\b", regex_escape(term_b), "\\b")
    has_a <- grepl(pat_a, s, perl = TRUE)
    has_b <- grepl(pat_b, s, perl = TRUE)

    if (mode == "same_abstract") {
        if (any(has_a) && any(has_b)) return(text) else return(character(0))
    }
    # default: same sentence (or fallback if both present in adjacent
    # sentences)
    same <- has_a & has_b
    if (any(same)) return(s[same])
    # adjacent fallback: window of 1
    n <- length(s)
    if (n >= 2) {
        for (i in seq_len(n - 1)) {
            if ((has_a[i] && has_b[i + 1]) || (has_b[i] && has_a[i + 1])) {
                return(paste(s[i], s[i + 1], sep = " "))
            }
        }
    }
    character(0)
}


regex_escape <- function(x) {
    gsub("([\\\\^$.|?*+()\\[\\]{}])", "\\\\\\1", x, perl = TRUE)
}
