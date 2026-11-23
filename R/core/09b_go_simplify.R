# ==============================================================================
# GO TERM SIMPLIFICATION FOR THE PATHWAY STEP (fgsea / ORA)
# ==============================================================================
#
# Opt-in (pathway.simplify_go, default false) for the proteomics and RNA pathway
# steps that run run_pathway_analysis(). Redundant GO terms are grouped by
# semantic similarity and one representative per group is kept in a SEPARATE
# simplified table; the full result tables and the pathway_results object that
# downstream targets consume are never modified.
#
# Restores the feature added in 4564b09 (simplify_go_results via rrvgo), whose
# definition was lost in merge 29ffe3e. Differences from that version: the full
# tables stay the primary output, similarity comes from the offline GOSemSim data
# built by .go_semdata() (R/core/09_enrichment.R), terms are grouped per GO
# ontology instead of assuming BP, and only terms below the significance cutoff
# are grouped.

#' GO semantic-similarity measures accepted by pathway.simplify_measure
#'
#' Wang uses only the GO graph (GO.db). The others are information-content
#' measures and also need the organism's OrgDb package.
#' @noRd
.go_simplify_measures <- c("Wang", "Resnik", "Lin", "Jiang", "Rel")

#' Resolve the GO simplification settings of a pathway config block
#'
#' The single reader of \code{pathway.simplify_go}, \code{simplify_threshold}
#' and \code{simplify_measure}; the pathway step and anything describing it go
#' through here rather than reading the keys themselves.
#'
#' @param pw_cfg The mode's \code{pathway} config block (may be NULL).
#' @param organism Organism name, used only to find the OrgDb that the
#'   information-content measures need. Ignored for Wang.
#' @return A list with \code{enabled} (logical), \code{threshold} (numeric),
#'   \code{measure} (character) and \code{orgdb} (OrgDb package name or NULL;
#'   always NULL for Wang).
#' @examples
#' resolve_go_simplify_settings(list(simplify_go = TRUE))
resolve_go_simplify_settings <- function(pw_cfg, organism = NULL) {
    pw_cfg <- pw_cfg %||% list()
    measure <- pw_cfg$simplify_measure %||% "Wang"
    orgdb <- NULL
    if (!identical(measure, "Wang") && !is.null(organism)) {
        org_info <- tryCatch(get_organism_info(organism), error = function(e) NULL)
        if (!is.null(org_info) && !is.null(org_info$orgdb) && !is.na(org_info$orgdb)) {
            orgdb <- org_info$orgdb
        }
    }
    list(
        enabled   = isTRUE(pw_cfg$simplify_go),
        threshold = pw_cfg$simplify_threshold %||% 0.7,
        measure   = measure,
        orgdb     = orgdb
    )
}

#' Validate the GO simplification keys of a pathway config block
#'
#' Called from the proteomics and RNA config validators.
#'
#' @param pw_cfg The mode's \code{pathway} config block (may be NULL).
#' @param prefix Key prefix used in error messages.
#' @return \code{invisible(TRUE)}; stops on an invalid value.
validate_go_simplify_config <- function(pw_cfg, prefix = "pathway") {
    if (is.null(pw_cfg)) return(invisible(TRUE))
    assert_scalar_bool(pw_cfg$simplify_go, paste0(prefix, "$simplify_go"), allow_null = TRUE)
    assert_scalar_num(pw_cfg$simplify_threshold, paste0(prefix, "$simplify_threshold"),
                      allow_null = TRUE, min_val = 0, max_val = 1)
    assert_one_of(pw_cfg$simplify_measure, paste0(prefix, "$simplify_measure"),
                  .go_simplify_measures, allow_null = TRUE)
    invisible(TRUE)
}

#' Look up the GO ontology (BP / MF / CC) of each GO identifier
#'
#' @param go_ids Character vector of GO identifiers.
#' @return Character vector named by \code{go_ids}; NA where GO.db does not
#'   know the term (obsolete or malformed) or GO.db is not installed.
#' @noRd
.go_term_ontology <- function(go_ids) {
    ont <- setNames(rep(NA_character_, length(go_ids)), go_ids)
    if (length(go_ids) == 0 || !requireNamespace("GO.db", quietly = TRUE)) return(ont)
    tryCatch({
        valid_keys <- AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID")
        known <- unique(go_ids[go_ids %in% valid_keys])
        if (length(known) > 0) {
            tab <- suppressMessages(AnnotationDbi::select(
                GO.db::GO.db, keys = known, columns = "ONTOLOGY", keytype = "GOID"))
            ont[] <- tab$ONTOLOGY[match(go_ids, tab$GOID)]
        }
    }, error = function(e) {
        message("    GO ontology lookup failed: ", e$message)
    })
    ont
}

#' Pairwise GO semantic similarity within one ontology
#'
#' @param go_ids GO identifiers, all from ontology \code{ont}.
#' @param ont "BP", "MF" or "CC".
#' @param measure One of \code{.go_simplify_measures}.
#' @param orgdb OrgDb package name for the information-content measures, or NULL.
#' @return Square numeric matrix named by the terms that have a similarity
#'   (terms with none are dropped, as rrvgo did), or NULL when the semantic data
#'   could not be built.
#' @noRd
.go_similarity_matrix <- function(go_ids, ont, measure, orgdb = NULL) {
    sem <- .go_semdata(ont, orgdb = if (identical(measure, "Wang")) NULL else orgdb)
    if (is.null(sem)) return(NULL)
    sim <- tryCatch(
        GOSemSim::termSim(go_ids, go_ids, semData = sem, method = measure),
        error = function(e) {
            message("    GO semantic similarity failed (", ont, ", ", measure, "): ", e$message)
            NULL
        }
    )
    if (is.null(sim)) return(NULL)
    sim <- as.matrix(sim)
    diag(sim) <- NA
    keep <- !apply(sim, 1, function(x) all(is.na(x)))
    sim <- sim[keep, keep, drop = FALSE]
    sim[is.na(sim)] <- 0
    diag(sim) <- 1
    sim
}

#' Group GO terms by similarity and pick one representative per group
#'
#' Complete-linkage hierarchical clustering of \code{1 - similarity}, cut at
#' height \code{threshold} (the convention of rrvgo::reduceSimMatrix(), which the
#' original implementation used): higher thresholds merge more terms. Within each
#' group the term with the lowest p-value is the representative; ties go to the
#' lower adjusted p-value, then to the earlier term.
#'
#' @param sim_matrix Square similarity matrix (values 0-1) named by term.
#' @param pval Named numeric vector of p-values for the same terms.
#' @param padj Named numeric vector of adjusted p-values (or NULL).
#' @param threshold Cut height, 0-1.
#' @return Data frame with columns \code{term}, \code{cluster} (integer) and
#'   \code{representative}, one row per term of \code{sim_matrix}, in its order.
#' @examples
#' sim <- matrix(c(1, 0.9, 0.1, 0.9, 1, 0.1, 0.1, 0.1, 1), 3,
#'               dimnames = list(c("GO:1", "GO:2", "GO:3"), c("GO:1", "GO:2", "GO:3")))
#' collapse_go_terms(sim, pval = c("GO:1" = 0.02, "GO:2" = 0.01, "GO:3" = 0.03))
collapse_go_terms <- function(sim_matrix, pval, padj = NULL, threshold = 0.7) {
    terms <- rownames(sim_matrix)
    if (length(terms) < 2) {
        return(data.frame(term = terms, cluster = seq_along(terms),
                          representative = terms, stringsAsFactors = FALSE))
    }
    hc <- stats::hclust(stats::as.dist(1 - sim_matrix), method = "complete")
    cluster <- stats::cutree(hc, h = threshold)[terms]

    p  <- as.numeric(pval[terms])
    pa <- if (is.null(padj)) rep(NA_real_, length(terms)) else as.numeric(padj[terms])
    rank_order <- order(p, pa, seq_along(terms), na.last = TRUE)
    rep_of <- tapply(terms[rank_order], cluster[rank_order], function(x) x[1])

    data.frame(term = terms,
               cluster = as.integer(cluster),
               representative = unname(rep_of[as.character(cluster)]),
               stringsAsFactors = FALSE)
}

#' Simplify the GO terms of one enrichment table
#'
#' Only rows whose \code{pathway} is a GO identifier take part; other rows pass
#' through unchanged. Of those, terms with adjusted p below \code{sig_cutoff} (at
#' most \code{max_terms}, lowest p first) are grouped per GO ontology with
#' \code{collapse_go_terms()}. The simplified table drops the non-representative
#' members of each group and keeps every other row of \code{results_df}.
#'
#' @param results_df One fgsea or ORA table from run_pathway_analysis().
#' @param settings Output of resolve_go_simplify_settings().
#' @param sig_cutoff Adjusted p-value below which terms are grouped (default 0.10,
#'   as in the original implementation).
#' @param max_terms Maximum number of terms grouped per table (default 500).
#' @return A list with \code{simplified} (data frame: the kept rows plus
#'   \code{go_cluster}, \code{n_terms_in_cluster}, \code{collapsed_terms}),
#'   \code{map} (data frame: one row per grouped term with its representative),
#'   \code{n_considered} (integer) and \code{status} (character).
simplify_go_results <- function(results_df, settings, sig_cutoff = 0.10, max_terms = 500) {
    empty <- function(status) {
        list(simplified = results_df, map = NULL, n_considered = 0L, status = status)
    }
    if (is.null(results_df) || nrow(results_df) == 0 || !"pathway" %in% colnames(results_df)) {
        return(empty("skipped: no GO terms"))
    }
    pval_col <- intersect(c("pval", "pvalue"), colnames(results_df))[1]
    if (is.na(pval_col) || !"padj" %in% colnames(results_df)) {
        return(empty("skipped: no p-value columns"))
    }

    ids  <- as.character(results_df$pathway)
    pval <- as.numeric(results_df[[pval_col]])
    padj <- as.numeric(results_df$padj)
    cand <- which(grepl("^GO:[0-9]{7}$", ids) & !is.na(padj) & padj < sig_cutoff)
    cand <- cand[order(pval[cand], padj[cand], cand)]
    if (length(cand) > max_terms) cand <- cand[seq_len(max_terms)]
    if (length(cand) < 2) {
        return(empty(sprintf("skipped: fewer than 2 GO terms with padj < %s", sig_cutoff)))
    }

    ont <- .go_term_ontology(ids[cand])
    maps <- list()
    sem_missing <- FALSE
    for (o in sort(unique(stats::na.omit(ont)))) {
        rows <- cand[!is.na(ont) & ont == o]
        if (length(rows) < 2) next
        sim <- .go_similarity_matrix(ids[rows], o, settings$measure, settings$orgdb)
        if (is.null(sim)) { sem_missing <- TRUE; next }
        if (nrow(sim) < 2) next
        rows <- rows[ids[rows] %in% rownames(sim)]
        collapsed <- collapse_go_terms(sim[ids[rows], ids[rows], drop = FALSE],
                                       pval = setNames(pval[rows], ids[rows]),
                                       padj = setNames(padj[rows], ids[rows]),
                                       threshold = settings$threshold)
        collapsed$row <- rows
        collapsed$ontology <- o
        collapsed$go_cluster <- paste0(o, "_", collapsed$cluster)
        maps[[o]] <- collapsed
    }
    if (length(maps) == 0) {
        return(empty(if (sem_missing) "skipped: GO semantic data unavailable"
                     else "skipped: fewer than 2 GO terms per ontology"))
    }

    m <- do.call(rbind, unname(maps))
    m <- m[order(m$row), , drop = FALSE]
    term_name <- if ("pathway_name" %in% colnames(results_df)) {
        as.character(results_df$pathway_name)
    } else {
        ids
    }
    name_of <- setNames(term_name, ids)
    map <- data.frame(
        pathway             = m$term,
        pathway_name        = unname(name_of[m$term]),
        ontology            = m$ontology,
        go_cluster          = m$go_cluster,
        representative      = m$representative,
        representative_name = unname(name_of[m$representative]),
        is_representative   = m$term == m$representative,
        pval                = pval[m$row],
        padj                = padj[m$row],
        stringsAsFactors    = FALSE
    )

    drop_rows <- m$row[m$term != m$representative]
    keep <- setdiff(seq_len(nrow(results_df)), drop_rows)
    simplified <- results_df[keep, , drop = FALSE]
    rep_rows <- match(ids[keep], map$pathway[map$is_representative])
    rep_cluster <- map$go_cluster[map$is_representative][rep_rows]
    simplified$go_cluster <- rep_cluster
    cluster_sizes <- table(map$go_cluster)
    simplified$n_terms_in_cluster <- as.integer(cluster_sizes[rep_cluster])
    simplified$collapsed_terms <- vapply(rep_cluster, function(cl) {
        if (is.na(cl)) return("")
        paste(map$pathway[map$go_cluster == cl & !map$is_representative], collapse = ";")
    }, character(1), USE.NAMES = FALSE)
    rownames(simplified) <- NULL

    list(simplified = simplified, map = map, n_considered = nrow(map), status = "simplified")
}

#' Simplify the GO tables of every contrast in a pathway result
#'
#' Tables without GO identifiers (KEGG, Reactome, custom GMTs with other IDs)
#' are not touched and not included. \code{pathway_results} is not modified.
#'
#' @param pathway_results Output of run_pathway_analysis().
#' @param settings Output of resolve_go_simplify_settings().
#' @param sig_cutoff,max_terms Passed to simplify_go_results().
#' @return A list with \code{tables} (named by "<contrast>|<analysis>", each a
#'   simplify_go_results() result plus \code{contrast} and \code{analysis}) and
#'   \code{summary} (one row per GO table).
simplify_pathway_results <- function(pathway_results, settings,
                                     sig_cutoff = 0.10, max_terms = 500) {
    tables <- list()
    summary_rows <- list()
    for (contrast_name in names(pathway_results)) {
        contrast_res <- pathway_results[[contrast_name]]
        for (analysis_name in names(contrast_res)) {
            res_df <- contrast_res[[analysis_name]]
            if (!is.data.frame(res_df) || nrow(res_df) == 0) next
            if (!any(grepl("^GO:[0-9]{7}$", as.character(res_df$pathway %||% character(0))))) next

            res <- simplify_go_results(res_df, settings, sig_cutoff = sig_cutoff,
                                       max_terms = max_terms)
            res$contrast <- contrast_name
            res$analysis <- analysis_name
            tables[[paste(contrast_name, analysis_name, sep = "|")]] <- res

            n_after <- nrow(res$simplified)
            summary_rows[[length(summary_rows) + 1]] <- data.frame(
                contrast = contrast_name, analysis = analysis_name,
                n_terms_full = nrow(res_df), n_go_terms_grouped = res$n_considered,
                n_terms_simplified = n_after, n_terms_collapsed = nrow(res_df) - n_after,
                sig_cutoff = sig_cutoff, max_terms = max_terms,
                threshold = settings$threshold, measure = settings$measure,
                status = res$status, stringsAsFactors = FALSE)
            message(sprintf("  GO simplification %s / %s: %s (%d -> %d rows)",
                            contrast_name, analysis_name, res$status, nrow(res_df), n_after))
        }
    }
    summary <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else NULL
    list(tables = tables, summary = summary)
}

#' Write simplified GO tables, term maps and the summary
#'
#' Files go to their own directory so that report globs over the full results
#' directory (e.g. \code{pathway_.*_ora.*\\.csv}) never pick them up:
#' \code{pathway_<contrast>_<analysis>_simplified.csv},
#' \code{pathway_<contrast>_<analysis>_simplified_map.csv} and
#' \code{go_simplification_summary.csv}.
#'
#' @param simplified Output of simplify_pathway_results().
#' @param output_dir Directory to write into (created if needed).
#' @return Character vector of written file paths (invisibly).
write_simplified_pathway_results <- function(simplified, output_dir) {
    if (is.null(simplified$summary)) return(invisible(character(0)))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    written <- character(0)
    for (res in simplified$tables) {
        if (!identical(res$status, "simplified")) next
        stem <- paste0("pathway_", gsub("[^a-zA-Z0-9_-]", "_", res$contrast), "_",
                       gsub("[^a-zA-Z0-9_-]", "_", res$analysis), "_simplified")
        f_tab <- file.path(output_dir, paste0(stem, ".csv"))
        f_map <- file.path(output_dir, paste0(stem, "_map.csv"))
        write.csv(res$simplified, f_tab, row.names = FALSE)
        write.csv(res$map, f_map, row.names = FALSE)
        written <- c(written, f_tab, f_map)
    }
    f_sum <- file.path(output_dir, "go_simplification_summary.csv")
    write.csv(simplified$summary, f_sum, row.names = FALSE)
    written <- c(written, f_sum)
    message("Saved ", length(written), " GO simplification file(s) to ", output_dir)
    invisible(written)
}

#' Methods sentence describing a GO simplification that ran
#'
#' Reads the summary written by write_simplified_pathway_results(), so the text
#' follows what ran rather than what the config asked for.
#'
#' @param summary_file Path to \code{go_simplification_summary.csv}.
#' @return A character scalar, or NULL when the file is absent or no table was
#'   simplified.
describe_go_simplification <- function(summary_file) {
    if (is.null(summary_file) || !file.exists(summary_file)) return(NULL)
    s <- tryCatch(read.csv(summary_file, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(s) || !"status" %in% names(s)) return(NULL)
    s <- s[s$status == "simplified", , drop = FALSE]
    if (nrow(s) == 0) return(NULL)
    sprintf(paste0(
        "For GO collections, a separate simplified table was also written: GO terms ",
        "with adjusted p < %s (at most %s per table) were grouped by semantic similarity ",
        "(GOSemSim, %s measure) using complete-linkage clustering of 1 - similarity cut ",
        "at %s, and the term with the lowest p-value in each group was kept. The full ",
        "enrichment tables are unchanged."),
        format(s$sig_cutoff[1]), format(s$max_terms[1]), s$measure[1], format(s$threshold[1]))
}

#' Run the optional GO simplification of a pathway step and write its files
#'
#' Entry point for the proteomics and RNA pathway steps. Does nothing unless
#' \code{pathway.simplify_go} is true. A failure is reported as a warning: the
#' full results are already written and are not affected.
#'
#' @param pathway_results Output of run_pathway_analysis().
#' @param pw_cfg The mode's \code{pathway} config block.
#' @param organism Organism name (for the OrgDb of IC-based measures).
#' @param output_dir Directory for the simplified files.
#' @return The simplify_pathway_results() result invisibly, or NULL when
#'   disabled or failed.
run_go_simplification <- function(pathway_results, pw_cfg, organism, output_dir) {
    settings <- resolve_go_simplify_settings(pw_cfg, organism)
    if (!settings$enabled) return(invisible(NULL))
    message("=== GO term simplification (threshold: ", settings$threshold,
            ", measure: ", settings$measure, ") ===")
    if (!identical(settings$measure, "Wang") && is.null(settings$orgdb)) {
        warning("GO simplification with the ", settings$measure, " measure needs an ",
                "organism OrgDb, and none is known for '", organism, "'. ",
                "Skipped; use simplify_measure: \"Wang\" or check annotation$organism.",
                call. = FALSE)
        return(invisible(NULL))
    }
    res <- tryCatch({
        out <- simplify_pathway_results(pathway_results, settings)
        write_simplified_pathway_results(out, output_dir)
        out
    }, error = function(e) {
        warning("GO simplification failed: ", conditionMessage(e),
                ". The full enrichment results are unaffected.", call. = FALSE)
        NULL
    })
    invisible(res)
}
