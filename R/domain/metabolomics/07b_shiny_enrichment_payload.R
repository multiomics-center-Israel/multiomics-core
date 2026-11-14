# R/domain/metabolomics/07b_shiny_enrichment_payload.R
#
# The metabolomics ORA and GSEA results in the shape the Shiny app reads for
# RNA-seq: the nine-component `enrichment` block, plus `enrichment_provenance`.
#
# This is a re-encoding, not an analysis. Every number comes from the result
# tables and every feature-level statement from the provenance recorded when
# the test ran (06r_enrichment_provenance.R). Nothing is recomputed, and a
# result whose provenance is not "ok" stops the build.
#
# Metabolomics differs from RNA-seq in one way that shapes everything here:
# the tests count COMPOUNDS, the app shows FEATURES, and several features can
# map to one compound. So four things are kept apart, never merged:
#   - the significant measured features that defined the ORA foreground
#     (`geneID` of an ORA row, and `features$significant` in the provenance);
#   - the compounds the ORA test counted (`Count`, `GeneRatio`, `BgRatio`);
#   - the pathways that are significant in an enrichment result (the display
#     tables under `enrichment$ora` / `enrichment$gsea`);
#   - exploratory membership: every measured feature of every pathway compound,
#     significant or not, ranked or not (`enrichment$pathway_membership`).
#
# The display tables hold the significant pathways only, as the RNA-seq ones do
# and as the app's validator requires. Every TESTED pathway, with its complete
# original statistics, is in `enrichment_provenance$<analysis>$<contrast>$all_tested`.
#
# build_shiny_payload_metabolomics() adds the two keys whenever the run has
# enrichment results (see metab_shiny_enrichment_keys()); there is no option.
#
# Reuses: %||%, .encode_genes(), validate_enrichment_payload() (RNA-seq payload),
#   build_metab_enrichment_manifest(), .METAB_ENRICHMENT_DB,
#   .METAB_ENRICHMENT_SIG_FDR (enrichment layout), validate_metab_enrichment_provenance(),
#   metab_pathway_feature_membership() (provenance).

.METAB_ENRICHMENT_PAYLOAD_SCHEMA <- 1L

# The one GSEA ranking metabolomics runs: log2 fold change. "fc" is the key the
# enrichment manifest and the app already use for it.
.METAB_ENRICHMENT_RANKING <- "fc"

# Where the non-directional ORA table of a contrast lives in the RNA-seq schema.
.METAB_ENRICHMENT_ORA_GROUP <- "contrasts_without_direction"


# ==== UNITS ====================================================================

#' The ORA and GSEA results of a run, one entry per analysis and contrast
#'
#' @param enrichment_res Result of \code{mod_metabolomics_enrichment()}.
#' @return List of units: list(analysis = "ORA"|"GSEA", method = "ora"|"gsea",
#'   contrast, res). Contrasts without a result table are left out.
#' @noRd
.metab_payload_units <- function(enrichment_res) {
    units <- list()
    for (method in c("ora", "gsea")) {
        by_contrast <- enrichment_res[[paste0(method, "_by_contrast")]] %||% list()
        for (contrast in names(by_contrast)) {
            res <- by_contrast[[contrast]]$res
            if (is.null(res) || !is.data.frame(res$table)) next
            units[[length(units) + 1L]] <- list(
                analysis = toupper(method), method = method,
                contrast = contrast, res = res)
        }
    }
    units
}

#' The pathways a unit passed to its test, with id and description
#'
#' @param provenance Provenance of one result.
#' @return The rows of \code{pathway_sources} that were given to the test (one
#'   per label).
#' @noRd
.metab_payload_tested_sources <- function(provenance) {
    src <- provenance$pathway_sources
    src[src$fate %in% c("tested", "not_returned"), , drop = FALSE]
}


# ==== PREFLIGHT ==================================================================

#' Check whether a run's ORA and GSEA results can become a Shiny payload
#'
#' Run this on real data before building the payload. It reports everything it
#' finds instead of stopping at the first problem. A \code{"blocker"} is a
#' reason the payload cannot be built correctly; \code{"info"} rows describe
#' the data and need no action.
#'
#' Blockers:
#' \itemize{
#'   \item a result without usable provenance;
#'   \item a feature id containing \code{"/"} (the app splits feature lists on it);
#'   \item duplicated or empty feature ids in the exported matrix;
#'   \item two tested pathways of one result sharing a pathway id, or a tested
#'     pathway without an id;
#'   \item a pathway id that names different pathways in different analyses or
#'     contrasts;
#'   \item a result row that provenance has no source for.
#' }
#'
#' @param enrichment_res Result of \code{mod_metabolomics_enrichment()}.
#' @param expr_features  Feature ids of the exported matrix, in its row order
#'   (\code{rownames(payload$expr_norm)}).
#' @param contrasts      Optional names of all DE contrasts, to report which
#'   have no ORA or GSEA result.
#' @return data.frame(check, level, detail), \code{level} being "ok", "info" or
#'   "blocker".
preflight_metab_enrichment_payload <- function(enrichment_res, expr_features,
                                               contrasts = NULL) {
    rows <- list()
    add <- function(check, level, detail = "") {
        rows[[length(rows) + 1L]] <<- data.frame(
            check = check, level = level, detail = detail, stringsAsFactors = FALSE)
    }
    some <- function(x, n = 5L) paste(utils::head(unique(x), n), collapse = ", ")

    expr_features <- as.character(expr_features)
    bad_expr <- is.na(expr_features) | !nzchar(expr_features)
    if (length(expr_features) == 0 || any(bad_expr)) {
        add("expr_feature_ids", "blocker",
            "the exported matrix has no feature ids, or some are empty")
    } else if (anyDuplicated(expr_features) > 0) {
        add("expr_feature_ids", "blocker", paste0(
            "duplicated feature ids in the exported matrix, e.g. ",
            some(expr_features[duplicated(expr_features)])))
    } else {
        add("expr_feature_ids", "ok")
    }

    units <- .metab_payload_units(enrichment_res)
    if (length(units) == 0) add("results", "info", "no ORA or GSEA result in this run")

    feature_ids <- expr_features
    id_labels   <- list()   # pathway id -> labels, over every unit
    for (u in units) {
        tag <- paste0(u$analysis, ":", u$contrast)
        problems <- validate_metab_enrichment_provenance(u$res$provenance)
        if (length(problems) > 0) {
            add(paste0("provenance:", tag), "blocker", paste(problems, collapse = "; "))
            next
        }
        add(paste0("provenance:", tag), "ok")
        p <- u$res$provenance
        feature_ids <- c(feature_ids, p$features$feature_id)

        src <- .metab_payload_tested_sources(p)
        ids <- src$pathway_id
        id_labels[[length(id_labels) + 1L]] <- data.frame(
            pathway_id = src$pathway_id, pathway = src$pathway, unit = rep(tag, nrow(src)),
            stringsAsFactors = FALSE)
        if (any(is.na(ids) | !nzchar(ids))) {
            add(paste0("pathway_ids:", tag), "blocker", "a tested pathway has no id")
        } else if (anyDuplicated(ids) > 0) {
            add(paste0("pathway_ids:", tag), "blocker", paste0(
                "two tested pathways share an id, e.g. ", some(ids[duplicated(ids)]),
                ". The app looks pathways up by id, so one would hide the other."))
        } else {
            add(paste0("pathway_ids:", tag), "ok")
        }
        unknown <- setdiff(as.character(u$res$table$pathway), src$pathway)
        if (length(unknown) > 0) {
            add(paste0("result_sources:", tag), "blocker", paste0(
                length(unknown), " result row(s) have no recorded source"))
        }

        fc <- p$feature_compounds
        if (identical(u$method, "ora")) {
            lost <- p$features$significant & is.na(p$features$compound_id)
            add(paste0("significant_features_without_compound:", tag), "info",
                paste0(sum(lost), " significant feature(s) have no usable compound ",
                       "id and are in no foreground"))
        }
        add(paste0("ambiguous_kegg_cells:", tag), "info",
            paste0(sum(fc$ambiguous), " feature(s) list more than one KEGG id; ",
                   "the first is used"))
    }

    # The block has ONE namespace of pathway ids for all analyses and contrasts
    # (pathway_membership is keyed by id alone). An id that names one pathway
    # in one unit and another pathway in another unit would merge the two.
    if (length(id_labels) > 0) {
        pairs <- unique(do.call(rbind, id_labels)[, c("pathway_id", "pathway")])
        pairs <- pairs[!is.na(pairs$pathway_id) & nzchar(pairs$pathway_id), , drop = FALSE]
        clash <- unique(pairs$pathway_id[duplicated(pairs$pathway_id)])
        if (length(clash) > 0) {
            add("pathway_ids_across_units", "blocker", paste0(
                "a pathway id names different pathways in different analyses or ",
                "contrasts, e.g. ", some(clash),
                ". The block keys pathways by id alone, so they would be merged."))
        } else {
            add("pathway_ids_across_units", "ok")
        }
    }

    slash <- unique(feature_ids[!is.na(feature_ids) & grepl("/", feature_ids, fixed = TRUE)])
    if (length(slash) > 0) {
        add("feature_ids_without_slash", "blocker", paste0(
            length(slash), " feature id(s) contain '/', e.g. ", some(slash),
            ". The app splits feature lists on '/', so these would be cut apart."))
    } else {
        add("feature_ids_without_slash", "ok")
    }

    outside <- setdiff(unique(feature_ids[!is.na(feature_ids)]), expr_features)
    add("features_outside_expr_norm", "info", paste0(
        length(outside), " feature(s) known to the enrichment are not rows of the ",
        "exported matrix; they are appended to gene_index and flagged"))

    if (!is.null(contrasts)) {
        have <- vapply(units, function(u) paste0(u$analysis, ":", u$contrast), character(1))
        want <- as.vector(outer(c("ORA", "GSEA"), as.character(contrasts), paste, sep = ":"))
        missing <- setdiff(want, have)
        add("units_without_result", "info",
            if (length(missing) == 0) "every contrast has an ORA and a GSEA result"
            else paste0(length(missing), " without a result: ",
                        paste(missing, collapse = ", ")))
    }

    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}


# ==== ONE UNIT ===================================================================

#' Re-encode one ORA result
#'
#' \code{Count}, \code{GeneRatio}, \code{BgRatio} and \code{FoldEnrichment} are
#' the test's own compound-based numbers. \code{geneID} lists the significant
#' measured features of the overlap compounds; \code{n_features} counts them,
#' and can exceed \code{Count} when a compound has several significant features.
#'
#' @param res    Result of \code{run_metabolomics_ora()} with "ok" provenance.
#' @param cutoff Adjusted-p threshold of the significant-only display table.
#' @return list(display, all_tested, membership, pathway_members, no_result).
#' @noRd
.metab_payload_ora_unit <- function(res, cutoff) {
    p   <- res$provenance
    tbl <- res$table
    src <- .metab_payload_tested_sources(p)
    k   <- match(as.character(tbl$pathway), src$pathway)

    membership <- metab_pathway_feature_membership(p)
    hit <- membership[membership$compound_in_overlap & membership$feature_significant, ,
                      drop = FALSE]
    features <- lapply(as.character(tbl$pathway), function(pw) {
        unique(hit$feature_id[hit$pathway == pw])
    })
    feature_ids <- vapply(features, paste, character(1), collapse = "/")
    n_features  <- lengths(features)
    significant <- !is.na(tbl$FDR) & tbl$FDR < cutoff

    description <- src$description[k]
    description[is.na(description) | !nzchar(description)] <- src$pathway_id[k][
        is.na(description) | !nzchar(description)]

    all_tested <- data.frame(pathway_id = src$pathway_id[k], description = description,
                             stringsAsFactors = FALSE)
    all_tested <- cbind(all_tested, tbl)
    all_tested$overlap_genes <- as.character(all_tested$overlap_genes)
    all_tested$significant   <- significant
    all_tested$feature_ids   <- feature_ids
    all_tested$n_features    <- as.integer(n_features)
    rownames(all_tested) <- NULL

    display <- data.frame(
        ID             = src$pathway_id[k],
        Description    = description,
        GeneRatio      = paste0(tbl$overlap, "/", tbl$fg_size),
        BgRatio        = paste0(tbl$pathway_size, "/", tbl$bg_size),
        pvalue         = tbl$raw_p,
        p.adjust       = tbl$FDR,
        geneID         = feature_ids,
        Count          = as.integer(tbl$overlap),
        FoldEnrichment = tbl$fold_enrichment,
        n_features     = as.integer(n_features),
        library        = as.character(tbl$library),
        stringsAsFactors = FALSE)[significant, , drop = FALSE]
    rownames(display) <- NULL

    members <- p$pathway_members
    members$pathway_id <- src$pathway_id[match(members$pathway, src$pathway)]

    list(display = display, all_tested = all_tested, membership = membership,
         pathway_members = members,
         no_result = .metab_payload_no_result(src, tbl))
}

#' Pathways that were given to a test and came back without a result
#'
#' fgsea keeps each member of a set once and then applies its size limits, so a
#' set that reached it with one compound listed twice is dropped and has no
#' row in the result. It was not tested, and it was not left out by the
#' pipeline either: it is listed here so that it does not simply vanish.
#' Fisher's test returns a row for every set it is given, so this is empty for
#' ORA.
#'
#' @param src Sources given to the test (\code{.metab_payload_tested_sources()}).
#' @param tbl The result table.
#' @return data.frame(pathway_id, description, pathway, library, n_listed, reason).
#' @noRd
.metab_payload_no_result <- function(src, tbl) {
    missing <- src[!(src$pathway %in% as.character(tbl$pathway)), , drop = FALSE]
    data.frame(pathway_id = missing$pathway_id, description = missing$description,
               pathway = missing$pathway, library = missing$library,
               n_listed = missing$n_listed,
               reason = rep("given_to_the_test_no_result_returned", nrow(missing)),
               stringsAsFactors = FALSE, row.names = NULL)
}

#' Re-encode one GSEA result
#'
#' The ranking and the leading edges are expressed in the representative
#' FEATURE of each compound, the one the test ranked. No other feature of a
#' compound appears in either.
#'
#' @param res    Result of \code{run_metabolomics_gsea()} with "ok" provenance.
#' @param cutoff Adjusted-p threshold of the significant-only display table.
#' @return list(display, all_tested, membership, pathway_members, ranking,
#'   leading_edge, no_result): \code{ranking} is list(features, score);
#'   \code{leading_edge} is a list of feature vectors named by pathway id, for
#'   the display rows.
#' @noRd
.metab_payload_gsea_unit <- function(res, cutoff) {
    p   <- res$provenance
    tbl <- res$table
    src <- .metab_payload_tested_sources(p)
    k   <- match(as.character(tbl$pathway), src$pathway)
    ranked <- p$ranked

    edges <- lapply(as.character(tbl$leading_edge), function(s) {
        if (is.na(s) || !nzchar(s)) return(character(0))
        compounds <- strsplit(s, ";", fixed = TRUE)[[1]]
        feats <- ranked$representative_feature[match(compounds, ranked$compound_id)]
        if (anyNA(feats)) {
            stop("[metabolomics enrichment payload] a leading-edge compound has no ",
                 "representative feature on record", call. = FALSE)
        }
        feats
    })
    significant <- !is.na(tbl$FDR) & tbl$FDR < cutoff

    description <- src$description[k]
    description[is.na(description) | !nzchar(description)] <- src$pathway_id[k][
        is.na(description) | !nzchar(description)]

    all_tested <- data.frame(pathway_id = src$pathway_id[k], description = description,
                             stringsAsFactors = FALSE)
    all_tested <- cbind(all_tested, tbl)
    all_tested$significant           <- significant
    all_tested$leading_edge_features <- vapply(edges, paste, character(1), collapse = "/")
    rownames(all_tested) <- NULL

    display <- data.frame(
        ID              = src$pathway_id[k],
        Description     = description,
        setSize         = as.integer(tbl$pathway_size),
        enrichmentScore = tbl$ES,
        NES             = tbl$NES,
        pvalue          = tbl$pval,
        p.adjust        = tbl$FDR,
        library         = as.character(tbl$library),
        stringsAsFactors = FALSE)[significant, , drop = FALSE]
    rownames(display) <- NULL

    members <- p$pathway_members
    members$pathway_id <- src$pathway_id[match(members$pathway, src$pathway)]

    list(display = display, all_tested = all_tested,
         membership = metab_pathway_feature_membership(p),
         pathway_members = members,
         ranking = list(features = ranked$representative_feature, score = ranked$logFC),
         leading_edge = stats::setNames(edges[significant], src$pathway_id[k][significant]),
         no_result = .metab_payload_no_result(src, tbl))
}


# ==== THE PAYLOAD ==================================================================

#' Build the Shiny enrichment block and its provenance from a metabolomics run
#'
#' Returns two things meant to sit side by side at the top level of the Shiny
#' payload:
#' \describe{
#'   \item{\code{enrichment}}{The nine components the app reads for RNA-seq:
#'     \code{available}, \code{gene_index}, \code{config}, \code{manifest},
#'     \code{ora}, \code{gsea}, \code{gsea_leading_edge}, \code{gsea_rankings},
#'     \code{pathway_membership}. Database \code{"combined"}, ranking
#'     \code{"fc"}; ids are FEATURE ids. Tables hold the significant pathways.}
#'   \item{\code{enrichment_provenance}}{Everything the block cannot hold: every
#'     tested pathway with its full statistics, compounds, per-feature roles,
#'     pathway members with multiplicity, and the state of each analysis.}
#' }
#' The two are joined by analysis ("ORA"/"GSEA"), contrast and pathway id.
#'
#' \code{gene_index} starts with \code{expr_features}, in that order, so its
#' first entries line up with the rows of the exported matrix. Features the
#' enrichment refers to that are not rows of that matrix follow, sorted, and
#' \code{enrichment_provenance$feature_index} says which are which.
#'
#' The significance threshold is the one the enrichment stage itself reports
#' with (\code{.METAB_ENRICHMENT_SIG_FDR}); there is no separate setting.
#'
#' An analysis that ran and found nothing significant has a manifest row with
#' status \code{"empty"} and a zero-row table. An analysis with no result has
#' no manifest row, and \code{enrichment_provenance$status} lists it as
#' \code{"no_result"}; the enrichment functions do not record why.
#'
#' @param enrichment_res Result of \code{mod_metabolomics_enrichment()}, or NULL.
#' @param expr_features  Feature ids of the exported matrix, in its row order.
#' @param contrasts      Optional names of all DE contrasts, so that contrasts
#'   without a result are listed in the status table.
#' @param config         Unused for now; accepted so that callers pass the same
#'   arguments as for RNA-seq.
#' @return list(enrichment, enrichment_provenance); both NULL when
#'   \code{enrichment_res} is NULL. Stops when the preflight finds a blocker.
build_metab_enrichment_payload <- function(enrichment_res, expr_features,
                                           contrasts = NULL, config = NULL) {
    if (is.null(enrichment_res)) {
        return(list(enrichment = NULL, enrichment_provenance = NULL))
    }
    expr_features <- as.character(expr_features)
    preflight <- preflight_metab_enrichment_payload(enrichment_res, expr_features,
                                                    contrasts)
    blockers <- preflight[preflight$level == "blocker", , drop = FALSE]
    if (nrow(blockers) > 0) {
        stop("[metabolomics enrichment payload] cannot be built:\n",
             paste0("  - ", blockers$check, ": ", blockers$detail, collapse = "\n"),
             call. = FALSE)
    }

    cutoff <- .METAB_ENRICHMENT_SIG_FDR
    db     <- .METAB_ENRICHMENT_DB
    rk     <- .METAB_ENRICHMENT_RANKING
    units  <- .metab_payload_units(enrichment_res)
    built  <- lapply(units, function(u) {
        if (identical(u$method, "ora")) .metab_payload_ora_unit(u$res, cutoff)
        else .metab_payload_gsea_unit(u$res, cutoff)
    })

    # ---- gene_index: the exported matrix first, then anything else referenced ----
    # Every feature the block or its provenance names, so that all of them can
    # be resolved and flagged: pathway members, ranked and leading-edge
    # features, and also features with no compound or outside every pathway.
    referenced <- unique(c(
        unlist(lapply(built, function(b) {
            c(b$membership$feature_id, b$ranking$features,
              unlist(b$leading_edge, use.names = FALSE))
        }), use.names = FALSE),
        unlist(lapply(units, function(u) {
            c(u$res$provenance$feature_compounds$feature_id,
              u$res$provenance$features$feature_id)
        }), use.names = FALSE)))
    referenced <- as.character(referenced)
    referenced <- referenced[!is.na(referenced)]
    gene_index <- c(expr_features, sort(setdiff(referenced, expr_features)))
    encode <- function(x, context) .encode_genes(x, gene_index, context = context)

    # ---- the nine components ----
    # Per-contrast lists are filled with `x[name] <- list(value)`, which keeps
    # an empty table or an empty leading-edge list as an entry of its own.
    ora_tabs <- list(); gsea_tabs <- list(); edge_tabs <- list(); rank_tabs <- list()
    membership <- list()
    for (i in seq_along(units)) {
        u <- units[[i]]; b <- built[[i]]
        if (identical(u$method, "ora")) {
            ora_tabs[u$contrast] <- list(b$display)
        } else {
            gsea_tabs[u$contrast] <- list(b$display)
            edge_tabs[u$contrast] <- list(lapply(
                b$leading_edge, encode, context = "GSEA leading-edge feature"))
            rank_tabs[u$contrast] <- list(list(
                idx   = encode(b$ranking$features, "ranking feature"),
                score = unname(as.numeric(b$ranking$score))))
        }
        # Exploratory membership: every measured feature of the pathway's
        # compounds. It does not depend on the contrast; results are merged.
        by_id <- split(b$membership$feature_id, b$membership$pathway_id)
        for (id in names(by_id)) {
            membership[[id]] <- union(membership[[id]], by_id[[id]])
        }
    }
    nest <- function(x, ...) {
        if (length(x) == 0) return(list())
        for (key in rev(c(...))) x <- stats::setNames(list(x), key)
        x
    }
    ora      <- nest(ora_tabs, db, .METAB_ENRICHMENT_ORA_GROUP)
    gsea     <- nest(gsea_tabs, db, rk)
    edges    <- nest(edge_tabs, db, rk)
    rankings <- nest(rank_tabs, rk)
    pathway_membership <- nest(
        lapply(membership, encode, context = "pathway membership feature"), db)

    # The manifest rows are the enrichment layout's own, so the app and the
    # files on disk describe the same units. Its storage_key is a directory
    # there; here it names the unit's place in this block instead.
    manifest <- build_metab_enrichment_manifest(lapply(units, function(u) {
        list(method = u$method, contrast = u$contrast, table = u$res$table,
             rel_dir = if (identical(u$method, "ora")) {
                 paste("ora", db, .METAB_ENRICHMENT_ORA_GROUP, u$contrast, sep = "/")
             } else {
                 paste("gsea", db, rk, u$contrast, sep = "/")
             })
    }))
    gsea_groups <- unique(manifest$group[manifest$analysis == "GSEA"])
    if (length(gsea_groups) > 0 && !identical(gsea_groups, rk)) {
        stop("[metabolomics enrichment payload] the enrichment manifest names the ",
             "GSEA ranking '", paste(gsea_groups, collapse = ", "), "', and this ",
             "block stores it as '", rk, "'", call. = FALSE)
    }

    libraries <- unique(unlist(lapply(units, function(u) {
        as.character(u$res$table$library)
    }), use.names = FALSE))
    enrichment <- list(
        available          = length(units) > 0,
        gene_index         = gene_index,
        config             = list(
            databases_run      = db,
            rankings           = if (length(rankings) > 0) rk else character(0),
            pvalue_cutoff      = cutoff,
            gsea_pvalue_cutoff = cutoff,
            padj_method        = "fdr",
            id_type            = "feature",
            libraries          = as.character(libraries[!is.na(libraries)])),
        manifest           = manifest,
        ora                = ora,
        gsea               = gsea,
        gsea_leading_edge  = edges,
        gsea_rankings      = rankings,
        pathway_membership = pathway_membership
    )
    validate_enrichment_payload(enrichment, strict = TRUE)

    # ---- provenance ----
    unit_provenance <- function(i) {
        u <- units[[i]]; b <- built[[i]]; p <- u$res$provenance
        sources <- p$pathway_sources
        sources$gmt_file <- basename(sources$gmt_file)
        c(list(all_tested = b$all_tested, no_result = b$no_result,
               membership = b$membership,
               pathway_members = b$pathway_members, pathway_sources = sources,
               features = p$features, diagnostics = p$diagnostics),
          if (identical(u$method, "ora")) {
              list(thresholds = p$thresholds, compounds = p$compounds)
          } else {
              list(ranked = p$ranked)
          })
    }
    by_analysis <- function(method) {
        idx <- which(vapply(units, function(u) identical(u$method, method), logical(1)))
        stats::setNames(lapply(idx, unit_provenance),
                        vapply(units[idx], function(u) u$contrast, character(1)))
    }

    all_contrasts <- unique(c(as.character(contrasts),
                              vapply(units, function(u) u$contrast, character(1))))
    status <- do.call(rbind, lapply(c("ORA", "GSEA"), function(analysis) {
        do.call(rbind, lapply(all_contrasts, function(contrast) {
            i <- which(vapply(units, function(u) {
                identical(u$analysis, analysis) && identical(u$contrast, contrast)
            }, logical(1)))
            if (length(i) == 0) {
                return(data.frame(analysis = analysis, contrast = contrast,
                                  state = "no_result", n_tested = NA_integer_,
                                  n_significant = NA_integer_,
                                  n_pathways_without_result = NA_integer_,
                                  provenance_status = NA_character_,
                                  stringsAsFactors = FALSE))
            }
            tested <- built[[i]]$all_tested
            n_sig  <- sum(tested$significant)
            data.frame(analysis = analysis, contrast = contrast,
                       state = if (n_sig > 0) "significant" else "no_significant_pathways",
                       n_tested = nrow(tested), n_significant = as.integer(n_sig),
                       n_pathways_without_result = nrow(built[[i]]$no_result),
                       provenance_status = units[[i]]$res$provenance$status,
                       stringsAsFactors = FALSE)
        }))
    }))
    if (is.null(status)) {
        status <- data.frame(analysis = character(0), contrast = character(0),
                             state = character(0), n_tested = integer(0),
                             n_significant = integer(0),
                             n_pathways_without_result = integer(0),
                             provenance_status = character(0), stringsAsFactors = FALSE)
    }

    reference_files <- unique(do.call(rbind, lapply(units, function(u) {
        ref <- u$res$provenance$reference_files
        data.frame(role = ref$role, file = basename(ref$path), md5 = ref$md5,
                   read = ref$read, stringsAsFactors = FALSE)
    })))
    if (!is.null(reference_files)) rownames(reference_files) <- NULL

    enrichment_provenance <- list(
        schema_version    = .METAB_ENRICHMENT_PAYLOAD_SCHEMA,
        id_type           = "feature",
        database          = db,
        ranking           = rk,
        significance      = list(column = "FDR", cutoff = cutoff),
        status            = status,
        feature_index     = data.frame(feature_id = gene_index,
                                       in_expr_norm = gene_index %in% expr_features,
                                       stringsAsFactors = FALSE),
        feature_compounds = if (length(units) > 0) {
            units[[1]]$res$provenance$feature_compounds
        },
        reference_files   = reference_files,
        ora               = by_analysis("ora"),
        gsea              = by_analysis("gsea"),
        preflight         = preflight
    )

    list(enrichment = enrichment, enrichment_provenance = enrichment_provenance)
}


#' The enrichment keys to add to a metabolomics Shiny payload
#'
#' What \code{build_shiny_payload_metabolomics()} adds, decided in one place:
#' \itemize{
#'   \item no enrichment in the run (\code{enrichment_res} is NULL): nothing;
#'   \item results that can be encoded: \code{enrichment} and
#'     \code{enrichment_provenance}. A run whose ORA and GSEA produced no
#'     result gets a block with \code{available = FALSE}, and each contrast
#'     listed as \code{no_result}; no statistic is made up for it;
#'   \item results that cannot be encoded correctly (the preflight reports a
#'     blocker, for instance results computed before provenance was recorded):
#'     no \code{enrichment} block, a warning, and an
#'     \code{enrichment_provenance} that says so:
#'     \code{list(schema_version, available = FALSE, reason = "blocked",
#'     preflight)}. The rest of the payload is unaffected.
#' }
#' The payload is always written: an enrichment problem must not cost the app
#' every other key. It is never dropped without a trace either.
#'
#' @param enrichment_res Result of \code{mod_metabolomics_enrichment()}, or NULL.
#' @param expr_features  Feature ids of the exported matrix, in its row order.
#' @param contrasts      Names of all DE contrasts (or NULL).
#' @param config         Passed on to \code{build_metab_enrichment_payload()}.
#' @return A named list of payload keys to set; empty when there is nothing to add.
metab_shiny_enrichment_keys <- function(enrichment_res, expr_features,
                                        contrasts = NULL, config = NULL) {
    if (is.null(enrichment_res)) return(list())

    preflight <- preflight_metab_enrichment_payload(enrichment_res, expr_features,
                                                    contrasts)
    blockers <- preflight[preflight$level == "blocker", , drop = FALSE]
    if (nrow(blockers) > 0) {
        # Raised so that it can never stop the payload: where warnings are
        # turned into errors it is dropped, and the reasons stay on record.
        .metab_prov_warn(
            "metabolomics Shiny payload: the ORA/GSEA enrichment block was not ",
            "added, because the results cannot be encoded correctly:\n",
            paste0("  - ", blockers$check, ": ", blockers$detail, collapse = "\n"),
            "\nEvery other key of the payload is written as usual.")
        return(list(enrichment_provenance = list(
            schema_version = .METAB_ENRICHMENT_PAYLOAD_SCHEMA,
            available      = FALSE,
            reason         = "blocked",
            preflight      = preflight)))
    }

    built <- build_metab_enrichment_payload(enrichment_res, expr_features,
                                            contrasts = contrasts, config = config)
    list(enrichment            = built$enrichment,
         enrichment_provenance = built$enrichment_provenance)
}


# ==== READING THE PAYLOAD ============================================================

#' Everything behind one pathway of one contrast
#'
#' Works for every tested pathway, significant or not. Reads only the payload,
#' so it can be used on a saved RDS without the pipeline's results.
#'
#' \code{state} is \code{"tested"}, or \code{"no_result"} for a pathway that was
#' given to the test and came back without a result (fgsea drops a set whose
#' members collapse below its size limit); \code{result} then has no
#' statistics. Otherwise \code{result} is the pathway's row with its complete
#' original statistics.
#' \code{compounds} has one row per member the library lists, duplicates and
#' unmeasured members included, with what the test did with it. \code{features}
#' has one row per measured feature of those compounds, with its role:
#' \itemize{
#'   \item ORA: \code{feature_significant}, \code{compound_in_universe},
#'     \code{compound_in_foreground}, \code{compound_in_overlap}.
#'   \item GSEA: \code{role}, \code{logFC}, \code{compound_ranked},
#'     \code{feature_is_representative}, \code{contributes_to_score},
#'     \code{in_leading_edge}.
#' }
#' A feature listed here did not necessarily enter the test: the flags say so.
#'
#' @param payload    A Shiny payload carrying \code{enrichment_provenance}.
#' @param analysis   "ORA" or "GSEA" (case-insensitive).
#' @param contrast   Contrast name.
#' @param pathway_id Pathway id (the \code{ID} column of the display tables,
#'   \code{pathway_id} of \code{all_tested}).
#' @return list(analysis, contrast, pathway_id, state, result, compounds,
#'   features).
metab_enrichment_pathway_detail <- function(payload, analysis, contrast, pathway_id) {
    prov <- payload$enrichment_provenance
    if (is.null(prov)) {
        stop("metab_enrichment_pathway_detail(): this payload has no ",
             "enrichment_provenance.", call. = FALSE)
    }
    if (isFALSE(prov$available)) {
        stop("metab_enrichment_pathway_detail(): this payload has no enrichment ",
             "block (reason: ", prov$reason %||% "not recorded", "). See ",
             "enrichment_provenance$preflight for what blocked it.", call. = FALSE)
    }
    key <- tolower(as.character(analysis))
    if (length(key) != 1L || !key %in% c("ora", "gsea")) {
        stop("metab_enrichment_pathway_detail(): analysis must be \"ORA\" or \"GSEA\".",
             call. = FALSE)
    }
    unit <- prov[[key]][[as.character(contrast)]]
    if (is.null(unit)) {
        st <- prov$status
        state <- st$state[st$analysis == toupper(key) & st$contrast == contrast]
        stop("metab_enrichment_pathway_detail(): there is no ", toupper(key),
             " result for contrast '", contrast, "'",
             if (length(state) == 1) paste0(" (state: ", state, ")"),
             ". Contrasts with one: ",
             paste(names(prov[[key]]), collapse = ", "), ".", call. = FALSE)
    }
    result <- unit$all_tested[unit$all_tested$pathway_id == pathway_id, , drop = FALSE]
    state  <- "tested"
    if (nrow(result) == 0) {
        # Given to the test, and no result came back (see `no_result`): there
        # are no statistics, but its compounds and features are on record.
        result <- unit$no_result[unit$no_result$pathway_id == pathway_id, , drop = FALSE]
        state  <- "no_result"
    }
    if (nrow(result) == 0) {
        stop("metab_enrichment_pathway_detail(): pathway '", pathway_id, "' was not ",
             "tested by ", toupper(key), " for '", contrast, "'. The tested pathways ",
             "are in enrichment_provenance$", key, "[[contrast]]$all_tested.",
             call. = FALSE)
    }

    members <- unit$pathway_members
    compounds <- members[!is.na(members$pathway_id) & members$pathway_id == pathway_id, ,
                         drop = FALSE]
    features <- unit$membership[unit$membership$pathway_id == pathway_id, , drop = FALSE]

    fc <- prov$feature_compounds
    at <- match(features$feature_id, fc$feature_id)
    features$id_source       <- fc$id_source[at]
    features$kegg_candidates <- fc$kegg_candidates[at]
    features$ambiguous       <- fc$ambiguous[at]
    features$in_expr_norm    <- prov$feature_index$in_expr_norm[
        match(features$feature_id, prov$feature_index$feature_id)]

    uf <- unit$features
    fa <- match(features$feature_id, uf$feature_id)
    if (identical(key, "gsea")) {
        features$role  <- uf$role[fa]
        features$logFC <- uf$logFC[fa]
        rk <- match(compounds$compound_id, unit$ranked$compound_id)
        compounds$rank                   <- unit$ranked$rank[rk]
        compounds$logFC                  <- unit$ranked$logFC[rk]
        compounds$representative_feature <- unit$ranked$representative_feature[rk]
    } else {
        features$in_de_table <- uf$in_de_table[fa]
        cp <- match(compounds$compound_id, unit$compounds$compound_id)
        compounds$n_features             <- unit$compounds$n_features[cp]
        compounds$n_significant_features <- unit$compounds$n_significant_features[cp]
    }
    rownames(result) <- NULL; rownames(compounds) <- NULL; rownames(features) <- NULL

    list(analysis = toupper(key), contrast = as.character(contrast),
         pathway_id = as.character(pathway_id), state = state,
         result = result, compounds = compounds, features = features)
}
