# R/domain/metabolomics/06r_enrichment_provenance.R
#
# Provenance of the metabolomics ORA and GSEA results.
#
# Both tests work on COMPOUNDS, and several measured features can map to one
# compound. The result tables therefore cannot say which features stood behind
# a compound, which feature represented it in the GSEA ranking, or which
# members of a pathway the test actually used. That information exists while
# the test runs and was dropped when it returned. The builders here keep it, as
# a `provenance` element beside the unchanged `table` (and `ranks`).
#
# Two rules:
#   1. Provenance never changes a statistic. The builders only read objects
#      the test already computed, and a failure in here cannot fail the test:
#      it yields a provenance with status "failed" and a warning.
#   2. Provenance records what was done; it does not work it out again. Where
#      a description has to be derived (aligning the members a GMT lists with
#      the members that matched), it is checked against the object the test
#      received, and a difference makes the status "inconsistent".
#
# `status` is always one of:
#   "ok"            built, and every consistency check passed;
#   "inconsistent"  built, but a check failed (see `diagnostics`);
#   "failed"        could not be built (see `error`).
# Code that needs provenance calls require_metab_enrichment_provenance(), which
# stops unless the status is "ok".
#
# Counting pathway members. A GMT set can list one compound twice, typically
# after HMDB ids are translated to KEGG. The two tests treat that differently,
# and `pathway_members` keeps both facts apart:
#   in_input_set            the occurrence was in the vector given to the test;
#   in_effective_statistic  the occurrence counted in the calculation.
# ORA counts every occurrence in the pathway size, so the two are equal. fgsea
# keeps each member once, so only the first occurrence is effective, and a set
# fgsea drops for its size has none. Neither behaviour is changed here.
#
# Reuses: %||%

.METAB_PROVENANCE_SCHEMA <- 1L


# ==== FAILURE HANDLING ========================================================

#' Mark a provenance step that could not be completed
#'
#' @param step Short name of the step.
#' @param e    The condition that stopped it.
#' @return An object of class \code{metab_provenance_error}: list(step, message).
metab_provenance_error <- function(step, e) {
    structure(list(step = step, message = conditionMessage(e)),
              class = "metab_provenance_error")
}

#' A provenance record for a build that failed
#'
#' @param method   "ora" or "gsea".
#' @param contrast Contrast name.
#' @param step     Step that failed.
#' @param message  What went wrong.
#' @return A provenance list with \code{status = "failed"}.
#' @noRd
.metab_provenance_failed <- function(method, contrast, step, message) {
    list(
        schema_version = .METAB_PROVENANCE_SCHEMA,
        method         = method,
        contrast       = contrast,
        status         = "failed",
        diagnostics    = data.frame(check = character(0), passed = logical(0),
                                    detail = character(0), stringsAsFactors = FALSE),
        error          = list(step = step, message = message)
    )
}

#' Build provenance without letting a failure reach the statistical result
#'
#' @param method   "ora" or "gsea".
#' @param contrast Contrast name.
#' @param expr     Expression that builds the provenance list.
#' @return The provenance, or a "failed" record (with a warning) on error.
#' @noRd
.metab_provenance_guard <- function(method, contrast, expr) {
    tryCatch(expr, error = function(e) {
        .metab_prov_warn("metabolomics ", toupper(method),
                         " provenance could not be built for '", contrast, "': ",
                         conditionMessage(e), ". The statistical result is ",
                         "unaffected; its provenance has status \"failed\".")
        .metab_provenance_failed(method, contrast, "build", conditionMessage(e))
    })
}

#' Warn about provenance without ever stopping the analysis
#'
#' Under \code{options(warn = 2)} a warning is raised as an error. Raised from
#' the provenance code, that error would end the ORA or GSEA call after its
#' statistics were computed. It is caught here instead: the warning is then
#' not shown, and the provenance record still carries the status and reason.
#'
#' @param ... Pieces of the message, pasted together.
#' @return NULL, invisibly.
#' @noRd
.metab_prov_warn <- function(...) {
    tryCatch(warning(..., call. = FALSE), error = function(e) NULL)
    invisible(NULL)
}

#' Checksum of one file, or NA
#'
#' Called from the ORA and GSEA functions beside each read of a reference
#' file, so it must never fail: anything that is not one existing file gives NA.
#'
#' @param path A file path (or NULL).
#' @return The md5 as a character scalar, or NA_character_.
#' @noRd
.metab_prov_md5 <- function(path) {
    tryCatch({
        p <- as.character(unlist(path, use.names = FALSE))
        if (length(p) != 1L || is.na(p) || !nzchar(p) || !file.exists(p) ||
            dir.exists(p)) {
            NA_character_
        } else {
            unname(tools::md5sum(p))
        }
    }, error = function(e) NA_character_)
}

#' One row of the diagnostics table
#'
#' @param check  Name of the check.
#' @param passed Whether it passed (anything but TRUE counts as failed).
#' @param detail Optional note for a failed check.
#' @return A one-row data.frame(check, passed, detail).
#' @noRd
.metab_prov_check <- function(check, passed, detail = "") {
    data.frame(check = check, passed = isTRUE(passed), detail = detail,
               stringsAsFactors = FALSE)
}

#' Set the status from the diagnostics and warn when a check failed
#'
#' @param prov   Provenance list without status/diagnostics.
#' @param checks List of one-row diagnostics data frames.
#' @return \code{prov} with \code{status}, \code{diagnostics} and \code{error}.
#' @noRd
.metab_prov_finish <- function(prov, checks) {
    diagnostics <- do.call(rbind, checks)
    rownames(diagnostics) <- NULL
    failed <- diagnostics$check[!diagnostics$passed]
    if (length(failed) > 0) {
        .metab_prov_warn("metabolomics ", toupper(prov$method), " provenance for '",
                         prov$contrast, "' does not agree with the result it ",
                         "describes (failed: ", paste(failed, collapse = ", "),
                         "). The statistical result is unaffected; its ",
                         "provenance has status \"inconsistent\".")
    }
    head <- list(
        schema_version = .METAB_PROVENANCE_SCHEMA,
        method         = prov$method,
        contrast       = prov$contrast,
        status         = if (length(failed) > 0) "inconsistent" else "ok",
        diagnostics    = diagnostics,
        error          = NULL
    )
    c(head, prov[setdiff(names(prov), c("method", "contrast"))])
}

#' Bind a list of data frames, or return an empty one
#'
#' @param rows  List of data frames (NULL entries are dropped).
#' @param empty Data frame to return when there is nothing to bind.
#' @return A data frame with fresh row names.
#' @noRd
.metab_prov_bind <- function(rows, empty) {
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) return(empty)
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}


# ==== FEATURES AND COMPOUNDS ===================================================

#' Describe how each measured feature got its compound id
#'
#' Called by \code{map_compounds_for_enrichment()} with the vectors it tracked
#' while choosing the ids. \code{id_source} is the step that made the final
#' choice, in that function's order of precedence: the \code{KEGG} column wins
#' over the HMDB -> KEGG mapping file, which wins over the name.
#'
#' \code{upstream_origin} says how far back the id can be traced.
#' \code{"not_tracked"} for an id taken from the \code{KEGG} column: the
#' enrichment stage reads that column as it finds it and cannot tell whether
#' the value came with the input annotation or was filled earlier from a
#' mapping file.
#'
#' @param feature_id  Feature ids, one per feature.
#' @param compound_id Final compound id per feature (NA when excluded).
#' @param valid       Logical: the feature has a usable compound id.
#' @param id_source   Step that set each id: "kegg_column", "hmdb_mapping",
#'   "name_column" or "matrix_rowname".
#' @param restricted  Logical: the id was dropped because it is not a KEGG id
#'   while the data carries KEGG ids.
#' @param kegg_cell   Raw content of the \code{KEGG} column (NA without one).
#' @return data.frame, one row per feature: feature_id, compound_id, id_source,
#'   upstream_origin, excluded_reason, kegg_cell, kegg_candidates,
#'   n_kegg_candidates, ambiguous.
describe_metab_feature_compounds <- function(feature_id, compound_id, valid,
                                             id_source, restricted, kegg_cell) {
    n <- length(compound_id)
    lens <- c(length(feature_id), length(valid), length(id_source),
              length(restricted), length(kegg_cell))
    if (is.null(feature_id) || any(lens != n)) {
        stop("describe_metab_feature_compounds(): the feature ids and the ",
             "tracked vectors do not have one entry per feature (", n,
             " compound ids; lengths ", paste(lens, collapse = ", "), ").")
    }

    # Every KEGG compound id in the cell, with the pattern the mapping uses.
    # The mapping takes the first; the rest are kept here so that a cell with
    # several ids is visible as ambiguous.
    cells <- as.character(kegg_cell)
    cells[is.na(cells)] <- ""
    candidates <- lapply(regmatches(cells, gregexpr("C[0-9]{5}", cells)), unique)
    n_candidates <- lengths(candidates)

    valid  <- as.logical(valid)
    source <- ifelse(valid, as.character(id_source), NA_character_)
    origin <- c(kegg_column    = "not_tracked",
                hmdb_mapping   = "enrichment_mapping_file",
                name_column    = "annotation",
                matrix_rowname = "annotation")[source]
    excluded <- ifelse(valid, NA_character_,
                       ifelse(as.logical(restricted), "not_kegg_annotated", "empty_id"))

    data.frame(
        feature_id        = as.character(feature_id),
        compound_id       = unname(ifelse(valid, as.character(compound_id), NA_character_)),
        id_source         = unname(source),
        upstream_origin   = unname(origin),
        excluded_reason   = unname(excluded),
        kegg_cell         = as.character(kegg_cell),
        kegg_candidates   = vapply(candidates, paste, character(1), collapse = ";"),
        n_kegg_candidates = as.integer(n_candidates),
        ambiguous         = n_candidates > 1L,
        stringsAsFactors  = FALSE,
        row.names         = NULL
    )
}

#' The feature record, or an error when it could not be built
#'
#' @param feature_compounds \code{map_compounds_for_enrichment()$feature_compounds}.
#' @return The data frame; stops otherwise.
#' @noRd
.metab_prov_feature_compounds <- function(feature_compounds) {
    if (inherits(feature_compounds, "metab_provenance_error")) {
        stop("the feature-to-compound record could not be built: ",
             feature_compounds$message)
    }
    if (!is.data.frame(feature_compounds)) {
        stop("the feature-to-compound record is missing")
    }
    feature_compounds
}


# ==== REFERENCE FILES ============================================================

#' The reference files a test was configured with, and their checksums
#'
#' \code{md5} is the checksum taken just before the file is read. For a GMT
#' library, \code{md5_after_read} is taken again just after it was parsed, so
#' the two bracket the read: when they agree, the parsed content is the content
#' with that checksum. \code{md5_at_end} is taken when the provenance is built.
#' Any difference among them means the file changed while the test ran, and the
#' provenance is then marked inconsistent: the result is real, but it cannot be
#' tied to one version of the file.
#'
#' The mapping table is read several times during a test (once for the
#' features, once per library), so it has no single read to bracket; its
#' checksum before the first read is compared with the one at the end.
#'
#' @param gmt_files        Resolved GMT paths, as the test received them.
#' @param mapping_file     Resolved HMDB -> KEGG path (or NULL).
#' @param gmt_capture      Per-library capture from the test's GMT loop, with
#'   the checksum of each library taken when it was read.
#' @param mapping_md5_read Checksum of the mapping file taken before its first
#'   read (NA when there is no such file).
#' @return data.frame(role, path, md5, md5_after_read, md5_at_end, read).
#'   \code{read} is FALSE, and \code{md5} NA, for a configured file the test
#'   did not read; \code{md5_after_read} is NA for the mapping table.
#' @noRd
.metab_prov_reference_files <- function(gmt_files, mapping_file, gmt_capture,
                                        mapping_md5_read) {
    gmt  <- as.character(unlist(gmt_files, use.names = FALSE))
    map  <- as.character(unlist(mapping_file, use.names = FALSE))
    path <- c(gmt, map)
    role <- c(rep("gmt", length(gmt)), rep("mapping", length(map)))
    ok   <- !is.na(path) & nzchar(path)
    path <- path[ok]
    role <- role[ok]

    read_file <- vapply(gmt_capture, function(g) as.character(g$file), character(1))
    read_md5  <- vapply(gmt_capture, function(g) {
        as.character(g$md5 %||% NA_character_)
    }, character(1))
    after_md5 <- vapply(gmt_capture, function(g) {
        as.character(g$md5_after_read %||% NA_character_)
    }, character(1))
    at  <- match(path, read_file)
    md5 <- ifelse(role == "gmt", read_md5[at],
                  as.character(mapping_md5_read %||% NA_character_)[1])
    read <- ifelse(role == "gmt", path %in% read_file, !is.na(md5))

    data.frame(role = role, path = path, md5 = as.character(md5),
               md5_after_read = as.character(ifelse(role == "gmt", after_md5[at],
                                                    NA_character_)),
               md5_at_end = vapply(path, .metab_prov_md5, character(1),
                                   USE.NAMES = FALSE),
               read = as.logical(read), stringsAsFactors = FALSE)
}

#' Whether every file that was read still has the content that was read
#'
#' @param reference_files Result of \code{.metab_prov_reference_files()}.
#' @return TRUE or FALSE.
#' @noRd
.metab_prov_references_unchanged <- function(reference_files) {
    used <- reference_files[reference_files$read, , drop = FALSE]
    # A library must have the same checksum before and after its read; the
    # mapping table has no after-read checksum to compare.
    bracketed <- ifelse(used$role == "gmt",
                        !is.na(used$md5_after_read) & used$md5 == used$md5_after_read,
                        TRUE)
    all(!is.na(used$md5) & !is.na(used$md5_at_end) & used$md5 == used$md5_at_end &
            bracketed)
}


# ==== PATHWAYS =====================================================================

#' The library entry a loaded pathway label came from
#'
#' Libraries are merged by label, and an equal label in a later library
#' replaces the earlier one, so the last library holding the label wins.
#'
#' @param gmt_capture Per-library capture from the test's GMT loop.
#' @param label       Pathway label.
#' @return list(library, pathway_id, listed, translated), or NULL.
#' @noRd
.metab_prov_winner <- function(gmt_capture, label) {
    for (i in rev(seq_along(gmt_capture))) {
        g <- gmt_capture[[i]]
        k <- match(label, g$labels)
        if (!is.na(k)) {
            return(list(library    = g$library,
                        pathway_id = g$ids[[k]],
                        listed     = as.character(unname(g$listed[[k]])),
                        translated = as.character(unname(g$translated[[k]]))))
        }
    }
    NULL
}

#' What became of every pathway line of every library
#'
#' One row per accepted GMT line, so that a line replaced by a later line with
#' the same id, or by an equal label in a later library, is still on record.
#' Two different ids of one file can also give the same label ("A" described as
#' "b - c", and "A - b" described as "c"). The merge then keeps the first of
#' them, and the others are recorded as replaced.
#'
#' @param gmt_capture Per-library capture from the test's GMT loop.
#' @param sets_loaded The merged sets, before the size filter.
#' @param keep        Logical, aligned with \code{sets_loaded}: passed the filter.
#' @param returned    Labels the test returned a result for; NULL when every
#'   kept set is tested (ORA).
#' @return data.frame(gmt_file, library, line, pathway_id, description, pathway,
#'   n_listed, fate).
#' @noRd
.metab_prov_pathway_sources <- function(gmt_capture, sets_loaded, keep, returned = NULL) {
    empty <- data.frame(gmt_file = character(0), library = character(0),
                        line = integer(0), pathway_id = character(0),
                        description = character(0), pathway = character(0),
                        n_listed = integer(0), fate = character(0),
                        stringsAsFactors = FALSE)
    all_labels <- lapply(gmt_capture, function(g) g$labels)
    rows <- list()
    for (i in seq_along(gmt_capture)) {
        g <- gmt_capture[[i]]
        later   <- unlist(all_labels[seq_along(all_labels) > i], use.names = FALSE)
        rec_ids <- vapply(g$records, function(r) r$pathway_id, character(1))
        for (j in seq_along(g$records)) {
            r     <- g$records[[j]]
            label <- g$labels[match(r$pathway_id, g$ids)]
            pos   <- match(label, names(sets_loaded))
            fate <- if (any(rec_ids[seq_along(rec_ids) > j] == r$pathway_id)) {
                "replaced_same_id_in_file"
            } else if (!identical(g$ids[[match(label, g$labels)]], r$pathway_id)) {
                # Another id of this file has the same label and comes first:
                # the merge looks sets up by label and finds that one.
                "replaced_same_label_in_file"
            } else if (label %in% later) {
                "replaced_same_label_later_library"
            } else if (is.na(pos)) {
                "not_loaded"
            } else if (!isTRUE(unname(keep[pos]))) {
                "too_few_matches"
            } else if (!is.null(returned) && !(label %in% returned)) {
                "not_returned"
            } else {
                "tested"
            }
            rows[[length(rows) + 1L]] <- data.frame(
                gmt_file = g$file, library = g$library, line = as.integer(r$line),
                pathway_id = r$pathway_id, description = r$description,
                pathway = label, n_listed = length(r$members), fate = fate,
                stringsAsFactors = FALSE)
        }
    }
    .metab_prov_bind(rows, empty)
}

#' Check that the captured libraries describe the sets that were loaded
#'
#' @param gmt_capture Per-library capture from the test's GMT loop.
#' @param sets_loaded The merged sets, before the size filter.
#' @return TRUE when every loaded set equals its winning library entry.
#' @noRd
.metab_prov_sources_match <- function(gmt_capture, sets_loaded) {
    all(vapply(seq_along(sets_loaded), function(i) {
        win <- .metab_prov_winner(gmt_capture, names(sets_loaded)[[i]])
        !is.null(win) &&
            identical(win$translated, as.character(unname(sets_loaded[[i]]))) &&
            identical(win$library, attr(sets_loaded[[i]], "library"))
    }, logical(1)))
}

#' The members of every set given to the test, one row per listed occurrence
#'
#' Duplicates are kept: a compound listed twice is two rows, numbered by
#' \code{occurrence}. The match to the universe repeats the test's own lookup,
#' and is then compared with the vector the test received.
#'
#' @param gmt_capture Per-library capture from the test's GMT loop.
#' @param sets_kept   The sets that passed the size filter, as listed.
#' @param sets_input  The same sets as given to the test (matched members only).
#' @param lookup      The test's lower-case lookup of its universe.
#' @return list(members, inputs_match): the data frame (pathway, position,
#'   member_as_listed, member_after_translation, compound_id, occurrence,
#'   in_input_set) and whether it reproduces \code{sets_input} exactly.
#' @noRd
.metab_prov_pathway_members <- function(gmt_capture, sets_kept, sets_input, lookup) {
    empty <- data.frame(pathway = character(0), position = integer(0),
                        member_as_listed = character(0),
                        member_after_translation = character(0),
                        compound_id = character(0), occurrence = integer(0),
                        in_input_set = logical(0), stringsAsFactors = FALSE)
    rows <- vector("list", length(sets_input))
    inputs_match <- TRUE
    for (i in seq_along(sets_input)) {
        pw         <- names(sets_input)[[i]]
        translated <- as.character(unname(sets_kept[[i]]))
        win        <- .metab_prov_winner(gmt_capture, pw)
        listed <- if (!is.null(win) && length(win$listed) == length(translated)) {
            win$listed
        } else {
            rep(NA_character_, length(translated))
        }
        matched  <- as.character(unname(lookup[tolower(translated)]))
        in_input <- !is.na(matched)
        if (!identical(matched[in_input], as.character(unname(sets_input[[i]])))) {
            inputs_match <- FALSE
        }
        occurrence <- rep(NA_integer_, length(translated))
        if (any(in_input)) {
            occurrence[in_input] <- as.integer(stats::ave(
                seq_len(sum(in_input)), matched[in_input], FUN = seq_along))
        }
        rows[[i]] <- data.frame(
            pathway = rep(pw, length(translated)),
            position = seq_along(translated),
            member_as_listed = listed, member_after_translation = translated,
            compound_id = matched, occurrence = occurrence,
            in_input_set = in_input, stringsAsFactors = FALSE)
    }
    list(members = .metab_prov_bind(rows, empty), inputs_match = inputs_match)
}


# ==== ORA ===========================================================================

#' Provenance of one ORA result
#'
#' Called at the end of \code{run_metabolomics_ora()} with the objects it
#' computed. Counts stay compound-based: \code{compounds} carries the universe
#' and the foreground, and \code{features} says, separately, whether a feature
#' is significant and whether its compound is in the universe or the
#' foreground. A feature that is not significant can belong to a foreground
#' compound, because another feature of that compound is.
#'
#' @param contrast          Contrast name.
#' @param feature_compounds \code{map_compounds_for_enrichment()$feature_compounds}.
#' @param universe          Background compounds (\code{bg_ids}).
#' @param foreground        Foreground compounds (\code{fg_ids}).
#' @param de_tbl            The contrast's DE table.
#' @param sig_mask          Logical per DE row: the feature is significant.
#' @param thresholds        list(pvalue_column, p_cutoff, logfc_cutoff).
#' @param gmt_capture       Per-library capture from the GMT loop.
#' @param sets_loaded       Merged sets before the size filter.
#' @param keep              Logical per loaded set: passed the filter.
#' @param sets_kept         Kept sets, as listed.
#' @param sets_input        Kept sets as given to Fisher's test.
#' @param lookup            Lower-case lookup of the universe (\code{bg_map}).
#' @param result_table      The result table.
#' @param gmt_files,mapping_file Resolved reference file paths.
#' @param mapping_md5_read  Checksum of the mapping file before its first read.
#' @return A provenance list (see the file header); never an error.
build_metab_ora_provenance <- function(contrast, feature_compounds,
                                       universe, foreground, de_tbl, sig_mask,
                                       thresholds, gmt_capture, sets_loaded, keep,
                                       sets_kept, sets_input, lookup, result_table,
                                       gmt_files, mapping_file,
                                       mapping_md5_read = NA_character_) {
    .metab_provenance_guard("ora", contrast, {
        fc      <- .metab_prov_feature_compounds(feature_compounds)
        refs    <- .metab_prov_reference_files(gmt_files, mapping_file, gmt_capture,
                                               mapping_md5_read)
        sources <- .metab_prov_pathway_sources(gmt_capture, sets_loaded, keep)
        mem     <- .metab_prov_pathway_members(gmt_capture, sets_kept, sets_input, lookup)

        # Fisher's test takes the pathway size as the length of the vector it
        # is given, so every occurrence counts.
        members <- mem$members
        members$in_effective_statistic <- members$in_input_set
        members$in_foreground <- members$in_input_set &
            members$compound_id %in% foreground

        # Features: every measured feature, and every feature of the DE table.
        de_feature   <- as.character(de_tbl$feature_id)
        all_features <- unique(c(fc$feature_id, de_feature))
        compound     <- fc$compound_id[match(all_features, fc$feature_id)]
        de_pos       <- match(all_features, de_feature)
        features <- data.frame(
            feature_id             = all_features,
            compound_id            = compound,
            in_de_table            = !is.na(de_pos),
            significant            = !is.na(de_pos) & as.logical(sig_mask)[de_pos],
            compound_in_universe   = !is.na(compound) & compound %in% universe,
            compound_in_foreground = !is.na(compound) & compound %in% foreground,
            stringsAsFactors       = FALSE)

        ids <- unique(c(as.character(universe), stats::na.omit(fc$compound_id)))
        sig_compounds <- features$compound_id[features$significant]
        compounds <- data.frame(
            compound_id            = ids,
            in_universe            = ids %in% universe,
            in_foreground          = ids %in% foreground,
            n_features             = as.integer(table(factor(fc$compound_id, levels = ids))),
            n_significant_features = as.integer(table(factor(sig_compounds, levels = ids))),
            stringsAsFactors       = FALSE)

        # ---- consistency with the result table ----
        per_pathway <- function(fun) {
            all(vapply(seq_len(nrow(result_table)), function(r) {
                rows <- members[members$pathway == result_table$pathway[[r]] &
                                    members$in_effective_statistic, , drop = FALSE]
                isTRUE(fun(result_table[r, , drop = FALSE], rows))
            }, logical(1)))
        }
        checks <- list(
            .metab_prov_check("universe_size",
                all(result_table$bg_size == length(universe))),
            .metab_prov_check("universe_is_measured_compounds",
                setequal(universe, stats::na.omit(fc$compound_id))),
            .metab_prov_check("foreground_size",
                all(result_table$fg_size == length(foreground))),
            .metab_prov_check("foreground_from_significant_features",
                setequal(foreground, unique(features$compound_id[
                    features$significant & features$compound_in_universe]))),
            .metab_prov_check("reference_files_unchanged",
                .metab_prov_references_unchanged(refs)),
            .metab_prov_check("library_sources",
                .metab_prov_sources_match(gmt_capture, sets_loaded) &&
                    !any(sources$fate == "not_loaded")),
            .metab_prov_check("pathway_inputs",
                mem$inputs_match &&
                    all(result_table$pathway %in% names(sets_input))),
            # Existing behaviour: the size counts every listed occurrence.
            .metab_prov_check("pathway_size",
                per_pathway(function(res, rows) res$pathway_size == nrow(rows))),
            # ... while the overlap counts each foreground compound once.
            .metab_prov_check("overlap",
                per_pathway(function(res, rows) {
                    res$overlap == length(unique(rows$compound_id[rows$in_foreground]))
                })),
            .metab_prov_check("overlap_compounds",
                per_pathway(function(res, rows) {
                    identical(as.character(res$overlap_genes),
                              paste(foreground[foreground %in% rows$compound_id],
                                    collapse = ";"))
                }))
        )

        .metab_prov_finish(list(
            method            = "ora",
            contrast          = contrast,
            reference_files   = refs,
            feature_compounds = fc,
            pathway_sources   = sources,
            pathway_members   = members,
            thresholds        = thresholds,
            compounds         = compounds,
            features          = features
        ), checks)
    })
}


# ==== GSEA ==========================================================================

#' Provenance of one GSEA result
#'
#' Called at the end of \code{run_metabolomics_gsea()} with the objects it
#' computed. \code{ranked} follows the final \code{ranks} vector entry by entry
#' and names the feature that represents each compound: the one the
#' deduplication kept, carried through it by \code{rank_feature}, not chosen
#' again here. \code{features} gives every feature its role, so that features
#' which did not enter the ranking stay visible without being counted.
#'
#' \code{pathway_members} separates the occurrences given to fgsea
#' (\code{in_input_set}) from the members fgsea used
#' (\code{in_effective_statistic}): fgsea keeps each member once and drops a
#' set that is then too small. \code{measured} marks members with measured
#' features whether or not they were ranked; it is for display only and never
#' implies a contribution to ES, NES or the leading edge.
#'
#' @param contrast          Contrast name.
#' @param feature_compounds \code{map_compounds_for_enrichment()$feature_compounds}.
#' @param de_tbl_mapped     The DE table with \code{compound_id}, before filtering.
#' @param de_tbl_ranked     The rows that entered the ranking (compound and logFC).
#' @param ranks             The final ranking (logFC named by compound).
#' @param rank_feature      Feature id named by compound, deduplicated with
#'   \code{ranks}; NULL when it could not be captured.
#' @param gmt_capture       Per-library capture from the GMT loop.
#' @param sets_loaded       Merged sets before the size filter.
#' @param keep              Logical per loaded set: passed the filter.
#' @param sets_kept         Kept sets, as listed.
#' @param sets_input        Kept sets as given to fgsea.
#' @param lookup            Lower-case lookup of the ranked compounds.
#' @param fgsea_res         fgsea's result.
#' @param result_table      The result table.
#' @param gmt_files,mapping_file Resolved reference file paths.
#' @param mapping_md5_read  Checksum of the mapping file before its first read.
#' @return A provenance list (see the file header); never an error.
build_metab_gsea_provenance <- function(contrast, feature_compounds, de_tbl_mapped,
                                        de_tbl_ranked, ranks, rank_feature,
                                        gmt_capture, sets_loaded, keep, sets_kept,
                                        sets_input, lookup, fgsea_res, result_table,
                                        gmt_files, mapping_file,
                                        mapping_md5_read = NA_character_) {
    .metab_provenance_guard("gsea", contrast, {
        fc       <- .metab_prov_feature_compounds(feature_compounds)
        refs     <- .metab_prov_reference_files(gmt_files, mapping_file, gmt_capture,
                                                mapping_md5_read)
        returned <- as.character(fgsea_res$pathway)
        sources  <- .metab_prov_pathway_sources(gmt_capture, sets_loaded, keep,
                                                returned = returned)
        mem      <- .metab_prov_pathway_members(gmt_capture, sets_kept, sets_input, lookup)

        # ---- the ranking and its representative features ----
        compound_ids <- names(ranks)
        if (is.null(rank_feature)) {
            stop("the feature behind each ranked compound was not captured")
        }
        if (anyDuplicated(compound_ids) > 0 || anyDuplicated(names(rank_feature)) > 0 ||
            !setequal(compound_ids, names(rank_feature))) {
            stop("the captured representative features do not line up with the ",
                 "ranked compounds")
        }
        representative <- unname(rank_feature[compound_ids])

        cand_compound <- as.character(de_tbl_ranked$compound_id)
        cand_logfc    <- de_tbl_ranked$logFC
        max_abs  <- tapply(abs(cand_logfc), cand_compound, max)
        at_max   <- abs(cand_logfc) == unname(max_abs[cand_compound])
        n_cand   <- tapply(cand_logfc, cand_compound, length)
        n_at_max <- tapply(at_max, cand_compound, sum)
        opposite <- tapply(seq_along(cand_logfc), cand_compound, function(i) {
            length(unique(sign(cand_logfc[i][at_max[i]]))) > 1L
        })
        ranked <- data.frame(
            rank                   = seq_along(ranks),
            compound_id            = compound_ids,
            logFC                  = unname(ranks),
            representative_feature = representative,
            n_features             = as.integer(n_cand[compound_ids]),
            tied_abs_logfc         = as.integer(n_at_max[compound_ids]) > 1L,
            tied_opposite_sign     = as.logical(opposite[compound_ids]),
            stringsAsFactors       = FALSE)

        # ---- every feature's role ----
        de_feature  <- as.character(de_tbl_mapped$feature_id)
        de_compound <- as.character(unname(de_tbl_mapped$compound_id))
        de_logfc    <- de_tbl_mapped$logFC
        has_compound <- !is.na(de_compound) & nzchar(de_compound)
        usable       <- has_compound & !is.na(de_logfc)
        is_rep       <- usable & de_feature == unname(rank_feature[de_compound])
        is_rep[is.na(is_rep)] <- FALSE
        role <- rep("not_representative", length(de_feature))
        role[!has_compound]                  <- "no_compound"
        role[has_compound & is.na(de_logfc)] <- "missing_logfc"
        role[is_rep]                         <- "representative"
        de_compound[!has_compound] <- NA_character_

        unranked <- fc[!(fc$feature_id %in% de_feature), , drop = FALSE]
        features <- data.frame(
            feature_id  = c(de_feature, unranked$feature_id),
            compound_id = c(de_compound, unranked$compound_id),
            in_de_table = c(rep(TRUE, length(de_feature)), rep(FALSE, nrow(unranked))),
            logFC       = c(de_logfc, rep(NA_real_, nrow(unranked))),
            role        = c(role, rep("not_in_de_table", nrow(unranked))),
            stringsAsFactors = FALSE)

        # ---- pathway members: given to fgsea, used by fgsea, merely measured ----
        members <- mem$members
        # fgsea matches each set to the ranking and keeps a member once, then
        # drops a set outside its size limits; such a set has no result.
        members$in_effective_statistic <- members$in_input_set &
            members$occurrence == 1L & members$pathway %in% returned
        members$in_effective_statistic[is.na(members$in_effective_statistic)] <- FALSE

        measured_ids    <- unique(stats::na.omit(fc$compound_id))
        measured_lookup <- stats::setNames(measured_ids, tolower(measured_ids))
        members$measured_compound_id <- as.character(unname(
            measured_lookup[tolower(members$member_after_translation)]))
        members$measured <- !is.na(members$measured_compound_id)

        leading <- stats::setNames(lapply(fgsea_res$leadingEdge, as.character), returned)
        members$in_leading_edge <- FALSE
        for (pw in unique(members$pathway)) {
            idx <- which(members$pathway == pw)
            members$in_leading_edge[idx] <- members$in_effective_statistic[idx] &
                members$compound_id[idx] %in% leading[[pw]]
        }

        # ---- consistency with the result ----
        rep_row <- match(paste(representative, compound_ids, sep = "\r"),
                         paste(as.character(de_tbl_ranked$feature_id), cand_compound,
                               sep = "\r"))
        per_pathway <- function(fun) {
            all(vapply(seq_len(nrow(result_table)), function(r) {
                pw   <- result_table$pathway[[r]]
                rows <- members[members$pathway == pw, , drop = FALSE]
                isTRUE(fun(result_table[r, , drop = FALSE], rows, leading[[pw]]))
            }, logical(1)))
        }
        checks <- list(
            .metab_prov_check("ranked_follows_ranks",
                identical(ranked$compound_id, names(ranks)) &&
                    identical(ranked$logFC, unname(ranks))),
            .metab_prov_check("ranked_universe",
                sum(usable) == nrow(de_tbl_ranked) &&
                    setequal(unique(cand_compound), compound_ids) &&
                    sum(is_rep) == length(ranks)),
            .metab_prov_check("representative_logfc",
                !anyNA(rep_row) && all(cand_logfc[rep_row] == unname(ranks))),
            .metab_prov_check("representative_is_largest_abs_logfc",
                all(abs(unname(ranks)) == unname(max_abs[compound_ids]))),
            .metab_prov_check("reference_files_unchanged",
                .metab_prov_references_unchanged(refs)),
            .metab_prov_check("library_sources",
                .metab_prov_sources_match(gmt_capture, sets_loaded) &&
                    !any(sources$fate == "not_loaded")),
            .metab_prov_check("pathway_inputs",
                mem$inputs_match && all(returned %in% names(sets_input)) &&
                    setequal(result_table$pathway, returned)),
            # Existing behaviour: fgsea's size counts each member once.
            .metab_prov_check("pathway_size",
                per_pathway(function(res, rows, le) {
                    res$pathway_size == sum(rows$in_effective_statistic)
                })),
            .metab_prov_check("leading_edge",
                per_pathway(function(res, rows, le) {
                    identical(as.character(res$leading_edge),
                              paste(le, collapse = ";")) &&
                        setequal(le, rows$compound_id[rows$in_leading_edge])
                }))
        )

        .metab_prov_finish(list(
            method            = "gsea",
            contrast          = contrast,
            reference_files   = refs,
            feature_compounds = fc,
            pathway_sources   = sources,
            pathway_members   = members,
            ranked            = ranked,
            features          = features
        ), checks)
    })
}


# ==== READING PROVENANCE ============================================================

#' What is wrong with a provenance record
#'
#' @param provenance The \code{provenance} element of an ORA or GSEA result.
#' @return Character vector of problems; empty when the record is usable.
validate_metab_enrichment_provenance <- function(provenance) {
    if (!is.list(provenance) || is.null(provenance$status)) {
        return("there is no provenance record")
    }
    if (!identical(provenance$schema_version, .METAB_PROVENANCE_SCHEMA)) {
        return(paste0("its schema version is ",
                      paste(format(provenance$schema_version), collapse = ""),
                      ", and this code reads version ", .METAB_PROVENANCE_SCHEMA))
    }
    if (identical(provenance$status, "failed")) {
        return(paste0("it could not be built: ",
                      provenance$error$message %||% "no message recorded"))
    }
    if (identical(provenance$status, "inconsistent")) {
        failed <- provenance$diagnostics$check[!provenance$diagnostics$passed]
        return(paste0("it does not agree with the result it describes (failed: ",
                      paste(failed, collapse = ", "), ")"))
    }
    if (!identical(provenance$status, "ok")) {
        return(paste0("its status is '", provenance$status, "'"))
    }
    needed <- c("reference_files", "feature_compounds", "pathway_sources",
                "pathway_members",
                switch(provenance$method %||% "",
                       ora  = c("thresholds", "compounds", "features"),
                       gsea = c("ranked", "features"),
                       "a method of 'ora' or 'gsea'"))
    missing <- setdiff(needed, names(provenance))
    if (length(missing) > 0) {
        return(paste0("it lacks: ", paste(missing, collapse = ", ")))
    }
    character(0)
}

#' Require usable provenance on an ORA or GSEA result
#'
#' For code that cannot do its job without provenance. A result whose
#' provenance is missing, failed or inconsistent stops here with the reason,
#' instead of being read as if it had none to give.
#'
#' @param res   Result of \code{run_metabolomics_ora()} or
#'   \code{run_metabolomics_gsea()}.
#' @param label Optional description of the result for the error message.
#' @return The provenance, invisibly; stops when it is not usable.
require_metab_enrichment_provenance <- function(res, label = NULL) {
    provenance <- if (is.list(res)) res$provenance else NULL
    problems <- validate_metab_enrichment_provenance(provenance)
    if (length(problems) > 0) {
        what <- label %||% paste0(
            "metabolomics ", toupper(provenance$method %||% "enrichment"),
            if (!is.null(provenance$contrast)) paste0(" '", provenance$contrast, "'"))
        stop("Provenance is required for ", what, ", but ",
             paste(problems, collapse = "; "),
             ". The statistical result itself is unaffected.", call. = FALSE)
    }
    invisible(provenance)
}

#' Features behind each pathway, with their role in the test
#'
#' Joins the members of every set given to the test to all measured features
#' of those compounds. It is meant for display, so it lists more than the test
#' used, and says so row by row:
#' \itemize{
#'   \item ORA: \code{compound_in_universe}, \code{feature_significant},
#'     \code{compound_in_foreground}, \code{compound_in_overlap}. The four are
#'     independent; the test counted compounds, not these rows.
#'   \item GSEA: \code{compound_ranked}, \code{feature_is_representative},
#'     \code{contributes_to_score} (both, in a set fgsea returned) and
#'     \code{in_leading_edge}. A feature that did not represent its compound,
#'     or a compound that was not ranked, never contributes and is never in the
#'     leading edge.
#' }
#'
#' @param provenance The \code{provenance} element of an ORA or GSEA result.
#' @return data.frame, one row per pathway x compound x feature: pathway,
#'   pathway_id, library, compound_id, feature_id, member_occurrences,
#'   in_effective_statistic, and the method's role flags.
metab_pathway_feature_membership <- function(provenance) {
    problems <- validate_metab_enrichment_provenance(provenance)
    if (length(problems) > 0) {
        stop("metab_pathway_feature_membership(): the provenance cannot be used: ",
             paste(problems, collapse = "; "), ".", call. = FALSE)
    }
    is_gsea <- identical(provenance$method, "gsea")
    members <- provenance$pathway_members
    feats   <- provenance$features

    # One compound per listed member: the one the test matched, or, for GSEA,
    # a measured compound the ranking did not contain.
    compound <- members$compound_id
    if (is_gsea) {
        compound <- ifelse(members$in_input_set, members$compound_id,
                           members$measured_compound_id)
    }
    members <- members[!is.na(compound), , drop = FALSE]
    compound <- compound[!is.na(compound)]

    key   <- paste(members$pathway, compound, sep = "\r")
    group <- factor(key, levels = unique(key))
    first <- !duplicated(key)
    any_by <- function(x) as.logical(tapply(x, group, any))
    pc <- data.frame(
        pathway                = members$pathway[first],
        compound_id            = compound[first],
        member_occurrences     = as.integer(table(group)),
        in_input_set           = any_by(members$in_input_set),
        in_effective_statistic = any_by(members$in_effective_statistic),
        stringsAsFactors       = FALSE)
    if (is_gsea) {
        pc$in_leading_edge <- any_by(members$in_leading_edge)
    } else {
        pc$in_foreground <- any_by(members$in_foreground)
    }

    src <- provenance$pathway_sources
    src <- src[src$fate %in% c("tested", "not_returned"), , drop = FALSE]
    src_row <- match(pc$pathway, src$pathway)

    # Every measured feature of each pathway compound.
    by_compound <- split(seq_len(nrow(feats)), feats$compound_id)
    hits <- unname(by_compound[pc$compound_id])
    p <- rep(seq_len(nrow(pc)), lengths(hits))
    f <- as.integer(unlist(hits, use.names = FALSE))

    out <- data.frame(
        pathway                = pc$pathway[p],
        pathway_id             = src$pathway_id[src_row][p],
        library                = src$library[src_row][p],
        compound_id            = pc$compound_id[p],
        feature_id             = feats$feature_id[f],
        member_occurrences     = pc$member_occurrences[p],
        in_effective_statistic = pc$in_effective_statistic[p],
        stringsAsFactors       = FALSE)
    if (is_gsea) {
        representative <- feats$role[f] == "representative"
        out$compound_ranked           <- pc$in_input_set[p]
        out$feature_is_representative <- representative
        out$contributes_to_score      <- pc$in_effective_statistic[p] & representative
        out$in_leading_edge           <- pc$in_leading_edge[p] & representative
    } else {
        out$compound_in_universe   <- pc$in_input_set[p]
        out$feature_significant    <- feats$significant[f]
        out$compound_in_foreground <- pc$in_foreground[p]
        out$compound_in_overlap    <- pc$in_effective_statistic[p] & pc$in_foreground[p]
    }
    rownames(out) <- NULL
    out
}
