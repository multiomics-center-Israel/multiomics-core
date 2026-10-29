# ============================================================
# Shiny Export — RNA-seq ENRICHMENT payload assembler (Stage 3B)
# ============================================================
# Builds the compact, portable `payload$enrichment` block from the Stage-3A
# enrichment primitives on `rna_pathway_res`. Representation transform only —
# NO biological result is recomputed (rankings, leading edge, membership, and
# result tables are carried through verbatim, only re-encoded as integer indices
# into a shared `gene_index`). See rnaseq-enrichment-shiny.md §14.
#
# This file is deliberately isolated from 06_shiny_export.R: Stage 3B builds and
# validates the block WITHOUT wiring it into build_shiny_payload_rnaseq() (that
# is Stage 3C). No plots/images/ggplot/pheatmap objects, no file paths, and no
# copies of expr_norm / de_stats / feature_annot enter this block.
# ============================================================

# Internal collection name (Stage-3A `group`) -> payload-facing ORA group name.
.ORA_GROUP_MAP <- c(
    contrasts              = "contrasts_with_direction",
    contrasts_wo_direction = "contrasts_without_direction",
    all_DE                 = "all_DE",
    partition              = "partition",
    binary_patterns        = "binary_patterns"
)
# ORA groups that hold exactly one table (stored directly, not keyed by contrast).
.ORA_SINGLE_GROUPS <- c("all_DE", "partition", "binary_patterns")


#' Encode a character vector of gene IDs to integer indices into gene_index
#'
#' Order- and duplicate-preserving. Any gene that cannot be resolved against
#' gene_index is a hard error (never silently dropped) — the union rule in
#' .build_gene_index() guarantees every referenced gene is present, so an
#' unresolved gene signals a real inconsistency.
#'
#' @param genes Character vector of gene IDs.
#' @param gene_index The shared gene universe (character).
#' @param context Short label used in the error message.
#' @return Integer vector of 1-based indices into gene_index.
#' @noRd
.encode_genes <- function(genes, gene_index, context = "gene reference") {
    if (length(genes) == 0) return(integer(0))
    idx <- match(as.character(genes), gene_index)
    if (anyNA(idx)) {
        missing <- unique(genes[is.na(idx)])
        stop(sprintf("[enrichment payload] %d unresolved %s not in gene_index (e.g. %s)",
                     length(missing), context,
                     paste(utils::head(missing, 5L), collapse = ", ")), call. = FALSE)
    }
    as.integer(idx)
}


#' Split a slash-separated gene string into a character vector (NA -> empty)
#' @noRd
.split_genes <- function(x) {
    if (length(x) != 1L || is.na(x) || !nzchar(x)) return(character(0))
    strsplit(x, "/", fixed = TRUE)[[1]]
}


#' Construct the shared, explicit gene_index (Stage 3B.2)
#'
#' Deterministic construction rule:
#'   gene_index = c( unique(gene_universe),
#'                   sort(unique( all enrichment-referenced genes NOT already in
#'                               gene_universe )) )
#' i.e. the measured expression universe first (in its native order, so the first
#' N indices align 1:1 with expr_norm rows), then every additional gene referenced
#' anywhere in enrichment (ranking names, GSEA leading edge, GSEA membership, ORA
#' geneID) appended in sorted order. Background/pathway genes that are NOT in
#' expr_norm are legitimately retained here (they are not dropped). The result is
#' validated to be unique.
#'
#' @param pathway_res rna_pathway_res (with Stage-3A siblings).
#' @param gene_universe Character vector of measured gene IDs (rownames(expr_norm));
#'   may be NULL/empty.
#' @return Character vector (the gene_index).
#' @noRd
.build_gene_index <- function(pathway_res, gene_universe = NULL) {
    parts <- list()
    add <- function(x) if (length(x)) parts[[length(parts) + 1L]] <<- as.character(x)

    # GSEA ranking names
    for (m in pathway_res$gsea_rankings %||% list())
        for (v in m) add(names(v))
    # GSEA membership genes
    for (db in pathway_res$pathway_membership %||% list())
        for (p in db) add(p)
    # GSEA leading edge (core_enrichment) + ORA geneID, from the result tables
    pw <- pathway_res$pathway_results %||% list()
    for (cont in setdiff(names(pw), "cluster_ora"))
        for (tab in pw[[cont]])
            if (!is.null(tab$core_enrichment))
                add(unlist(lapply(tab$core_enrichment, .split_genes), use.names = FALSE))
    if (!is.null(pw$cluster_ora))
        for (tab in pw$cluster_ora)
            if (!is.null(tab$geneID))
                add(unlist(lapply(tab$geneID, .split_genes), use.names = FALSE))

    referenced <- if (length(parts)) unique(unlist(parts, use.names = FALSE)) else character(0)
    universe   <- unique(as.character(gene_universe %||% character(0)))
    extra      <- sort(setdiff(referenced, universe))
    gene_index <- c(universe, extra)

    if (anyDuplicated(gene_index))
        stop("[enrichment payload] gene_index contains duplicates", call. = FALSE)
    gene_index
}


#' Restructure ORA tables into the nested payload schema via enrichment_index
#'
#' Uses the Stage-3A index (never parses concatenated keys). Preserves every
#' significant row and the `geneID` column. GO simplify tables are nested in
#' parallel; KEGG (has_simplify = FALSE) never gets a fabricated simplify entry.
#' @noRd
.restructure_ora <- function(pathway_res) {
    idx <- pathway_res$enrichment_index
    co  <- pathway_res$pathway_results$cluster_ora
    ora <- list()
    ora_idx <- idx[idx$analysis == "ORA", , drop = FALSE]
    for (i in seq_len(nrow(ora_idx))) {
        r    <- ora_idx[i, ]
        db   <- r$database
        pg   <- .ORA_GROUP_MAP[[r$group]] %||% r$group
        tab  <- co[[r$storage_key]]
        if (is.null(tab)) next
        if (pg %in% .ORA_SINGLE_GROUPS) {
            ora[[db]][[pg]] <- tab
        } else {
            ora[[db]][[pg]][[r$item]] <- tab
        }
        # Simplified GO result (only when the index says one is stored).
        if (!is.na(r$simplify_key)) {
            stab <- co[[r$simplify_key]]
            if (!is.null(stab)) {
                if (pg %in% .ORA_SINGLE_GROUPS) {
                    ora[[db]]$simplify[[pg]] <- stab
                } else {
                    ora[[db]]$simplify[[pg]][[r$item]] <- stab
                }
            }
        }
    }
    ora
}


#' Restructure GSEA tables + extract compact leading edge via enrichment_index
#'
#' Returns list(gsea, gsea_leading_edge). The heavy `core_enrichment` gene string
#' is encoded to integer indices (order-preserving) and REMOVED from each table
#' only AFTER a successful lossless encode, together with fgsea's comma-joined
#' `leadingEdge` copy of the same genes; the short `leading_edge` summary
#' string is kept. All significant rows are preserved (no top-N). The fgsea
#' engine stores every tested set, so tables are cut to the rows at
#' `sig_cutoff` here, which leaves the payload as it was with the
#' clusterProfiler engine (whose tables hold only those rows already).
#' @param sig_cutoff Adjusted-p cutoff (the configured GSEA cutoff).
#' @noRd
.restructure_gsea <- function(pathway_res, gene_index, sig_cutoff = 0.05) {
    idx <- pathway_res$enrichment_index
    pw  <- pathway_res$pathway_results
    gsea <- list(); le <- list()
    gidx <- idx[idx$analysis == "GSEA", , drop = FALSE]
    for (i in seq_len(nrow(gidx))) {
        r        <- gidx[i, ]
        db       <- r$database; rk <- r$group; contrast <- r$item
        tab      <- pw[[r$container]][[r$storage_key]]
        if (is.null(tab)) next
        if ("padj" %in% names(tab)) tab <- .gsea_significant_rows(tab, sig_cutoff)
        if ("core_enrichment" %in% names(tab)) {
            ids <- as.character(tab$ID)
            enc <- lapply(tab$core_enrichment, function(s)
                .encode_genes(.split_genes(s), gene_index,
                              context = sprintf("GSEA leading-edge gene (%s/%s/%s)", db, rk, contrast)))
            names(enc) <- ids
            le[[db]][[rk]][[contrast]] <- enc
            tab$core_enrichment <- NULL   # strip only after successful encode
            tab$leadingEdge     <- NULL
        }
        gsea[[db]][[rk]][[contrast]] <- tab
    }
    list(gsea = gsea, gsea_leading_edge = le)
}


#' Encode the exact GSEA rankings to compact indexed storage (Stage 3B.4)
#'
#' One entry per ranking_method x contrast (db-independent — never duplicated by
#' database). Each entry is list(idx = int[], score = num[]); round-trips exactly
#' to the original named numeric vector via setNames(score, gene_index[idx]).
#' @noRd
.encode_gsea_rankings <- function(pathway_res, gene_index) {
    out <- list()
    for (rk in names(pathway_res$gsea_rankings %||% list())) {
        for (contrast in names(pathway_res$gsea_rankings[[rk]])) {
            v <- pathway_res$gsea_rankings[[rk]][[contrast]]
            out[[rk]][[contrast]] <- list(
                idx   = .encode_genes(names(v), gene_index,
                                      context = sprintf("ranking gene (%s/%s)", rk, contrast)),
                score = unname(as.numeric(v)))
        }
    }
    out
}


#' Encode GSEA pathway membership to compact indexed storage (Stage 3B.5)
#'
#' database -> pathway_id -> int[] into gene_index. Round-trips exactly to the
#' original character membership via gene_index[idx].
#' @noRd
.encode_pathway_membership <- function(pathway_res, gene_index) {
    out <- list()
    for (db in names(pathway_res$pathway_membership %||% list())) {
        for (pid in names(pathway_res$pathway_membership[[db]])) {
            out[[db]][[pid]] <- .encode_genes(
                pathway_res$pathway_membership[[db]][[pid]], gene_index,
                context = sprintf("membership gene (%s/%s)", db, pid))
        }
    }
    out
}


#' Minimal enrichment config snapshot for the Shiny layer (Stage 3B.9)
#' @noRd
.build_enrichment_config <- function(pathway_res, config = NULL) {
    enr <- config$modes$rna$enrichment %||% list()
    cuts <- resolve_rna_gsea_settings(config)
    list(
        # Configured run set -> lets Shiny distinguish "not run" (absent) databases.
        databases_run = as.character(enr$databases %||%
            unique(pathway_res$enrichment_manifest$database)),
        # Ranking methods present (internal keys; Shiny maps to display labels).
        rankings       = names(pathway_res$gsea_rankings %||% list()),
        # The cutoffs the module ran ORA and filtered GSEA with (one resolver).
        pvalue_cutoff  = cuts$ora_cutoff,
        gsea_pvalue_cutoff = cuts$sig_cutoff,
        padj_method    = enr$padj_method %||% "fdr"
    )
}


#' Build the compact `payload$enrichment` block from rna_pathway_res
#'
#' Stage 3B assembler. Consumes the Stage-3A enrichment primitives and produces
#' the approved compact schema. Returns NULL when the input is absent or is not
#' the offline enrichment shape (e.g. the online workflow, or disabled/empty
#' enrichment) — the natural place for the `enrichment = NULL` decision; the
#' actual wiring into build_shiny_payload_rnaseq() is Stage 3C.
#'
#' @param pathway_res rna_pathway_res, with Stage-3A siblings
#'   (enrichment_manifest, enrichment_index, gsea_rankings, pathway_membership).
#' @param gene_universe Character vector of measured gene IDs (rownames(expr_norm)).
#' @param config Full pipeline config (for the minimal config snapshot); optional.
#' @return The `enrichment` list, or NULL.
#' @export
build_enrichment_payload <- function(pathway_res, gene_universe = NULL, config = NULL) {
    # ---- offline-shape / availability detection -> NULL fallback ----
    if (is.null(pathway_res)) return(NULL)
    if (is.null(pathway_res$enrichment_index) ||
        is.null(pathway_res$gsea_rankings) ||
        is.null(pathway_res$pathway_results) ||
        length(pathway_res$pathway_results) == 0) {
        return(NULL)
    }

    gene_index <- .build_gene_index(pathway_res, gene_universe)

    enr_config <- .build_enrichment_config(pathway_res, config)
    gsea_split <- .restructure_gsea(pathway_res, gene_index,
                                    sig_cutoff = enr_config$gsea_pvalue_cutoff)

    enrichment <- list(
        available          = TRUE,
        gene_index         = gene_index,
        config             = enr_config,
        manifest           = pathway_res$enrichment_manifest %||% .empty_enrichment_manifest(),
        ora                = .restructure_ora(pathway_res),
        gsea               = gsea_split$gsea,
        gsea_leading_edge  = gsea_split$gsea_leading_edge,
        gsea_rankings      = .encode_gsea_rankings(pathway_res, gene_index),
        pathway_membership = .encode_pathway_membership(pathway_res, gene_index)
    )
    # Self-validate the assembled block (structure + compact-representation
    # invariants) before it leaves the assembler.
    validate_enrichment_payload(enrichment, strict = TRUE)
    enrichment
}


#' Payload ORA group node name for a manifest (group, subgroup) — subgroup-aware
#'
#' §17 R4. The R2 manifest merges the two contrast collections under one
#' `group = "contrasts"`; the with-/without-direction distinction is carried by
#' `subgroup` (up/down/up_down vs all). This maps a manifest (group, subgroup) to
#' the payload `ora[[db]]` node key: contrasts + subgroup "all" ->
#' `contrasts_without_direction`; contrasts + up/down/up_down ->
#' `contrasts_with_direction`; every other group via `.ORA_GROUP_MAP`.
#' @noRd
.ora_payload_group_for <- function(group, subgroup) {
    if (identical(as.character(group), "contrasts"))
        return(if (identical(as.character(subgroup), "all"))
                   "contrasts_without_direction" else "contrasts_with_direction")
    .ORA_GROUP_MAP[[as.character(group)]] %||% as.character(group)
}

#' Cluster filter for one ORA manifest selection (generic, omics-agnostic)
#'
#' §17 R4. Returns the `Cluster` value to keep, or NA meaning "the whole table"
#' (a compare / up_down / all / all_DE view). Fully determined by the R2 manifest
#' columns — no new column, no biology:
#'   group_type "direction": subgroup up|down -> that Cluster; up_down|all -> whole
#'   group_type partition|pattern: item "compare" -> whole; else Cluster == item
#'   otherwise (all_DE / NA): whole
#' @noRd
.ora_selection_filter <- function(group_type, subgroup, item) {
    gt <- as.character(group_type); sg <- as.character(subgroup); it <- as.character(item)
    if (!is.na(gt) && gt == "direction")
        return(if (!is.na(sg) && sg %in% c("up", "down")) sg else NA_character_)
    if (!is.na(gt) && gt %in% c("partition", "pattern"))
        return(if (identical(it, "compare")) NA_character_ else it)
    NA_character_
}

#' Resolve one ORA manifest selection to its (filtered) compact-payload table
#'
#' §17 R4 — the executable resolution contract that lets the three-state selector
#' turn a manifest row (db, group, item, subgroup, group_type) into the exact ORA
#' table WITHOUT any duplicated per-subgroup storage: UP/DOWN filter the ONE
#' `contrasts_with_direction[[contrast]]` table on `Cluster`; UP+DOWN returns that
#' whole table; All returns the whole `contrasts_without_direction[[contrast]]`
#' table (distinct node — never a union); a single cluster/pattern filters the ONE
#' `partition`/`binary_patterns` table; Compare returns that whole table; all_DE
#' returns its single table. Generic and omics-agnostic (manifest columns only).
#'
#' @param enrichment The compact enrichment block.
#' @param db,group,item,subgroup,group_type Manifest-row fields (group_type may be
#'   NULL -> inferred from `group`).
#' @return The selected data.frame (possibly 0-row for an evaluated-but-empty
#'   subgroup), or NULL when the underlying node is absent (an entirely-empty
#'   evaluated unit that stored no table).
#' @export
resolve_ora_selection <- function(enrichment, db, group, item, subgroup = NA_character_,
                                  group_type = NULL) {
    ora_db <- enrichment$ora[[db]]
    if (is.null(ora_db)) return(NULL)
    pg   <- .ora_payload_group_for(group, subgroup)
    node <- if (pg %in% .ORA_SINGLE_GROUPS) ora_db[[pg]] else ora_db[[pg]][[item]]
    if (is.null(node)) return(NULL)
    # Infer group_type from `group` when the caller omits it (NULL) or passes NA,
    # so a partition/pattern selection never silently degrades to the whole table.
    gt <- if (is.null(group_type) || is.na(group_type))
              switch(as.character(group),
                     contrasts       = "direction",
                     partition       = "partition",
                     binary_patterns = "pattern",
                     NA_character_)
          else group_type
    filt <- .ora_selection_filter(gt, subgroup, item)
    if (is.na(filt) || !("Cluster" %in% names(node))) return(node)
    node[as.character(node$Cluster) == filt, , drop = FALSE]
}

#' Validate the compact enrichment payload structure
#'
#' Deeper, enrichment-specific structural validation kept OUT of the generic
#' Shiny contract (which only treats `enrichment` as an optional list). NULL is
#' valid (enrichment absent / not run). Checks the top-level schema, gene_index
#' uniqueness, the manifest columns (incl. the `status` three-state field and the
#' §17 R2 `subgroup`/`group_type`/`group_label` hierarchy), the §17 R4 resolution
#' invariant (every ORA manifest row resolves through resolve_ora_selection() to a
#' row count consistent with its three-state `status`/`n_significant`, proving
#' UP/DOWN/UP+DOWN/All + Compare all resolve from the compact tables with NO
#' duplication), and the compact-representation invariant that GSEA tables never
#' carry the heavy `core_enrichment` string (it lives in `gsea_leading_edge`).
#'
#' @param enrichment The enrichment block, or NULL.
#' @param strict If TRUE, violations stop(); otherwise warning().
#' @return invisibly TRUE.
#' @export
validate_enrichment_payload <- function(enrichment, strict = TRUE) {
    signal <- function(msg) if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
    if (is.null(enrichment)) return(invisible(TRUE))   # absent enrichment is valid

    if (!is.list(enrichment)) {
        signal("enrichment must be a list or NULL"); return(invisible(FALSE))
    }
    req <- c("available", "gene_index", "config", "manifest", "ora", "gsea",
             "gsea_leading_edge", "gsea_rankings", "pathway_membership")
    miss <- setdiff(req, names(enrichment))
    if (length(miss)) signal(paste("enrichment missing keys:", paste(miss, collapse = ", ")))

    if (!is.character(enrichment$gene_index) || anyDuplicated(enrichment$gene_index) != 0)
        signal("enrichment$gene_index must be a unique character vector")

    if (!is.data.frame(enrichment$manifest)) {
        signal("enrichment$manifest must be a data.frame")
    } else {
        m <- enrichment$manifest
        mcols <- c("analysis", "database", "group", "item", "subgroup",
                   "group_type", "group_label", "evaluated",
                   "status", "n_significant", "has_simplify", "storage_key")
        if (!all(mcols %in% names(m))) {
            signal("enrichment$manifest missing required columns")
        } else if (nrow(m) > 0) {
            # §17 R2 hierarchy vocabulary (NA allowed for GSEA / all_DE / single groups).
            ok_sg <- is.na(m$subgroup) | m$subgroup %in% c("up", "down", "up_down", "all")
            if (!all(ok_sg))
                signal("enrichment$manifest$subgroup outside {NA, up, down, up_down, all}")
            ok_gt <- is.na(m$group_type) | m$group_type %in% c("direction", "partition", "pattern")
            if (!all(ok_gt))
                signal("enrichment$manifest$group_type outside {NA, direction, partition, pattern}")

            # §17 R4 resolution invariant: every ORA manifest row must resolve
            # through resolve_ora_selection() to a row count consistent with its
            # three-state status — proving UP/DOWN/UP+DOWN/All + Compare all come
            # from the compact whole-table nodes (no duplicated per-subgroup tables).
            om <- m[m$analysis == "ORA", , drop = FALSE]
            for (i in seq_len(nrow(om))) {
                r   <- om[i, ]
                sel <- resolve_ora_selection(enrichment, r$database, r$group, r$item,
                                             r$subgroup, r$group_type)
                if (identical(r$status, "significant")) {
                    if (is.null(sel) || nrow(sel) != r$n_significant)
                        signal(sprintf(
                            "ORA row [%s/%s/%s/subgroup=%s] does not resolve to n_significant=%s rows",
                            r$database, r$group, r$item, r$subgroup, r$n_significant))
                } else if (identical(r$status, "empty")) {
                    # empty is fine either way: node absent (whole unit empty) or 0 filtered rows.
                    if (!is.null(sel) && nrow(sel) != 0L)
                        signal(sprintf(
                            "empty ORA row [%s/%s/%s/subgroup=%s] resolves to %d rows (expected 0)",
                            r$database, r$group, r$item, r$subgroup, nrow(sel)))
                }
            }
        }
    }

    # Compact-representation invariant: GSEA result tables must NOT carry the
    # heavy core_enrichment gene string (it is stored in gsea_leading_edge).
    for (db in enrichment$gsea %||% list())
        for (rk in db)
            for (tab in rk)
                if (is.data.frame(tab) && "core_enrichment" %in% names(tab))
                    signal("enrichment$gsea tables must not contain core_enrichment")

    invisible(TRUE)
}
