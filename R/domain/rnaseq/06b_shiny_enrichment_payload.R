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
#' only AFTER a successful lossless encode; the short `leading_edge` summary
#' string is kept. All significant rows are preserved (no top-N).
#' @noRd
.restructure_gsea <- function(pathway_res, gene_index) {
    idx <- pathway_res$enrichment_index
    pw  <- pathway_res$pathway_results
    gsea <- list(); le <- list()
    gidx <- idx[idx$analysis == "GSEA", , drop = FALSE]
    for (i in seq_len(nrow(gidx))) {
        r        <- gidx[i, ]
        db       <- r$database; rk <- r$group; contrast <- r$item
        tab      <- pw[[r$container]][[r$storage_key]]
        if (is.null(tab)) next
        if ("core_enrichment" %in% names(tab)) {
            ids <- as.character(tab$ID)
            enc <- lapply(tab$core_enrichment, function(s)
                .encode_genes(.split_genes(s), gene_index,
                              context = sprintf("GSEA leading-edge gene (%s/%s/%s)", db, rk, contrast)))
            names(enc) <- ids
            le[[db]][[rk]][[contrast]] <- enc
            tab$core_enrichment <- NULL   # strip only after successful encode
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
    pcut <- enr$pvalue_cutoff %||% 0.05
    list(
        # Configured run set -> lets Shiny distinguish "not run" (absent) databases.
        databases_run = as.character(enr$databases %||%
            unique(pathway_res$enrichment_manifest$database)),
        # Ranking methods present (internal keys; Shiny maps to display labels).
        rankings       = names(pathway_res$gsea_rankings %||% list()),
        pvalue_cutoff  = pcut,
        gsea_pvalue_cutoff = enr$gsea_pvalue_cutoff %||% pcut,
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

    gsea_split <- .restructure_gsea(pathway_res, gene_index)

    enrichment <- list(
        available          = TRUE,
        gene_index         = gene_index,
        config             = .build_enrichment_config(pathway_res, config),
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


#' Validate the compact enrichment payload structure
#'
#' Deeper, enrichment-specific structural validation kept OUT of the generic
#' Shiny contract (which only treats `enrichment` as an optional list). NULL is
#' valid (enrichment absent / not run). Checks the top-level schema, gene_index
#' uniqueness, the manifest columns (incl. the `status` three-state field), and
#' the compact-representation invariant that GSEA tables never carry the heavy
#' `core_enrichment` string (it lives in `gsea_leading_edge`).
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
        mcols <- c("analysis", "database", "group", "item", "evaluated",
                   "status", "n_significant", "has_simplify", "storage_key")
        if (!all(mcols %in% names(enrichment$manifest)))
            signal("enrichment$manifest missing required columns")
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
