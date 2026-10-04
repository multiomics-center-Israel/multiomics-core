# =============================================================================
# Loadings GSEA: pathway enrichment on DIABLO loadings and MOFA2 weights
# =============================================================================
#
# The ORA path in 07_enrichment.R takes the top 50 features of a component by
# absolute loading and asks which pathways they crowd. That throws the sign
# away and makes the answer depend on an arbitrary cut. GSEA instead ranks every
# feature of the view by its signed loading and asks whether a pathway's members
# sit towards one end of that ranking, so each pathway gets a direction.
#
# The gene sets are the ones the ORA path already tests against, reached the
# same way: KEGG through clusterProfiler for organisms with a KEGG code and an
# OrgDb, the configured GMTs otherwise, and KEGG compound pathways for
# metabolites. Only the test changes.
#
# Which test runs is a config choice, `modes.multiomics.enrichment.loadings
# .method`: "gsea" (the default) or "ora" (the earlier behaviour, unchanged).


#' Which enrichment test the loadings get
#'
#' @param config Full config object.
#' @return "gsea" or "ora". Absent means "gsea"; the config validator rejects
#'   any other value, so the warning here only fires for a config that
#'   bypassed it.
#' @examples
#' loadings_enrichment_method(list())  # "gsea"
loadings_enrichment_method <- function(config) {
    method <- config$modes$multiomics$enrichment$loadings$method %||% "gsea"
    method <- tolower(trimws(as.character(method)[1]))
    if (!method %in% c("gsea", "ora")) {
        warning("modes.multiomics.enrichment.loadings.method must be \"gsea\" or ",
                "\"ora\", not \"", method, "\"; using \"gsea\".", call. = FALSE)
        method <- "gsea"
    }
    method
}


#' Record which test produced the loadings enrichment in this directory
#'
#' The report reads the figures by glob, and both tests write into the same two
#' directories. Without a record, a run that switched from one test to the other
#' would have the report show the previous test's files as this run's.
#'
#' @param out_dir The loadings-enrichment output directory.
#' @param method "gsea" or "ora".
#' @return Invisibly, the path written.
record_loadings_method <- function(out_dir, method) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(out_dir, "loadings_enrichment_method.txt")
    writeLines(method, path)
    invisible(path)
}


#' Size bounds and seed for the loadings GSEA
#'
#' @param config Full config object.
#' @return List with `min_size`, `max_size` and `seed`.
#' @keywords internal
.loadings_gsea_settings <- function(config) {
    lcfg <- config$modes$multiomics$enrichment$loadings %||% list()
    list(
        min_size = as.integer(lcfg$min_size %||% 10L),
        max_size = as.integer(lcfg$max_size %||% 500L),
        seed     = as.integer(config$params$seed %||% 1L)
    )
}


#' Resolve GENE_N feature ids without changing the vector's length
#'
#' The aligned counterpart of \code{resolve_gene_n_ids()}, which drops the ids
#' it cannot resolve. That is right for a feature list, and wrong here: each id
#' carries a loading, and dropping entries would shift every later id onto the
#' wrong value. Same mapping, same accepted spellings (bare `GENE_12` and the
#' MOFA2 view-suffixed `GENE_12_proteomics`); an unresolvable id becomes NA in
#' place.
#'
#' @param feature_ids Character vector of feature ids.
#' @param harmonization_res Harmonization result carrying
#'   \code{gene_protein_mapping}.
#' @param omics_type "transcriptomics" or "proteomics"; anything else is
#'   returned untouched.
#' @return Character vector the same length as \code{feature_ids}.
#' @keywords internal
.resolve_feature_ids_aligned <- function(feature_ids, harmonization_res, omics_type) {
    ids <- as.character(feature_ids)
    gpm <- harmonization_res$gene_protein_mapping
    id_col <- c(transcriptomics = "gene_id", proteomics = "protein_id")[omics_type]
    if (is.null(gpm) || is.na(id_col) || !id_col %in% names(gpm)) return(ids)

    is_gene_n <- grepl("^GENE_\\d+(_.+)?$", ids)
    if (!any(is_gene_n)) return(ids)

    idx <- as.integer(sub("^GENE_(\\d+)(_.+)?$", "\\1", ids[is_gene_n]))
    idx[idx < 1 | idx > nrow(gpm)] <- NA_integer_
    ids[is_gene_n] <- as.character(gpm[[id_col]])[idx]
    ids
}


#' Collapse signed loadings onto one value per key
#'
#' Several features can land on one gene or one compound -- protein isoforms,
#' paralogues, metabolites annotated to the same KEGG compound -- and GSEA would
#' otherwise score that key more than once. The rule is the one
#' \code{rank_compounds_for_gsea()} uses: the largest absolute value wins, a
#' `+x` / `-x` tie goes to the positive value, and the result is ordered
#' explicitly, so it does not depend on the order features arrived in.
#'
#' @param values Numeric vector of signed loadings or weights.
#' @param keys Character vector of the same length: the gene or compound id
#'   each value belongs to. NA or empty keys are dropped.
#' @return Named numeric vector, unique names, sorted decreasing;
#'   \code{numeric(0)} when nothing is usable.
#' @examples
#' collapse_loading_ranks(c(0.4, -0.9, 0.2), c("g1", "g1", "g2"))
#' # g2 = 0.2, g1 = -0.9
collapse_loading_ranks <- function(values, keys) {
    values <- as.numeric(values)
    keys <- as.character(keys)
    keep <- is.finite(values) & !is.na(keys) & nzchar(keys)
    if (!any(keep)) return(numeric(0))

    ranks <- stats::setNames(values[keep], keys[keep])
    ranks <- ranks[order(-abs(ranks), names(ranks), -ranks)]
    ranks <- ranks[!duplicated(names(ranks))]
    ranks[order(-ranks, names(ranks))]
}


#' One ranking per DIABLO block and component
#'
#' `block.plsda` loadings are dense -- every feature of a block has one -- so
#' the ranking covers the whole block, not a selection. The outcome block `Y`
#' is not an omics layer and is skipped.
#'
#' @param diablo_results Output of \code{run_diablo_integration()}: uses
#'   `loadings` and `sample_scores` (one features-by-components and one
#'   samples-by-components matrix per block).
#' @return Named list of entries, each a list with `label`, `integration`,
#'   `omics`, `axis`, `values` (named by feature) and `scores` (named by
#'   sample, or NULL).
diablo_loading_entries <- function(diablo_results) {
    loadings <- diablo_results$loadings
    scores <- diablo_results$sample_scores
    entries <- list()
    if (is.null(loadings)) return(entries)

    for (om in setdiff(names(loadings), "Y")) {
        m <- loadings[[om]]
        if (is.null(m) || NROW(m) == 0) next
        m <- as.matrix(m)
        sc <- scores[[om]]
        for (comp in colnames(m)) {
            s <- NULL
            if (!is.null(sc) && comp %in% colnames(sc)) {
                s <- stats::setNames(as.numeric(sc[, comp]), rownames(sc))
            }
            label <- paste0("DIABLO_", om, "_", comp)
            entries[[label]] <- list(
                label = label, integration = "DIABLO", omics = om, axis = comp,
                values = stats::setNames(as.numeric(m[, comp]), rownames(m)),
                scores = s)
        }
    }
    entries
}


#' One ranking per MOFA2 view and factor
#'
#' The first `n_factors` factors, as the ORA path takes.
#'
#' @param mofa_results Output of \code{run_mofa_integration()}: uses `weights`
#'   (one features-by-factors matrix per view) and `factors` (a data frame with
#'   `sample_id` and one column per factor, or a samples-by-factors matrix).
#' @param n_factors Number of leading factors to rank.
#' @return Named list of entries in the shape \code{diablo_loading_entries()}
#'   returns.
mofa_weight_entries <- function(mofa_results, n_factors = 3L) {
    weights <- mofa_results$weights
    factors <- mofa_results$factors
    entries <- list()
    if (is.null(weights)) return(entries)

    factor_scores <- function(fname) {
        if (is.null(factors)) return(NULL)
        if (is.data.frame(factors) && "sample_id" %in% names(factors) &&
            fname %in% names(factors)) {
            return(stats::setNames(as.numeric(factors[[fname]]), factors$sample_id))
        }
        if (is.matrix(factors) && fname %in% colnames(factors)) {
            return(stats::setNames(as.numeric(factors[, fname]), rownames(factors)))
        }
        NULL
    }

    for (view in names(weights)) {
        w <- weights[[view]]
        if (is.null(w) || NROW(w) == 0) next
        w <- as.matrix(w)
        for (k in seq_len(min(ncol(w), n_factors))) {
            fname <- colnames(w)[k]
            label <- paste0("MOFA_", view, "_", fname)
            entries[[label]] <- list(
                label = label, integration = "MOFA2", omics = view, axis = fname,
                values = stats::setNames(as.numeric(w[, k]), rownames(w)),
                scores = factor_scores(fname))
        }
    }
    entries
}


#' Sample-to-condition lookup from the harmonized MAE
#'
#' @param harmonization_res Harmonization result carrying `mae`.
#' @param config Full config object, for the condition column.
#' @return Named character vector, condition per sample; NULL when the MAE or
#'   the column is missing.
loadings_sample_conditions <- function(harmonization_res, config) {
    mae <- harmonization_res$mae
    if (is.null(mae)) return(NULL)
    cond_col <- config$modes$multiomics$condition_column %||%
                config$design$condition_column %||% "condition"
    cd <- as.data.frame(SummarizedExperiment::colData(mae))
    if (!cond_col %in% names(cd)) return(NULL)
    stats::setNames(as.character(cd[[cond_col]]), rownames(cd))
}


#' Which condition scores higher on a component or factor
#'
#' The sign of a DIABLO component or MOFA2 factor is arbitrary: refitting can
#' flip it, and with it every NES. A positive NES means the pathway's features
#' load in the same direction as the samples that score high on the axis, so
#' the NES is only readable once it is said which condition those samples
#' belong to. Group means of the sample scores decide it; ties break on the
#' condition name so the answer does not depend on sample order.
#'
#' @param scores Named numeric vector of sample scores on the axis.
#' @param conditions Named character vector, condition per sample.
#' @param axis Axis label for the note, e.g. "Factor1".
#' @return List with `higher`, `lower` (condition names, NA when fewer than two
#'   conditions have scores) and `note`, a sentence for the plot subtitle.
#' @examples
#' axis_orientation(c(s1 = 2, s2 = 1, s3 = -1, s4 = -2),
#'                  c(s1 = "B", s2 = "B", s3 = "A", s4 = "A"), "comp1")$higher  # "B"
axis_orientation <- function(scores, conditions, axis) {
    none <- list(higher = NA_character_, lower = NA_character_,
                 note = paste0("Condition means on ", axis,
                               " not available; the NES sign is relative to ",
                               "the axis orientation."))
    if (is.null(scores) || is.null(conditions)) return(none)

    common <- intersect(names(scores), names(conditions))
    s <- scores[common]
    g <- conditions[common]
    ok <- is.finite(s) & !is.na(g) & nzchar(g)
    if (length(unique(g[ok])) < 2) return(none)

    m <- tapply(s[ok], g[ok], mean)
    m <- m[order(-m, names(m))]
    higher <- names(m)[1]
    lower <- names(m)[length(m)]
    list(higher = higher, lower = lower,
         note = sprintf("%s scores higher than %s on %s; positive NES points towards %s",
                        higher, lower, axis, higher))
}


#' GSEA of one gene view's loadings
#'
#' The gene-set source follows \code{enrich_feature_list()}: KEGG through
#' clusterProfiler, keyed on NCBI gene ids, when the organism has a KEGG code
#' and an OrgDb; otherwise the view's configured GMTs, one collection per file
#' as \code{load_gene_sets()} builds them, scored by the core fgsea in
#' \code{run_pathway_analysis()}. Every call is seeded. A configured GMT that
#' does not exist is skipped with a warning; the others are still scored.
#'
#' @param values Named numeric vector, signed loading per feature id.
#' @param omics_type "transcriptomics" or "proteomics".
#' @param harmonization_res Harmonization result.
#' @param config Full config object.
#' @param kegg_org KEGG organism code, or NULL.
#' @param org_db OrgDb object, or NULL.
#' @param min_size,max_size Gene-set size bounds.
#' @param seed Seed for the stochastic scoring.
#' @param exclude_classes KEGG BRITE classes to drop after scoring, or NULL.
#' @return Data frame with `pathway`, `ID`, `NES`, `ES`, `pvalue`, `padj`,
#'   `setSize`, `leadingEdge`, `database`, `method`; NULL when nothing scored,
#'   including when none of the configured GMTs exists.
gene_loadings_gsea <- function(values, omics_type, harmonization_res, config,
                               kegg_org, org_db, min_size, max_size, seed,
                               exclude_classes = NULL) {
    keys <- .resolve_feature_ids_aligned(names(values), harmonization_res, omics_type)

    if (!is.null(kegg_org) && !is.null(org_db)) {
        ids <- unique(keys[!is.na(keys)])
        id_map <- tryCatch(
            map_feature_ids_to_entrez(
                de_tables = list(loadings = data.frame(feature_id = ids,
                                                       stringsAsFactors = FALSE)),
                omics_type = omics_type,
                harmonization_res = harmonization_res,
                org_db = org_db),
            error = function(e) {
                message("    ID mapping failed: ", conditionMessage(e))
                NULL
            })
        if (is.null(id_map) || nrow(id_map) == 0) return(NULL)

        ranks <- collapse_loading_ranks(values, id_map$ENTREZID[match(keys, id_map$feature_id)])
        if (length(ranks) < min_size) {
            message("    Too few features mapped to NCBI gene ids (", length(ranks), ")")
            return(NULL)
        }
        message("    GSEA (KEGG, clusterProfiler): ", length(ranks), " ranked genes")

        res <- tryCatch(
            withr::with_seed(seed, clusterProfiler::gseKEGG(
                geneList = ranks, organism = kegg_org, keyType = "ncbi-geneid",
                minGSSize = min_size, maxGSSize = max_size,
                pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE,
                seed = FALSE)),
            error = function(e) {
                message("    clusterProfiler::gseKEGG failed: ", conditionMessage(e))
                NULL
            })
        if (is.null(res)) return(NULL)
        df <- as.data.frame(res)
        if (nrow(df) == 0) return(NULL)

        out <- data.frame(
            pathway = df$Description, ID = df$ID, NES = df$NES,
            ES = df$enrichmentScore, pvalue = df$pvalue, padj = df$p.adjust,
            setSize = df$setSize, leadingEdge = gsub("/", ",", df$core_enrichment),
            database = "KEGG", method = "gsea", stringsAsFactors = FALSE)
    } else {
        cfg_key <- c(transcriptomics = "rna", proteomics = "proteomics")[omics_type]
        if (is.na(cfg_key) || is.null(config)) return(NULL)
        gmt_path <- unlist(config$modes[[cfg_key]]$pathway$gmt_file, use.names = FALSE)
        if (length(gmt_path) == 0 || !any(nzchar(gmt_path))) return(NULL)

        # A missing file costs only its own collection. None left is NULL here:
        # load_gene_sets() given no GMT would fall back to generating gene sets
        # for a non-model organism, which is not what this view was configured
        # to test against.
        gmt_abs <- resolve_input_path(config, gmt_path[nzchar(gmt_path)])
        missing_gmt <- !file.exists(gmt_abs)
        if (any(missing_gmt)) {
            warning("Loadings GSEA (GMT): ", omics_type, " gmt_file not found, ",
                    "skipped: ", paste(gmt_abs[missing_gmt], collapse = ", "),
                    call. = FALSE)
        }
        if (all(missing_gmt)) return(NULL)

        gene_sets <- load_gene_sets(config$global$organism,
                                    pathway_database = character(0),
                                    gmt_file = gmt_abs[!missing_gmt])
        if (length(gene_sets) == 0) return(NULL)

        ranks <- collapse_loading_ranks(values, keys)
        if (length(ranks) < min_size) {
            message("    Too few ranked features for GSEA (", length(ranks), ")")
            return(NULL)
        }
        message("    GSEA (GMT, fgsea): ", length(ranks), " ranked features, ",
                length(gene_sets), " collection(s)")

        tabs <- run_pathway_analysis(
            de_tables = list(loadings = data.frame(FeatureID = names(ranks),
                                                   stat = unname(ranks),
                                                   stringsAsFactors = FALSE)),
            gene_sets = gene_sets, method = "fgsea",
            min_size = min_size, max_size = max_size, seed = seed)$loadings
        if (length(tabs) == 0) return(NULL)

        out <- .rbind_fill(lapply(tabs, function(df) {
            label <- if ("pathway_name" %in% names(df)) df$pathway_name else df$pathway
            label[is.na(label) | !nzchar(label)] <- df$pathway[is.na(label) | !nzchar(label)]
            data.frame(
                pathway = label, ID = df$pathway, NES = df$NES, ES = df$ES,
                pvalue = df$pval, padj = df$padj, setSize = df$size,
                leadingEdge = df$leadingEdge, database = df$database,
                method = "fgsea", stringsAsFactors = FALSE)
        }))
        if (is.null(out) || nrow(out) == 0) return(NULL)
    }

    if (length(unlist(exclude_classes)) > 0) {
        out <- out[keep_kegg_pathways(out$ID, exclude = exclude_classes,
                                      kegg_org = kegg_org,
                                      label = "loadings GSEA pathways"), , drop = FALSE]
    }
    if (nrow(out) == 0) return(NULL)
    out
}


#' GSEA of one metabolomics view's loadings
#'
#' Metabolites reach KEGG compound ids through \code{map_metabolite_ids_to_kegg()},
#' the mapper the ORA path uses, and are scored by \code{run_compound_gsea()}
#' with the loading as the ranking statistic.
#'
#' @param values Named numeric vector, signed loading per feature id.
#' @param harmonization_res Harmonization result.
#' @param cache_dir Directory for the compound-pathway cache.
#' @param min_size,max_size Pathway size bounds.
#' @param seed Seed for the stochastic scoring.
#' @param exclude_classes KEGG BRITE classes to drop after scoring, or NULL.
#' @param cpd_pathways Compound-pathway associations; fetched through
#'   \code{get_kegg_compound_pathways()} when NULL.
#' @return Data frame of scored pathways, or NULL.
metabolite_loadings_gsea <- function(values, harmonization_res, cache_dir,
                                     min_size, max_size, seed,
                                     exclude_classes = NULL, cpd_pathways = NULL) {
    ids <- names(values)
    id_map <- tryCatch(
        map_metabolite_ids_to_kegg(
            de_tables = list(loadings = data.frame(feature_id = ids,
                                                   stringsAsFactors = FALSE)),
            harmonization_res = harmonization_res),
        error = function(e) {
            message("    Metabolite ID mapping failed: ", conditionMessage(e))
            NULL
        })
    if (is.null(id_map) || nrow(id_map) == 0) {
        message("    Could not map metabolomics features to KEGG compound IDs")
        return(NULL)
    }

    de_mapped <- merge(data.frame(feature_id = ids, statistic = unname(values),
                                  stringsAsFactors = FALSE),
                       id_map, by = "feature_id")
    de_mapped$KEGG_ID <- de_mapped$KEGG_CPD

    if (is.null(cpd_pathways)) cpd_pathways <- get_kegg_compound_pathways(cache_dir)
    run_compound_gsea(de_mapped, cache_dir = cache_dir, min_gs = min_size,
                      max_gs = max_size, seed = seed,
                      exclude_classes = exclude_classes, cpd_pathways = cpd_pathways)
}


#' NES bar chart for one component or factor
#'
#' @param gsea_df One entry's GSEA table: `pathway`, `ID`, `NES`, `pvalue`,
#'   `padj`.
#' @param title Plot title, usually the entry label.
#' @param subtitle Orientation note from \code{axis_orientation()}, or NULL.
#' @param top_n Number of pathways, taken by p-value.
#' @param nes_guide Position of the dashed guide lines, at -nes_guide and
#'   +nes_guide.
#' @return A ggplot, or NULL when there is nothing to draw.
plot_loadings_gsea_nes <- function(gsea_df, title, subtitle = NULL, top_n = 15,
                                   nes_guide = 2) {
    if (is.null(gsea_df) || nrow(gsea_df) == 0) return(NULL)
    df <- gsea_df[is.finite(gsea_df$NES) & is.finite(gsea_df$pvalue), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)

    df <- df[order(df$pvalue, df$ID), , drop = FALSE]
    df <- df[seq_len(min(top_n, nrow(df))), , drop = FALSE]

    lab <- ifelse(nchar(df$pathway) > 45, paste0(substr(df$pathway, 1, 42), "..."),
                  df$pathway)
    lab <- make.unique(lab, sep = " ")
    df$label <- factor(lab, levels = lab[order(df$NES, lab)])
    df$fdr_text <- paste("FDR", formatC(df$padj, format = "g", digits = 2))
    df$text_hjust <- ifelse(df$NES >= 0, -0.1, 1.1)
    x_lim <- max(abs(df$NES), nes_guide) * 1.45

    ggplot2::ggplot(df, ggplot2::aes(x = NES, y = label)) +
        ggplot2::geom_col(fill = "steelblue", width = 0.7) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey50") +
        ggplot2::geom_vline(xintercept = c(-nes_guide, nes_guide),
                            colour = "red", linetype = "dashed") +
        ggplot2::geom_text(ggplot2::aes(label = fdr_text, hjust = text_hjust),
                           size = 2.6, colour = "grey30") +
        ggplot2::scale_x_continuous(limits = c(-x_lim, x_lim)) +
        ggplot2::labs(
            title = paste("Loadings GSEA:", title), subtitle = subtitle,
            x = "NES", y = NULL,
            caption = sprintf(paste0("Top %d sets by p-value. Dashed red lines at ",
                                     "NES = -%s and +%s are a visual guide, not a ",
                                     "significance test."),
                              nrow(df), nes_guide, nes_guide)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = 7),
            plot.subtitle = ggplot2::element_text(size = 8),
            plot.caption = ggplot2::element_text(size = 7, colour = "grey30"))
}


#' Remove a previous run's loadings GSEA files from one directory
#'
#' Scoped to the files this path writes, so the ORA outputs beside them are
#' left alone.
#'
#' @param dir Output directory.
#' @param combined_name File name of the combined table.
#' @return Invisibly, the paths removed.
#' @keywords internal
.clear_loadings_gsea_outputs <- function(dir, combined_name) {
    old <- c(list.files(dir, pattern = "_gsea\\.csv$|_gsea_nes\\.png$", full.names = TRUE),
             file.path(dir, combined_name))
    old <- old[file.exists(old)]
    unlink(old)
    invisible(old)
}


#' Remove a previous run's loadings ORA files from one directory
#'
#' The ORA counterpart of \code{.clear_loadings_gsea_outputs()}: the
#' `<label>_enrichment.csv` / `.png` files and the combined table, and nothing
#' else, so GSEA outputs and the KEGG caches beside them are left alone.
#'
#' @param dir Output directory.
#' @param combined_name File name of the combined ORA table.
#' @return Invisibly, the paths removed.
#' @keywords internal
.clear_loadings_ora_outputs <- function(dir, combined_name) {
    old <- c(list.files(dir, pattern = "_enrichment\\.(csv|png)$", full.names = TRUE),
             file.path(dir, combined_name))
    old <- old[file.exists(old)]
    unlink(old)
    invisible(old)
}


#' Remove the ORA files of both integrations before an ORA run
#'
#' @param out_dir The loadings-enrichment output directory.
#' @return Invisibly, the paths removed.
#' @keywords internal
.clear_loadings_ora_all <- function(out_dir) {
    invisible(c(
        .clear_loadings_ora_outputs(file.path(out_dir, "diablo_loadings"),
                                    "diablo_loadings_enrichment_all.csv"),
        .clear_loadings_ora_outputs(file.path(out_dir, "mofa_loadings"),
                                    "mofa_weights_enrichment_all.csv")))
}


#' Remove everything the loadings enrichment step writes
#'
#' For a run that skips or fails the step: the method record and both tests'
#' outputs for both integrations go, so the report finds nothing rather than
#' a previous run presented as this one. Caches and anything else in the
#' directories are left alone.
#'
#' @param out_dir The loadings-enrichment output directory.
#' @return Invisibly, the paths removed.
clear_loadings_enrichment_outputs <- function(out_dir) {
    record <- file.path(out_dir, "loadings_enrichment_method.txt")
    removed <- c(
        .clear_loadings_gsea_outputs(file.path(out_dir, "diablo_loadings"),
                                     "diablo_loadings_gsea_all.csv"),
        .clear_loadings_gsea_outputs(file.path(out_dir, "mofa_loadings"),
                                     "mofa_weights_gsea_all.csv"),
        .clear_loadings_ora_all(out_dir),
        record[file.exists(record)])
    unlink(record)
    invisible(removed)
}


#' Score, write and plot a set of loadings rankings
#'
#' @param entries Entries from \code{diablo_loading_entries()} or
#'   \code{mofa_weight_entries()}.
#' @param ctx List of shared inputs: `harmonization_res`, `config`, `kegg_org`,
#'   `org_db`, `min_size`, `max_size`, `seed`, `exclude_classes`,
#'   `cpd_pathways`, `conditions`, `cache_dir`.
#' @param out_dir Directory for `<label>_gsea.csv`, `<label>_gsea_nes.png` and
#'   the combined table.
#' @param combined_name File name of the combined table.
#' @return Combined data frame of every entry's results, or NULL.
run_loadings_gsea_entries <- function(entries, ctx, out_dir, combined_name) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    .clear_loadings_gsea_outputs(out_dir, combined_name)

    all_results <- list()
    for (e in entries) {
        message("  ", e$label, ": ", length(e$values), " ranked features")
        res <- if (identical(e$omics, "metabolomics") && is.null(ctx$cpd_pathways)) {
            # The one shared lookup already failed; NULL here would make
            # metabolite_loadings_gsea() fetch again, once per entry.
            message("    KEGG compound pathways unavailable; skipping ", e$label)
            NULL
        } else if (identical(e$omics, "metabolomics")) {
            metabolite_loadings_gsea(
                e$values, ctx$harmonization_res, cache_dir = ctx$cache_dir,
                min_size = ctx$min_size, max_size = ctx$max_size, seed = ctx$seed,
                exclude_classes = ctx$exclude_classes, cpd_pathways = ctx$cpd_pathways)
        } else if (e$omics %in% c("transcriptomics", "proteomics")) {
            gene_loadings_gsea(
                e$values, e$omics, ctx$harmonization_res, ctx$config,
                kegg_org = ctx$kegg_org, org_db = ctx$org_db,
                min_size = ctx$min_size, max_size = ctx$max_size, seed = ctx$seed,
                exclude_classes = ctx$exclude_classes)
        } else {
            NULL
        }
        if (is.null(res) || nrow(res) == 0) {
            message("    No pathways scored for ", e$label)
            next
        }

        orient <- axis_orientation(e$scores, ctx$conditions, e$axis)
        res$integration <- e$integration
        res$omics <- e$omics
        res$axis <- e$axis
        res$higher_scoring_group <- orient$higher
        res$lower_scoring_group <- orient$lower
        res <- res[order(res$pvalue, res$ID), , drop = FALSE]
        rownames(res) <- NULL

        write.csv(res, file.path(out_dir, paste0(e$label, "_gsea.csv")), row.names = FALSE)
        p <- plot_loadings_gsea_nes(res, e$label, subtitle = orient$note)
        if (!is.null(p)) {
            ggplot2::ggsave(file.path(out_dir, paste0(e$label, "_gsea_nes.png")),
                            plot = p, width = 9, height = 6, dpi = 150)
        }
        all_results[[e$label]] <- res
    }

    combined <- .rbind_fill(all_results)
    if (is.null(combined) || nrow(combined) == 0) return(NULL)
    rownames(combined) <- NULL
    write.csv(combined, file.path(out_dir, combined_name), row.names = FALSE)
    combined
}


#' GSEA on DIABLO loadings and MOFA2 weights
#'
#' The GSEA branch of \code{run_loadings_enrichment()}; same inputs, same two
#' output directories, same return shape.
#'
#' @param integration_res Output from mod_multiomics_integration().
#' @param harmonization_res Output from mod_multiomics_harmonization().
#' @param config Full config object.
#' @param out_dir Loadings-enrichment output directory.
#' @return List with `diablo` and `mofa` combined tables (either may be NULL).
run_loadings_gsea <- function(integration_res, harmonization_res, config, out_dir) {
    # Cleared up front, for both integrations: the method record will say
    # "gsea", so a previous run's files left behind by an integration that is
    # now absent, or by a failure further down, would be reported as current.
    diablo_dir <- file.path(out_dir, "diablo_loadings")
    mofa_dir <- file.path(out_dir, "mofa_loadings")
    .clear_loadings_gsea_outputs(diablo_dir, "diablo_loadings_gsea_all.csv")
    .clear_loadings_gsea_outputs(mofa_dir, "mofa_weights_gsea_all.csv")

    organism <- config$global$organism
    settings <- .loadings_gsea_settings(config)

    diablo_entries <- if (!is.null(integration_res$diablo_results)) {
        diablo_loading_entries(integration_res$diablo_results)
    } else list()
    mofa_entries <- if (!is.null(integration_res$mofa_results)) {
        mofa_weight_entries(integration_res$mofa_results)
    } else list()

    # One compound-pathway read for both methods, and none at all when no
    # metabolomics view is ranked.
    has_metab <- any(vapply(c(diablo_entries, mofa_entries),
                            function(e) identical(e$omics, "metabolomics"), logical(1)))
    cpd_pathways <- if (has_metab) {
        tryCatch(get_kegg_compound_pathways(out_dir), error = function(e) {
            message("  KEGG compound pathways unavailable: ", conditionMessage(e))
            NULL
        })
    } else NULL

    # run_compound_gsea() reads the class table from this cache only (it never
    # fetches), so without it a class exclusion would fail open for metabolites
    # while the gene views apply it.
    exclude_classes <- .excluded_pathway_classes(config)
    if (has_metab && !is.null(cpd_pathways) && length(unlist(exclude_classes)) > 0) {
        tryCatch(kegg_pathway_categories(cache_dir = out_dir), error = function(e) {
            message("  KEGG pathway classification unavailable: ", conditionMessage(e))
            NULL
        })
    }

    ctx <- list(
        harmonization_res = harmonization_res, config = config,
        kegg_org = get_kegg_organism(organism), org_db = get_organism_db(organism),
        min_size = settings$min_size, max_size = settings$max_size,
        seed = settings$seed, exclude_classes = exclude_classes,
        cpd_pathways = cpd_pathways,
        conditions = loadings_sample_conditions(harmonization_res, config),
        cache_dir = out_dir)

    results <- list()
    if (length(diablo_entries) > 0) {
        message("Running GSEA on DIABLO loadings...")
        results$diablo <- tryCatch(
            run_loadings_gsea_entries(diablo_entries, ctx, diablo_dir,
                                      "diablo_loadings_gsea_all.csv"),
            error = function(e) {
                message("  DIABLO loadings GSEA failed: ", conditionMessage(e))
                NULL
            })
    }
    if (length(mofa_entries) > 0) {
        message("Running GSEA on MOFA2 weights...")
        results$mofa <- tryCatch(
            run_loadings_gsea_entries(mofa_entries, ctx, mofa_dir,
                                      "mofa_weights_gsea_all.csv"),
            error = function(e) {
                message("  MOFA2 weights GSEA failed: ", conditionMessage(e))
                NULL
            })
    }

    message("Loadings GSEA complete")
    results
}
