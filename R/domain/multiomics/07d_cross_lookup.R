#' Cross-omics pathway lookup: one layer's top pathways, read in another
#'
#' The meta-analysis heatmaps rank pathways on how many layers support them,
#' so a layer whose own results are weak -- no feature past its DE threshold,
#' nothing below FDR -- rarely sets the rows. This view starts from one layer
#' instead: its top pathways by its own evidence, whether or not they are
#' significant, and what the other layer shows for the same pathways. It runs
#' in both directions.
#'
#' Descriptive only. Every p-value shown is the one the layer's own enrichment
#' produced; nothing is re-tested or re-adjusted.
#'
#' Identity and method choice are the meta-analysis's own:
#' \code{.layer_contribution()} decides which rows of a layer are read -- its
#' rank-based (GSEA) rows where it has them, otherwise its ORA rows -- and keys
#' them with \code{pathway_join_key()}. The target is always read in the same
#' contrast as the source row it sits beside, so two NES are never compared
#' across different comparisons.


#' Best row per pathway within the chosen rows of one layer
#'
#' The collapse is the one \code{merge_pathway_pvalues()} makes: the smallest
#' raw p-value per join key, across contrasts and collections -- or, with
#' \code{by_contrast}, per join key within each contrast, which is how a target
#' is read so it can be matched to the source row's own contrast. Ties go to the
#' larger |NES|, then to the join key and the contrast, so the order never
#' depends on the order rows arrived in.
#'
#' @param df Enrichment data frame for one layer.
#' @param contrib Its \code{.layer_contribution()} result, which supplies the
#'   join keys, the raw p-values and the rows the layer contributes.
#' @param keep Logical vector, one per row of \code{df}: the rows to read.
#'   Defaults to the contributing rows; the target's ORA values beside its
#'   chosen method pass the ORA rows instead.
#' @param by_contrast Collapse per (pathway, contrast) rather than per pathway.
#' @return Data frame with one row per join key (or per join key and contrast)
#'   -- \code{norm_id}, \code{label}, \code{method}, \code{NES}, \code{p},
#'   \code{padj}, \code{n_measured}, \code{contrast}, \code{contrast_key} --
#'   sorted best first; NULL when no row is usable. \code{n_measured} is the
#'   number of pathway members present in the layer's ranked list (fgsea
#'   `size`, GSEA `setSize`) and is NA for ORA rows, whose set-size columns
#'   count something else. \code{contrast_key} is
#'   \code{normalize_contrast_key()} of the contrast, NA where there is none.
#' @keywords internal
.lookup_layer_stats <- function(df, contrib, keep = contrib$keep,
                                by_contrast = FALSE) {
    if (!is.data.frame(df) || nrow(df) == 0 || is.null(contrib) ||
        !is.null(contrib$problem)) return(NULL)
    ok <- keep & !is.na(contrib$keys) & !is.na(contrib$pvals)
    if (!any(ok)) return(NULL)

    # First non-missing value per row across the named columns, numeric.
    num_col <- function(cols) {
        vals <- rep(NA_real_, nrow(df))
        for (col in intersect(cols, names(df))) {
            v <- df[[col]]
            if (!is.numeric(v)) v <- suppressWarnings(as.numeric(as.character(v)))
            fill <- is.na(vals) & !is.na(v)
            vals[fill] <- v[fill]
        }
        vals
    }

    method <- if ("method" %in% names(df)) {
        tolower(trimws(as.character(df$method)))
    } else {
        rep("unspecified", nrow(df))
    }
    n_measured <- num_col(c("size", "setSize"))
    n_measured[!(method %in% .RANK_BASED_METHODS)] <- NA_real_
    contrast <- if ("contrast" %in% names(df)) as.character(df$contrast)
                else rep(NA_character_, nrow(df))

    tab <- data.frame(
        norm_id      = contrib$keys,
        label        = pathway_display_label(df),
        method       = method,
        NES          = num_col("NES"),
        p            = contrib$pvals,
        padj         = .ora_adjusted_p_values(df),
        n_measured   = n_measured,
        contrast     = contrast,
        contrast_key = normalize_contrast_key(contrast),
        stringsAsFactors = FALSE
    )[ok, , drop = FALSE]

    tab <- tab[order(tab$p, -abs(tab$NES), tab$norm_id, tab$contrast_key,
                     na.last = TRUE), , drop = FALSE]
    group <- if (by_contrast) .lookup_pair_key(tab) else tab$norm_id
    tab <- tab[!duplicated(group), , drop = FALSE]
    rownames(tab) <- NULL
    tab
}


#' Pathway-and-contrast key used to match a source row to its target row
#'
#' A row with no contrast only ever matches another row with no contrast: a
#' source row from one comparison is never paired with a target row from
#' another, or from a table that does not say which comparison it is.
#'
#' @param tab Rows from \code{.lookup_layer_stats()}.
#' @return Character vector, one key per row.
#' @keywords internal
.lookup_pair_key <- function(tab) {
    paste(tab$norm_id, ifelse(is.na(tab$contrast_key), "<none>", tab$contrast_key),
          sep = "\r")
}


#' Look up one layer's top pathways in another layer
#'
#' The source layer is ranked on its own raw p-value, from the rows
#' \code{.layer_contribution()} chooses for it -- rank-based where it has them,
#' which is what "ranked by GSEA" means here. Two rankings are returned:
#'
#' \itemize{
#'   \item \code{all}: the source's top \code{top_n} pathways, whether or not the
#'     target has a result for them. Many can be absent from the target -- a
#'     gene-based pathway may contain no measured compound -- and saying so is
#'     part of the answer.
#'   \item \code{tested_in_target}: the source's top \code{top_n} among the
#'     pathways the target has a result for in the same contrast.
#' }
#'
#' Each source pathway keeps its best row across contrasts, and the target is
#' read in that row's contrast (matched on \code{normalize_contrast_key()}): a
#' target with no result in that contrast is reported as such, never filled in
#' from another contrast. \code{source_contrast} and \code{target_contrast}
#' say which comparison each value comes from. For the target the chosen
#' method's values are reported and, where the target also carries ORA rows,
#' its ORA p and adjusted p in the same contrast beside them.
#'
#' @param pathway_tables Named list of per-omics enrichment data frames -- the
#'   tables the meta-analysis merged, rank-based rows included.
#' @param from Name of the layer to rank.
#' @param to Name of the layer to read the ranked pathways in.
#' @param top_n Number of pathways per ranking.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return Data frame, one row per (ranking, pathway), or NULL when the source
#'   layer has no usable rows. Columns: ranking, rank, norm_id, pathway,
#'   source_layer, source_method, source_contrast, source_NES, source_p,
#'   source_padj, source_n_measured, target_layer, target_status,
#'   target_method, target_contrast, target_NES, target_p, target_padj,
#'   target_n_measured, target_ora_p, target_ora_padj.
#' @examples
#' tabs <- list(
#'   proteomics   = data.frame(ID = c("map00010", "map00020"), pval = c(0.01, 0.2),
#'                             NES = c(1.9, -1.1), padj = c(0.1, 0.5),
#'                             method = "fgsea"),
#'   metabolomics = data.frame(ID = "map00010", pvalue = 0.3, NES = -0.8,
#'                             padj = 0.6, method = "fgsea"))
#' build_cross_omics_lookup(tabs, "proteomics", "metabolomics", top_n = 2)
build_cross_omics_lookup <- function(pathway_tables, from, to, top_n = 15,
                                     kegg_org = NULL) {
    src_df <- pathway_tables[[from]]
    if (!is.data.frame(src_df) || nrow(src_df) == 0) return(NULL)
    src <- .lookup_layer_stats(src_df, .layer_contribution(src_df, kegg_org))
    if (is.null(src)) return(NULL)

    tgt <- NULL
    tgt_ora <- NULL
    tgt_df <- pathway_tables[[to]]
    if (is.data.frame(tgt_df) && nrow(tgt_df) > 0) {
        tgt_cb <- .layer_contribution(tgt_df, kegg_org)
        tgt <- .lookup_layer_stats(tgt_df, tgt_cb, by_contrast = TRUE)
        if (!is.null(tgt) && "method" %in% names(tgt_df)) {
            is_ora <- !is.na(tgt_df$method) &
                tolower(trimws(as.character(tgt_df$method))) == "ora"
            tgt_ora <- .lookup_layer_stats(tgt_df, tgt_cb, keep = is_ora,
                                           by_contrast = TRUE)
        }
    }

    top_n <- suppressWarnings(as.integer(top_n))
    if (length(top_n) != 1 || is.na(top_n) || top_n < 1) top_n <- 15L

    tested <- if (is.null(tgt)) rep(FALSE, nrow(src))
              else .lookup_pair_key(src) %in% .lookup_pair_key(tgt)
    blocks <- list(
        all              = utils::head(src, top_n),
        tested_in_target = utils::head(src[tested, , drop = FALSE], top_n)
    )

    out <- lapply(names(blocks), function(b) {
        .assemble_lookup_block(blocks[[b]], b, from, to, tgt, tgt_ora)
    })
    out <- Filter(Negate(is.null), out)
    if (length(out) == 0) return(NULL)
    out <- do.call(rbind, out)
    rownames(out) <- NULL
    out
}


#' One ranking of a lookup, with the target's values joined on
#'
#' Target values are joined on pathway and contrast together
#' (\code{.lookup_pair_key()}), so each target value is the one for the source
#' row's own contrast.
#'
#' @param blk Source rows for this ranking, best first, from
#'   \code{.lookup_layer_stats()}.
#' @param ranking Name of the ranking ("all" or "tested_in_target").
#' @param from,to Layer names.
#' @param tgt Target stats per pathway and contrast for its chosen method, or
#'   NULL.
#' @param tgt_ora Target stats per pathway and contrast for its ORA rows, or
#'   NULL.
#' @return Data frame for this ranking, or NULL when \code{blk} is empty.
#' @keywords internal
.assemble_lookup_block <- function(blk, ranking, from, to, tgt, tgt_ora) {
    if (is.null(blk) || nrow(blk) == 0) return(NULL)
    n <- nrow(blk)

    key <- .lookup_pair_key(blk)
    t_idx <- if (is.null(tgt)) rep(NA_integer_, n) else match(key, .lookup_pair_key(tgt))
    o_idx <- if (is.null(tgt_ora)) rep(NA_integer_, n) else match(key, .lookup_pair_key(tgt_ora))
    pick <- function(tab, idx, col, na) if (is.null(tab)) rep(na, n) else tab[[col]][idx]

    tested <- !is.na(t_idx)

    label <- blk$label
    tgt_label <- pick(tgt, t_idx, "label", NA_character_)
    label[is.na(label)] <- tgt_label[is.na(label)]
    label[is.na(label)] <- blk$norm_id[is.na(label)]

    data.frame(
        ranking           = ranking,
        rank              = seq_len(n),
        norm_id           = blk$norm_id,
        pathway           = label,
        source_layer      = from,
        source_method     = blk$method,
        source_contrast   = blk$contrast,
        source_NES        = blk$NES,
        source_p          = blk$p,
        source_padj       = blk$padj,
        source_n_measured = blk$n_measured,
        target_layer      = to,
        target_status     = ifelse(tested, "tested",
                                   ifelse(is.na(blk$contrast_key),
                                          "not in target results",
                                          "not in target results for this contrast")),
        target_method     = pick(tgt, t_idx, "method", NA_character_),
        target_contrast   = pick(tgt, t_idx, "contrast", NA_character_),
        target_NES        = pick(tgt, t_idx, "NES", NA_real_),
        target_p          = pick(tgt, t_idx, "p", NA_real_),
        target_padj       = pick(tgt, t_idx, "padj", NA_real_),
        target_n_measured = pick(tgt, t_idx, "n_measured", NA_real_),
        target_ora_p      = pick(tgt_ora, o_idx, "p", NA_real_),
        target_ora_padj   = pick(tgt_ora, o_idx, "padj", NA_real_),
        stringsAsFactors  = FALSE
    )
}


#' Heatmap of one ranking of a lookup
#'
#' One row per pathway, in rank order, and two columns: the ranked layer and the
#' layer it was read in. Fill is the NES on a diverging scale whose centre is a
#' mid grey, so a row reads as "same direction" or "opposite" at a glance; each
#' cell also prints its NES, with a star where the raw p-value is below 0.05.
#'
#' A cell with no NES is left unfilled, outlined only: "ORA" where the layer has
#' only an ORA result, a dash where it has no result for the pathway at all.
#' The two greys used to be confusable -- a NES near zero sat at a near-white
#' midpoint that read as the missing-value grey, and only its number told them
#' apart.
#'
#' @param blk Rows of one ranking from \code{build_cross_omics_lookup()}, in
#'   rank order.
#' @return A ggplot object.
.lookup_heatmap <- function(blk) {
    from <- blk$source_layer[1]
    to <- blk$target_layer[1]
    labels <- disambiguate_pathway_labels(truncate_pathway_label(blk$pathway, 45),
                                          blk$norm_id)
    src_col <- paste0(from, "\n(ranked)")
    tgt_col <- paste0(to, "\n(looked up)")

    cell_text <- function(nes, p, tested) {
        ifelse(!tested, "\u2013",
               ifelse(is.na(nes), "ORA",
                      paste0(sprintf("%.1f", nes),
                             ifelse(!is.na(p) & p < 0.05, "*", ""))))
    }
    long <- rbind(
        data.frame(layer = src_col, pathway = labels, NES = blk$source_NES,
                   label = cell_text(blk$source_NES, blk$source_p,
                                     rep(TRUE, nrow(blk))),
                   stringsAsFactors = FALSE),
        data.frame(layer = tgt_col, pathway = labels, NES = blk$target_NES,
                   label = cell_text(blk$target_NES, blk$target_p,
                                     blk$target_status == "tested"),
                   stringsAsFactors = FALSE)
    )
    long$layer <- factor(long$layer, levels = c(src_col, tgt_col))
    long$pathway <- factor(long$pathway, levels = rev(labels))
    lim <- max(1, abs(long$NES), na.rm = TRUE)
    method <- unique(stats::na.omit(blk$source_method))

    ggplot2::ggplot(long, ggplot2::aes(x = layer, y = pathway, fill = NES)) +
        ggplot2::geom_tile(colour = "grey70") +
        ggplot2::geom_text(ggplot2::aes(label = label), size = 2.6,
                           colour = "grey15") +
        ggplot2::scale_fill_gradient2(
            name = "NES", low = "#2166ac", mid = "#a6a6a6", high = "#b2182b",
            midpoint = 0, limits = c(-lim, lim), na.value = NA) +
        ggplot2::scale_x_discrete(position = "top") +
        ggplot2::labs(
            title = sprintf("Top %s pathways and their %s results", from, to),
            subtitle = sprintf("Ranked by %s raw p-value (%s). Nominal p-values.",
                               from, paste(method, collapse = "/")),
            # Wrapped by hand: one line of this length is cut off at this
            # figure's width.
            caption = paste0("Colour and number: NES, grey at zero; ",
                             "* raw p < 0.05. Both layers read in the same contrast.\n",
                             "An unfilled cell has no NES: ORA only, or a dash ",
                             "for no result in that layer and contrast."),
            x = NULL, y = NULL) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(size = 7),
            axis.text.x = ggplot2::element_text(size = 9),
            plot.title = ggplot2::element_text(hjust = 0.5),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
            plot.caption = ggplot2::element_text(size = 8, colour = "grey30"))
}


#' Draw one ranking of a lookup as a heatmap
#'
#' @param lookup Result of \code{build_cross_omics_lookup()}.
#' @param ranking Which ranking to draw ("all" or "tested_in_target").
#' @param out_path PNG path to write.
#' @return Invisibly, \code{out_path}, or NULL when the ranking has no rows.
plot_cross_omics_lookup <- function(lookup, ranking = "all", out_path) {
    blk <- lookup[lookup$ranking == ranking, , drop = FALSE]
    if (nrow(blk) == 0) return(invisible(NULL))
    blk <- blk[order(blk$rank), , drop = FALSE]

    ggplot2::ggsave(out_path, plot = .lookup_heatmap(blk), width = 7,
                    height = max(4, 2.2 + nrow(blk) * 0.3), dpi = 300)
    invisible(out_path)
}


#' Write every direction's lookup table and figures
#'
#' One table per ordered pair of layers,
#' \code{cross_lookup_<from>_to_<to>.tsv}, holding both rankings, and one
#' figure per ranking: \code{cross_lookup_<from>_to_<to>.png} for the source's
#' top pathways and \code{..._tested.png} for those the target tested.
#'
#' @param pathway_tables Named list of per-omics enrichment data frames, as
#'   merged by the meta-analysis.
#' @param omics Layer names, in order.
#' @param out_dir Directory to write into.
#' @param top_n Number of pathways per ranking.
#' @param kegg_org Active KEGG organism code for the run, or NULL.
#' @return Named list, one element per direction written, each a named
#'   character vector of the paths written (tsv, all, tested_in_target).
write_cross_omics_lookups <- function(pathway_tables, omics, out_dir,
                                      top_n = 15, kegg_org = NULL) {
    written <- list()
    for (from in omics) {
        for (to in setdiff(omics, from)) {
            lk <- build_cross_omics_lookup(pathway_tables, from, to, top_n, kegg_org)
            if (is.null(lk) || nrow(lk) == 0) next

            stem <- file.path(out_dir, paste0("cross_lookup_", .collection_slug(from),
                                              "_to_", .collection_slug(to)))
            paths <- c(tsv = paste0(stem, ".tsv"))
            utils::write.table(lk, paths[["tsv"]], sep = "\t", quote = FALSE,
                               row.names = FALSE)

            for (rk in c("all", "tested_in_target")) {
                png_path <- paste0(stem, if (rk == "all") "" else "_tested", ".png")
                drawn <- tryCatch(
                    plot_cross_omics_lookup(lk, rk, png_path),
                    error = function(e) {
                        message("  Cross-omics lookup figure failed (", from, " to ",
                                to, ", ", rk, "): ", e$message)
                        NULL
                    })
                if (!is.null(drawn) && file.exists(png_path)) paths[[rk]] <- png_path
            }
            written[[paste0(from, "_to_", to)]] <- paths
        }
    }
    if (length(written) > 0) {
        message("  Cross-omics lookups written: ", paste(names(written), collapse = ", "))
    }
    written
}


#' Remove lookup tables and figures a previous run left in this directory
#'
#' The report finds these by glob, so a direction or ranking this run does not
#' produce -- a layer dropped, the lookup switched off -- would otherwise be
#' shown as this run's. Same reasoning, and same narrow scope, as
#' \code{.clear_collection_heatmaps()}.
#'
#' @param out_dir Directory the current run is about to write into.
#' @return Invisibly, the paths removed.
#' @keywords internal
.clear_cross_lookup_outputs <- function(out_dir) {
    if (is.null(out_dir) || !dir.exists(out_dir)) return(invisible(character(0)))
    stale <- list.files(out_dir, pattern = "^cross_lookup_.*\\.(tsv|png)$",
                        full.names = TRUE)
    if (length(stale) > 0) unlink(stale)
    invisible(stale)
}


#' Lookup settings from the config, with defaults
#'
#' Off unless the config turns it on -- it is exploratory, and writes a table and
#' two figures per ordered pair of layers; \code{top_n} falls back to 15 when
#' absent or unusable. The config validator fills the same defaults, and they
#' are repeated here so a config that bypassed it behaves the same.
#'
#' @param config Full config object.
#' @return List with \code{enabled} (logical) and \code{top_n} (integer).
#' @keywords internal
.cross_lookup_config <- function(config) {
    cfg <- config$modes$multiomics$enrichment$cross_lookup %||% list()
    enabled <- cfg$enabled %||% FALSE
    top_n <- suppressWarnings(as.integer(cfg$top_n %||% 15L))
    if (length(top_n) != 1 || is.na(top_n) || top_n < 1) top_n <- 15L
    list(enabled = isTRUE(enabled), top_n = top_n)
}
