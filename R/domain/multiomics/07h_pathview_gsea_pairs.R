#' KEGG maps for pathways chosen by GSEA and by enzyme-metabolite pairs
#'
#' The Multi-ORA renderers draw only pathways that ORA called enriched, and ORA
#' sees only the short list of DE features. A pathway that rank-based GSEA calls
#' enriched, or that holds a changed enzyme paired with a measured metabolite,
#' then has no map at all -- even though the report discusses it elsewhere.
#' This file selects those pathways and draws them, in their own PDF, so each
#' map's reason for being there stays stated beside it.


#' KEGG global and overview maps (01100-01320)
#'
#' pathview cannot overlay data on these, and every metabolic enzyme sits in
#' 01100, so a pair table would otherwise ask for it on almost every row.
#'
#' @keywords internal
.GLOBAL_KEGG_MAP_PATTERN <- "^01[1-3][0-9]{2}$"


#' Collect the KEGG pathways of changed enzymes in the enzyme-metabolite table
#'
#' A row qualifies when its enzyme is a DE hit; the metabolite side is not
#' required, because a pathway with a changed enzyme is worth seeing whether or
#' not its paired metabolite moved. A row may list several pathways separated
#' by ";". Global maps are dropped (see \code{.GLOBAL_KEGG_MAP_PATTERN}).
#'
#' Pure: reads no files and draws nothing.
#'
#' @param pairs Data frame from \code{build_enzyme_metabolite_pairs()}, with
#'   `contrast`, `enzyme_hit`, `pathway_id` and optionally `enzyme_padj`,
#'   `enzyme_p` and `gene_symbol`.
#' @return Named list keyed by canonical contrast key, each element a list of
#'   `label` (first raw contrast spelling), `pathways` (bare KEGG map numbers,
#'   ranked by the best enzyme p-value on them) and `genes` (named by pathway,
#'   the changed enzymes' symbols joined with ", "). Empty list when nothing
#'   qualifies.
pair_pathways_by_contrast <- function(pairs) {
    need <- c("contrast", "enzyme_hit", "pathway_id")
    if (is.null(pairs) || !is.data.frame(pairs) || nrow(pairs) == 0 ||
        !all(need %in% names(pairs))) {
        return(list())
    }

    hit <- as.logical(pairs$enzyme_hit)
    raw_contrast <- as.character(pairs$contrast)
    keep <- !is.na(hit) & hit & !is.na(pairs$pathway_id) &
        nzchar(as.character(pairs$pathway_id)) &
        !is.na(raw_contrast) & nzchar(trimws(raw_contrast))
    if (!any(keep)) return(list())
    p <- pairs[keep, , drop = FALSE]
    raw_contrast <- raw_contrast[keep]

    # Rank on the enzyme's own evidence; padj where present, else raw p.
    score <- if ("enzyme_padj" %in% names(p) && any(!is.na(p$enzyme_padj))) {
        suppressWarnings(as.numeric(p$enzyme_padj))
    } else if ("enzyme_p" %in% names(p)) {
        suppressWarnings(as.numeric(p$enzyme_p))
    } else {
        rep(NA_real_, nrow(p))
    }
    genes <- if ("gene_symbol" %in% names(p)) {
        as.character(p$gene_symbol)
    } else {
        rep(NA_character_, nrow(p))
    }

    ids <- strsplit(as.character(p$pathway_id), ";", fixed = TRUE)
    long <- data.frame(
        row = rep(seq_len(nrow(p)), lengths(ids)),
        pathway = normalize_kegg_pathway_id(trimws(unlist(ids))),
        stringsAsFactors = FALSE
    )
    long <- long[nzchar(long$pathway) &
                 grepl("^[0-9]{5}$", long$pathway) &
                 !grepl(.GLOBAL_KEGG_MAP_PATTERN, long$pathway), , drop = FALSE]
    if (nrow(long) == 0) return(list())
    long$ckey <- normalize_contrast_key(raw_contrast[long$row])
    long$score <- score[long$row]
    long$gene <- genes[long$row]
    long$label <- raw_contrast[long$row]

    out <- list()
    for (k in unique(long$ckey)) {
        d <- long[long$ckey == k, , drop = FALSE]
        best <- tapply(ifelse(is.na(d$score), Inf, d$score), d$pathway, min)
        ord <- order(as.numeric(best), names(best))
        pw <- names(best)[ord]
        gene_txt <- vapply(pw, function(id) {
            g <- sort(unique(stats::na.omit(d$gene[d$pathway == id])))
            paste(g, collapse = ", ")
        }, character(1))
        out[[k]] <- list(label = d$label[1], pathways = pw,
                         genes = stats::setNames(gene_txt, pw))
    }
    out
}


#' Decide which pathways get a GSEA/pair map, and why
#'
#' GSEA rows come from this run's per-omics enrichment frames, KEGG accessions
#' only, rank-based methods only, adjusted p below the shared
#' \code{.PATHVIEW_THRESHOLDS$fdr_alpha}, capped at `top_n` per contrast. Pair
#' pathways are all kept: each one names a specific changed enzyme the reader
#' was pointed at. A pathway reached both ways is listed once with both reasons.
#'
#' Pure: reads no files and draws nothing.
#'
#' @param per_omics_enrichment Named list of per-omics enrichment frames
#'   (`multiomics_cross_enrichment$per_omics`).
#' @param enzyme_pairs Enzyme-metabolite pair table, or NULL.
#' @param kegg_org KEGG organism code of the run (e.g. "rno").
#' @param top_n Maximum GSEA pathways per contrast.
#' @return Data frame with one row per contrast and pathway: `ckey`, `contrast`,
#'   `pathway` (bare map number), `source` ("GSEA", "enzyme-metabolite pair" or
#'   both joined by " + "), `gsea_padj` and `pair_enzymes`. Zero rows when
#'   nothing qualifies.
select_gsea_pair_pathways <- function(per_omics_enrichment, enzyme_pairs,
                                      kegg_org, top_n = 10) {
    empty <- data.frame(ckey = character(0), contrast = character(0),
                        pathway = character(0), source = character(0),
                        gsea_padj = numeric(0), pair_enzymes = character(0),
                        stringsAsFactors = FALSE)

    gsea <- if (length(per_omics_enrichment) > 0) {
        suppressMessages(.kegg_hits_by_contrast(per_omics_enrichment, kegg_org,
                                                methods = c("fgsea", "gsea")))
    } else {
        list()
    }
    pairs <- pair_pathways_by_contrast(enzyme_pairs)

    rows <- list()
    for (k in union(names(gsea), names(pairs))) {
        g <- gsea[[k]]
        g_ids <- if (is.null(g)) character(0) else {
            ids <- g$pathways[!grepl(.GLOBAL_KEGG_MAP_PATTERN, g$pathways)]
            utils::head(ids, top_n)
        }
        p_ids <- if (is.null(pairs[[k]])) character(0) else pairs[[k]]$pathways
        ids <- c(g_ids, setdiff(p_ids, g_ids))
        if (length(ids) == 0) next
        in_g <- ids %in% g_ids
        in_p <- ids %in% p_ids
        rows[[k]] <- data.frame(
            ckey = k,
            contrast = if (!is.null(g)) g$label else pairs[[k]]$label,
            pathway = ids,
            source = ifelse(in_g & in_p, "GSEA + enzyme-metabolite pair",
                            ifelse(in_g, "GSEA", "enzyme-metabolite pair")),
            gsea_padj = if (is.null(g)) NA_real_ else unname(g$scores[ids]),
            pair_enzymes = if (is.null(pairs[[k]])) NA_character_ else
                unname(pairs[[k]]$genes[ids]),
            stringsAsFactors = FALSE
        )
    }
    if (length(rows) == 0) return(empty)
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}


#' Render KEGG maps for GSEA- and enzyme-pair-selected pathways
#'
#' Gene nodes carry proteomics (and transcriptomics, when present) log2FC and
#' compound nodes metabolomics log2FC, each from the contrast the pathway was
#' selected in, and each filtered by \code{filter_changed_features()} -- the
#' same node rule, and so the same caption, as the other pathview renderers.
#'
#' Writes `pathview_gsea_and_pairs.pdf` (one captioned page per map) and
#' `pathview_gsea_and_pairs_selection.tsv` (why each map is there) into
#' `out_dir`, and the PNGs into `out_dir/pathview` with suffix `gsea_pair_*`.
#'
#' @param per_omics_enrichment Named list of per-omics enrichment frames.
#' @param enzyme_pairs Enzyme-metabolite pair table, or NULL.
#' @param de_results Named list of DE results per omics.
#' @param harmonization_res Harmonization result.
#' @param config Full config.
#' @param out_dir Multi-ORA output directory.
#' @return Path to the compiled PDF, or NULL when nothing was drawn.
generate_gsea_pair_pathview <- function(per_omics_enrichment, enzyme_pairs,
                                        de_results, harmonization_res,
                                        config, out_dir) {
    if (!requireNamespace("pathview", quietly = TRUE)) return(NULL)

    organism <- config$global$organism %||% "human"
    kegg_org <- get_kegg_organism(organism)
    org_db <- get_organism_db(organism)
    if (is.null(kegg_org) || is.null(org_db)) return(NULL)

    pv_cfg <- config$modes$multiomics$enrichment$pathview %||% list()
    sel <- select_gsea_pair_pathways(per_omics_enrichment, enzyme_pairs,
                                     kegg_org, top_n = pv_cfg$top_n %||% 10)
    excl <- .excluded_pathway_classes(config)
    if (nrow(sel) > 0 && length(excl) > 0) {
        sel <- sel[keep_kegg_pathways(sel$pathway, exclude = excl,
                                      kegg_org = kegg_org,
                                      label = "GSEA/pair pathview maps"), ,
                   drop = FALSE]
    }
    if (nrow(sel) == 0) {
        message("  GSEA/pair pathview: no GSEA-enriched or enzyme-pair KEGG ",
                "pathways to draw")
        return(NULL)
    }
    message("  GSEA/pair pathview: ", nrow(sel), " pathway map(s) selected")

    gene_tables <- list()
    gene_maps <- list()
    for (om in intersect(c("transcriptomics", "proteomics"), names(de_results))) {
        tabs <- tryCatch(extract_de_tables(de_results[[om]], om, harmonization_res),
                         error = function(e) NULL)
        if (is.null(tabs) || length(tabs) == 0) next
        id_map <- tryCatch(map_feature_ids_to_entrez(tabs, om, harmonization_res,
                                                     org_db),
                           error = function(e) NULL)
        if (is.null(id_map) || nrow(id_map) == 0) next
        gene_tables[[om]] <- tabs
        gene_maps[[om]] <- id_map
    }
    metab_tables <- NULL
    metab_map <- NULL
    if ("metabolomics" %in% names(de_results)) {
        metab_tables <- tryCatch(
            extract_de_tables(de_results$metabolomics, "metabolomics",
                              harmonization_res),
            error = function(e) NULL)
        if (!is.null(metab_tables)) {
            metab_map <- tryCatch(map_metabolite_ids_to_kegg(metab_tables,
                                                             harmonization_res),
                                  error = function(e) NULL)
        }
    }

    if (!exists("bods", envir = globalenv())) {
        utils::data("bods", package = "pathview", envir = globalenv())
    }
    pv_dir <- file.path(out_dir, "pathview")
    dir.create(pv_dir, recursive = TRUE, showWarnings = FALSE)
    pv_dir <- normalizePath(pv_dir, winslash = "/", mustWork = FALSE)

    sel$drawn <- FALSE
    page_labels <- character(0)
    generated <- withr::with_dir(pv_dir, {
        made <- character(0)
        for (k in unique(sel$ckey)) {
            gene_fc <- list()
            for (om in names(gene_tables)) {
                df <- filter_changed_features(
                    .de_table_for_contrast(gene_tables[[om]], k))
                if (is.null(df) || nrow(df) == 0) next
                m <- merge(df, gene_maps[[om]], by = "feature_id")
                fc <- tapply(m$log2fc, m$ENTREZID, mean, na.rm = TRUE)
                gene_fc[[om]] <- stats::setNames(as.numeric(fc), names(fc))
            }
            gene_data <- NULL
            if (length(gene_fc) == 1) {
                gene_data <- gene_fc[[1]]
            } else if (length(gene_fc) == 2) {
                ids <- unique(unlist(lapply(gene_fc, names)))
                gene_data <- matrix(NA_real_, length(ids), 2,
                                    dimnames = list(ids, names(gene_fc)))
                for (om in names(gene_fc)) gene_data[names(gene_fc[[om]]), om] <- gene_fc[[om]]
            }

            cpd_data <- NULL
            mdf <- filter_changed_features(.de_table_for_contrast(metab_tables, k))
            if (!is.null(mdf) && nrow(mdf) > 0 && !is.null(metab_map) &&
                nrow(metab_map) > 0) {
                md <- merge(mdf, metab_map, by = "feature_id")
                ok <- is.finite(md$log2fc)
                if (any(ok)) {
                    fc <- tapply(md$log2fc[ok], md$KEGG_CPD[ok], mean, na.rm = TRUE)
                    cpd_data <- stats::setNames(as.numeric(fc), names(fc))
                }
            }
            if (is.null(gene_data) && is.null(cpd_data)) {
                message("  GSEA/pair pathview: no feature passed the node ",
                        "thresholds for contrast ", k, "; skipping it")
                next
            }

            out_suffix <- paste0("gsea_pair_", .contrast_out_key(k))
            two_layers <- !is.null(dim(gene_data)) && ncol(gene_data) > 1
            for (i in which(sel$ckey == k)) {
                pid <- sel$pathway[i]
                tryCatch({
                    pathview::pathview(gene.data = gene_data, cpd.data = cpd_data,
                                       pathway.id = pid, species = kegg_org,
                                       out.suffix = out_suffix, kegg.dir = pv_dir,
                                       keys.align = "y", match.data = TRUE,
                                       multi.state = two_layers,
                                       same.layer = FALSE)
                    f <- c(paste0(kegg_org, pid, ".", out_suffix, ".multi.png"),
                           paste0(kegg_org, pid, ".", out_suffix, ".png"))
                    f <- f[file.exists(f)]
                    if (length(f) > 0) {
                        made <- c(made, file.path(pv_dir, f[1]))
                        sel$drawn[i] <- TRUE
                        page_labels <- c(page_labels, paste0(
                            sel$contrast[i], "  |  ", kegg_org, pid,
                            "  |  selected by: ", sel$source[i]))
                        message("    GSEA/pair pathview: ", kegg_org, pid,
                                " (", sel$source[i], ")")
                    }
                }, error = function(e) {
                    message("    GSEA/pair pathview failed for ", pid, ": ",
                            conditionMessage(e))
                })
            }
        }
        made
    })

    # The selection is written even when a render failed, so a missing map can
    # be told apart from a pathway that was never selected.
    sel$pathway <- paste0(kegg_org, sel$pathway)
    utils::write.table(sel[, setdiff(names(sel), "ckey")],
                       file.path(out_dir, "pathview_gsea_and_pairs_selection.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE, na = "")

    if (length(generated) == 0) return(NULL)
    pdf_path <- .compile_pathview_pdf(
        generated, file.path(out_dir, "pathview_gsea_and_pairs.pdf"), page_labels)
    if (!is.null(pdf_path)) {
        message("  GSEA/pair pathview: ", length(generated), " maps -> ",
                basename(pdf_path))
    }
    pdf_path
}
