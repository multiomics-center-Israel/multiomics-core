# R/domain/metabolomics/06f_pathway_volcano.R
#
# Volcano plot with an enrichment-pathway selector: pick a significant GSEA or
# QEA pathway and its measured metabolites are colored on the volcano.
#
#   build_pathway_feature_members() — pathway label -> measured feature ids
#   select_overlay_pathways()       — significant pathways, in selector order
#   plot_volcano_pathway_overlay()  — plotly volcano with the pathway dropdown
#
# Reuses: read_gmt_list, translate_gmt_hmdb_to_kegg, make_pathway_labels,
#   map_compounds_for_enrichment, resolve_input_path (06_enrichment.R, core).


#' Map enrichment pathways to the measured features they contain
#'
#' Reads the configured GMT file(s) the way GSEA and QEA do (HMDB to KEGG
#' translation, "ID - Description" labels), then maps each pathway's compounds
#' to feature ids through the same per-feature compound map. Compound IDs are
#' matched case-insensitively, and every feature carrying a member compound is
#' included, so duplicate annotations of one compound all count.
#'
#' @param pre    Preprocessing list (needs \code{row_data} and \code{expr_raw}).
#' @param config Full pipeline config (reads \code{modes$metabolomics$enrichment}).
#' @return list(members, feature_compound): \code{members} is a named list,
#'   pathway label -> character vector of feature ids, holding only pathways
#'   with at least one measured feature; \code{feature_compound} is a named
#'   character vector, feature id -> compound id. Both are empty when no GMT
#'   or no compound mapping is available.
build_pathway_feature_members <- function(pre, config) {
    empty   <- list(members = list(), feature_compound = character(0))
    enr_cfg <- config$modes$metabolomics$enrichment %||% list()

    # as.character(): no gmt_file resolves to NULL, and file.exists(NULL) errors.
    gmt_files <- as.character(resolve_input_path(config, enr_cfg$gmt_file))
    gmt_files <- gmt_files[!is.na(gmt_files) & nzchar(gmt_files) & file.exists(gmt_files)]
    if (length(gmt_files) == 0 || is.null(pre$row_data) || is.null(pre$expr_raw)) {
        return(empty)
    }
    mapping_file <- resolve_input_path(config, enr_cfg$mapping_file)

    mapped <- map_compounds_for_enrichment(pre$row_data, pre$expr_raw, mapping_file)
    fmap   <- mapped$feature_map
    fmap   <- fmap[!is.na(fmap) & nzchar(fmap)]
    if (length(fmap) == 0) return(empty)

    # compound id (lower case) -> feature ids carrying it
    by_compound <- split(names(fmap), tolower(unname(fmap)))

    members <- list()
    for (gf in gmt_files) {
        parsed <- read_gmt_list(gf, include_descriptions = TRUE)
        sets   <- translate_gmt_hmdb_to_kegg(parsed$sets, mapping_file)
        if (length(sets) == 0) next
        names(sets) <- make_pathway_labels(names(sets), parsed$descriptions)
        for (pw in names(sets)) {
            feats <- unique(unlist(by_compound[tolower(sets[[pw]])], use.names = FALSE))
            if (length(feats) > 0) members[[pw]] <- feats
        }
    }

    list(members = members, feature_compound = fmap)
}


#' Significant GSEA and QEA pathways for the volcano selector
#'
#' Keeps pathways with FDR below \code{fdr_cutoff} that have at least one
#' measured feature in \code{members}. GSEA pathways come first, ordered by
#' |NES| (highest first, ties by FDR), then QEA pathways ordered by FDR
#' (lowest first).
#'
#' @param gsea_tbl   GSEA table (pathway, FDR, NES, optional leading_edge), or NULL.
#' @param qea_tbl    QEA table (pathway, FDR), or NULL.
#' @param members    Named list from \code{build_pathway_feature_members()$members}.
#' @param fdr_cutoff Significance cutoff on FDR (default 0.05).
#' @return data.frame with columns method, pathway, FDR, NES, leading_edge,
#'   n_measured, label; zero rows when nothing passes.
select_overlay_pathways <- function(gsea_tbl = NULL, qea_tbl = NULL,
                                    members = list(), fdr_cutoff = 0.05) {
    empty <- data.frame(
        method = character(0), pathway = character(0), FDR = numeric(0),
        NES = numeric(0), leading_edge = character(0), n_measured = integer(0),
        label = character(0), stringsAsFactors = FALSE
    )

    pick <- function(tbl, method) {
        if (is.null(tbl) || nrow(tbl) == 0 ||
            !all(c("pathway", "FDR") %in% colnames(tbl))) {
            return(empty)
        }
        fdr  <- as.numeric(tbl$FDR)
        keep <- !is.na(fdr) & fdr < fdr_cutoff &
            as.character(tbl$pathway) %in% names(members)
        tbl <- tbl[keep, , drop = FALSE]
        if (nrow(tbl) == 0) return(empty)

        fdr <- as.numeric(tbl$FDR)
        nes <- if ("NES" %in% colnames(tbl)) as.numeric(tbl$NES) else rep(NA_real_, nrow(tbl))
        ord <- if (method == "GSEA") order(-abs(nes), fdr) else order(fdr)
        tbl <- tbl[ord, , drop = FALSE]
        fdr <- fdr[ord]
        nes <- nes[ord]

        pw <- as.character(tbl$pathway)
        le <- if ("leading_edge" %in% colnames(tbl)) {
            as.character(tbl$leading_edge)
        } else {
            rep(NA_character_, nrow(tbl))
        }
        stats_txt <- if (method == "GSEA") {
            sprintf("NES %.2f, FDR %s", nes, signif(fdr, 2))
        } else {
            sprintf("FDR %s", signif(fdr, 2))
        }
        data.frame(
            method       = method,
            pathway      = pw,
            FDR          = fdr,
            NES          = nes,
            leading_edge = le,
            n_measured   = vapply(pw, function(p) length(members[[p]]), integer(1),
                                  USE.NAMES = FALSE),
            label        = sprintf("%s | %s (%s)", method, pw, stats_txt),
            stringsAsFactors = FALSE
        )
    }

    out <- rbind(pick(gsea_tbl, "GSEA"), pick(qea_tbl, "QEA"))
    rownames(out) <- NULL
    out
}


#' Interactive volcano with a pathway selector
#'
#' Draws every feature in grey, plus one hidden trace per pathway holding that
#' pathway's measured features in red. A dropdown shows one pathway at a time,
#' in the order of \code{pathways}. With no pathways the dropdown holds a
#' single entry that does nothing, so nothing can be selected. Pathways whose
#' members all lack a fold change or p-value are left out of the dropdown.
#'
#' @param de_tbl     Contrast table with feature_id, logFC, the \code{p_col}
#'   column and optionally Name.
#' @param pathways   data.frame from \code{select_overlay_pathways()}.
#' @param members    Named list, pathway label -> feature ids.
#' @param feature_compound Named character, feature id -> compound id; used to
#'   report GSEA leading-edge membership on hover. NULL skips that line.
#' @param p_col      Column plotted on the y axis ("adj.P.Val" or "P.Value").
#' @param p_cutoff   Significance cutoff drawn as a horizontal line.
#' @param lfc_cutoff |log2 FC| cutoff drawn as vertical lines.
#' @param fdr_cutoff Pathway FDR cutoff, quoted in the empty-selector text.
#' @return A plotly htmlwidget, or NULL when plotly or the needed columns are missing.
plot_volcano_pathway_overlay <- function(de_tbl, pathways, members,
                                         feature_compound = NULL,
                                         p_col = "adj.P.Val", p_cutoff = 0.05,
                                         lfc_cutoff = log2(1.5),
                                         fdr_cutoff = 0.05) {
    if (!requireNamespace("plotly", quietly = TRUE)) return(NULL)
    if (is.null(de_tbl) || !all(c("feature_id", "logFC", p_col) %in% colnames(de_tbl))) {
        return(NULL)
    }
    df <- de_tbl[!is.na(de_tbl$logFC) & !is.na(de_tbl[[p_col]]), , drop = FALSE]
    if (nrow(df) == 0) return(NULL)

    fid  <- as.character(df$feature_id)
    name <- if ("Name" %in% colnames(df)) as.character(df$Name) else fid
    no_name <- is.na(name) | !nzchar(name)
    name[no_name] <- fid[no_name]
    y     <- -log10(pmax(as.numeric(df[[p_col]]), 1e-300))
    hover <- sprintf("Metabolite: %s<br>log2(FC): %.3f<br>%s: %s",
                     name, df$logFC, p_col, signif(as.numeric(df[[p_col]]), 3))

    # Selector entries with at least one plotted feature
    if (is.null(pathways)) pathways <- select_overlay_pathways(members = members)
    member_idx <- lapply(pathways$pathway, function(pw) which(fid %in% members[[pw]]))
    has_pts    <- lengths(member_idx) > 0
    pathways   <- pathways[has_pts, , drop = FALSE]
    member_idx <- member_idx[has_pts]
    n_pw <- nrow(pathways)

    fig <- plotly::plot_ly(
        x = df$logFC, y = y, type = "scatter", mode = "markers",
        text = hover, hoverinfo = "text",
        marker = list(color = "rgba(160,160,160,0.5)", size = 6),
        name = "All metabolites", showlegend = FALSE
    )

    for (k in seq_len(n_pw)) {
        idx <- member_idx[[k]]
        le_txt <- ""
        le_raw <- pathways$leading_edge[k]
        if (identical(pathways$method[k], "GSEA") && !is.null(feature_compound) &&
            !is.na(le_raw) && nzchar(le_raw)) {
            le <- tolower(strsplit(le_raw, ";", fixed = TRUE)[[1]])
            in_le <- tolower(unname(feature_compound[fid[idx]])) %in% le
            le_txt <- ifelse(in_le, "<br>Leading edge: yes", "<br>Leading edge: no")
        }
        fig <- plotly::add_trace(
            fig,
            x = df$logFC[idx], y = y[idx], type = "scatter", mode = "markers",
            text = paste0(hover[idx], "<br>Pathway: ", pathways$pathway[k], le_txt),
            hoverinfo = "text",
            marker = list(color = "#b2182b", size = 9,
                          line = list(color = "black", width = 0.5)),
            name = pathways$pathway[k], visible = (k == 1),
            showlegend = FALSE, inherit = FALSE
        )
    }

    short_label <- function(x, max_len = 90) {
        ifelse(nchar(x) > max_len, paste0(substr(x, 1, max_len - 3), "..."), x)
    }
    buttons <- if (n_pw > 0) {
        lapply(seq_len(n_pw), function(k) {
            list(method = "restyle",
                 args   = list("visible", as.list(c(TRUE, seq_len(n_pw) == k))),
                 label  = short_label(pathways$label[k]))
        })
    } else {
        list(list(method = "skip", args = list(),
                  label = sprintf("No GSEA or QEA pathway with FDR < %s", fdr_cutoff)))
    }

    dashed <- list(color = "grey", width = 1, dash = "dash")
    shapes <- list(
        list(type = "line", xref = "paper", x0 = 0, x1 = 1,
             y0 = -log10(p_cutoff), y1 = -log10(p_cutoff), line = dashed),
        list(type = "line", yref = "paper", y0 = 0, y1 = 1,
             x0 = -lfc_cutoff, x1 = -lfc_cutoff, line = dashed),
        list(type = "line", yref = "paper", y0 = 0, y1 = 1,
             x0 = lfc_cutoff, x1 = lfc_cutoff, line = dashed)
    )

    plotly::layout(
        fig,
        xaxis  = list(title = "log2(FC)", zeroline = FALSE),
        yaxis  = list(title = paste0("-log10(", p_col, ")")),
        shapes = shapes,
        margin = list(t = 80),
        updatemenus = list(list(
            type = "dropdown", active = 0,
            x = 0, xanchor = "left", y = 1.15, yanchor = "top",
            buttons = buttons
        ))
    )
}
