#' GO hierarchy of the layers' GSEA results
#'
#' Draws, for each gene-based layer, contrast, GO ontology and NES sign, the
#' part of the GO graph above the GO terms that the layer's own fgsea step
#' found at padj below the configured cutoff. No enrichment is computed here:
#' the seeds and their values are read from the existing results
#' (prot_pathway_res, rna_pathway_res), and the hierarchy comes from GO.db.
#' Output: multiomics/go_hierarchy/<layer>/<contrast>/<ONT>_<direction>.*,
#' an index the report reads, and <mode>/Datasets/go_dag_edges.tsv.

# Layers that carry gene-level GO fgsea results, and the mode folder each
# writes its Datasets/ export into.
.go_hierarchy_layer_modes <- c(proteomics = "proteomics", transcriptomics = "rna")

#' Resolve the GO hierarchy settings
#'
#' The one reader of \code{modes.multiomics.go_hierarchy}: the targets, the
#' writer and the report all call it, so a default is decided in one place.
#'
#' @param config Full config list (or one whose \code{modes$multiomics} is the
#'   multiomics block).
#' @return list(enabled, ontologies (subset of "BP", "MF", "CC"), fdr,
#'   max_seed_terms, max_nodes).
go_hierarchy_config <- function(config) {
    gh <- (config$modes$multiomics %||% list())$go_hierarchy %||% list()
    onts <- unlist(gh$ontologies %||% "GO_BP", use.names = FALSE)
    list(
        enabled        = isTRUE(gh$enabled),
        ontologies     = unique(sub("^GO_", "", toupper(as.character(onts)))),
        fdr            = as.numeric(gh$fdr %||% 0.05),
        max_seed_terms = as.integer(gh$max_seed_terms %||% 30L),
        max_nodes      = as.integer(gh$max_nodes %||% 120L)
    )
}

#' Collect a layer's GO fgsea tables in one long table
#'
#' Reads the contrast-keyed result that run_pathway_analysis() returns
#' (proteomics: the list itself; RNA online path: its \code{pathway_results}
#' slot). Only fgsea tables of GO collections are taken: \code{GO_BP},
#' \code{GO_MF}, \code{GO_CC}, or \code{GO} (the mixed collection built for
#' non-model organisms, split here by each term's GO.db ontology). Other
#' collections, ORA tables and the clusterProfiler GSEA of the offline RNA
#' path are not fgsea-on-GO results and are skipped.
#'
#' @param pathway_res A layer's pathway result.
#' @param ontology_of Function mapping GO ids to "BP"/"MF"/"CC" (NA when
#'   unknown); used only for the mixed "GO" collection.
#' @return data.frame with \code{contrast}, \code{ontology}, \code{go_id},
#'   \code{term}, \code{NES}, \code{padj}, \code{pval}, \code{size},
#'   \code{database}; zero rows when the layer has none.
collect_go_fgsea_tables <- function(pathway_res,
                                    ontology_of = function(ids) go_term_info(ids)$ontology) {
    empty <- data.frame(contrast = character(0), ontology = character(0),
                        go_id = character(0), term = character(0),
                        NES = numeric(0), padj = numeric(0), pval = numeric(0),
                        size = integer(0), database = character(0),
                        stringsAsFactors = FALSE)
    if (is.null(pathway_res)) return(empty)
    contrast_list <- if (is.list(pathway_res$pathway_results)) {
        pathway_res$pathway_results
    } else {
        pathway_res
    }
    if (!is.list(contrast_list) || is.data.frame(contrast_list)) return(empty)

    rows <- list()
    for (cn in names(contrast_list)) {
        res <- contrast_list[[cn]]
        if (!is.list(res) || is.data.frame(res)) next
        for (key in names(res)) {
            df <- res[[key]]
            if (!is.data.frame(df) || nrow(df) == 0) next
            is_fgsea <- if ("method" %in% colnames(df)) {
                identical(unique(as.character(df$method)), "fgsea")
            } else {
                grepl("_fgsea$", key)
            }
            db <- if ("database" %in% colnames(df)) as.character(df$database[1]) else
                sub("_fgsea$", "", key)
            if (!is_fgsea || !grepl("^GO(_BP|_MF|_CC)?$", db)) next
            if (!all(c("pathway", "padj", "NES") %in% colnames(df))) next
            go_ids <- as.character(df$pathway)
            keep <- grepl("^GO:[0-9]{7}$", go_ids)
            if (!any(keep)) next
            df <- df[keep, , drop = FALSE]
            go_ids <- go_ids[keep]
            ont <- if (db == "GO") ontology_of(go_ids) else rep(sub("^GO_", "", db), length(go_ids))
            term <- if ("pathway_name" %in% colnames(df)) as.character(df$pathway_name) else go_ids
            term[is.na(term) | !nzchar(term)] <- go_ids[is.na(term) | !nzchar(term)]
            rows[[length(rows) + 1]] <- data.frame(
                contrast = cn, ontology = ont, go_id = go_ids, term = term,
                NES = as.numeric(df$NES), padj = as.numeric(df$padj),
                pval = if ("pval" %in% colnames(df)) as.numeric(df$pval) else NA_real_,
                size = if ("size" %in% colnames(df)) as.integer(df$size) else NA_integer_,
                database = db, stringsAsFactors = FALSE)
        }
    }
    if (length(rows) == 0) return(empty)
    out <- do.call(rbind, rows)
    out <- out[!is.na(out$ontology), , drop = FALSE]
    rownames(out) <- NULL
    out
}

#' Pick the seed terms of one contrast and ontology, split by NES sign
#'
#' Seeds are terms with \code{padj < fdr} and a finite, non-zero NES. Within
#' each sign they are ordered by padj, then raw p, then GO id, and capped at
#' \code{max_seed_terms}. If one term appears twice (the same GO term from two
#' collections), its smallest padj is kept.
#'
#' @param tbl Rows of collect_go_fgsea_tables() for one contrast and ontology.
#' @param fdr Adjusted p-value cutoff.
#' @param max_seed_terms Most seeds per direction.
#' @return list(positive, negative), each list(seeds = data.frame,
#'   n_significant = number of terms passing before the cap).
select_go_dag_seeds <- function(tbl, fdr = 0.05, max_seed_terms = 30L) {
    sig <- tbl[!is.na(tbl$padj) & tbl$padj < fdr &
                   is.finite(tbl$NES) & tbl$NES != 0, , drop = FALSE]
    sig <- sig[order(sig$padj, sig$pval, sig$go_id, na.last = TRUE), , drop = FALSE]
    sig <- sig[!duplicated(sig$go_id), , drop = FALSE]
    pick <- function(rows) {
        list(seeds = utils::head(rows[, c("go_id", "term", "NES", "padj", "pval", "size"),
                                      drop = FALSE], max_seed_terms),
             n_significant = nrow(rows))
    }
    list(positive = pick(sig[sig$NES > 0, , drop = FALSE]),
         negative = pick(sig[sig$NES < 0, , drop = FALSE]))
}

#' Split a contrast name into its two arms
#'
#' Recognises "A_vs_B", "A vs. B" / "Avs.B" (proteomics strips spaces),
#' "A vs B" and "A - B". Anything else returns NULL, and callers fall back to
#' naming the NES sign.
#'
#' @param contrast Contrast name.
#' @return character(2) c(numerator, denominator), or NULL.
split_contrast_arms <- function(contrast) {
    for (sep in c("_vs_", "vs\\.", " vs ", " - ")) {
        parts <- trimws(strsplit(contrast, sep, perl = TRUE)[[1]])
        if (length(parts) == 2 && all(nzchar(parts))) return(parts)
    }
    NULL
}

#' Reader-facing label of one NES direction
#'
#' Positive NES means the set sits towards the top of the layer's ranking,
#' which every gene-based layer orders as numerator over denominator.
#'
#' @param contrast Contrast name.
#' @param direction "positive" or "negative".
#' @return E.g. "higher in Perfused", or "positive NES" when the contrast name
#'   cannot be split.
go_direction_label <- function(contrast, direction = c("positive", "negative")) {
    direction <- match.arg(direction)
    arms <- split_contrast_arms(contrast)
    if (is.null(arms)) return(paste(direction, "NES"))
    paste("higher in", if (direction == "positive") arms[1] else arms[2])
}

#' Build every GO hierarchy figure of the run
#'
#' Pure apart from GO.db reads: returns the graphs and plots, writes nothing.
#'
#' @param enrichment_results Named list of layer pathway results, as in the
#'   multiomics_enrichment_results target; only the gene-based layers in
#'   \code{.go_hierarchy_layer_modes} are read.
#' @param config Full config list.
#' @param parent_tables Optional named list (by ontology) of parent tables; when
#'   absent each is read once from GO.db with go_parent_table().
#' @param term_info Optional go_term_info()-shaped table; when absent, from GO.db.
#' @return NULL when disabled, else list(settings, units): \code{units} is a
#'   list of one entry per layer x contrast x ontology x direction (plus one
#'   per layer that has no GO fgsea result) with \code{layer}, \code{contrast},
#'   \code{ontology}, \code{direction}, \code{direction_label},
#'   \code{n_significant}, \code{dag}, \code{plot} and \code{note}.
build_go_hierarchy_dags <- function(enrichment_results, config,
                                    parent_tables = NULL, term_info = NULL) {
    gh <- go_hierarchy_config(config)
    if (!gh$enabled) return(NULL)
    if (is.null(parent_tables)) parent_tables <- list()
    units <- list()
    add_unit <- function(...) units[[length(units) + 1]] <<- list(...)

    for (layer in names(.go_hierarchy_layer_modes)) {
        res <- enrichment_results[[layer]]
        if (is.null(res)) next
        tbl <- collect_go_fgsea_tables(res)
        tbl <- tbl[tbl$ontology %in% gh$ontologies, , drop = FALSE]
        if (nrow(tbl) == 0) {
            add_unit(layer = layer, contrast = NA_character_, ontology = NA_character_,
                     direction = NA_character_, direction_label = NA_character_,
                     n_significant = 0L, dag = NULL, plot = NULL,
                     note = paste("No fgsea result on GO",
                                  paste(gh$ontologies, collapse = "/"),
                                  "gene sets for this layer."))
            next
        }
        for (cn in unique(tbl$contrast)) {
            for (ont in intersect(gh$ontologies, unique(tbl$ontology[tbl$contrast == cn]))) {
                sub_tbl <- tbl[tbl$contrast == cn & tbl$ontology == ont, , drop = FALSE]
                seeds <- select_go_dag_seeds(sub_tbl, gh$fdr, gh$max_seed_terms)
                for (dirn in c("positive", "negative")) {
                    s <- seeds[[dirn]]
                    label <- go_direction_label(cn, dirn)
                    if (nrow(s$seeds) == 0) {
                        add_unit(layer = layer, contrast = cn, ontology = ont,
                                 direction = dirn, direction_label = label,
                                 n_significant = 0L, dag = NULL, plot = NULL,
                                 note = sprintf("No GO %s term with padj < %s and %s NES.",
                                                ont, format(gh$fdr), dirn))
                        next
                    }
                    if (is.null(parent_tables[[ont]])) parent_tables[[ont]] <- go_parent_table(ont)
                    dag <- build_go_dag(s$seeds, parent_tables[[ont]],
                                        term_info = term_info, max_nodes = gh$max_nodes)
                    note <- go_hierarchy_note(s$n_significant, dag$n_seeds_requested,
                                              dag$n_seeds_shown, gh$max_nodes)
                    plot <- plot_go_dag(
                        dag, fdr = gh$fdr,
                        title = sprintf("GO %s terms, %s, %s (%s)", ont,
                                        gsub("_", " ", cn), label, layer),
                        subtitle = paste0("fgsea padj < ", format(gh$fdr),
                                          ". ", note))
                    add_unit(layer = layer, contrast = cn, ontology = ont,
                             direction = dirn, direction_label = label,
                             n_significant = s$n_significant, dag = dag,
                             plot = plot, note = note)
                }
            }
        }
    }
    list(settings = gh, units = units)
}

#' One sentence saying how many significant terms a figure shows
#'
#' @param n_significant Terms passing the cutoff in this direction.
#' @param n_requested Seeds after the max_seed_terms cap.
#' @param n_shown Seeds drawn after the node limit.
#' @param max_nodes The node limit.
#' @return Character scalar.
go_hierarchy_note <- function(n_significant, n_requested, n_shown, max_nodes) {
    txt <- sprintf("%d of %d significant terms shown as seeds", n_shown, n_significant)
    if (n_requested < n_significant) {
        txt <- paste0(txt, " (top ", n_requested, " by padj")
        txt <- paste0(txt, if (n_shown < n_requested)
            sprintf(", then reduced to keep the graph within %d nodes)", max_nodes)
            else ")")
    } else if (n_shown < n_requested) {
        txt <- paste0(txt, sprintf(" (reduced by padj to keep the graph within %d nodes)",
                                   max_nodes))
    }
    paste0(txt, ".")
}

#' Figure legend for one GO hierarchy figure
#'
#' Kept here rather than in the template so its worst-case length is tested
#' (test-go-hierarchy.R) under the report's 100-word legend limit.
#'
#' @param ontology Ontology code, e.g. "BP".
#' @param fdr Adjusted p cutoff for seed terms.
#' @param direction "positive" or "negative".
#' @param direction_label From go_direction_label().
#' @param note From go_hierarchy_note().
#' @return Character scalar.
go_hierarchy_legend <- function(ontology, fdr, direction, direction_label, note) {
    paste(
        sprintf("GO %s terms with fgsea adjusted p < %s and %s NES",
                ontology, format(fdr), direction),
        sprintf("(%s), with their ancestors.", direction_label),
        "Filled nodes: seed terms, shaded by adjusted p (darker is",
        "smaller). Hollow grey nodes: context terms, not significant",
        "or not tested. Solid lines: is_a; dashed: part_of. More",
        "general terms are higher.", note)
}

#' Remove GO hierarchy outputs of an earlier run
#'
#' The report finds the figures through the index, so a run that draws nothing
#' must not leave the previous one's files to be read as current.
#'
#' @param out_dir The go_hierarchy output directory.
#' @param run_dir The run directory (holding the mode folders).
#' @return Invisibly, the paths removed.
clear_go_hierarchy_outputs <- function(out_dir, run_dir) {
    stale <- character(0)
    if (!is.null(out_dir) && dir.exists(out_dir)) stale <- out_dir
    if (!is.null(run_dir)) {
        tsv <- file.path(run_dir, unname(.go_hierarchy_layer_modes), "Datasets",
                         "go_dag_edges.tsv")
        stale <- c(stale, tsv[file.exists(tsv)])
    }
    if (length(stale) > 0) unlink(stale, recursive = TRUE)
    invisible(stale)
}

#' Write the GO hierarchy figures, tables and index
#'
#' Per unit with a graph: <ONT>_<direction>.png and .pdf, _nodes.csv and
#' _edges.csv under go_hierarchy/<layer>/<contrast>/. One index,
#' go_hierarchy_index.csv, lists every unit (with or without a figure), with
#' paths relative to \code{out_dir}. Per layer, the union of drawn edges goes
#' to <run_dir>/<mode>/Datasets/go_dag_edges.tsv.
#'
#' @param built Output of build_go_hierarchy_dags(), or NULL.
#' @param out_dir The go_hierarchy output directory.
#' @param run_dir The run directory.
#' @param term_info Optional go_term_info()-shaped table for the edge export's
#'   parent names; when absent, from GO.db.
#' @return Character vector of the files written (empty when \code{built} is
#'   NULL).
write_go_hierarchy_outputs <- function(built, out_dir, run_dir, term_info = NULL) {
    clear_go_hierarchy_outputs(out_dir, run_dir)
    if (is.null(built)) return(character(0))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    safe <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
    written <- character(0)
    index_rows <- list()
    edges_by_layer <- list()

    units <- built$units
    if (!is.null(built$error)) {
        units <- c(units, list(list(
            layer = NA_character_, contrast = NA_character_,
            ontology = NA_character_, direction = NA_character_,
            direction_label = NA_character_, n_significant = 0L, dag = NULL,
            note = paste("Not drawn:", built$error))))
    }
    for (u in units) {
        files <- c(png = NA_character_, pdf = NA_character_,
                   nodes_csv = NA_character_, edges_csv = NA_character_)
        n_nodes <- 0L; n_edges <- 0L; n_shown <- 0L
        if (!is.null(u$dag)) {
            rel_dir <- file.path(safe(u$layer), safe(u$contrast))
            dir.create(file.path(out_dir, rel_dir), recursive = TRUE, showWarnings = FALSE)
            stem <- file.path(rel_dir, paste0(u$ontology, "_", u$direction))
            files <- c(png = paste0(stem, ".png"), pdf = paste0(stem, ".pdf"),
                       nodes_csv = paste0(stem, "_nodes.csv"),
                       edges_csv = paste0(stem, "_edges.csv"))
            nodes <- u$dag$nodes
            edges <- u$dag$edges
            n_nodes <- nrow(nodes); n_edges <- nrow(edges); n_shown <- u$dag$n_seeds_shown
            sz <- go_dag_figure_size(n_nodes)
            ggplot2::ggsave(file.path(out_dir, files[["png"]]), u$plot,
                            width = sz[["width"]], height = sz[["height"]], dpi = 150)
            ggplot2::ggsave(file.path(out_dir, files[["pdf"]]), u$plot,
                            width = sz[["width"]], height = sz[["height"]])
            utils::write.csv(nodes, file.path(out_dir, files[["nodes_csv"]]), row.names = FALSE)
            utils::write.csv(edges, file.path(out_dir, files[["edges_csv"]]), row.names = FALSE)
            written <- c(written, file.path(out_dir, files))
            edges_by_layer[[u$layer]] <- rbind(edges_by_layer[[u$layer]], edges)
        }
        index_rows[[length(index_rows) + 1]] <- data.frame(
            layer = u$layer, contrast = u$contrast, ontology = u$ontology,
            direction = u$direction, direction_label = u$direction_label,
            n_significant = u$n_significant, n_seeds_shown = n_shown,
            n_nodes = n_nodes, n_edges = n_edges, note = u$note,
            png = files[["png"]], pdf = files[["pdf"]],
            nodes_csv = files[["nodes_csv"]], edges_csv = files[["edges_csv"]],
            stringsAsFactors = FALSE)
    }

    index <- if (length(index_rows) > 0) do.call(rbind, index_rows) else
        data.frame(layer = character(0), contrast = character(0),
                   ontology = character(0), direction = character(0),
                   direction_label = character(0), n_significant = integer(0),
                   n_seeds_shown = integer(0), n_nodes = integer(0),
                   n_edges = integer(0), note = character(0), png = character(0),
                   pdf = character(0), nodes_csv = character(0),
                   edges_csv = character(0))
    index_file <- file.path(out_dir, "go_hierarchy_index.csv")
    utils::write.csv(index, index_file, row.names = FALSE)
    written <- c(written, index_file)

    for (layer in names(edges_by_layer)) {
        tsv <- go_dag_edges_export(edges_by_layer[[layer]], term_info = term_info)
        ds_dir <- file.path(run_dir, .go_hierarchy_layer_modes[[layer]], "Datasets")
        dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
        tsv_file <- file.path(ds_dir, "go_dag_edges.tsv")
        utils::write.table(tsv, tsv_file, sep = "\t", quote = FALSE,
                           row.names = FALSE, na = "")
        written <- c(written, tsv_file)
    }
    written
}

#' The portal's GO edge table
#'
#' @param edges Row-bound edge tables (\code{child}, \code{parent},
#'   \code{relation}) of one layer's figures.
#' @param term_info Optional go_term_info()-shaped table for the parents; when
#'   absent, from GO.db.
#' @return data.frame with \code{child}, \code{parent}, \code{relation},
#'   \code{parent_name}, \code{parent_ontology}, one row per distinct link,
#'   ordered by child then parent.
go_dag_edges_export <- function(edges, term_info = NULL) {
    edges <- unique(edges[, c("child", "parent", "relation"), drop = FALSE])
    edges <- edges[order(edges$child, edges$parent, edges$relation), , drop = FALSE]
    if (is.null(term_info)) term_info <- go_term_info(edges$parent)
    idx <- match(edges$parent, term_info$go_id)
    edges$parent_name <- term_info$term[idx]
    edges$parent_ontology <- term_info$ontology[idx]
    rownames(edges) <- NULL
    edges
}

#' Read the GO hierarchy index for the report
#'
#' @param out_dir The go_hierarchy output directory.
#' @return The index data.frame with absolute figure paths, or NULL when the
#'   index is absent.
read_go_hierarchy_index <- function(out_dir) {
    f <- file.path(out_dir, "go_hierarchy_index.csv")
    if (!file.exists(f)) return(NULL)
    idx <- utils::read.csv(f, stringsAsFactors = FALSE)
    for (col in c("png", "pdf", "nodes_csv", "edges_csv")) {
        has <- !is.na(idx[[col]]) & nzchar(idx[[col]])
        idx[[col]][has] <- file.path(out_dir, idx[[col]][has])
    }
    idx
}
