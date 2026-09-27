# R/domain/de_integration/00_inputs.R
#
# What a DE-integration layer holds and which parts of it a run uses: the
# column names of each contrast, the contrasts a table carries, how requested
# contrasts resolve and pair across layers, and the input files to track. A
# layer is a finished DE table -- one of this pipeline's own exports, or any
# table with a column map. 01_read_layer.R builds the standardized frames every
# later step works on.

#' Column names for one contrast of a DE-integration layer
#'
#' The pipeline's own exports are described by \code{get_contrast_cols()}, the
#' shared naming contract, so this reads exactly the columns the producers
#' write. A generic table names its columns in the layer's \code{columns} map;
#' that map is also where an external table's own spelling of its observed
#' counts belongs (\code{n_obs_num}, \code{n_obs_den}).
#'
#' @param format One of \code{de_integration_formats()}.
#' @param contrast Contrast label as it appears in the column suffixes.
#' @param columns The layer's \code{columns} map (generic format only).
#' @param id_col The layer's \code{id_col} override (native formats only):
#'   the feature id column, for a run that exported its ids under another
#'   name. NULL keeps the format's own id column.
#' @return List with \code{id}, \code{symbol}, \code{gene_id},
#'   \code{description}, \code{log2fc} (candidates, in preference order),
#'   \code{linear_ratio}, \code{linear_fc}, \code{pvalue}, \code{padj},
#'   \code{hit}, \code{hit_kind} ("pass" or "updown") and
#'   \code{n_obs_prefix} (the prefix of the per-group observed-count columns);
#'   NULL for anything the format does not carry.
#' @examples
#' de_layer_columns("proteomics_summary", "A_vs_B")$pvalue   # "pvalue.imputs.A_vs_B"
#' de_layer_columns("proteomics_final", "A_vs_B", id_col = "Protein.Group")$id
de_layer_columns <- function(format, contrast, columns = NULL, id_col = NULL) {
    if (identical(format, "generic")) {
        cols <- columns %||% list()
        return(list(
            id = cols$id, symbol = cols$symbol, gene_id = cols$gene_id,
            description = cols$description,
            log2fc = cols$log2fc, linear_ratio = NULL, linear_fc = cols$linear_fc,
            pvalue = cols$pvalue, padj = cols$padj,
            hit = cols$hit, hit_kind = "pass",
            n_obs_prefix = NULL,
            n_obs_num = cols$n_obs_num, n_obs_den = cols$n_obs_den
        ))
    }

    omics <- sub("_.*$", "", format)
    final <- endsWith(format, "_final")
    if (identical(omics, "proteomics")) {
        cc <- get_contrast_cols(contrast, mode = "proteomics")
        list(
            id = id_col %||% "FeatureID", symbol = "Genes", gene_id = NULL,
            description = "First.Protein.Description",
            log2fc = cc$log2fc,
            # Spelled like the other columns: the proteomics exports strip spaces
            # from contrast names (normalize_contrast_name()).
            linear_ratio = paste0("linearRatio.imputs.", normalize_contrast_name(contrast)),
            linear_fc = cc$fc, pvalue = cc$p, padj = cc$padj,
            hit = if (final) cc$updown else cc$pass,
            hit_kind = if (final) "updown" else "pass",
            # Written by build_group_raw_stats_proteomics() when the group
            # summaries are exported (modes.proteomics.excel.group_cv).
            n_obs_prefix = "N.observed."
        )
    } else {
        cc <- get_contrast_cols(contrast, mode = "rna")
        list(
            # The RNA exports rename FeatureID to Gene (05_outputs_legacy.R).
            id = id_col %||% c("Gene", "FeatureID"), symbol = NULL,
            gene_id = id_col %||% c("Gene", "FeatureID"),
            description = NULL,
            # log2FoldChange.<c> is what exports written before log2FC.<c> carry.
            log2fc = c(cc$log2fc, paste0("log2FoldChange.", contrast)),
            linear_ratio = NULL, linear_fc = cc$fc, pvalue = cc$p, padj = cc$padj,
            hit = if (final) cc$updown else cc$pass,
            hit_kind = if (final) "updown" else "pass",
            n_obs_prefix = NULL
        )
    }
}


#' Contrasts a DE-integration layer's table holds
#'
#' Read off the p-value columns, through \code{.de_summary_candidates()}, the
#' same helper that decides which columns a pre-computed summary resolves to.
#'
#' @param cn Column names of the layer's table.
#' @param format The layer's format.
#' @param generic_contrast For a generic table: the label of its one contrast.
#' @return Character vector of contrast labels.
list_layer_contrasts <- function(cn, format, generic_contrast = "contrast") {
    if (identical(format, "generic")) return(generic_contrast)
    prefix <- if (startsWith(format, "proteomics")) "pvalue.imputs" else "pvalue"
    unique(.de_summary_candidates(cn, prefix)$contrast)
}


#' Pick the table's contrast that a requested label names
#'
#' The exact label, else the one contrast whose \code{normalize_contrast_key()}
#' matches, so "A_vs_B", "A vs. B" and "A - B" agree. Nothing else: a label
#' someone wrote into \code{comparisons} that the table does not hold is an
#' error, never quietly swapped for whatever the table does hold. (Inferred
#' comparisons name the tables' own labels, so they always resolve here; the
#' pairing of two differently named lone contrasts is decided, and announced,
#' in \code{resolve_dei_comparisons()}.)
#'
#' @param requested Contrast label asked for.
#' @param available Contrast labels the table holds.
#' @param layer Layer name, for the error message.
#' @return The matching label from \code{available}.
resolve_layer_contrast <- function(requested, available, layer) {
    if (requested %in% available) return(requested)
    by_key <- available[normalize_contrast_key(available) ==
                        normalize_contrast_key(requested)]
    if (length(by_key) == 1) return(by_key)
    why <- if (length(by_key) > 1) {
        paste0("several contrasts match it by key (", paste(by_key, collapse = ", "), ")")
    } else {
        "no contrast matches it"
    }
    stop("Layer '", layer, "': contrast '", requested, "' cannot be resolved: ", why,
         ". The table holds: ", paste(available, collapse = ", "),
         ". Name one of these in modes.de_integration.comparisons.", call. = FALSE)
}


#' Comparisons to run, each naming one contrast per layer
#'
#' Configured comparisons are taken as written (their contrasts are resolved
#' against each table later). Without any, contrasts are paired across layers
#' by \code{normalize_contrast_key()}: every key two or more layers share
#' becomes a comparison. When every layer holds exactly one contrast and they
#' share no key, those single contrasts are paired, with a message -- two runs
#' of one experiment often spell their contrast differently.
#'
#' @param cfg The validated \code{modes$de_integration} section.
#' @param layer_contrasts Named list: layer name -> contrasts its table holds.
#' @return List of comparisons, each \code{list(name, members, flip)} with
#'   \code{members} a named character vector (layer -> contrast).
resolve_dei_comparisons <- function(cfg, layer_contrasts) {
    if (length(cfg$comparisons %||% list()) > 0) {
        return(lapply(cfg$comparisons, function(cp) list(
            name = cp$name,
            members = unlist(cp$members),
            flip = as.character(cp$flip %||% character(0)))))
    }

    keys <- lapply(layer_contrasts, normalize_contrast_key)
    # Two contrasts of one layer on one key cannot be told apart by name, and
    # leaving that layer out of the comparison would drop its data unsaid.
    for (ly in names(keys)) {
        dup <- unique(keys[[ly]][duplicated(keys[[ly]])])
        if (length(dup) > 0) {
            same <- layer_contrasts[[ly]][keys[[ly]] %in% dup]
            stop("Layer '", ly, "' holds contrasts that name the same comparison (",
                 paste(same, collapse = ", "), "), so they cannot be paired by name. ",
                 "Name the pairs in modes.de_integration.comparisons.", call. = FALSE)
        }
    }
    all_keys <- sort(unique(unlist(keys)))
    comps <- list()
    for (k in all_keys) {
        members <- unlist(lapply(names(layer_contrasts), function(ly) {
            hit <- layer_contrasts[[ly]][keys[[ly]] == k]
            if (length(hit) == 1) stats::setNames(hit, ly) else NULL
        }))
        if (length(members) >= 2) {
            # Named after the first layer's own spelling; the key only decided
            # which contrasts belong together.
            comps[[length(comps) + 1]] <- list(name = unname(members[1]),
                                               members = members,
                                               flip = character(0))
        }
    }
    if (length(comps) > 0) return(comps)

    if (all(lengths(layer_contrasts) == 1)) {
        members <- unlist(layer_contrasts)
        message("  No contrast name is shared across layers; pairing each layer's ",
                "only contrast: ", paste(names(members), members, sep = " = ",
                                         collapse = ", "), ".")
        return(list(list(name = "comparison_1", members = members,
                         flip = character(0))))
    }
    stop("No contrast is shared by two or more layers, and at least one layer ",
         "holds several. Name the pairs in modes.de_integration.comparisons.",
         call. = FALSE)
}


#' Every input file a DE-integration run reads
#'
#' Tracked by a file target, so editing any of them reruns what depends on it.
#' A missing file stops the run here, naming the layer and key.
#'
#' @param config Full config.
#' @return Character vector of resolved paths.
de_integration_input_files <- function(config) {
    cfg <- config$modes$de_integration
    paths <- character(0)
    for (ly in cfg$layers) {
        keys <- list(path = ly$path, annotation_file = ly$annotation_file,
                     observed.matrix = ly$observed$matrix,
                     observed.samplesheet = ly$observed$samplesheet,
                     observed.contrasts_file = ly$observed$contrasts_file)
        for (k in names(keys)) {
            if (is.null(keys[[k]])) next
            p <- resolve_input_path(config, keys[[k]])
            if (!file.exists(p)) {
                stop("Layer '", ly$name, "': ", k, " not found: ", p, call. = FALSE)
            }
            paths <- c(paths, p)
        }
    }
    unique(paths)
}
