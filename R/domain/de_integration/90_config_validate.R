# R/domain/de_integration/90_config_validate.R
#
# Config contract for the DE-integration mode: layers of finished DE tables
# (proteomics or RNA-seq, the pipeline's own exports or any table with a column
# map), compared without the sample-level matrices the multiomics mode needs.

#' Formats a DE-integration layer can be read from
#'
#' @return Character vector of the accepted \code{format} values.
de_integration_formats <- function() {
    c("proteomics_summary", "proteomics_final", "rnaseq_summary", "rnaseq_final",
      "generic")
}


#' Validate the DE-integration config and fill its defaults
#'
#' Checks every layer and comparison and fills the defaults the readers rely
#' on, so a mistake is reported once, at config time, naming the key to fix --
#' not as an empty table three targets later.
#'
#' @param cfg The \code{modes$de_integration} section of the config.
#' @return The same section with defaults filled in.
#' @examples
#' cfg <- list(layers = list(
#'     list(name = "cells", omics_type = "proteomics",
#'          format = "proteomics_summary", path = "cells_summary.tsv"),
#'     list(name = "media", omics_type = "proteomics",
#'          format = "proteomics_summary", path = "media_summary.tsv")))
#' validate_de_integration_config(cfg)$hits$p_cutoff   # 0.05
validate_de_integration_config <- function(cfg) {
    where <- "modes.de_integration"
    layers <- cfg$layers
    if (!is.list(layers) || length(layers) < 2) {
        stop(where, ".layers needs at least two layers to compare; found ",
             length(layers), ".", call. = FALSE)
    }

    formats <- de_integration_formats()
    layer_names <- character(0)
    for (i in seq_along(layers)) {
        ly <- layers[[i]]
        at <- sprintf("%s.layers[[%d]]", where, i)
        nm <- ly$name
        if (!is.character(nm) || length(nm) != 1 || !grepl("^[a-z][a-z0-9_]*$", nm)) {
            stop(at, ".name must be one lower-case word (letters, digits, '_', ",
                 "starting with a letter); it becomes a column prefix. Got: ",
                 paste(format(nm), collapse = " "), call. = FALSE)
        }
        if (nm %in% layer_names) {
            stop(at, ".name '", nm, "' is used by another layer; names must be unique.",
                 call. = FALSE)
        }
        layer_names <- c(layer_names, nm)

        if (!identical(length(ly$omics_type), 1L) ||
            !ly$omics_type %in% c("proteomics", "rnaseq")) {
            stop(at, " ('", nm, "').omics_type must be \"proteomics\" or \"rnaseq\".",
                 call. = FALSE)
        }
        if (!identical(length(ly$format), 1L) || !ly$format %in% formats) {
            stop(at, " ('", nm, "').format must be one of: ",
                 paste(formats, collapse = ", "), ".", call. = FALSE)
        }
        if (!startsWith(ly$format, "generic") &&
            !startsWith(ly$format, ly$omics_type)) {
            stop(at, " ('", nm, "'): format \"", ly$format, "\" is not a ",
                 ly$omics_type, " export. Use a ", ly$omics_type,
                 "_* format or \"generic\" with a column map.", call. = FALSE)
        }
        if (!is.character(ly$path) || length(ly$path) != 1 || !nzchar(ly$path)) {
            stop(at, " ('", nm, "').path must be the path to the DE table.",
                 call. = FALSE)
        }
        if (identical(ly$format, "generic")) {
            cols <- ly$columns %||% list()
            missing <- c("id", "pvalue")[!c("id", "pvalue") %in% names(cols)]
            if (is.null(cols$log2fc) && is.null(cols$linear_fc)) {
                missing <- c(missing, "log2fc (or linear_fc)")
            }
            if (length(missing) > 0) {
                stop(at, " ('", nm, "') is \"generic\" and its columns map lacks: ",
                     paste(missing, collapse = ", "), ".", call. = FALSE)
            }
        }
        if (!is.null(ly$observed)) {
            obs <- ly$observed
            need <- c("matrix", "samplesheet", "sample_col", "contrasts_file")
            gap <- need[!need %in% names(obs)]
            if (length(gap) > 0) {
                stop(at, " ('", nm, "').observed needs ", paste(gap, collapse = ", "),
                     " to count observed values per group.", call. = FALSE)
            }
        }
        layers[[i]]$label <- ly$label %||% nm
    }
    cfg$layers <- layers

    comps <- cfg$comparisons %||% list()
    comp_names <- character(0)
    for (j in seq_along(comps)) {
        cp <- comps[[j]]
        at <- sprintf("%s.comparisons[[%d]]", where, j)
        if (!is.character(cp$name) || length(cp$name) != 1 || !nzchar(cp$name)) {
            stop(at, ".name is required.", call. = FALSE)
        }
        if (cp$name %in% comp_names) {
            stop(at, ".name '", cp$name, "' is used twice.", call. = FALSE)
        }
        comp_names <- c(comp_names, cp$name)
        members <- cp$members
        if (!is.list(members) || length(members) < 2) {
            stop(at, " ('", cp$name, "').members must name a contrast for at least ",
                 "two layers, e.g. {cells: A_vs_B, media: A_vs_B}.", call. = FALSE)
        }
        unknown <- setdiff(names(members), layer_names)
        if (length(unknown) > 0) {
            stop(at, " ('", cp$name, "').members names unknown layer(s): ",
                 paste(unknown, collapse = ", "), ". Layers: ",
                 paste(layer_names, collapse = ", "), ".", call. = FALSE)
        }
        flip <- unlist(cp$flip %||% character(0), use.names = FALSE)
        bad_flip <- setdiff(flip, names(members))
        if (length(bad_flip) > 0) {
            stop(at, " ('", cp$name, "').flip names layer(s) not in its members: ",
                 paste(bad_flip, collapse = ", "), ".", call. = FALSE)
        }
        comps[[j]]$flip <- as.character(flip)
    }
    cfg$comparisons <- comps

    hits <- cfg$hits %||% list()
    hits$use_table_flag   <- hits$use_table_flag %||% TRUE
    hits$p_cutoff         <- hits$p_cutoff %||% 0.05
    hits$use_adjusted     <- hits$use_adjusted %||% TRUE
    hits$linear_fc_cutoff <- hits$linear_fc_cutoff %||% 1.5
    if (!is.numeric(hits$p_cutoff) || hits$p_cutoff <= 0 || hits$p_cutoff > 1) {
        stop(where, ".hits.p_cutoff must be a number in (0, 1].", call. = FALSE)
    }
    if (!is.numeric(hits$linear_fc_cutoff) || hits$linear_fc_cutoff < 1) {
        stop(where, ".hits.linear_fc_cutoff is a linear fold change and must be ",
             ">= 1 (1.5 means 1.5-fold either way).", call. = FALSE)
    }
    cfg$hits <- hits

    conc <- cfg$concordance %||% list()
    conc$well_observed_min <- conc$well_observed_min %||% 2
    cfg$concordance <- conc

    cfg
}
