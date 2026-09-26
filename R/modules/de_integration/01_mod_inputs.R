# R/modules/de_integration/01_mod_inputs.R

#' Load every DE-integration layer for the comparisons it takes part in
#'
#' Each table is read once. Its contrasts are listed first, so the comparisons
#' can be resolved -- configured, or paired by contrast name -- and then each
#' layer is standardized for exactly the contrasts those comparisons ask of it.
#' A layer no comparison uses is reported and left out.
#'
#' @param config Full, validated config.
#' @return List with \code{layers} (named by layer; see \code{read_de_layer()}),
#'   \code{comparisons} (see \code{resolve_dei_comparisons()}) and
#'   \code{summary} (one row per layer and contrast).
mod_dei_load_layers <- function(config) {
    cfg <- config$modes$de_integration
    names(cfg$layers) <- vapply(cfg$layers, function(ly) ly$name, character(1))

    raw <- lapply(cfg$layers, function(ly) {
        read_table_auto(resolve_input_path(config, ly$path))
    })
    layer_contrasts <- lapply(cfg$layers, function(ly) {
        list_layer_contrasts(names(raw[[ly$name]]), ly$format, ly$contrast %||% "contrast")
    })
    comparisons <- resolve_dei_comparisons(cfg, layer_contrasts)

    layers <- list()
    for (ly in cfg$layers) {
        wanted <- unlist(lapply(comparisons, function(cp) unname(cp$members[ly$name])))
        wanted <- wanted[!is.na(wanted)]
        if (length(wanted) == 0) {
            message("  Layer '", ly$name, "' is in no comparison; left out.")
            next
        }
        message("  Reading layer '", ly$name, "' (", ly$format, "): ",
                paste(unique(wanted), collapse = ", "))
        layers[[ly$name]] <- read_de_layer(
            ly, wanted, config, hits_default = cfg$hits,
            well_observed_min = cfg$concordance$well_observed_min,
            df = raw[[ly$name]])
    }

    list(layers = layers, comparisons = comparisons,
         summary = summarize_de_layers(layers))
}


#' One row per layer and contrast: what was read and how
#'
#' @param layers Named list of layers from \code{read_de_layer()}.
#' @return Data frame of the layers' provenance rows, or an empty frame.
summarize_de_layers <- function(layers) {
    rows <- lapply(layers, function(ly) ly$provenance)
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) return(data.frame())
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}


#' Write the layer summary and the resolved comparisons
#'
#' @param dei_layers Output of \code{mod_dei_load_layers()}.
#' @param out_dir The mode's output directory.
#' @return Paths of the files written.
write_dei_layer_summary <- function(dei_layers, out_dir) {
    dir <- file.path(out_dir, "layers")
    comps <- do.call(rbind, lapply(dei_layers$comparisons, function(cp) data.frame(
        comparison = cp$name, layer = names(cp$members),
        contrast = unname(cp$members), flipped = names(cp$members) %in% cp$flip,
        stringsAsFactors = FALSE)))
    c(save_tsv(dei_layers$summary, dir, "layer_summary.tsv"),
      save_tsv(comps, dir, "comparisons.tsv"))
}
