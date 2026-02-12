#' Load proteomics inputs
#' @param config list as returned by load_config()
#' @return list (protein, sample_map, meta, contrasts, engine, ...)
load_proteomics_inputs <- function(config) {
    load_omics_inputs(config, mode = "proteomics")
}

#' Validate proteomics inputs
validate_proteomics_inputs <- function(inputs, cfg) {
    protein <- inputs$protein
    sample_map <- inputs$sample_map
    meta <- inputs$metadata

    id_cols <- cfg$id_columns
    eff_cols <- cfg$effects

    check_has_cols(protein, id_cols$protein_id, df_name = "protein")

    if (!is.null(sample_map)) {
        # Validate sample map columns and consistency
        check_has_cols(sample_map, c(id_cols$map_from, id_cols$map_to), df_name = "sample_map")

        annot_cols <- unlist(id_cols$protein_annot %||% character(0))
        non_sample_cols <- unique(c(id_cols$protein_id, annot_cols))
        protein_sample_cols <- setdiff(colnames(protein), non_sample_cols)

        check_all_in(protein_sample_cols, sample_map[[id_cols$map_from]], "protein cols", "sample_map raw")
        check_all_in(sample_map[[id_cols$map_to]], meta[[eff_cols$samples]], "sample_map unified", "meta samples")
    } else {
        # No sample map: protein column names must already match metadata sample IDs
        annot_cols <- unlist(id_cols$protein_annot %||% character(0))
        non_sample_cols <- unique(c(id_cols$protein_id, annot_cols))
        protein_sample_cols <- setdiff(colnames(protein), non_sample_cols)

        meta_samples <- meta[[eff_cols$samples]]
        matched <- intersect(protein_sample_cols, meta_samples)
        if (length(matched) == 0) {
            stop("No sample map provided and protein column names do not match metadata sample IDs. ",
                 "Provide a sample_map file or ensure column names match.")
        }
        message(sprintf("  No sample map: %d/%d protein columns matched metadata sample IDs directly.",
                        length(matched), length(protein_sample_cols)))
    }

    invisible(TRUE)
}

#' Generic loader helper (could be in core/01_io or locally here if only used by omics loaders)
load_omics_inputs <- function(config, mode = c("proteomics", "rna")) {
    mode <- match.arg(mode)
    cfg <- config$modes[[mode]]
    if (is.null(cfg)) stop("No config for mode ", mode)

    files <- cfg$files
    inputs <- list()

    for (nm in names(files)) {
        rel <- files[[nm]]
        if (is.null(rel) || !nzchar(rel)) next
        abs <- resolve_raw_path(config, rel)
        if (!file.exists(abs)) stop("File not found: ", abs)
        inputs[[nm]] <- read_table_auto(abs)
    }

    if (!is.null(cfg$engine)) inputs$engine <- cfg$engine
    inputs
}
