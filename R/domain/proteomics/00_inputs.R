#' Load proteomics inputs
#' @param config list as returned by load_config()
#' @return list (protein, sample_map, meta, contrasts, engine, ...)
load_proteomics_inputs <- function(config) {
    inputs <- load_omics_inputs(config, mode = "proteomics")

    cfg <- config$modes$proteomics
    is_preprocessed <- identical(cfg$input$format, "preprocessed")

    if (is_preprocessed && is.null(inputs$protein) && !is.null(inputs$preprocessed_protein)) {
        message("[load_proteomics_inputs] Using preprocessed_protein as protein matrix")
        inputs$protein <- inputs$preprocessed_protein
        inputs$source_type <- "preprocessed"
    }

    inputs
}

#' Validate proteomics inputs
validate_proteomics_inputs <- function(inputs, cfg) {
    protein <- inputs$protein
    sample_map <- inputs$sample_map
    meta <- inputs$metadata

    id_cols <- cfg$id_columns
    eff_cols <- cfg$effects

    # For preprocessed inputs, only validate protein ID column exists
    if (identical(inputs$source_type, "preprocessed") ||
        identical(cfg$input$format, "preprocessed")) {
        if (!id_cols$protein_id %in% colnames(protein)) {
            stop(sprintf(
                "Preprocessed protein table missing ID column '%s'. Available: %s",
                id_cols$protein_id, paste(colnames(protein), collapse = ", ")
            ))
        }
        return(invisible(TRUE))
    }

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

# load_omics_inputs and validate_contrasts_content live in R/core/01_io.R
