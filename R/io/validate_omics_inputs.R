#' Check that required columns exist in a data frame
#'
#' @param df A data.frame / tibble
#' @param required Character vector of required column names
#' @param df_name Optional name of the data frame (for error messages)
check_has_cols <- function(df, required, df_name = deparse(substitute(df))) {
  required <- unlist(required)
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(
      "In table '", df_name, "': missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

#' Check that all values in x are present in y
#'
#' @param x Vector of values that must be contained in y
#' @param y Reference vector
#' @param label_x Name/label for x (used in error message)
#' @param label_y Name/label for y (used in error message)
check_all_in <- function(x, y, label_x = "x", label_y = "y") {
  missing <- setdiff(x, y)
  if (length(missing) > 0) {
    stop(
      "Values in ", label_x, " not found in ", label_y, ": ",
      paste(head(missing, 10), collapse = ", "),
      if (length(missing) > 10) " ... (truncated)" else "",
      call. = FALSE
    )
  }
}

#' Generic omics input validation
#'
#' This function loads inputs for a given omics mode (e.g. proteomics, rna)
#' using \code{load_omics_inputs()}, and then dispatches to a mode-specific
#' validation function.
#'
#' @param config List as returned by \code{load_config()}
#' @param mode   Character, e.g. "proteomics", "rna"
#'
#' @return (invisibly) the loaded inputs (a named list)
validate_omics_inputs <- function(config,
                                  mode = c("proteomics", "rna",
                                           "metabolomics", "lipidomics")) {
  mode <- match.arg(mode)
  
  # Load inputs using the generic loader
  inputs <- load_omics_inputs(config, mode = mode)
  cfg    <- config$modes[[mode]]
  
  # Dispatch to mode-specific validator
  switch(
    mode,
    proteomics  = validate_proteomics_inputs(inputs,  cfg),
    rna         = validate_rna_inputs(inputs,         cfg),
    metabolomics = validate_metabolomics_inputs(inputs, cfg),
    lipidomics   = validate_lipidomics_inputs(inputs,   cfg)
  )
  
  invisible(inputs)
}

#' Validate proteomics inputs (wide matrix: samples in columns)
#'
#' Assumes:
#' - protein: one row per protein, one column per sample (wide format)
#'            with non-sample columns defined in cfg$id_columns
#' - sample_map: maps from raw sample names (map_from) to unified sample IDs (map_to)
#' - metadata: one row per sample, with sample ID = effects$samples
#'
#' @param inputs List returned by load_proteomics_inputs()
#' @param cfg    config$modes$proteomics
#'
#' @return TRUE (invisibly) if validation passes; otherwise stops with an error
validate_proteomics_inputs <- function(inputs, cfg) {
  protein    <- inputs$protein
  sample_map <- inputs$sample_map
  meta       <- inputs$metadata
  contrasts  <- inputs$contrasts
  
  id_cols  <- cfg$id_columns
  eff_cols <- cfg$effects
  
  # 1) Column existence checks ----------------------------------
  
  # protein: must have the protein ID column
  check_has_cols(
    protein,
    id_cols$protein_id,
    df_name = "protein"
  )
  
  # sample_map: must have mapping columns
  check_has_cols(
    sample_map,
    c(id_cols$map_from, id_cols$map_to),
    df_name = "sample_map"
  )
  
  # metadata: must have sample label column + effect columns (color/shape/label)
  meta_required <- unique(c(
    eff_cols$samples,
    eff_cols$color,
    eff_cols$shape
  ))
  meta_required <- meta_required[!is.na(meta_required)]
  
  check_has_cols(
    meta,
    meta_required,
    df_name = "metadata"
  )
  
  # 2) Determine sample IDs from protein matrix -----------------
  
  # Non-sample columns: protein_id + optional annotation columns
  annot_cols <- id_cols$protein_annot
  if (is.null(annot_cols)) {
    annot_cols <- character(0)
  } else {
    annot_cols <- unlist(annot_cols)
  }
  
  non_sample_cols <- unique(c(id_cols$protein_id, annot_cols))
  
  # Wide format: all remaining columns are assumed to be samples
  protein_sample_cols <- setdiff(colnames(protein), non_sample_cols)
  
  # 3) Cross-table consistency checks ---------------------------
  
  # Here we assume:
  # - protein column names correspond to the *raw* sample names
  # - those raw names are stored in sample_map[[map_from]]
  #
  # So: all protein sample columns must exist in sample_map map_from
  check_all_in(
    x = protein_sample_cols,
    y = sample_map[[id_cols$map_from]],
    label_x = "protein sample columns",
    label_y = "sample_map map_from"
  )
  
  # And all unified sample IDs (map_to) must exist in metadata samples column
  meta_sample_ids <- meta[[eff_cols$samples]]
  
  check_all_in(
    x = sample_map[[id_cols$map_to]],
    y = meta_sample_ids,
    label_x = "sample_map map_to",
    label_y = "metadata label"
  )
  
  invisible(TRUE)
}


#' Validate RNA inputs
#'
#' Placeholder / skeleton. To be filled once RNA files and id_columns
#' are fully defined in the config.
#'
#' @param inputs List returned by \code{load_rna_inputs()}
#' @param cfg    The \code{config$modes$rna} sub-list
validate_rna_inputs <- function(inputs, cfg) {
  # Example structure; adjust once your RNA inputs format is finalized
  # counts     <- inputs$counts
  # sample_map <- inputs$sample_map
  # meta       <- inputs$metadata
  # contrasts  <- inputs$contrasts
  #
  # id_cols  <- cfg$id_columns
  # eff_cols <- cfg$effects
  #
  # check_has_cols(counts, id_cols$gene_id, df_name = "counts")
  #
  # (Add consistency checks similar to proteomics)
  
  invisible(TRUE)
}

#' Validate metabolomics inputs
#'
#' Placeholder function. Implement when metabolomics file structure is defined.
validate_metabolomics_inputs <- function(inputs, cfg) {
  # Implement when metabolomics inputs are added to the config
  invisible(TRUE)
}

#' Validate lipidomics inputs
#'
#' Placeholder function. Implement when lipidomics file structure is defined.
validate_lipidomics_inputs <- function(inputs, cfg) {
  # Implement when lipidomics inputs are added to the config
  invisible(TRUE)
}
