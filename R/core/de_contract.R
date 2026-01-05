#' Extract DE feature IDs from a DE result object
#'
#' Omics-agnostic helper that returns feature identifiers passing
#' differential expression thresholds as defined in config.
#'
#' @param de_res DE result object containing a summary table
#' @param config Global config
#'
#' @return character vector of feature IDs
get_de_features <- function(de_res, config) {
  stopifnot(is.list(de_res), !is.null(de_res$summary_df))
  df <- de_res$summary_df
  stopifnot(is.data.frame(df))
  
  # Accept legacy + new column names
  id_candidates <- c("feature_id", "FeatureID", "Protein", "protein_id")
  id_col <- id_candidates[id_candidates %in% colnames(df)][1]
  
  if (is.na(id_col) || is.null(id_col)) {
    stop("DE summary missing id column. Expected one of: ",
         paste(id_candidates, collapse = ", "))
  }
  
  # Standardize internally
  df$feature_id <- df[[id_col]]
  
  # ---- your existing logic below, but use df$feature_id ----
  # Example: any_contrast = pass_any == 1
  if (!"pass_any_contrast" %in% colnames(df)) {
    stop("DE summary missing 'pass_any_contrast'. (expected from summarize_limma_mult_imputation)")
  }
  
  feats <- df$feature_id[!is.na(df$pass_any_contrast) & df$pass_any_contrast == 1]
  unique(as.character(feats))
}

#' Standardized per-contrast column names in the DE summary table
#'
#' This project uses a "wide" summary format where each contrast creates
#' multiple columns (fc/p/padj/pass/upDown/manual_cutoffs).
#'
#' @param contrast Character scalar. Contrast name (e.g., "T_vs_C").
#' @return Named list with standardized column names.
get_contrast_cols <- function(contrast) {
  stopifnot(is.character(contrast), length(contrast) == 1, nzchar(contrast))
  
  list(
    fc     = paste0("linearFC.imputs.", contrast),
    p      = paste0("pvalue.imputs.",   contrast),
    padj   = paste0("padj.imputs.",     contrast),
    pass   = paste0("pass.imputs.",     contrast),
    updown = paste0("upDown.imputs.",   contrast),
    manual = paste0("manual_cutoffs.",  contrast)
  )
}


