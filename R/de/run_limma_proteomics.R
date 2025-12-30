#' Run limma differential analysis for proteomics with contrast support
#'
#' Fits a limma linear model on an imputed proteomics expression matrix and
#' returns per-contrast result tables with feature annotations.
#'
#' Key behaviors:
#' - aligns metadata rows to expression matrix columns by sample ID
#' - builds a ~0 + group design matrix (default group column: "Condition")
#' - applies user-defined contrasts (Numerator - Denominator) from `contrasts_df`
#' - aligns feature annotations from `prot_tbl` to matrix rownames using a configurable ID column
#' - reorders `topTable()` output to match expression row order to prevent silent misalignment
#'
#' @param expr_imp Numeric matrix (features × samples) of imputed expression values.
#' @param meta data.frame with one row per sample.
#' @param contrasts_df data.frame with columns: Contrast_name, Factor, Numerator, Denominator.
#' @param prot_tbl data.frame with feature annotations (must include the feature ID column).
#' @param cfg Full config list (expects `modes$proteomics$effects` and `modes$proteomics$id_columns`).
#'
#' @return A list with aligned metadata, design matrix, contrasts, fitted model and per-contrast DE tables.
#' @export
run_limma_proteomics <- function(expr_imp, meta, contrasts_df, prot_tbl, cfg) {
  p_cfg <- cfg$modes$proteomics
  
  sample_col <- p_cfg$effects$samples %||% "SampleID"
  group_col  <- p_cfg$effects$color   %||% "Condition"
  
  protein_id_col <- p_cfg$id_columns$protein_id %||% "Protein.Group"
  default_annot <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description")
  annot_cols <- unique(c(protein_id_col, p_cfg$id_columns$protein_annot %||% default_annot))
  
  assert_numeric_matrix(expr_imp, "expr_imp")
  meta_aligned <- align_meta_to_expr(expr_imp, meta, sample_col)
  
  meta_aligned[[group_col]] <- factor(meta_aligned[[group_col]])
  design <- model.matrix(stats::as.formula(paste0("~ 0 + ", group_col)), data = meta_aligned)
  colnames(design) <- levels(meta_aligned[[group_col]])
  
  stopifnot(all(contrasts_df$Factor == group_col))
  contrast_formulas <- setNames(
    paste(contrasts_df$Numerator, contrasts_df$Denominator, sep = " - "),
    contrasts_df$Contrast_name
  )

  contrast_matrix <- limma::makeContrasts(contrasts = contrast_formulas, levels = design)
  colnames(contrast_matrix) <- names(contrast_formulas)
  
  fit2 <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr_imp, design), contrast_matrix))
  
  ann <- align_annotations_to_expr(expr_imp, prot_tbl, protein_id_col, annot_cols)
  feature_id <- ann[[protein_id_col]]
  annot_out <- setdiff(annot_cols, protein_id_col)
  
  de_tables <- lapply(colnames(contrast_matrix), function(cn) {
    de <- limma::topTable(fit2, coef = cn, adjust.method="BH", sort.by="none", number=Inf)
    de <- align_de_to_expr(de, expr_imp, contrast_name = cn)
    data.frame(
      FeatureID = feature_id,
      Contrast  = cn,
      ann[, annot_out, drop = FALSE],
      de[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")],
      check.names = FALSE
    )
  })
  
  # IMPORTANT: names(de_tables) are used by summarize_limma_mult_imputation()
  names(de_tables) <- colnames(contrast_matrix)
 
  list(
    meta_aligned = meta_aligned,
    design = design,
    contrast_formulas = contrast_formulas,
    contrast_matrix = contrast_matrix,
    fit2 = fit2,
    de_tables = de_tables
  )
}
