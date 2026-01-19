# R/pipeline/mod_proteomics_de.R

#' Proteomics DE module (orchestration only)
#'
#' Runs DE according to cfg$de$method (currently supports "limma").
#' No file I/O here.
#'
#' @param pre     output of preprocess_proteomics()
#' @param inputs  output of load_proteomics_inputs()
#' @param config  full config
#' @param verbose logical
#' @return list(method, runs, runs_de_tables, summary_df)
mod_proteomics_de <- function(pre, inputs, config, verbose = FALSE) {
  assert_pre_contract(pre, stage = "proteomics")
  
  cfg <- config$modes$proteomics
  
  # ---- choose DE method ----
  method <- cfg$de$method %||% "limma"
  
  if (identical(method, "limma")) {
    
    # 1) imputations
    imputations <- make_imputations_proteomics(
      expr_mat = pre$expr_filt,
      cfg      = config,
      verbose  = verbose
    )
    
    # 2) run limma on imputations
    runs <- run_limma_multimp(
      imputations  = imputations,
      meta         = pre$meta,
      contrasts_df = inputs$contrasts,
      prot_tbl     = inputs$protein,
      cfg          = config,
      verbose      = verbose
    )
    
    # 3) collect DE tables + summarize
    runs_de_tables <- lapply(runs, function(x) x$de_tables)
    
    summary_df <- summarize_limma_mult_imputation(
      runs_de_tables = runs_de_tables,
      config         = config
    )
    
    return(list(
      method        = "limma",
      imputations   = imputations,    
      runs          = runs,
      runs_de_tables = runs_de_tables,
      summary_df    = summary_df
    ))
  }
  
  stop("Unsupported proteomics DE method: ", method)
}
#' Run limma differential analysis for proteomics with contrast support
#'
#' Fits a limma linear model on an imputed proteomics expression matrix and
#' returns per-contrast result tables with feature annotations.
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
  meta_aligned <- align_meta_to_expr(expr_imp, meta, p_cfg)
  
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
    expr_imp = expr_imp,
    meta_aligned = meta_aligned,
    design = design,
    contrast_formulas = contrast_formulas,
    contrast_matrix = contrast_matrix,
    fit2 = fit2,
    de_tables = de_tables
  )
}


#' Run limma proteomics DE on a precomputed list of imputed datasets
#'
#' Use this function when imputations were already generated and you want to:
#' - avoid recomputing imputations (e.g., when trying a different DE method), or
#' - decouple imputation from DE in a `{targets}` pipeline.
#'
#' @return A list of length `length(imputations)`. Each element is the output of
#'         `run_limma_proteomics()` for the corresponding imputed dataset.

run_limma_multimp <- function(imputations, meta, contrasts_df, prot_tbl, cfg,
                              verbose = FALSE) {
  
  validate_proteomics_imputations(imputations = imputations,
                                  meta = meta, cfg = cfg)
  
  runs <- vector("list", length(imputations))
  
  for (i in seq_along(imputations)) {
    if (isTRUE(verbose)) {
      message(sprintf("Limma on imputations: %d / %d", i, length(imputations)))
    }
    
    runs[[i]] <- run_limma_proteomics(
      expr_imp     = imputations[[i]],
      meta         = meta,
      contrasts_df = contrasts_df,
      prot_tbl     = prot_tbl,
      cfg          = cfg
    )
  }
  
  runs
}



