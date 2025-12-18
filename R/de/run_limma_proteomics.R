run_limma_proteomics <- function(expr_imp, meta, contrasts_df, prot_tbl, cfg) {
  stopifnot(is.matrix(expr_imp))
  stopifnot(all(c("SampleID","Condition") %in% colnames(meta)))
  stopifnot(all(c("Contrast_name","Factor","Numerator","Denominator") %in% colnames(contrasts_df)))
  
  # 1) align meta to expr columns
  meta_aligned <- meta[match(colnames(expr_imp), meta$SampleID), , drop = FALSE]
  stopifnot(all(meta_aligned$SampleID == colnames(expr_imp)))
  
  # 2) design (Condition only)  ~0+cond  (matches old script)
  meta_aligned$Condition <- factor(meta_aligned$Condition)
  design <- model.matrix(~ 0 + Condition, data = meta_aligned)
  colnames(design) <- levels(meta_aligned$Condition)
  
  # 3) contrasts (Numerator - Denominator) + keep Contrast_name (matches old script)
  stopifnot(all(contrasts_df$Factor == "Condition"))
  contrast_formulas <- paste(contrasts_df$Numerator, contrasts_df$Denominator, sep = " - ")
  names(contrast_formulas) <- contrasts_df$Contrast_name
  
  contrast.matrix <- limma::makeContrasts(contrasts = contrast_formulas, levels = design)
  colnames(contrast.matrix) <- names(contrast_formulas)  # <-- match old script Contrast_name
  
  # 4) fit
  fit  <- limma::lmFit(expr_imp, design)
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2)
  
  # 5) Build per-contrast annotated topTables (match old: sort.by="none", number=Inf)
  protein_id_col <- cfg$modes$proteomics$id_columns$protein_id
  stopifnot(protein_id_col %in% colnames(prot_tbl))
  idx <- match(rownames(expr_imp), prot_tbl[[protein_id_col]])
  
  if (anyNA(idx)) {
    stop(sprintf("prot_tbl missing %d proteins (e.g. %s)",
                 sum(is.na(idx)), paste(head(rownames(expr_imp)[is.na(idx)], 3), collapse = ", ")))
  }
  
  
  ann <- prot_tbl[idx, c("Protein.Group","Protein.Names","Genes","First.Protein.Description"), drop = FALSE]
  
  de_tables <- lapply(colnames(contrast.matrix), function(cn) {
    de <- limma::topTable(fit2, coef = cn, adjust.method = "BH", sort.by = "none", number = Inf)
    
    de$FeatureID <- ann$Protein.Group
    de$Contrast  <- cn
    
    de_out <- cbind(
      FeatureID = de$FeatureID,
      Contrast  = de$Contrast,
      ann[, c("Protein.Names","Genes","First.Protein.Description"), drop = FALSE],
      de[, c("logFC","AveExpr","t","P.Value","adj.P.Val","B")]
    )
    
    de_out
  })
  names(de_tables) <- colnames(contrast.matrix)

  # return only objects (no files yet)
  list(
    meta_aligned     = meta_aligned,
    design           = design,
    contrast_formulas= contrast_formulas,
    contrast.matrix  = contrast.matrix,
    fit2             = fit2,
    de_tables = de_tables
  )
}
