#' ANOVA + Tukey HSD for proteomics (multi-group DE)
#'
#' For >2 groups, runs one-way ANOVA per protein, then Tukey HSD for
#' pairwise comparisons. Returns results in the same structure as
#' run_limma_proteomics() for pipeline compatibility.
#'
#' @param expr_imp numeric matrix (proteins x samples), imputed
#' @param meta     data.frame with sample metadata
#' @param contrasts_df data.frame with Contrast_name, Factor, Numerator, Denominator
#' @param prot_tbl data.frame with protein annotations
#' @param cfg      full config list
#' @return list with de_tables (named list of per-contrast data.frames)
run_anova_posthoc <- function(expr_imp, meta, contrasts_df, prot_tbl, cfg) {
    p_cfg <- cfg$modes$proteomics

    sample_col <- p_cfg$effects$samples %||% "SampleID"
    group_col  <- p_cfg$effects$color %||% "Condition"
    protein_id_col <- p_cfg$id_columns$protein_id %||% "Protein.Group"
    p_adjust_method <- p_cfg$de$p_adjust_method %||% "BH"

    default_annot <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description")
    annot_cols <- unique(c(protein_id_col, p_cfg$id_columns$protein_annot %||% default_annot))
    annot_cols <- as.character(unlist(annot_cols))
    annot_cols <- annot_cols[nzchar(annot_cols)]

    assert_numeric_matrix(expr_imp, "expr_imp")
    meta_aligned <- align_meta_to_expr(expr_imp, meta, p_cfg)

    group <- as.character(meta_aligned[[group_col]])
    names(group) <- meta_aligned[[sample_col]]

    ann <- align_annotations_to_expr(expr_imp, prot_tbl, protein_id_col, annot_cols)
    feature_id <- ann[[protein_id_col]]
    annot_out <- setdiff(annot_cols, protein_id_col)

    de_table_cfg <- p_cfg$de_table %||% list()
    target_id_col <- de_table_cfg$id_col %||% "FeatureID"

   
    # Pre-compute factor + contrast lookup keys once (outside any loop)
    grp_factor <- factor(group)
    n_proteins <- nrow(expr_imp)
    n_contrasts <- nrow(contrasts_df)
    contrast_names <- contrasts_df$Contrast_name
    stopifnot(length(contrast_names) > 0)
    
    # Pre-compute Tukey row keys for each contrast (forward + reverse direction)
    contrast_fwd <- paste0(contrasts_df$Numerator, "-", contrasts_df$Denominator)
    contrast_rev <- paste0(contrasts_df$Denominator, "-", contrasts_df$Numerator)
    
    # Pre-allocate result containers
    anova_pvals <- numeric(n_proteins)
    logFC_mat <- matrix(NA_real_, nrow = n_proteins, ncol = n_contrasts)
    pval_mat <- matrix(NA_real_, nrow = n_proteins, ncol = n_contrasts)
    
    # Single per-protein loop: fit aov + Tukey ONCE, extract all contrasts
    for (i in seq_len(n_proteins)) {
      vals <- as.numeric(expr_imp[i, ])
      fit <- tryCatch(
        stats::aov(vals ~ grp_factor),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        anova_pvals[i] <- NA_real_
        next
      }
      
      # Omnibus ANOVA p-value
      sf <- summary(fit)
      anova_pvals[i] <- sf[[1]][["Pr(>F)"]][1]
      
      # Tukey HSD: all pairwise comparisons in one shot
      tukey <- tryCatch(stats::TukeyHSD(fit), error = function(e) NULL)
      if (is.null(tukey)) next
      
      # tukey is a named list with one element per predictor; we have one predictor
      tukey_df <- as.data.frame(tukey[[1]])
      tukey_rownames <- rownames(tukey_df)
      
      # Extract each requested contrast from the single Tukey table
      for (ci in seq_len(n_contrasts)) {
        if (contrast_fwd[ci] %in% tukey_rownames) {
          logFC_mat[i, ci] <- tukey_df[contrast_fwd[ci], "diff"]
          pval_mat[i, ci] <- tukey_df[contrast_fwd[ci], "p adj"]
        } else if (contrast_rev[ci] %in% tukey_rownames) {
          # Reverse direction: flip the sign of diff
          logFC_mat[i, ci] <- -tukey_df[contrast_rev[ci], "diff"]
          pval_mat[i, ci] <- tukey_df[contrast_rev[ci], "p adj"]
        }
        # else: leave as NA (group level missing from Tukey output)
      }
    }
    
    # Pre-compute average expression once (vectorized)
    ave_expr <- rowMeans(expr_imp, na.rm = TRUE)
    
    # Build per-contrast tables from pre-computed matrices (no nested loops)
    de_tables <- lapply(seq_along(contrast_names), function(ci) {
      cn <- contrast_names[ci]
      logFC_vals <- logFC_mat[, ci]
      p_vals <- pval_mat[, ci]
      adj_p_vals <- stats::p.adjust(p_vals, method = p_adjust_method)
      
      df_out <- data.frame(
        TEMP_ID_COL = feature_id,
        Contrast = cn,
        ann[, annot_out, drop = FALSE],
        logFC = logFC_vals,
        AveExpr = ave_expr,
        t = NA_real_,
        P.Value = p_vals,
        adj.P.Val = adj_p_vals,
        B = NA_real_,
        check.names = FALSE
      )
      names(df_out)[names(df_out) == "TEMP_ID_COL"] <- target_id_col
      df_out
    })
    
    names(de_tables) <- contrast_names

    list(
      expr_imp = expr_imp,
      meta_aligned = meta_aligned,
      de_tables = de_tables,
      anova_pvals = anova_pvals
    )
}
