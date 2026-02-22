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

    assert_numeric_matrix(expr_imp, "expr_imp")
    meta_aligned <- align_meta_to_expr(expr_imp, meta, p_cfg)

    group <- as.character(meta_aligned[[group_col]])
    names(group) <- meta_aligned[[sample_col]]

    ann <- align_annotations_to_expr(expr_imp, prot_tbl, protein_id_col, annot_cols)
    feature_id <- ann[[protein_id_col]]
    annot_out <- setdiff(annot_cols, protein_id_col)

    de_table_cfg <- p_cfg$de_table %||% list()
    target_id_col <- de_table_cfg$id_col %||% "FeatureID"

    # Per-protein ANOVA
    n_proteins <- nrow(expr_imp)
    anova_pvals <- numeric(n_proteins)

    for (i in seq_len(n_proteins)) {
        vals <- as.numeric(expr_imp[i, ])
        grp <- group
        fit <- tryCatch(
            stats::aov(vals ~ factor(grp)),
            error = function(e) NULL
        )
        if (is.null(fit)) {
            anova_pvals[i] <- NA_real_
        } else {
            sf <- summary(fit)
            anova_pvals[i] <- sf[[1]][["Pr(>F)"]][1]
        }
    }

    # Build per-contrast tables using Tukey HSD results
    contrast_names <- contrasts_df$Contrast_name
    stopifnot(length(contrast_names) > 0)

    de_tables <- lapply(seq_along(contrast_names), function(ci) {
        cn <- contrast_names[ci]
        num <- contrasts_df$Numerator[ci]
        den <- contrasts_df$Denominator[ci]

        logFC_vals <- numeric(n_proteins)
        p_vals <- numeric(n_proteins)

        for (i in seq_len(n_proteins)) {
            vals <- as.numeric(expr_imp[i, ])
            grp <- factor(group)
            fit <- tryCatch(stats::aov(vals ~ grp), error = function(e) NULL)
            if (is.null(fit)) {
                logFC_vals[i] <- NA_real_
                p_vals[i] <- NA_real_
                next
            }
            tukey <- tryCatch(stats::TukeyHSD(fit), error = function(e) NULL)
            if (is.null(tukey)) {
                logFC_vals[i] <- NA_real_
                p_vals[i] <- NA_real_
                next
            }
            tukey_df <- as.data.frame(tukey$grp)
            # Tukey names rows like "B-A"; we need "num-den"
            target_row <- paste0(num, "-", den)
            rev_row <- paste0(den, "-", num)

            if (target_row %in% rownames(tukey_df)) {
                logFC_vals[i] <- tukey_df[target_row, "diff"]
                p_vals[i] <- tukey_df[target_row, "p adj"]
            } else if (rev_row %in% rownames(tukey_df)) {
                logFC_vals[i] <- -tukey_df[rev_row, "diff"]
                p_vals[i] <- tukey_df[rev_row, "p adj"]
            } else {
                logFC_vals[i] <- NA_real_
                p_vals[i] <- NA_real_
            }
        }

        adj_p_vals <- stats::p.adjust(p_vals, method = p_adjust_method)

        df_out <- data.frame(
            TEMP_ID_COL = feature_id,
            Contrast = cn,
            ann[, annot_out, drop = FALSE],
            logFC = logFC_vals,
            AveExpr = rowMeans(expr_imp, na.rm = TRUE),
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
        de_tables = de_tables
    )
}
