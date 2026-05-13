#' t-test / Welch test DE for proteomics
#'
#' Per-protein two-sample t-test (equal or unequal variance) with optional
#' paired test support. Returns results in the same structure as
#' run_limma_proteomics() for pipeline compatibility.
#'
#' @param expr_imp numeric matrix (proteins x samples), imputed
#' @param meta     data.frame with sample metadata
#' @param contrasts_df data.frame with Contrast_name, Factor, Numerator, Denominator
#' @param prot_tbl data.frame with protein annotations
#' @param cfg      full config list
#' @param var_equal logical, TRUE for Student's t-test, FALSE for Welch
#' @return list with de_tables (named list of per-contrast data.frames)
run_ttest_de <- function(expr_imp, meta, contrasts_df, prot_tbl, cfg, var_equal = TRUE) {
    p_cfg <- cfg$modes$proteomics

    sample_col <- p_cfg$effects$samples %||% "SampleID"
    group_col  <- p_cfg$effects$color %||% "Condition"
    protein_id_col <- p_cfg$id_columns$protein_id %||% "Protein.Group"
    p_adjust_method <- p_cfg$de$p_adjust_method %||% "BH"
    paired <- isTRUE(p_cfg$de$paired)
    pairing_col <- p_cfg$de$pairing_col

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

    n_proteins <- nrow(expr_imp)
    contrast_names <- contrasts_df$Contrast_name
    stopifnot(length(contrast_names) > 0)

    de_tables <- lapply(seq_along(contrast_names), function(ci) {
        cn <- contrast_names[ci]
        num <- contrasts_df$Numerator[ci]
        den <- contrasts_df$Denominator[ci]

        idx_num <- which(group == num)
        idx_den <- which(group == den)

        if (length(idx_num) < 2 || length(idx_den) < 2) {
            stop(sprintf("Contrast '%s': need at least 2 samples per group (got %d vs %d).",
                         cn, length(idx_num), length(idx_den)))
        }

        # For paired tests, align samples by pairing column
        is_paired <- paired && !is.null(pairing_col) && pairing_col %in% colnames(meta_aligned)
        if (is_paired) {
            pair_num <- meta_aligned[[pairing_col]][idx_num]
            pair_den <- meta_aligned[[pairing_col]][idx_den]
            common_pairs <- intersect(pair_num, pair_den)
            if (length(common_pairs) < 2) {
                stop(sprintf("Paired test for '%s': fewer than 2 common pairs found.", cn))
            }
            # Reorder to match pairs
            idx_num <- idx_num[match(common_pairs, pair_num)]
            idx_den <- idx_den[match(common_pairs, pair_den)]
        }

        logFC_vals <- numeric(n_proteins)
        p_vals <- numeric(n_proteins)
        t_vals <- numeric(n_proteins)

        for (i in seq_len(n_proteins)) {
            x_num <- as.numeric(expr_imp[i, idx_num])
            x_den <- as.numeric(expr_imp[i, idx_den])

            tt <- tryCatch(
                stats::t.test(x_num, x_den, var.equal = var_equal, paired = is_paired),
                error = function(e) NULL
            )

            if (is.null(tt)) {
                logFC_vals[i] <- mean(x_num, na.rm = TRUE) - mean(x_den, na.rm = TRUE)
                p_vals[i] <- NA_real_
                t_vals[i] <- NA_real_
            } else {
                logFC_vals[i] <- mean(x_num, na.rm = TRUE) - mean(x_den, na.rm = TRUE)
                p_vals[i] <- tt$p.value
                t_vals[i] <- tt$statistic
            }
        }

        adj_p_vals <- stats::p.adjust(p_vals, method = p_adjust_method)

        df_out <- data.frame(
            TEMP_ID_COL = feature_id,
            Contrast = cn,
            ann[, annot_out, drop = FALSE],
            logFC = logFC_vals,
            AveExpr = rowMeans(expr_imp, na.rm = TRUE),
            t = t_vals,
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
