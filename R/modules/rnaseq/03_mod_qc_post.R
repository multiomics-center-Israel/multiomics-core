#' RNA-seq QC post-DE module
#'
#' Writes QC artifact plots that depend on DE results (Volcano, MA, Heatmap).
#' Mirrors Proteomics logic using core plotting functions.
#'
#' @param pre        preprocessed RNA object
#' @param de_res     DE results list
#' @param config     full YAML config list
#' @param out_dir    output run directory
#' @return list(files)
mod_rnaseq_qc_post <- function(pre, de_res, config, out_dir) {
    assert_pre_contract(pre, stage = "rna")
    dirs <- create_legacy_output_dirs(out_dir, create = TRUE)
    out_qc_post <- file.path(dirs$diagnostic_plots, "DE_plots")
    ensure_dir(out_qc_post)

    mat <- pre$expr_work
    meta <- pre$meta
    cfg <- config$modes$rna

    # Extract primary color (handle array config for multi-color PCA)
    color_config <- cfg$effects$color %||% "Group"
    eff_color <- as.character(color_config[[1]])

    files <- character(0)
    plots <- list()

    # 1. Ensure summary_df exists
    sumdf <- de_res$summary_df
    if (is.null(sumdf)) {
        sumdf <- tryCatch(build_rnaseq_summary_df(de_res$tables, config$modes$rna$de), error = function(e) NULL)
        if (is.null(sumdf)) {
            warning("RNA QC-post: Could not build summary_df. Skipping QC post plots.")
            return(list(files = character(0)))
        }
    }

    # 2. Extract per-contrast tables using the Shim
    # (Helper defined below)
    tables <- get_rna_de_tables_qc_post(sumdf, cfg, mat)

    # 3. Generate Volcano and MA plots
    id_col_cfg <- cfg$de$id_col %||% "Gene" # or FeatureID

    for (cn in names(tables)) {
        de_tbl <- tables[[cn]]
        if (is.null(de_tbl) || nrow(de_tbl) == 0) {
            message(sprintf("[QC_post] Skipping contrast '%s': empty DE table.", cn))
            next
        }

        # Validation for P-values
        if (all(is.na(de_tbl$adj.P.Val))) {
            warning(sprintf("[QC_post] Skipping contrast '%s': all adj.P.Val are NA.", cn))
            next
        }

        # --- Volcano ---
        # "Y-axis must extend beyond 5–10 for real DE" -> Checking this programmatically is hard without stopping valid runs.
        # But we can define the threshold
        max_nlp <- max(-log10(de_tbl$adj.P.Val), na.rm = TRUE)
        if (!is.infinite(max_nlp) && max_nlp < 1.3) { # < 1.3 means no p < 0.05
            message(sprintf("[QC_post] Note: Contrast '%s' has low significance (min adj.P.Val > 0.05). Volcano might look flat.", cn))
        }

        # Emit one volcano per p-value type (padj + raw pval)
        for (ptype in c("padj", "pval")) {
            p_volcano <- plot_volcano(
                de_tbl,
                cfg = cfg,
                title = paste0("Volcano: ", cn,
                               " (", if (ptype == "padj") "adj.P.Val" else "P.Value", ")"),
                pvalue_type = ptype
            )
            f_volcano <- file.path(out_qc_post, sprintf("volcano_%s_%s.png", cn, ptype))
            ggplot2::ggsave(f_volcano, plot = p_volcano, width = 8, height = 6, dpi = 150)
            files <- c(files, f_volcano)
            plots[[paste0("volcano_", cn, "_", ptype)]] <- p_volcano
        }
        # Backwards-compat alias
        f_legacy <- file.path(out_qc_post, sprintf("volcano_%s.png", cn))
        f_padj   <- file.path(out_qc_post, sprintf("volcano_%s_padj.png", cn))
        if (file.exists(f_padj) && !file.exists(f_legacy)) {
            file.copy(f_padj, f_legacy, overwrite = TRUE)
            files <- c(files, f_legacy)
            plots[[paste0("volcano_", cn)]] <- plots[[paste0("volcano_", cn, "_padj")]]
        }

        # --- MA Plot ---
        # Requires 'A' / 'AveExpr'. Our Shim ensures this exists.
        if (all(is.na(de_tbl$AveExpr))) {
            warning(sprintf("[QC_post] Skipping MA for '%s': AveExpr all NA.", cn))
            next
        }

        p_ma <- plot_ma(
            de_tbl,
            cfg = cfg,
            title = paste0("MA: ", cn),
            use_adj = TRUE
        )

        # "Clamp Y-range to avoid outlier compression"
        # Since we can't modify core plot_ma, we check if we can modify the plot object (ggplot)
        # We can add coordinate limits to the ggplot object
        y_lims <- stats::quantile(de_tbl$logFC, c(0.01, 0.99), na.rm = TRUE)
        # Expand slightly
        y_span <- diff(y_lims)
        y_lims <- y_lims + c(-0.1, 0.1) * y_span

        # If the range is super small (e.g. 0), default to some reasonable range
        if (diff(y_lims) < 0.5) y_lims <- c(-1, 1)

        p_ma <- p_ma + ggplot2::coord_cartesian(ylim = y_lims)

        f_ma <- file.path(out_qc_post, sprintf("ma_%s.png", cn))
        ggplot2::ggsave(f_ma, plot = p_ma, width = 8, height = 6, dpi = 150)
        files <- c(files, f_ma)
        plots[[paste0("ma_", cn)]] <- p_ma
    }

    # 4. Global Top DE Heatmap (Legacy/Bonus feature)
    de_ids <- get_rna_de_features_qc(sumdf)
    de_ids <- de_ids[!is.na(de_ids)]

    # Get max top DE features from config (default 100)
    qc_post_cfg <- cfg$qc_post %||% list()
    max_top_de <- qc_post_cfg$max_top_de_features %||% 100

    if (length(de_ids) >= 2) {
        # choose top N by min padj across contrasts
        padj_cols <- grep("^padj\\.", names(sumdf), value = TRUE)
        min_padj <- apply(sumdf[, padj_cols, drop = FALSE], 1, function(x) {
            if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
        })

        valid_mask <- !is.na(min_padj)
        valid_count <- sum(valid_mask)

        if (valid_count == 0) {
            warning("[QC_post] No features with valid padj values for top DE heatmap")
        } else {
            valid_padj <- min_padj[valid_mask]
            valid_ids <- sumdf$FeatureID[valid_mask]
            ord <- order(valid_padj)

            top_n <- min(max_top_de, length(ord))
            top_ids <- valid_ids[ord[seq_len(top_n)]]

            message(sprintf(
                "[QC_post] Selecting %d/%d features for top DE heatmap (requested: %d)",
                top_n, valid_count, max_top_de
            ))

            m <- mat[intersect(top_ids, rownames(mat)), , drop = FALSE]
            if (nrow(m) >= 2) {
                # z-score by gene (for visualization only)
                m_z <- t(scale(t(m)))

                # Replace NaN/Inf with 0 (NaN from SD=0, Inf from near-zero SD)
                m_z[!is.finite(m_z)] <- 0

                # Build column annotation, validating columns exist in meta
                annot_cols <- cfg$effects$heatmap_annotations %||% NULL
                if (!is.null(annot_cols)) {
                    # Filter to columns that actually exist in meta
                    valid_annot_cols <- intersect(annot_cols, colnames(meta))
                    if (length(valid_annot_cols) == 0) {
                        warning(sprintf(
                            "[QC_post] heatmap_annotations columns [%s] not found in metadata. Using color column.",
                            paste(annot_cols, collapse = ", ")
                        ))
                        ann <- data.frame(Condition = meta[[eff_color]])
                        rownames(ann) <- meta[[cfg$effects$samples]]
                    } else {
                        if (length(valid_annot_cols) < length(annot_cols)) {
                            missing <- setdiff(annot_cols, valid_annot_cols)
                            warning(sprintf(
                                "[QC_post] Annotation columns not found in metadata: %s",
                                paste(missing, collapse = ", ")
                            ))
                        }
                        ann <- meta[, valid_annot_cols, drop = FALSE]
                        rownames(ann) <- meta[[cfg$effects$samples]]
                    }
                } else {
                    ann <- data.frame(Condition = meta[[eff_color]])
                    rownames(ann) <- meta[[cfg$effects$samples]]
                }

                # Align annotation to matrix columns (only keep matching samples)
                common_samples <- intersect(colnames(m_z), rownames(ann))
                if (length(common_samples) == 0) {
                    warning("[QC_post] No matching sample IDs between expression matrix and annotation. Skipping heatmap.")
                } else {
                    ann <- ann[common_samples, , drop = FALSE]
                    m_z <- m_z[, common_samples, drop = FALSE]

                    hm_png <- file.path(dirs$diagnostic_plots, "heatmap_top_DE.png")

                    tryCatch(
                        {
                            hm_obj <- plot_heatmap_core(
                                expr_mat = m_z,
                                annotation_col = ann,
                                title = sprintf("Top DE (%d) - Z-score", nrow(m_z)),
                                scale_rows = FALSE, # already z-scored
                                cluster_cols = TRUE,
                                max_rows = max_top_de
                            )
                            save_heatmap_to_file(hm_obj, hm_png, width = 2000, height = 1400, res = 150)
                            files <- c(files, hm_png)
                        },
                        error = function(e) {
                            warning(sprintf("[QC_post] heatmap_top_DE failed: %s", conditionMessage(e)))
                        }
                    )
                }
            }
        } # Close valid_count > 0 check
    }

    list(files = unique(files), plots = plots)
}

# ---- Helpers ----

#' Shim: Extract & Normalise per-contrast DE tables from summary_df
#' Returns list: contrast -> data.frame(FeatureID, logFC, P.Value, adj.P.Val, AveExpr)
get_rna_de_tables_qc_post <- function(summary_df, cfg, expr_mat) {
    stopifnot(is.data.frame(summary_df))
    # RNA summary always has FeatureID (fixed previously)
    src_id_col <- "FeatureID"

    # FIX 2: Updated to handle new column naming (removed ".imputs.")
    # RNA summary now uses: "padj.<Contrast>", "pvalue.<Contrast>", "linearFC.<Contrast>"

    cols <- colnames(summary_df)
    padj_cols <- grep("^padj\\.", cols, value = TRUE)

    if (length(padj_cols) == 0) {
        return(list())
    }

    # Extract contrast names from padj columns
    contrasts <- sub("^padj\\.", "", padj_cols)

    out <- setNames(vector("list", length(contrasts)), contrasts)

    for (cn in contrasts) {
        # RNA column names after FIX 2
        p_col_adj <- paste0("padj.", cn)
        p_col_raw <- paste0("pvalue.", cn)
        fc_col <- paste0("linearFC.", cn)

        if (!all(c(p_col_adj, fc_col) %in% colnames(summary_df))) next

        # Convert linearFC to log2FC for plotting
        lin_fc <- summary_df[[fc_col]]
        log_fc <- log2(abs(lin_fc)) * sign(lin_fc)

        # Build Table
        de_tbl <- data.frame(
            FeatureID = summary_df[[src_id_col]],
            logFC = log_fc,
            P.Value = if (p_col_raw %in% colnames(summary_df)) summary_df[[p_col_raw]] else summary_df[[p_col_adj]],
            adj.P.Val = summary_df[[p_col_adj]],
            stringsAsFactors = FALSE
        )

        # Add 'A' (AveExpr) from expression matrix
        # This is CRITICAL for plot_ma
        de_tbl <- add_A_from_expr(de_tbl, expr_mat, id_col = "FeatureID")
        # Ensure column is named 'AveExpr' for plot_ma (add_A_from_expr adds '.A')
        names(de_tbl)[names(de_tbl) == ".A"] <- "AveExpr"

        out[[cn]] <- de_tbl
    }
    out
}

get_rna_de_features_qc <- function(summary_df) {
    if (is.null(summary_df) || nrow(summary_df) == 0) {
        return(character(0))
    }

    # Check for pass_any_contrast column (handles both boolean TRUE and numeric 1)
    if ("pass_any_contrast" %in% colnames(summary_df)) {
        pass_vals <- summary_df$pass_any_contrast
        # Handle both boolean TRUE and numeric 1
        is_pass <- !is.na(pass_vals) & (pass_vals == TRUE | pass_vals == 1)
        return(summary_df$FeatureID[which(is_pass)])
    }

    # Fallback 1: RNA-seq style _pass columns (e.g., "ContrastA_pass")
    pass_cols <- grep("_pass$", colnames(summary_df), value = TRUE)
    if (length(pass_cols) > 0) {
        # Convert to numeric matrix for rowSums (handles TRUE/FALSE and 1/0)
        pass_mat <- as.matrix(summary_df[, pass_cols, drop = FALSE])
        mode(pass_mat) <- "numeric"
        row_sums <- rowSums(pass_mat, na.rm = TRUE)
        return(summary_df$FeatureID[which(row_sums > 0)])
    }

    # Fallback 2: Proteomics style sum.pass columns (unlikely for RNA)
    pass_cols <- grep("^sum\\.pass\\.", colnames(summary_df), value = TRUE)
    if (length(pass_cols) > 0) {
        row_sums <- rowSums(summary_df[, pass_cols, drop = FALSE], na.rm = TRUE)
        return(summary_df$FeatureID[which(row_sums > 0)])
    }

    character(0)
}
