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
    out_qc_post <- file.path(dirs$diagnostic_plots, "QC_post")
    ensure_dir(out_qc_post)

    mat <- pre$expr_work
    meta <- pre$meta
    cfg <- config$modes$rna
    eff_color <- cfg$effects$color %||% "Group"

    files <- character(0)
    plots <- list()

    # 1. Ensure summary_df exists
    sumdf <- de_res$summary_df
    if (is.null(sumdf)) {
        sumdf <- tryCatch(build_rnaseq_summary_df(de_res$tables, config), error = function(e) NULL)
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

        # Passing use_adj=TRUE implies core looks for adj.P.Val
        p_volcano <- plot_volcano(
            de_tbl,
            cfg = cfg,
            title = paste0("Volcano: ", cn),
            use_adj = TRUE
        )
        f_volcano <- file.path(out_qc_post, sprintf("volcano_%s.png", cn))
        ggplot2::ggsave(f_volcano, plot = p_volcano, width = 8, height = 6, dpi = 150)
        files <- c(files, f_volcano)
        plots[[paste0("volcano_", cn)]] <- p_volcano

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
    # Recalculate de_ids logic
    de_ids <- get_rna_de_features_qc(sumdf)
    de_ids <- de_ids[!is.na(de_ids)]

    if (length(de_ids) >= 2) {
        # choose top N by min padj across contrasts
        padj_cols <- grep("^padj\\.", names(sumdf), value = TRUE)
        min_padj <- apply(sumdf[, padj_cols, drop = FALSE], 1, function(x) {
            if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
        })
        ord <- order(min_padj, na.last = NA)
        top_n <- min(2000L, length(ord))
        top_ids <- sumdf$FeatureID[ord][seq_len(top_n)]

        m <- mat[intersect(top_ids, rownames(mat)), , drop = FALSE]
        if (nrow(m) > 2) {
            # z-score by gene (for visualization only)
            m_z <- t(scale(t(m)))
            m_z <- m_z[complete.cases(m_z), , drop = FALSE]

            ann <- data.frame(Condition = meta[[eff_color]])
            rownames(ann) <- rownames(meta)

            # Align annotation to matrix columns
            ann <- ann[colnames(m_z), , drop = FALSE]

            hm_png <- file.path(dirs$diagnostic_plots, "heatmap_top_DE.png")

            grDevices::png(hm_png, width = 1400, height = 1200, res = 150)
            pheatmap::pheatmap(
                m_z,
                show_rownames = FALSE,
                annotation_col = ann,
                main = sprintf("Top DE (%d) - Z-score", nrow(m_z)),
                cluster_cols = TRUE
            )
            grDevices::dev.off()
            files <- c(files, hm_png)
        }
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

    # Identify contrasts from columns like padj.imputs.<Contrast>
    # Note: RNA might use "padj.<Contrast>" directly if not using "imputs" naming convention everywhere
    # But current RNA summary (domain/rnaseq/04_de_summary) uses "padj.<name>"

    # Check column patterns for RNA summary DF:
    # It constructs cols like: "padj.<Contrast>", "pvalue.<Contrast>", "log2FoldChange.<Contrast>"
    # Or strict Proteomics compatibility might have imposed "padj.imputs."?
    # Let's check the patterns present.

    # Try generic pattern detection
    cols <- colnames(summary_df)
    padj_cols <- grep("^padj\\.", cols, value = TRUE)
    padj_cols <- padj_cols[!grep("^padj\\.imputs\\.", padj_cols)] # exclude if mix

    # If using Proteomics style "padj.imputs.", adapt
    padj_imp <- grep("^padj\\.imputs\\.", cols, value = TRUE)

    if (length(padj_imp) > 0) {
        # It's using imputation style naming (unlikely for pure RNA but possible if shimmed)
        contrasts <- sub("^padj\\.imputs\\.", "", padj_imp)
        style <- "imputs"
    } else {
        # RNA native
        contrasts <- sub("^padj\\.", "", padj_cols)
        style <- "native"
    }

    if (length(contrasts) == 0) {
        return(list())
    }

    out <- setNames(vector("list", length(contrasts)), contrasts)

    for (cn in contrasts) {
        if (style == "imputs") {
            p_col_adj <- paste0("padj.imputs.", cn)
            p_col_raw <- paste0("pvalue.imputs.", cn)
            fc_col <- paste0("linearFC.imputs.", cn) # Proteomics uses linearFC
        } else {
            p_col_adj <- paste0("padj.", cn)
            p_col_raw <- paste0("pvalue.", cn)
            fc_col <- paste0("log2FoldChange.", cn) # RNA native
        }

        if (!all(c(p_col_adj, fc_col) %in% colnames(summary_df))) next

        # Extract base values
        if (style == "imputs") {
            lin_fc <- summary_df[[fc_col]]
            log_fc <- log2(abs(lin_fc)) * sign(lin_fc)
        } else {
            log_fc <- summary_df[[fc_col]]
        }

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
    if ("pass_any_contrast" %in% colnames(summary_df)) {
        return(summary_df$FeatureID[which(summary_df$pass_any_contrast == 1)])
    }
    pass_cols <- grep("^sum\\.pass\\.", colnames(summary_df), value = TRUE)
    if (length(pass_cols) > 0) {
        row_sums <- rowSums(summary_df[, pass_cols, drop = FALSE], na.rm = TRUE)
        return(summary_df$FeatureID[which(row_sums > 0)])
    }
    character(0)
}
