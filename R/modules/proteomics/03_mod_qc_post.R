#' Proteomics QC post-DE module
#'
#' Writes QC artifacts that depend on DE results.
#' Supports DE source:
#' - "summary": builds per-contrast tables from summary_df (wide)
#' - "table1" : uses first imputation run tables (legacy-ish)
#'
#' @param pre        preprocessed proteomics object
#' @param de_res     DE results list
#' @param config     full YAML config list
#' @param out_dir    output run directory
#' @param de_source  "summary" or "table1"
#' @return list(plots, files)
mod_proteomics_qc_post <- function(pre, de_res, config, out_dir, de_source = c("summary", "table1")) {
    de_source <- match.arg(de_source)

    cfg <- config$modes$proteomics

    # Check if QC post is enabled
    if (isFALSE(cfg$qc_post$enabled)) {
        message("QC post-DE disabled in config.")
        return(list(plots = list(), files = character(0)))
    }

    dirs <- create_legacy_output_dirs(out_dir)
    out_qc_post <- file.path(dirs$diagnostic_plots, "DE_plots")
    ensure_dir(out_qc_post)

    plots <- list()
    files <- character(0)

    # Read plot/output toggles from config
    do_volcano <- isTRUE(cfg$qc_post$plots$volcano %||% TRUE)
    do_ma      <- isTRUE(cfg$qc_post$plots$ma %||% TRUE)
    do_write   <- isTRUE(cfg$qc_post$outputs$write_de_tables %||% TRUE)

    # Write the summary if available and enabled
    if (do_write && !is.null(de_res$summary_df)) {
        f_all <- file.path(out_qc_post, "de_summary_all.tsv")
        save_tsv_path(de_res$summary_df, f_all)
        files <- c(files, f_all)
    }

    # unified: get tables per contrast (list)
    tables <- get_de_tables_qc_post(de_res, cfg, de_source = de_source)


    # per-contrast plots
    for (cn in names(tables)) {
        de_tbl <- tables[[cn]]
        if (is.null(de_tbl)) next

        # Volcano — emit one variant per p-value type so the report can show both
        if (do_volcano) {
            for (ptype in c("padj", "pval")) {
                p_volcano <- plot_volcano(de_tbl,
                    cfg = cfg,
                    title = paste0("Volcano: ", cn, " (", de_source,
                                   ", ", if (ptype == "padj") "adj.P.Val" else "P.Value", ")"),
                    pvalue_type = ptype
                )
                f_volcano <- file.path(out_qc_post,
                                       sprintf("volcano_%s_%s_%s.png", cn, de_source, ptype))
                ggplot2::ggsave(f_volcano, plot = p_volcano, width = 8, height = 6, dpi = 150)
                plots[[paste0("volcano_", cn, "_", ptype)]] <- p_volcano
                files <- c(files, f_volcano)
            }
            # Backwards-compat alias for downstream consumers (Excel exporter etc.)
            f_volcano_legacy <- file.path(out_qc_post,
                                          sprintf("volcano_%s_%s.png", cn, de_source))
            f_volcano_padj   <- file.path(out_qc_post,
                                          sprintf("volcano_%s_%s_padj.png", cn, de_source))
            if (file.exists(f_volcano_padj) && !file.exists(f_volcano_legacy)) {
                file.copy(f_volcano_padj, f_volcano_legacy, overwrite = TRUE)
                files <- c(files, f_volcano_legacy)
            }
        }

        # MA
        if (do_ma) {
            id_col <- cfg$de_table$id_col %||% "FeatureID"
            de_tbl_ma <- de_tbl
            if (!("AveExpr" %in% colnames(de_tbl_ma))) {
                de_tbl_ma <- add_A_from_expr(de_tbl_ma, pre$expr_imp_single, id_col = "FeatureID")
            }

            p_ma <- plot_ma(de_tbl_ma,
                cfg = cfg,
                title = paste0("MA: ", cn, " (", de_source, ")"),
                use_adj = TRUE
            )
            f_ma <- file.path(out_qc_post, sprintf("ma_%s_%s.png", cn, de_source))
            ggplot2::ggsave(f_ma, plot = p_ma, width = 8, height = 6, dpi = 150)

            plots[[paste0("ma_", cn)]] <- p_ma
            files <- c(files, f_ma)
        }
    }

    # ---------- Top DE heatmaps ----------
    heatmap_sizes <- cfg$qc$top_de_heatmap_sizes %||% c(25, 50, 100)
    for (cn in names(tables)) {
        de_tbl <- tables[[cn]]
        if (is.null(de_tbl)) next
        for (n_top in heatmap_sizes) {
            ph <- plot_top_de_heatmap(de_tbl, pre$expr_imp_single, pre$meta, cfg, n_top, cn)
            if (!is.null(ph)) {
                f_hm <- file.path(out_qc_post, sprintf("heatmap_top%d_%s.png", n_top, cn))
                save_heatmap_to_file(ph, f_hm, width = 1200, height = 800 + n_top * 3, res = 150)
                files <- c(files, f_hm)
                plots[[sprintf("heatmap_top%d_%s", n_top, cn)]] <- ph
            }
        }
    }

    list(plots = plots, files = unique(files))
}


# ---- helpers for QC_post ----

get_de_tables_qc_post <- function(de_res, cfg, de_source = c("summary", "table1"), use_adj = TRUE) {
    de_source <- match.arg(de_source)

    if (de_source == "summary") {
        if (is.null(de_res$summary_df)) stop("QC_post: de_res$summary_df is missing.")
        return(qc_post_tables_from_summary(de_res$summary_df, cfg, use_adj = use_adj))
    }

    # table1: list(contrast -> limma topTable) already, so return as-is
    if (is.null(de_res$runs_de_tables) || length(de_res$runs_de_tables) < 1) {
        stop("QC_post: de_source='table1' requested but de_res$runs_de_tables is missing.")
    }

    tabs <- de_res$runs_de_tables[[1]]
    for (cn in names(tabs)) {
        if (!("FeatureID" %in% colnames(tabs[[cn]]))) {
            if (is.null(rownames(tabs[[cn]]))) stop("QC_post(table1): missing FeatureID and rownames are NULL for contrast ", cn)
            tabs[[cn]]$FeatureID <- rownames(tabs[[cn]])
        }
    }
    return(tabs)
}


# returns named list: contrast -> de_tbl (with logFC + P.Value/adj.P.Val)
qc_post_tables_from_summary <- function(summary_df, cfg, use_adj = TRUE) {
    stopifnot(is.data.frame(summary_df))

    src_id_col <- cfg$de_table$id_col %||% "FeatureID"
    if (!(src_id_col %in% colnames(summary_df))) {
        stop("summary_df missing id col: ", src_id_col)
    }

    contrasts <- sub("^padj\\.imputs\\.", "", grep("^padj\\.imputs\\.", colnames(summary_df), value = TRUE))
    if (length(contrasts) == 0) stop("No contrasts found in summary_df (expected columns like padj.imputs.<contrast>).")

    out <- setNames(vector("list", length(contrasts)), contrasts)

    for (cn in contrasts) {
        p_col <- if (isTRUE(use_adj)) paste0("padj.imputs.", cn) else paste0("pvalue.imputs.", cn)
        fc_col <- paste0("linearFC.imputs.", cn) # this is linear fold-change in your summary
        if (!(p_col %in% colnames(summary_df))) stop("summary_df missing: ", p_col)
        if (!(fc_col %in% colnames(summary_df))) stop("summary_df missing: ", fc_col)

        de_tbl <- data.frame(
            FeatureID = summary_df[[src_id_col]],
            logFC = signed_fc_to_log2(summary_df[[fc_col]]),
            P.Value = as.numeric(summary_df[[paste0("pvalue.imputs.", cn)]]),
            adj.P.Val = as.numeric(summary_df[[paste0("padj.imputs.", cn)]]),
            stringsAsFactors = FALSE
        )

        # optional annotations if exist
        for (a in intersect(c("Protein.Names", "Genes", "First.Protein.Description"), colnames(summary_df))) {
            de_tbl[[a]] <- summary_df[[a]]
        }

        out[[cn]] <- de_tbl
    }

    out
}
