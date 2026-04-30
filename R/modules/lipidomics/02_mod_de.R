# R/modules/lipidomics/02_mod_de.R
#
# Lipidomics DE module: runs differential analysis, generates volcano/MA plots,
# and saves result tables.


#' Lipidomics differential expression module
#'
#' @param pre     List from preprocess_lipidomics().
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @return list conforming to the DE contract
mod_lipidomics_de <- function(pre, config, out_dir) {
    stage <- "lipidomics"
    assert_pre_contract(pre, stage = stage)

    cfg <- config$modes$lipidomics
    de_cfg <- cfg$de %||% list()

    if (isTRUE(de_cfg$skip) || identical(de_cfg$enabled, FALSE)) {
        message("lipidomics DE: disabled in config — skipping")
        return(list(summary_df = NULL, method = NULL))
    }

    dirs <- create_legacy_output_dirs(out_dir)
    out_qc <- dirs$diagnostic_plots
    out_ds <- dirs$datasets

    de_res <- run_lipidomics_de(pre, config)
    assert_de_contract(de_res, stage = stage)

    plots <- list()
    files <- character(0)

    # ---- Build feature name lookup ----
    name_map <- NULL
    if (!is.null(pre$row_data) && "Name" %in% colnames(pre$row_data)) {
        name_map <- stats::setNames(
            as.character(pre$row_data$Name),
            as.character(pre$row_data$feature_id)
        )
    }

    # ---- Save summary table ----
    summary_out <- de_res$summary_df
    if (!is.null(name_map)) {
        summary_out$Name <- name_map[as.character(summary_out$feature_id)]
        cn <- colnames(summary_out)
        summary_out <- summary_out[, c("feature_id", "Name",
                                        setdiff(cn, c("feature_id", "Name")))]
    }

    # Add lipid class to summary
    if (!is.null(pre$row_data) && "lipid_class" %in% colnames(pre$row_data)) {
        class_map <- stats::setNames(
            as.character(pre$row_data$lipid_class),
            as.character(pre$row_data$feature_id)
        )
        summary_out$lipid_class <- class_map[as.character(summary_out$feature_id)]
        cn <- colnames(summary_out)
        # Move lipid_class right after Name
        name_pos <- which(cn == "Name")
        if (length(name_pos) == 0) name_pos <- 1
        reorder <- c(cn[seq_len(name_pos)], "lipid_class",
                     setdiff(cn[(name_pos + 1):length(cn)], "lipid_class"))
        summary_out <- summary_out[, reorder]
    }

    f_summary <- save_tsv(summary_out, out_ds, "de_summary.tsv")
    files <- c(files, f_summary)

    # ---- Per-contrast outputs ----
    contrasts <- de_cfg$contrasts
    if (is.list(contrasts)) contrasts <- unlist(contrasts)

    for (ctr in contrasts) {
        ctr_label <- make_contrast_label(ctr)
        ctr_tbl <- extract_contrast_table(de_res$summary_df, ctr_label)

        if (!is.null(name_map)) {
            ctr_tbl$Name <- name_map[as.character(ctr_tbl$feature_id)]
            cn <- colnames(ctr_tbl)
            ctr_tbl <- ctr_tbl[, c("feature_id", "Name",
                                    setdiff(cn, c("feature_id", "Name")))]
        }

        f_ctr <- save_tsv(ctr_tbl, out_ds, paste0("de_", ctr_label, ".tsv"))
        files <- c(files, f_ctr)

        # Volcano — emit one variant per p-value type (padj + raw pval)
        for (ptype in c("padj", "pval")) {
            f_volcano <- file.path(out_qc, paste0("volcano_", ctr_label, "_", ptype, ".png"))
            p_volcano <- tryCatch({
                pv <- plot_volcano(ctr_tbl, cfg,
                                   title = paste0("Volcano: ", ctr,
                                                  " (", if (ptype == "padj") "adj.P.Val" else "P.Value", ")"),
                                   pvalue_type = ptype)
                ggplot2::ggsave(f_volcano, pv, width = 8, height = 6, dpi = 300)
                pv
            }, error = function(e) {
                warning("Volcano plot (", ptype, ") failed for ", ctr_label, ": ", e$message)
                NULL
            })
            if (!is.null(p_volcano)) {
                files <- c(files, f_volcano)
                plots[[paste0("volcano_", ctr_label, "_", ptype)]] <- p_volcano
            }
        }
        # Backwards-compat alias
        f_legacy <- file.path(out_qc, paste0("volcano_", ctr_label, ".png"))
        f_padj   <- file.path(out_qc, paste0("volcano_", ctr_label, "_padj.png"))
        if (file.exists(f_padj) && !file.exists(f_legacy)) {
            file.copy(f_padj, f_legacy, overwrite = TRUE)
            files <- c(files, f_legacy)
            plots[[paste0("volcano_", ctr_label)]] <- plots[[paste0("volcano_", ctr_label, "_padj")]]
        }

        # MA
        f_ma <- file.path(out_qc, paste0("MA_", ctr_label, ".png"))
        p_ma <- tryCatch({
            pm <- plot_ma(ctr_tbl, cfg, title = paste0("MA: ", ctr))
            ggplot2::ggsave(f_ma, pm, width = 8, height = 6, dpi = 300)
            pm
        }, error = function(e) {
            warning("MA plot failed for ", ctr_label, ": ", e$message)
            NULL
        })
        if (!is.null(p_ma)) {
            files <- c(files, f_ma)
            plots[[paste0("ma_", ctr_label)]] <- p_ma
        }

        # P-value histogram
        f_phist <- file.path(out_qc, paste0("pvalue_hist_", ctr_label, ".png"))
        p_phist <- tryCatch({
            ph <- plot_pvalue_histogram(ctr_tbl, title = paste0("P-value: ", ctr))
            ggplot2::ggsave(f_phist, ph, width = 7, height = 5, dpi = 300)
            ph
        }, error = function(e) {
            warning("P-value histogram failed: ", e$message)
            NULL
        })
        if (!is.null(p_phist)) {
            files <- c(files, f_phist)
            plots[[paste0("pval_hist_", ctr_label)]] <- p_phist
        }
    }

    de_res$plots <- plots
    de_res$files <- unique(files)
    de_res
}
