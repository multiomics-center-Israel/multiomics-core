# R/modules/metabolomics/02_mod_de.R
#
# Metabolomics DE module: runs differential analysis, generates QC plots
# (volcano, MA, p-value histogram), and saves result tables.
#
# Reuses core plotting: plot_volcano, plot_ma
# Reuses domain: run_metabolomics_de, extract_contrast_table


#' Metabolomics differential expression module
#'
#' @param pre     List from preprocess_metabolomics().
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @return list conforming to the DE contract: summary_df, method, de_tables,
#'   de_model, plots, files
mod_metabolomics_de <- function(pre, config, inputs, out_dir) {
    stage <- "metabolomics"
    assert_pre_contract(pre, stage = stage)

    cfg <- config$modes$metabolomics
    de_cfg <- cfg$de %||% list()

    # Skip if DE is disabled
    if (isTRUE(de_cfg$skip) || identical(de_cfg$enabled, FALSE)) {
        message("metabolomics DE: disabled in config — skipping")
        return(list(summary_df = NULL, method = NULL))
    }

    dirs <- create_legacy_output_dirs(out_dir)
    out_qc <- file.path(dirs$diagnostic_plots, "QC_post")
    ensure_dir(out_qc)
    out_ds <- dirs$datasets

    # ---- Run DE or load pre-computed ----
    de_table_files <- cfg$files$de_table
    has_precomputed <- !is.null(de_table_files) && length(de_table_files) > 0 &&
                       all(nzchar(unlist(de_table_files)))

    if (has_precomputed) {
        message("metabolomics DE: loading pre-computed DE tables")
        de_res <- load_precomputed_metabolomics_de(config)
    } else {
        de_res <- run_metabolomics_de(pre, config, inputs$contrasts)
    }
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

    # ---- Save summary table (with Name column) ----
    summary_out <- de_res$summary_df
    if (!is.null(name_map)) {
        summary_out$Name <- name_map[as.character(summary_out$feature_id)]
        # Move Name to second column (after feature_id)
        cn <- colnames(summary_out)
        summary_out <- summary_out[, c("feature_id", "Name", setdiff(cn, c("feature_id", "Name")))]
    }
    f_summary <- save_tsv(summary_out, out_ds, "de_summary.tsv")
    files <- c(files, f_summary)

    # ---- Per-contrast outputs ----
    # Labels come from de_tables names (set by run_metabolomics_de or precomputed loader)
    contrast_labels <- names(de_res$de_tables)

    for (ctr_label in contrast_labels) {
        ctr_tbl <- extract_contrast_table(de_res$summary_df, ctr_label)

        # Add Name column
        if (!is.null(name_map)) {
            ctr_tbl$Name <- name_map[as.character(ctr_tbl$feature_id)]
            cn <- colnames(ctr_tbl)
            ctr_tbl <- ctr_tbl[, c("feature_id", "Name", setdiff(cn, c("feature_id", "Name")))]
        }

        # Save per-contrast table
        f_ctr <- save_tsv(ctr_tbl, out_ds, paste0("de_", ctr_label, ".tsv"))
        files <- c(files, f_ctr)

        # Volcano plot (using core plot_volcano)
        f_volcano <- file.path(out_qc, paste0("volcano_", ctr_label, ".png"))
        p_volcano <- tryCatch({
            pv <- plot_volcano(ctr_tbl, cfg, title = paste0("Volcano: ", ctr_label))
            ggplot2::ggsave(f_volcano, pv, width = 8, height = 6, dpi = 300)
            pv
        }, error = function(e) {
            warning("Volcano plot failed for ", ctr_label, ": ", e$message)
            NULL
        })
        if (!is.null(p_volcano)) {
            files <- c(files, f_volcano)
            plots[[paste0("volcano_", ctr_label)]] <- p_volcano
        }

        # MA plot (using core plot_ma)
        f_ma <- file.path(out_qc, paste0("MA_", ctr_label, ".png"))
        p_ma <- tryCatch({
            pm <- plot_ma(ctr_tbl, cfg, title = paste0("MA: ", ctr_label))
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
            ph <- plot_pvalue_histogram(ctr_tbl, title = paste0("P-value: ", ctr_label))
            ggplot2::ggsave(f_phist, ph, width = 7, height = 5, dpi = 300)
            ph
        }, error = function(e) {
            warning("P-value histogram failed for ", ctr_label, ": ", e$message)
            NULL
        })
        if (!is.null(p_phist)) {
            files <- c(files, f_phist)
            plots[[paste0("pval_hist_", ctr_label)]] <- p_phist
        }
    }

    # Return full DE result with module outputs
    de_res$plots <- plots
    de_res$files <- unique(files)
    de_res
}


# ---- module-local plot helpers ----------------------------------------------

#' P-value distribution histogram
#'
#' @param de_tbl data.frame with P.Value column.
#' @param title  Plot title.
#' @return ggplot object.
plot_pvalue_histogram <- function(de_tbl, title = "P-value Distribution") {
    stopifnot("P.Value" %in% colnames(de_tbl))

    ggplot2::ggplot(de_tbl, ggplot2::aes(x = P.Value)) +
        ggplot2::geom_histogram(bins = 50, fill = "steelblue",
                                 color = "white", boundary = 0) +
        ggplot2::labs(title = title, x = "Raw p-value", y = "Count") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 13)
        )
}
