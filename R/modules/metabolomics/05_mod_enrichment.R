# R/modules/metabolomics/05_mod_enrichment.R
#
# Metabolomics enrichment module: runs QEA and ssGSEA, generates plots,
# and saves result tables.
#
# Reuses domain: run_metabolomics_qea, run_metabolomics_ssgsea,
#   plot_enrichment_barplot, plot_ssgsea_boxplots


#' Metabolomics enrichment module
#'
#' Runs QEA (MetaboAnalystR/globaltest) and ssGSEA (GSVA) enrichment analysis.
#' Each is independently optional — if one fails or is unavailable, the other
#' still runs.
#'
#' @param pre     List from preprocess_metabolomics().
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @return list(qea, ssgsea, plots, files) where qea/ssgsea are sub-results.
mod_metabolomics_enrichment <- function(pre, config, out_dir) {
    assert_pre_contract(pre, stage = "metabolomics")

    enr_cfg <- config$modes$metabolomics$enrichment %||% list()

    if (!isTRUE(enr_cfg$run_enrichment)) {
        message("metabolomics enrichment: disabled — skipping")
        return(NULL)
    }

    dirs    <- create_legacy_output_dirs(out_dir)
    out_enr <- dirs$enrichment
    out_ds  <- dirs$datasets

    plots <- list()
    files <- character(0)

    # ---- QEA ----
    qea_res <- tryCatch(
        run_metabolomics_qea(pre, config),
        error = function(e) {
            warning("metabolomics QEA failed: ", e$message)
            NULL
        }
    )

    if (!is.null(qea_res) && !is.null(qea_res$table)) {
        f_qea <- save_tsv(qea_res$table, out_ds, "enrichment_qea_results.tsv")
        files <- c(files, f_qea)

        # QEA barplot
        f_qea_plot <- file.path(out_enr, "enrichment_qea_barplot.png")
        p_qea <- tryCatch({
            pq <- plot_enrichment_barplot(qea_res$table, top_n = 20,
                                           title = "Pathway Enrichment (QEA)")
            if (!is.null(pq)) {
                ggplot2::ggsave(f_qea_plot, pq, width = 12, height = 8, dpi = 300)
            }
            pq
        }, error = function(e) {
            warning("QEA barplot failed: ", e$message)
            NULL
        })
        if (!is.null(p_qea)) {
            files <- c(files, f_qea_plot)
            plots$qea_barplot <- p_qea
        }

        # Per-library barplots
        if (!is.null(qea_res$table$library)) {
            for (lib_name in unique(qea_res$table$library)) {
                lib_df <- qea_res$table[qea_res$table$library == lib_name, ]
                lib_label <- toupper(gsub("_pathway$", "", lib_name))
                f_lib <- file.path(out_enr,
                                    paste0("enrichment_qea_", tolower(lib_label), ".png"))
                p_lib <- tryCatch({
                    pl <- plot_enrichment_barplot(lib_df, top_n = 20,
                                                  title = paste0(lib_label, " Enrichment"))
                    if (!is.null(pl)) {
                        ggplot2::ggsave(f_lib, pl, width = 12, height = 8, dpi = 300)
                    }
                    pl
                }, error = function(e) NULL)
                if (!is.null(p_lib)) {
                    files <- c(files, f_lib)
                    plots[[paste0("qea_", tolower(lib_label))]] <- p_lib
                }
            }
        }
    }

    # ---- ssGSEA ----
    ssgsea_res <- tryCatch(
        run_metabolomics_ssgsea(pre, config),
        error = function(e) {
            warning("metabolomics ssGSEA failed: ", e$message)
            NULL
        }
    )

    if (!is.null(ssgsea_res) && !is.null(ssgsea_res$table)) {
        f_ssgsea <- save_tsv(ssgsea_res$table, out_ds, "enrichment_ssgsea_results.tsv")
        files <- c(files, f_ssgsea)

        # ssGSEA barplot
        f_ssgsea_plot <- file.path(out_enr, "enrichment_ssgsea_barplot.png")
        p_ssgsea <- tryCatch({
            ps <- plot_enrichment_barplot(ssgsea_res$table, top_n = 20,
                                           title = "ssGSEA Pathway Enrichment")
            if (!is.null(ps)) {
                ggplot2::ggsave(f_ssgsea_plot, ps, width = 12, height = 8, dpi = 300)
            }
            ps
        }, error = function(e) {
            warning("ssGSEA barplot failed: ", e$message)
            NULL
        })
        if (!is.null(p_ssgsea)) {
            files <- c(files, f_ssgsea_plot)
            plots$ssgsea_barplot <- p_ssgsea
        }

        # ssGSEA boxplots (if Wilcoxon results available)
        if (!is.null(ssgsea_res$scores) &&
            any(ssgsea_res$table$significant, na.rm = TRUE)) {
            sample_col <- config$modes$metabolomics$effects$samples %||% "sample_id"
            cond_col <- enr_cfg$condition_column %||%
                        config$modes$metabolomics$de$condition_column %||%
                        config$modes$metabolomics$effects$color %||% "sample_type"
            conditions <- factor(pre$meta[[cond_col]][
                match(colnames(ssgsea_res$scores), pre$meta[[sample_col]])
            ])

            f_box <- file.path(out_enr, "enrichment_ssgsea_boxplots.png")
            p_box <- tryCatch({
                n_sig <- sum(ssgsea_res$table$significant, na.rm = TRUE)
                plot_h <- max(6, ceiling(min(n_sig, 12) / 3) * 3.5)
                pb <- plot_ssgsea_boxplots(ssgsea_res$scores, conditions,
                                            ssgsea_res$table)
                if (!is.null(pb)) {
                    ggplot2::ggsave(f_box, pb, width = 12, height = plot_h,
                                     dpi = 300)
                }
                pb
            }, error = function(e) {
                warning("ssGSEA boxplots failed: ", e$message)
                NULL
            })
            if (!is.null(p_box)) {
                files <- c(files, f_box)
                plots$ssgsea_boxplots <- p_box
            }
        }
    }

    list(
        qea    = qea_res,
        ssgsea = ssgsea_res,
        plots  = plots,
        files  = unique(files)
    )
}
