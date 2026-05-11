# R/modules/metabolomics/05_mod_enrichment.R
#
# Metabolomics enrichment module: runs QEA, ssGSEA, ORA, and GSEA,
# generates plots, and saves result tables.
#
# Reuses domain: run_metabolomics_qea, run_metabolomics_ssgsea,
#   run_metabolomics_ora, run_metabolomics_gsea,
#   plot_enrichment_barplot, plot_ssgsea_boxplots,
#   plot_gsea_nes_barplot, plot_ora_dotplot


#' Metabolomics enrichment module
#'
#' Runs up to four enrichment methods (each independently optional):
#'   QEA (globaltest), ssGSEA (GSVA), ORA (Fisher), GSEA (fgsea).
#'
#' @param pre     List from preprocess_metabolomics().
#' @param de_res  DE results (needed by ORA and GSEA; NULL if unavailable).
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @return list(qea, ssgsea, ora, gsea, plots, files).
mod_metabolomics_enrichment <- function(pre, de_res = NULL, config, out_dir) {
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

        # QEA lollipop plot
        f_qea_loll <- file.path(out_enr, "enrichment_qea_lollipop.png")
        p_qea_loll <- tryCatch({
            pl <- plot_qea_lollipop(qea_res$table, top_n = 20)
            if (!is.null(pl)) {
                ggplot2::ggsave(f_qea_loll, pl, width = 12, height = 8, dpi = 300)
            }
            pl
        }, error = function(e) NULL)
        if (!is.null(p_qea_loll)) {
            files <- c(files, f_qea_loll)
            plots$qea_lollipop <- p_qea_loll
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

        # ssGSEA lollipop plot
        f_ssgsea_loll <- file.path(out_enr, "enrichment_ssgsea_lollipop.png")
        p_ssgsea_loll <- tryCatch({
            pl <- plot_ssgsea_lollipop(ssgsea_res$table, top_n = 20)
            if (!is.null(pl)) {
                ggplot2::ggsave(f_ssgsea_loll, pl, width = 12, height = 8, dpi = 300)
            }
            pl
        }, error = function(e) NULL)
        if (!is.null(p_ssgsea_loll)) {
            files <- c(files, f_ssgsea_loll)
            plots$ssgsea_lollipop <- p_ssgsea_loll
        }
    }

    # ---- ORA ----
    ora_res <- tryCatch(
        run_metabolomics_ora(pre, de_res, config),
        error = function(e) {
            warning("metabolomics ORA failed: ", e$message)
            NULL
        }
    )

    if (!is.null(ora_res) && !is.null(ora_res$table)) {
        f_ora <- save_tsv(ora_res$table, out_ds, "enrichment_ora_results.tsv")
        files <- c(files, f_ora)

        # ORA dot plot
        f_ora_plot <- file.path(out_enr, "enrichment_ora_dotplot.png")
        p_ora <- tryCatch({
            po <- plot_ora_dotplot(ora_res$table, top_n = 20,
                                    title = "ORA — Pathway Over-Representation")
            if (!is.null(po)) {
                ggplot2::ggsave(f_ora_plot, po, width = 12, height = 8, dpi = 300)
            }
            po
        }, error = function(e) {
            warning("ORA dotplot failed: ", e$message)
            NULL
        })
        if (!is.null(p_ora)) {
            files <- c(files, f_ora_plot)
            plots$ora_dotplot <- p_ora
        }

        # ORA barplot (reuse enrichment barplot)
        f_ora_bar <- file.path(out_enr, "enrichment_ora_barplot.png")
        p_ora_bar <- tryCatch({
            pb <- plot_enrichment_barplot(ora_res$table, top_n = 20,
                                           title = "ORA — Pathway Enrichment")
            if (!is.null(pb)) {
                ggplot2::ggsave(f_ora_bar, pb, width = 12, height = 8, dpi = 300)
            }
            pb
        }, error = function(e) NULL)
        if (!is.null(p_ora_bar)) {
            files <- c(files, f_ora_bar)
            plots$ora_barplot <- p_ora_bar
        }

        # ORA lollipop plot
        f_ora_loll <- file.path(out_enr, "enrichment_ora_lollipop.png")
        p_ora_loll <- tryCatch({
            pl <- plot_ora_lollipop(ora_res$table, top_n = 20)
            if (!is.null(pl)) {
                ggplot2::ggsave(f_ora_loll, pl, width = 12, height = 8, dpi = 300)
            }
            pl
        }, error = function(e) NULL)
        if (!is.null(p_ora_loll)) {
            files <- c(files, f_ora_loll)
            plots$ora_lollipop <- p_ora_loll
        }
    }

    # ---- GSEA ----
    gsea_res <- tryCatch(
        run_metabolomics_gsea(pre, de_res, config),
        error = function(e) {
            warning("metabolomics GSEA failed: ", e$message)
            NULL
        }
    )

    if (!is.null(gsea_res) && !is.null(gsea_res$table)) {
        f_gsea <- save_tsv(gsea_res$table, out_ds, "enrichment_gsea_results.tsv")
        files <- c(files, f_gsea)

        # GSEA NES barplot
        f_gsea_plot <- file.path(out_enr, "enrichment_gsea_nes_barplot.png")
        p_gsea <- tryCatch({
            pg <- plot_gsea_nes_barplot(gsea_res$table, top_n = 20,
                                         title = "GSEA — Normalized Enrichment Scores")
            if (!is.null(pg)) {
                ggplot2::ggsave(f_gsea_plot, pg, width = 12, height = 8, dpi = 300)
            }
            pg
        }, error = function(e) {
            warning("GSEA NES barplot failed: ", e$message)
            NULL
        })
        if (!is.null(p_gsea)) {
            files <- c(files, f_gsea_plot)
            plots$gsea_nes_barplot <- p_gsea
        }

        # GSEA FDR barplot (reuse enrichment barplot)
        f_gsea_bar <- file.path(out_enr, "enrichment_gsea_barplot.png")
        p_gsea_bar <- tryCatch({
            pb <- plot_enrichment_barplot(gsea_res$table, top_n = 20,
                                           title = "GSEA — Pathway Enrichment")
            if (!is.null(pb)) {
                ggplot2::ggsave(f_gsea_bar, pb, width = 12, height = 8, dpi = 300)
            }
            pb
        }, error = function(e) NULL)
        if (!is.null(p_gsea_bar)) {
            files <- c(files, f_gsea_bar)
            plots$gsea_barplot <- p_gsea_bar
        }

        # GSEA lollipop plot
        f_gsea_loll <- file.path(out_enr, "enrichment_gsea_lollipop.png")
        p_gsea_loll <- tryCatch({
            pl <- plot_gsea_lollipop(gsea_res$table, top_n = 20)
            if (!is.null(pl)) {
                ggplot2::ggsave(f_gsea_loll, pl, width = 12, height = 8, dpi = 300)
            }
            pl
        }, error = function(e) NULL)
        if (!is.null(p_gsea_loll)) {
            files <- c(files, f_gsea_loll)
            plots$gsea_lollipop <- p_gsea_loll
        }
    }

    # ---- Mummichog ----
    mummichog_res <- NULL
    mummi_cfg <- enr_cfg$mummichog %||% list()
    if (isTRUE(mummi_cfg$enabled)) {
        mummichog_res <- tryCatch({
            mummi_out <- run_mummichog_all(pre, de_res, config, out_dir)
            if (!is.null(mummi_out) && length(mummi_out$results) > 0) {
                out_mummi <- file.path(out_dir, "mummichog_enrichment")

                for (org_label in names(mummi_out$results)) {
                    res <- mummi_out$results[[org_label]]
                    if (is.null(res) || is.null(res$table)) next

                    # Save TSV
                    f_tsv <- save_tsv(res$table, out_mummi,
                                      paste0("mummichog_", org_label, "_results.tsv"))
                    files <- c(files, f_tsv)

                    # Barplot
                    f_bar <- file.path(out_mummi,
                                       paste0("mummichog_", org_label, "_barplot.png"))
                    p_bar <- tryCatch({
                        pb <- plot_mummichog_barplot(
                            res$table, top_n = 20,
                            title = paste0("Mummichog — ", toupper(org_label)))
                        if (!is.null(pb))
                            ggplot2::ggsave(f_bar, pb, width = 12, height = 8, dpi = 300)
                        pb
                    }, error = function(e) NULL)
                    if (!is.null(p_bar)) {
                        files <- c(files, f_bar)
                        plots[[paste0("mummichog_", org_label, "_barplot")]] <- p_bar
                    }

                    # Lollipop
                    f_loll <- file.path(out_mummi,
                                        paste0("mummichog_", org_label, "_lollipop.png"))
                    p_loll <- tryCatch({
                        pl <- plot_mummichog_lollipop(
                            res$table, top_n = 20,
                            title = paste0("Mummichog — ", toupper(org_label)))
                        if (!is.null(pl))
                            ggplot2::ggsave(f_loll, pl, width = 12, height = 8, dpi = 300)
                        pl
                    }, error = function(e) NULL)
                    if (!is.null(p_loll)) {
                        files <- c(files, f_loll)
                        plots[[paste0("mummichog_", org_label, "_lollipop")]] <- p_loll
                    }
                }

                # Cross-organism heatmap
                if (length(mummi_out$results) > 1) {
                    f_heat <- file.path(out_mummi, "mummichog_organism_heatmap.png")
                    p_heat <- tryCatch({
                        ph <- plot_mummichog_organism_heatmap(mummi_out$results)
                        if (!is.null(ph))
                            ggplot2::ggsave(f_heat, ph, width = 14, height = 10, dpi = 300)
                        ph
                    }, error = function(e) NULL)
                    if (!is.null(p_heat)) {
                        files <- c(files, f_heat)
                        plots$mummichog_organism_heatmap <- p_heat
                    }
                }
            }
            mummi_out
        }, error = function(e) {
            warning("mummichog enrichment failed: ", e$message)
            NULL
        })
    }

    list(
        qea       = qea_res,
        ssgsea    = ssgsea_res,
        ora       = ora_res,
        gsea      = gsea_res,
        mummichog = mummichog_res,
        plots     = plots,
        files     = unique(files)
    )
}
