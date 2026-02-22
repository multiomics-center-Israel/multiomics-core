# R/modules/metabolomics/03_mod_feature_selection.R
#
# Metabolomics feature selection module: runs RF importance and PLS-DA VIP,
# generates plots, and saves result tables.
#
# Reuses domain: run_metabolomics_rf, plot_rf_importance,
#   run_metabolomics_plsda, plot_plsda_scores, plot_plsda_vip


#' Metabolomics feature selection module
#'
#' Runs RF and PLS-DA (each independently optional) and returns both
#' sub-results.  If one method is disabled or unavailable the other still runs.
#'
#' @param pre     List from preprocess_metabolomics().
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @return list(rf, plsda, plots, files) where rf/plsda may be NULL.
#'   Returns NULL only if both methods are disabled/unavailable.
mod_metabolomics_feature_selection <- function(pre, config, out_dir) {
    assert_pre_contract(pre, stage = "metabolomics")

    dirs   <- create_legacy_output_dirs(out_dir)
    out_qc <- dirs$diagnostic_plots
    out_ds <- dirs$datasets

    plots <- list()
    files <- character(0)

    # ---- Random Forest ----
    rf_res <- tryCatch(
        run_metabolomics_rf(pre, config),
        error = function(e) {
            warning("metabolomics RF failed: ", e$message)
            NULL
        }
    )

    if (!is.null(rf_res)) {
        # Save importance table
        f_imp <- save_tsv(rf_res$importance_df, out_ds, "rf_importance.tsv")
        files <- c(files, f_imp)

        # Importance bar plot
        rf_cfg <- config$modes$metabolomics$rf %||% list()
        top_n  <- rf_cfg$top_n %||% 20

        f_rf_plot <- file.path(out_qc, "rf_importance.png")
        p_rf <- tryCatch({
            pp <- plot_rf_importance(rf_res$importance_df, top_n = top_n)
            ggplot2::ggsave(f_rf_plot, pp, width = 8, height = 6, dpi = 300)
            pp
        }, error = function(e) {
            warning("RF plot failed: ", e$message)
            NULL
        })
        if (!is.null(p_rf)) {
            files <- c(files, f_rf_plot)
            plots$rf_importance <- p_rf
        }
    }

    # ---- PLS-DA ----
    plsda_res <- tryCatch(
        run_metabolomics_plsda(pre, config),
        error = function(e) {
            warning("metabolomics PLS-DA failed: ", e$message)
            NULL
        }
    )

    if (!is.null(plsda_res)) {
        plsda_cfg <- config$modes$metabolomics$plsda %||% list()
        colors    <- plsda_cfg$colors %||% NULL
        vip_top_n <- plsda_cfg$vip_top_n %||% 15

        # Save VIP table
        f_vip <- save_tsv(plsda_res$vip_df, out_ds, "plsda_vip_scores.tsv")
        files <- c(files, f_vip)

        # Scores plot
        f_scores <- file.path(out_qc, "plsda_scores.png")
        p_scores <- tryCatch({
            ps <- plot_plsda_scores(plsda_res, colors = colors)
            ggplot2::ggsave(f_scores, ps, width = 8, height = 7, dpi = 300)
            ps
        }, error = function(e) {
            warning("PLS-DA scores plot failed: ", e$message)
            NULL
        })
        if (!is.null(p_scores)) {
            files <- c(files, f_scores)
            plots$plsda_scores <- p_scores
        }

        # VIP plot
        f_vip_plot <- file.path(out_qc, "plsda_vip_scores.png")
        p_vip <- tryCatch({
            pv <- plot_plsda_vip(plsda_res, top_n = vip_top_n, colors = colors)
            ggplot2::ggsave(f_vip_plot, pv, width = 9, height = 7, dpi = 300)
            pv
        }, error = function(e) {
            warning("PLS-DA VIP plot failed: ", e$message)
            NULL
        })
        if (!is.null(p_vip)) {
            files <- c(files, f_vip_plot)
            plots$plsda_vip <- p_vip
        }
    }

    # Both disabled / unavailable
    if (is.null(rf_res) && is.null(plsda_res)) {
        message("feature selection: both RF and PLS-DA disabled/unavailable — skipping")
        return(NULL)
    }

    list(
        rf    = rf_res,
        plsda = plsda_res,
        plots = plots,
        files = unique(files)
    )
}
