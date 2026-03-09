# R/modules/lipidomics/07_mod_biomarker.R
#
# Lipidomics biomarker discovery module: ROC, AUC, Cohen's d, power.


#' Lipidomics biomarker discovery module
#'
#' @param pre     List from preprocess_lipidomics().
#' @param de_res  List from mod_lipidomics_de().
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @return list(biomarker_df, roc_curves, mv_roc, plots, files)
mod_lipidomics_biomarker <- function(pre, de_res, config, out_dir) {
    assert_pre_contract(pre, stage = "lipidomics")

    cfg <- config$modes$lipidomics
    bm_cfg <- cfg$biomarker %||% list()

    if (identical(bm_cfg$enabled, FALSE)) {
        message("lipidomics biomarker: disabled in config -- skipping")
        return(NULL)
    }

    dirs   <- create_legacy_output_dirs(out_dir)
    out_qc <- dirs$diagnostic_plots
    out_ds <- dirs$datasets

    plots <- list()
    files <- character(0)

    # ---- Biomarker metrics ----
    biomarker_df <- tryCatch(
        compute_biomarker_metrics(pre, config),
        error = function(e) {
            warning("Biomarker metrics failed: ", e$message)
            NULL
        }
    )

    if (!is.null(biomarker_df)) {
        f_bm <- save_tsv(biomarker_df, out_ds, "biomarker_discovery.tsv")
        files <- c(files, f_bm)

        top_n <- bm_cfg$top_n %||% 20

        # Biomarker table plot
        f_bm_plot <- file.path(out_qc, "biomarker_table.png")
        p_bm <- tryCatch({
            pp <- plot_biomarker_table(biomarker_df, top_n = top_n)
            ggplot2::ggsave(f_bm_plot, pp, width = 10, height = 8, dpi = 300)
            pp
        }, error = function(e) {
            warning("Biomarker table plot failed: ", e$message)
            NULL
        })
        if (!is.null(p_bm)) {
            files <- c(files, f_bm_plot)
            plots$biomarker_table <- p_bm
        }
    }

    # ---- ROC curves (top lipids) ----
    roc_top_n <- bm_cfg$roc_top_n %||% 10
    roc_curves <- tryCatch(
        compute_roc_curves(pre, config, top_n = roc_top_n),
        error = function(e) {
            warning("ROC curves failed: ", e$message)
            NULL
        }
    )

    if (!is.null(roc_curves)) {
        f_roc <- file.path(out_qc, "roc_curves_top_lipids.png")
        p_roc <- tryCatch({
            pr <- plot_roc_curves(roc_curves, top_n = roc_top_n)
            ggplot2::ggsave(f_roc, pr, width = 8, height = 7, dpi = 300)
            pr
        }, error = function(e) {
            warning("ROC curve plot failed: ", e$message)
            NULL
        })
        if (!is.null(p_roc)) {
            files <- c(files, f_roc)
            plots$roc_curves <- p_roc
        }
    }

    # ---- Multivariate ROC (RF-based) ----
    mv_roc <- tryCatch(
        compute_multivariate_roc(pre, config),
        error = function(e) {
            warning("Multivariate ROC failed: ", e$message)
            NULL
        }
    )

    if (!is.null(mv_roc)) {
        f_mv_roc <- file.path(out_qc, "roc_multivariate_rf.png")
        p_mv_roc <- tryCatch({
            pm <- plot_multivariate_roc(mv_roc)
            ggplot2::ggsave(f_mv_roc, pm, width = 7, height = 6, dpi = 300)
            pm
        }, error = function(e) {
            warning("Multivariate ROC plot failed: ", e$message)
            NULL
        })
        if (!is.null(p_mv_roc)) {
            files <- c(files, f_mv_roc)
            plots$roc_multivariate <- p_mv_roc
        }
    }

    list(
        biomarker_df = biomarker_df,
        roc_curves   = roc_curves,
        mv_roc       = mv_roc,
        plots        = plots,
        files        = unique(files)
    )
}
