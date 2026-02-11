# R/modules/metabolomics/03_mod_rf.R
#
# Metabolomics random forest module: runs RF feature importance,
# generates importance bar plot, and saves result tables.
#
# Reuses domain: run_metabolomics_rf, plot_rf_importance


#' Metabolomics random forest feature importance module
#'
#' @param pre     List from preprocess_metabolomics().
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @return list(importance_df, model, method, plots, files) or NULL.
mod_metabolomics_rf <- function(pre, config, out_dir) {
    assert_pre_contract(pre, stage = "metabolomics")

    rf_res <- run_metabolomics_rf(pre, config)
    if (is.null(rf_res)) return(NULL)

    dirs   <- create_legacy_output_dirs(out_dir)
    out_qc <- dirs$diagnostic_plots
    out_ds <- dirs$datasets

    plots <- list()
    files <- character(0)

    # Save importance table
    f_imp <- save_tsv(rf_res$importance_df, out_ds, "rf_importance.tsv")
    files <- c(files, f_imp)

    # Importance bar plot
    rf_cfg <- config$modes$metabolomics$rf %||% list()
    top_n  <- rf_cfg$top_n %||% 20

    f_plot <- file.path(out_qc, "rf_importance.png")
    p <- tryCatch({
        pp <- plot_rf_importance(rf_res$importance_df, top_n = top_n)
        ggplot2::ggsave(f_plot, pp, width = 8, height = 6, dpi = 300)
        pp
    }, error = function(e) {
        warning("RF plot failed: ", e$message)
        NULL
    })
    if (!is.null(p)) {
        files <- c(files, f_plot)
        plots$rf_importance <- p
    }

    rf_res$plots <- plots
    rf_res$files <- unique(files)
    rf_res
}
