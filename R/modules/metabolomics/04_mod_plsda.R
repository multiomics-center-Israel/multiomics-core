# R/modules/metabolomics/04_mod_plsda.R
#
# Metabolomics PLS-DA module: runs PLS-DA, generates scores + VIP plots,
# and saves VIP tables.
#
# Reuses domain: run_metabolomics_plsda, plot_plsda_scores, plot_plsda_vip


#' Metabolomics PLS-DA module
#'
#' @param pre     List from preprocess_metabolomics().
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @return list(model, vip_scores, explained_variance, vip_df, plots, files)
#'   or NULL.
mod_metabolomics_plsda <- function(pre, config, out_dir) {
    assert_pre_contract(pre, stage = "metabolomics")

    plsda_res <- run_metabolomics_plsda(pre, config)
    if (is.null(plsda_res)) return(NULL)

    dirs   <- create_legacy_output_dirs(out_dir)
    out_qc <- dirs$diagnostic_plots
    out_ds <- dirs$datasets

    plsda_cfg <- config$modes$metabolomics$plsda %||% list()
    colors    <- plsda_cfg$colors %||% NULL
    vip_top_n <- plsda_cfg$vip_top_n %||% 15

    plots <- list()
    files <- character(0)

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

    plsda_res$plots <- plots
    plsda_res$files <- unique(files)
    plsda_res
}
