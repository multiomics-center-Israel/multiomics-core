# R/modules/lipidomics/09_mod_qc_enhanced.R
#
# Enhanced QC module: PCA loading/scree, PLS-DA CV/confusion, PERMANOVA,
# lipid class correlation, category donut chart.


#' Lipidomics enhanced QC module
#'
#' Generates additional QC plots beyond the basic QC module:
#' PCA loading and scree plots, PLS-DA cross-validation and confusion matrix,
#' lipid class correlation heatmap, category donut chart, and PERMANOVA.
#'
#' @param pre             List from preprocess_lipidomics().
#' @param qc_res          List from mod_lipidomics_qc_pre() (for PCA objects).
#' @param feature_sel_res List from mod_lipidomics_feature_selection() (for PLS-DA).
#' @param config          Full config.
#' @param out_dir         Output directory for this mode.
#' @return list(plots, files, permanova)
mod_lipidomics_qc_enhanced <- function(pre, qc_res, feature_sel_res,
                                        config, out_dir) {
    assert_pre_contract(pre, stage = "lipidomics")

    dirs   <- create_legacy_output_dirs(out_dir)
    out_qc <- dirs$diagnostic_plots
    out_ds <- dirs$datasets
    cfg    <- config$modes$lipidomics

    plots <- list()
    files <- character(0)

    # ---- PCA Loading Plot ----
    pca_obj <- qc_res$objects$norm_log_counts_pca
    if (!is.null(pca_obj)) {
        f_loading <- file.path(out_qc, "pca_loading_plot.png")
        p_loading <- tryCatch({
            pl <- plot_pca_loading(pca_obj, pcs = c(1, 2), top_n = 20)
            ggplot2::ggsave(f_loading, pl, width = 9, height = 8, dpi = 300)
            pl
        }, error = function(e) {
            warning("PCA loading plot failed: ", e$message)
            NULL
        })
        if (!is.null(p_loading)) {
            files <- c(files, f_loading)
            plots$pca_loading <- p_loading
        }

        # ---- PCA Scree Plot ----
        f_scree <- file.path(out_qc, "pca_scree_plot.png")
        p_scree <- tryCatch({
            ps <- plot_pca_scree(pca_obj)
            ggplot2::ggsave(f_scree, ps, width = 8, height = 5, dpi = 300)
            ps
        }, error = function(e) {
            warning("PCA scree plot failed: ", e$message)
            NULL
        })
        if (!is.null(p_scree)) {
            files <- c(files, f_scree)
            plots$pca_scree <- p_scree
        }
    }

    # ---- PLS-DA Cross-Validation ----
    plsda_res <- if (!is.null(feature_sel_res)) feature_sel_res$plsda else NULL
    if (!is.null(plsda_res) && !is.null(plsda_res$X)) {
        cv_result <- tryCatch(
            compute_plsda_cv(plsda_res$X, plsda_res$condition),
            error = function(e) {
                warning("PLS-DA cross-validation failed: ", e$message)
                NULL
            }
        )
        if (!is.null(cv_result)) {
            f_cv <- file.path(out_qc, "plsda_cross_validation.png")
            p_cv <- tryCatch({
                pc <- plot_plsda_cv(cv_result)
                ggplot2::ggsave(f_cv, pc, width = 8, height = 5, dpi = 300)
                pc
            }, error = function(e) {
                warning("PLS-DA CV plot failed: ", e$message)
                NULL
            })
            if (!is.null(p_cv)) {
                files <- c(files, f_cv)
                plots$plsda_cv <- p_cv
            }
        }

        # ---- PLS-DA Confusion Matrix ----
        confusion <- tryCatch(
            compute_plsda_confusion(plsda_res),
            error = function(e) {
                warning("PLS-DA confusion matrix failed: ", e$message)
                NULL
            }
        )
        if (!is.null(confusion)) {
            f_cm <- file.path(out_qc, "plsda_confusion_matrix.png")
            p_cm <- tryCatch({
                pcm <- plot_plsda_confusion(confusion)
                ggplot2::ggsave(f_cm, pcm, width = 6, height = 5, dpi = 300)
                pcm
            }, error = function(e) {
                warning("PLS-DA confusion plot failed: ", e$message)
                NULL
            })
            if (!is.null(p_cm)) {
                files <- c(files, f_cm)
                plots$plsda_confusion <- p_cm
            }
        }
    }

    # ---- Lipid Class Correlation Heatmap ----
    f_class_cor <- file.path(out_qc, "lipid_class_correlation.png")
    p_class_cor <- tryCatch({
        pc <- plot_lipid_class_correlation(pre, config)
        # pheatmap saves directly; check if file was made
        if (!is.null(pc)) {
            grDevices::png(f_class_cor, width = 800, height = 700, res = 150)
            print(pc)
            grDevices::dev.off()
        }
        pc
    }, error = function(e) {
        warning("Lipid class correlation failed: ", e$message)
        NULL
    })
    if (!is.null(p_class_cor) && file.exists(f_class_cor)) {
        files <- c(files, f_class_cor)
        plots$class_correlation <- p_class_cor
    }

    # ---- Category Donut Charts (one per condition) ----
    cat_data <- tryCatch(
        compute_lipid_category_composition(pre, config),
        error = function(e) { NULL }
    )
    if (!is.null(cat_data)) {
        groups <- colnames(cat_data$cat_norm)
        donut_plots <- list()
        for (grp in groups) {
            f_donut_grp <- file.path(out_qc, paste0("lipid_category_donut_", grp, ".png"))
            p_donut_grp <- tryCatch({
                pd <- plot_lipid_category_donut(cat_data, group = grp)
                ggplot2::ggsave(f_donut_grp, pd, width = 9, height = 9, dpi = 300)
                pd
            }, error = function(e) {
                warning("Category donut chart failed for ", grp, ": ", e$message)
                NULL
            })
            if (!is.null(p_donut_grp)) {
                files <- c(files, f_donut_grp)
                donut_plots[[grp]] <- p_donut_grp
            }
        }
        # Store per-group donuts; also keep first as legacy key for template
        plots$category_donuts <- donut_plots
        if (length(donut_plots) > 0) {
            plots$category_donut <- donut_plots[[1]]
        }
    }

    # ---- PERMANOVA ----
    permanova <- tryCatch(
        run_permanova(pre, config),
        error = function(e) {
            warning("PERMANOVA failed: ", e$message)
            NULL
        }
    )
    if (!is.null(permanova)) {
        perm_df <- data.frame(
            statistic = c("F", "R2", "p_value"),
            value = c(permanova$F_statistic, permanova$R2, permanova$p_value),
            stringsAsFactors = FALSE
        )
        f_perm <- save_tsv(perm_df, out_ds, "permanova_result.tsv")
        files <- c(files, f_perm)
    }

    list(
        plots     = plots,
        files     = unique(files),
        permanova = permanova
    )
}
