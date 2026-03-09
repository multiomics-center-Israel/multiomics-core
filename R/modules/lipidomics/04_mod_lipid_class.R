# R/modules/lipidomics/04_mod_lipid_class.R
#
# Lipidomics-specific module: lipid class composition, chain profiling,
# and class-level enrichment.


#' Lipidomics class analysis module
#'
#' @param pre     List from preprocess_lipidomics().
#' @param de_res  DE result list (or NULL).
#' @param config  Full config.
#' @param out_dir Output directory for this mode.
#' @return list(class_comp, chain_sat, chain_len, class_ora, plots, files)
mod_lipidomics_class_analysis <- function(pre, de_res = NULL, config, out_dir) {
    assert_pre_contract(pre, stage = "lipidomics")

    dirs   <- create_legacy_output_dirs(out_dir)
    out_qc <- dirs$diagnostic_plots
    out_ds <- dirs$datasets

    plots <- list()
    files <- character(0)

    # ---- Lipid class composition ----
    class_comp <- tryCatch(
        compute_lipid_class_composition(pre, config),
        error = function(e) {
            warning("Lipid class composition failed: ", e$message)
            NULL
        }
    )

    if (!is.null(class_comp)) {
        # Save class composition tables
        abs_df <- as.data.frame(class_comp$class_abs)
        abs_df$lipid_class <- rownames(abs_df)
        f_abs <- save_tsv(abs_df, out_ds, "lipid_class_composition_absolute.tsv")
        files <- c(files, f_abs)

        norm_df <- as.data.frame(class_comp$class_norm)
        norm_df$lipid_class <- rownames(norm_df)
        f_norm <- save_tsv(norm_df, out_ds, "lipid_class_composition_normalized.tsv")
        files <- c(files, f_norm)

        # Plots
        f_abs_plot <- file.path(out_qc, "lipid_class_absolute.png")
        p_abs <- tryCatch({
            pa <- plot_lipid_class_barplot(class_comp, type = "absolute")
            ggplot2::ggsave(f_abs_plot, pa, width = 10, height = 7, dpi = 300)
            pa
        }, error = function(e) NULL)
        if (!is.null(p_abs)) {
            files <- c(files, f_abs_plot)
            plots$class_absolute <- p_abs
        }

        f_norm_plot <- file.path(out_qc, "lipid_class_normalized.png")
        p_norm <- tryCatch({
            pn <- plot_lipid_class_barplot(class_comp, type = "normalized")
            ggplot2::ggsave(f_norm_plot, pn, width = 10, height = 7, dpi = 300)
            pn
        }, error = function(e) NULL)
        if (!is.null(p_norm)) {
            files <- c(files, f_norm_plot)
            plots$class_normalized <- p_norm
        }
    }

    # ---- Chain saturation ----
    chain_sat <- tryCatch(
        compute_chain_saturation(pre, config),
        error = function(e) {
            warning("Chain saturation failed: ", e$message)
            NULL
        }
    )

    if (!is.null(chain_sat)) {
        f_sat <- save_tsv(chain_sat, out_ds, "chain_saturation.tsv")
        files <- c(files, f_sat)

        f_sat_plot <- file.path(out_qc, "chain_saturation.png")
        p_sat <- tryCatch({
            ps <- plot_chain_saturation(chain_sat)
            ggplot2::ggsave(f_sat_plot, ps, width = 8, height = 6, dpi = 300)
            ps
        }, error = function(e) NULL)
        if (!is.null(p_sat)) {
            files <- c(files, f_sat_plot)
            plots$chain_saturation <- p_sat
        }
    }

    # ---- Chain length ----
    chain_len <- tryCatch(
        compute_chain_length_distribution(pre, config),
        error = function(e) {
            warning("Chain length failed: ", e$message)
            NULL
        }
    )

    if (!is.null(chain_len)) {
        f_len <- save_tsv(chain_len, out_ds, "chain_length_distribution.tsv")
        files <- c(files, f_len)

        f_len_plot <- file.path(out_qc, "chain_length_distribution.png")
        p_len <- tryCatch({
            pl <- plot_chain_length(chain_len)
            ggplot2::ggsave(f_len_plot, pl, width = 8, height = 6, dpi = 300)
            pl
        }, error = function(e) NULL)
        if (!is.null(p_len)) {
            files <- c(files, f_len_plot)
            plots$chain_length <- p_len
        }
    }

    # ---- Lipid class ORA ----
    class_ora <- tryCatch(
        lipid_class_ora(de_res, pre$row_data),
        error = function(e) {
            warning("Lipid class ORA failed: ", e$message)
            NULL
        }
    )

    if (!is.null(class_ora)) {
        f_ora <- save_tsv(class_ora, out_ds, "lipid_class_ora.tsv")
        files <- c(files, f_ora)

        f_ora_plot <- file.path(out_qc, "lipid_class_ora.png")
        p_ora <- tryCatch({
            po <- plot_lipid_class_ora(class_ora)
            if (!is.null(po)) {
                ggplot2::ggsave(f_ora_plot, po, width = 10, height = 7, dpi = 300)
            }
            po
        }, error = function(e) NULL)
        if (!is.null(p_ora)) {
            files <- c(files, f_ora_plot)
            plots$class_ora <- p_ora
        }
    }

    list(
        class_comp = class_comp,
        chain_sat  = chain_sat,
        chain_len  = chain_len,
        class_ora  = class_ora,
        plots      = plots,
        files      = unique(files)
    )
}
