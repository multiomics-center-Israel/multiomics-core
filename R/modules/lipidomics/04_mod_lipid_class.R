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
        # Wider plot if faceted by bond type
        sat_width <- if ("bond_type" %in% colnames(chain_sat) &&
                         length(unique(chain_sat$bond_type)) > 1) 14 else 8
        p_sat <- tryCatch({
            ps <- plot_chain_saturation(chain_sat)
            ggplot2::ggsave(f_sat_plot, ps, width = sat_width, height = 6, dpi = 300)
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
        # Wider plot if faceted by bond type
        len_width <- if ("bond_type" %in% colnames(chain_len) &&
                         length(unique(chain_len$bond_type)) > 1) 14 else 8
        p_len <- tryCatch({
            pl <- plot_chain_length(chain_len)
            ggplot2::ggsave(f_len_plot, pl, width = len_width, height = 6, dpi = 300)
            pl
        }, error = function(e) NULL)
        if (!is.null(p_len)) {
            files <- c(files, f_len_plot)
            plots$chain_length <- p_len
        }
    }

    # ---- Per-class chain profiles ----
    chain_profiles <- tryCatch(
        compute_per_class_chain_profiles(pre, config),
        error = function(e) {
            warning("Per-class chain profiles failed: ", e$message)
            NULL
        }
    )

    if (!is.null(chain_profiles) && length(chain_profiles) > 0) {
        for (cls in names(chain_profiles)) {
            cls_df <- chain_profiles[[cls]]
            plot_key <- paste0("chain_profile_", cls)
            f_cls_plot <- file.path(out_qc, paste0("chain_profile_", cls, ".png"))
            p_cls <- tryCatch({
                pc <- plot_per_class_chain_profile(cls_df, cls)
                ggplot2::ggsave(f_cls_plot, pc, width = 10, height = 6, dpi = 300)
                pc
            }, error = function(e) NULL)
            if (!is.null(p_cls)) {
                files <- c(files, f_cls_plot)
                plots[[plot_key]] <- p_cls
            }
        }
    }

    # ---- DE-aware structural bubble (per contrast) ----
    if (!is.null(de_res) && !is.null(de_res$summary_df)) {
        cfg_li      <- config$modes$lipidomics %||% list()
        de_cfg      <- cfg_li$de %||% list()
        sig_cutoff  <- de_cfg$p_cutoff %||% 0.05
        contrasts   <- de_cfg$contrasts
        if (is.list(contrasts)) contrasts <- unlist(contrasts)

        if (!is.null(de_cfg$logfc_cutoff)) {
            logfc_cut <- de_cfg$logfc_cutoff
        } else {
            logfc_cut <- log2(de_cfg$linear_fc_cutoff %||% 1.5)
        }

        for (ctr in contrasts %||% character(0)) {
            ctr_label <- make_contrast_label(ctr)
            ctr_tbl   <- tryCatch(
                extract_contrast_table(de_res$summary_df, ctr_label),
                error = function(e) NULL
            )
            if (is.null(ctr_tbl)) next

            p_bub <- tryCatch(
                plot_de_chain_bubble(
                    ctr_tbl, pre$row_data,
                    sig_cutoff   = sig_cutoff,
                    logfc_cutoff = logfc_cut,
                    title = paste0("Significant lipids by chain length × DB: ",
                                   ctr)
                ),
                error = function(e) {
                    warning("plot_de_chain_bubble failed for ", ctr_label,
                            ": ", e$message)
                    NULL
                }
            )
            if (!is.null(p_bub)) {
                f_bub <- file.path(out_qc,
                                   paste0("de_chain_bubble_", ctr_label, ".png"))
                tryCatch(
                    ggplot2::ggsave(f_bub, p_bub, width = 8, height = 6,
                                    dpi = 300),
                    error = function(e) NULL
                )
                files <- c(files, f_bub)
                plots[[paste0("de_chain_bubble_", ctr_label)]] <- p_bub
            }
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
        class_comp     = class_comp,
        chain_sat      = chain_sat,
        chain_len      = chain_len,
        chain_profiles = chain_profiles,
        class_ora      = class_ora,
        plots          = plots,
        files          = unique(files)
    )
}
